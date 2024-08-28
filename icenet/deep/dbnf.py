# Block Neural Autoregressive Flow (BNAF)
# Generative (Neural Likelihood) Functions
#
# https://arxiv.org/abs/1904.04676
# https://github.com/nicola-decao/BNAF (MIT license)
#
#
# m.mieskolainen@imperial.ac.uk, 2024


import os
import torch
import torch.nn.functional as F

import numpy as np
from tqdm import tqdm
from torch.utils import data


from icenet.deep  import optimize

from . bnaf import *
from icenet.tools import aux
from icenet.tools import aux_torch


def compute_log_p_x(model, x):
    """ Model log-density value log[pdf(x), where x is the data vector]
    
    log p(x) = log p(z) + sum_{k=1}^K log|det J_{f_k}|
    
    Args:
        model : model object
        x     : N minibatch vectors
    Returns:
        log-likelihood value
    """

    # Evaluate the non-diagonal and the diagonal part
    y, log_diag = model(x)

    # Sum of log-likelihoods (log product) pushed through the flow, evaluated for the unit Gaussian
    log_p_y = torch.distributions.Normal(torch.zeros_like(y), torch.ones_like(y)).log_prob(y).sum(dim=-1)
    
    return log_p_y + log_diag


def get_pdf(model, x) :
    """ Evaluate learned density (pdf) at point x
    
    Args:
        model :  model object
        x     :  input vector(s)
    Returns:
        pdf value

    Examples:
        > x = torch.tensor([[1.0, 2.0]])
        > l = get_pdf(model,x)
    """
    return (torch.exp(compute_log_p_x(model=model, x=x))).detach().cpu().numpy()


def predict(X, models, return_prob=True, EPS=1E-9):
    """
    2-class density ratio pdf(x,S) / pdf(x,B) for each vector x.
    
    Args:
        param       : input parameters
        X           : pytorch tensor of vectors
        models      : list of model objects
        return_prob : return pdf(S) / (pdf(S) + pdf(B)), else pdf(S) / pdf(B)
    
    Returns:
        likelihood ratio (or alternatively probability)
    """
    
    print(__name__ + f'.predict: Computing density (likelihood) ratio for N = {X.shape[0]} events | return_prob = {return_prob}')
    
    bgk_pdf = get_pdf(models[0], X)
    sgn_pdf = get_pdf(models[1], X)
    
    if return_prob:
        out = sgn_pdf / np.clip(sgn_pdf + bgk_pdf, EPS, None)
    else:
        out = sgn_pdf / np.clip(bgk_pdf, EPS, None)
    
    out[~np.isfinite(out)] = 0
    
    return out


class Dataset(torch.utils.data.Dataset):

    def __init__(self, X, W):
        """ Initialization """
        self.X = X
        self.W = W

    def __len__(self):
        """ Return the total number of samples """
        return self.X.shape[0]

    def __getitem__(self, index):
        """ Generates one sample of data """
        # Use ellipsis ... to index over scalar [,:], vector [,:,:], tensor [,:,:,..,:] indices
        return self.X[index,...], self.W[index]


def train(model, optimizer, scheduler, trn_x, val_x,
          trn_weights, val_weights, param, modeldir, save_name):
    """ Train the model density.
    
    Args:
        model       : initialized model object
        optimizer   : optimizer object
        scheduler   : optimization scheduler
        trn_x       : training vectors
        val_x       : validation vectors
        trn_weights : training weights
        val_weights : validation weights
        param       : parameters
        modeldir    : directory to save the model
    """
    
    model, device = optimize.model_to_cuda(model, param['device'])

    # TensorboardX
    if 'tensorboard' in param and param['tensorboard']:
        from tensorboardX import SummaryWriter
        writer = SummaryWriter(os.path.join('tmp/tensorboard/', save_name))

    if trn_weights is None:
        trn_weights = torch.ones(trn_x.shape[0], dtype=torch.float32)
    
    if val_weights is None:
        val_weights = torch.ones(val_x.shape[0], dtype=torch.float32)
    
    # Datasets
    training_set   = Dataset(trn_x, trn_weights)
    validation_set = Dataset(val_x, val_weights)
    
    # N.B. We use 'sampler' with 'BatchSampler', which loads a set of events using multiple event indices (faster) than the default
    # one which takes events one-by-one and concatenates the results (slow).
    params_train = {'batch_size'  : None,
                    'num_workers' : param['num_workers'],
                    'sampler'     : torch.utils.data.BatchSampler(
                        torch.utils.data.RandomSampler(training_set), param['opt_param']['batch_size'], drop_last=False
                    ),
                    'pin_memory'  : True}
    
    params_validate = {'batch_size'  : None,
                    'num_workers' : param['num_workers'],
                    'sampler'     : torch.utils.data.BatchSampler(
                        torch.utils.data.RandomSampler(validation_set), param['eval_batch_size'], drop_last=False
                    ),
                    'pin_memory'  : True}
    
    training_loader   = torch.utils.data.DataLoader(training_set,   **params_train)
    validation_loader = torch.utils.data.DataLoader(validation_set, **params_validate)
    
    # Loss function
    """
    Note:
        log-likelihood functions can be weighted linearly, due to
        \\prod_i p_i(\\theta; x_i)**w_i ==\\log==> \\sum_i w_i \\log p_i(\\theta; x_i)
    """
    def lossfunc(model, x, weights):
        w = weights / torch.sum(weights, dim=0)
        lossvec = compute_log_p_x(model, x)
        return -(lossvec * w).sum(dim=0) # log-likelihood
    
    trn_losses = []
    val_losses = []
    
    # Training loop
    for epoch in tqdm(range(param['opt_param']['start_epoch'], param['opt_param']['start_epoch'] + param['opt_param']['epochs']), ncols = 88):
        
        model.train() # !

        train_loss  = []

        for batch_x, batch_weights in training_loader:

            # Transfer to GPU
            batch_x       = batch_x.to(device, dtype=torch.float32, non_blocking=True)
            batch_weights = batch_weights.to(device, dtype=torch.float32, non_blocking=True)

            loss = lossfunc(model=model, x=batch_x, weights=batch_weights)
            
            # Zero gradients, calculate loss, calculate gradients and update parameters
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=param['opt_param']['clip_norm'])
            optimizer.step()

            train_loss.append(loss)

        train_loss = torch.stack(train_loss).mean()
        optimizer.swap()

        # ----------------------------------------------------------------
        # Compute validation loss
        
        if epoch == 0 or (epoch % param['evalmode']) == 0:
        
            model.eval() # !

            validation_loss = []
            
            with torch.no_grad():
            
                for batch_x, batch_weights in validation_loader:

                    # Transfer to GPU
                    batch_x       = batch_x.to(device, dtype=torch.float32, non_blocking=True)
                    batch_weights = batch_weights.to(device, dtype=torch.float32, non_blocking=True)
                    
                    loss = lossfunc(model=model, x=batch_x, weights=batch_weights)
                    validation_loss.append(loss)
            
            validation_loss = torch.stack(validation_loss).mean()
            optimizer.swap()
        
        # ----------------------------------------------------------------
        
        # Save metrics
        trn_losses.append(train_loss.item())
        val_losses.append(validation_loss.item())
        
        # Step scheduler
        stop = scheduler.step(validation_loss)
        
        print('Epoch {:3d}/{:3d} | Train: loss: {:4.3f} | Validation: loss: {:4.3f} | lr: {:0.3E}'.format(
            epoch,
            param['opt_param']['start_epoch'] + param['opt_param']['epochs'],
            train_loss.item(),
            validation_loss.item(),
            scheduler.get_last_lr()[0])
        )
        
        # Save
        filename = f'{modeldir}/{save_name}_{epoch}.pth'
        aux_torch.save_torch_model(model     = model,
                                   optimizer = optimizer,
                                   epoch     = epoch,
                                   losses    = {'trn_losses': trn_losses, 'val_losses': val_losses},
                                   filename  = filename)()

        if 'tensorboard' in param and param['tensorboard']:
            writer.add_scalar('lr', scheduler.get_last_lr()[0], epoch)
            writer.add_scalar('loss/validation', validation_loss.item(), epoch)
            writer.add_scalar('loss/train', train_loss.item(), epoch)
        
        if stop:
            break
        

def create_model(param, verbose=False, rngseed=0):
    """ Construct the network object.
    
    Args:
        param : parameters
    Returns:
        model : model object
    """

    # For random permutations
    np.random.seed(rngseed)
    
    flows = []
    for f in range(param['flows']):

        # Layers of the flow block
        layers = []
        for _ in range(param['layers']):
            layers.append(MaskedWeight(param['n_dims'] * param['hidden_dim'],
                                       param['n_dims'] * param['hidden_dim'], dim=param['n_dims']))
            layers.append(Tanh())
        
        # Add this flow block
        flows.append(
            BNAF(*([MaskedWeight(param['n_dims'], param['n_dims'] * param['hidden_dim'], dim=param['n_dims']), Tanh()] + \
                   layers + \
                   [MaskedWeight(param['n_dims'] * param['hidden_dim'], param['n_dims'], dim=param['n_dims'])]),\
                 res=param['residual'] if f < param['flows'] - 1 else None
            )
        )

        # Flow permutations
        if f < param['flows'] - 1:
            flows.append(Permutation(param['n_dims'], param['perm']))

    # Create the model
    model  = Sequential(*flows)
    
    params = sum((p != 0).sum() if len(p.shape) > 1 else torch.tensor(p.shape).item()
                 for p in model.parameters()).item()

    # Print model information    
    print('{}'.format(model))
    print('Parameters={}, n_dims={}'.format(params, param['n_dims']))

    return model


def load_models(param, modelnames, modeldir, device):
    """ Load models from files
    """
    
    models = []
    for i in range(len(modelnames)):
        print(__name__ + f'.load_models: Loading model[{i}] from {modelnames[i]}')

        model      = create_model(param['model_param'], verbose=False)
        param['opt_param']['start_epoch'] = 0

        filename   = aux.create_model_filename(path=modeldir, label=modelnames[i], \
            epoch=param['readmode'], filetype='.pth')
        checkpoint = torch.load(filename, map_location='cpu')
        
        model.load_state_dict(checkpoint['model'])
        model, device = optimize.model_to_cuda(model, device_type=device)
        
        model.eval() # Turn on eval mode!

        models.append(model)

    return models, device
