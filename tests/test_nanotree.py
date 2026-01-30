# Fast root file processor tests (playground)
# 
# m.mieskolainen@imperial.ac.uk, 2022

import sys
import pprint
import pickle
import uproot
import time
import awkward as ak
import numpy as np

sys.path.append(".")
from icenet.tools import iceroot
from icenet.tools import io
from icenet.algo import analytic

#path     = '/home/user/travis-stash/input/icedqcd'
#path     = '/vols/cms/mc3909'
path = '/vols/cms/khl216/bparkProductionAll_V1p3'

#datasets = 'bparkProductionV2/HiddenValley_vector_m_10_ctau_10_xiO_1_xiL_1_privateMC_11X_NANOAODSIM_v2_generationForBParking/output_1.root'
datasets = 'hiddenValleyGridPack_vector_m_2_ctau_10_xiO_1_xiL_1/data_1.root'
key      = 'Events;1'

ids = ['nsv', 'sv_.*', 'cpf_.*', 'Jet_.*', 'muonSV_.*']

#ids = ['nsv', 'sv_deltaR', 'sv_mass']
#ids = ['sv_ptrel', 'sv_deta', 'sv_dphi', 'sv_deltaR', 'sv_mass', 'sv_chi2']
#ids = None

entry_stop = 1000
rootfile   = io.glob_expand_files(datasets=datasets, datapath=path)

for library in ['ak']:

    t       = time.time()
    out,ids = iceroot.load_tree(rootfile=rootfile, tree=key, entry_stop=entry_stop, ids=ids, library=library)
    elapsed = time.time() - t


    #print(ak.fields(out.nsv))
    #count   = (np.isfinite(out['sv']['dxysig']))

    #print(count)
    #exit()
 
    print(f'Opening with <{library}> took: {elapsed} sec')
	
    idx = out.sv.chi2 > 0.0
    print(idx)
    out['sv'] = out['sv'][idx]
 	
    print(len(out))
    print(out)

    N = 20

    # Add in a new record
    out['MODEL_m'] = ak.Array(np.random.rand(len(out)))
    print(len(out))
    print(out)
	
    #out = out[:, out.cpf.jetIdx > 0]

    print(ak.fields(out))
	
    X = out

    def loop(X):
        for i in range(N):
            nsv = X[i].nsv
            print(f'\n[Event {i}]:')
            print(f'ak:')
            print(f'X[{i}]: {X[i]}')
            print(f'MODEL_m: {X[i].MODEL_m}')


            jetIdx =  X[i].cpf.jetIdx
            mask   = (X[i].cpf.jetIdx != -1)

            print(mask)

            # Alternative indexing
            print(X[i]['cpf']['jetIdx'])
			
            p = X[i].cpf[mask]
            print(f'  px: {p.px} py: {p.py} pz: {p.pz} trackSip2dVal: {p.trackSip2dVal} trackSip3dVal: {p.trackSip3dVal}')

            # print(f'list:')
            # print(f'{X_list[i]}')

            for j in range(nsv):
	            print(f'sv[{j}].mass: {X[i].sv[j].mass}')
                #print(f'{X_list[i]}')
        

    print(f'---------------------------------------------------')
	
	
    print(f'Testing ak under different representations')

    X    = out[:N]
    print(f'{len(out)}')
	
	#X_list  = out[:N].tolist()
    loop(X)    

    muonSV_phi = np.sign(X['muonSV']['y'])*np.arccos(X['muonSV']['x']/np.sqrt((X['muonSV']['x'])**2 + (X['muonSV']['y'])**2))
    muonSV_theta = np.arccos(X['muonSV']['z']/np.sqrt((X['muonSV']['x'])**2 + (X['muonSV']['y'])**2+(X['muonSV']['z'])**2))
    muonSV_eta = -np.log(np.tan(muonSV_theta/2))

    muonSV_dphi = ak.fill_none(ak.pad_none(muonSV_phi[:,1:2],1),0)-ak.fill_none(ak.pad_none(muonSV_phi[:,0:1],1),0)
    muonSV_dphi = analytic.phi_phasewrap(muonSV_dphi)
    muonSV_deta = ak.fill_none(ak.pad_none(muonSV_eta[:,1:2],1),0)-ak.fill_none(ak.pad_none(muonSV_eta[:,0:1],1),0)

    muonSV_deltaR = np.sqrt(muonSV_dphi**2 + muonSV_deta**2)

    print("muonSV_phi = ", muonSV_phi)
    print("muonSV_deltaR = ", muonSV_deltaR)
    print("ak.firsts(muonSV_phi) = ", ak.Array(ak.firsts(muonSV_phi)))
    print("muonSV_phi[...,:1] = ", muonSV_phi[...,:1])
    print("muonSV_phi[...,:2] = ", muonSV_phi[...,:2])
    print("muonSV_phi[...,0:1] = ", muonSV_phi[...,0:1])
    print("muonSV_phi[...,1:2] = ", muonSV_phi[...,1:2])

    muonSV_phi_constructed = ak.Array([[1.23], [0.531, -3.14], [0.1,0.2,0.3]])
    muonSV_eta_constructed = ak.Array([[1.01], [0.31, 0.052], [0.2,0.3,0.4]])

    muonSV_dphi = []
    muonSV_deta = []
    for i in range(len(muonSV_phi)):
        print("i = ", i)
        dphi = []
        deta = []
        for j in range(len(muonSV_phi[i])):
            print("j = ", j)
            if j == 0:
                dphi.append(muonSV_phi[i][j])
                deta.append(muonSV_eta[i][j])
            else:
                dphi.append(muonSV_phi[i][j] - muonSV_phi[i][0])
                deta.append(muonSV_eta[i][j] - muonSV_eta[i][0])
        muonSV_dphi.append(dphi)
        muonSV_deta.append(deta)

    muonSV_dphi = ak.Array(muonSV_dphi)
    muonSV_dphi = analytic.phi_phasewrap(muonSV_dphi)
    muonSV_deta = ak.Array(muonSV_deta)

    print("muonSV_dphi = ", muonSV_dphi)
    print("muonSV_deta = ", muonSV_deta)
    print("muonSV_deltaR = ", np.sqrt(muonSV_dphi**2 + muonSV_deta**2))

    def dRmSV_compute(x, y, z):
  
        N = len(x)
        deltaR_values = []

        phi = np.sign(y)*np.arccos(x/np.sqrt(x**2 + y**2))
        theta = np.arccos(z/np.sqrt(x**2 + y**2+ z**2))
        eta = -np.log(np.tan(theta/2))

        for i in range(N):
            for j in range(i+1, N):
                dphi = analytic.phi_phasewrap(phi[j] - phi[i])
                deta = eta[j] - eta[i]

                deltaR_value = np.sqrt(dphi**2 + deta**2)

                deltaR_values.append(deltaR_value)
  
        return np.array(deltaR_values)

    dtype = np.float32
    N_max = 8
    null_value = -999.0
    num_of_events = N

    dRmSV = np.full((int(num_of_events), int((N_max**2 - N_max)/2)), null_value, dtype=dtype)

    for i in range(num_of_events):
        vec = dRmSV_compute(X['muonSV']['x'][i], X['muonSV']['y'][i], X['muonSV']['z'][i])
        print("vec = ", vec)
        vec = np.sort(vec)
        print("Sorted vec = ", vec)
        dRmSV[i, 0:len(vec)] = vec 
        
    print("dRmSV = ", dRmSV)

    

    '''
    muonSV_phi_0_constructed = ak.Array([[1.23], [0.531], [0.1]])
    print("muonSV_phi_constructed - muonSV_phi_0_constructed = ", muonSV_phi_constructed - muonSV_phi_0_constructed)
    '''

    muonSV_dphi_constructed = ak.fill_none(ak.pad_none(muonSV_phi_constructed[:,1:2],1),0)-ak.fill_none(ak.pad_none(muonSV_phi_constructed[:,0:1],1),0)
    muonSV_dphi_constructed = analytic.phi_phasewrap(muonSV_dphi_constructed)
    muonSV_deta_constructed = ak.fill_none(ak.pad_none(muonSV_eta_constructed[...,1:2],1),0)-ak.fill_none(ak.pad_none(muonSV_eta_constructed[...,0:1],1),0)

    print("muonSV_phi_constructed[,0:1] = ", muonSV_phi_constructed[:,0:1])
    print("muonSV_phi_constructed[,1:2] = ", muonSV_phi_constructed[:,1:2])
    print("muonSV_phi_constructed[...,1:2]-muonSV_phi_constructed[...,0:1] = ", ak.fill_none(ak.pad_none(muonSV_phi_constructed[...,1:2],1),0)-ak.fill_none(ak.pad_none(muonSV_phi_constructed[...,0:1],1),0))
    print("muonSV_eta_constructed[...,1:2]-muonSV_eta_constructed[...,0:1] = ", ak.fill_none(ak.pad_none(muonSV_eta_constructed[...,1:2],1),0)-ak.fill_none(ak.pad_none(muonSV_eta_constructed[...,0:1],1),0))
    print("muonSV_dR = ", np.sqrt(muonSV_dphi_constructed**2 + muonSV_deta_constructed**2))

    muonSV_phi_0 = np.sign(X['muonSV', 'y'][0])*np.arccos(X['muonSV', 'y'][0]/np.sqrt((X['muonSV', 'x'][0])**2 + (X['muonSV', 'y'][0])**2))
    muonSV_theta_0 = np.arccos(X['muonSV', 'z'][0]/np.sqrt((X['muonSV', 'x'][0])**2 + (X['muonSV', 'y'][0])**2+(X['muonSV', 'z'][0])**2))
    muonSV_eta_0 = -np.log(np.tan(muonSV_theta_0/2))

    muonSV_dphi = analytic.phi_phasewrap(muonSV_phi - muonSV_phi_0)
    muonSV_deta = muonSV_eta - muonSV_eta_0
    deltaR_muonSV = muonSV_dphi**2 + muonSV_deta**2

    print('muonSV_phi_0 = ', muonSV_phi_0)
    print('deltaR_muonSV = ', deltaR_muonSV)

    #muonSV_phi = np.sign(data.x['muonSV', 'y'])*np.arccos(data.x['muonSV', 'y']/np.sqrt((data.x['muonSV', 'x'])**2 + (data.x['muonSV', 'y'])**2))
    #muonSV_theta = np.arccos(data.x['muonSV', 'z']/np.sqrt((data.x['muonSV', 'x'])**2 + (data.x['muonSV', 'y'])**2+(data.x['muonSV', 'z'])**2))
    #muonSV_eta = -np.log(np.tan(muonSV_theta/2))

    #muonSV_phi_0 = np.sign(data.x['muonSV', 'y'][0])*np.arccos(data.x['muonSV', 'y'][0]/np.sqrt((data.x['muonSV', 'x'][0])**2 + (data.x['muonSV', 'y'][0])**2))
	#muonSV_theta_0 = np.arccos(data.x['muonSV', 'z'][0]/np.sqrt((data.x['muonSV', 'x'][0])**2 + (data.x['muonSV', 'y'][0])**2+(data.x['muonSV', 'z'][0])**2))
	#muonSV_eta_0 = -np.log(np.tan(muonSV_theta_0/2))

    #Add dR between the first two muon SVs for events with at least two muon SVs
    '''
    if data.x['nmuonSV'] >= 2:
        muonSV_phi_0 = np.sign(data.x['muonSV', 'y'][0])*np.arccos(data.x['muonSV', 'y'][0]/np.sqrt((data.x['muonSV', 'x'][0])**2 + (data.x['muonSV', 'y'][0])**2))
        muonSV_theta_0 = np.arccos(data.x['muonSV', 'z'][0]/np.sqrt((data.x['muonSV', 'x'][0])**2 + (data.x['muonSV', 'y'][0])**2+(data.x['muonSV', 'z'][0])**2))
        muonSV_eta_0 = -np.log(np.tan(muonSV_theta_0/2))
        muonSV_phi_1 = np.sign(data.x['muonSV', 'y'][1])*np.arccos(data.x['muonSV', 'y'][1]/np.sqrt((data.x['muonSV', 'x'][1])**2 + (data.x['muonSV', 'y'][1])**2))
        muonSV_theta_1 = np.arccos(data.x['muonSV', 'z'][1]/np.sqrt((data.x['muonSV', 'x'][1])**2 + (data.x['muonSV', 'y'][1])**2+(data.x['muonSV', 'z'][1])**2))
        muonSV_eta_1 = -np.log(np.tan(muonSV_theta_1/2))
        data.x['deltaR_muonSV'] = np.sqrt((muonSV_eta_0 - muonSV_eta_1)**2 + (analytic.phi_phasewrap(muonSV_phi_0 - muonSV_phi_1))**2)
    '''
    '''
    muonSV_phi = np.sign(data.x['muonSV', 'y'])*np.arccos(data.x['muonSV', 'y']/np.sqrt((data.x['muonSV', 'x'])**2 + (data.x['muonSV', 'y'])**2))
    muonSV_theta = np.arccos(data.x['muonSV', 'z']/np.sqrt((data.x['muonSV', 'x'])**2 + (data.x['muonSV', 'y'])**2+(data.x['muonSV', 'z'])**2))
    muonSV_eta = -np.log(np.tan(muonSV_theta/2))

    muonSV_phi_0 = np.sign(data.x['muonSV', 'y'][0])*np.arccos(data.x['muonSV', 'y']/np.sqrt((data.x['muonSV', 'x'])**2 + (data.x['muonSV', 'y'])**2))
    '''

    print(f'---------------------------------------------------')
    print(f'Testing nested object property based selection')

    # Put a requirement on the objects
    cut_str = ['X.nsv >= 1',
	           'ak.sum(X.sv.dxysig >= 5,        -1)',
	           'ak.sum(X.Jet.pt    > 40.0,      -1)',
	           'ak.sum(np.abs(X.Jet.eta) < 2.0, -1)',
	           'ak.sum(np.logical_and(X.Jet.pt  > 150.0, np.abs(X.Jet.eta)  < 2.4), -1) > 0']

    cuts = []
    for i in range(len(cut_str)):
        cuts.append(eval(cut_str[i]))
        print(f'cuts[{i}] = {cuts[i]}')
	
