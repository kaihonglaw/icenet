#!/bin/sh
#
# Generate dynamic YAML files
#
# Use * or other glob wildcards for filenames
#
# Run with: source runme.sh

#source $HOME/setconda.sh
__conda_setup="$('/vols/cms/khl216/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/vols/cms/khl216/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/vols/cms/khl216/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/vols/cms/khl216/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
conda activate icenet

ICEPATH="/vols/cms/khl216/icenet_new_models/icenet"
cd $ICEPATH
echo "$(pwd)"
source $ICEPATH/setenv.sh

# Training
#python configs/dqcd/include/ymlgen.py --process 'QCD'    --filerange '[0-50]'      --outputfile configs/dqcd/include/QCD_newmodels.yml
python configs/dqcd/include/ymlgen.py --process 'scenarioA' --filerange '[0-20]'      --outputfile configs/dqcd/include/scenarioA_all_model_points.yml
python configs/dqcd/include/ymlgen.py --process 'data-D' --filerange '[0-215]'      --outputfile configs/dqcd/include/data-D_newmodels.yml 

# Deployment
#python configs/dqcd/include/ymlgen.py --process 'QCD'    --filerange '[51-100000]' --outputfile configs/dqcd/include/QCD_newmodels_deploy.yml
python configs/dqcd/include/ymlgen.py --process 'scenarioA' --filerange '[21-100000]' --outputfile configs/dqcd/include/scenarioA_all_model_points_deploy.yml
python configs/dqcd/include/ymlgen.py --process 'data-D' --filerange '[216-1074]'      --outputfile configs/dqcd/include/data-D_newmodels_deploy.yml



