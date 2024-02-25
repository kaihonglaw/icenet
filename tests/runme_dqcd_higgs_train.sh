#!/bin/sh
#
# Execute training and evaluation for the DQCD analysis
#
# Run with: source runme.sh

source $HOME/setconda.sh
conda activate icenet

ICEPATH="/vols/cms/mmieskol/icenet"
cd $ICEPATH
echo "$(pwd)"
source $ICEPATH/setenv.sh

CONFIG="tune0.yml"
#DATAPATH="/home/user/travis-stash/input/icedqcd"
DATAPATH="/vols/cms/mc3909"

CONDITIONAL=1

#mkdir "figs/dqcd/config-[$CONFIG]" -p # for output ascii dump

if [ ${maxevents+x} ]; then MAX="--maxevents $maxevents"; else MAX=""; fi

# Use * or other glob wildcards for filenames

# Generate steering YAML for QCD
python configs/dqcd/include/ymlgen.py --process 'QCD'        --filerange '[0-10]'

# Higgs
python configs/dqcd/include/ymlgen.py --process 'higgs'      --filerange '*'

python analysis/dqcd.py --runmode genesis  $MAX --inputmap mc_map__higgs_all.yml --config $CONFIG --datapath $DATAPATH
python analysis/dqcd.py --runmode train    $MAX --inputmap mc_map__higgs_all.yml --modeltag higgs_all --config $CONFIG --datapath $DATAPATH --use_conditional $CONDITIONAL
python analysis/dqcd.py --runmode eval     $MAX --inputmap mc_map__higgs_all.yml --modeltag higgs_all --config $CONFIG --datapath $DATAPATH --use_conditional $CONDITIONAL
python analysis/dqcd.py --runmode optimize $MAX --inputmap mc_map__higgs_all.yml --modeltag higgs_all --config $CONFIG --datapath $DATAPATH --use_conditional $CONDITIONAL
