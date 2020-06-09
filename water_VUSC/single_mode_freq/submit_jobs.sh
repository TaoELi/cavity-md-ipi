#!/bin/bash

# This script submit all necessary jobs accordingly

RUN=run_NVEfreq

for freq in 3550 3950 3350 4150 4950 3150 3750 4350 4550 4750 5150 5350
do
    DIR=Freq_$freq
    mkdir $DIR
    cp data.lmp in.lmp init.pdb input_equlibrate.xml input_traj.xml.bak photon_params.json run_diff.sh $DIR
    cd $DIR
    echo "move in $DIR"
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.$freq/" input_equlibrate.xml
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.$freq/" input_traj.xml.bak
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.$freq/" in.lmp
    sed -i "s/3550/$freq/" photon_params.json
    mv run_diff.sh NVE_$freq.sh
    sbatch NVE_$freq.sh
    cd ..
done
