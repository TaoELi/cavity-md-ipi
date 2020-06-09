#!/bin/bash

# This script submit all necessary jobs accordingly

RUN=run_NVE

for E0 in 1e-4 3e-4 5e-4 0e-4 2e-4 4e-4 6e-4 7e-4 8e-4
do
    DIR=E0_$E0
    mkdir $DIR
    cp data.lmp in.lmp init.pdb input_equlibrate.xml input_traj.xml.bak photon_params.json run_diff.sh $DIR
    cd $DIR
    echo "move in $DIR"
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.E0_$E0/" input_equlibrate.xml
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.E0_$E0/" input_traj.xml.bak
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.E0_$E0/" in.lmp
    sed -i "s/0.0005/$E0/" photon_params.json
    mv run_diff.sh NVE_$E0.sh
    sbatch NVE_$E0.sh
    cd ..
done
