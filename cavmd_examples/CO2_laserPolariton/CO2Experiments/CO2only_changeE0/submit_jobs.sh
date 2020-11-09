#!/bin/bash

# This script submit all necessary jobs accordingly

RUN=only_CO2

for E0 in 0e-4 5e-5 1e-4 2e-4 4e-4 5e-4
do
    DIR=E0_$E0
    mkdir $DIR
    cp data.lmp in.lmp init.xyz input_equlibrate.xml input_traj.xml.bak photon_params.json run_diff.sh $DIR
    cd $DIR
    echo "move in $DIR"
    sed -i "s/mesitylene-pimd.1/taoli.$RUN.E0_$E0/" input_equlibrate.xml
    sed -i "s/mesitylene-pimd.1/taoli.$RUN.E0_$E0/" input_traj.xml.bak
    sed -i "s/mesitylene-pimd.1/taoli.$RUN.E0_$E0/" in.lmp
    sed -i "s/0.0005/$E0/" photon_params.json
    mv run_diff.sh Oco2_$E0.sh
    sbatch Oco2_$E0.sh
    cd ..
done
