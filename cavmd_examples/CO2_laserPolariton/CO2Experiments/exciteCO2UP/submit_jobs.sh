#!/bin/bash

# This script submit all necessary jobs accordingly

RUN=ExciteCO2UP

for Amp in 6e-4 6e-3 #8e-4 1e-3 2e-3 3e-3 4e-3 5e-3
do
    DIR=Amp_$Amp
    mkdir $DIR
    cp data.lmp in.lmp photon_params.json run_diff.sh input_traj.xml.bak  $DIR
    cd $DIR
    echo "move in $DIR"
    sed -i "s/mesitylene-pimd.1/co2-nemd.$RUN.Amp_$Amp/" input_traj.xml.bak
    sed -i "s/mesitylene-pimd.1/co2-nemd.$RUN.Amp_$Amp/" in.lmp
    sed -i "s/0.05/$Amp/" photon_params.json
    
    mv run_diff.sh NECO2_$Amp.sh
    sbatch NECO2_$Amp.sh
    cd ..
done
