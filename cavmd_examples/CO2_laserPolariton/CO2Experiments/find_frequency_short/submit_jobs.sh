#!/bin/bash

# This script submit all necessary jobs accordingly

RUN=find_freq

for Amp in 9e-2 9.5e-1 1e-1 1.1e-1 1.2e-1 2.5e-2 3.5e-2 4.5e-2 5.5e-2 6.5e-2 7.5e-2 8.5e-2 7e-2 8e-2 9e-2 #1e-2 2e-3 4e-3 5e-3 7e-3 8e-3 9e-3 2e-2 4e-2 5e-2 3e-4 6e-4 1e-3 3e-3 6e-3  3e-2 6e-2
do
    DIR=Amp_$Amp
    mkdir $DIR
    cp data.lmp in.lmp init.xyz input_equlibrate.xml photon_params.json run_diff.sh $DIR
    cd $DIR
    echo "move in $DIR"
    sed -i "s/mesitylene-pimd.1/taoli.$RUN.Amp_$Amp/" input_equlibrate.xml
    sed -i "s/mesitylene-pimd.1/taoli.$RUN.Amp_$Amp/" in.lmp
    sed -i "s/0.05/$Amp/" photon_params.json
    mv run_diff.sh Amp_$Amp.sh
    sbatch Amp_$Amp.sh
    cd ..
done
