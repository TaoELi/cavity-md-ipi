#!/bin/bash

# This script submit all necessary jobs accordingly

RUN=run_CO2_freq

E0=2e-4

for freq in 2025 2075 2125 2175 2225 2275 2325 2375 2425 2475 2450 2525 2550 2000 2050 2100 2150 2200 2250 2300 2350  2400 2500
do
    DIR=Freq_$freq
    mkdir $DIR
    cp data.lmp in.lmp init.xyz input_equlibrate.xml input_traj.xml.bak photon_params.json run_diff.sh $DIR
    cd $DIR
    echo "move in $DIR"
    sed -i "s/mesitylene-pimd.1/taoli.$RUN.freq_$freq/" input_equlibrate.xml
    sed -i "s/mesitylene-pimd.1/taoli.$RUN.freq_$freq/" input_traj.xml.bak
    sed -i "s/mesitylene-pimd.1/taoli.$RUN.freq_$freq/" in.lmp
    sed -i "s/0.0005/$E0/" photon_params.json
    sed -i "s/2380/$freq/" photon_params.json
    mv run_diff.sh freq_$freq.sh
    sbatch freq_$freq.sh
    cd ..
done
