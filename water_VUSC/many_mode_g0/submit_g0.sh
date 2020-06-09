#!/bin/bash

# This script is used to submit jobs for water simulation but supports N modes field :)

omega0=3550.0

RUN=Mode_G0_20ps
N=4
for E0 in 0e-4 1e-4 2e-4 3e-4 4e-4 5e-4 6e-4
do
    dir=E0_$E0
    mkdir $dir
    cp data.lmp in.lmp photon_params.json init.pdb input_equlibrate.xml input_traj.xml.bak run_diff.sh $dir
    cd $dir
    # change photon_params to E0 = $E0
    sed -i "s/\"n_modes\" : 1/\"n_modes\" : $N/" photon_params.json
    
    if [ $N -eq 1 ]; then
	# If N = 1; we run the previous simulation
	echo "now perform N = $N"
    else
	echo "now perform N = $N not 1"
	# change fundamental freq to omeag * 2.0 / N
	omega=$(echo - | awk "{print 3550*2.0/$N}")	
	sed -i "s/3550.0/$omega/" photon_params.json
	# change coupling strength to 4e-4 * 2.0 / N
	E0_now=$(echo - | awk "{print $E0*2.0/$N}")	
	sed -i "s/0.0002/$E0_now/" photon_params.json

	echo "fundamental freq $omega fundamental strength $E0_now"
	
	cat photon_params.json
	# For N > 1, I need to change init.pdb configuration to include more photons
	for i in $( seq 3 1 $(($N*2)) )
	do
	    num=$((651+$N*2-$i))
	    added="ATOM    $num    L   1     1       0.000   0.000   0.000  0.00  0.00            0"
	    sed -i "653i $added" init.pdb
	done
    fi
    
    # apart from these, the job submission flowchart is the same as previous, namely change ID
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.E0_$E0/" input_equlibrate.xml
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.E0_$E0/" input_traj.xml.bak
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.E0_$E0/" in.lmp
    mv run_diff.sh E020ps_$E0.sh
    sbatch E020ps_$E0.sh

    cd ..
done
