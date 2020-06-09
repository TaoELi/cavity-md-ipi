#!/bin/bash

# This script is used to submit jobs for water simulation but supports N modes field :)

omega0=3550.0

RUN=Mode20
for N in 1 2 4 6 8 10 12 14 16 18 20
do
    dir=N_$N
    mkdir $dir
    cp data.lmp in.lmp photon_params.json init.pdb input_equlibrate.xml input_traj.xml.bak run_diff.sh $dir
    cd $dir
    # change photon_params to n_modes = $N
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
	E0=$(echo - | awk "{print 0.0002*2.0/$N}")	
	sed -i "s/0.0002/$E0/" photon_params.json

	echo "fundamental freq $omega fundamental strength $E0"
	
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
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.N_$N/" input_equlibrate.xml
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.N_$N/" input_traj.xml.bak
    sed -i "s/h2o-pimd.1.E0_5e-4/h2o-pimd.1.$RUN.N_$N/" in.lmp
    mv run_diff.sh Mode20_$N.sh
    sbatch Mode20_$N.sh

    cd ..
done
