#!/bin/bash -l

#SBATCH -q shared
#SBATCH --mem=8GB
#SBATCH -t 48:00:00
#SBATCH -L SCRATCH,project
#SBATCH -C haswell
#SBATCH -A m3138

module load python3
module load lammps/2018.12.12-hsw
export PYTHONPATH=$PYTHONPATH:~/lib/python3.7/site-packages/
export PATH=~/bin:$PATH


FILE_EQU=RESTART_traj_0
CHECKPOINT=equlibrate.checkpoint
if [ ! -f "$FILE_EQU" ]; then
    # Run for equlibrium
    # first check if the equlibrating job is not finished
    if [ -f "$CHECKPOINT" ]; then
	echo "Continuing equlibrating job 150 ps"
	i-pi $CHECKPOINT &> log &
	sleep 60s
	srun -n 1 -c 1 --cpu-bind=cores lmp_cori < in.lmp
	sleep 60s
    else
	echo "Performing equlibrating job 150 ps"
	i-pi input_equlibrate.xml &> log &
	sleep 60s
	srun -n 1 -c 1 --cpu-bind=cores lmp_cori < in.lmp
	sleep 60s
    fi
    mv $CHECKPOINT $FILE_EQU || exit 1
else
    echo "Skip equlibrating job"
fi

cp $FILE_EQU simu_0.checkpoint || exit 1

PATTERN="<step>40000</step>"
# Run for different slices for properties calculation
if [ ! -f "$FILE_EQU" ]; then
    echo "NO INPUT FILE FOR TRAJ SIMULATION"
else
    for i in {1..40}
    do
	echo "Dealing with $i slice"
	CHECKPOINT=simu_$i.checkpoint
	if grep -q $PATTERN $CHECKPOINT; then
	    echo "found checkpoint finished"
	    echo "Skip $i-th sequential job"
	else
	    echo "not found checkpoint finished"
	    a=$(wc -c < $CHECKPOINT)
	    if [ ! -f "$CHECKPOINT" ] || [ $a -le 1 ]; then
		echo "Performing $i-th simulation 20 ps"
		cp input_traj.xml.bak input_traj_$i.xml
		sed -i "s/'simu'/'simu_$i'/g" input_traj_$i.xml
		sed -i "s/RESTART/simu_$(($i-1)).checkpoint/g" input_traj_$i.xml
		i-pi input_traj_$i.xml &> log &
		sleep 60s
		srun -n 1 -c 1 --cpu-bind=cores lmp_cori < in.lmp
		sleep 60s
	    else
		echo "Continuing $i-th simulation 20 ps"
		i-pi $CHECKPOINT &> log &
		sleep 60s
		srun -n 1 -c 1 --cpu-bind=cores lmp_cori < in.lmp
		sleep 60s
	    fi
	    cp $CHECKPOINT RESTART_traj_$i || exit 1
	fi
    done
fi
