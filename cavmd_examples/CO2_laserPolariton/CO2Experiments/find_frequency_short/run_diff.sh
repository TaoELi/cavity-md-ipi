#!/bin/bash -l

#SBATCH -q shared
#SBATCH --mem=8GB
#SBATCH -t 2:00:00
#SBATCH -L SCRATCH,project
#SBATCH -C haswell
#SBATCH -A m3138

module load python3
module load lammps/2018.12.12-hsw
export PYTHONPATH=$PYTHONPATH:~/lib/python3.7/site-packages/
export PATH=~/bin:$PATH


FILE_EQU=RESTART_traj_0
CHECKPOINT=simu.checkpoint
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

