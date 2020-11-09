#!/bin/bash -l

#SBATCH -q shared
#SBATCH --mem=8GB
#SBATCH -t 48:00:00
#SBATCH -L SCRATCH,project
#SBATCH -C haswell
#SBATCH -A m3138

module load python3
module load lammps
export PYTHONPATH=$PYTHONPATH:~/lib/python3.7/site-packages/
export PATH=~/bin:$PATH

cp photon_params.json photon_params.json.bak

PATTERN="<step>200000</step>"
PATTERN_INIT="<step>40000</step>"
# This simulation starts from each equilibrium traj and rerun it by sending an external
# pluse to the system. 
for i in {1..40}
do
    echo "Dealing with $i trajectory"
    CHECKPOINT=simu_$i.checkpoint
    if grep -q $PATTERN $CHECKPOINT; then
        echo "found traj finished"
        echo "Skip $i-th sequential job"
    else
        echo "not found checkpoint finished"
        if [ ! -f "$CHECKPOINT" ]; then
	   echo "even not found checkpoint existed!"
	   echo "Performing $i-th simulation 100 ps..."
	   cp input_traj.xml.bak input_traj_$i.xml
	   sed -i "s/'simu'/'simu_$i'/g" input_traj_$i.xml
	   cp ../../CO2only_changeFreq/Freq_2e-4/simu_$(($i-1)).checkpoint init_$(($i-1)).checkpoint || exit 1 
	   sed -i "s/RESTART/init_$(($i-1)).checkpoint/g" input_traj_$i.xml
	   # if this eq (init) checkpoint file has not finished; exit simulation
	   if [ $i -ge 2 ]; then
	       if grep -q $PATTERN_INIT init_$(($i-1)).checkpoint; then
		   echo "this eq (init) checkfile is already finished"
	       else
		   rm init_$(($i-1)).checkpoint
		   echo "this eq (init) checkfile has not finished, exiting..."
		   exit 1
	       fi
	   fi
	   cp photon_params.json.bak photon_params.json
	   # change to a random phase for the incoming pulse
	   phase=$(seq 0 .0001 6.2832 | shuf | head -n1)
	   echo "phase is $phase"
	   sed -i "s/3.14/$phase/g" photon_params.json
	   
	   i-pi input_traj_$i.xml &> log_$i &
	   sleep 60s
	   srun -n 1 -c 1 --cpu-bind=cores lmp_cori < in.lmp
	   sleep 60s
        else
           echo "found checkpoint but not finished, meaning that last time the simulation stops at No.$i traj"
	   # under this situation, in order to avoid bad dynamics, we delete this simulation and redo it for safety.
	   echo "Redo $i-th simulation 100 ps" 
	   rm simu_$i.xc.xyz* simu_$i.out 
	   i-pi input_traj_$i.xml &> log_$i &
	   sleep 60s
	   srun -n 1 -c 1 --cpu-bind=cores lmp_cori < in.lmp
	   sleep 60s
        fi
    fi
done
