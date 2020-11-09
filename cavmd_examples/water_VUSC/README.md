# CavMD for VSC and V-USC of liquid water

The folders here contain all necessary files to generate figures in the following paper:

- Li, T. E., Subotnik, J. E., & Nitzan, A. (2020). Cavity molecular dynamics simulations of liquid water under vibrational ultrastrong coupling. Proceedings of the National Academy of Sciences, 117(31), 18324–18331. https://doi.org/10.1073/pnas.2009272117

## NERSC users

The simulation was run on the government computer NERSC.

If using NERSC, go to each folder (e.g., single_mode_g0/, where most published results can be obtained here) and run <pre><code>./submit_xx.sh </code></pre> and the jobs will be automatically submitted in NERSC. The whole simulation may take longer than 48 hours. If you find that your job is killed by server while the job has not completely finished, please resubmit with command ./submit_xx.sh and the job can continue from the last checkpoint. Actually, rough results (which are close to publication results by have some noise) can be obtained by running a few trajectories. If users just want to recover the results in https://doi.org/10.1073/pnas.2009272117, running a few trajectories are OK.

After simulation, we should analysis 80 * 20 ps (or just a few 20 ps) equilibrium trajectories. Command line <pre><code>python collect_all_data_N.py folder/subfolders </code></pre> to obtain all necessary molecular properties. Go to each folder, run
<pre><code>module load texlive
python plot_which_you_are_interested_in.py </code></pre> and you will obtain the published figures.

## Non-NERSC users

For users who have no access to NERSC but have access to other linux clusters, the following commands may need to be modified. Taking folder single_mode_g0/ for example,
in **submit_jobs.sh**, the command
<pre><code>sbatch NVE_$E0.sh </code></pre>
may need to be modified if the job in the cluster is not submitted through sbatch but with other programs (e.g., PBS).
In **run_diff.sh**, please change the resource allocation definitions (here a single CPU for 48 hours is requested):
<pre><code>#!/bin/bash -l
#SBATCH -q shared
#SBATCH --mem=8GB
#SBATCH -t 48:00:00
#SBATCH -L SCRATCH,project
#SBATCH -C haswell
#SBATCH -A m3138
</code></pre>
Second, please load enough packages so that the modified i-pi and LAMMPS can be performed directly in the script.
<pre><code>module load python3
#module load python/2.7-anaconda-2019.07
module load lammps
#export PYTHONPATH=$PYTHONPATH:~/lib/python2.7/site-packages/
export PATH=~/bin:$PATH
</code></pre>
Lastly, the command to call LAMMPS need to be modified to that is suitable for different clusters:
<pre><code>srun -n 1 -c 1 --cpu-bind=cores lmp_cori < in.lmp</code></pre>
