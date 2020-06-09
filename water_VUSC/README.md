# HOW TO run cavity MD for liquid water?

1. "Install" the modified i-pi on your personal computer or goverment supercomputer (e.g., nersc); see [here](http://ipi-code.org/resources/documentation/) for a guide. The simplest way is to add <pre><code> source ~/which_path/cavity-md-ipi/i-pi/env.sh </code></pre> to your .bashrc and then source it.

2. Ensure that LAMMPS is also installed and you can use it with command <pre><code>lmp < in.lmp </code></pre>

3. Make sure that you can run simple test jobs with i-pi by calling LAMMPS.

3.1. If using NERSC, go to each folder (e.g., single_mode_g0/) and run <pre><code>./submit_xx.sh </code></pre> and the jobs will be automatically submitted in NERSC. The whole simulation may take longer than 48 hours. If you find that your job is killed by server while the job has not completely finished, please resubmit with command ./submit_xx.sh and the job can continue from the last checkpoint.

3.2. If using personal computer, please go to folder test/ and run <pre><code> i-pi input_traj_1.xml &> log & </code></pre> then <pre><code>lmp < in.lmp </code></pre> After a while (~ 1h), you will finish one-trajectory simulation.

4. After simulation, we should analysis 80 * 20 ps equilibrium trajectories. Command line <pre><code>python collect_*.py folder/subfolders </code></pre> to obtain all necessary molecular properties.

5. Go to each folder, run <pre><code>python plot_which_you_are_interested_in.py </code></pre> and you will obtain the published figures.
