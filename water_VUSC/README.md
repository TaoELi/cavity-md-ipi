# HOW TO run cavity MD for liquid water?

1. "Install" the modified i-pi on your personal computer or goverment supercomputer (e.g., nersc); see [here](http://ipi-code.org/resources/documentation/) for a guide. One can do it by 
<pre><code>
cd which_path_you_download/cavity-md-ipi/i-pi/
python2 setup.py build
python2 setup.py install --prefix=~
</code></pre> 
which will allow the installation of i-pi at $HOME/bin/i-pi.

2. Ensure that LAMMPS is also installed and you can use it with command <pre><code>lmp < in.lmp </code></pre>

3. Make sure that you can run simple test jobs with i-pi by calling LAMMPS.

3.1. If using NERSC, go to each folder (e.g., single_mode_g0/) and run <pre><code>./submit_xx.sh </code></pre> and the jobs will be automatically submitted in NERSC. The whole simulation may take longer than 48 hours. If you find that your job is killed by server while the job has not completely finished, please resubmit with command ./submit_xx.sh and the job can continue from the last checkpoint.

3.2. If using personal computer, please go to folder test/ and run <pre><code> i-pi input_traj_1.xml &> log & </code></pre> then <pre><code>lmp < in.lmp </code></pre> After a while (~ 1h), you will finish one-trajectory simulation.

4. After simulation, we should analysis 80 * 20 ps equilibrium trajectories. Command line <pre><code>python collect_*.py folder/subfolders </code></pre> to obtain all necessary molecular properties.

5. Go to each folder, run <pre><code>python plot_which_you_are_interested_in.py </code></pre> and you will obtain the published figures.

# What is different from conventional i-pi simulation?

The difference is small. Assuming that you have prepared all files for runing conventional nuclear dynamics with i-pi + LAMMPS, then please add photon modes at the end of the init.pdb file like this:

<pre><code>
ATOM    649    L   1     1       0.000   0.000   0.000  0.00  0.00            0
ATOM    650    L   1     1       0.000   0.000   0.000  0.00  0.00            0
</code></pre>

Here, I add two photons (labelled by L). Typically the number of photons should be even due to the two polarization  directions (which are set as x and y directions here) of each single cavity mode, where the cavity direction is set at z-drection. For the LAMMPS data file (data.lmp), please leave it alone and do not do modifications, because we will not throw the information of photons to LAMMPS. 

Second, please prepare a photon_params.json file to control the parameters of cavity mode:
<pre><code>
{
  "apply_photon" : true,
  "eff_mass" : 1.0,
  "freqs_cm" : 3550.0,
  "E0" : 0.0005
}
</code></pre>
These information are mandatory for cavity MD simulations:
- "apply_photon" denotes whether including cavity effects or not.
- "eff_mass" denotes the effective mass of photons, which is taken as 1.0 a.u. (atomic units) for convenience.
- "freqs_cm" denotes the frequency of the fundamental photon mode in units of wave number.
- "E0" denotes <img src="https://latex.codecogs.com/svg.latex?\tilde{\varepsilon}" /> (effective coupling strength in units of a.u.) for the fundamental photon mode; see [here](https://arxiv.org/abs/2004.04888) for details.

# Multiple Rabi splitting

Besides, there are some optional parameters which allow more complicated cavity MD simulations. For example, if the photon_params.json is 
<pre><code>
{
  "apply_photon" : true,
  "eff_mass" : 1.0,
  "freqs_cm" : 1775,
  "E0" : 0.0002,
  "n_modes" : 4
}
</code></pre>
- The new option "n_modes" will include 4 (four different cavity modes with spacing "freqs_cm") * 2 (two polarization directions) photons, which will allow simulating multiple Rabi splitting, i.e., different cavity modes forms Rabi splittings with different vibrational normal modes of molecules.  Of course, in init.pdb, please add 4 * 2 photons at the end. All photon modes are still at x or y-direction. Note that here I do not include the photons with in-plane wave vectors rather than 0.

Remember that "E0" denotes the effective coupling strength for the fundamental (or with smallest frequency) cavity mode.

# Going beyond water simulation

If you want to go beyond the cavity MD simulations of liquid water and try other molecules, please change init.pdb, in.lmp, and data.lmp to the system you are interested in and also add an additional control to photon_params.json:
<pre><code>
{
  "apply_photon" : true,
  "eff_mass" : 1.0,
  "freqs_cm" : 1775,
  "E0" : 0.0002,
  "n_modes" : 4,
  "charge_array": [0.33, -0.33, ...]
}
</code></pre>
Here, the option "charge_array" will redefine the partial charge of each atom in the order of configurations. The number of partial charges should match the number of atoms in the simulations.
