# Running  CavMD/ReaxFF using i-pi+LAMMMPS

In CavMD simulations, the cavity photons interact with the dipole derivatives (classically speaking, the partial charges of nuclei) at every time step. In (nonreactive) classical MD, the partial charges of nuclei are always constants. When chemical reactions are simulated, however, bonding breaking/formation implies  changeable dipole derivatives (partial charges). Hence, when using CavMD simulations to study chemical reactions inside a cavity, we need to update dipole derivatives (partial charges) at every time step.

This additional requirement differentiates between (reactive) CavMD with conventional MD. When chemical reactions are implemented on the level of reactive force fields (ReaxFF), the .cpp and .h files in this folder are needed for the communication between i-pi and LAMMPS/ReaxFF.

## Practical installation procedure

1. Please copy these two files (.cpp and .h) to the LAMMPS source code, e.g., lammps-stable_3Mar2020 [https://lammps.sandia.gov/tars/lammps-3Mar20.tar.gz].:
<pre><code>cp fix_cavph.cpp fix_cavph.h lammps-patch_15Jun2020/src/USER-MISC/
 </code></pre>

2. In the LAMMPS folder, recompile the LAMMPS source code and generate the new **lmp** file
 <pre><code>make clean
 cmake -C ../cmake/presets/most.cmake -C ../cmake/presets/nolib.cmake -D PKG_GPU=off ../cmake
 make -j 18 </code></pre>

3. Copy this binary file **lmp**  to the folder where your $PATH environment recognizes, e.g.,
 <pre><code> cp lmp ~/local/bin/ </code></pre>
 where we have assumed that $PATH contains ~/local/bin/.

4. Finally, ensure **lmp** is callable:
 <pre><code> which lmp </code></pre>
 If everything works fine, you will find the output looks like "~/loca/bin/lmp".

## Practical usage procedure

The above installation procedure creates a new fix in LAMMPS: **cavphipi**. This fix resembles the original **ipi** fix in LAMMPS but will also send the partial charges of nuclei to the i-pi code at every time step, which is necessary for a ReaxFF simulation.

In order to run a CavMD/ReaxFF simulation on the i-pi/LAMMPS framework, we need to do the following three steps:

1. In lmp.in (LAMMPS input file), please change the **fix ipi** to **fix cavphipi** to allow the functionality of the new communicator.

2. In photon_params.json (the CavMD input file), please add a new line **"update_charge": true** to allow i-pi read the partial charges from **cavphipi**.

3. Please also prepare the other LAMMPS/ReaxFF files that are necessary for a simulation outside the cavity.

These steps will allow to run CavMD/ReaxFF on the ipi+LAMMPS infrastructure. Please also check the tutorial for more details.  
