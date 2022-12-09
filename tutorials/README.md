## Tutorials

Three tutorials are available for exploring different features of CavMD:

### Polariton spectrum under collective vibrational strong coupling

See [Rabi_splitting/](Rabi_splitting/). This tutorial uses force-field CavMD to simulation the polariton spectrum of liquid water and can be easily done with a personal laptop.

### Chemical reaction for a single molecule under vibrational strong coupling

See [single_molecule_reaction_vsc/](single_molecule_reaction_vsc/). This tutorial runs first-principles CavMD and can be easily done with a personal laptop.

### Large-scale CavMD QM/MM simulations under collective vibrational strong coupling

See [large_scale_qmmm_reaction_vsc/](large_scale_qmmm_reaction_vsc/). This tutorial runs large-scale QM/MM CavMD and can be done on a Linux cluster. **However, because currently the QM/MM feature requires to slightly modify the Q-Chem source code, and the Q-Chem package is not open-source, please send me an email (taoli@sas.upenn.edu) if you need this tutorial!**

**We are currently working on a more efficient implementation of QM/MM CavMD. In the near future, one may run the simulations with the commercial version of Q-Chem (or other electronic structure packages) + our modified i-pi code without touching any source code. Stay tuned!**

**Before checking the tutorials, please install i-pi as follows.**

## Installation of CavMD

Install the modified i-pi (in folder ../i-pi-master-py3/) package on your personal Linux computer or goverment supercomputer (e.g., NERSC); see [here](http://ipi-code.org/resources/documentation/) for a guide. Here, the simplest installation method is provided.

(Before installation) I strongly suggest the readers to install Anaconda or Miniconda (https://docs.conda.io/en/latest/miniconda.html) before installing i-pi for better experience. Python 3 (up to 3.9) is needed for simulation.

  Please open a terminal
```bash
cd which_path_you_download/cavity-md-ipi/i-pi-master-py3/
python3 setup.py build
python3 setup.py install
```
 copy the following command to the end of your ~/.bashrc file
```bash
source which_path_you_download/cavity-md-ipi/i-pi-master-py3/env.sh
```
in the terminal, run
```bash
 source ~/.bashrc
 which i-pi
 ```
 If the output looks similar as
 ```bash
which_path_you_download/cavity-md-ipi/i-pi-master-py3/bin/i-pi
 ```
 it means we have correctly installed the i-pi package.
