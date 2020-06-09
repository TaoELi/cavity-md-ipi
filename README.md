# Cavity Molecular Dynamics Simulations

## What is cavity MD?

Cavity MD is an approach to simulate coupled photon-nuclear dynamics for realistic molecules in optical microcavities. It aims to accurately describe vibrational strong (VSC) or ultrastrong (V-USC) coupling, i.e., when a cavity mode is resonantly coupled to a vibrational normal mode of molecules and forms a collective Rabi splitting in the infrared (IR) spectrum. For more details, please check [a recent publication](https://arxiv.org/abs/2004.04888).

This approach is still under developing and more features will be reported and updated in the near future.

## What can cavity MD do NOW?

- Computation of Rabi splitting in the IR spectrum for realistic molecular systems from VSC to V-USC without fitting experimental parameters.

- Expolaring other possible modifications on molecular properties under VSC or V-USC.

## Implementation

The implementation of cavity MD relies on conventional MD packages. Currently, I develop cavity MD on top of the [i-pi](http://ipi-code.org/) package, which provides a user-friendly interface for nuclear dynamics. Because the [i-pi](http://ipi-code.org/) package is written in python, it enables a fast realization of new algorithms and ideas when the cavity is considered.

In the future, when the development of cavity MD is stablized, I will also implement this approach in other MD packages for a faster and computationally cheaper realization. If you are a MD developer and has questions on the implementation of cavity MD in your package, please feel free to contact me: t.e.li@outlook.com or taoli@sas.upenn.edu.

## How to use cavity MD to do simulation NOW?

Please check the README in folder "water_VUSC" for a step-by-step introduction on [simulating VSC and V-USC for liquid water](https://arxiv.org/abs/2004.04888). If you are experienced in MD, it will be very easy to replace the force field and configuration of liquid water to other molecules and do your own simulations. Again, for any questions, feel free to contact me: t.e.li@outlook.com or taoli@sas.upenn.edu.

I will also keep updating new examples for other applications.

## Citations

- If you find cavity MD useful for your research, please cite:

Li, T. E., Nitzan, A., & Subotnik, J. E. (2020). Cavity molecular dynamics simulations of liquid water under vibrational ultrastrong coupling. [http://arxiv.org/abs/2004.04888](http://arxiv.org/abs/2004.04888)

- If you directly use the provided code in your simulations, please also cite the original [i-pi](http://ipi-code.org/) package which provides a fancy and user-friendly interface for simulating nuclear dynamics:

Kapil, V., Rossi, M., Marsalek, O., Petraglia, R., Litman, Y., Spura, T., … Ceriotti, M. (2019). i-PI 2.0: A universal force engine for advanced molecular simulations. Computer Physics Communications, 236, 214–223. [https://doi.org/10.1016/j.cpc.2018.09.020](https://doi.org/10.1016/j.cpc.2018.09.020)

