# CP2K

CP2K is a quantum chemistry and solid state physics software package that can perform atomistic simulations of solid state, liquid, molecular, periodic, material, crystal, and biological systems. CP2K provides a general framework for different modeling methods such as DFT using the mixed Gaussian and plane waves approaches GPW and GAPW. Supported theory levels include DFTB, LDA, GGA, MP2, RPA, semi-empirical methods (AM1, PM3, PM6, RM1, MNDO, ...), and classical force fields (AMBER, CHARMM, ...). CP2K can do simulations of molecular dynamics, metadynamics, Monte Carlo, Ehrenfest dynamics, vibrational analysis, core level spectroscopy, energy minimization, and transition state optimization using NEB or dimer method.

CP2K is written in Fortran 2003 and can be run efficiently in parallel using a combination of multi-threading, MPI, and CUDA.

## Links

* [CP2K.org](https://www.cp2k.org) for showcases of scientific work, tutorials, exercises, presentation slides, etc.
* [The manual](https://manual.cp2k.org/) with descriptions of all the keywords for the CP2K input file
* [The dashboard](https://dashboard.cp2k.org) to get an overview of the currently tested architectures
* [The Google group](https://groups.google.com/group/cp2k) to get help if you could not find an answer in one of the previous links

## Directory organization

* `INSTALL.md`: How to build and setup CP2K
* `arch`: Collection of definitions for different architectures and compilers
* `src`: The source code
* `tests`: Inputs for regression tests
* `tools`: Mixed collection of useful scripts related to cp2k
* `data`: Simulation parameters e.g. basis sets and pseudopotentials

Additional directories created during build process:

* `lib`: Libraries built during compilation
* `obj`: Objects and other intermediate compilation-time files
* `exe`: Where the executables will be located
