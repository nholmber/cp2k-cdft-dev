# CP2K CDFT development branch 
This repository contains the constrained density functional theory (CDFT) development branch for [**CP2K**](https://www.cp2k.org/ "CP2K Project"). CP2K is a freely available quantum chemistry package to perform atomistic simulations of solid state, liquid, molecular, periodic, material, crystal, and biological systems. 

The theoretical basis, features, and usage of the CDFT/CP2K implementation will in the future be documented at the [Wiki](https://github.com/nholmber/cp2k-cdft-dev/wiki). In the meanwhile, please refer to the publication: 

Holmberg, Nico; Laasonen, Kari, *Efficient Constrained Density Functional Theory Implementation for Simulation of Condensed Phase Electron Transfer Reactions*, J. Chem. Theory Comput., Just Accepted Manuscript (2016), doi: [10.1021/acs.jctc.6b01085](https://dx.doi.org/10.1021/acs.jctc.6b01085 "Online Version of Publication").


## Available branches

1. `main`: latest stable build (ideally a clone of the official cp2k [`trunk`](https://github.com/cp2k/cp2k) branch, with all CDFT features pushed upstream via svn) 
	* currently based on cp2k version `svn:15309` with some backport buxfixes from newer versions, and **NOT** yet available from the [official SVN mirror]("https://sourceforge.net/projects/cp2k/?source=navbar")  
	* code cleanup (debug prints, etc.) and manual documentation still in complete 

2. unstable development branches (might not compile, incomplete implementation of new features, code not cleaned) 
	* `spin-constraint`: will add possibility to define combined charge + spin constraints
	* `trunk-port` (TODO): implement necessary modifications to merge `main` branch to upstream `trunk` 

## Installation

1. Clone the latest stable build in a directory of your choice `git clone https://github.com/nholmber/cp2k-cdft-dev`
2. Follow the installation instructions in the file [`INSTALL`](INSTALL), which is also available at the official [CP2K project page]("https://www.cp2k.org/howto:compile") (note the difference in versions) 

## Running a CDFT calculation

An official tutorial is in preparation and will be added to the [Wiki](https://github.com/nholmber/cp2k-cdft-dev/wiki). A minimal how-to guide is included in the folder [`example-inputs/`](example_inputs/). A `HTML` keyword manual is also available in [`manual/`](manual/)

## Bugs, ideas and contributing

The implemented CDFT features in the `main` branch have been tested, but some unforseen bugs might still be present. If you encounter problems please open a new issue in the [issue tracker](https://github.com/nholmber/cp2k-cdft-dev/issues) and include a minimal input file for reproducing the bug. New development ideas and contributions are also welcomed.

