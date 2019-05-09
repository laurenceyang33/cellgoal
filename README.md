# Cellgoal project

Online Resource for 
**Estimating Cellular Goals from High-Dimensional Biological Data**

## Requirements for ADMM solver
1. Fortran compiler.
	- tested on PGI 18.4 and 18.10: required for OpenACC-based parallelization
	- tested on gfortran 8.2.0
1. cmake
1. make
1. Linear solver:
	1. umfpack: provided by suitesparse package in Linux
	1. MA57: please obtain from www.hsl.rl.ac.uk/catalogue/ma57.html

## Requirements to run test suite
1. Python >= 2.7
1. cobrapy
	- to load metabolic models
	- to generate synthetic fluxes
	- to validate the estimated goal reaction
	
## Installation
1. Install MA57 linear solver
   1. it needs metis-4, not metis-5 that is commonly in package managers. Get metis-4.0.3 from [Karypis lab](http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD).

### ADMM solver
1. cd admm
1. mkdir build
1. cd build
1. cmake ..
1. make

### Test suite
1. python setup.py develop
