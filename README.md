# Cellgoal project

Online Resource for 
**Estimating Cellular Goals from High-Dimensional Biological Data**

## Requirements for ADMM solver
1. Fortran compiler.
	- tested on PGI 18.4: required for OpenACC-based parallelization
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
