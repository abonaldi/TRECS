# T-RECS

T-RECS is a software package composed of three fortran codes
- sampler_continuum to create radio continuum simulated catalogues
- sampler_hi to create HI simulated catalogues
- wrapper to collate different redshift slices of a simulation and project the coordinates to a chosen field of view in the sky. 

and python scripts
- xmatch_hi.py to associate conterparts between continuum and hi outputs
- xmatch_clustering.py to associate counterparts between continuum, hi or cross-matched catalogues and DM cosmological simulation, to add clustering
- wrapper_input.py to prepare the input file for wrapper
- addkeys_vizier.py to make the catalogue format compliant with vizier


To build the package for the first time:

1) download the archive 

https://www.dropbox.com/s/crkzwho0hqc565g/TRECS_Inputs.tgz?dl=0

(size 8 GB)
which contains input files with all the models used by the sampler code. Once extracted, 

tar -xvf TRECS_Inputs.tgz

move the folder without renaming it in the same folder where the T-RECS code folder is. Your folder should now contain: 

TRECS  TRECS_Inputs

2) check that you have a fortran compiler and these libraries are installed: GSL, Lapack, Healpix, cfitsio

3) edit  The Makefile for sampler_continuum, sampler_hi and  wrapper with your paths and compiler options

4) run "make" in the sampler_continuum, sampler_hi and wrapper folders. The executables are saved in the bin folder. 


the python scripts require astropy and astropy.io to run. 