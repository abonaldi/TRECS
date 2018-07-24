# T-RECS

T-RECS is a software package composed of two fortran codes, sampler and wrapper, and some ancillary python scripts. 

To build the package for the first time:

1) download the archive 

which contains input files with all the models used by the sampler code. Once extracted, copy the folder without renaming it in the same folder where the T-RECS folder is. 


2) check that you have a fortran compiler and these libraries are installed: GSL, Lapack, Healpix, cfitsio

3) edit  sampler/Makefile and  wrapper/Makefile with your paths and compiler options

4) run "make" from both the sampler and the wrapper folders. the executables are saved in the bin folder. 


the python scripts require astropy and astropy.io to run. 