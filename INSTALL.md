# T-RECS

T-RECS is a software package composed of three fortran executables
- trecs_sampler_continuum to create radio continuum simulated catalogues
- trecs_sampler_hi to create HI simulated catalogues
- trecs_wrapper to collate different redshift slices of a simulation and
  project the coordinates to a chosen field of view in the sky. 

and some ancillary python scripts for catalogue cross-matching.

All the executables above are conveniently wrapped into a bash script called 'trecs' that
accepts different flags for customised usage.

To build the package for the first time:

1) download the archive 

   https://www.dropbox.com/s/3u4wtk1fxps6fwg/TRECS_Inputs.zip?dl=0

   (size 8.68 GB)
   which contains input files with all the models used by the sampler code.
   Once downloaded, extract it with 

   gunzip TRECS_Inputs.zip

   the location of the extracted folder will need to be added to your parameter file when running.

2) check that you have a fortran compiler and the following libraries and packages installed:
   - fortran: GSL, Lapack, Healpix, cfitsio
   - python: NumPy, AstroPy, SKLearn

3) edit the input makefile `make.inc` (copying it from examples/make.inc within the root repository directory) with your paths and compiler options

4) run `make`. The executables are saved in the bin folder at the PREFIX provided in `make.inc` 

Everything should be installed correctly, an example of parameter file and of the frequency list file
are available in example/parameter_file.ini and example/frequency_list.dat files respectively
(within the root directory).
