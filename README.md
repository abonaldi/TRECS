# T-RECS

This is the code used to produce the Tiered Radio Extragalactic Continuum Simulation (T-RECS, Bonaldi et al. 2018). 
It can be run to produce radio sources catalogues with user-defined frequencies, area and depth. 

Once installed (see the INSTALL file in this folder) the code runs as follows:

1) run sampler from the code/bin folder. The calling sequence is "sampler file.inp", where file.inp contains all the configuration parameters (examples are in the code/bin folder). While running, the code produces catalogues per redshift bin and per galaxy sub-population, on a user-specified output folder. 
For demanding simulations, the code can be easily parallelised by running several instances of the code each processing a different redshift interval. The redshift interval is controlled by the the z_min, z_max parameters that can be set in the input file (default is z_min=0, z_max=8, which means no parallelization). Those do not need to be related to the redshift bins intrinsically defined in the code. 

2) create input files for the wrapper code, which will position the field of view on the user-specified central coordinates and collate all the output catalogues in just two, one for all AGNs and one for all SFGs. This input file can be automatically created with the python/wrapper_input.py script 

3) run wrapper from the code/bin folder with the input files generated, which are separate for SFGs and AGNs. This is not computationally demanding and can be run multiple times, using the same catalogue inputs, to project the simulated sky onto different fields. 


