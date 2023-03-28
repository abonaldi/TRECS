# T-RECS

This is the code used to produce the Tiered Radio Extragalactic Continuum Simulation
(T-RECS, Bonaldi et al. 2018). 
It can be run to produce radio sources catalogues with user-defined frequencies, area and depth. 

Once the code is installed (see the INSTALL.md file in this folder) and the installation
path has been added to your search path (i.e. export PATH=${PATH}:/path/to/prefix/bin)
the code runs by calling :

```bash
$ trecs [OPTIONAL FLAGS] -p/--params /path/to/parameter_file.ini
```

Available options can be printed on screen by calling

```bash
$ trecs -h/--help
```

A dummy parameter file is available at example/parameter_file.ini as well as a dummy
frequency list file at example/frequency_list.

According to the optional flags with which the master script is called different possible
results will be produced.

## Optional arguments

`-c`/`--continuum`	will run the continuum simulation

`-i`/`--HI`        	will run the HI simulation

`-x`/`--xmatch`    	cross-matches the continuum and HI simulations
	     		(this requires the two flags -c AND -i to have been used,
	       	 	not necessarily on the same run, as long as the output
		 	paths in the parameter file are compatible)

`-C`/`--clustering`	will add clustering properties based on the coordinates of
			sub-haloes in a lightcone built from the P-Millenium simullation

`-w`/`--wrap [tag]` 	will wrap raw catalogues in a single fits file (and eventually
	     		apply a rotation to the coordinates towards some required
			central latitude and longitude. This is required as the above
			options produce catalogues per redshift bin and per galaxy sub-population,
			on a user-specified output folder.

`-h`/`--help`		displays a help message with usage instructions

**NOTE 1)** For demanding simulations, the code can be easily parallelised by running
several instances of the master script code each processing a different redshift interval.
The redshift interval is controlled by the the z_min, z_max parameters that can be set in the
input parameter file (default is z_min=0, z_max=8, which means no parallelization).
The aforementioned parameters do not need to be related to the redshift bins intrinsically
defined in the code.
The resulting catalogues can be then wrapped in a second moment by calling the trecs master scrip
with only the wrap flag.

**NOTE 2)** the wrapping part of the code will position the field of view on the user-specified
central coordinates and collate all the output catalogues in just one. 
This is not computationally demanding and can be run multiple times, using the same catalogue inputs,
to project the simulated sky onto different fields. 


