CC = gcc
F90 = gfortran
F90FLAGS = -O3 -fopenmp -std=gnu
CFITSIO_DIR = /opt/cfitsio/4.2.0
HEALPIX_DIR = /opt/healpix/3.82
LAPACK_DIR = /opt/lapack/3.11.0
I_CFITSIO =  -I$(CFITSIO_DIR)/include/
L_CFITSIO =  -L$(CFITSIO_DIR)/lib -lcfitsio
I_HEALPIX = -I$(HEALPIX_DIR)/include
L_HEALPIX =  -L$(HEALPIX_DIR)/lib -lhealpix -lgif
L_GSL = -L$(GSL_LIB) -lgsl -lgslcblas
L_LAPACK = -L$(LAPACK_DIR) -llapack
LDFLAGS = -lgcc

INCLUDE := $(I_CFITSIO) $(I_HEALPIX) 
LINK := $(LINK) $(L_HEALPIX) $(L_CFITSIO) $(L_LAPACK) $(LDFLAGS) 

OBJ := ../modules/random_modules.o ../modules/sampler_io.o  ../modules/sampler_modules.o sampler.o

default: all
all : $(OBJ)
	$(F90) $(F90FLAGS) -o ../bin/sampler_continuum $(OBJ) $(LINK) 

%.o : %.f90
	$(F90) $(F90FLAGS) $(INCLUDE) -o $@ -c $<

clean: 
	rm -f *.o *.mod $(OBJ) ../modules/*.mod ../bin/sampler_continuum
