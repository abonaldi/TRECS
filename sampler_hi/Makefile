

F90 = gfortran
CC = cc
F90FLAGS = -O3 -openmp
I_CFITSIO =  -I/home/a.bonaldi/data-cold-for-backup/software/cfitsio/include/
L_CFITSIO =  -L/home/a.bonaldi/data-cold-for-backup/software/cfitsio/lib -lcfitsio
I_HEALPIX = -I/home/a.bonaldi/data-cold-for-backup/software/Healpix_3.31/include
L_HEALPIX =  -L/home/a.bonaldi/data-cold-for-backup/software/Healpix_3.31/lib -lhealpix -lgif
L_GSL= -L//home/a.bonaldi/data-cold-for-backup/software/libtyler/ -lgsl -lgslcblas
L_LAPACK=-L/home/bonaldi/libraries/ -llapack 
LDFLAGS = 

INCLUDE := $(I_CFITSIO) $(I_HEALPIX) 
LINK := $(LINK) $(L_HEALPIX) $(L_CFITSIO)  $(LDFLAGS) 

OBJ := ../modules/random_modules.o ../modules/sampler_io.o  ../modules/sampler_modules.o sampler.o ../../liblapack.so

default: all
all : $(OBJ)
	$(F90) $(F90FLAGS) -o ../bin/sampler_hi $(OBJ) $(LINK) 

%.o : %.f90
	$(F90) $(F90FLAGS) $(INCLUDE) -o $@ -c $<

clean: 
	rm *.o ../modules/*.o *.mod
