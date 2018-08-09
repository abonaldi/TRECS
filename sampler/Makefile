

F90 = gfortran
CC = cc
F90FLAGS = -O3 -openmp
I_CFITSIO =  -I/home/a.bonaldi/software/cfitsio/include/
L_CFITSIO =  -L/home/a.bonaldi/software/cfitsio/lib -lcfitsio
I_HEALPIX = -I/home/a.bonaldi/software/Healpix_3.31/include
L_HEALPIX =  -L/home/a.bonaldi/software/Healpix_3.31/lib -lhealpix -lgif
L_GSL= -L/data/tbourke/bogan/tlb/progs/star-namaka/lib/ -lgsl -lgslcblas
L_LAPACK=-L/home/bonaldi/libraries/ -llapack 
LDFLAGS = 

INCLUDE := $(I_CFITSIO) $(I_HEALPIX) 
LINK := $(LINK) $(L_HEALPIX) $(L_CFITSIO)  $(LDFLAGS) 


OBJ := sampler_io.o random_modules.o sampler_modules.o sampler.o ../../liblapack.so
default: all
all : $(OBJ)
	$(F90) $(F90FLAGS) -o ../bin/sampler $(OBJ) $(LINK) 

%.o : %.f90
	$(F90) $(F90FLAGS) $(INCLUDE) -o $@ -c $<

clean: rm *.o
