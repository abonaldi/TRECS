######################################################
# T-RECS: The Tiered Radio Extragalactic Continuum   #
# Top Level Makefile                                 #
#                                                    #
# DO NOT modify this Makefile                        #
# Note: Modify the make.inc file to suit your system #
######################################################

TOPSRCDIR := $(PWD)
include $(TOPSRCDIR)/make.inc

F90FLAGS += $(OMPFLAG)

L_GSL     = -L$(GSL_LIB) -lsgl
I_CFITSIO = -I$(CFITSIO_INC)
L_CFITSIO = -L$(CFITSIO_LIB) -lcfitsio
I_HEALPIX = -I$(HEALPIX_INC)
L_HEALPIX = -L$(HEALPIX_LIB) -lhealpix -lgif
L_LAPACK  = -L$(LAPACK_DIR) -llapack

INCLUDE := $(I_CFITSIO) $(I_HEALPIX) 
LINK := $(LINK) $(L_HEALPIX) $(L_CFITSIO) $(L_LAPACK) $(LDFLAGS) 

###############################################################################################
# Paths within the top directory

SRCDIR = $(TOPSRCDIR)/src
MODDIR = $(PREFIX)/mod
ifdef F90RMOD
F90FLAGS += $(F90RMOD) $(MODDIR)
endif
PYTHONDIR = $(TOPSRCDIR)/python

###############################################################################################
# Define objects

OBJ_MOD = $(BUILDDIR)/random_modules.o $(BUILDDIR)/sampler_io.o $(BUILDDIR)/sampler_modules.o
OBJ_CON = $(BUILDDIR)/sampler_continuum.o
OBJ_HyI = $(BUILDDIR)/sampler_hi.o
OBJ_WRP = $(BUILDDIR)/wrapper.o

###############################################################################################
# Rules to build library

default: all

all: mkdirs modules continuum hi xmatch_hi xmatch_clustering wrapper ending

mkdirs:
	mkdir -p $(BUILDDIR) $(MODDIR) $(PREFIX)/bin 

modules: $(OBJ_MOD)

continuum: $(OBJ_MOD) $(OBJ_CON)
	$(F90) $(F90FLAGS) -o $(PREFIX)/bin/trecs_sampler_continuum $(OBJ_MOD) $(OBJ_CON) $(LINK)

hi: $(OBJ_MOD) $(OBJ_HyI)
	$(F90) $(F90FLAGS) -o $(PREFIX)/bin/trecs_sampler_hi $(OBJ_MOD) $(OBJ_HyI) $(LINK)

xmatch_hi: $(PYTHONDIR)/xmatch_hi.py
	sed "s|@PYTHONSHEBANG@|! $(shell which python)|" $(PYTHONDIR)/xmatch_hi.py \
	> $(PREFIX)/bin/trecs_xmatch_hi && \
	chmod u+x $(PREFIX)/bin/trecs_xmatch_hi

xmatch_clustering: $(PYTHONDIR)/xmatch_clustering.py
	sed "s|@PYTHONSHEBANG@|! $(shell which python)|" $(PYTHONDIR)/xmatch_clustering.py \
	> $(PREFIX)/bin/trecs_xmatch_clustering && \
	chmod u+x $(PREFIX)/bin/trecs_xmatch_clustering

wrapper: $(OBJ_MOD) $(OBJ_WRP)
	$(F90) $(F90FLAGS) -o $(PREFIX)/bin/trecs_wrapper $(OBJ_MOD) $(OBJ_WRP) $(LINK)

###############################################################################################
# Generic rules

$(BUILDDIR)/%.o : $(SRCDIR)/%.f90
	$(F90) $(F90FLAGS) $(INCLUDE) -o $@ -c $<

###############################################################################################
# Phony rules

.PHONY: clean ending

ending:
	cp $(SRCDIR)/trecs.sh $(PREFIX)/bin/trecs && chmod u+x $(PREFIX)/bin/trecs
	$(info )
	$(info Success!!)
	$(info )
	$(info Library built in $(PREFIX)/bin)
	$(info To use it remember to add it to your search path)
	$(info now removing build directory)
	rm -rf $(BUILDDIR)

clean:
	$(info cleaning up)
	rm -f $(BUILDDIR)/*.o $(MODDIR)/*.mod *~ *# $(PREFIX)/bin/trecs*

###############################################################################################
