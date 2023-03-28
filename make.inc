####################################################
# T-RECS: The Tiered Radio Extragalactic Continuum #
# Make include file                                #
#                                                  #
# Modify this file to suit your system             #
# Note: DO NOT modify the general Makefile         #
####################################################

SHELL = /bin/sh

# PREFIX sets the root directory for installation of the binaries
#
PREFIX = /opt/trecs/2.0

#  CC is the C compiler, normally invoked with options CFLAGS.
#
CC = gcc
CFLAGS = -O3 -std=C99

# Modify the F90 and F90FLAGS definitions to the desired compiler
# and desired compiler options for your machine.
#
F90 = gfortran
F90FLAGS = -O3

# Provide the flag in your compiler used to redirect
# .mod builds in directory different than current
# (e.g. gfortran -> '-J', ifort -> '-module')
# If commented or left empty modules will be compiled
# in the current directory
#
F90RMOD = -J

# OpenMP flag for your system
#
OMPFLAG = -fopenmp

# Eventual additional link flags
#
LDFLAGS = -lgcc

###############################################################################################
# Dependencies

# GSL library
#
GSL_DIR =
GSL_LIB ?= $(GSL_LIB)
L_GSL = -L$(GSL_LIB) -lgsl -lgslcblas

# CFITSIO library
#
CFITSIO_DIR = /opt/cfitsio/4.2.0
CFITSIO_INC = $(CFITSIO_DIR)/include
CFITSIO_LIB = $(CFITSIO_DIR)/lib

# HEALPIX library
#
HEALPIX_DIR = /opt/healpix/3.82
HEALPIX_INC = $(HEALPIX_DIR)/include
HEALPIX_LIB = $(HEALPIX_DIR)/lib

# LAPACK library
#
LAPACK_DIR = /opt/lapack/3.11.0
LAPACK_LIB = $(LAPACK_DIR)

###############################################################################################

BUILDDIR = /tmp/TRECS_build

###############################################################################################


