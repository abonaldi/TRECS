#!/bin/bash

#############################################################################################
# Argument parsing

docontinuum=false
dohi=false
doxmatch=false
dowrap=false
params=
while [ $# -gt 0 ]
do
    case $1 in
	-c|--continuum)
	    docontinuum=true
	    ;;
	-i|--HI)
	    dohi=true
	    ;;
	-x|--xmatch)
	    doxmatch=true
	    ;;
	-w|--wrap)
	    dowrap=true
	    ;;
	-p|--params)
	    params=$2
	    if [ ! -f "$params" ] ; then
		echo "Error: provided parameter file does not exist or is not a file"
		exit 1
	    fi
	    shift
	    ;;
	-h|--help)
	    cat <<EOF >&2
********************************************************************
********************************************************************
****      T-RECS, The Tiered Radio Extragalactic Continuum      ****
****              Welcome to the master script!                 ****
****                        RRRAWRRR!                           ****
********************************************************************
********************************************************************

Usage:
    $(basename $0) [OPTIONAL(s)] -p|--params [FILE]

Mandatory arguments:

    -p, --params /path/to/parameter_file
         provide the path to a T-RECS parameter file

Optional arguments:

    -c, --continuum
        whether to generate the continuum properties 
        (default=false)

    -i, --HI
        whether to generate the neutral hydrogen properties 
        (default=false)

    -x, --xmatch
        whether to cross-match continuum properties with
        HI properties in case both have been generated 
        (default=false)

    -w, --wrap
    	whether to wrap all the resulting catalogues into one
	(NOTE: the un-wrapped catalogues will be removed)

    -h, --help
        display this help message and exit

EOF
	    exit 0
	    ;;
	*)
	    echo "Unknown option $1"
	    exit 1
	    ;;
    esac
    shift
done

if [ ! -f "$params" ] ; then
    echo "Error: parameter file not provided! See usage by running"
    echo "$ $0 --help"
    exit 1
fi

#############################################################################################
# Checking python dependencies before running

# necessary for both Wrapper and X-matcher:
if [ $dowrap = true || $doxmatch = true ]; then

   # AstroPy
   python -c 'import astropy' 2>/dev/null
   if [ $? != 0 ]; then
       echo "Error: TRECS needs the AstroPy package of python to run" 1>&2
       exit 1
   fi
   
fi

# necessary only for the X-matcher:
if [ $doxmatch = true ]; then

    # SKLearn
    python -c 'import sklearn' 2>/dev/null
    if [ $? != 0 ]; then
	echo "Error: TRECS needs the SKLearn package of python to run" 1>&2
	exit 1
    fi

    # NumPy
    python -c 'import numpy' 2>/dev/null
    if [ $? != 0 ]; then
	echo "Error: TRECS needs the NumPy package of python to run" 1>&2
	exit 1
    fi

fi

#############################################################################################
# Continuum on-demand

if [ $docontinuum = true ]; then
    trecs_sampler_continuum $params
fi

#############################################################################################
# HI on-demand

if [ $dohi = true ]; then
    trecs_sampler_hi $params
fi

#############################################################################################
# Cross-match on-demand

if [ $doxmatch = true ]; then
    # trecs_sampler_xmatch $params
    echo
fi

#############################################################################################
# here do wrapping

if [ $dowrap = true ]; then
    # trecs_wrapper
    echo
fi

#############################################################################################
