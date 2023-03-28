#!/bin/bash

#############################################################################################
# Argument parsing

docontinuum=false
dohi=false
doxmatch=false
doclustering=false
dowrap=false
wraptag=none
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
	-C|--clustering)
	    doclustering=true
	    ;;
	-w|--wrap)
	    dowrap=true
	    wraptag=$2
	    case $wraptag in
		"-"*|"")
		    wraptag=none
		    ;;
		HI|HI_clustered)
		    shift
		    ;;
		continuum|continuum_clustered)
		    shift
		    ;;
		HI_continuum|HI_continuum_clustered)
		    shift
		    ;;
		*)
		    echo "Error: provided wrap-tag '$wraptag' is not valid" 1>&2
		    exit 1
		    ;;
	    esac
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

    -C, --clustering
    	whether to add clustering properties to a catalogue
	NOTE that this will be possible only for catalogues with
	size <= 5 degrees. Which catalogue to run this on
	is controlled through the parameter file

    -w, --wrap [tag]
    	whether to wrap all the resulting catalogues into one.
	Argument 'tag' is optional, if none is passed it will
	either be read from the parameter file or set to its 
	default value (i.e. 'continuum')

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

# necessary only for the X-matcher:
if [ $doxmatch = true ]; then

    # AstroPy
    python -c 'import astropy' 2>/dev/null
    if [ $? != 0 ]; then
	echo "Error: TRECS needs the AstroPy package of python to run X-matching" 1>&2
	exit 1
    fi

    # SKLearn
    python -c 'import sklearn' 2>/dev/null
    if [ $? != 0 ]; then
	echo "Error: TRECS needs the SKLearn package of python to run X-matching" 1>&2
	exit 1
    fi

    # NumPy
    python -c 'import numpy' 2>/dev/null
    if [ $? != 0 ]; then
	echo "Error: TRECS needs the NumPy package of python to run X-matching" 1>&2
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
# Cross-match HI on-demand

if [ $doxmatch = true ]; then
    trecs_xmatch_hi $params
fi

#############################################################################################
# Add clustering on-demand

if [ $doclustering = true ]; then
    trecs_xmatch_clustering $params
fi

#############################################################################################
# here do wrapping

if [ $dowrap = true ]; then
    if [ $wraptag = none ]; then
        trecs_wrapper --params $params
    else
	trecs_wrapper --params $params --tag $wraptag
    fi
fi

#############################################################################################
