#!/bin/bash

# initiate and purge
. /etc/profile.d/modules.sh
module purge

# PYTHON enviroment
module add python/2.7.16

# The computing libraries
module add gsl/2.5
module add openmpi/3.1.4 
module add fftw/3.3.8_openmpi_3.1.4

# Default SPARX format: HDF
module add hdf5/1.10.5_openmpi_3.1.4

# Additional dependencies
# module add intel/2015

# The image output format library
module add cfitsio
#module add miriad/2011

# redirect sparx to cluster version
HOSTNAME=`hostname`
CLUSTERNAME=${HOSTNAME:0:2}
SPARXVERSION='sparx-'$CLUSTERNAME
SPARX_VERSION='sparx_'$CLUSTERNAME
PRESPARX='presparx-'$CLUSTERNAME

alias sparx=$SPARXVERSION
alias presparx=$PRESPARX

