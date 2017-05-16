#!/bin/bash

. /etc/profile.d/modules.sh
module purge

module add icc/16.0
module add cfitsio
module add openmpi/2.0.1_ic16.0
module add fftw/3.3.4_openmpi_2.0.1_ic16.0
module add python/2.7.11
module add hdf5/1.8.16_openmpi_2.0.1_ic16.0 
module add gsl/2.1_ic16.0
module add miriad



HOSTNAME=`hostname`
CLUSTERNAME=${HOSTNAME:0:2}
SPARXVERSION='sparx-'$CLUSTERNAME
SPARX_VERSION='sparx_'$CLUSTERNAME
PRESPARX='presparx-'$CLUSTERNAME
# SPARX PATH
PATH=$PATH:$HOME/opt/sparx/bin
export PYTHONPATH=$PYTHONPATH:$HOME/opt/sparx/lib/python2.7/site-packages
alias sparx=$SPARXVERSION
alias presparx=$PRESPARX
