#!/bin/bash

. /etc/profile.d/modules.sh
module purge

module add cfitsio
module add python/2.7.11
module add miriad

# module add openmpi/2.0.1_ic16.0
# module add fftw/3.3.4_openmpi_2.0.1_ic16.0
# module add hdf5/1.8.16_openmpi_2.0.1_ic16.0 
# module add gsl/2.1_ic16.0
# module add icc/16.0

module add hdf5/1.8.16_openmpi_1.10.2_ic15.0
module add openmpi/1.10.2_ic15.0
module add icc/15.0
module add gsl/2.1_ic15.0
module add fftw/3.3.4_openmpi_1.10.2_ic15.0

# redirect sparx to cluster version
HOSTNAME=`hostname`
CLUSTERNAME=${HOSTNAME:0:2}
SPARXVERSION='sparx-'$CLUSTERNAME
SPARX_VERSION='sparx_'$CLUSTERNAME
PRESPARX='presparx-'$CLUSTERNAME

alias sparx=$SPARXVERSION
alias presparx=$PRESPARX

