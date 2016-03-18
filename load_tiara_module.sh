#!/bin/bash

. /etc/profile.d/modules.sh
module purge

HOSTNAME=`hostname`
CLUSTERNAME=${HOSTNAME:0:2}
SPARXVERSION='sparx-'$CLUSTERNAME
SPARX_VERSION='sparx_'$CLUSTERNAME

module add icc
module add torque
module add cfitsio
module add openmpi
module add fftw

if [ "$HOSTNAME" == "gwhpc" -o "$CLUSTERNAME" == "oc" -o "$CLUSTERNAME" == "tc" ];then 
    module add python/2.7.2
    module add miriad/2011.7_gfortran
    module add HDFView/2.1
    module add HDF/5-1.8.10_ic13.0_openmpi_1.6.3
    module add gsl/1.13_ic11.0
    module add git

elif [ "$HOSTNAME" == "ibhpc" -o "$CLUSTERNAME" == "px" ];then 
    module add python/2.7.2
    module add miriad/2011.7_gfortran
    module add HDF/5-1.8.10_ic13.0_openmpi_1.6.3
    module add gsl/1.13_ic11.0
    module add git

elif [ "$HOSTNAME" == "ashpc" -o "$CLUSTERNAME" == "xl" ];then
    module add python/2.7.11
    module add hdf5/1.8.16_openmpi_1.10.2_ic15.0
    module add gsl/1.16_ic15.0
    
fi


