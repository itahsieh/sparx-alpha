#!/bin/bash

. /etc/profile.d/modules.sh
module purge

HOSTNAME=`hostname`
CLUSTERNAME=${HOSTNAME:0:2}
SPARXVERSION='sparx-'$CLUSTERNAME
SPARX_VERSION='sparx_'$CLUSTERNAME


module add git
module add icc
module add torque
module add cfitsio
module add openmpi
module add gsl
module add fftw

if [ "$HOSTNAME" == "gwhpc" -o "$CLUSTERNAME" == "oc" -o "$CLUSTERNAME" == "tc" ];then 
    module add python/2.5.6
    module add HDFView/2.1
    module add HDF/5-1.8.10_ic13.0_openmpi_1.6.3
    
elif [ "$HOSTNAME" == "ibhpc" -o "$CLUSTERNAME" == "px" ];then 
    module add python/2.5.6
    module add miriad/2011.7_gfortran
    module add HDF/5-1.8.10_ic13.0_openmpi_1.6.3

elif [ "$HOSTNAME" == "ashpc" -o "$CLUSTERNAME" == "xl" ];then
    module add python/2.7.11
    
fi


