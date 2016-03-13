#!/bin/bash

. /etc/profile.d/modules.sh
module purge

HOSTNAME=`hostname`
CLUSTERNAME=${HOSTNAME:0:2}
SPARXVERSION='sparx-'$CLUSTERNAME
SPARX_VERSION='sparx_'$CLUSTERNAME


if [ "$HOSTNAME" == "gwhpc" -o "$CLUSTERNAME" == "oc" -o "$CLUSTERNAME" == "tc" ];then 
    module add torque
    module add git
    module add autoconf/2.68
    module add automake/1.11.1
    module add python/2.5.6
    module add HDFView/2.1
    module add cfitsio/3.181
    module add miriad/2011.7_gfortran
    module add HDF/5-1.8.10_ic13.0_openmpi_1.6.3
    module add gsl/1.5_ic11.0
    module add libtool/2.4
    module add fftw/3.3_ic11.0

elif [ "$HOSTNAME" == "ibhpc" -o "$CLUSTERNAME" == "px" ];then 
    module add torque/4.2.2
    module add python/2.5.6
    module add cfitsio/3.181
    module add miriad/2011.7_gfortran
    module add HDF/5-1.8.10_ic13.0_openmpi_1.6.3
    module add gsl/1.5_ic11.0
    module add libtool/2.4
    module add fftw/3.3_ic11.0

elif [ "$HOSTNAME" == "ashpc" -o "$CLUSTERNAME" == "xl" ];then
    module add cfitsio/3.380
    module add openmpi/1.10.2_ic15.0
    module add fftw/3.3.4_openmpi_1.10.2_ic15.0

fi


