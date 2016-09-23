#!/bin/bash

###########################
# Swith testing cases     #
# The available cases are #
# AMC/P2A/LTE/Q-PARA      #
###########################
CASE=${1} 

# the color label
LIGHTBLUE='\033[1;34m'
YELLOW='\033[1;33m'
LIGHTCYAN='\033[1;36m'
RED='\033[0;31m'
NC='\033[0m' # No Color



cd unit_tests
# create "tmp" directory only if it doesn't exist. then enter it
mkdir -p tmp; cd tmp
# remove the testing log file
rm -f test.log

SAVE_LOG='| tee -a test.log'

case $CASE in

# preprocessing
  "P2A")
        source ../benchmark/2002_p2a_benchmark/test $SAVE_LOG
        ;;
  "SHU1D")
        \cp ../../preprocessor/presparx/Shu1D/* ./
        presparx -o model -e
        ;;
  "DISK_SPH2D")
        \cp ../../preprocessor/presparx/Disk_sph2d/* ./
        presparx -o model -v -p
        ;;
  "DISK_CYL2D")
        \cp ../../preprocessor/presparx/Disk_cyl2d/* ./
        presparx -o model -v -p
        ;;

# AMC solver
  "AMC")
###########################
# testing preprocessor    #
###########################
        #source ../grid/sph1d/test $SAVE_LOG
###########################
# testing AMC             #
###########################
        source ../amc/test $SAVE_LOG
###########################
# testing postprocessor   #
###########################
        #source ../telsim/coldens/test  $SAVE_LOG
        #source ../telsim/cont/test     $SAVE_LOG
        #source ../telsim/zeeman/test   $SAVE_LOG
        #source ../telsim/overlap/test  $SAVE_LOG
        #source ../telsim/lte/test      $SAVE_LOG
        
# imaging / postprocessing
        ;;
  "LINE")
        source ../telsim/line/test       $SAVE_LOG
        ;;
  "LTE")
        source ../grid/sph1d/test       $SAVE_LOG
        source ../telsim/lte/test       $SAVE_LOG
        ;;
  "ZEEMAN")
        source ../telsim/zeeman/test       $SAVE_LOG
        ;;
  "CONT")
        source ../telsim/cont/test       $SAVE_LOG
        ;;
  "COLDENS")
        source ../telsim/coldens/test   $SAVE_LOG
        ;;
  "CONTRIBUTION")
        #source ../telsim/contribution/sph1d/test     $SAVE_LOG
        #source ../telsim/contribution/sph3d/test     $SAVE_LOG
        source ../telsim/contribution/cyl3d/test      $SAVE_LOG
        ;;
 

        
# Algorithm testing
  "QMC")
        source ../algorithm/AMC_accuracy/test       $SAVE_LOG
        #gfortran ../algorithm/AMC_accuracy/pops_error.f90 -o pops_error
        #./pops_error
        #gnuplot ../algorithm/AMC_accuracy/fit  $SAVE_LOG
        ;;
  "SOR")
        #source ../algorithm/ALI_convergency/test    $SAVE_LOG
        gfortran ../algorithm/ALI_convergency/pops_error.f90 -o pops_error
        ./pops_error
        #gnuplot ../algorithm/ALI_convergency/plot  $SAVE_LOG
        ;;

#  parallelization and queuing system 
  "Q-PARA")
        source ../benchmark/2002_p2a_benchmark/test    $SAVE_LOG 
        ;;        

# defualt : exit        
  *)
        exit 1
        ;;
esac

cd ../..

