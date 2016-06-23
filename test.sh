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

case $CASE in
  "AMC")
###########################
# testing preprocessor    #
###########################
        #source ../grid/sph1d/test       | tee -a test.log
        \cp ../../preprocessor/presparx/Shu1D/* ./
        presparx -o model
###########################
# testing AMC             #
###########################
        source ../amc/test      | tee -a test.log
###########################
# testing postprocessor   #
###########################
        #source ../telsim/coldens/test  | tee -a test.log
        #source ../telsim/cont/test     | tee -a test.log
        #source ../telsim/zeeman/test   | tee -a test.log
        #source ../telsim/overlap/test  | tee -a test.log
        #source ../telsim/lte/test      | tee -a test.log
        ;;
  "P2A")
        source ../benchmark/2002_p2a_benchmark/test     | tee -a test.log
        ;;
  "LTE")
        source ../grid/sph1d/test       | tee -a test.log
        source ../telsim/lte/test       | tee -a test.log
        ;;
  "Q-PARA")
        source ../benchmark/2002_p2a_benchmark/test     | tee -a test.log  
        ;;
  "CONTRIBUTION")
        #source ../telsim/contribution/sph1d/test       | tee -a test.log
        #source ../telsim/contribution/sph3d/test       | tee -a test.log
        source ../telsim/contribution/cyl3d/test       | tee -a test.log
        ;;
  "LINE")
        source ../telsim/line/test       | tee -a test.log
        ;;
  "QMC")
        source ../algorithm/AMC_accuracy/test       | tee -a test.log
        #gfortran ../algorithm/AMC_accuracy/pops_error.f90 -o pops_error
        #./pops_error
        #gnuplot ../algorithm/AMC_accuracy/fit  | tee -a test.log
        ;;
  "SOR")
        #source ../algorithm/ALI_convergency/test       | tee -a test.log
        gfortran ../algorithm/ALI_convergency/pops_error.f90 -o pops_error
        ./pops_error
        #gnuplot ../algorithm/ALI_convergency/plot  | tee -a test.log
        ;;
  "SHU1D")
        \cp ../../preprocessor/presparx/Shu1D/* ./
        presparx -o model
        ;;
  *)
        exit 1
        ;;
esac

cd ../..

