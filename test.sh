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

# Benchmark case
    "P2A")
        #source ../benchmark/2002_p2a_benchmark/test $SAVE_LOG
        \cp ../../preprocessor/presparx/P2A/* ./
        presparx -o model -p
        
        
        # generating non-LTE level population
        printf "${YELLOW}RUNNING AMC${NC}\n"
        POPSFILE='pops_sparx'
        rm -rf $POPSFILE $POPSFILE*
        sparx run task_amc \
        source=model \
        out=$POPSFILE \
        trace='True' \
        dat='True'
        
        \cp ../benchmark/2002_p2a_benchmark/model_1.d ./
        \cp ../benchmark/2002_p2a_benchmark/pops_ratran.dat ./
        gnuplot ../benchmark/2002_p2a_benchmark/plot
        
        ;;
    "P2B")
        #source ../benchmark/2002_p2a_benchmark/test $SAVE_LOG
        \cp ../../preprocessor/presparx/P2B/* ./
        presparx -o model -p
        
        
        # generating non-LTE level population
        printf "${YELLOW}RUNNING AMC${NC}\n"
        POPSFILE='pops_sparx'
        rm -rf $POPSFILE $POPSFILE*
        sparx run task_amc \
        source=model \
        out=$POPSFILE \
        trace='True' \
        dat='True'
        
        \cp ../benchmark/2002_p2b_benchmark/model_2.d ./
        \cp ../benchmark/2002_p2b_benchmark/pops_ratran.dat ./
        gnuplot ../benchmark/2002_p2b_benchmark/plot
        
        ;;
# preprocessing
    "SHU1D")
        \cp ../../preprocessor/presparx/Shu_sph1d/* ./
        presparx -o model -epv
        ;;
    "AGB1D")
        \cp ../../preprocessor/presparx/AGB_sph1d/* ./
        presparx -o model -epv
        ;;
    "DISK_SPH2D")
        \cp ../../preprocessor/presparx/Disk_sph2d/* ./
        presparx -o model -epv
        ;;
    "DISK_CYL2D")
        \cp ../../preprocessor/presparx/Disk_cyl2d/* ./
        presparx -o model -epv
        ;;
    "N1333")
        \cp ../../preprocessor/presparx/N1333I4A/* ./
        presparx -o model 
        ;;
    "COMET2D")
        \cp ../../preprocessor/presparx/comet2D/* ./
        presparx -o model -epv
        ;;
        
    "ZEUS")
        \cp ../../preprocessor/presparx/Zeus/* ./
        presparx -o model -cepv
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
        ;;
# imaging / postprocessing
    "LINE")
        source ../telsim/line/test       $SAVE_LOG
        ;;
    "LTE")
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
        #source ../telsim/contribution/cyl3d/test      $SAVE_LOG
        source ../telsim/contribution/test      $SAVE_LOG
        ;;
    "SOURCE")
        source ../telsim/OuterSource/test     $SAVE_LOGG
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
    "VELO_INTERP")
        source ../algorithm/VELO_INTERP_accuracy/test    $SAVE_LOG
        #gfortran ../algorithm/VELO_INTERP_accuracy/pops_error.f90 -o pops_error
        #./pops_error
        #gnuplot ../algorithm/VELO_INTERP_accuracy/fit  $SAVE_LOG
        ;;

#  parallelization and queuing system 
    "Q-PARA")
        qsub ../parallel_job/test
        while [ 1 ]; do clear;qstat -u $USER;tail -n 20 history.log;sleep 5;done
        ;;

        
# visualization
    "LINECTB")
        source ../visual/linectb/sph1d/test       $SAVE_LOG
        ;;
    "CONTCTB")
        source ../visual/contctb/sph1d/test       $SAVE_LOG        
        ;;
    "MODEL2VTK")
        source ../visual/model2vtk/sph1d/test       $SAVE_LOG   
        ;;
        
    "POPS")
        source ../pops2ascii/test       $SAVE_LOG   
        ;;
# defualt : exit        
  *)
        exit 1
        ;;
esac

cd ../..

