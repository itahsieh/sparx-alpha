#!/bin/bash

# configure SPARX command
sparx_version=`which sparx`

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

# testing preprocessor
source ../grid/sph1d/test | tee -a test.log

# testing AMC
source ../amc/test  | tee -a test.log

# testing postprocessor
source ../telsim/line/test     | tee -a test.log
source ../telsim/coldens/test  | tee -a test.log
source ../telsim/cont/test     | tee -a test.log
source ../telsim/zeeman/test   | tee -a test.log
source ../telsim/overlap/test  | tee -a test.log
source ../telsim/lte/test      | tee -a test.log


cd ../..

