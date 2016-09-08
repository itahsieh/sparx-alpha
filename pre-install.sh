#!/bin/bash

TMP_DIR=$HOME/tmp
OPT_DIR=$HOME/opt

mkdir -p $TMP_DIR
cd $TMP_DIR

# OPENMPI
wget https://www.open-mpi.org/software/ompi/v2.0/downloads/openmpi-2.0.0.tar.gz  --no-check-certificate -O- | tar xz
cd openmpi-2.0.0 
./configure --prefix=$OPT_DIR/openmpi-2.0.0
make all install
cd ..

#FITSIO
wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio3390.tar.gz -O- | tar xz
cd cfitsio
./configure --prefix=$OPT_DIR/cfitsio
make shared
make install
cd ..

# HDF5
wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.17.tar -O- | tar x
cd hdf5-1.8.17 
CC=$OPT_DIR/openmpi-2.0.0/bin/mpicc ./configure --enable-parallel --prefix=$OPT_DIR/hdf5-1.8.17
make 
make check 
make install
cd ..
 




cd ..

rm -rf $TMP_DIR
