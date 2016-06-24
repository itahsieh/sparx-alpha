#! /bin/bash

# INSTALL SPARX FOR UBUNTU
sudo apt-get install \
libpython2.7-dev \
libgsl0-dev \
libtool \
libfftw3-dev \
libcfitsio3-dev \
libhdf5-openmpi-dev \
libx11-dev \
python-matplotlib \
python-sympy \
python-tables

CODENAME=`lsb_release -c | awk 'END {print $NF}'`
if   [ $CODENAME == 'xenial' ]; then
  HDF_INCLUDE='--with-include=/usr/include/hdf5/openmpi'
  MPI_INCLUDE='--with-include=/usr/lib/openmpi/include'
  HDF_LIB='--with-lib=/usr/lib/x86_64-linux-gnu/hdf5/openmpi'
elif [ $CODENAME == 'trusty' ]; then
  HDF_INCLUDE='--with-include=/usr/include/openmpi'
fi


destination=$HOME/opt/sparx
rm -rf build/* $destination

python setup.py install --prefix=$destination \
--with-include=/usr/include \
--with-lib=/usr/lib \
$HDF_INCLUDE \
$MPI_INCLUDE \
$HDF_LIB


if ! grep -q "# SPARX PATH" ~/.bashrc; then
  echo "# SPARX PATH" >> ~/.bashrc
  echo 'PATH=$PATH:$HOME/opt/sparx/bin' >> ~/.bashrc
  echo 'export PYTHONPATH=$HOME/opt/sparx/lib/python2.7/site-packages:$PYTHONPATH' >> ~/.bashrc
fi