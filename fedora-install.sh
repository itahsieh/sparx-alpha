#! /bin/bash

# INSTALL SPARX FOR FEDORA 24
sudo yum install \
redhat-rpm-config \
python-devel \
gsl-devel \
libtool \
fftw-devel \
cfitsio-devel \
hdf5-openmpi-devel \
libX11-devel \
python2-matplotlib \
python2-scipy \
sympy \
python-tables

HDF_INCLUDE='--with-include=/usr/include/hdf5/openmpi'
HDF_LIB='--with-lib=/usr/lib64/openmpi/lib'
HDF_INCLUDE='--with-include=/usr/include/openmpi-x86_64'
FITSIO_INCLUDE='--with-include=/usr/include/cfitsio'

destination=$HOME/opt/sparx
rm -rf build/* $destination

python setup.py install --prefix=$destination \
$HDF_INCLUDE \
$FITSIO_INCLUDE \
$HDF_LIB


if ! grep -q "# SPARX PATH" ~/.bashrc; then
  echo "# SPARX PATH" >> ~/.bashrc
  echo 'PATH=$PATH:$HOME/opt/sparx/bin' >> ~/.bashrc
  echo 'export PYTHONPATH=$HOME/opt/sparx/lib64/python2.7/site-packages:$PYTHONPATH' >> ~/.bashrc
fi
