#! /bin/bash


sudo yum install \
redhat-rpm-config \
python-devel \
numpy \
gsl-devel \
libtool \
fftw-devel \
cfitsio-devel \
hdf5-openmpi-devel \
libX11-devel \
python2-matplotlib \
sympy \
python-tables

FEDORA_VERSION=`grep -o '[0-9]*' /etc/fedora-release`


if [ $FEDORA_VERSION == '22' ]; then 
  # INSTALL SPARX FOR FEDORA 22
  sudo yum install scipy 
 
  # setup bash environment
  if ! grep -q "# SPARX ENVIRONMENT" ~/.bashrc; then
    echo "# SPARX ENVIRONMENT" >> ~/.bashrc
    echo 'export LIBRARY_PATH=$LIBRARY_PATH:/usr/lib64/openmpi/lib' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/openmpi/lib' >> ~/.bashrc
  fi

elif [ $FEDORA_VERSION == '24' ]; then
  # INSTALL SPARX FOR FEDORA 24
  sudo yum install python2-scipy 

fi

HDF_LIB='--with-lib=/usr/lib64/openmpi/lib'
HDF_INCLUDE='--with-include=/usr/include/openmpi-x86_64'
FITSIO_INCLUDE='--with-include=/usr/include/cfitsio'


destination=$HOME/opt/sparx
rm -rf build/* $destination

python setup.py install --prefix=$destination \
$HDF_INCLUDE $HDF_LIB \
$FITSIO_INCLUDE


if ! grep -q "# SPARX PATH" ~/.bashrc; then
  echo "# SPARX PATH" >> ~/.bashrc
  echo 'PATH=$PATH:$HOME/opt/sparx/bin' >> ~/.bashrc
  echo 'export PYTHONPATH=$HOME/opt/sparx/lib64/python2.7/site-packages:$PYTHONPATH' >> ~/.bashrc
fi
