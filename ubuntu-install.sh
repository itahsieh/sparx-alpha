#! /bin/bash

CODENAME=`lsb_release -c | awk 'END {print $NF}'`

# INSTALL SPARX DEPENDENCIES ON UBUNTU
if [ $CODENAME == 'bionic' ] || [ $CODENAME == 'xenial' ] || [ $CODENAME == 'trusty' ]; then
	echo 'UBUNTU VERSION:'$CODENAME
else
	echo 'unsupported ubuntu version:'$CODENAME
	exit 1
fi

sudo apt-get install \
	libpython2.7-dev \
	libgsl0-dev \
	libtool \
	libfftw3-dev \
	libhdf5-openmpi-dev \
	libx11-dev \
	python-matplotlib \
	python-sympy \
	python-tables

if [ $CODENAME == 'bionic' ]; then
	sudo apt-get install \
		libcfitsio-dev \
		libopenmpi-dev
elif [ $CODENAME == 'xenial' ] || [ $CODENAME == 'trusty' ]; then
	sudo apt-get install libcfitsio3-dev
fi


# ENVIRONMENT VARIABLES
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/usr/include/hdf5/openmpi:/usr/lib/openmpi/include
export LIBRARY_PATH=$LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/hdf5/openmpi

if [ $CODENAME == 'bionic' ]; then
  MPI_INCLUDE='--with-include=/usr/lib/x86_64-linux-gnu/openmpi/include'
elif [ $CODENAME == 'xenial' ]; then
  HDF_INCLUDE='--with-include=/usr/include/hdf5/openmpi'
  MPI_INCLUDE='--with-include=/usr/lib/openmpi/include'
  HDF_LIB='--with-lib=/usr/lib/x86_64-linux-gnu/hdf5/openmpi'
elif [ $CODENAME == 'trusty' ]; then
  HDF_INCLUDE='--with-include=/usr/include/openmpi'
fi


# SPARX INSTALLATION
destination=$HOME/opt/sparx
rm -rf build/* $destination

python setup.py install --prefix=$destination \
--with-include=/usr/include \
--with-lib=/usr/lib \
$HDF_INCLUDE \
$MPI_INCLUDE \
$HDF_LIB


# BASHRC APPENDIND
if ! grep -q "# SPARX PATH" ~/.bashrc; then
  echo "# SPARX PATH" >> ~/.bashrc
  echo 'PATH=$PATH:$HOME/opt/sparx/bin' >> ~/.bashrc
  echo 'export PYTHONPATH=$HOME/opt/sparx/lib/python2.7/site-packages:$PYTHONPATH' >> ~/.bashrc
fi
