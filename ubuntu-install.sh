#! /bin/bash

# INSTALL SPARX FOR UBUNTU
sudo apt-get install \
libpython2.7-dev \
libhdf5-openmpi-7-dev \
libgsl0-dev \
libtool \
libfftw3-dev \
libcfitsio3-dev 

destination=$HOME/opt/sparx
rm -rf build/* $destination
python setup.py install --prefix=$destination 

if ! grep -q "# SPARX PATH" ~/.bashrc; then
  echo "# SPARX PATH" >> ~/.bashrc
  echo "PATH=$PATH:$HOME/opt/sparx/bin" >> ~/.bashrc
  echo "export PYTHONPATH=$HOME/opt/sparx/lib/python2.7/site-packages:$PYTHONPATH" >> ~/.bashrc
fi
