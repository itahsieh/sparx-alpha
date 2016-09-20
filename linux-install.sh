#! /bin/bash

# check bashrc to source SPARX's module
BASHRC=$HOME/.bashrc
SPARX_ENV=$HOME/.SPARX_dependent
if ! grep -q "# SPARX ENVIROMENT" $BASHRC ; then
  cp SPARX_ENV $SPARX_ENV
  echo "# SPARX ENVIROMENT" >> $BASHRC
  echo "source $SPARX_ENV"  >> $BASHRC
  source $SPARX_ENV
fi

destination=$HOME/opt/sparx
rm -rf build/* $destination

python setup.py install \
--prefix=$destination \
--with-include=$OPENMPI_PATH 

if ! grep -q "# SPARX PATH" $BASHRC; then
  echo "# SPARX PATH" >> $BASHRC
  echo 'PATH=$PATH:'$destination/bin >> $BASHRC
  echo 'export PYTHONPATH='$destination'/lib/python2.7/site-packages:$PYTHONPATH' >> $BASHRC
fi
