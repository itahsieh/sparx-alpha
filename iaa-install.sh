#! /bin/bash

# check bashrc to source SPARX's module
if ! grep -q "# SPARX ENVIROMENT" ~/.bashrc ; then
  echo "# SPARX ENVIROMENT" >> ~/.bashrc
  echo "source /asiaa/home/ithsieh/.SPARX_dependent" \
        >> ~/.bashrc
  source $HOME/.bashrc
fi

destination=$HOME/opt/sparx
rm -rf build/* $destination
python setup.py install \
--prefix=$destination 

if ! grep -q "# SPARX PATH" ~/.bashrc; then
  echo "# SPARX PATH" >> ~/.bashrc
  echo "PATH=$PATH:$HOME/opt/sparx/bin" >> ~/.bashrc
  echo "export PYTHONPATH=$HOME/opt/sparx/lib/python2.7/site-packages:$PYTHONPATH" >> ~/.bashrc
fi
