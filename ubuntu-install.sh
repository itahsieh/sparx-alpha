#! /bin/bash
destination=$HOME/opt/sparx
rm -rf $HOME/pysparx/build/* $destination
python setup.py install --prefix=$destination 


