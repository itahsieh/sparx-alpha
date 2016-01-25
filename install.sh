#! /bin/bash

dev=1

if [ $dev == 1 ]; then
  destination=$HOME/opt/sparx-dev
else
  destination=$HOME/opt/sparx
fi

#rm -rf /home/users/schung/usr/lib/python2.5/site-packages/sparx* /home/users/schung/usr/bin/sparx* #$HOME/opt/sparx
rm -rf $HOME/pysparx/build/* $destination

python setup.py install \
--prefix=$destination --lam \
--with-include=$FFTW_HOME/include/ \
--with-include=$HDF5_HOME/include/ \
--with-include=$LAM_HOME/include/ \
--with-include=$CFITSIO_HOME/include/ \
--with-lib=$GSL_HOME/lib \
--with-lib=$FFTW_HOME/lib \
--with-lib=$HDF5_HOME/lib \
--with-lib=$LAM_HOME/lib \
--with-lib=$CFITSIO_HOME/lib

if [ $dev == 1 ]; then
  sed -e 's/import sparx._sparx as _sparx/import sparxdev._sparx as _sparx/g' -e 's/from sparx.tasks/from sparxdev.tasks/g' $destination/bin/sparx > $destination/bin/sparx-dev
  chmod 755 $destination/bin/sparx-dev
  rm $destination/bin/sparx
  sed -e 's/from sparx import _sparx/from sparxdev import _sparx/g' $HOME/pysparx/lib/sparx/tasks.py > $destination/lib/python2.5/site-packages/sparxdev/tasks.py
fi

