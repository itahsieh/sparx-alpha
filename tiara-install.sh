#! /bin/bash

LOAD_MODULE_FILE=~/.load_sparx_module

# LOAD MODULE AND DEFINE SPARXVERSION
source $LOAD_MODULE_FILE
printf "${YELLOW}LOAD MODULE${NC}\n"

destination=$HOME/opt/$SPARXVERSION

# the color label
LIGHTBLUE='\033[1;34m'
YELLOW='\033[1;33m'
LIGHTCYAN='\033[1;36m'
RED='\033[0;31m'
NC='\033[0m' # No Color


# CHECK if LOAD_MODULE_FILE exits
if [ ! -e $LOAD_MODULE_FILE ];then
  cp load_tiara_module.sh $LOAD_MODULE_FILE
  printf "${LIGHTBLUE}COPY load_tiara_module.sh TO $LOAD_MODULE_FILE${NC}\n"
fi


# INSTALLATION
rm -rf build/* $destination
python setup.py install \
--prefix=$destination \
--version=$CLUSTERNAME \
--with-include=$FFTW_HOME/include/ \
--with-include=$HDF5_HOME/include/ \
--with-include=$OPENMPI_HOME/include/ \
--with-include=$CFITSIO_HOME/include/ \
--with-lib=$GSL_HOME/lib \
--with-lib=$FFTW_HOME/lib \
--with-lib=$HDF5_HOME/lib \
--with-lib=$LAM_HOME/lib \
--with-lib=$CFITSIO_HOME/lib
printf "${LIGHTCYAN}BUILDING IS DONE!${NC}\n"


# REPLACE sparx_module in FILE to 'sparx_CLUSTERNAME' 
FILE='__init__.py'
if [ $CLUSTERNAME == 'xl' -o  $HOSTNAME == 'ashpc' ];then
  PYTHON_NAME='python2.7'
else
  PYTHON_NAME='python2.5'
fi
ABS_INTI_PATH="$destination/lib/$PYTHON_NAME/site-packages/$SPARX_VERSION/$FILE"
sed -e "s/SPARX_VERSION/$SPARX_VERSION/g" \
        lib/sparx/$FILE > $ABS_INTI_PATH
printf "${LIGHTBLUE}CREATE ${FILE}${NC}\n"


# REPLACE the first 'sparx' occurrence in bin/sparx TO 'sparx-CLUSTERNAME'
EXECUTABLE=$destination/bin/$SPARXVERSION
sed -e "s/SPARX_VERSION/$SPARX_VERSION/g" \
        bin/sparx > $EXECUTABLE
chmod 755 $EXECUTABLE
rm $destination/bin/sparx
printf "${RED}CREATE bin/${SPARXVERSION}${NC}\n"


mv $destination/lib/$PYTHON_NAME/site-packages/sparx/_sparx.so \
        $destination/lib/$PYTHON_NAME/site-packages/$SPARX_VERSION/_sparx.so
rm -rf $destination/lib/$PYTHON_NAME/site-packages/sparx
        
        
# check LOAD_MODULE_FILE to set SPARXVERSION's PATH
if ! grep -q "# $SPARXVERSION PATH" $LOAD_MODULE_FILE; then
  echo "# $SPARXVERSION PATH" \
        >> $LOAD_MODULE_FILE
  echo 'PATH=$PATH:$HOME/opt/'$SPARXVERSION'/bin' \
        >> $LOAD_MODULE_FILE
  echo 'export PYTHONPATH=$PYTHONPATH:'$HOME'/opt/'$SPARXVERSION'/lib/python2.5/site-packages/' \
        >> $LOAD_MODULE_FILE
  echo 'alias sparx=$SPARXVERSION' \
        >> $LOAD_MODULE_FILE
  printf "${LIGHTCYAN}UPDATE $LOAD_MODULE_FILE${NC}\n"
fi


# check bashrc to source SPARX's module
if ! grep -q "# SPARX ENVIROMENT" ~/.bashrc ; then
  echo "# SPARX ENVIROMENT" >> ~/.bashrc
  echo "source $LOAD_MODULE_FILE" >> ~/.bashrc
  printf "${RED}UPDATE ~/.bashrc${NC}\n"
fi



