#! /bin/bash
# distinatin of the binary path
destination=$HOME/opt/sparx


# the color label
LIGHTBLUE='\033[1;34m'
YELLOW='\033[1;33m'
LIGHTCYAN='\033[1;36m'
RED='\033[0;31m'
NC='\033[0m' # No Color


# CHECK if LOAD_MODULE_FILE exits
LOAD_MODULE_FILE=~/.load_sparx_module
if [ ! -e $LOAD_MODULE_FILE ];then
  cp load_tiara_module.sh $LOAD_MODULE_FILE
  echo 'PATH=$PATH:'$destination'/bin' | tee -a $LOAD_MODULE_FILE
  echo 'export PYTHONPATH=$PYTHONPATH:'$destination'/lib/python2.7/site-packages' | tee -a $LOAD_MODULE_FILE
  printf "${LIGHTBLUE}COPY load_tiara_module.sh TO $LOAD_MODULE_FILE ${NC}\n"
else
  printf "${LIGHTBLUE}PLEASE manually COPY load_tiara_module.sh TO $LOAD_MODULE_FILE ${NC}\n"
fi

# LOAD MODULE AND DEFINE SPARXVERSION
source $LOAD_MODULE_FILE
printf "${YELLOW}LOAD MODULE${NC}\n"

# check bashrc to source SPARX's module
if ! grep -q "# SPARX ENVIROMENT" ~/.bashrc ; then
  echo '###### SPARX START ######' >> ~/.bashrc
  echo '# SPARX PATH' >> ~/.bashrc
  echo 'if [ ! -z "$CLUSTERNAME" ]' >> ~/.bashrc
  echo 'then' >> ~/.bashrc
  echo '  # SPARX ENVIROMENT' >> ~/.bashrc
  echo '  source $LOAD_MODULE_FILE' >> ~/.bashrc
  echo 'fi' >> ~/.bashrc
  echo '###### SPARX END ######' >> ~/.bashrc
  printf "${RED}UPDATE ~/.bashrc${NC}\n"
else
  echo '###### SPARX START ######'
  echo '# SPARX PATH'
  echo 'if [ ! -z "$CLUSTERNAME" ]'
  echo 'then'
  echo '  # SPARX ENVIROMENT'
  echo '  source ~/.load_sparx_module'
  echo 'fi'
  echo '###### SPARX END ######'
  printf "${LIGHTBLUE}manually UPDATE above script into ~/.bashrc, if needed${NC}\n"
fi



# INSTALLATION
rm -rf build/*
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
--with-lib=$OPENMPI_HOME/lib \
--with-lib=$CFITSIO_HOME/lib \
--with-lib=$MIR/lib/linux \
&& printf "${LIGHTCYAN}BUILDING IS DONE!${NC}\n"

if [ $? -eq 1 ]; then
  exit $?
fi

# REPLACE the first 'sparx' occurrence in bin/sparx TO 'sparx-CLUSTERNAME'
EXECUTABLE=$destination/bin/$SPARXVERSION
sed -e "s/import sparx/import $SPARX_VERSION/" \
    -e "s/from sparx/from $SPARX_VERSION/" \
        $destination/bin/sparx > $EXECUTABLE
chmod 755 $EXECUTABLE
rm $destination/bin/sparx
printf "${RED}CREATE bin/${SPARXVERSION}${NC}\n"


# REPLACE the first 'sparx' occurrence in bin/presparx TO 'sparx-CLUSTERNAME'
EXECUTABLE=$destination/bin/pre$SPARXVERSION
sed -e "s/import sparx/import $SPARX_VERSION/" \
    -e "s/from sparx/from $SPARX_VERSION/" \
        bin/presparx > $EXECUTABLE
chmod 755 $EXECUTABLE
rm $destination/bin/presparx
printf "${RED}CREATE bin/pre${SPARXVERSION}${NC}\n"



# REPLACE sparx_module in FILE to 'sparx_CLUSTERNAME' 
SPARX_PYTHONPATH=$destination/lib/python2.7/site-packages
for FILE in `ls lib/sparx | grep .py`;do
  sed -e "s/import sparx/import $SPARX_VERSION/" \
      -e "s/from sparx/from $SPARX_VERSION/" \
          lib/sparx/$FILE > $SPARX_PYTHONPATH/$SPARX_VERSION/$FILE
done
printf "${LIGHTBLUE}CREATE $SPARX_VERSION/*.py${NC}\n"





