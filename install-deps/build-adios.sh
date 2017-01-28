#!/bin/bash

# Edit to suite your system
PREFIX=$HOME/gkylsoft/adios
PHDF5_INC_DIR=$HOME/software/hdf5-par/include
PHDF5_LIB_DIR=$HOME/software/hdf5-par/lib
PHDF5_LIBS=hdf5

#curl http://users.nccs.gov/~pnorbert/adios-1.10.0.tar.gz > adios-1.10.0.tar.gz
#gunzip adios-1.10.0.tar.gz
#tar -xvf adios-1.10.0.tar
cd adios-1.10.0
CFLAGS="-fPIC" ./configure --prefix=$PREFIX --disable-fortran
make install


#CFLAGS="-fPIC" ./configure --prefix=$PREFIX --disable-fortran --with-phdf5-incdir=$PHDF5_INC_DIR --with-phdf5-libdir=$PHDF5_LIB_DIR
