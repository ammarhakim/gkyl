#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/adios-1.13.0

# delete old checkout and builds
rm -rf adios-1.13.0.tar* adios-1.13.0

curl -L http://users.nccs.gov/~pnorbert/adios-1.13.0.tar.gz > adios-1.13.0.tar.gz
gunzip adios-1.13.0.tar.gz
tar -xvf adios-1.13.0.tar
cd adios-1.13.0
./configure --prefix=$PREFIX --disable-fortran CFLAGS="-fPIC" CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --enable-shared=no
make install

# soft-link
ln -sf $PREFIX $GKYLSOFT/adios
