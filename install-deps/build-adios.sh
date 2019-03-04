#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/adios-1.13.1

# delete old checkout and builds
rm -rf adios-1.13.1.tar* adios-1.13.1

curl -L http://users.nccs.gov/~pnorbert/adios-1.13.1.tar.gz > adios-1.13.1.tar.gz
gunzip adios-1.13.1.tar.gz
tar -xvf adios-1.13.1.tar
cd adios-1.13.1
./configure --prefix=$PREFIX --disable-fortran --without-netcdf CFLAGS="-fPIC" CC=$MPICC CXX=$MPICXX MPICC=$MPICC MPICXX=$MPICXX --enable-shared=no
make install

# soft-link
ln -sf $PREFIX $GKYLSOFT/adios
