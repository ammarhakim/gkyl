#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/adios-1.13.1

# Tar file is checked in to repo
gunzip adios-1.13.1.tar.gz
tar -xvf adios-1.13.1.tar
cd adios-1.13.1
./configure --prefix=$PREFIX --disable-fortran --without-netcdf CFLAGS="-fPIC" CC=$MPICC CXX=$MPICXX MPICC=$MPICC MPICXX=$MPICXX --enable-shared=no
make install

# soft-link
ln -sfn $PREFIX $GKYLSOFT/adios
