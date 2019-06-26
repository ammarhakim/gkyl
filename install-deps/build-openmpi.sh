#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/openmpi-3.1.2

# delete old checkout and builds
rm -rf openmpi-3.1.2.tar* openmpi-3.1.2

curl -L https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.2.tar.gz > openmpi-3.1.2.tar.gz
gunzip -c openmpi-3.1.2.tar.gz | tar xf -
cd openmpi-3.1.2
./configure --prefix=$PREFIX --enable-mpi-fortran=none CC=$CC CXX=$CXX
make all install

# soft-link 
ln -sf $PREFIX $GKYLSOFT/openmpi
