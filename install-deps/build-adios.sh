#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/adios2-2.9.2

# delete old checkout and builds
rm -rf adios2-2.9.2.tar* adios2-2.9.2

curl -L https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.9.2.tar.gz > adios2-2.9.2.tar.gz
gunzip adios2-2.9.2.tar.gz
mkdir adios2-2.9.2
tar -xvf adios2-2.9.2.tar -C adios2-2.9.2 --strip-components 1

cd adios2-2.9.2
mkdir build
cd build

CC=$MPICC CXX=$MPICXX cmake ../ -DCMAKE_C_FLAGS="-g -O3 -fPIC" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX -DBUILD_TESTING=OFF -DADIOS2_BUILD_EXAMPLES=OFF -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_HDF5=OFF -DADIOS2_USE_Python=OFF -DBUILD_SHARED_LIBS=ON 

make -j 32 VERBOSE=1
make install

# soft-link
ln -sfn $PREFIX $GKYLSOFT/adios2
