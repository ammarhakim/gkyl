#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/fftw-3.3.10

# delete old checkout and builds
rm -rf fftw-3.3.10.tar* fftw-3.3.10

curl -L http://www.fftw.org/fftw-3.3.10.tar.gz > fftw-3.3.10.tar.gz
gunzip fftw-3.3.10.tar.gz
tar -xvf fftw-3.3.10.tar
cd fftw-3.3.10
./configure --prefix=$PREFIX --disable-fortran CFLAGS="-fPIC" CC=$MPICC CXX=$MPICXX MPICC=$MPICC MPICXX=$MPICXX --enable-mpi
make
make install

# soft-link
ln -sfn $PREFIX $GKYLSOFT/fftw
