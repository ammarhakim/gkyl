#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/openmpi-3.0.0

# delete old checkout and builds
rm -rf openmpi-3.0.0.tar* openmpi-3.0.0

curl -L https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz > openmpi-3.0.0.tar.gz
gunzip -c openmpi-3.0.0.tar.gz | tar xf -
cd openmpi-3.0.0
./configure --prefix=$PREFIX --enable-mpi-fortran=none
make all install

# write information to build-opts
echo "# MPI library locations

MPICC=$PREFIX/bin/mpicc
MPICXX=$PREFIX/bin/mpicxx
" >> build-opts.sh
