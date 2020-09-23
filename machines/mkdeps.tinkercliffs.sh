#!/bin/bash

module reset
module load foss
#module load gcc/8.2.0
#module load openmpi/gcc/64/4.0.4

# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='$HOME/tinkercliffs/gkylsoft'
cd install-deps
./mkdeps.sh CC=gcc CXX=g++ MPICC=mpicc MPICXX=mpicxx --prefix=$GKYLSOFT --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no --build-luarocks=yes
