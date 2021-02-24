#!/bin/bash

module reset
#module load foss # gcc
module load intel/2019b # intel


# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='$HOME/tinkercliffs/gkylsoft'
cd install-deps
# gcc
#./mkdeps.sh CC=gcc CXX=g++ MPICC=mpicc MPICXX=mpicxx --prefix=$GKYLSOFT --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no --build-luarocks=yes 
# intel
./mkdeps.sh CC=icc CXX=icpc MPICC=mpicc MPICXX=mpicxx --prefix=$GKYLSOFT --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no --build-luarocks=yes
