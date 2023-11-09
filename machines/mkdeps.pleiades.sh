#!/bin/bash
module load mpi-hpe/mpt
module load comp-intel

CC=icc
CXX=icpc
MPICC=icc 
MPICXX=icpc 
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
cd install-deps
export GKYLSOFT=$HOME/gkylsoft
./mkdeps.sh CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --prefix=$GKYLSOFT --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
