#!/bin/bash
module load mpi-sgi/mpt
module load comp-intel

CC=icc
CXX=icpc
MPICC=icc 
MPICXX=icpc 
export GKYLSOFT='~/gkylsoft'
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
cd install-deps
./mkdeps.sh CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
