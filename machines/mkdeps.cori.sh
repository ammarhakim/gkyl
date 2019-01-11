#!/bin/bash
module load cray-mpich
CC=icc
CXX=icpc
MPICC=cc
MPICXX=CC
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$(readlink -f ../gkylsoft)
cd install-deps
./mkdeps.sh --prefix=$GKYLSOFT CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --build-luajit=yes --build-adios=no --build-eigen=no --build-openmpi=no
