#!/bin/bash
module load cray-python/3.6.1.1
CC=icc
CXX=icpc
MPICC=cc
MPICXX=CC
export GKYLSOFT=$HOME/gkylsoft
cd install-deps
./mkdeps.sh CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --prefix=$GKYLSOFT --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
