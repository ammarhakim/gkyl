#!/bin/bash
module load cray-mpich
CC=icc
CXX=icpc
MPICC=cc
MPICXX=CC
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --build-luajit=yes --build-adios=no --build-eigen=no --build-openmpi=no
