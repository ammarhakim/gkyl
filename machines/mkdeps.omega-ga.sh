#!/bin/bash
##  Add the following lines to your .bashrc file:
#   export GCCBASE=/nfs-projects/opt/theta/gcc-9-toolkit
#   export PATH=$GCCBASE/bin:$PATH
#   export LD_LIBRARY_PATH=$GCCBASE/lib:$GCCBASE/lib64:$LD_LIBRARY_PATH

CC=mpicc
CXX=mpicxxc
MPICC=mpicc
MPICXX=mpicxx
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
