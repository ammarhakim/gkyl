#!/bin/bash
##  Add the following lines to your .bashrc file:
#   export GCCBASE=/nfs-projects/opt/theta/gcc-9-toolkit
#   export PATH=$GCCBASE/bin:$PATH
#   export LD_LIBRARY_PATH=$GCCBASE/lib:$GCCBASE/lib64:$LD_LIBRARY_PATH

module load cray-python/3.6.1.1
CC=mpicc
CXX=mpicxxc
MPICC=mpicc
MPICXX=mpicxx
export GKYLSOFT=$HOME/gkylsoft
export MACHINE_NAME='theta-ga'
cd install-deps
./mkdeps.sh CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --prefix=$GKYLSOFT --build-gkylzero=yes --build-luajit=yes --build-adios=yes --build-openmpi=no
