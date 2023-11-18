#!/bin/bash

# Edit the paths and options in the following command to suit your system
module load intel/2021.1.2
module load openmpi/intel-2021.1/4.1.0
module load anaconda3/2020.11

# Build directory
OUT=build
# Install location
GKYLSOFT=$HOME/gkylsoft

# Compile flags (set optimization/debug flags here)
CC=mpicc
CXX=mpicxx
CXXFLAGS='-O3,-std=c++17'

# LuaJIT options
LUAJIT_INC_DIR=$GKYLSOFT/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$GKYLSOFT/luajit/lib
LUAJIT_SHARE_DIR=$GKYLSOFT/luajit/share/luajit-2.1.0-beta3

## MPI options
MPICC=mpicc
MPICXX=mpicxx
ENABLE_MPI="--enable-mpi"
I_MPI_ROOT=/usr/local/openmpi/4.1.0/intel20211
MPI_INC_DIR=$I_MPI_ROOT/include
MPI_LIB_DIR=$I_MPI_ROOT/lib64
MPI_LINK_LIBS="mpi"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$GKYLSOFT/adios/include
ADIOS_LIB_DIR=$GKYLSOFT/adios/lib

OPENBLAS_INC_DIR=$GKYLSOFT/OpenBLAS/include
OPENBLAS_LIB_DIR=$GKYLSOFT/OpenBLAS/lib

SUPERLU_INC_DIR=$GKYLSOFT/superlu/include
SUPERLU_LIB_DIR=$GKYLSOFT/superlu/lib

GKYLZERO_INC_DIR=$GKYLSOFT/gkylzero/include
GKYLZERO_LIB_DIR=$GKYLSOFT/gkylzero/lib

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT -p $GKYLSOFT --prefix=$GKYLSOFT/gkyl --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --enable-gkylzero --gkylzero-inc-dir=$GKYLZERO_INC_DIR --gkylzero-lib-dir=$GKYLZERO_LIB_DIR --enable-superlu --superlu-inc-dir=$SUPERLU_INC_DIR --superlu-lib-dir=$SUPERLU_LIB_DIR --enable-openblas --openblas-inc-dir=$OPENBLAS_INC_DIR --openblas-lib-dir=$OPENBLAS_LIB_DIR configure"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd
