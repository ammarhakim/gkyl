#!/bin/bash

# Edit the paths and options in the following command to suit your system
module load PrgEnv-gnu/8.3.3
module load nvhpc-mixed/21.11
#module load cray-mpich/8.1.15
module load python/3.9-anaconda-2021.11
module load craype-x86-milan

# Build directory
OUT=build
# Install location
PREFIX=$HOME/gkylsoft/gkyl

# Compile flags (set optimization/debug flags here)
CC=cc
CXX=CC
MPICC=cc
MPICXX=CC
CXXFLAGS='-O3,-std=c++17'

# LuaJIT options
LUAJIT_INC_DIR=$HOME/gkylsoft/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$HOME/gkylsoft/luajit/lib
LUAJIT_SHARE_DIR=$HOME/gkylsoft/luajit/share/luajit-2.1.0-beta3

## MPI options
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR= #$MPI_DIR/include
MPI_LIB_DIR= #$MPI_DIR/lib
MPI_LINK_LIBS="mpi" #"mpich"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$HOME/gkylsoft/adios/include
ADIOS_LIB_DIR=$HOME/gkylsoft/adios/lib

# EIGEN options
EIGEN_INC_DIR=$HOME/gkylsoft/eigen3/include/eigen3

OPENBLAS_INC_DIR=$HOME/gkylsoft/OpenBLAS/include
OPENBLAS_LIB_DIR=$HOME/gkylsoft/OpenBLAS/lib

SUPERLU_INC_DIR=$HOME/gkylsoft/superlu/include
SUPERLU_LIB_DIR=$HOME/gkylsoft/superlu/lib

GKYLZERO_INC_DIR=$HOME/gkylsoft/gkylzero/include
GKYLZERO_LIB_DIR=$HOME/gkylsoft/gkylzero/lib

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --enable-gkylzero --gkylzero-inc-dir=$GKYLZERO_INC_DIR --gkylzero-lib-dir=$GKYLZERO_LIB_DIR --enable-superlu --superlu-inc-dir=$SUPERLU_INC_DIR --superlu-lib-dir=$SUPERLU_LIB_DIR --enable-openblas --openblas-inc-dir=$OPENBLAS_INC_DIR --openblas-lib-dir=$OPENBLAS_LIB_DIR configure"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd