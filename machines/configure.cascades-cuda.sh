#!/bin/bash

# Edit the paths and options in the following command to suit your system
module purge
module load gcc/7.3.0
module load openmpi/3.0.0
module load cuda/10.1.168 

# Build directory
OUT=build
# Install location
PREFIX=$HOME/cascades/gkylsoft/gkyl-cuda

# Compile flags (set optimization/debug flags here)
CC=gcc
CXX=g++
CXXFLAGS='-O3,-std=c++17'

# LuaJIT options
LUAJIT_INC_DIR=$HOME/cascades/gkylsoft/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$HOME/cascades/gkylsoft/luajit/lib
LUAJIT_SHARE_DIR=$HOME/cascades/gkylsoft/luajit/share/luajit-2.1.0-beta3

## MPI options
MPICC=mpicc
MPICXX=mpicxx
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=$OMPI_INC
MPI_LIB_DIR=$OMPI_LIB
MPI_LINK_LIBS="mpi"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$HOME/cascades/gkylsoft/adios/include
ADIOS_LIB_DIR=$HOME/cascades/gkylsoft/adios/lib

# EIGEN options
EIGEN_INC_DIR=$HOME/cascades/gkylsoft/eigen3/include/eigen3

# CUDA options
CUTOOLS_INC_DIR=$CUDA_INC
CUTOOLS_LIB_DIR=$CUDA_LIB
CUTOOLS_LINK_LIBS=cudart

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --eigen-inc-dir=$EIGEN_INC_DIR --cuda-inc-dir=$CUTOOLS_INC_DIR --cuda-lib-dir=$CUTOOLS_LIB_DIR configure --disable-sqlite"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd

