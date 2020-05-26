#!/bin/bash

# Edit the paths and options in the following command to suit your system
#module load gcc/7.3.0
#module load cuda
module load gcc
module load openmpi
module unload gcc
export PATH=~/gkylsoft/gcc8/bin:/usr/local/cuda-10.2/bin:/usr/local/cuda-10.2/NsightCompute-2019.1:${PATH}
export LD_LIBRARY_PATH=~/gkylsoft/gcc8/lib64:/usr/local/cuda-10.2/lib64:${LD_LIBRARY_PATH}

# Build directory
OUT=build
# Install location
PREFIX=$HOME/gkylsoft/gkyl

# Compile flags (set optimization/debug flags here)
CC=gcc
CXX=g++
CXXFLAGS='-O3,-std=c++17'

# LuaJIT options
LUAJIT_INC_DIR=$HOME/gkylsoft/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$HOME/gkylsoft/luajit/lib
LUAJIT_SHARE_DIR=$HOME/gkylsoft/luajit/share/luajit-2.1.0-beta3

## MPI options
MPICC=mpicc
MPICXX=mpicxx
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=/usr/pppl/gcc/9.3-pkgs/openmpi-4.0.3/include
MPI_LIB_DIR=/usr/pppl/gcc/9.3-pkgs/openmpi-4.0.3/lib
#MPI_INC_DIR=/usr/pppl/gcc/7.3-pkgs/openmpi-3.0.0/include
#MPI_LIB_DIR=/usr/pppl/gcc/7.3-pkgs/openmpi-3.0.0/lib
MPI_LINK_LIBS="mpi"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$HOME/gkylsoft/adios/include
ADIOS_LIB_DIR=$HOME/gkylsoft/adios/lib

# EIGEN options
EIGEN_INC_DIR=$HOME/gkylsoft/eigen3/include/eigen3

# CUDA options
CUTOOLS_INC_DIR=/usr/local/cuda-10.2/include
CUTOOLS_LIB_DIR=/usr/local/cuda-10.2/lib64
#CUTOOLS_INC_DIR=/usr/pppl/cuda/10.1.243/include
#CUTOOLS_LIB_DIR=/usr/pppl/cuda/10.1.243/lib64
CUTOOLS_LINK_LIBS=cudart

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --cuda-inc-dir=$CUTOOLS_INC_DIR --cuda-lib-dir=$CUTOOLS_LIB_DIR --cuda-link-libs=$CUTOOLS_LINK_LIBS configure"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd
