#!/bin/bash

# Edit the paths and options in the following command to suit your system
# The following modules should be loaded to configure on Traverse.
# module load cudatoolkit
# module load rh/devtoolset/8
# module load openmpi/cuda-10.2/devtoolset-8/4.0.3/64
# module load git/2.18

# Build directory
OUT=build
# Install location
GKYLSOFT=$HOME/gkylsoft
PREFIX=$GKYLSOFT/gkyl

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
MPI_INC_DIR=/usr/local/openmpi/cuda-10.2/4.0.3/devtoolset-8/ppc64le/include
MPI_LIB_DIR=/usr/local/openmpi/cuda-10.2/4.0.3/devtoolset-8/ppc64le/lib64
MPI_LINK_LIBS="mpi"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$GKYLSOFT/adios2/include
ADIOS_LIB_DIR=$GKYLSOFT/adios2/lib64
ADIOS_LINK_LIBS="adios2_c_mpi"

# EIGEN options
EIGEN_INC_DIR=$GKYLSOFT/eigen3/include/eigen3

# CUDA options
ENABLE_CUDA="--enable-cuda"
CUTOOLS_INC_DIR=$CPATH
CUTOOLS_LIB_DIR=$LIBRARY_PATH
CUTOOLS_LINK_LIBS="cudart,cublas"

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT -p $GKYLSOFT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --adios-link-libs=$ADIOS_LINK_LIBS $ENABLE_CUDA --cuda-inc-dir=$CUTOOLS_INC_DIR --cuda-lib-dir=$CUTOOLS_LIB_DIR --cuda-link-libs=$CUTOOLS_LINK_LIBS configure"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd
