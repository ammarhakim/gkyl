#!/bin/bash

# Edit the paths and options in the following command to suit your system
# The following modules should be loaded to configure on Saga.
# module load fosscuda/2020b
# module load Eigen/3.3.8-GCCcore-10.2.0
# module load Python/3.8.6-GCCcore-10.2.0

# Build directory
OUT=build
# Install location
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
if [ -z "$GKYLSOFT" ]
  then
    export GKYLSOFT=$(readlink -f ../gkylsoft)
fi
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
MPI_INC_DIR=/cluster/software/OpenMPI/4.0.5-gcccuda-2020b/include
MPI_LIB_DIR=/cluster/software/OpenMPI/4.0.5-gcccuda-2020b/lib
MPI_LINK_LIBS="mpi"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$GKYLSOFT/adios2/include
ADIOS_LIB_DIR=$GKYLSOFT/adios2/lib64
ADIOS_LINK_LIBS="adios2_c_mpi"

# CUDA options
CUTOOLS_INC_DIR=$CPATH
CUTOOLS_LIB_DIR=$LIBRARY_PATH
CUTOOLS_LINK_LIBS="cudart,cublas"

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX -p $GKYLSOFT --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --adios-link-libs=$ADIOS_LINK_LIBS --cuda-inc-dir=$CUTOOLS_INC_DIR --cuda-lib-dir=$CUTOOLS_LIB_DIR --cuda-link-libs=$CUTOOLS_LINK_LIBS --disable-cuda configure"
echo $cmd
$cmd
