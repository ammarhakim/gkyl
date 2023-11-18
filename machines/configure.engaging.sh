#!/bin/bash

# Edit the paths and options in the following command to suit your system
module load gcc/9.3.0
module load openmpi/4.0.5

# Build directory
OUT=build
# if we are in machines directory, go up a directory before 
# executing commands in this script
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi

# Install location... by default, puts gkylsoft directory
# on same level as gkyl directory (where we are now, presumably)
export GKYLSOFT=$(readlink -f ../gkylsoft)
PREFIX=$GKYLSOFT/gkyl

# Compile flags (set optimization/debug flags here)
CC=gcc
CXX=g++
CXXFLAGS='-O3,-std=c++17'

# LuaJIT options
LUAJIT_INC_DIR=$GKYLSOFT/luajit/include/luajit-2.1/
LUAJIT_LIB_DIR=$GKYLSOFT/luajit/lib/
LUAJIT_SHARE_DIR=$GKYLSOFT/luajit/share/luajit-2.1.0-beta3/

## MPI options
MPICC=mpicc
MPICXX=mpicxx
ENABLE_MPI="--enable-mpi"
# no need to specify these paths, loading the openmpi module 
# ensures that everything can be found
MPI_INC_DIR=
MPI_LIB_DIR=
MPI_LINK_LIBS="mpi"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$GKYLSOFT/adios/include/
ADIOS_LIB_DIR=$GKYLSOFT/adios/lib/

# CUDA options
ENABLE_CUDA="--disable-cuda" # disable CUDA, even if nvcc is found

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR $ENABLE_CUDA configure"
echo $cmd
$cmd
