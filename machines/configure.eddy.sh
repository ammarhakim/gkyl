#!/bin/bash

# Edit the paths and options in the following command to suit your system
module load intel
module unload intel-mpi
module load openmpi/intel-18.0/3.0.0/64 

# Build directory
OUT=build
# Install location
PREFIX=$HOME/gkylsoft/gkyl

# Compile flags (set optimization/debug flags here)
CC=icc
CXX=icpc
CXXFLAGS='-O3,-std=c++17'

# LuaJIT options
LUAJIT_INC_DIR=$HOME/gkylsoft/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$HOME/gkylsoft/luajit/lib
LUAJIT_SHARE_DIR=$HOME/gkylsoft/luajit/share/luajit-2.1.0-beta3

## MPI options
MPICC=mpicc
MPICXX=mpicxx
ENABLE_MPI="--enable-mpi"
I_MPI_ROOT=/usr/local/openmpi/3.0.0/intel180/x86_64
MPI_INC_DIR=$I_MPI_ROOT/include
MPI_LIB_DIR=$I_MPI_ROOT/lib64
MPI_LINK_LIBS="mpi"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$HOME/gkylsoft/adios/include
ADIOS_LIB_DIR=$HOME/gkylsoft/adios/lib

# EIGEN options
EIGEN_INC_DIR=$HOME/gkylsoft/eigen3/include/eigen3

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR configure"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd
