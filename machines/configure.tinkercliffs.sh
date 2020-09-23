#!/bin/bash

# Edit the paths and options in the following command to suit your system
module reset
module load foss
#module load gcc/8.2.0
#module load openmpi/gcc/64/4.0.4

# Build directory
OUT=build
# Install location
PREFIX=$HOME/tinkercliffs/gkylsoft/gkyl

# Compile flags (set optimization/debug flags here)
CC=gcc
CXX=g++
CXXFLAGS='-O3,-std=c++17'

# LuaJIT options
LUAJIT_INC_DIR=$HOME/tinkercliffs/gkylsoft/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$HOME/tinkercliffs/gkylsoft/luajit/lib
LUAJIT_SHARE_DIR=$HOME/tinkercliffs/gkylsoft/luajit/share/luajit-2.1.0-beta3

## MPI options
MPICC=mpicc
MPICXX=mpicxx
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=/apps/easybuild/software/tinkercliffs-rome/OpenMPI/4.0.3-GCC-9.3.0/bin
MPI_LIB_DIR=/apps/easybuild/software/tinkercliffs-rome/OpenMPI/4.0.3-GCC-9.3.0/lib
MPI_LINK_LIBS="mpi"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$HOME/tinkercliffs/gkylsoft/adios/include
ADIOS_LIB_DIR=$HOME/tinkercliffs/gkylsoft/adios/lib

# EIGEN options
EIGEN_INC_DIR=$HOME/tinkercliffs/gkylsoft/eigen3/include/eigen3

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --eigen-inc-dir=$EIGEN_INC_DIR configure --disable-sqlite"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd

