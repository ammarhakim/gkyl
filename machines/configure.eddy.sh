#!/bin/bash

# Edit the paths and options in the following command to suit your system
module load intel
module load intel-mpi

# Build directory
OUT=build-par
# Install location
PREFIX=$HOME/gkylsoft/gkyl

# Compile flags (set optimization/debug flags here)
CC=icc
CXX=icpc
CXXFLAGS='-O3,-std=c++14,-axMIC-AVX512'

# LuaJIT options
LUAJIT_INC_DIR=$HOME/gkylsoft/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$HOME/gkylsoft/luajit/lib
LUAJIT_SHARE_DIR=$HOME/gkylsoft/luajit/share/luajit-2.1.0-beta3

## MPI options
MPICC=mpicc
MPICXX=mpicxx
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=$TACC_IMPI_INC
MPI_LIB_DIR=$TACC_IMPI_LIB
MPI_LINK_LIBS="mpi,mpicxx"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$HOME/gkylsoft/adios/include
ADIOS_LIB_DIR=$HOME/gkylsoft/adios/lib

# EIGEN options
EIGEN_INC_DIR=$HOME/gkylsoft/eigen/include/eigen3

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR configure"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd
