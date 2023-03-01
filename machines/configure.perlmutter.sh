#!/bin/bash

# Edit the paths and options in the following command to suit your system
module load cray-mpich
module unload adios
module unload zlib
module unload darshan
CC=cc
CXX=CC
MPICC=cc
MPICXX=CC

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

# Build directory
OUT=build

# Compile flags (set optimization/debug flags here)
CXXFLAGS='-O3,-std=c++17,-fPIC'

# MPI options
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=$MPICH_DIR/include
MPI_LIB_DIR=$MPICH_DIR/lib
MPI_LINK_LIBS="mpich"

# ADIOS options
ENABLE_ADIOS="--enable-adios" 
ADIOS_INC_DIR=$GKYLSOFT/adios/include
ADIOS_LIB_DIR=$GKYLSOFT/adios/lib

# EIGEN options
EIGEN_INC_DIR=$GKYLSOFT/eigen3/include/eigen3

# LuaJIT options
LUAJIT_INC_DIR=$GKYLSOFT/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$GKYLSOFT/luajit/lib
LUAJIT_SHARE_DIR=$GKYLSOFT/luajit/share/luajit-2.1.0-beta3

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR  $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --eigen-inc-dir=$EIGEN_INC_DIR configure"
echo $cmd
$cmd
