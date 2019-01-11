#!/bin/bash

# Edit the paths and options in the following command to suit your system
module load cray-mpich
module load adios
module load eigen3
module load zlib
CC=icc
CXX=icpc
MPICC=cc
MPICXX=CC

# Build directory
OUT=build
# Install location
GKYLSOFT='/global/common/software/m2116/nmandell/gkylsoft'
PREFIX=$GKYLSOFT/gkyl

# Compile flags (set optimization/debug flags here)
CXXFLAGS='-O3,-std=c++14,-fast'

# LuaJIT options
LUAJIT_INC_DIR=$GKYLSOFT/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$GKYLSOFT/luajit/lib
LUAJIT_SHARE_DIR=$GKYLSOFT/luajit/share/luajit-2.1.0-beta3

# MPI options
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=$MPICH_DIR/include
MPI_LIB_DIR=$MPICH_DIR/lib
MPI_LINK_LIBS="mpich"

# ADIOS options
ENABLE_ADIOS="--enable-adios" 
ADIOS_INC_DIR=$ADIOS_DIR/include
ADIOS_LIB_DIR="$ADIOS_DIR/lib,$ZLIB_DIR,/usr/common/software/minixml/2.9/hsw//intel/lib"
ADIOS_LINK_LIBS="z,mxml"

# EIGEN options
EIGEN_INC_DIR=$EIGEN3_DIR/include/eigen3

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR  $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --adios-link-libs=$ADIOS_LINK_LIBS --eigen-inc-dir=$EIGEN_INC_DIR configure"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd
