#!/bin/bash

# Edit the paths and options in the following command to suit your system
module load cray-mpich
module load adios
module load eigen
module load sz
module load zfp
module load zlib
CC=icc
CXX=icpc
MPICC=cc
MPICXX=CC

# Build directory
OUT=build
# Install location
PREFIX=$HOME/gkylsoft/gkyl

# Compile flags (set optimization/debug flags here)
CXXFLAGS='-O3,-std=c++11'

# LuaJIT options
LUAJIT_INC_DIR=$HOME/gkylsoft/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$HOME/gkylsoft/luajit/lib
LUAJIT_SHARE_DIR=$HOME/gkylsoft/luajit/share/luajit-2.1.0-beta3

# MPI options
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=$MPICH_DIR/include
MPI_LIB_DIR=$MPICH_DIR/lib
MPI_LINK_LIBS="mpich"

# ADIOS options
ENABLE_ADIOS="--enable-adios" 
ADIOS_INC_DIR=$ADIOS_INC
ADIOS_LIB_DIR="$ADIOS_LIB,$SZ_LIB,$ZFP_LIB,$ZLIB_LIBDIR"
ADIOS_LINK_LIBS="zfp,sz,z"

# EIGEN options
EIGEN_INC_DIR=$EIGEN_DIR/include/eigen3

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR  $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --adios-link-libs=$ADIOS_LINK_LIBS --eigen-inc-dir=$EIGEN_INC_DIR configure"
cd ..
echo $cmd
$cmd
