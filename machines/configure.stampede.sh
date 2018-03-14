#!/bin/bash

# Edit the paths and options in the following command to suit your system

# Build directory
OUT=build-par
# Install location
PREFIX=$HOME/gkylsoft/gkyl

# Compile flags (set optimization/debug flags here)
export CC='icc'
export CXX='icpc'
CXXFLAGS='-O3,-xCORE-AVX2,-axMIC-AVX512,-qopt-zmm-usage=high'

# LuaJIT options
LUAJIT_INC_DIR=$HOME/gkylsoft/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$HOME/gkylsoft/luajit/lib
LUAJIT_SHARE_DIR=$HOME/gkylsoft/luajit/share/luajit-2.1.0-beta3

## MPI options
export MPICC='mpicc'
export MPICXX='mpicxx'
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=/opt/intel/compilers_and_libraries_2018.0.128/linux/mpi/intel64/include
MPI_LIB_DIR=/opt/intel/compilers_and_libraries_2018.0.128/linux/mpi/intel64/lib
MPI_LINK_LIBS="mpi,mpicxx"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$HOME/gkylsoft/adios-1.13.0/include
ADIOS_LIB_DIR=$HOME/gkylsoft/adios-1.13.0/lib

# You probably do not need to modify the command itself
cmd="./waf --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR configure"
echo $cmd
$cmd
