#!/bin/bash

# Edit the paths and options in the following command to suit your system
module reset
#module load foss # gcc
module load intel/2019b # intel

# Build directory
OUT=build
# Install location
GKYLSOFT=$HOME/tinkercliffs/gkylsoft
PREFIX=$GKYLSOFT/gkyl

# Compile flags (set optimization/debug flags here) for GCC
#CC=gcc
#CXX=g++
#CXXFLAGS='-O3,-std=c++17,-mtune=znver2,-march=znver2,-mavx2'

# Compile flags (set optimization/debug flags here) for Intel
CC=icc
CXX=icpc
CXXFLAGS='-O3,-std=c++17,-ffreestanding,-march=core-avx2'

# LuaJIT options
LUAJIT_INC_DIR=$GKYLSOFT/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$GKYLSOFT/luajit/lib
LUAJIT_SHARE_DIR=$GKYLSOFT/luajit/share/luajit-2.1.0-beta3

## MPI options
MPICC=mpicc
MPICXX=mpicxx
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=/apps/easybuild/software/tinkercliffs-rome/impi/2018.5.288-iccifort-2019.5.281/bin64/
MPI_LIB_DIR=/apps/easybuild/software/tinkercliffs-rome/impi/2018.5.288-iccifort-2019.5.281/lib64/
#MPI_INC_DIR=/apps/easybuild/software/tinkercliffs-rome/OpenMPI/4.0.3-GCC-9.3.0/bin
#MPI_LIB_DIR=/apps/easybuild/software/tinkercliffs-rome/OpenMPI/4.0.3-GCC-9.3.0/lib
MPI_LINK_LIBS="mpi"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$GKYLSOFT/adios2/include
ADIOS_LIB_DIR=$GKYLSOFT/adios2/lib64
ADIOS_LINK_LIBS="adios2_c_mpi"

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT -p $GKYLSOFT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --adios-link-libs=$ADIOS_LINK_LIBS configure --disable-sqlite"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd

