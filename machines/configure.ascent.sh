#!/bin/bash

#.MF: It's possible this config (and/or the mkdeps) file for Ascent needs to be revised.
#.    Compilation fails at some medium to large size kernels, with error message:
#.      fatal error: error writing to /tmp/cc37dAA0.s: Cannot allocate memory.
#.    To overcome this problem, compile with ./waf build install -j1
#.    However, this eventually leads the to the following errors:
#.      [578/602] Compiling DataStruct/CartFieldDeviceImpl.cu
#.      /autofs/nccsopen-svm1_sw/ascent/gcc/8.1.1/include/c++/8.1.1/type_traits(347): error: identifier "__ieee128" is undefined
#.
#.      /autofs/nccsopen-svm1_sw/ascent/gcc/8.1.1/include/c++/8.1.1/bits/std_abs.h(101): error: identifier "__ieee128" is undefined
#.
#.      /autofs/nccsopen-svm1_sw/ascent/gcc/8.1.1/include/c++/8.1.1/bits/std_abs.h(102): error: identifier "__ieee128" is undefined
#.
#.      3 errors detected in the compilation of "/tmp/tmpxft_000043b2_00000000-6_CartFieldDeviceImpl.cpp1.ii".
#.
#.      Waf: Leaving directory `/autofs/nccsopen-svm1_home/manaf/gkeyll/gkyl/build'
#.      Build failed
#.       -> task in 'datastruct_cu' failed with exit status 1 (run with -v to display more information)

# Edit the paths and options in the following command to suit your system
module unload xl/16.1.1-5
module load gcc/8.1.1
module load spectrum-mpi/10.3.1.2-20200121
module load adios/1.13.1-py2
module load cuda/10.1.243

# Build directory
OUT=build
# Install location
PREFIX=$HOME/gkylsoft/gkyl

# Compile flags (set optimization/debug flags here)
CC=mpicc
CXX=mpiCC
CXXFLAGS='-O3,-std=c++17'

# LuaJIT options
LUAJIT_INC_DIR=$HOME/gkylsoft/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$HOME/gkylsoft/luajit/lib
LUAJIT_SHARE_DIR=$HOME/gkylsoft/luajit/share/luajit-2.1.0-beta3

## MPI options
MPICC=mpicc
MPICXX=mpiCC
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=$MPI_ROOT/include64
MPI_LIB_DIR=$MPI_ROOT/lib64
MPI_LINK_LIBS="mpi_ibm"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$OLCF_ADIOS_ROOT/include
ADIOS_LIB_DIR=$OLCF_ADIOS_ROOT/lib

# EIGEN options
EIGEN_INC_DIR=$HOME/gkylsoft/eigen3/include/eigen3

# CUDA options
CUTOOLS_INC_DIR=$CUDA_DIR/include/
CUTOOLS_LIB_DIR=$CUDA_DIR/lib64/
CUTOOLS_LINK_LIBS=cudart

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --cuda-inc-dir=$CUTOOLS_INC_DIR --cuda-lib-dir=$CUTOOLS_LIB_DIR configure"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd
