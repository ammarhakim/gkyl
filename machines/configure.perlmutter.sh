#!/bin/bash

# Edit the paths and options in the following command to suit your system
module load PrgEnv-gnu/8.3.3
module load cray-mpich/8.1.22
module load python/3.9-anaconda-2021.11
module load cudatoolkit/11.7
module load nccl/2.14.3
module unload darshan

# Build directory
OUT=build
# Install location
GKYLSOFT=$HOME/gkylsoft
PREFIX=$GKYLSOFT/gkyl

# Compile flags (set optimization/debug flags here)
CC=cc
CXX=CC
MPICC=cc
MPICXX=CC
CXXFLAGS='-O3,-std=c++17'

# LuaJIT options
LUAJIT_INC_DIR=$HOME/gkylsoft/luajit/include/luajit-2.1
LUAJIT_LIB_DIR=$HOME/gkylsoft/luajit/lib
LUAJIT_SHARE_DIR=$HOME/gkylsoft/luajit/share/luajit-2.1.0-beta3

## MPI options
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=
MPI_LIB_DIR=
MPI_LINK_LIBS="mpi" #"mpich"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$GKYLSOFT/adios/include
ADIOS_LIB_DIR=$GKYLSOFT/adios/lib

# EIGEN options
EIGEN_INC_DIR=$GKYLSOFT/eigen3/include/eigen3

OPENBLAS_INC_DIR=$GKYLSOFT/OpenBLAS/include
OPENBLAS_LIB_DIR=$GKYLSOFT/OpenBLAS/lib

SUPERLU_INC_DIR=$GKYLSOFT/superlu/include
SUPERLU_LIB_DIR=$GKYLSOFT/superlu/lib

GKYLZERO_INC_DIR=$GKYLSOFT/gkylzero/include
GKYLZERO_LIB_DIR=$GKYLSOFT/gkylzero/lib

# CUDA options
ENABLE_CUDA="--enable-cuda"
CUTOOLS_INC_DIR=$CPATH
CUTOOLS_LIB_DIR=$LD_LIBRARY_PATH
CUTOOLS_LINK_LIBS="cudart"

# NCCL
ENABLE_NCCL="--enable-nccl"
NCCL_INC_DIR=/global/common/software/nersc/pm-2022q3/sw/nccl/2.14.3-1_cuda11.7/include
NCCL_LIB_DIR=/global/common/software/nersc/pm-2022q3/sw/nccl/2.14.3-1_cuda11.7/lib
NCCL_LINK_LIBS="nccl"

# You probably do not need to modify the command itself
cmd="./waf CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --out=$OUT -p $GKYLSOFT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR --eigen-inc-dir=$EIGEN_INC_DIR --enable-gkylzero --gkylzero-inc-dir=$GKYLZERO_INC_DIR --gkylzero-lib-dir=$GKYLZERO_LIB_DIR --enable-superlu --superlu-inc-dir=$SUPERLU_INC_DIR --superlu-lib-dir=$SUPERLU_LIB_DIR --enable-openblas --openblas-inc-dir=$OPENBLAS_INC_DIR --openblas-lib-dir=$OPENBLAS_LIB_DIR $ENABLE_CUDA --cuda-inc-dir=$CUTOOLS_INC_DIR --cuda-lib-dir=$CUTOOLS_LIB_DIR --cuda-link-libs=$CUTOOLS_LINK_LIBS $ENABLE_NCCL --nccl-inc-dir=$NCCL_INC_DIR --nccl-lib-dir=$NCCL_LIB_DIR --nccl-link-libs=$NCCL_LINK_LIBS configure"
# if we are in machines directory, go up a directory before executing cmd
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
echo $cmd
$cmd
