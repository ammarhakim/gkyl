#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/adios2-2.9.1
PYTHON_DIR=$HOME/anaconda3/envs/pgkylEnv/

# delete old checkout and builds
rm -rf adios2-2.9.1.tar* adios2-2.9.1

curl -L https://github.com/ornladios/ADIOS2/archive/refs/tags/v2.9.1.tar.gz > adios2-2.9.1.tar.gz
gunzip adios2-2.9.1.tar.gz
mkdir adios2-2.9.1
tar -xvf adios2-2.9.1.tar -C adios2-2.9.1 --strip-components 1

cd adios2-2.9.1
mkdir build
cd build

CC=$MPICC CXX=$MPICXX cmake ../ -DCMAKE_C_FLAGS="-g -O3 -fPIC" -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX -DBUILD_TESTING=OFF -DADIOS2_BUILD_EXAMPLES=OFF -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_HDF5=OFF -DADIOS2_USE_Python=ON -DPYTHON_LIBRARY=$PYTHON_DIR/lib/python3.11/ -DPYTHON_INCLUDE_DIR=$PYTHON_DIR/include/python3.11/ -DPYTHON_EXECUTABLE=$PYTHON_DIR/bin/python -DPython_NumPy_INCLUDE_DIR=$PYTHON_DIR/lib/python3.11/site-packages/numpy/ -DPythonModule_mpi4py_PATH=$PYTHON_DIR/lib/python3.11/site-packages/mpi4py/

make -j 32 VERBOSE=1
make install

# soft-link
ln -sfn $PREFIX $GKYLSOFT/adios2
