#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/eigen-3.3.7

# delete old checkout and builds
rm -rf eigen-3.3.7.tar* eigen-eigen-*

curl -L https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.tar.gz > eigen-3.3.7.tar.gz
gunzip eigen-3.3.7.tar.gz
tar -xvf eigen-3.3.7.tar

mkdir -p $PREFIX/include/eigen3
cp -R eigen-*/Eigen $PREFIX/include/eigen3
cp -R eigen-*/signature_of_eigen3_matrix_library $PREFIX/include/eigen3
cp -R eigen-*/unsupported $PREFIX/include/eigen3

rm -rf eigen-3.3.7.tar* eigen-eigen-*

# softlink to make finding easier
ln -sf $PREFIX $GKYLSOFT/eigen3

