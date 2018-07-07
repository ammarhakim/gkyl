#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/eigen-3.3.4

# delete old checkout and builds
rm -rf eigen-3.3.4.tar* eigen-eigen-*

curl -L http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz > eigen-3.3.4.tar.gz
gunzip eigen-3.3.4.tar.gz
tar -xvf eigen-3.3.4.tar

mkdir -p $PREFIX/include/eigen3
cp -R eigen-eigen*/Eigen $PREFIX/include/eigen3
cp -R eigen-eigen*/signature_of_eigen3_matrix_library $PREFIX/include/eigen3
cp -R eigen-eigen*/unsupported $PREFIX/include/eigen3

rm -rf eigen-3.3.4.tar* eigen-eigen-*

# softlink to make finding easier
ln -sf $PREFIX $GKYLSOFT/eigen3

