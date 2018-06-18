#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/eigen-3.3.4

# delete old checkout and builds
rm -rf eigen-3.3.4.tar* eigen-eigen-*

curl -L http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz > eigen-3.3.4.tar.gz
gunzip eigen-3.3.4.tar.gz
tar -xvf eigen-3.3.4.tar
mv eigen-eigen-* $PREFIX

# softlink to make finding easier
ln -sf $PREFIX $GKYLSOFT/eigen3

# build tests and docs (not necessary)
#cd $PREFIX
#mkdir build-dir; cd build-dir
#cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX
#make install

