#!/bin/sh
curl -L http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz > eigen-3.3.4.tar.gz
gunzip eigen-3.3.4.tar.gz
tar -xvf eigen-3.3.4.tar

cd eigen-eigen-*
mkdir build-dir; cd build-dir
cmake ../ -DCMAKE_INSTALL_PREFIX=$PREFIX
make install
