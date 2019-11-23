#!/bin/bash

source ./build-opts.sh

## This script builds ZeroMQ

# Install prefix
PREFIX=$GKYLSOFT/zeromq-4.2.2

# delete old checkout and builds
rm -rf zeromq-4.2.2 zeromq-4.2.2.tar*

curl -L https://github.com/zeromq/libzmq/releases/download/v4.2.2/zeromq-4.2.2.tar.gz > zeromq-4.2.2.tar.gz
gunzip zeromq-4.2.2.tar.gz
tar -xvf zeromq-4.2.2.tar
cd zeromq-4.2.2
./configure --prefix=$PREFIX
make install

# softlink to make finding easier
ln -sf $PREFIX $GKYLSOFT/zeromq

