#!/bin/bash

source ./build-opts.sh

## This script builds CZMQ (C wrapper around ZMQ)

# Install prefix
PREFIX=$GKYLSOFT/czmq-4.0.2

# delete old checkout and builds
rm -rf czmq-4.0.2 czmq-4.0.2.tar*

curl -L https://github.com/zeromq/czmq/releases/download/v4.0.2/czmq-4.0.2.tar.gz > czmq-4.0.2.tar.gz
gunzip czmq-4.0.2.tar.gz
tar -xvf czmq-4.0.2.tar
cd czmq-4.0.2
export CFLAGS=-I$$GKYLSOFT/zeromq-4.2.2/include
export LDFLAGS=-L$GKYLSOFT/zeromq-4.2.2/lib
export PKG_CONFIG_PATH=$GKYLSOFT/zeromq-4.2.2/lib/pkgconfig
./configure --prefix=$PREFIX
make install

# softlink to make finding easier
ln -sf $PREFIX $GKYLSOFT/czmq
