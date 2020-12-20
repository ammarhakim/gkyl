#!/bin/bash

source ./build-opts.sh

## This script builds OpenRESTY branch of LuaJIT

# Install prefix
PREFIX=$GKYLSOFT/luajit-2.1.0-beta3-openresty

# delete old checkout and builds
rm -rf luajit2

git clone https://github.com/openresty/luajit2.git
cd luajit2
# IF NOT ON MACOSX, DOWNGRADE LUAJIT TO A VERSION THAT WORKS WITH ADIOS ON CLUSTERS.
# THIS VERSION IS FROM JUNE 15th, 2020 AND CAN BE FOUND HERE
# https://github.com/openresty/luajit2/commit/a44f53acf53603e7d9b88352de035b1804be4e88
# CURRENTLY DEBUGGING WHAT CHANGE TO LUAJIT IS LEADING TO NOT BEING ABLE TO USE ADIOS
if [ -z "$MACOSX_DEPLOYMENT_TARGET" ]
then 
   git checkout a44f53acf53603e7d9b88352de035b1804be4e88
fi
make PREFIX=$PREFIX CC=$CC CFLAGS=-fPIC
make XCFLAGS=-DLUAJIT_ENABLE_GC64 install PREFIX=$PREFIX

# softlink to make finding easier
ln -sf $PREFIX $GKYLSOFT/luajit

# soft-link executable name "lua". This allows running various tools
# (luarocks) needing lua executable to run
ln -sf $PREFIX/bin/luajit-2.1.0-beta3 $PREFIX/bin/lua

