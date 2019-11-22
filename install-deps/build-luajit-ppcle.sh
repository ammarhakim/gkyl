#!/bin/bash

source ./build-opts.sh

## This script builds OpenRESTY branch of LuaJIT

# Install prefix
PREFIX=$GKYLSOFT/luajit-2.1.0-ppcle

# delete old checkout and builds
rm -rf luajit2

git clone https://github.com/siddhesh/LuaJIT.git
cd LuaJIT
make PREFIX=$PREFIX CC=$CC
make XCFLAGS=-DLUAJIT_ENABLE_GC64 install PREFIX=$PREFIX

# softlink to make finding easier
ln -sf $PREFIX $GKYLSOFT/luajit

# soft-link executable name "lua". This allows running various tools
# (luarocks) needing lua executable to run
ln -sf $PREFIX/bin/luajit-2.1.0-ppcle $PREFIX/bin/lua
