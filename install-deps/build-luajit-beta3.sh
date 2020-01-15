#!/bin/bash

source ./build-opts.sh

## If "git checkout v2.1' does not work just download the ZIP file for
## LuaJIT 2.1beta3 and install that instead.

# Install prefix
PREFIX=$GKYLSOFT/luajit-2.1.0-beta3

# delete old checkout and builds
rm -rf luajit-2.0

git clone http://luajit.org/git/luajit-2.0.git
cd luajit-2.0
git checkout v2.1
make PREFIX=$PREFIX CC=$CC
make XCFLAGS=-DLUAJIT_ENABLE_GC64 install PREFIX=$PREFIX

