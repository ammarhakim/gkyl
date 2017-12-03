#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/luarocks-2.4.3

# LuaJIT install information
LUA_BIN=$GKYLSOFT/luajit/bin
LUA_INCLUDE=$GKYLSOFT/luajit/include/luajit-2.1
LUA_LIB=$GKYLSOFT/luajit/lib

# delete old checkout and builds
rm -rf luarocks-2.4.3 luarocks-2.4.3.tar*

curl -L https://luarocks.org/releases/luarocks-2.4.3.tar.gz > luarocks-2.4.3.tar.gz
tar zxpf luarocks-2.4.3.tar.gz
cd luarocks-2.4.3
cmd="./configure --prefix=$PREFIX --with-lua-bin=$LUA_BIN --with-lua-include=$LUA_INCLUDE --with-lua-lib=$LUA_LIB "
echo $cmd
$cmd
make bootstrap

# softlink to make finding easier
ln -sf $PREFIX $GKYLSOFT/luarocks
