#!/bin/bash

source ./build-opts.sh

## This script builds Gkylzero

# Install prefix
PREFIX=$GKYLSOFT/gkylzero

# delete old checkout and builds
rm -rf gkylzero

git clone https://github.com/ammarhakim/gkylzero.git
cd gkylzero
if [ "$HOSTNAME" = login1.frontera.tacc.utexas.edu ]; then
    ./machines/mkdeps.frontera.sh
    ./machines/configure.frontera.sh
else
    if [ $(uname) == 'Darwin' ]; then
        ./machines/mkdeps.macos.sh
        ./machines/configure.macos.sh
    else
        if [ $USE_GPU ]; then
            ./machines/mkdeps.linux.sh
            ./machines/configure.linux.dev.sh
        else 
            ./machines/mkdeps.linux.sh
            ./machines/configure.linux.host.sh
        fi
    fi
fi
make -j
make install

