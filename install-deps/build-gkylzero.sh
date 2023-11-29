#!/bin/bash

source ./build-opts.sh

## This script builds Gkylzero

# Install prefix
PREFIX=$GKYLSOFT/gkylzero

# delete old checkout and builds
rm -rf gkylzero

git clone https://github.com/ammarhakim/gkylzero.git
cd gkylzero
if [ "$MACHINE_NAME" = 'frontera' ]; then
    PREFIX=$GKYLSOFT ./machines/mkdeps.frontera.sh
    PREFIX=$GKYLSOFT ./machines/configure.frontera.sh
elif [ "$MACHINE_NAME" = 'stellar-intel' ]; then
    PREFIX=$GKYLSOFT ./machines/mkdeps.stellar-intel.sh
    PREFIX=$GKYLSOFT ./machines/configure.stellar-intel.sh  
elif [ "$MACHINE_NAME" = 'stellar-amd' ]; then
    PREFIX=$GKYLSOFT ./machines/mkdeps.stellar-amd.sh
    PREFIX=$GKYLSOFT ./machines/configure.stellar-amd.sh  
elif [ "$MACHINE_NAME" = 'perlmutter-cpu' ]; then
    PREFIX=$GKYLSOFT ./machines/mkdeps.perlmutter.cpu.sh
    PREFIX=$GKYLSOFT ./machines/configure.perlmutter.cpu.sh     
elif [ "$MACHINE_NAME" = 'perlmutter-gpu' ]; then
    PREFIX=$GKYLSOFT ./machines/mkdeps.perlmutter.gpu.sh
    PREFIX=$GKYLSOFT ./machines/configure.perlmutter.gpu.sh   
elif [ "$MACHINE_NAME" = 'della-gpu' ]; then
    PREFIX=$GKYLSOFT ./machines/mkdeps.della-gpu.sh
    PREFIX=$GKYLSOFT ./machines/configure.della-gpu.sh  
elif [ "$MACHINE_NAME" = 'macos' ]; then
    PREFIX=$GKYLSOFT ./machines/mkdeps.macos.sh
    PREFIX=$GKYLSOFT ./machines/configure.macos.sh  
elif [ "$MACHINE_NAME" = 'linux' ]; then
    PREFIX=$GKYLSOFT ./machines/mkdeps.linux.sh
    PREFIX=$GKYLSOFT ./machines/configure.linux.cpu.sh  
elif [ "$MACHINE_NAME" = 'linux-gpu' ]; then
    PREFIX=$GKYLSOFT ./machines/mkdeps.linux.sh
    PREFIX=$GKYLSOFT ./machines/configure.linux.gpu.sh  
else
    echo 'Machine option not implemented in build-gkylzero.sh'
    exit 1
fi
make -j16
make install

