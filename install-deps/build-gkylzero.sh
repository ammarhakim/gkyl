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
    ./machines/mkdeps.frontera.sh
    ./machines/configure.frontera.sh
elif [ "$MACHINE_NAME" = 'stellar-intel' ]; then
    ./machines/mkdeps.stellar-intel.sh
    ./machines/configure.stellar-intel.sh  
elif [ "$MACHINE_NAME" = 'stellar-amd' ]; then
    ./machines/mkdeps.stellar-amd.sh
    ./machines/configure.stellar-amd.sh  
elif [ "$MACHINE_NAME" = 'perlmutter-cpu' ]; then
    ./machines/mkdeps.perlmutter-cpu.sh
    ./machines/configure.perlmutter-cpu.sh     
elif [ "$MACHINE_NAME" = 'perlmutter-gpu' ]; then
    ./machines/mkdeps.perlmutter-gpu.sh
    ./machines/configure.perlmutter-gpu.sh   
elif [ "$MACHINE_NAME" = 'della-gpu' ]; then
    ./machines/mkdeps.della-gpu.sh
    ./machines/configure.della-gpu.sh  
elif [ "$MACHINE_NAME" = 'macos' ]; then
    ./machines/mkdeps.macos.sh
    ./machines/configure.macos.sh  
elif [ "$MACHINE_NAME" = 'linux' ]; then
    ./machines/mkdeps.linux.sh
    ./machines/configure.linux.host.sh  
elif [ "$MACHINE_NAME" = 'linux-dev' ]; then
    ./machines/mkdeps.linux.sh
    ./machines/configure.linux.dev.sh  
fi
make -j
make install

