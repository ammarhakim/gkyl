#!/bin/bash

source ./build-opts.sh

# Install prefix
PREFIX=$GKYLSOFT/spectra

# delete old checkout and builds
rm -rf spectra.tar* spectra-master.tar*

curl -L https://github.com/yixuan/spectra/tarball/master > spectra-master.tar.gz
gunzip spectra-master.tar.gz
tar -xvf spectra-master.tar

mv yixuan-spectra-* $PREFIX

rm -rf spectra-master.tar* 

# softlink to make finding easier
ln -sf $PREFIX $GKYLSOFT/spectra

