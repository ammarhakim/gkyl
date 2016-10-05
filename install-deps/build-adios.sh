#!/bin/bash

# Edit to suite your system
PREFIX=$HOME/gkylsoft/adios
# You may have to explicitly set the compilers to use, in which case,
# uncomment the following lines
#CC=gcc
#CXX=g++
#MPICC=mpicc
#MPICXX=mpicxx

wget http://users.nccs.gov/~pnorbert/adios-1.10.0.tar.gz
gunzip adios-1.10.0.tar.gz
tar -xvf adios-1.10.0.tar
cd adios-1.10.0
./configure --prefix=$PREFIX --disable-fortran
make install
