#!/bin/bash
conda config --set anaconda_upload yes

cd openmpi
conda build .
cd ..

cd luajit
conda build .
cd ..

cd adios
conda build .
cd ..

cd adiospy
conda build .
cd ..

cd gkyl
conda build .
cd ..

conda build purge
