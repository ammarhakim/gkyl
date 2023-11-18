# The following modules should be loaded to build on Adroit.
# module load intel
# module load openmpi/cuda-10.2/intel-19.1/4.0.3/64
# module load cudatoolkit/10.2
# module load rh/devtoolset/8
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --build-luajit=yes --build-adios=yes --build-openmpi=no
