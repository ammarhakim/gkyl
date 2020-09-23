# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/gkylsoft'
cd install-deps
# first build OpenMPI
./mkdeps.sh CC=gcc CXX=g++ --build-openmpi=yes
# now build rest of packages
./mkdeps.sh CC=gcc CXX=g++ MPICC=$GKYLSOFT/openmpi-3.1.2/bin/mpicc MPICXX=$GKYLSOFT/openmpi-3.1.2/bin/mpicxx --build-luajit=yes --build-adios=yes --build-eigen=yes