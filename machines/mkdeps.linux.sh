# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft
export MACHINE_NAME='linux'
cd install-deps
# first build OpenMPI
./mkdeps.sh CC=gcc CXX=g++ --build-openmpi=yes --prefix=$GKYLSOFT
# now build rest of packages
./mkdeps.sh CC=gcc CXX=g++ MPICC=$GKYLSOFT/openmpi-4.0.5/bin/mpicc MPICXX=$GKYLSOFT/openmpi-4.0.5/bin/mpicxx --build-gkylzero=yes --build-luajit=yes --build-adios=yes
