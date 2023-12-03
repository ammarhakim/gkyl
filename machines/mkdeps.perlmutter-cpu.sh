module load PrgEnv-gnu/8.3.3
module load cray-mpich/8.1.25
module unload darshan
CC=cc
CXX=CC
MPICC=cc
MPICXX=CC
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft
export MACHINE_NAME='perlmutter-cpu'
cd install-deps
./mkdeps.sh CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --prefix=$GKYLSOFT --build-gkylzero=yes --build-luajit=yes --build-adios=yes --build-openmpi=no
