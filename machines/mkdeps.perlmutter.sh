module load PrgEnv-gnu/8.3.3
module load cray-mpich/8.1.22
module load python/3.9-anaconda-2021.11
module load cudatoolkit/11.7
module load nccl/2.15.5-ofi
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
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
