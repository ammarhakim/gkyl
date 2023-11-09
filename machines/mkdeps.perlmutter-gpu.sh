module load python/3.9-anaconda-2021.11
module load openmpi/5.0.0rc12
module load cudatoolkit/12.0
module load nccl/2.18.3-cu12
module unload darshan
CC=cc
CXX=CC
MPICC=mpicc
MPICXX=mpicxx
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/perlmutter/gkeyll/code/gpu/gkylsoft
export MACHINE_NAME='perlmutter-gpu'
cd install-deps
./mkdeps.sh CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --prefix=$GKYLSOFT --build-gkylzero=yes --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
