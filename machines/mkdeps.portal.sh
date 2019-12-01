module load gcc
module load openmpi
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
