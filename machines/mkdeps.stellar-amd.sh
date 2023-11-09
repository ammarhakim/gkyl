module load gcc/8
module load cudatoolkit/12.0
module load openmpi/cuda-11.1/gcc/4.1.1
module load anaconda3/2020.11
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft-amd
export MACHINE_NAME='stellar-amd'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --prefix=$GKYLSOFT --build-gkylzero=yes --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
