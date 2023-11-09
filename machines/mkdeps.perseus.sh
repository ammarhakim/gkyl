# The following modules should always be loaded in bashrc (or equivalent) to build on perseus.

#module load intel/19.1/64/19.1.1.217
#module load openmpi/intel-17.0/1.10.2/64 

# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --prefix=$GKYLSOFT --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
