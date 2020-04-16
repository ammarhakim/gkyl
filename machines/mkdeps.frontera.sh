module swap intel/19.0.5 intel/18.0.5
module swap impi/19.0.5 impi/18.0.5
export PATH=$PATH:/usr/sbin/
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
