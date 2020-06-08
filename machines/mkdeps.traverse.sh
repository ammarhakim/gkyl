# The following modules should be loaded to build on Traverse.
# module load cudatoolkit
# module load rh/devtoolset/8
# module load openmpi/devtoolset-8/4.0.3rc1/64
# module load git/2.18
# In addition, make sure /sbin/ is in your path to link luajit
# export PATH=$PATH:/sbin/

# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
