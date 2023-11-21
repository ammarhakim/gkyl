# The following modules should be loaded to build on Traverse.
# module load cudatoolkit
# module load rh/devtoolset/8
# module load openmpi/cuda-10.2/devtoolset-8/4.0.3/64
# module load git/2.18
# In addition, make sure /sbin/ is in your path to link luajit
# export PATH=$PATH:/sbin/

# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft
export MACHINE_NAME='traverse'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --prefix=$GKYLSOFT --build-gkylzero=yes --build-luajit=yes --build-adios=yes --build-openmpi=no
