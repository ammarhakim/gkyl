# The following modules should be loaded to build on Saga.
# module load fosscuda/2020b
# module load Eigen/3.3.8-GCCcore-10.2.0
# module load Python/3.8.6-GCCcore-10.2.0

# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
if [ -z "$GKYLSOFT" ]
  then
    export GKYLSOFT=$(readlink -f ../gkylsoft)
fi
cd install-deps
./mkdeps.sh --prefix=$GKYLSOFT CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --build-luajit=yes --build-adios=yes --build-eigen=no --build-openmpi=no
