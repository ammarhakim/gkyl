module purge
module load intel/17.0
module load openmpi/3.1.2
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/cascades/gkylsoft'
cd install-deps
./mkdeps.sh CC=icc CXX=icpc MPICC=mpicc MPICXX=mpicxx --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
