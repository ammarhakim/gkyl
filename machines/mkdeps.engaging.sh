module use /home/software/psfc/modulefiles
module add psfc/config
module load intel/2017-01
module load impi/2017-01
module load psfc/adios/1.13.1 
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=icc CXX=icpc MPICC=mpicc MPICXX=mpicxx --build-luajit=yes --build-adios=no --build-eigen=yes --build-openmpi=no
