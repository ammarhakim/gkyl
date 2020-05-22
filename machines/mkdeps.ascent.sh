#.MF: It's possible this mkdeps (and/or the config) file for Ascent needs to be revised.
#.    Compilation fails at some medium to large size kernels, with error message:
#.    fatal error: error writing to /tmp/cc37dAA0.s: Cannot allocate memory.

module unload xl/16.1.1-5
module load gcc/8.1.1
module load spectrum-mpi/10.3.1.2-20200121
module load adios/1.13.1-py2
# If we are in machines directory, go up a directory.
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpiCC MPICC=mpicc MPICXX=mpiCC --build-luajit=yes --build-adios=no --build-eigen=yes --build-openmpi=no
