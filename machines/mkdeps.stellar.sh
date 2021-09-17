module load intel/2021.1.2
module load openmpi/intel-2021.1/4.1.0
module load fftw/intel-2021.1/openmpi-4.1.0/3.3.9
module load anaconda3/2020.11
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no --build-fftw=no
