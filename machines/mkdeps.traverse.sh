module load pgi/19.5/64
module load openmpi/pgi-19.5/3.1.4/64
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/gkylsoft'
cd install-deps
# first build LuaJIT with GCC
./mkdeps.sh CC=gcc --build-luajit-ppcle=yes --build-adios=no --build-eigen=no --build-openmpi=no
# build rest of the code
./mkdeps.sh CC=pgcc CXX=pgc++ MPICC=mpicc MPICXX=mpicxx --build-luajit-ppcle=no --build-adios=yes --build-eigen=yes --build-openmpi=no
