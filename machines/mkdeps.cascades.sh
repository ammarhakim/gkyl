module purge
module load gcc/7.3.0
module load openmpi/3.0.0
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft
export MACHINE_NAME='cascades'
cd install-deps
./mkdeps.sh CC=gcc CXX=g++ MPICC=mpicc MPICXX=mpicxx --prefix=$GKYLSOFT --build-gkylzero=yes --build-luajit=yes --build-adios=yes --build-openmpi=no --build-luarocks=yes
