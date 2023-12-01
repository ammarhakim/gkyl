module load gcc/9.3.0
module load openmpi/4.0.5
# if we are in machines directory, go up a directory before 
# executing commands in this script
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi

export GKYLSOFT=$HOME/gkylsoft
export MACHINE_NAME='engaging'
cd install-deps
./mkdeps.sh --prefix=$GKYLSOFT CC=gcc CXX=g++ MPICC=mpicc MPICXX=mpicxx --build-gkylzero=yes --build-luajit=yes --build-adios=yes --build-openmpi=no
