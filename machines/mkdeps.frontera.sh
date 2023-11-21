export PATH=$PATH:/usr/sbin/
# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft
export MACHINE_NAME='frontera'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --prefix=$GKYLSOFT --build-gkylzero=yes --build-gkylzero=yes --build-luajit=yes --build-adios=yes --build-openmpi=no
