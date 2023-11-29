# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft/
export MACHINE_NAME='puma'
cd install-deps
# now build rest of packages
./mkdeps.sh CC=icc CXX=icpc MPICC=mpicc MPICXX=mpicxx --build-gkylzero=yes --build-luajit=yes --build-adios=yes
