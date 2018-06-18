module swap intel/18.0.0
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --build-luajit=yes --build-adios=yes --build-eigen=yes --build-openmpi=no
