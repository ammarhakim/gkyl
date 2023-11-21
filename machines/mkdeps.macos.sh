# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft
export MACHINE_NAME='macos'
cd install-deps
# 1) Build OpenMPI.
./mkdeps.sh CC=clang CXX=clang++ --build-openmpi=yes --build-luajit=no  --prefix=$GKYLSOFT

# 2) Build gkylzero and adios.
./mkdeps.sh CC=clang CXX=clang++ MPICC=$GKYLSOFT/openmpi-4.0.5/bin/mpicc MPICXX=$GKYLSOFT/openmpi-4.0.5/bin/mpicxx --prefix=$GKYLSOFT --build-gkylzero=yes --build-adios=yes --build-luajit=no --build-openmpi=no

# 3) Build luajit. Some Macs (newer ones) have to comment out the next 2 lines.
# Get OSX version (XX.XX.X) and remove patch version so that only XX.X is
# assigned to MACOSX_DEPLOYMENT_TARGET.
vers=`sw_vers -productVersion`
export MACOSX_DEPLOYMENT_TARGET=${vers%.*}
./mkdeps.sh CC=clang CXX=clang++ MPICC=$GKYLSOFT/openmpi-4.0.5/bin/mpicc MPICXX=$GKYLSOFT/openmpi-4.0.5/bin/mpicxx --prefix=$GKYLSOFT --build-luajit=yes --build-openmpi=no 
