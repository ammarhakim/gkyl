# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft
export MACHINE_NAME='macos'
# get OSX version (XX.XX.X)
vers=`sw_vers -productVersion`
# remove patch version, so that only XX.XX
export MACOSX_DEPLOYMENT_TARGET=${vers%.*}
cd install-deps
# first build OpenMPI
./mkdeps.sh CC=clang CXX=clang++ --build-openmpi=yes --prefix=$GKYLSOFT
# now build rest of packages
./mkdeps.sh CC=clang CXX=clang++ MPICC=$GKYLSOFT/openmpi-4.0.5/bin/mpicc MPICXX=$GKYLSOFT/openmpi-4.0.5/bin/mpicxx --prefix=$GKYLSOFT --build-gkylzero=yes --build-luajit=yes --build-adios=yes
