#.Install file for PPPL's Portal cluster.
#.As of May 30 2020 one must build on the CentOS 7 nodes:
#.  ssh ppplusername@portalc7.pppl.gov
#.and use the following modules
#.  module load gcc
#.  module load openmpi
#.  module swap gcc gcc/8.4.0
#.  module load cuda/10.2
#.  module load git
#.Also ensure that /sbin/ in in your PATH.

# if we are in machines directory, go up a directory
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT=$HOME/gkylsoft
export MACHINE_NAME='portal'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpicxx MPICC=mpicc MPICXX=mpicxx --prefix=$GKYLSOFT --build-gkylzero=yes --build-luajit=yes --build-adios=yes --build-openmpi=no
