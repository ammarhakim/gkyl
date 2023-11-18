#.MF: It's possible this mkdeps (and/or the config) file for Ascent needs to be revised.
#.    Compilation fails at some medium to large size kernels, with error message:
#.      fatal error: error writing to /tmp/cc37dAA0.s: Cannot allocate memory.
#.    To overcome this problem, compile with ./waf build install -j1
#.    However, this eventually leads the to the following errors:
#.      [578/602] Compiling DataStruct/CartFieldDeviceImpl.cu
#.      /autofs/nccsopen-svm1_sw/ascent/gcc/8.1.1/include/c++/8.1.1/type_traits(347): error: identifier "__ieee128" is undefined
#.      
#.      /autofs/nccsopen-svm1_sw/ascent/gcc/8.1.1/include/c++/8.1.1/bits/std_abs.h(101): error: identifier "__ieee128" is undefined
#.      
#.      /autofs/nccsopen-svm1_sw/ascent/gcc/8.1.1/include/c++/8.1.1/bits/std_abs.h(102): error: identifier "__ieee128" is undefined
#.      
#.      3 errors detected in the compilation of "/tmp/tmpxft_000043b2_00000000-6_CartFieldDeviceImpl.cpp1.ii".
#.      
#.      Waf: Leaving directory `/autofs/nccsopen-svm1_home/manaf/gkeyll/gkyl/build'
#.      Build failed
#.       -> task in 'datastruct_cu' failed with exit status 1 (run with -v to display more information)

module unload xl/16.1.1-5
module load gcc/8.1.1
module load spectrum-mpi/10.3.1.2-20200121
module load adios/1.13.1-py2
# If we are in machines directory, go up a directory.
if [ `dirname "$0"` == "." ] 
  then
    cd ..
fi
export GKYLSOFT='~/gkylsoft'
cd install-deps
./mkdeps.sh CC=mpicc CXX=mpiCC MPICC=mpicc MPICXX=mpiCC --build-luajit=yes --build-adios=no --build-openmpi=no
