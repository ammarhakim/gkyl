#!/bin/sh

# Defaults
PREFIX=$HOME/gkylsoft

# default build options
CC=gcc
CXX=g++
MPICC=mpicc
MPICXX=mpicxx

BUILD_LUAJIT=yes
BUILD_LUAROCKS=no
BUILD_ADIOS=no
BUILD_OPENMPI=no
BUILD_EIGEN=no

# ----------------------------------------------------------------------------
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------------

# Help

show_help() {
cat <<EOF

./mkdeps.sh CC=cc CXX=cxx MPICC=mpicc MPICXX=mpicxx

Build Gkyl dependencies. By default, only LuaJIT is built.

CC 
CXX
MPICC
MPICXX                      C, C++, MPI C and MPI C++ compilers to use

-h
--help                      This help.
--prefix=DIR                Prefix where dependencies should be installed.
                            Default is $HOME/gkylsoft

The following flags specify which libraries to build. By default, only
LuaJIT is built. If you build libraries that depend on MPI please
specify the MPI C and C++ compilers to use.

--build-luajit              [yes] Should we build LuaJIT?
--build-luarocks            [no] Should we build Luarocks?
--build-adios               [no] Should we build ADIOS?
--build-openmpi             [no] Should we build OpenMPI?
--build-eigen               [no] Should we build Eigen?
EOF
}

# Helper functions

find_program() {
   prog=`command -v "$1" 2>/dev/null`
   if [ -n "$prog" ]
   then
      dirname "$prog"
   fi
}

die() {
   echo "$*"
   echo
   echo "Dependency builds failed."
   echo
   exit 1
}

# ----------------------------------------------------------------------------
# MAIN PROGRAM
# ----------------------------------------------------------------------------

# Parse options

while [ -n "$1" ]
do
   value="`echo $1 | sed 's/[^=]*.\(.*\)/\1/'`"
   key="`echo $1 | sed 's/=.*//'`"
   if `echo "$value" | grep "~" >/dev/null 2>/dev/null`
   then
      echo
      echo '*WARNING*: the "~" sign is not expanded in flags.'
      echo 'If you mean the home directory, use $HOME instead.'
      echo
   fi
   case "$key" in
   -h)
      show_help
      exit 0
      ;;
   --help)
      show_help
      exit 0
      ;;
   CC)
      [ -n "$value" ] || die "Missing value in flag $key."
      CC="$value"
      ;;
   CXX)
      [ -n "$value" ] || die "Missing value in flag $key."
      CXX="$value"
      ;;
   MPICC)
      [ -n "$value" ] || die "Missing value in flag $key."
      MPICC="$value"
      ;;
   MPICXX)
      [ -n "$value" ] || die "Missing value in flag $key."
      MPICXX="$value"
      ;;   
   --prefix)
      [ -n "$value" ] || die "Missing value in flag $key."
      PREFIX="$value"
      ;;
   --build-openmpi)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_OPENMPI="$value"
      ;;
   --build-eigen)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_EIGEN="$value"
      ;;
   --build-luajit)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_LUAJIT="$value"
      ;;
   --build-adios)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_ADIOS="$value"
      ;;
   --build-luarocks)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_LUAROCKS="$value"
      ;;
   *)
      die "Error: Unknown flag: $1"
      ;;
   esac
   shift
done

# if mpicc doesn't work (because it doesn't exist or it's not in path), try to use installed openmpi version
if ! [ -x "$(command -v $MPICC)" ]
then
    MPICC=$PREFIX/openmpi-3.0.0/bin/mpicc
    MPICXX=$PREFIX/openmpi-3.0.0/bin/mpicxx
fi
# if mpicc still doesn't work, force to install openmpi
if ! [ -x "$(command -v $MPICC)" ] 
then
    BUILD_OPENMPI="yes"
fi

# Write out build options for scripts to use
cat <<EOF1 > build-opts.sh
# Generated automatically! Do not edit

# Installation directory
GKYLSOFT=$PREFIX
# Various compilers
CC=$CC
CXX=$CXX
MPICC=$MPICC
MPICXX=$MPICXX

EOF1

build_openmpi() {
    if [ "$BUILD_OPENMPI" = "yes" ]
    then
	echo "Building OpenMPI"
	./build-openmpi.sh
    fi
}

build_eigen() {
    if [[ "$BUILD_EIGEN" = "yes" || ! -f $PREFIX/eigen3/include/eigen3/Eigen/Core ]]
    then
	echo "Building EIGEN"
	./build-eigen.sh
    fi
}

build_luajit() {
    if [[ "$BUILD_LUAJIT" = "yes" || ! -f $PREFIX/luajit/include/luajit-2.1/lua.hpp ]]
    then    
	echo "Building LuaJIT"
	./build-luajit-beta3.sh
    fi
}

build_adios() {
    if [[ "$BUILD_ADIOS" = "yes" || ! -f $PREFIX/adios/include/adios.h ]]
    then    
	echo "Building ADIOS"
	./build-adios.sh
    fi
}

build_luarocks() {
    if [ "$BUILD_LUAROCKS" = "yes" ]
    then    
	echo "Building Luarocks"
	./build-luarocks.sh
    fi
}

echo "Installations will be in $PREFIX"

build_openmpi
build_luajit
build_luarocks
build_adios
build_eigen
