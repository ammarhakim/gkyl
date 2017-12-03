#!/bin/sh

# Defaults
PREFIX=$HOME/gkylsoft

# default build options
BUILD_LUAJIT=yes
BUILD_ADIOS=yes
BUILD_OPENMPI=no
BUILD_EIGEN=no
BUILD_PETSC=no
BUILD_LUAROCKS=no

# ----------------------------------------------------------------------------
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------------

# Help

show_help() {
cat <<EOF
Build Gkyl dependencies. Default libraries built are: LuaJIT and ADIOS.

--help                      This help.
--prefix=DIR                Prefix where dependencies should be installed.
                            Default is $HOME/gkylsoft
--build-openmpi             [no] Should we build OpenMPI?
--build-eigen               [no] Should we build Eigen?
--build-petsc               [no] Should we build Petsc?
--build-luajit              [yes] Should we build LuaJIT?
--build-adios               [yes] Should we build ADIOS?
--build-luarocks            [no] Should we build Luarocks?
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
   echo "configure failed."
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
   --help)
      show_help
      exit 0
      ;;
   --prefix)
      [ -n "$value" ] || die "Missing value in flag $key."
      PREFIX="$value"
      PREFIX_SET=yes
      ;;
   --build-openmpi)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_OPENMPI="$value"
      ;;
   --build-eigen)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_EIGEN="$value"
      ;;
   --build-petsc)
      [ -n "$value" ] || die "Missing value in flag $key."
      BUILD_PETSC="$value"
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

echo "# Configuration information. Generated automatically. Do not edit!

# location where Gkyl dependencies will be installed
GKYLSOFT=$PREFIX
" > build-opts.sh

build_openmpi() {
    if [ "$BUILD_OPENMPI" = "yes" ]
    then
	echo "Building OpenMPI"
	./build-openmpi.sh
    else
	# write information to build-opts
	echo "
# MPI library locations
MPICC=mpicc
MPICXX=mpicxx
" >> build-opts.sh
    fi
}

build_eigen() {
    if [ "$BUILD_EIGEN" = "yes" ]
    then
	echo "Building EIGEN"
	./build-eigen.sh
    fi
}

build_luajit() {
    if [ "$BUILD_LUAJIT" = "yes" ]
    then    
	echo "Building LuaJIT"
	./build-luajit-beta3.sh
    fi
}

build_adios() {
    if [ "$BUILD_ADIOS" = "yes" ]
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

echo "Installations will be in  $PREFIX"

build_openmpi
build_luajit
build_luarocks
build_adios
build_eigen
