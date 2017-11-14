# Build directory
OUT=build-par

# Compile flags (set optimization/debug flags here)
CXXFLAGS='-O3,-Wall'

# LuaJIT options
LUAJIT_INC_DIR=$PREFIX/include/luajit-2.1
LUAJIT_LIB_DIR=$PREFIX/lib
LUAJIT_SHARE_DIR=$PREFIX/share/luajit-2.1.0-beta3

## MPI options
ENABLE_MPI="--enable-mpi"
MPI_INC_DIR=$PREFIX/include
MPI_LIB_DIR=$PREFIX/lib
MPI_LINK_LIBS="mpi"

# ADIOS options
ENABLE_ADIOS="--enable-adios" # set to blank to disable ADIOS
ADIOS_INC_DIR=$PREFIX/include
ADIOS_LIB_DIR=$PREFIX/lib

# You probably do not need to modify the command itself
cmd="./waf --out=$OUT --prefix=$PREFIX --cxxflags=$CXXFLAGS --luajit-inc-dir=$LUAJIT_INC_DIR --luajit-lib-dir=$LUAJIT_LIB_DIR --luajit-share-dir=$LUAJIT_SHARE_DIR  $ENABLE_MPI --mpi-inc-dir=$MPI_INC_DIR --mpi-lib-dir=$MPI_LIB_DIR --mpi-link-libs=$MPI_LINK_LIBS $ENABLE_ADIOS --adios-inc-dir=$ADIOS_INC_DIR --adios-lib-dir=$ADIOS_LIB_DIR configure build install"
echo $cmd
$cmd
