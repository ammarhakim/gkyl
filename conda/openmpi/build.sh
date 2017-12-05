curl https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz > openmpi-3.0.0.tar.gz
gunzip -c openmpi-3.0.0.tar.gz | tar xf -
cd openmpi-3.0.0
CFLAGS="-fPIC" ./configure --prefix=$PREFIX --enable-mpi-fortran=none
make all install
