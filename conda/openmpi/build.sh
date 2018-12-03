curl -L https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.2.tar.gz > openmpi-3.1.2.tar.gz
gunzip -c openmpi-3.1.2.tar.gz | tar xf -
cd openmpi-3.1.2
./configure --prefix=$PREFIX --enable-mpi-fortran=none CC=$CC CXX=$CXX
make all install
