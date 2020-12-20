curl -L https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz > openmpi-4.0.5.tar.gz
gunzip -c openmpi-4.0.5.tar.gz | tar xf -
cd openmpi-4.0.5
./configure --prefix=$PREFIX --enable-mpi-fortran=none CC=$CC CXX=$CXX
make all install
