curl https://www.open-mpi.org/software/ompi/v2.1/downloads/openmpi-2.1.1.tar.gz > openmpi-2.1.1.tar.gz
gunzip -c openmpi-2.1.1.tar.gz | tar xf -
cd openmpi-2.1.1
./configure --prefix=$PREFIX 
make all install
