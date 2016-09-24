wget http://users.nccs.gov/~pnorbert/adios-1.10.0.tar.gz
gunzip adios-1.10.0.tar.gz
tar -xvf adios-1.10.0.tar
cd adios-1.10.0
./configure --prefix=$HOME/gkylsoft/adios --disable-fortran
make install
