curl http://users.nccs.gov/~pnorbert/adios-1.11.0.tar.gz > adios-1.11.0.tar.gz
gunzip adios-1.11.0.tar.gz
tar -xvf adios-1.11.0.tar
cd adios-1.11.0
CFLAGS="-fPIC" ./configure --prefix=$PREFIX --disable-fortran
make install
