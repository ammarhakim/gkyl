curl -L http://users.nccs.gov/~pnorbert/adios-1.13.1.tar.gz > adios-1.13.1.tar.gz
gunzip adios-1.13.1.tar.gz
tar -xvf adios-1.13.1.tar
cd adios-1.13.1
./configure --prefix=$PREFIX --disable-fortran --without-netcdf CFLAGS="-fPIC" CC=$CC CXX=$CXX MPICC=$MPICC MPICXX=$MPICXX --enable-shared=no
make install
