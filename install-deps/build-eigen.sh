wget http://bitbucket.org/eigen/eigen/get/3.2.9.tar.bz2
bunzip2 3.2.9.tar.bz2
tar -xvf 3.2.9.tar
cd eigen-eigen-*
mkdir build-dir; cd build-dir
cmake ../ -DCMAKE_INSTALL_PREFIX=$HOME/gkylsoft
make install

