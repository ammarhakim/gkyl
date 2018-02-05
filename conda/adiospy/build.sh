curl -L http://users.nccs.gov/~pnorbert/adios-1.13.0.tar.gz > adios-1.13.0.tar.gz
gunzip adios-1.13.0.tar.gz
tar -xvf adios-1.13.0.tar
cd adios-1.13.0/wrappers/numpy
# removing '-rlt' from one line in the Makefile that breaks the build on Macs
sed 's/ -lrt//g' Makefile > tmp
mv tmp Makefile
make python
$PYTHON setup.py install --prefix=$PREFIX
