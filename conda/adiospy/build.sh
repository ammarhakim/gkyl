curl http://users.nccs.gov/~pnorbert/adios-1.11.0.tar.gz > adios-1.11.0.tar.gz
gunzip adios-1.11.0.tar.gz
tar -xvf adios-1.11.0.tar
cd adios-1.11.0/wrappers/numpy
# removing '-rlt' from one line in the Makefile that breaks the build on Macs
sed -e '33s/.*/        python setup.py build_ext/' -i ''  Makefile 
make python
$PYTHON setup.py install --prefix=$PREFIX
