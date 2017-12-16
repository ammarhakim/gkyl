Depending on your system, building dependencies can be complicated.
On a Mac or Linux machine you can simply run the mkdeps.sh script, but
first, please cheack details by running

./mkdeps.sh -h

On most supercomputers you will likely need to use the system
recommended compilers and MPI libraries. In this case, you should pass
the appropriate compilers to mkdeps.sh, and then build libraries not
provided by the system. In practice, this likely means LuaJIT, ADIOS
and, optionally LuaRocks and Eigen.