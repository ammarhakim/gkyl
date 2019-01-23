Depending on your system, building dependencies can be complicated.
On a Mac (read below too) or Linux machine you can simply run the
mkdeps.sh script, but first, please cheack details by running

./mkdeps.sh -h

On most supercomputers you will likely need to use the system
recommended compilers and MPI libraries. In this case, you should pass
the appropriate compilers to mkdeps.sh, and then build libraries not
provided by the system. In practice, this likely means LuaJIT, ADIOS,
Eigen and, optionally LuaRocks.

After upgrading to MacOS 10.13 or 10.14, you may have trouble building
LuaJIT. One way around this is to declare the following environment
variable before running mkdeps.sh:

export MACOSX_DEPLOYMENT_TARGET=10.X

where 'X' is 13 (High Sierra) or 14 (Mojave).
