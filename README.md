# About

This is the Gkyl code. The name is pronounced as in the book "The
Strange Case of Dr. Jekyll and Mr. Hyde". Gkyl is written in LuaJIT,
with time-critical parts written in C++. Gkyl is developed at
Princeton Plasma Physics Laboratory (PPPL).

# License

Gkyl can be used freely for research at universities, national
laboratories and other non-profit institutions. Nightly zip-balls will
be provided. Access to the source-code repository is restricted to
PPPL and our research partners.

# Dependencies

You must build and install the dependencies yourself, or use existing
builds for your system. Most supercomputer centers have optimized,
pre-built libraries. On these systems, you will only need to install
LuaJIT yourself, and perhaps Eigen. All other dependencies should
already exist.

Build instructions for dependencies are provided in the build sections
below. Gkyl depends on the following tools and packages:

- A C/C++ compiler supporting C++11 extensions.
- LuaJIT
- HDF5
- MPI
- Eigen

Optionally, you will need

- Petsc for linear and non-linear solvers.

# Building Gkyl and its dependencies

## Building dependencies

Building dependencies can be complicated. Some example scripts are
provided in the install-deps directory. You can use these as-is or
modify them as needed to suit your system. If you build your own
dependencies and want Waf to find them automatically, install them in
the $HOME/gkylsoft directory. Obviously, the dependencies need to be
built only once.

## Building Gkyl

Once you have all dependencies installed, you can build Gkyl
itself. Gkyl uses the Waf build system. You do NOT need to install waf
as it is included with the distribution. To build the gkyl executable
do:

  ./waf configure build

The builds are created in a 'build' directory. The executable is
build/gkyl. It takes a single parameter, the name of the LuaJIT script
to run. Note that, in general, you **can not** run the executable from
the build directory. (The reason for this is that besides the
executable a bunch of LuaJIT files are also needed to run most
simulations). Hence, the code **must** be installed using the
following command (the prefix can be whatever you please, of course):

  ./waf install --prefix=$HOME/gkylsoft

If you need to clean up a build do:

  ./waf distclean

If you need to uninstall do:

  ./waf uninstall

