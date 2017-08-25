# About

This is the Gkyl code. The name is pronounced as in the book "The
Strange Case of Dr. Jekyll and Mr. Hyde". Gkyl is written in LuaJIT,
with some wrapper and time-critical parts written in C/C++. Gkyl is
developed at Princeton Plasma Physics Laboratory (PPPL).

# License

Gkyl can be used freely for research at universities, national
laboratories and other non-profit institutions. Release zip-balls will
be provided. Access to the source-code repository is restricted to
those who need to modify the code. In practice, this means people at
PPPL and our research partners who have jointly funded projects with
us.

# Dependencies

You must build and install the dependencies yourself, or use existing
builds for your system. Most supercomputer centers have optimized,
pre-built libraries for most dependencies. On these systems, you will
probably only need to install LuaJIT.

Build instructions for dependencies are provided in the build sections
below. Gkyl depends on the following tools and packages:

- A C/C++ compiler
- LuaJIT
- MPI
- ADIO IO library

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
as it is included with the distribution. To configure gkyl for
building do

- Copy the configure-par.sh-in script to configure-par.sh

- Modify this script to reflect the locations of various libraries on
  your machine. In particular, if you are using pre-built libraries
  you will likely need to change the information about MPI and ADIOS.

- Run the configure-par.sh script

Once the configuration is complete, run the following command to build
and install:

~~~~~~~
./waf build install
~~~~~~~

The builds are created in a 'build-par' directory. The executable is
build-par/gkyl. It takes a single parameter, the name of the LuaJIT script
to run. Note that, in general, you **can not** run the executable from
the build directory. (The reason for this is that besides the
executable a bunch of LuaJIT files are also needed to run most
simulations). Hence, the code **must** be installed as shown above.

If you need to clean up a build do:

~~~~~~~
./waf distclean
~~~~~~~

If you need to uninstall do:

~~~~~~~
./waf uninstall
~~~~~~~

If you want to zip up the code for copying to another machine, do:

~~~~~~~
./waf dist
~~~~~~~

This will create a zip file with the code. The repo hidden files and
generated files will not be zipped up.

## Note on building LuaJIT

To avoid linking to the shared LuaJIT library, please delete the
installed LuaJIT shared libraries.

LuaJIT builds easily on most machines with standard GCC
compiler. Often, you may run into problems on older gcc as they do not
include the log2 and exp2 functions unless c99 standard is enabled. To
get around this, modify the src/Makefile in LuaJIT. To do this, change
the line:

~~~~~~~
CC= $(DEFAULT_CC)
~~~~~~~

to:

~~~~~~~
CC= $(DEFAULT_CC) -std=gnu99
~~~~~~~
