# About

This is the Gkeyll code. The name is pronounced as in the book "The
Strange Case of Dr. Jekyll and Mr. Hyde" which is required reading for
all members of the Gkeyll Team. Gkeyll is written in a combination of
LuaJIT and C/C++.  Gkeyll is developed at Princeton Plasma Physics
Laboratory (PPPL) and is copyrighted 2016-2024 by Ammar Hakim and the
Gkeyll Team.

Documentation for the code is available at http://gkeyll.rtfd.io.

# Installation

Building gkyl requires

- A modern C/C++ compiler.
- NVIDIA Cuda compiler (nvcc) if using (NVIDIA) GPUs.
- The NCCL multi-GPU communication library.
- gkylzero (https://github.com/ammarhakim/gkylzero) compiled and
  installed.

The following instructions assume that these tools are present.

Installing gkyl consists of three steps: building or indicating,
dependencies, configuring, and building the executable. If we have built
gkyl in the computer of interest before, we likely saved scripts
(machine files) that simplify this process. If we haven't, you'll
need to either build new machine files or perform each of the
installation steps manually.                                 

## On a known computer using machine files

We provide a set of "machine files" to ease the build process.
These are stored in the `machines/` directory. For example, to
build on Perlmutter please run
```
./machines/mkdeps.perlmutter.sh
./machines/configure.perlmutter.sh
```
The first of these will install whichever dependencies are needed
(e.g. LuaJIT). The default is for these installations take place
in ```$HOME/gkylsoft```; if a different install directory is desired,
specify it via the ```--prefix``` argument. The second of these
steps tells our ```waf``` build system where to find the
dependencies and where to place the gkyl executable
(```$HOME/gkylsoft``` by default, but a non-default directory can
be specified by changing GKYLSOFT). These two steps only need to be
done once, unless one wishes to change the dependencies.

Next, manually load the modules listed at the top of the
```machines/configure.<machine name>.sh``` file. For example,
for Perlmutter do:
```
module load PrgEnv-gnu/8.3.3
module load cray-mpich/8.1.22
module load python/3.9-anaconda-2021.11
module load cudatoolkit/11.7
module load nccl/2.15.5-ofi
module unload darshan
```
Finally, build and install the gkyl executable using
```
make install
```
The result will be a `gkyl` executable located in the
`$HOME/gkylsoft/gkyl/bin/` directory.

## On a new computer (no machine files available).

For systems that do not already have corresponding files in the
`machines/` directory, we encourage you to author machine files
for your machine following the existing ones as guides.
Instructions can be found in `machines/README.md`.

## Testing the build.

As a preliminary test, just to make sure the `gkyl` executable is
ok, you can do
```
$HOME/gkylsoft/gkyl/bin/gkyl -v
```
This will print some version information and the libraries `gkyl`
 was built with. Since gkyl is a parallel code, and some clusters
don't allow simply calling the `gkyl` executable (especially on the
login node), you may have to use `mpirun`, `mpiexec` or `srun`
(see your cluster's documentation) to run gkyl with, for example,
```
srun -n 1 $HOME/gkylsoft/gkyl/bin/gkyl -v
```

You can run a regression test as a first simulation. For example,
to run the Vlasov-Maxwell 2x2v Weibel regression test on a CPU, do
```
cd Regression/vm-weibel/
srun -n 1 $HOME/gkylsoft/gkyl/bin/gkyl rt-weibel-2x2v-p2.lua
```
and to run it on a GPU you may use
```
srun -n 1 $HOME/gkylsoft/gkyl/bin/gkyl -g rt-weibel-2x2v-p2.lua
```

You can run the full suite of unit tests using
```
cd Regression/
$HOME/gkylsoft/gkyl/bin/gkyl runregression config
$HOME/gkylsoft/gkyl/bin/gkyl runregression rununit
```
