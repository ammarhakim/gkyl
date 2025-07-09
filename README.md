**This version of Gkeyll is NO LONGER SUPPORTED AND SHOULD NOT BE USED!!. It remains here only as an archive of past work**

# About

This is the Gkeyll code. The name is pronounced as in the book "The
Strange Case of Dr. Jekyll and Mr. Hyde" which is required reading for
all members of the Gkeyll Team. Gkeyll is written in a combination of
LuaJIT and C/C++.  Gkeyll is developed at Princeton Plasma Physics
Laboratory (PPPL) and is copyrighted 2016-2023 by Ammar Hakim and the
Gkeyll Team.

Documentation for the code is available at http://gkeyll.rtfd.io.

# Installation

Building gkyl requires

- a modern C/C++ compiler.
- Python 3 (for use in the `waf` build system and post-processing).
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
Finally, build the gkyl executable using
```
./waf build install
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

# Diagnostic tools

The `postgkyl` python package has been developed for plotting diagnostic
files from Gkeyll.  It can be installed via `conda` using

```
conda install -c gkyl -c conda-forge postgkyl
```

For more information about `postgkyl` and how to use it, please see
https://gkeyll.readthedocs.io/en/latest/postgkyl/main.html.

# Code contribution and formatting guidelines

All contributions to the code that improve the code via new
functionality and/or refactoring of existing functionality are
welcomes. Please strive for excellence in your programming and follow
carefully the rest of the code structure while doing so.


## Formatting guidelines

Formatting guidelines given below are meant to reduce the thought
given to minor (but asthetically important) issues. There are as many
opionions on how to format code as there are developers. Hence, in
Gkeyll these guidelines have been determined by the lead developer of
the code and are not open for further discussion.

- **Do not** modify existing code alignment or comments unless the code is
  wrong or the comment is incorrect, or if the formatting is
  egregiously bad.
  
- **Do not** align multiple consecutive statements with = signs.

- **Do not** mix tabs and spaces. Uses **spaces consistently**

- **Leave single space** between LHS and RHS expressions.

- You **may or may not** leave spaces between operators.

- You **may or may not** leave spaces after a comma in a function
  call.

- **Do not** comment obvious pieces of code.

- **Comment** function call signatures for user-facing functions.

# License

**Gkeyll can be used freely for research at universities, national
laboratories and other research institutions. 
If you want to use Gkeyll in a commercial environment,
please ask us first.**

We follow an *open-source but closed development model*. Even though
read access to the code is available to everyone, write access to the
source-code repository is restricted to those who need to modify the
code. In practice, this means researchers at PPPL and our partner
institutions. In particular, this means that for write access you
either need to have jointly funded projects or jointly supervised
graduate students or postdocs with Princeton University/PPPL.

In general, we allow users to "fork" the code to make their own
modifications. However, we would appreciate if you would work with us
to merge your features back into the main-line (if those features are
useful to the larger Gkeyll team). You can submit a "pull request" and
we will try our best to merge your changes into the
mainline. Contributed code should compile and have sufficient
unit/regression tests.

# Authors

Gkeyll is developed at the Princeton Plasma Physics Laboratory (PPPL),
a Department of Energy (DOE) national lab, managed by Princeton
University. Funding for the code comes from Department of Energy,
Airforce Office of Scientific Research, Advanced Projects Agency -
Energy, National Science Foundation and NASA.

The institutions involved in Gkeyll development are PPPL, Princeton
University, Virginia Tech, University of Maryland and MIT.

The CEO and Algorithm Alchemist of the project is Ammar Hakim.

The lead physicists for the project are Greg Hammett, Jason TenBarge
and Ammar Hakim.

The major contributors to the code are: Noah Mandell, Manaure (Mana)
Francisquez, Petr Cagas, James (Jimmy) Juno, Liang Wang and Tess
Bernard.

```
sed -i '' -e "s/[[:space:]]* =/ =/g"
```
