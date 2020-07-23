# About

This is the Gkeyll code. The name is pronounced as in the book "The
Strange Case of Dr. Jekyll and Mr. Hyde" which is required reading for
all members of the Gkeyll Team. Gkeyll is written in a combination of
LuaJIT and C++.  Gkeyll is developed at Princeton Plasma Physics
Laboratory (PPPL) and is copyrighted 2016-2020 by Ammar Hakim and the
Gkeyll Team.

Documentation for the code is available at http://gkeyll.rtfd.io.

# Getting dependencies and building the code

Building Gkeyll requires a modern C/C++ compiler and Python 3 (for use in the `waf` build system and post-processing). The following instructions assume that these tools are present.

For systems on which Gkeyll has been built before, the code can be built in three steps using scripts found in the `machines/` directory.
1. Install dependencies using a `mkdeps` script from the `machines/` directory:
```
./machines/mkdeps.[SYSTEM].sh
```
where `[SYSTEM]` should be replaced by the name of the system you are building on, such as `macosx` or `eddy`. By default, installations will be made in `~/gkylsoft/`. 

2. Configure `waf` using a `configure` script from the `machines/` directory: 
```
./machines/configure.[SYSTEM].sh
```

**Steps 1 and 2 should only need to be done on the first build, unless one wishes to change the dependencies.**

3. Build the code using
```
./waf build install
```

## Building on non-native systems.

For systems that do not already have corresponding files in the `machines/` directory, we encourage you to add files for your machine. Instructions can be found in `machines/README.md`.

## Testing the build.

As a preliminary test, just to make sure the `gkyl` executable is ok, you can do
```
~/gkylsoft/gkyl/bin/gkyl -v
```
This will print some version information and the libraries `gkyl` was built with.

You can run the full suite of unit tests using
```
cd Regression/

~/gkylsoft/gkyl/bin/gkyl runregression config

~/gkylsoft/gkyl/bin/gkyl runregression rununit
```

# Diagnostic tools

The `postgkyl` python package has been developed for plotting diagnostic files from Gkeyll. 
It can be installed via `conda` using

```
conda install -c gkyl postgkyl
```

For more information about `postgkyl` and how to use it, please see
https://gkeyll.readthedocs.io/en/latest/postgkyl/usage.html.

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
