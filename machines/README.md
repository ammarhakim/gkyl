There are two types of scripts in this directory for each machine
(usually a supercomputer center) that gkyl is supported on:

mkdeps.[machine].sh 
configure.[machine].sh

The dependencies should be build before configuring gkyl.

These scripts can be run from either within the machines (this)
directory or the main gkyl directory. Both scripts should should only
need to be run once, after which gkyl can be built with the usual

./waf build install

from the main gkyl directory.

To develop a build on an new machine, copy the mkdeps and configure
scripts from an existing machine, and make the following changes:

In mkdeps.[machine].sh:

-- For each of the depedencies (BUILD...), if the machine has existing
   installations (via modules, etc), set the corresponding BUILD flag
   to no

-- Change the compilers to those supported on the machine

In configure.[machine].sh :

-- As above, change the compilers to those supported on the machine

-- For each of the dependencies that were not built in the mkdeps
   file, specify paths and library names needed by the scripts in
   waf_tools. On many systems, this can be aided by loading modules.

