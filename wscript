## -*- python -*-
# Top-level build-script for Gkyl
##    _______     ___
## + 6 @ |||| # P ||| +

import datetime
import os
import platform
import sys

APPNAME = 'gkyl'
VER = "0.1"

now = datetime.datetime.now()
VERSION = VER + "-"+now.strftime("%Y-%m-%d")

top = '.'
out = 'build'

# extra flags to pass to linker
EXTRA_LINK_FLAGS = []

def options(opt):
    opt.load('compiler_c compiler_cxx') 
    opt.load('gkyl luajit mpi adios',
             tooldir='waf_tools')

def configure(conf):
    r"""Configure Gkyl build"""

    # load tools
    conf.load('compiler_c compiler_cxx')
    conf.check_gkyl()
    conf.check_luajit()
    conf.check_mpi()
    conf.check_adios()

    # standard install location for dependencies
    gkydepsDir = os.path.expandvars('$HOME/gkylsoft')

    # add current build directory to pick up config header
    conf.env.append_value('INCLUDES', ['.'])
    
    # load options for math and dynamic library
    conf.env.LIB_M = ['m']
    conf.env.LIB_DL = ['dl']

    conf.env.append_value("RPATH", conf.env.LIBDIR)

    # write out configuration info into header
    conf.write_config_header('gkylconfig.h')

def build(bld):
    # recurse down directories and build C++ code
    bld.recurse("Lib") 
    bld.recurse("Comm")
    bld.recurse("Unit")
    bld.recurse("Updater")
    bld.recurse("DataStruct")
    bld.recurse("Eq")

    # build executable
    buildExec(bld)    

    ### install LuaJIT code

    # - xsys
    xsys_dir = bld.path.find_dir('xsys')
    bld.install_files(
        "${PREFIX}/bin/xsys",
        xsys_dir.ant_glob('**/*.lua'),
        cwd=xsys_dir, relative_trick=True)
    
    # - Unit
    Unit_dir = bld.path.find_dir('Unit')
    bld.install_files(
        "${PREFIX}/bin/Unit",
        ["Unit/unit.lua", "Unit/init.lua"],
        cwd=Unit_dir, relative_trick=True)

    # - Lib
    Lib_dir = bld.path.find_dir('Lib')
    bld.install_files(
        "${PREFIX}/bin/Lib",
        Lib_dir.ant_glob('**/*.lua'),
        cwd=Lib_dir, relative_trick=True)

    # - Grid
    Grid_dir = bld.path.find_dir('Grid')
    bld.install_files(
        "${PREFIX}/bin/Grid",
        Grid_dir.ant_glob('**/*.lua'),
        cwd=Grid_dir, relative_trick=True)

    # - DataStruct
    DataStruct_dir = bld.path.find_dir('DataStruct')
    bld.install_files(
        "${PREFIX}/bin/DataStruct",
        DataStruct_dir.ant_glob('**/*.lua'),
        cwd=DataStruct_dir, relative_trick=True)

    # - Eq
    Eq_dir = bld.path.find_dir('Eq')
    bld.install_files(
        "${PREFIX}/bin/Eq",
        Eq_dir.ant_glob('**/*.lua'),
        cwd=Eq_dir, relative_trick=True)

    # - Updater
    Updater_dir = bld.path.find_dir('Updater')
    bld.install_files(
        "${PREFIX}/bin/Updater",
        Updater_dir.ant_glob('**/*.lua'),
        cwd=Updater_dir, relative_trick=True)

    # - App
    App_dir = bld.path.find_dir('App')
    bld.install_files(
        "${PREFIX}/bin/App",
        App_dir.ant_glob('**/*.lua'),
        cwd=App_dir, relative_trick=True)

    # - Comm
    Comm_dir = bld.path.find_dir('Comm')
    bld.install_files(
        "${PREFIX}/bin/Comm",
        Comm_dir.ant_glob('**/*.lua'),
        cwd=Comm_dir, relative_trick=True)

    # - Io
    Io_dir = bld.path.find_dir('Io')
    bld.install_files(
        "${PREFIX}/bin/Io",
        Io_dir.ant_glob('**/*.lua'),
        cwd=Io_dir, relative_trick=True)

    # - Basis
    Basis_dir = bld.path.find_dir('Basis')
    bld.install_files(
        "${PREFIX}/bin/Basis",
        Basis_dir.ant_glob('**/*.lua'),
        cwd=Basis_dir, relative_trick=True)

def buildExec(bld):
    r"""Build top-level executable"""
    if platform.system() == 'Darwin' and platform.machine() == 'x86_64':
        # we need to append special flags to get stuff to work on a 64 bit Mac
        EXTRA_LINK_FLAGS.append('-pagezero_size 10000 -image_base 100000000')

    # slightly modify Linux linker that thinks that he is smart and is
    # not exporting the symbols that are only used in the Lua part
    if platform.system() == 'Linux':
        bld.env.LINKFLAGS_cstlib = ['-Wl,-Bstatic,-E']
        bld.env.LINKFLAGS_cxxstlib = ['-Wl,-Bstatic,-E']
        bld.env.STLIB_MARKER = '-Wl,-Bstatic,-E'

    # build gkyl executable
    bld.program(
        source ='gkyl.cxx', target='gkyl',
        includes = 'Unit Lib Comm',
        use='lib datastruct eq unit comm updater basis LUAJIT ADIOS MPI M DL',
        linkflags = EXTRA_LINK_FLAGS,
        rpath = bld.env.RPATH,
        lib = 'pthread'
    )

def dist(ctx):
    ctx.algo = "zip" # use ZIP instead of tar.bz2
    ctx.excl = " **/.waf* **/*~ **/*.pyc **/*.swp **/.lock-w* configure-par.sh **/.hg **/.hgignore install-deps/build-opts.sh install-deps/luajit-2.0 install-deps/eigen-eigen-* install-deps/adios-1.* install-deps/luarocks-2.4.3* install-deps/openmpi-* build build-par build-ser"
