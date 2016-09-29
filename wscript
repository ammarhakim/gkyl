## -*- python -*-
# Top-level build-script for Gkyl
##    _______     ___
## + 6 @ |||| # P ||| +

import os
import sys
sys.path.insert(0, './waf_tools')
import commands

APPNAME = 'gkyl'
VERSION = '1.0'

top = '.'
out = 'build'

# extra flags to pass to linker
EXTRA_LINK_FLAGS = []

def options(opt):
    opt.load('compiler_c compiler_cxx') 
    opt.load('gkyl')
    opt.load('luajit') 
    opt.load('mpi')   

def configure(conf):
    r"""Configure Gkyl build"""

    # load tools
    conf.load('compiler_c compiler_cxx')
    conf.check_gkyl()
    conf.check_luajit()
    conf.check_mpi()

    # standard install location for dependencies
    gkydepsDir = os.path.expandvars('$HOME/gkylsoft')

    # add current build directory to pick up config header
    conf.env.append_value('INCLUDES', ['.'])
    
    # load options for LuaJIT
    conf.env.INCLUDES_LUAJIT = [gkydepsDir+'/luajit/include/luajit-2.1']
    conf.env.STLIB_LUAJIT = ['luajit-5.1']
    conf.env.STLIBPATH_LUAJIT = [gkydepsDir+'/luajit/lib']

    # load options for math and dynamic library
    conf.env.LIB_M = ['m']
    conf.env.LIB_DL = ['dl']

    # write out configuration info into header
    conf.write_config_header('gkylconfig.h')

def build(bld):
    ### recurse down directories and build C++ code
    bld.recurse("Unit")

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

    # Install wrapper shell script
    bld.install_files("${PREFIX}/bin", "xgkyl", chmod=0755)

def buildExec(bld):
    r"""Build top-level executable"""
    uname = os.uname()
    if uname[0] == 'Darwin' and uname[4] == 'x86_64':
        # we need to append special flags to get stuff to work on a Mac
        EXTRA_LINK_FLAGS.append('-pagezero_size 10000 -image_base 100000000')

    # build gkyl executable
    bld.program(
        source='gkyl.cxx', target='gkyl',
        includes = 'Unit',
        use='gkunit LUAJIT M DL MPI',
        linkflags = EXTRA_LINK_FLAGS
    )

def dist(ctx):
    ctx.algo = "zip" # use ZIP instead of tar.bz2
    ctx.excl = " **/.waf-1* **/*~ **/*.pyc **/*.swp **/.lock-w* **/.hg **/.hgignore build install-deps/luajit-2.0 install-deps/eigen-eigen-*"
