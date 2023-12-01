# -*- python -*-
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

from waflib import TaskGen
from waflib.Options import options

def options(opt):
    opt.load('compiler_c compiler_cxx') 
    opt.load('gkyl luajit mpi adios sqlite3 cutools nccl gkylzero superlu openblas',
             tooldir='waf_tools')

def configure(conf):
    r"""Configure Gkyl build"""

    # load tools
    conf.load('compiler_c compiler_cxx')
    conf.load('cutools', tooldir='waf_tools')
    conf.check_gkyl()
    conf.check_luajit()
    conf.check_mpi()
    conf.check_adios()
    conf.check_sqlite3()
    conf.check_cutools()
    conf.check_nccl()
    conf.check_superlu()
    conf.check_openblas()
    conf.check_gkylzero()

    # standard install location for dependencies
    gkydepsDir = os.path.expandvars('$HOME/gkylsoft')

    # add current build directory to pick up config header
    conf.env.append_value('INCLUDES', ['.'])
    
    # load options for math and dynamic library
    conf.env.LIB_M = ['m']
    conf.env.LIB_DL = ['dl']

    # write out configuration info into header
    conf.write_config_header('gkylconfig.h')


from waflib import Task
class GitTip(Task.Task):
    always_run = True # need to force running every time
    
    run_str = r'echo \#define GKYL_GIT_CHANGESET  \"`git describe --abbrev=12 --always --dirty=+`\" > ${TGT}'

def build(bld):

    gitTip = GitTip(env=bld.env)
    gitTip.set_outputs(bld.path.find_or_declare('gkylgittip.h'))
    bld.add_to_group(gitTip)

    if bld.jobs > 16:
      bld.jobs = 16
    
    # recurse down directories and build C++ code
    bld.recurse("Comm")
    bld.recurse("Eq")
    bld.recurse("Grid")
    bld.recurse("Lib")
    bld.recurse("Proto")
    bld.recurse("Tool")
    bld.recurse("Unit")
    bld.recurse("Updater")
    if bld.env['CUTOOLS_FOUND']:
        bld.recurse("Cuda")

    # Sometimes there is an issue with an existing build of sqlite on
    # a Linux machine. In that case, sqlite support can be
    # disabled. WARNING: This means the regression testing system
    # won't work.
    if bld.env['USE_SQLITE']:
        bld.recurse("sqlite3")

    # build executable
    buildExec(bld)

    ### install LuaJIT code

    # - xsys
    xsys_dir = bld.path.find_dir('xsys')
    bld.install_files(
        "${PREFIX}/bin/xsys",
        xsys_dir.ant_glob('**/*.lua'),
        cwd=xsys_dir, relative_trick=True)

    # - sci
    sci_dir = bld.path.find_dir('sci')
    bld.install_files(
        "${PREFIX}/bin/sci",
        sci_dir.ant_glob('**/*.lua'),
        cwd=sci_dir, relative_trick=True)
    
    # - Unit
    Unit_dir = bld.path.find_dir('Unit')
    bld.install_files(
        "${PREFIX}/bin/Unit",
        Unit_dir.ant_glob('**/*.lua'),
        cwd=Unit_dir, relative_trick=True)

    # - Regression
    Regression_dir = bld.path.find_dir('Regression')
    bld.install_files(
        "${PREFIX}/bin/Regression",
        Regression_dir.ant_glob('**/*.lua'),
        cwd=Regression_dir, relative_trick=True)

    # - Tool
    Tool_dir = bld.path.find_dir('Tool')
    bld.install_files(
        "${PREFIX}/bin/Tool",
        Tool_dir.ant_glob('**/*.lua'),
        cwd=Tool_dir, relative_trick=True)
    bld.install_files(
        "${PREFIX}/bin/Tool",
        Tool_dir.ant_glob('**/*.md'),
        cwd=Tool_dir, relative_trick=True)

    # - Sqlite3
    Sqlite_dir = bld.path.find_dir('sqlite3')
    bld.install_files(
        "${PREFIX}/bin/sqlite3",
        Sqlite_dir.ant_glob('**/*.lua'),
        cwd=Sqlite_dir, relative_trick=True)    

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

    # - Cuda
    Cuda_dir = bld.path.find_dir('Cuda')
    bld.install_files(
        "${PREFIX}/bin/Cuda",
        Cuda_dir.ant_glob('**/*.lua'),
        cwd=Cuda_dir, relative_trick=True)

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

    # - Proto/
    Proto_dir = bld.path.find_dir('Proto')
    bld.install_files(
        "${PREFIX}/bin/Proto",
        Proto_dir.ant_glob('**/*.lua'),
        cwd=Proto_dir, relative_trick=True)

def appendToList(target, val):
    if type(val) == list:
        target.extend(val)
    else:
        target.append(val)
        
def buildExec(bld):
    r"""Build top-level executable"""

    # Link flags on Linux
    if platform.system() == 'Linux':
        bld.env.LINKFLAGS_cstlib = ['-Wl,-Bstatic,-E']
        bld.env.LINKFLAGS_cxxstlib = ['-Wl,-Bstatic,-E']
        bld.env.STLIB_MARKER = '-Wl,-Bstatic,-E'
        bld.env.SHLIB_MARKER = '-Wl,-Bdynamic,--no-as-needed'

    # list of objects to use
    useList = 'lib datastruct eq unit comm updater tool proto basis grid gkylzero LUAJIT ADIOS MPI SUPERLU OPENBLAS M DL'
    if bld.env['USE_SQLITE']:
        useList = ' sqlite3 ' + useList
    if bld.env['CUTOOLS_FOUND']:
        useListCuda = ' cuda CUTOOLS unit_cu updater_cu NCCL comm_cu '
        useList = useListCuda + useList

    # set RPATH
    fullRpath = []
    appendToList(fullRpath, bld.env.RPATH) # user specified RPATH
    appendToList(fullRpath, bld.env.LIBDIR)
    appendToList(fullRpath, bld.env.LIBPATH_LUAJIT)
    appendToList(fullRpath, bld.env.LIBPATH_gkylzero)
    appendToList(fullRpath, bld.env.LIBPATH_ADIOS)


    # build gkyl executable

    source = ['gkyl.cxx']
    if bld.env['CUTOOLS_FOUND']:
        bld(
            target='culink',
            use = useListCuda,
            features = 'cxx culink',
            name = 'culink'
        )
        bld.add_group() 
        appendToList(source, 'culink.o')

    bld.program(
        source = source, target='gkyl',
        includes = 'Unit Lib Comm',
        use = useList,
        linkflags = EXTRA_LINK_FLAGS,
        rpath = fullRpath,
        lib = 'pthread ' + bld.env.EXTRALIBS,
    )

def dist(ctx):
    ctx.algo = "zip" # use ZIP instead of tar.bz2
    ctx.excl = " **/.waf* **/*~ **/*.pyc **/*.swp **/.lock-w* configure-par.sh **/.hg **/.hgignore install-deps/build-opts.sh install-deps/luajit-2.0 install-deps/adios2-2.* install-deps/luarocks-2.4.3* install-deps/openmpi-* build build-par build-ser"
