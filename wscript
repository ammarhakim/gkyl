## -*- python -*-
# Top-level build-script for Gkyl

import os

APPNAME = 'gkyl'
VERSION = '1.0'

top = '.'
out = 'build'

# extra flags to pass to linker
EXTRA_LINK_FLAGS = []

def options(opt):
    opt.load('compiler_c compiler_cxx')

def configure(conf):
    r"""Configure Gkyl build"""

    # load tools
    conf.load('compiler_c compiler_cxx')

    # standard install location for dependencies
    gkydepsDir = os.path.expandvars('$HOME/software')
    
    # load options for LuaJIT
    conf.env.INCLUDES_LUAJIT = [gkydepsDir+'/luajit/include/luajit-2.1']
    conf.env.STLIB_LUAJIT = ['luajit-5.1']
    conf.env.STLIBPATH_LUAJIT = [gkydepsDir+'/luajit/lib']

    # load options for math and dynamic library
    conf.env.LIB_M = ['m']
    conf.env.LIB_DL = ['dl']

def build(bld):
    # recurse down directories

    # build executable
    buildExec(bld)

def buildExec(bld):
    r"""Build top-level executable"""
    uname = os.uname()
    if uname[0] == 'Darwin' and uname[4] == 'x86_64':
        # we need to append special flags to get stuff to work on a Mac
        EXTRA_LINK_FLAGS.append('-pagezero_size 10000 -image_base 100000000')

    # build idjit executable
    bld.program(
        source='gkyl.cxx', target='gkyl',
        use='LUAJIT M DL',
        linkflags = EXTRA_LINK_FLAGS
    )

def dist(ctx):
    ctx.algo = "zip" # use ZIP instead of tar.bz2
    ctx.excl = " **/.waf-1* **/*~ **/*.pyc **/*.swp **/.lock-w* **/.hg **/.hgignore build install-deps/luajit-2.0"
