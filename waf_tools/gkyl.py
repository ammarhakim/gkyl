"""
Top-level GKYL paths
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('-p', type='string', help='Path to Gkyl dependency directory', dest='gkylDepsDir',
                   default=os.path.expandvars('$HOME/gkylsoft'))
    opt.add_option('--cxxflags', type='string', help='Compiler flags', dest='gkcxxflags',
                   default='-O3,-Wall,-std=c++14')
    opt.add_option('--prefix', type='string', help='Install path', dest='prefix',
                   default=os.path.expandvars('$HOME/gkylsoft/gkyl'))

@conf
def check_gkyl(conf):
    conf.start_msg("Setting dependency path:")
    conf.end_msg(conf.options.gkylDepsDir)

    conf.start_msg("Setting prefix:")
    conf.end_msg(conf.options.prefix)
    conf.env.PREFIX = conf.options.prefix

    conf.env.append_value('CXXFLAGS', conf.options.gkcxxflags.split(','))
    conf.env.append_value('CFLAGS', conf.options.gkcxxflags.split(','))
    conf.start_msg("Checking if CXXFLAGS work")
    conf.check_cxx(fragment="""#include<stdio.h>\nint main(){return 0;}\n""", execute=True)
    conf.end_msg("Flags work (%s)" % conf.options.gkcxxflags)
    
    return 1

def detect(conf):
    return detect_gkyl(conf)
