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
    opt.add_option('--cflags', type='string', help='Compiler flags', dest='gkcflags',
                   default='-O3,-Wall')
    opt.add_option('--debug', help='Debug flags', dest='gkdebug',
                   action='store_true', default=False)
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
    conf.env.append_value('CFLAGS', conf.options.gkcflags.split(','))
    conf.env.append_value('LDFLAGS', conf.options.gkcflags.split(','))
    if conf.options.gkdebug:
      conf.env.append_value('CXXFLAGS', '-g')
      conf.env.append_value('CFLAGS', '-g')
     
    #conf.start_msg("Checking if CXXFLAGS work")
    #conf.check_cxx(fragment="""#include<stdio.h>\nint main(){return 0;}\n""", execute=True)
    #conf.end_msg("Flags work (%s)" % conf.options.gkcxxflags)

    #conf.start_msg("Checking if CFLAGS work")
    #conf.check_cxx(fragment="""#include<stdio.h>\nint main(){return 0;}\n""", execute=True)
    #conf.end_msg("Flags work (%s)" % conf.options.gkcflags)
    
    return 1

def detect(conf):
    return detect_gkyl(conf)
