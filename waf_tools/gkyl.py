"""
Top-level GKYL paths
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('-p', type='string', help='Path to Gkyl dependency directory', dest='gkylDepDir',
                   default=os.path.expandvars('$HOME/gkylsoft'))

@conf
def check_gkyl(conf):
    conf.start_msg("Setting dependency path:")
    conf.end_msg(conf.options.gkylDepDir)
    return 1

def detect(conf):
    return detect_gkyl(conf)
