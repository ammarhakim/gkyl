"""
Determine if we should use sqlite3
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-sqlite', help=('Enable sqlite3'),
                   dest='enable_sqlite', action='store_true',
                   default=True)
    opt.add_option('--disable-sqlite', help=('Disable sqlite'),
                   dest='enable_sqlite', action='store_false')

@conf
def check_sqlite3(conf):
    opt = conf.options

    conf.start_msg('Checking for Sqlite3')
    
    conf.env['USE_SQLITE'] = True
    if not conf.options.enable_sqlite:
        conf.env['USE_SQLITE'] = False
        conf.end_msg("Not using Sqlite3")
    else:
        conf.define("USING_SQLITE3", 1)
        conf.end_msg("Using Sqlite3")
    
    return 1

def detect(conf):
    return detect_eigen(conf)
