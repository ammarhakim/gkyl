"""
Detect Lapack includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--lapack-lib-dir', type='string', help='Path to LAPACK libraries', dest='lapackLibDir')

@conf
def check_lapack(conf):
    opt = conf.options
    conf.env['LAPACK_FOUND'] = False
    
    if conf.options.lapackLibDir:
        conf.env.LIBPATH_LAPACK = conf.options.lapackLibDir
    else:
        conf.env.LIBPATH_LAPACK = [conf.options.gkylDepsDir+'/lapack/lib']

    if COMPILER_CXX == 'g++':
        conf.env.LIB_LAPACK = ["lapack"]
    elif COMPILER_CXX == 'clang':
        pass
        
    conf.env['LAPACK_FOUND'] = True
    
    return 1

def detect(conf):
    return detect_lapack(conf)
