"""
Detect Lapack includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--lapack-lib-dir', type='string', help='Path to LAPACK libraries', dest='lapackLibDir')
    opt.add_option('--lapack-link-libs', type='string', help='List of LAPACK libraries to link',
                   dest='lapackLinkLibs', default='lapack,blas')

@conf
def check_lapack(conf):
    opt = conf.options
    conf.env['LAPACK_FOUND'] = False
    
    if conf.options.lapackLibDir:
        conf.env.LIBPATH_LAPACK = conf.options.lapackLibDir
    else:
        conf.env.LIBPATH_LAPACK = [conf.options.gkylDepsDir+'/lapack/lib']

    libList = conf.options.lapackLinkLibs
    conf.env.LIB_LAPACK = libList.split(',')

    conf.start_msg('Checking for LAPACK')
    # the header_name below is completely random: this only checks for libraries
    conf.check(header_name='stdio.h', features='cxx cxxprogram', use="LAPACK", mandatory=True)
    conf.end_msg("Found LAPACK")
        
    conf.env['LAPACK_FOUND'] = True
    
    return 1

def detect(conf):
    return detect_lapack(conf)
