"""
Detect OPENBLAS includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-openblas', help=('Enable OPENBLAS'),
                   dest='enable_openblas', action='store_true',
                   default=True)
    opt.add_option('--openblas-inc-dir', type='string', help='Path to OPENBLAS includes', dest='openblasIncDir')
    opt.add_option('--openblas-lib-dir', type='string', help='Path to OPENBLAS libraries', dest='openblasLibDir')
    opt.add_option('--openblas-link-libs', type='string', help='List of libraries to link with OPENBLAS',
                   dest='openblasLinkLibs')

@conf
def check_openblas(conf):
    opt = conf.options
    conf.env['OPENBLAS_FOUND'] = False

    if not conf.options.enable_openblas:
        return
    
    if conf.options.openblasIncDir:
        conf.env.INCLUDES_OPENBLAS = conf.options.openblasIncDir
    else:
        conf.env.INCLUDES_OPENBLAS = conf.options.gkylDepsDir+'/OpenBLAS/include'

    if conf.options.openblasLibDir:
        #conf.env.LIBPATH_OPENBLAS = conf.options.openblasLibDir
        conf.env.LIBPATH_OPENBLAS = conf.options.openblasLibDir.split(',')
    else:
        conf.env.LIBPATH_OPENBLAS = conf.options.gkylDepsDir+'/OpenBLAS/lib'

    conf.env.LIB_OPENBLAS = ["openblas"]

    if conf.options.openblasLinkLibs:
        libList = conf.options.openblasLinkLibs
        conf.env.append_value('LIB_OPENBLAS', libList.split(','))
         
    conf.start_msg('Checking for OPENBLAS')
    conf.check(header_name='lapack.h', features='cxx cxxprogram', use="OPENBLAS", mandatory=True)
    conf.end_msg("Found OPENBLAS")
    conf.env['OPENBLAS_FOUND'] = True
    return 1

def detect(conf):
    return detect_openblas(conf)
