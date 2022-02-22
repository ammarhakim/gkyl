"""
Detect gkylzero includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-gkylzero', help=('Enable gkylzero'),
                   dest='enable_gkylzero', action='store_true',
                   default=True)
    opt.add_option('--gkylzero-inc-dir', type='string', help='Path to gkylzero includes', dest='gkylzeroIncDir')
    opt.add_option('--gkylzero-lib-dir', type='string', help='Path to gkylzero libraries', dest='gkylzeroLibDir')
    opt.add_option('--gkylzero-link-libs', type='string', help='List of libraries to link with gkylzero',
                   dest='gkylzeroLinkLibs', default='gkylzero')

@conf
def check_gkylzero(conf):
    opt = conf.options
    conf.env['gkylzero_FOUND'] = False

    if not conf.options.enable_gkylzero:
        return
    
    if conf.options.gkylzeroIncDir:
        conf.env.INCLUDES_gkylzero = conf.options.gkylzeroIncDir
    else:
        conf.env.INCLUDES_gkylzero = conf.options.gkylDepsDir+'/gkylzero/include'

    if conf.options.gkylzeroLibDir:
        conf.env.LIBPATH_gkylzero = conf.options.gkylzeroLibDir.split(',')
    else:
        conf.env.LIBPATH_gkylzero = conf.options.gkylDepsDir+'/gkylzero/lib'

    libList = conf.options.gkylzeroLinkLibs
    conf.env.LIB_gkylzero = libList.split(',')
         
    conf.start_msg('Checking for gkylzero')
    conf.check(header_name='gkylzero.h', features='cxx cxxprogram', use="gkylzero SUPERLU OPENBLAS", mandatory=True)
    conf.end_msg("Found gkylzero")
    conf.env['gkylzero_FOUND'] = True
    return 1

def detect(conf):
    return detect_gkylzero(conf)
