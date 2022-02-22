"""
Detect SUPERLU includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-superlu', help=('Enable SUPERLU'),
                   dest='enable_superlu', action='store_true',
                   default=True)
    opt.add_option('--superlu-inc-dir', type='string', help='Path to SUPERLU includes', dest='superluIncDir')
    opt.add_option('--superlu-lib-dir', type='string', help='Path to SUPERLU libraries', dest='superluLibDir')
    opt.add_option('--superlu-link-libs', type='string', help='List of libraries to link with SUPERLU',
                   dest='superluLinkLibs')

@conf
def check_superlu(conf):
    opt = conf.options
    conf.env['SUPERLU_FOUND'] = False

    if not conf.options.enable_superlu:
        return
    
    if conf.options.superluIncDir:
        conf.env.INCLUDES_SUPERLU = conf.options.superluIncDir
    else:
        conf.env.INCLUDES_SUPERLU = conf.options.gkylDepsDir+'/superlu/include'

    if conf.options.superluLibDir:
        #conf.env.STLIBPATH_SUPERLU = conf.options.superluLibDir
        conf.env.STLIBPATH_SUPERLU = conf.options.superluLibDir.split(',')
    else:
        conf.env.STLIBPATH_SUPERLU = conf.options.gkylDepsDir+'/superlu/lib'

    conf.env.STLIB_SUPERLU = ["superlu"]

    if conf.options.superluLinkLibs:
        libList = conf.options.superluLinkLibs
        conf.env.append_value('STLIB_SUPERLU', libList.split(','))
         
    conf.start_msg('Checking for SUPERLU')
    conf.check(header_name='slu_ddefs.h', features='cxx cxxprogram', use="SUPERLU", mandatory=True)
    conf.end_msg("Found SUPERLU")
    conf.env['SUPERLU_FOUND'] = True
    return 1

def detect(conf):
    return detect_superlu(conf)
