"""
Detect NCCL includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-nccl', help=('Enable NCCL build'),
                   dest='enable_nccl', action='store_true',
                   default=False)
    opt.add_option('--disable-nccl', help=('Disable NCCL build'),
                   dest='enable_nccl', action='store_false')
    opt.add_option('--nccl-inc-dir', type='string', help='Path to NCCL includes', dest='ncclIncDir')
    opt.add_option('--nccl-lib-dir', type='string', help='Path to NCCL libraries', dest='ncclLibDir')
    opt.add_option('--nccl-link-libs', type='string', help='List of NCCL libraries to link',
                   dest='ncclLinkLibs', default='nccl')

@conf
def check_nccl(conf):
    opt = conf.options
    conf.env['NCCL_FOUND'] = False
    
    if not conf.options.enable_nccl:
        return
    if conf.options.ncclIncDir:
        conf.env.INCLUDES_NCCL = conf.options.ncclIncDir
    else:
        conf.env.INCLUDES_NCCL = conf.options.gkylDepsDir+"/nccl/include"

    if conf.options.ncclLibDir:
        conf.env.LIBPATH_NCCL = conf.options.ncclLibDir
    else:
        conf.env.LIBPATH_NCCL = conf.options.gkylDepsDir+"/nccl/lib"

    libList = conf.options.ncclLinkLibs
    conf.env.LIB_NCCL = libList.split(',')
    
    conf.start_msg('Checking for NCCL')
    conf.check(header_name='nccl.h', features='cxx cxxprogram', use="NCCL", mandatory=True)
    conf.end_msg("Found NCCL")

    conf.env['NCCL_FOUND'] = True
    return 1

def detect(conf):
    return detect_nccl(conf)
