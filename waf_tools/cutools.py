"""
Detect CUDA compiler, includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--cutools-inc-dir', type='string', help='Path to CUTOOLS includes', dest='cutoolsIncDir')
    opt.add_option('--cutools-lib-dir', type='string', help='Path to CUTOOLS libraries', dest='cutoolsLibDir')
    opt.add_option('--cutools-link-libs', type='string', help='List of CUTOOLS libraries to link',
                   dest='cutoolsLinkLibs', default='cutools')

@conf
def check_cutools(conf):
    opt = conf.options
    conf.env['CUTOOLS_FOUND'] = False
    
    if conf.options.cutoolsIncDir:
        conf.env.INCLUDES_CUTOOLS = conf.options.cutoolsIncDir

    if conf.options.cutoolsLibDir:
        conf.env.LIBPATH_CUTOOLS = conf.options.cutoolsLibDir

    libList = conf.options.cutoolsLinkLibs
    conf.env.LIB_CUTOOLS = libList.split(',')
    
    conf.start_msg('Checking for NVCC compiler')
    try:
        conf.find_program('nvcc', var='nvcc', mandatory=True)
        conf.end_msg("Found NVCC")
        conf.env['CUTOOLS_FOUND'] = True
    except:
        conf.end_msg("Not found NVCC", "YELLOW")
        conf.env['CUTOOLS_FOUND'] = False
        
    return 1

def detect(conf):
    return detect_cutools(conf)
