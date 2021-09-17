"""
Detect FFTW includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-fftw', help=('Enable FFTW'),
                   dest='enable_fftw', action='store_true',
                   default=True)
    opt.add_option('--disable-fftw', help=('Disable FFTW'),
                   dest='enable_fftw', action='store_false')
    opt.add_option('--fftw-inc-dir', type='string', help='Path to FFTW includes', dest='fftwIncDir')
    opt.add_option('--fftw-lib-dir', type='string', help='Path to FFTW libraries', dest='fftwLibDir')
    opt.add_option('--fftw-link-libs', type='string', help='List of FFTW libraries to link',
                   dest='fftwLinkLibs', default='fftw3')

@conf
def check_fftw(conf):
    opt = conf.options
    conf.env['FFTW_FOUND'] = False
    
    if not conf.options.enable_fftw:
        return
    if conf.options.fftwIncDir:
        conf.env.INCLUDES_FFTW = conf.options.fftwIncDir
    else:
        conf.env.INCLUDES_FFTW = conf.options.gkylDepsDir#+"/openmpi-3.1.2/include"

    if conf.options.fftwLibDir:
        conf.env.LIBPATH_FFTW = conf.options.fftwLibDir
    else:
        conf.env.LIBPATH_FFTW = conf.options.gkylDepsDir#+"/openmpi-3.1.2/lib"

    libList = conf.options.fftwLinkLibs
    conf.env.LIB_FFTW = libList.split(',')
    
    conf.start_msg('Checking for FFTW')
    conf.check(header_name='fftw3.h', features='cxx cxxprogram', use="FFTW", mandatory=True)
    conf.end_msg("Found FFTW")

    conf.env['FFTW_FOUND'] = True
    return 1

def detect(conf):
    return detect_fftw(conf)
