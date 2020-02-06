"""
Detect SPECTRA includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-spectra', help=('Enable SPECTRA'),
                   dest='enable_spectra', action='store_true',
                   default=True)
    opt.add_option('--disable-spectra', help=('Disable SPECTRA'),
                   dest='enable_spectra', action='store_false')
    opt.add_option('--spectra-inc-dir', type='string', help='Path to SPECTRA includes', dest='spectraIncDir')

@conf
def check_spectra(conf):
    opt = conf.options
    conf.env['SPECTRA_FOUND'] = False

    if not conf.options.enable_spectra:
        return
    
    if conf.options.spectraIncDir:
        conf.env.INCLUDES_SPECTRA = conf.options.spectraIncDir
    else:
        conf.env.INCLUDES_SPECTRA = conf.options.gkylDepsDir+'/spectra/include'

    conf.start_msg('Checking for SPECTRA')
    conf.check(header_name='Spectra/SymEigsSolver.h', features='cxx cxxprogram', use='SPECTRA', mandatory=False)
    conf.end_msg("Found SPECTRA")
    conf.env['SPECTRA_FOUND'] = True
    return 1

def detect(conf):
    return detect_spectra(conf)
