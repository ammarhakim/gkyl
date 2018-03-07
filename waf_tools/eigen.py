"""
Detect EIGEN includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--disable-eigen', help=('Disable EIGEN'),
                   dest='disable_eigen', action='store_true',
                   default=False)
    opt.add_option('--eigen-inc-dir', type='string', help='Path to EIGEN includes', dest='eigenIncDir')

@conf
def check_eigen(conf):
    opt = conf.options
    conf.env['EIGEN_FOUND'] = False

    if conf.options.disable_eigen:
        return
    
    if conf.options.eigenIncDir:
        conf.env.INCLUDES_EIGEN = conf.options.eigenIncDir
    else:
        conf.env.INCLUDES_EIGEN = conf.options.gkylDepsDir+'/eigen-3.3.4/include/eigen3'

    #conf.env.append_value('CXXFLAGS', '-D_NOMPI')
    #conf.env.append_value('CFLAGS', '-D_NOMPI')
    #conf.env.STLIB_EIGEN = ["eigen"]

    conf.start_msg('Checking for EIGEN')
    conf.check(header_name='Eigen/Core', features='cxx cxxprogram', use='EIGEN', mandatory=True)
    conf.end_msg("Found EIGEN")
    conf.env['EIGEN_FOUND'] = True
    return 1

def detect(conf):
    return detect_eigen(conf)
