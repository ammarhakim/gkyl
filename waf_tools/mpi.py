"""
Detect MPI includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-mpi', help=('Enable parallel build'),
                   dest='enable_mpi', action='store_true',
                   default=False)
    opt.add_option('--mpi-inc-dir', type='string', help='Path to MPI includes', dest='mpiIncDir')
    opt.add_option('--mpi-lib-dir', type='string', help='Path to MPI libraries', dest='mpiLibDir')

@conf
def check_mpi(conf):
    opt = conf.options
    conf.env['MPI_FOUND'] = False
    
    if not conf.options.enable_mpi:
	return
    if conf.options.mpiIncDir:
	conf.env.INCLUDES_MPI = conf.options.mpiIncDir
    else:
        conf.fatal("Please specify MPI include directories by --mpi-inc-dir")
    if conf.options.mpiLibDir:
	conf.env.LIBPATH_MPI = conf.options.mpiLibDir
        conf.env.LIB_MPI = ["mpi", "mpi_cxx"]
    else:
        conf.fatal("Please specify MPI library directories by --mpi-lib-dir")
        
    conf.start_msg('Checking for MPI')
    conf.check(header_name='mpi.h', features='cxx cxxprogram', use="MPI", mandatory=True)
    conf.end_msg("Found MPI")

    conf.env['MPI_FOUND'] = True
    return 1

def detect(conf):
    return detect_mpi(conf)
