"""
Detect ADIOS includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--enable-adios', help=('Enable ADIOS'),
                   dest='enable_adios', action='store_true',
                   default=True)
    opt.add_option('--adios-inc-dir', type='string', help='Path to ADIOS includes', dest='adiosIncDir')
    opt.add_option('--adios-lib-dir', type='string', help='Path to ADIOS libraries', dest='adiosLibDir')
    opt.add_option('--adios-link-libs', type='string', help='List of libraries to link with ADIOS',
                   dest='adiosLinkLibs')

@conf
def check_adios(conf):
    opt = conf.options
    conf.env['ADIOS_FOUND'] = False

    if not conf.options.enable_adios:
        return
    
    if conf.options.adiosIncDir:
        conf.env.INCLUDES_ADIOS = conf.options.adiosIncDir
    else:
        conf.env.INCLUDES_ADIOS = conf.options.gkylDepsDir+'/adios2/include'

    if conf.options.adiosLibDir:
        #conf.env.LIBPATH_ADIOS = conf.options.adiosLibDir
        conf.env.LIBPATH_ADIOS = conf.options.adiosLibDir.split(',')
    else:
        conf.env.LIBPATH_ADIOS = conf.options.gkylDepsDir+'/adios2/lib'

    if conf.options.enable_mpi:
        conf.env.LIB_ADIOS = ["adios2_c", "adios2_c_mpi"]
        conf.env.append_value('CFLAGS', '-DADIOS2_USE_MPI -isystem')
        conf.env.append_value('CXXFLAGS', '-DADIOS2_USE_MPI')
    else:
        conf.env.append_value('CXXFLAGS', '-D_NOMPI')
        conf.env.append_value('CFLAGS', '-D_NOMPI')
        conf.env.LIB_ADIOS = ["adios_nompi"]

    if conf.options.adiosLinkLibs:
        libList = conf.options.adiosLinkLibs
        conf.env.append_value('LIB_ADIOS', libList.split(','))
         
    conf.start_msg('Checking for ADIOS')
    conf.check(header_name='adios2_c.h', features='cxx cxxprogram', use="ADIOS MPI", mandatory=True)
    conf.end_msg("Found ADIOS")
    conf.env['ADIOS_FOUND'] = True
    return 1

def detect(conf):
    return detect_adios(conf)
