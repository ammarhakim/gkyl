"""
Detect CUDA compiler, includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf
from waflib import Task
from waflib.TaskGen import extension
from waflib.Tools import ccroot, c_preproc

class cutools(Task.Task):
        run_str = '${NVCC} -c -dc -O3 --shared ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${CPPPATH_ST:INCPATHS} ${DEFINES_ST:DEFINES} ${CXX_SRC_F} ${SRC} -o ${OUT}/${TGT}'
        color   = 'GREEN'
        ext_in  = ['.h']
        vars    = ['CCDEPS']
        scan    = c_preproc.scan
        shell   = False

@extension('.cu', '.cuda')
def c_hook(self, node):
    self.create_compiled_task('cutools', node)
    return self

def options(opt):
    opt.add_option('--cuda-inc-dir', type='string', help='Path to CUTOOLS includes', dest='cuIncDir')
    opt.add_option('--cuda-lib-dir', type='string', help='Path to CUTOOLS libraries', dest='cuLibDir')
    opt.add_option('--cuda-link-libs', type='string', help='List of CUTOOLS libraries to link',
                   dest='cuLinkLibs', default='cudart')

@conf
def check_cutools(conf):
    opt = conf.options
    conf.env['CUTOOLS_FOUND'] = False
    
    if conf.options.cuIncDir:
        conf.env.INCLUDES_CUTOOLS = conf.options.cuIncDir.split(":")

    if conf.options.cuLibDir:
        conf.env.LIBPATH_CUTOOLS = conf.options.cuLibDir.split(":")
    conf.env.OUT = conf.options.out

    libList = conf.options.cuLinkLibs
    conf.env.LIB_CUTOOLS = libList.split(',')
    
    conf.start_msg('Checking for NVCC compiler')
    try:
        conf.find_program('nvcc', var='NVCC', mandatory=True)
        conf.end_msg("Found NVCC")
        conf.check(header_name='cuda.h', features='cxx cxxprogram', use="CUTOOLS", mandatory=True)
        conf.check(header_name='cuda_runtime.h', features='cxx cxxprogram', use="CUTOOLS", mandatory=True)
        conf.end_msg("Linking to libraries work")
        conf.env['CUTOOLS_FOUND'] = True
    except:
        conf.end_msg("Not found NVCC", "YELLOW")
        conf.env['CUTOOLS_FOUND'] = False
        
    return 1

def detect(conf):
    return detect_cutools(conf)
