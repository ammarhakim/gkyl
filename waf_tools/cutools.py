"""
Detect CUDA compiler, includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf
from waflib import Task
from waflib.TaskGen import extension, feature, after_method, before_method
from waflib.Tools import ccroot, c_preproc
from waflib.Tools.cxx import cxxprogram

class nvcc(Task.Task):
        run_str = '${NVCC} -x cu -c -dc -O3 -std=c++14 -lineinfo --maxrregcount=255 -arch=sm_70 --compiler-options="-fPIC" ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${CPPPATH_ST:INCPATHS} ${DEFINES_ST:DEFINES} ${CXX_SRC_F} ${SRC} -o ${OUT}/${TGT}'
        color   = 'GREEN'
        ext_in  = ['.h']
        vars    = ['CCDEPS']
        scan    = c_preproc.scan
        shell   = False

@extension('.cu', '.cuda')
def hook(self, node):
    self.create_compiled_task('nvcc', node)
    return self

# this allows cuda compilation of non *.cu files if features=nvcc is passed
@feature('nvcc')
@before_method('process_source')
def change_to_cu(self):
    self.source = self.to_nodes(getattr(self, 'source', []))
    for x in self.source:
       self.create_compiled_task('nvcc', x)

    self.source = []

# enable features=cushlib, so that cu files linked into shared library
@before_method('process_rule')
@feature('cushlib')
def update_pattern_cushlib(self):
    self.env['cushlib_PATTERN'] = self.env['cshlib_PATTERN']

class cushlib(ccroot.link_task):
    run_str = '${NVCC} -shared -arch=sm_70 ${CXXLNK_SRC_F} ${SRC} ${CXXLNK_TGT_F} ${TGT[0].abspath()}'
    inst_to = '${LIBDIR}'

# enable features=culink, so that linking is done with nvcc -dlink
@before_method('process_rule')
@feature('culink')
def update_pattern_culink(self):
    self.env['culink_PATTERN'] = "%s.o"

class culink(ccroot.link_task):
    run_str = '${NVCC} -dlink -arch=sm_70 ${CXXLNK_SRC_F} ${SRC} ${CXXLNK_TGT_F} ${TGT[0].abspath()}'

def options(opt):
    opt.add_option('--enable-cuda', help=('Enable CUDA GPU build (if NVCC is found)'),
                   dest='enable_cuda', action='store_true',
                   default=False)
    opt.add_option('--disable-cuda', help=('Disable CUDA GPU build (even if NVCC can be found)'),
                   dest='enable_cuda', action='store_false')
    opt.add_option('--cuda-inc-dir', type='string', help='Path to CUTOOLS includes', dest='cuIncDir')
    opt.add_option('--cuda-lib-dir', type='string', help='Path to CUTOOLS libraries', dest='cuLibDir')
    opt.add_option('--cuda-link-libs', type='string', help='List of CUTOOLS libraries to link',
                   dest='cuLinkLibs', default='cudart')

@conf
def check_cutools(conf):
    opt = conf.options
    conf.env['CUTOOLS_FOUND'] = False
    
    if not conf.options.enable_cuda:
        return
    
    conf.start_msg('Checking for NVCC compiler')
    try:
        conf.find_program('nvcc', var='NVCC', mandatory=True)
        conf.end_msg("Found NVCC")
        if conf.options.cuIncDir:
            conf.env.INCLUDES_CUTOOLS = conf.options.cuIncDir.split(":")

        if conf.options.cuLibDir:
            conf.env.LIBPATH_CUTOOLS = conf.options.cuLibDir.split(":")
        conf.env.OUT = conf.options.out

        libList = conf.options.cuLinkLibs
        conf.env.LIB_CUTOOLS = libList.split(',')
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
