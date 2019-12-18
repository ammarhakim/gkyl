"""
Detect CUDA compiler, includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf
from waflib import Task
from waflib.TaskGen import extension, feature, after_method, before_method
from waflib.Tools import ccroot, c_preproc

class cuda(Task.Task):
        run_str = '${NVCC} -c -dc -O3 -D USE_CUDA_H=0 ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${CPPPATH_ST:INCPATHS} ${DEFINES_ST:DEFINES} ${CXX_SRC_F} ${SRC} -o ${OUT}/${TGT}'
        color   = 'GREEN'
        ext_in  = ['.h']
        vars    = ['CCDEPS']
        scan    = c_preproc.scan
        shell   = False

class cudacpp(Task.Task):
        run_str = '${NVCC} -x cu -c -dc -O3 -D USE_CUDA_H=0 ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${CPPPATH_ST:INCPATHS} ${DEFINES_ST:DEFINES} ${CXX_SRC_F} ${SRC} -o ${OUT}/${TGT}'
        color   = 'GREEN'
        ext_in  = ['.h']
        vars    = ['CCDEPS']
        scan    = c_preproc.scan
        shell   = False

@extension('.cu', '.cuda')
def hook(self, node):
    self.create_compiled_task('cuda', node)
    return self

# this allows cuda compilation of non *.cu files if features=nvcc is passed
@feature('nvcc')
@before_method('process_source', 'process_rule')
def change_to_cu(self):
    self.source = self.to_nodes(getattr(self, 'source', []))
    for x in self.source:
       self.create_compiled_task('cudacpp', x)

    self.source = []

# enable features=culink, so that linking is done with nvcc
@after_method('process_source')
@feature('culink')
def call_apply_link(self):
    pass

class culink(ccroot.link_task):
    run_str = '${NVCC} -O3 -dlink ${SRC} ${CXXLNK_TGT_F} ${TGT} ${FRAMEWORKPATH_ST:FRAMEWORKPATH} ${FRAMEWORK_ST:FRAMEWORK} ${ARCH_ST:ARCH} ${STLIBPATH_ST:STLIBPATH} ${STLIB_ST:STLIB} ${LIBPATH_ST:LIBPATH} ${LIB_ST:LIB}'
    vars    = ['LINKDEPS']
    ext_out = ['.bin']
    inst_to = '${BINDIR}'

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
