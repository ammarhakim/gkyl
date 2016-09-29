"""
Detect LuaJIT includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--luajit-inc-dir', type='string', help='Path to LUAJIT includes', dest='luaJitIncDir')
    opt.add_option('--luajit-lib-dir', type='string', help='Path to LUAJIT libraries', dest='luaJitLibDir')

@conf
def check_luajit(conf):
    opt = conf.options
    conf.env['LIB_LUAJIT'] = ''
    conf.env['LUAJIT_FOUND'] = False
    
    if conf.options.luaJitIncDir:
	conf.env.INCLUDES_LUAJIT = conf.options.luaJitIncDir
    else:
        conf.env.INCLUDES_LUAJIT = [conf.options.gkylDepDir+'/luajit/include/luajit-2.1']
        
    if conf.options.luaJitLibDir:
	conf.env.LIBPATH_LUAJIT = conf.options.luaJitLibDir
    else:
        conf.env.STLIBPATH_LUAJIT = [conf.options.gkylDepDir+'/luajit/lib']

    conf.env.STLIB_LUAJIT = ["luajit-5.1"]
        
    conf.start_msg('Checking for LUAJIT include')
    conf.check(header_name='lua.hpp', features='cxx cxxprogram', use="LUAJIT", mandatory=True)
    conf.end_msg("Found LuaJIT")
    return 1

def detect(conf):
    return detect_luajit(conf)
