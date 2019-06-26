"""
Detect ZMQ includes and libraries.
"""

import os, glob, types
from waflib.Configure import conf

def options(opt):
    opt.add_option('--zmq-inc-dir', type='string', help='Path to ZMQ includes', dest='zmqIncDir')
    opt.add_option('--zmq-lib-dir', type='string', help='Path to ZMQ libraries', dest='zmqLibDir')

@conf
def check_zmq(conf):
    opt = conf.options
    conf.env['ZMQ_FOUND'] = False

    # include directory
    if conf.options.zmqIncDir:
        conf.env.INCLUDES_ZMQ = conf.options.zmqIncDir
    else:
        conf.env.INCLUDES_ZMQ = [conf.options.gkylDepsDir+'/zeromq/include/']

    # lib directory
    if conf.options.zmqLibDir:
        conf.env.STLIBPATH_ZMQ = conf.options.zmqLibDir
        conf.env.STLIB_ZMQ = ["zmq"]
    else:
        conf.env.STLIBPATH_ZMQ = [conf.options.gkylDepsDir+'/zeromq/lib']
        conf.env.STLIB_ZMQ = ["zmq"]
         
    conf.start_msg('Checking for ZMQ')
    try:
        conf.check(header_name='zmq.h', features='cxx cxxprogram', use='ZMQ', mandatory=True)
        conf.end_msg("Found ZMQ")
        conf.env['ZMQ_FOUND'] = True
    except:
       conf.env['ZMQ_FOUND'] = False
       conf.end_msg("ZMQ not found", "YELLOW")
              
    return 1

def detect(conf):
    return detect_zmq(conf)
