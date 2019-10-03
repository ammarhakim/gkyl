// Gkyl ------------------------------------------------------------------------
//
// Class representing the top Gkeyll application (executable)
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

// library includes
#include <lua.hpp>

// std include
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <tuple>

#ifdef HAVE_MPI_H
# include <mpi.h>
# include <GkMpiFuncs.h>
#endif

#ifdef HAVE_ADIOS_H
# include <adios.h>
# include <adios_read.h>
#endif

#ifdef HAVE_EIGEN_CORE
#include <Eigen/Core>
#endif

#ifdef HAVE_ZMQ_H
#include <zmq.h>
#endif

// Gkyl/Lib includes
#include <lfs.h>
#include <whereami.h>
#include <base64.h>

/**
 * Top-level object representing the complete simulation.
 */
class Gkyl {
  public:
    /**
     * @param luaExpr Lua expression to execute before running input-file or tool
     * @param inpFile Name of input file
     * @param args    Commands/arguments to feed to app or tool
     */
    Gkyl(const std::string& luaExpr, const std::string& inpFile, const std::list<std::string>& args);

    /** Run simulation */
    int run();
    
  private:
    bool hasInpFile;
};

Gkyl::Gkyl(const std::string& luaExpr, const std::string& inpFile, const std::list<std::string>& args) {
  
}

int Gkyl::run() {
  return 0;
}
