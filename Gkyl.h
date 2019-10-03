// Gkyl ------------------------------------------------------------------------
//
// Class representing the top Gkeyll application (executable)
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

// std include
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <utility>
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

#include <lua.hpp>

// Gkyl/Lib includes
#include <lfs.h>
#include <whereami.h>
#include <base64.h>

// Find location of executable
std::string
findExecPath() {
  auto len = wai_getExecutablePath(NULL, 0, NULL);
  char *path = (char*) malloc(len+1);
  int dirname_len; wai_getExecutablePath(path, len, &dirname_len);
  path[len] = '\0';
  std::string execPath(path); // = std::string(path, dirname_len);
  free(path);
  return execPath;
}

// Read contents of input file. This uses MPI broadcast to get file
// across processors so as not to hit the file system with huge number
// of file open calls at once.
std::string readInputFile(const std::string& inpFile) {
  int size;
  std::string buffer;

  // read contents on rank 0
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
  if (rank == 0) {
    std::ifstream inpf(inpFile.c_str()); 
    inpf.seekg(0, std::ios::end);
    size = inpf.tellg();
    buffer = std::string(size, ' ');
    inpf.seekg(0);
    inpf.read(&buffer[0], size);
  }

  // send to all other ranks
  MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank > 0)
    buffer = std::string(size, ' ');
  MPI_Bcast(&buffer[0], size, MPI_CHAR, 0, MPI_COMM_WORLD);

  return buffer;
}

/**
 * Top-level object representing complete simulation.
 */
class Gkyl {
  public:
    /**
     * @param luaExpr Lua expression to execute before running input-file or tool
     * @param inpFile Name of input file
     * @param args    Commands/arguments to feed to app or tool
     */
    Gkyl(const std::string& luaExpr, const std::string& inpFile, const std::list<std::string>& args);

    /** Show list of resitered tools */
    void showToolList();

    /** Run simulation */
    int run();
    
  private:
    /* True if we have input file */    
    bool hasInpFile;
    /* Lua expression to execute */
    std::string luaExpr;    
    /* Name of input file */
    std::string inpFile;
    /* List of app/tool arguments */
    std::list<std::string> args;
    /* List of tools */
    std::map<std::string, std::pair<std::string, std::string> > toolList;
    /* Path to executable */
    std::string execPath;
    /* Input file contents */
    std::string inpFileContent;
};

Gkyl::Gkyl(const std::string& luaExpr, const std::string& inpFile, const std::list<std::string>& args)
  : hasInpFile(true), luaExpr(luaExpr), inpFile(inpFile), args(args), execPath(findExecPath())
{
  toolList = {
    { "help", { "help.lua", "    Gkeyll help system" } },
    { "examples", {"examples.lua", "Example input files"} },
    { "queryrdb", {"queryrdb.lua", "Query/modify regression test DB"} }
  };

  if (inpFile.empty())
    hasInpFile = false;
  else
    inpFileContent = readInputFile(inpFile);
}

void Gkyl::showToolList() {
  std::cout << "Following tools are available. Query tool help for more information." << std::endl;
  for (auto t : toolList)
    std::cout << " " << t.first << "  " << t.second.second << std::endl;
}

int Gkyl::run() {
  
  return 0;
}
