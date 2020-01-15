// Gkyl ------------------------------------------------------------------------
//
// Class representing the top Gkeyll executable
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
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

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

#ifdef HAVE_CUDA_H
#include <cuda.h>
#include <cuda_runtime.h>
#endif

#include <lua.hpp>

// Gkyl/Lib includes
#include <lfs.h>
#include <whereami.h>
#include <base64.h>
#include <gkylgittip.h>

// Find location of executable
std::string
findExecPath() {
  int len = wai_getExecutablePath(NULL, 0, NULL);
  char *path = (char*) malloc(len+1);
  int dirname_len; wai_getExecutablePath(path, len, &dirname_len);
  path[len] = '\0';
  std::string execPath(path, dirname_len);
  free(path);
  return execPath;
}

// Read contents of input file. This uses MPI broadcast to get file
// across processors so that in parallel we do not to hit system with
// huge number of file open calls at once.
std::string readInputFile(const std::string& inpFile) {
  // check if file exists
  int fileExists = 1;
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    std::ifstream inpf(inpFile.c_str());
    if (!inpf.good())
      fileExists = 0;
  }
  // send information about file existence to all other ranks
  MPI_Bcast(&fileExists, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (fileExists == 0)
    // all ranks must throw exception, or the code will hang
    throw std::runtime_error("File " + inpFile + " does not exist");

  int size;  
  std::string buffer;  
  // read contents on rank 0  
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
    /** dtor */
    ~Gkyl();
    /** Run simulation */
    int run();
    /** Return list of registered tools */
    std::vector<std::pair<std::string, std::string>> getToolList() const;

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
    std::map<std::string, std::pair<std::string, std::string>> toolList;
    /* Path to executable */
    std::string execPath;
    /* Input file contents */
    std::string inpFileContent;
    /* Lua state object */
    lua_State *L;

    /* Create string with top-level Lua definitions */
    std::string createTopLevelDefs() const;
    /* Run Lua code in string */
    void runLua(const std::string& lua);
};

Gkyl::Gkyl(const std::string& luaExpr, const std::string& inpFileNm, const std::list<std::string>& args)
  : hasInpFile(!inpFileNm.empty()), luaExpr(luaExpr), inpFile(inpFileNm), args(args), execPath(findExecPath())
{
  toolList = {
    { "help", { "help.lua", "Gkeyll help system" } },
    { "examples", {"examples.lua", "Example input files"} },
    { "queryrdb", {"queryrdb.lua", "Query/modify regression test DB"} },
    { "runregression", {"runregression.lua", "Run regression/unit tests"} }
  };

  if (hasInpFile) {
    // check if tool specified
    auto tool = toolList.find(inpFile);
    if (toolList.end() != tool) 
      inpFile = execPath + "/Tool/" + tool->second.first;
    inpFileContent = readInputFile(inpFile);
  }

  // initialize Lua (I want to replace this with sol2 library)
  L =  luaL_newstate();
  if (NULL == L)
    throw std::runtime_error("Unable to create a new Lua interpreter state.");
  lua_gc(L, LUA_GCSTOP, 0);  // stop GC during initialization
  luaL_openlibs(L);  // open standard libraries
  luaopen_lfs(L); // open lua file-system library
  lua_gc(L, LUA_GCRESTART, -1); // restart GC

  // load variable definitions into Lua
  runLua(createTopLevelDefs());
}

Gkyl::~Gkyl() {
  lua_close(L);
}

int Gkyl::run() {
  runLua(luaExpr);
  if (hasInpFile)
    runLua(inpFileContent);
  return 0;
}

std::vector<std::pair<std::string, std::string> >
Gkyl::getToolList() const {
  std::vector<std::pair<std::string, std::string>> tlist;
  for (auto t : toolList) tlist.push_back({t.first, t.second.second});
  return tlist;
}

std::string Gkyl::createTopLevelDefs() const {
  std::ostringstream varDefs;

  // find executable location and modify package paths
  varDefs << "package.path = package.path .. \";"
          << execPath << "/?.lua;"
          << execPath << "/Lib/?.lua;" // we need add Lib to allow using external libraries
          << execPath << "/?/init.lua" << "\"" << std::endl;

  varDefs << "package.cpath = package.cpath .. \";"
          << execPath << "/../lib/lib?.so;"
          << execPath << "/../lib/lib?.dylib;"
          << "\"" << std::endl;

  // info about build
  varDefs << "GKYL_EXEC_PATH = \"" << execPath << "\"" << std::endl;
  varDefs << "GKYL_EXEC = \"" << execPath << "/gkyl\"" << std::endl;
  varDefs << "GKYL_GIT_CHANGESET = \"" << GKYL_GIT_CHANGESET << "\"" << std::endl;
  varDefs << "GKYL_BUILD_DATE = \"" << __DATE__ << " " << __TIME__ << "\"" << std::endl;
  
#ifdef HAVE_MPI_H
  varDefs << "GKYL_HAVE_MPI = true" << std::endl;
#else
  varDefs << "GKYL_HAVE_MPI = false" << std::endl;
#endif

#ifdef HAVE_ADIOS_H
  varDefs << "GKYL_HAVE_ADIOS = true" << std::endl;
#else
  varDefs << "GKYL_HAVE_ADIOS = false" << std::endl;
#endif

#ifdef HAVE_EIGEN_CORE
  varDefs << "GKYL_HAVE_EIGEN = true" << std::endl;
#else
  varDefs << "GKYL_HAVE_EIGEN = false" << std::endl;
#endif

#ifdef HAVE_ZMQ_H
  varDefs << "GKYL_HAVE_ZMQ = true" << std::endl;
  int zmajor, zminor, zpatch;
  zmq_version (&zmajor, &zminor, &zpatch);
  varDefs << "GKYL_ZMQ_VERSION = { ";
  varDefs << "major = " << zmajor << ", ";
  varDefs << "minor = " << zminor << ", ";
  varDefs << "patch = " << zpatch << ", ";
  varDefs << "}" << std::endl;
#else
  varDefs << "GKYL_HAVE_ZMQ = false" << std::endl;
#endif

#ifdef USING_SQLITE3
  varDefs << "GKYL_HAVE_SQLITE3 = true" << std::endl;
#else
  varDefs << "GKYL_HAVE_SQLITE3 = false" << std::endl;
#endif

#ifdef HAVE_CUDA_H
  varDefs << "GKYL_HAVE_CUDA = true" << std::endl;
  varDefs << "GKYL_DEFAULT_NUM_THREADS = 256" << std::endl;  
  int cuDriverVersion;
  cudaDriverGetVersion(&cuDriverVersion);
  varDefs << "GKYL_CUDA_DRIVER_VERSION = " << cuDriverVersion << std::endl;
#else
  varDefs << "GKYL_HAVE_CUDA = false" << std::endl;
#endif
  
  // numeric limits
  varDefs << "GKYL_MIN_DOUBLE = " << std::numeric_limits<double>::min() << std::endl;
  varDefs << "GKYL_MIN_FLOAT = " << std::numeric_limits<float>::min() << std::endl;
  varDefs << "GKYL_MAX_DOUBLE = " << std::numeric_limits<double>::max() << std::endl;
  varDefs << "GKYL_MAX_FLOAT = " << std::numeric_limits<float>::max() << std::endl;
  varDefs << "GKYL_EPSILON = " << std::numeric_limits<double>::epsilon() << std::endl;

  // set some JIT parameters to fiddle around with optimizations
  varDefs << "jit.opt.start('callunroll=40', 'loopunroll=80', 'maxmcode=40960', 'maxtrace=8000', 'maxrecord=16000', 'minstitch=3')"
          << std::endl;

  // output prefix
  std::string snm(inpFile);
  auto const trunc = inpFile.find_last_of(".", snm.size());
  if (std::string::npos != trunc)
    snm.erase(trunc, snm.size());
  varDefs << "GKYL_OUT_PREFIX = '" << snm << "'" << std::endl;

  // we need to base64 encode file contents to avoid issues with Lua
  // loadstr method getting confused with embedded strings etc
  std::string inpfEncoded = base64_encode(
    reinterpret_cast<const unsigned char*>(inpFileContent.c_str()), inpFileContent.length());

  varDefs << "GKYL_INP_FILE_CONTENTS = \"" << inpfEncoded  << "\" " << std::endl;

  // append list of unparsed command-line parameters so app/tool has
  // access to them
  varDefs << "GKYL_COMMANDS = {}" << std::endl;
  do {
    int i = 1;
    for (auto a : args)
      varDefs << "GKYL_COMMANDS[" << i++ << "] = \"" << a << "\"" << std::endl;
  } while (0);

  varDefs << "GKYL_TOOLS = {}" << std::endl;
  // append list of tools
  for (auto tool : toolList)
    varDefs << "GKYL_TOOLS." << tool.first << " = { \""
            << tool.second.first  << "\", \"" << tool.second.second
            << "\" }" << std::endl;

  return varDefs.str();
}

void Gkyl::runLua(const std::string& lua) {
  if (luaL_loadstring(L, lua.c_str()) || lua_pcall(L, 0, LUA_MULTRET, 0)) {
    // some error occured
    const char* ret = lua_tostring(L, -1);
    throw std::runtime_error(std::string("*** LOAD ERROR ***\n ").append(ret));
  }  
}
