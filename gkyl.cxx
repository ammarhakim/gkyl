// Gkyl ------------------------------------------------------------------------
//
// Top-level entry point into Gkyl
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <gkylconfig.h>
#include <gkylhgtip.h>

// library includes
#include <lua.hpp>

// std include
#include <fenv.h> 
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdlib.h>
#include <string>

#if defined(__clang__)
// nothing to include
#elif defined(__powerpc__)
// nothing to include
#elif defined(__GNUC__) || defined(__GNUG__)
# include <xmmintrin.h>
#endif

# include <mpi.h>
# include <GkMpiFuncs.h>

#ifdef HAVE_ADIOS_H
# include <adios.h>
# include <adios_read.h>
#endif

#ifdef HAVE_EIGEN_CORE
#include <Eigen/Core>
#endif

// Gkyl includes
#include <lfs.h>
#include <whereami.h>
#include <base64.h>

// A simple logger for parallel simulations
class Logger {
  public:
    Logger() : rank(0) {
#ifdef HAVE_MPI_H
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif      
    }
    void log(const char *m) const {
      if (rank == 0)
        std::cout << m << std::endl;
    }
  private:
    int rank;
};

// Class to store information about various tools
class GkToolInfo {
  public:
    GkToolInfo(const std::string& script, const std::string& description)
      : script(script), description(description) {
    }

    // Name of script and short description
    std::string script, description;
};

// Finish simulation
int finish(int err) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (err != 0) {
    MPI_Abort(MPI_COMM_WORLD, err);
  }
  else {
    MPI_Finalize();
  }

  return err;
}

// Find location of executable
std::string
findExecPath() {
  int len = wai_getExecutablePath(NULL, 0, NULL);
  char *path = (char*) malloc(len+1);
  int dirname_len;
  wai_getExecutablePath(path, len, &dirname_len);
  path[len] = '\0';
  std::string execPath = std::string(path, dirname_len);
  free(path);
  return execPath;
}

// Create top-level variable definitions
std::string
createTopLevelDefs(int argc, char **argv, const std::string& inpFile,
  const std::map<std::string, GkToolInfo>& toolList) {
  // load compile-time constants into Lua compiler so they become
  // available to scripts
  std::ostringstream varDefs;

  // find executable location and modify package paths
  auto execPath = findExecPath();
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
  varDefs << "GKYL_HG_CHANGESET = \"" << GKYL_HG_CHANGESET << "\"" << std::endl;
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
  
  // numeric limits
  varDefs << "GKYL_MIN_DOUBLE = " << std::numeric_limits<double>::min() << std::endl;
  varDefs << "GKYL_MIN_FLOAT = " << std::numeric_limits<float>::min() << std::endl;
  varDefs << "GKYL_MAX_DOUBLE = " << std::numeric_limits<double>::max() << std::endl;
  varDefs << "GKYL_MAX_FLOAT = " << std::numeric_limits<float>::max() << std::endl;
  varDefs << "GKYL_EPSILON = " << std::numeric_limits<double>::epsilon() << std::endl;
  
  // set some JIT parameters to fiddle around with optimizations
  varDefs << "jit.opt.start('callunroll=40', 'loopunroll=80', 'maxmcode=40960', 'maxtrace=8000', 'maxrecord=16000', 'minstitch=3')"
          << std::endl;

  // This is change from commit e2448fa. Rerolling it back because of unforseen consequences
  // std::string myInpFile(inpFile);

  // // figure out output prefix. This is complicated by the fact the
  // // input file (or tool file) may be specified via full path
  // auto const slashPos = myInpFile.find_last_of('/');
  // if (std::string::npos != slashPos)
  //   myInpFile = myInpFile.substr(slashPos+1);

  // std::string snm(myInpFile);
  // auto const trunc = myInpFile.find_last_of(".", snm.size());
  // if (std::string::npos != trunc)
  std::string snm(inpFile);
  auto const trunc = inpFile.find_last_of(".", snm.size());
  if (std::string::npos != trunc)
    snm.erase(trunc, snm.size());
  varDefs << "GKYL_OUT_PREFIX = '" << snm << "'" << std::endl;

  // read complete input file and make it available to Lua
  std::ifstream inpf(inpFile.c_str()); 
  inpf.seekg(0, std::ios::end);
  size_t size = inpf.tellg();
  std::string buffer(size, ' ');
  inpf.seekg(0);
  inpf.read(&buffer[0], size);

  // we need to base64 encode file contents to avoid issues with Lua
  // loadstr method getting confused with embedded strings etc
  std::string inpfEncoded = base64_encode(
    reinterpret_cast<const unsigned char*>(buffer.c_str()), buffer.length());
  varDefs << "GKYL_INP_FILE_CONTENTS = \"" << inpfEncoded  << "\" " << std::endl;

  varDefs << "GKYL_COMMANDS = {}" << std::endl;
  // just append list of commands 
  for (unsigned i=2; i<argc; ++i)
    varDefs << "GKYL_COMMANDS[" << i-1 << "] = \"" << argv[i]  << "\"" << std::endl;

  varDefs << "GKYL_TOOLS = {}" << std::endl;
  // append list of tools
  for (auto tool : toolList)
    varDefs << "GKYL_TOOLS." << tool.first << " = { \""
            << tool.second.script  << "\", \"" << tool.second.description
            << "\" }" << std::endl;

  // std::cout << varDefs.str() << std::endl;
    
  return varDefs.str();
}

// These dummy calls are a hack to force the linker to pull in symbols
// from static libraries. There is a better way to do this, I just do
// not know it yet.
#ifdef HAVE_ADIOS_H
int _adios_init(const char *cf, MPI_Comm comm) { return adios_init(cf, comm); }
ADIOS_FILE*
_adios_read_open_file(const char *fname, enum ADIOS_READ_METHOD method, MPI_Comm comm)
{ return adios_read_open_file(fname, method, comm); }
#endif

void
showUsage(const std::map<std::string, GkToolInfo>& tlist) {
  Logger logger;
  
  std::cout << "Usage: gkyl <tool> ..." << std::endl;
  std::cout << "Usage: gkyl LUA-SCRIPT <cmd> ..." << std::endl;
  std::cout << std::endl;
  std::cout << "Supported tools are" << std::endl;
  for (auto tool : tlist)
    std::cout << " " << tool.first << ": " << tool.second.description << std::endl;
}

int
main(int argc, char **argv) {
  // This prevents denormalized floats from occuring in
  // code. Otherwise, the code become horribly slow in some rare (but
  // not impossible to reproduce) situations.
#if defined(__clang__)
  fesetenv(FE_DFL_DISABLE_SSE_DENORMS_ENV);
#elif defined(__powerpc__)
  // not sure what the POWER calls are for denormalized floats
#elif defined(__GNUC__) || defined(__GNUG__)
  _MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);  
#endif

  MPI_Init(&argc, &argv);
  
  Logger logger;

  // list of tools
  std::map<std::string, GkToolInfo> toolList;
  // function to make inserting into list easier
  auto insertInfo = [](std::map<std::string, GkToolInfo>& tl,
    const std::string& nm, const std::string& script, const std::string& descr ) {
    tl.insert(
      std::pair<std::string, GkToolInfo>(nm, GkToolInfo(script, descr)));
  };
  // register various tools  
  insertInfo(toolList, "help", "help.lua", "Gkeyll help system");
  insertInfo(toolList, "examples", "examples.lua", "Example input files");
  insertInfo(toolList, "runtests", "runtests.lua", "Run unit/regression tests");

  if (argc < 2) {
    showUsage(toolList);
    return finish(0);
  }

  std::string inpFile(argv[1]);
  // check if we should use a tool
  auto tool = toolList.find(argv[1]);
  if (toolList.end() != tool) {
    auto execPath = findExecPath(); // path of tool is wrt to executable location
    inpFile =  execPath + "/Tool/" + tool->second.script;
  }

  // check if file exists  
  std::ifstream f(inpFile.c_str());
  if (!f.good()) {
    std::cerr << "Unable to open input file '" << inpFile << "'" << std::endl;
    showUsage(toolList);
    return finish(1);
  }
  f.close();

  // initialize LuaJIT and load libraries
  lua_State *L = luaL_newstate();
  if (NULL == L) {
    // NOTE: we need to use cerr and not 'logger' when something goes
    // wrong in top-level executable. Otherwise error message won't
    // appear anywhere.
    std::cerr << "Unable to create a new LuaJIT interpreter state. Quitting" << std::endl;
    return finish(1);
  }
  lua_gc(L, LUA_GCSTOP, 0);  // stop GC during initialization
  luaL_openlibs(L);  // open standard libraries
  luaopen_lfs(L); // open lua file-system library
  lua_gc(L, LUA_GCRESTART, -1); // restart GC
  
  std::string topDefs = createTopLevelDefs(argc, argv, inpFile, toolList);

  // load variable definitions etc
  if (luaL_loadstring(L, topDefs.c_str()) || lua_pcall(L, 0, LUA_MULTRET, 0)) {
    // some error occured
    const char* ret = lua_tostring(L, -1);
    std::cerr << "*** LOAD ERROR *** " << ret << std::endl;
    lua_close(L);
    return finish(1);    
  }

  // load and run input file
  if (luaL_loadfile(L, inpFile.c_str()) || lua_pcall(L, 0, LUA_MULTRET, 0)) {
    // some error occured
    const char* ret = lua_tostring(L, -1);
    std::cerr << "*** ERROR *** " << ret << std::endl;
    lua_close(L);
    return finish(1);
  }
  lua_close(L);

  return finish(0);
}
