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
#include <sstream>
#include <stdlib.h>

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
createTopLevelDefs(int argc, char **argv) {
  // load compile-time constants into Lua compiler so they become
  // available to scripts
  std::ostringstream varDefs;

  // find executable location and modify package paths
  std::string execPath = findExecPath();
  varDefs << "package.path = package.path .. \";"
          << execPath << "/?.lua;"
          << execPath << "/Lib/?.lua;" // we need add Lib to allow using external libraries
          << execPath << "/?/init.lua" << "\"" << std::endl;

  varDefs << "package.cpath = package.cpath .. \";"
          << execPath << "/../lib/lib?.so;"
          << execPath << "/../lib/lib?.dylib;"
          << "\"" << std::endl;

  // info about build  
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
  varDefs << "jit.opt.start('callunroll=20', 'loopunroll=60')" << std::endl;

  std::string inpFile(argv[1]);
  std::string snm(argv[1]);
  unsigned trunc = inpFile.find_last_of(".", snm.size());
  if (trunc > 0)
    snm.erase(trunc, snm.size());
  varDefs << "GKYL_OUT_PREFIX = '" << snm << "'" << std::endl;

  varDefs << "GKYL_COMMANDS = {}" << std::endl;
  // just append list of commands 
  for (unsigned i=2; i<argc; ++i)
    varDefs << "GKYL_COMMANDS[" << i-1 << "] = \"" << argv[i]  << "\"" << std::endl;

  //std::cout << varDefs.str() << std::endl;
    
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
showUsage() {
  Logger logger;
  
  logger.log("Usage: gkyl LUA-SCRIPT <cmd> ...");
  logger.log(" List of <cmd> is optional");
  logger.log("See particular App documentation for supported commands.\n");
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

  if (argc < 2) {
    showUsage();
    return finish(0);
  }

  // check if file exists  
  std::string inpFile(argv[1]);
  std::ifstream f(inpFile.c_str());
  if (!f.good()) {
    std::cerr << "Unable to open input file '" << inpFile << "'" << std::endl;
    showUsage();
    return finish(1);
  }
  f.close();

  // initialize LuaJIT and load libraries
  lua_State *L = luaL_newstate();
  if (L==NULL) {
    // NOTE: we need to use cerr and not 'logger' when something goes
    // wrong in top-level executable. Otherwise error message won't
    // appear anywhere.
    std::cerr << "Unable to create a new LuaJIT interpreter state. Quitting" << std::endl;
    return finish(1);
  }
  lua_gc(L, LUA_GCSTOP, 0);  // stop GC during initialization
  luaL_openlibs(L);  // open standard libraries
  luaopen_lfs(L); // open lue file-system library
  lua_gc(L, LUA_GCRESTART, -1); // restart GC
  
  std::string topDefs = createTopLevelDefs(argc, argv);

  // load variable definitions etc
  if (luaL_loadstring(L, topDefs.c_str()) || lua_pcall(L, 0, LUA_MULTRET, 0)) {
    // some error occured
    const char* ret = lua_tostring(L, -1);
    std::cerr << "*** LOAD ERROR *** " << ret << std::endl;
    lua_close(L);
    return finish(1);    
  }

  // load and run input file
  if (luaL_loadfile(L, argv[1]) || lua_pcall(L, 0, LUA_MULTRET, 0)) {
    // some error occured
    const char* ret = lua_tostring(L, -1);
    std::cerr << "*** ERROR *** " << ret << std::endl;
    lua_close(L);
    return finish(1);
  }
  lua_close(L);  

  return finish(0);
}
