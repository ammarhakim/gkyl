// Gkyl ------------------------------------------------------------------------
//
// Top-level entry point into Gkyl
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <gkylconfig.h>

// library includes
#include <lua.hpp>

// std include
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>

#ifdef HAVE_MPI_H
# include <mpi.h>
# include <GkMpiFuncs.h>
#endif

#ifdef HAVE_ADIOS_H
# include <adios.h>
#endif

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
#ifdef HAVE_MPI_H
  if (err != 0)
    MPI_Abort(MPI_COMM_WORLD, err);
  else 
    MPI_Finalize();
#endif
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
  return std::string(path, dirname_len); // path just of directory
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
          << execPath << "/?/init.lua" << "\"";
  
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

  return varDefs.str();
}

// These dummy calls are a hack to force the linker to pull in symbols
// from static libraries. There is a better way to do this, I just do
// not know it yet.
#ifdef HAVE_ADIOS_H
int _adios_init(const char *cf, MPI_Comm comm) { return adios_init(cf, comm); }
#endif
#ifdef HAVE_MPI_H
MPI_Comm _get_MPI_COMM_WORLD() { return get_MPI_COMM_WORLD(); }
#endif

int
main(int argc, char **argv) {
#ifdef HAVE_MPI_H
  MPI_Init(&argc, &argv);
#endif
  Logger logger;

  if (argc != 2) {
    logger.log("Usage: gkyl LUA-SCRIPT");
    return finish(0);
  }

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

  // check if file exists  
  std::string inpFile(argv[1]);
  std::ifstream f(inpFile.c_str());
  if (!f.good()) {
    std::cerr << "Unable to open input file '" << inpFile << "'" << std::endl;
    std::cerr << "Usage: gkyl LUA-SCRIPT" << std::endl;
    return finish(1);
  }
  f.close();

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
