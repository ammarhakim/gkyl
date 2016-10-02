// Gkyl ------------------------------------------------------------------------
//
// Top-level entry point into Gkyl
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <gkylconfig.h>

// luajit includes
#include <lua.hpp>

// std include
#include <iostream>
#include <vector>

#ifdef HAVE_MPI_H
# include <mpi.h>
# include <GkMpiFuncs.h>
#endif

#ifdef HAVE_ADIOS_H
# include <adios.h>
#endif

/* A simple logger for parallel simulations */
class DumDumLogger {
  public:
    DumDumLogger() : rank(0) {
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

int finish(int err) {
#ifdef HAVE_MPI_H
  MPI_Finalize();
#endif  
  return err;  
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
  DumDumLogger logger;

  if (argc != 2) {
    logger.log("Usage: gkyl LUA-SCRIPT");
    return finish(0);
  }

  // initialize LuaJIT and load libraries
  lua_State *L = luaL_newstate();
  if (L==NULL) {
    logger.log("Unable to create a new LuaJIT interpreter state. Quitting");
    return finish(0);
  }
  lua_gc(L, LUA_GCSTOP, 0);  // stop collector during initialization
  luaL_openlibs(L);  // open standard libraries
  lua_gc(L, LUA_GCRESTART, -1);

  if (luaL_loadfile(L, argv[1]) || lua_pcall(L, 0, LUA_MULTRET, 0)) {
    // some error occured
    const char* ret = lua_tostring(L, -1);
    logger.log("*** ERROR ***"); logger.log(ret);
    lua_close(L);
    return finish(0);
  }
  lua_close(L);
  return finish(0);
}
