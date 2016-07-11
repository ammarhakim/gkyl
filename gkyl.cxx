/** Top-level entry point into Gkyl */
//    _______     ___
// + 6 @ |||| # P ||| +

// luajit includes
#include <lua.hpp>

// std include
#include <iostream>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

int
main(int argc, char **argv)
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif
  
  if (argc != 2)
  {
    std::cout << "Usage: gkyl LUA-SCRIPT" << std::endl;
    return 1;
  }
  
  lua_State *L = luaL_newstate();
  if (L==NULL)
  {
    std::cerr << "Unable to create a new LuaJIT interpreter state. Quitting" << std::endl;
    return 1;
  }
  lua_gc(L, LUA_GCSTOP, 0);  // stop collector during initialization
  luaL_openlibs(L);  // open libraries
  lua_gc(L, LUA_GCRESTART, -1);

  if (luaL_loadfile(L, argv[1]) || lua_pcall(L, 0, LUA_MULTRET, 0))
  { // some error occured
    const char* ret = lua_tostring(L, -1);
    std::cerr << "ERROR: " << ret << std::endl;
    lua_close(L);
    return 1;
  }

  const char* ret = lua_tostring(L, -1); // value returned by script
  std::cout << ret << std::endl;
  lua_close(L);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif  
  return 0;
}
