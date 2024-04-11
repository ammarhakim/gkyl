// Gkyl ------------------------------------------------------------------------
//
// Top-level entry point into Gkyl
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <gkylzero.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#endif

// Input arguments to gkyl executable. You must call release_app_args
// when done with this struct
struct app_args {
  bool use_gpu; // should this be run on GPU?
  bool step_mode; // run for fixed number of steps? (for valgrind/cuda-memcheck)
  bool trace_mem; // should we trace memory allocation/deallocations?
  int num_steps; // number of steps
  bool is_restart; // is this a restarted sim?

  char *echunk; // chunk of lua code to execute
  
  int num_opt_args; // number of optional arguments
  char **opt_args; // optional arguments
};

static int
calc_output_prefix_len(const char *fn)
{
  const char *suff = strrchr(fn, '.');
  return strlen(fn) - (suff ? strlen(suff) : 0);
}

static const char *get_fname(const char *fn) { return strrchr(fn, '/'); }

// show usage
static void
show_usage()
{
  printf("This is the Gkeyll code. See gkeyll.rtfd.io for details.\n");
  printf("Type 'gkyl man' for help.\n\n");

  printf("gkyl [OPTIONS] app/tool [APP-OPTIONS]\n");
  printf("Available options are\n");
  printf("  -e chunk   Execute string 'chunk'\n");
  printf("  -t         Show list of registered tools\n");
  printf("  -v         Show version information\n");
  printf("  -g         Run on NVIDIA GPU (if available and built with CUDA)\n\n");
  printf("  -m         Run memory tracer\n");    

  printf("Most app input files take the following commands:\n");
  printf("  run        Run simulation. Default if nothing specified\n");
  printf("  init       Only initialize simulation but do not run it\n");
  printf("  restart    Restart simulation from last saved restart frame\n\n");

  printf("Individual tools may take other options and commands. See their specific help.\n");
}

static struct app_args*
parse_app_args(int argc, char **argv)
{
  bool use_gpu = false;
  bool trace_mem = false;

  struct app_args *args = gkyl_malloc(sizeof(*args));

  args->use_gpu = false;
  args->step_mode = false;
  args->num_steps = INT_MAX;
  args->is_restart = false;
  args->num_opt_args = 0;
  args->echunk = 0;

  int c;
  while ((c = getopt(argc, argv, "+hvtme:g")) != -1) {
    switch (c)
    {
      case 'h':
        show_usage();
        exit(-1);
        break;

      case 't':
        // TODO
        break;

      case 'e':
        args->echunk = gkyl_malloc(strlen(optarg)+1);
        strcpy(args->echunk, optarg);
        break;
        
      case 'g':
        use_gpu = true;
        break;

      case 'm':
        trace_mem = true;
        break;

      case '?':
        break;
    }
  }

  args->use_gpu = use_gpu;
  args->trace_mem = trace_mem;

  args->num_opt_args = 0;
  // collect remaining options into a list
  for (int oind=optind; oind < argc; ++oind) args->num_opt_args += 1;
  args->opt_args = gkyl_malloc(sizeof(char*)*args->num_opt_args);

  for (int i=0, oind=optind; oind < argc; ++oind, ++i) {
    args->opt_args[i] = gkyl_malloc(strlen(argv[oind]+1));
    strcpy(args->opt_args[i], argv[oind]);
  }
  return args;
}

static void
release_opt_args(struct app_args *args)
{
  for (int i=0; i<args->num_opt_args; ++i)
    gkyl_free(args->opt_args[i]);
  gkyl_free(args->opt_args);
  if (args->echunk)
    gkyl_free(args->echunk);
  gkyl_free(args);
}

int
main(int argc, char **argv)
{
  struct app_args *app_args = parse_app_args(argc, argv);
  
#ifdef GKYL_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  if (app_args->trace_mem) {
    gkyl_cu_dev_mem_debug_set(true);
    gkyl_mem_debug_set(true);
  }
  
  lua_State *L = luaL_newstate();
  lua_gc(L, LUA_GCSTOP, 0);
  luaL_openlibs(L);
  gkyl_zero_lw_openlibs(L);
  gkyl_vlasov_lw_openlibs(L);
  gkyl_moment_lw_openlibs(L);
  lua_gc(L, LUA_GCRESTART, -1);

#ifdef GKYL_HAVE_MPI
  lua_pushboolean(L, 1);
  lua_setglobal(L, "GKYL_HAVE_MPI");
  
  lua_pushlightuserdata(L, MPI_COMM_WORLD);
  lua_setglobal(L, "GKYL_MPI_COMM");
#else
  lua_pushboolean(L, 0);
  lua_setglobal(L, "GKYL_HAVE_MPI");

  lua_pushboolean(L, 0);
  lua_setglobal(L, "GKYL_MPI_COMM");
#endif

  // run Lua code (if it exists) before running input file
  if (app_args->echunk)
    glua_run_lua(L, app_args->echunk, strlen(app_args->echunk), 0);

  if (app_args->num_opt_args > 0) {
    // read input file and run it
    const char *inp_name = app_args->opt_args[0];
     if (gkyl_check_file_exists(inp_name)) {
      // set file prefix as a global (THIS WILL BE REPLACED!!)
      const char *suff = get_fname(inp_name);
      const char *suff1 = suff ? suff+1 : inp_name;
      lua_pushlstring(L, suff1, calc_output_prefix_len(suff1));
      lua_setglobal(L, "GKYL_OUT_PREFIX");
      
      int64_t sz;
      char *buff = gkyl_load_file(inp_name, &sz);
      glua_run_lua(L, buff, sz, stderr);
      gkyl_free(buff);
    }
  }
  
  lua_close(L);  
  
  release_opt_args(app_args);

#ifdef GKYL_HAVE_MPI
  MPI_Finalize();
#endif
  
  return 0;
}
