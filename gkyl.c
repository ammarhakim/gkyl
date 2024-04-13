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
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <whereami.h>

#include <gkylzero.h>

#ifdef GKYL_HAVE_MPI
#include <mpi.h>
#endif

// Tool description struct
struct tool_description {
  const char *tool_name;
  const char *tool_lua;
  const char *tool_help;
};

// List of available Tools
static struct tool_description tool_list[] = {
  {"man", "man.lua", "Gkeyll online manual"},
  {"woman", "man.lua", "Gkeyll online manual (Woe without man)"},
  {"queryrdb", "queryrdb.lua", "Query/modify regression test DB"},
  {"runregression", "runregression.lua", "Run regression/unit tests"},
  {"comparefiles", "comparefiles.lua", "Compare two BP files"},
  {"exacteulerrp", "exacteulerrp.lua", "Exact Euler Riemann problem solver"},
  {"multimomlinear", "multimomlinear.lua",
   "Linear dispersion solver for multi-moment, multifluid equations"},
  {"eqdskreader", "eqdskreader.lua",
   "Read eqdsk file, writing to ADIOS-BP files"},
#ifdef GKYL_HAVE_NCCL
  {"deviceinfo", "deviceinfo.lua", "Information about device"},
#endif
  {0, 0}
};

static int max2(int a, int b) { return a>b ? a : b; }

// Show list of available Tools
static void
show_tool_list(void)
{
  printf("Following tools are available. Query tool help for more information.\n\n");
  for (int i=0; tool_list[i].tool_name != 0; ++i)
    printf("%s %s\n", tool_list[i].tool_name, tool_list[i].tool_help);
}

// Returns tool Lua script name given tool name. Returns 0 if Tool
// does no exist
static const char *
get_tool_from_name(const char *nm)
{
  for (int i=0; tool_list[i].tool_name != 0; ++i)
    if (strcmp(tool_list[i].tool_name, nm) == 0)
      return tool_list[i].tool_lua;
  return 0;
}

static char *
find_exec_path(void)
{
  int len = wai_getExecutablePath(NULL, 0, NULL);
  char *path = gkyl_malloc(len+1);
  int dirname_len; wai_getExecutablePath(path, len, &dirname_len);
  path[dirname_len] = '\0'; // only directory path is returned
  return path;
}

// Input arguments to gkyl executable. You must call release_app_args
// when done with this struct
struct app_args {
  bool use_gpu; // should this be run on GPU?
  bool step_mode; // run for fixed number of steps? (for valgrind/cuda-memcheck)
  bool trace_mem; // should we trace memory allocation/deallocations?
  int num_steps; // number of steps
  bool is_restart; // is this a restarted sim?
  bool use_mpi; // should we use MPI?

  char *echunk; // chunk of lua code to execute
  
  int num_opt_args; // number of optional arguments
  char **opt_args; // optional arguments

  char *exec_path; // location of executable
};

static void
release_opt_args(struct app_args *args)
{
  for (int i=0; i<args->num_opt_args; ++i)
    gkyl_free(args->opt_args[i]);
  gkyl_free(args->opt_args);
  
  if (args->echunk)
    gkyl_free(args->echunk);

  gkyl_free(args->exec_path);
  gkyl_free(args);
}

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
  printf("  -S         Do not initialize MPI\n");  

  printf("Most app input files take the following commands:\n");
  printf("  run        Run simulation. Default if nothing specified\n");
  printf("  init       Only initialize simulation but do not run it\n");
  printf("  restart    Restart simulation from last saved restart frame\n\n");

  printf("Individual tools may take other options and commands. See their specific help.\n");
}

static struct app_args*
parse_app_args(int argc, char **argv)
{
  struct app_args *args = gkyl_malloc(sizeof(*args));

  args->use_gpu = false;
  args->step_mode = false;

#ifdef GKYL_HAVE_MPI  
  args->use_mpi = true;
#else
  args->use_mpi = false;
#endif
  
  args->trace_mem = false;
  args->is_restart = false;
  args->num_steps = INT_MAX;
  args->num_opt_args = 0;
  args->echunk = 0;

  int c;
  while ((c = getopt(argc, argv, "+hvtmSe:g")) != -1) {
    switch (c)
    {
      case 'h':
        show_usage();
        exit(-1);
        break;

      case 't':
        show_tool_list();
        exit(1);
        break;

      case 'e':
        args->echunk = gkyl_malloc(strlen(optarg)+1);
        strcpy(args->echunk, optarg);
        break;
        
      case 'g':
        args->use_gpu = true;
        break;

      case 'S':
        args->use_mpi = false;
        break;
        
      case 'm':
        args->trace_mem = true;
        break;

      case '?':
        break;
    }
  }

  args->num_opt_args = 0;
  // collect remaining options into a list
  for (int oind=optind; oind < argc; ++oind) args->num_opt_args += 1;
  args->opt_args = gkyl_malloc(sizeof(char*)*args->num_opt_args);

  for (int i=0, oind=optind; oind < argc; ++oind, ++i) {
    args->opt_args[i] = gkyl_malloc(strlen(argv[oind])+1);
    strcpy(args->opt_args[i], argv[oind]);
  }

  args->exec_path = find_exec_path();
  
  return args;
}

int
main(int argc, char **argv)
{
  struct app_args *app_args = parse_app_args(argc, argv);
  
  if (app_args->use_mpi)
    MPI_Init(&argc, &argv);

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

  if (app_args->use_mpi) {
    lua_pushboolean(L, true);
    lua_setglobal(L, "GKYL_HAVE_MPI");
    
    lua_pushlightuserdata(L, MPI_COMM_WORLD);
    lua_setglobal(L, "GKYL_MPI_COMM");
  }
  else {
    lua_pushboolean(L, false);
    lua_setglobal(L, "GKYL_HAVE_MPI");
    
    lua_pushboolean(L, false);
    lua_setglobal(L, "GKYL_MPI_COMM");
  }

#ifdef GKYL_HAVE_NCCL
    lua_pushboolean(L, true);
    lua_setglobal(L, "GKYL_HAVE_CUDA");
#else
    lua_pushboolean(L, false);
    lua_setglobal(L, "GKYL_HAVE_CUDA");
#endif  

  lua_pushstring(L, app_args->exec_path);
  lua_setglobal(L, "GKYL_EXEC_PATH");

  // push extra arguments into a Lua table to Tools and App can get
  // them
  lua_newtable(L);
  for (int i=1; i<app_args->num_opt_args; ++i) {
    lua_pushinteger(L, i);
    lua_pushstring(L, app_args->opt_args[i]);
    lua_rawset(L, -3);
  }
  lua_setglobal(L, "GKYL_COMMANDS");

  // set package paths so we find installed libraries
  do {
    const char *fmt = "package.path = package.path .. \";%s/?.lua;%s/Lib/?.lua;%s/?/init.lua\"";
    size_t len = gkyl_calc_strlen(fmt, app_args->exec_path, app_args->exec_path, app_args->exec_path);

    char *str = gkyl_malloc(len+1);
    snprintf(str, len+1, fmt, app_args->exec_path, app_args->exec_path, app_args->exec_path);
    glua_run_lua(L, str, strlen(str), 0);
    
    gkyl_free(str);
  } while (0);

  // run Lua code (if it exists) before running input file
  if (app_args->echunk)
    glua_run_lua(L, app_args->echunk, strlen(app_args->echunk), 0);

  if (app_args->num_opt_args > 0) {

    const char *inp_name = app_args->opt_args[0];
    if (gkyl_check_file_exists(inp_name)) {

      const char *suff = get_fname(inp_name);
      const char *suff1 = suff ? suff+1 : inp_name;
      lua_pushlstring(L, suff1, calc_output_prefix_len(suff1));
      lua_setglobal(L, "GKYL_OUT_PREFIX");
      
      int64_t sz = 0;
      char *buff = gkyl_load_file(inp_name, &sz);
      glua_run_lua(L, buff, sz, stderr);
      gkyl_free(buff);
    }
    else {
      // check if it is a Tool, and run it if so
      const char *tlua = get_tool_from_name(app_args->opt_args[0]);
      if (tlua) {
        const char *fmt = "%s/Tool/%s";
        size_t len = gkyl_calc_strlen(fmt, app_args->exec_path, tlua);
        char *tool_name = gkyl_malloc(len+1);
        snprintf(tool_name, len+1, fmt, app_args->exec_path, tlua);
        
        int64_t sz = 0;
        char *buff = gkyl_load_file(tool_name, &sz);
        glua_run_lua(L, buff, sz, stderr);
        gkyl_free(buff);
        
        gkyl_free(tool_name);
      }
    }
  }
  
  lua_close(L);  
  
  if (app_args->use_mpi)
    MPI_Finalize();

  release_opt_args(app_args);
  
  return 0;
}
