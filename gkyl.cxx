// Gkyl ------------------------------------------------------------------------
//
// Top-level entry point into Gkyl
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

// gkyl includes
#include <gkylconfig.h>
#include <gkylhgtip.h>
#include <Gkyl.h>

#ifdef HAVE_MPI_H
# include <mpi.h>
#endif

// std includes
#include <fenv.h> 
#include <iostream>
#include <list>
#include <string>
#include <unistd.h>
#include <vector>

// includes for denormalized floats
#if defined(__clang__)
// nothing to include
#elif defined(__powerpc__)
// nothing to include
#elif defined(__GNUC__) || defined(__GNUG__)
# include <xmmintrin.h>
#endif

// Finish simulation
int finish(int err) {
  int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (0 != err)
    MPI_Abort(MPI_COMM_WORLD, err);
  else 
    MPI_Finalize();
  return err;
}

// show usage
void showUsage() {
  std::cout << "This is the Gkeyll code. See gkeyll.rtfd.io for details." << std::endl;
  std::cout << std::endl;

  std::cout << "gkyl [OPTIONS] app/tool [OPTIONS]" << std::endl;
  std::cout << "Available options are" << std::endl;
  std::cout << "  -e chunk   Execute string 'chunk'" << std::endl;
  std::cout << "  -v         Show version information" << std::endl;
  std::cout << std::endl;

  std::cout << "Most app input files take the following commands:" << std::endl;
  std::cout << "  run        Run simulation. Default if nothing specified" << std::endl;
  std::cout << "  init       Only initialize simulation but do not run it" << std::endl;
  std::cout << "  restart    Restart simulation from last saved restart frame" << std::endl;
  std::cout << std::endl;

  std::cout << "Individual tools may take other options and commands. See their specific help." << std::endl;
}

// show version information
void showVersion() {
  std::cout << "Changeset: " << GKYL_HG_CHANGESET << std::endl;
  std::cout << "Built on: " << __DATE__ << " " << __TIME__ << std::endl;
  std::cout << "Built with:";
#ifdef HAVE_MPI_H
  std::cout << " MPI";
#endif      
#ifdef HAVE_ADIOS_H
  std::cout << " Adios";
#endif    
#ifdef HAVE_EIGEN_CORE
  std::cout << " Eigen";
#endif
#ifdef USING_SQLITE3
  std::cout << " Sqlite3";
#endif    
#ifdef HAVE_ZMQ_H
  std::cout << " ZeroMQ";
#endif
  std::cout << std::endl;
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

  std::string luaExpr;
  // parse options: I am not using any heavy-duty tools as basic
  // requirements are very simple. Sadly, even "sophisticated" tools
  // do not support very special need of gkeyll in which app or tool
  // can have complex command parsers of their own
  char c;
  while ((c = getopt(argc, argv, "hve:")) != -1)
    switch (c)
    {
      case 'h':
          showUsage();
          return finish(0);
          break;

      case 'v':
          showVersion();
          return finish(0);

      case 'e':
          luaExpr.append(optarg);
          break;

      case '?':
          return finish(0);
    }

  std::list<std::string> argRest;
  // collect remaining options into a list
  for (; optind < argc; ++optind)
    argRest.push_back(argv[optind]);

  std::string inpFile;
  // get input file name or tool (if any)
  if (argRest.size() > 0) {
    inpFile = argRest.front();
    argRest.pop_front();
  }

  Gkyl app(luaExpr, inpFile, argRest);
  return app.run();
}
