// Gkyl ------------------------------------------------------------------------
//
// Kernel driver tool: runs various kernels for timing and profiling
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

// gkyl includes
#include <gkylconfig.h>

#include <GkylBasisTypes.h>
#include <GkylKernelDriver.h>
#include <GkylRange.h>
#include <GkylVlasovKernelDriver.h>

#include <cxxopts.hpp>

// std includes
#include <chrono>
#include <iostream>
#include <memory>
#include <string>

int
main(int argc, char **argv) {
  cxxopts::Options options("kerneldriver", "Run kernel timers");
  options.add_options()
    ("h,help", "Show usage")
    ("e,equation", "Equation system, one of 'vlasov', 'gyrokinetic'", cxxopts::value<std::string>()->default_value("vlasov"))
    ("b,basis", "Basis functions, one of 'ms', 'mo'", cxxopts::value<std::string>()->default_value("ms"))
    ("c,cdim", "Config space dimensions", cxxopts::value<unsigned>()->default_value("1"))
    ("v,vdim", "Velocity space dimensions", cxxopts::value<unsigned>()->default_value("1"))
    ("p,polyOrder", "Polynomial order", cxxopts::value<unsigned>()->default_value("1"))
    ("n,nloop", "Number of times to run kernels", cxxopts::value<unsigned>()->default_value("100"))
    ("f,faceupdates", "Face based updates (default: false)", cxxopts::value<bool>()->default_value("false"))
    ("d,cellupdates", "Cell based updates", cxxopts::value<bool>()->default_value("true"))
    ;

  auto result = options.parse(argc, argv);

  if (result.count("help")) {
    std::cout << options.help() << std::endl;
    return 0;
  }

  auto cdim = result["cdim"].as<unsigned>();
  auto vdim = result["vdim"].as<unsigned>();
  auto polyOrder = result["polyOrder"].as<unsigned>();
  auto nloop = result["nloop"].as<unsigned>();
  auto basis = result["basis"].as<std::string>();

  unsigned basisType = Gkyl::G_SERENDIPITY_C;
  if (basis == "mo")
    basisType = Gkyl::G_MAX_ORDER_C;
  else if (basis == "ms")
    basisType = Gkyl::G_SERENDIPITY_C;
  else {
    std::cout << "Basis '" << basis << "' not recognized. Must be 'ms' or 'mo'" << std::endl;
    return 0;
  }

  std::unique_ptr<Gkyl::KernelDriver> driver;  
  auto eqn = result["equation"].as<std::string>();
  // construct driver based on equation system specified
  if (eqn == "vlasov") {
    driver = std::make_unique<Gkyl::VlasovKernelDriver>(cdim, vdim, polyOrder, basisType);
  }
  else if (eqn == "gyrokinetic") {
    std::cout << "Equation system '" << eqn << "' NYI!" << std::endl;
    return 0;    
  }
  else {
    std::cout << "Equation system '" << eqn << "' not recognized" << std::endl;
    return 0;
  }

  if (result["cellupdates"].as<bool>()) {
    std::cout << "Running cell-based update kernels ..." << std::endl;
    
    auto t1 = std::chrono::high_resolution_clock::now();
    driver->cellBasedUpdate(nloop);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto tm = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    std::cout << "... completed in " << tm.count() << " sec" << std::endl;
  }
  if (result["faceupdates"].as<bool>()) {
    std::cout << "Running face-based update kernels ..." << std::endl;

    auto t1 = std::chrono::high_resolution_clock::now();
    driver->faceBasedUpdate(nloop);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto tm = std::chrono::duration_cast<std::chrono::duration<double>>(t2-t1);
    std::cout << "... completed in " << tm.count() << " sec" << std::endl;
  }
  
  return 0;
}
