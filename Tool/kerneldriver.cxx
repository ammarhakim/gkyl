// Gkyl ------------------------------------------------------------------------
//
// Kernel driver tool: runs various kernels for timing and profiling
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

// gkyl includes
#include <gkylconfig.h>

#include <GkBasisTypes.h>
#include <GkKernelDriver.h>
#include <GkVlasovKernelDriver.h>

// std includes
#include <iostream>
#include <string>

int
main(int argc, char **argv) {
  Gkyl::KernelDriver *driver = new Gkyl::VlasovKernelDriver(3, 3, 1, Gkyl::G_SERENDIPITY_C);
  driver->cellBasedUpdate(1000000);
  
  return 0;
}
