// Gkyl ------------------------------------------------------------------------
//
// Kernel driver tool: runs various kernels for timing and profiling
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

// gkyl includes
#include <gkylconfig.h>
#include <GkKernelDriver.h>
#include <GkBasisTypes.h>

// std includes
#include <iostream>
#include <string>


int
main(int argc, char **argv) {
  Gkyl::KernelDriver *driver = new Gkyl::VlasovKernelDriver(1, 3, 2, Gkyl::G_SERENDIPITY_C);
  driver->cellUpdate(100);
  
  return 0;
}
