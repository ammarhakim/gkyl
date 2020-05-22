// Gkyl ------------------------------------------------------------------------
//
// Top-level entry point into Gkyl
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

// gkyl includes
#include <gkylconfig.h>

// std includes
#include <iostream>
#include <string>

// include kernels declarations
#include <VlasovTmplModDecl.h>

int
main(int argc, char **argv) {
  Gkyl::VlasovModDecl<1,2,1,Gkyl::G_SERENDIPITY_C> vKernels1x1vp1Ser;
  
  return 0;
}
