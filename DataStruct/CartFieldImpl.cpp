// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for CartField
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <cmath>
#include <CartFieldImpl.h>

void gkylCartFieldAccumulate(unsigned nv, double fact, const double *inp, double *out) {
  for (unsigned n=0; n<nv; ++n)
    out[n] += fact*inp[n];
}

void gkylCartFieldAssign(unsigned nv, double fact, const double *inp, double *out) {
  for (unsigned n=0; n<nv; ++n)
    out[n] = fact*inp[n];
}


