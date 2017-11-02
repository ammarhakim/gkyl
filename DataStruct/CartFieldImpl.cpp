// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for CartField
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <cmath>
#include <CartFieldImpl.h>

void gkylCartFieldAccumulate(unsigned s, unsigned nv, double fact, const double *inp, double *out) {
  for (unsigned n=s; n<(s+nv); ++n)
    out[n] += fact*inp[n];
}

void gkylCartFieldAssign(unsigned s, unsigned nv, double fact, const double *inp, double *out) {
  for (unsigned n=s; n<(s+nv); ++n)
    out[n] = fact*inp[n];
}


