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

void gkylCartFieldScale(unsigned s, unsigned nv, double fact, double *out) {
  for (unsigned n=s; n<(s+nv); ++n)
    out[n] *= fact;
}

void gkylCartFieldAbs(unsigned s, unsigned nv, double *out) {
  for (unsigned n=s; n<(s+nv); ++n)
    out[n] = fabs(out[n]);
}

void copyFromField(double *data, double *f, unsigned numComponents, unsigned c) {
  for (unsigned k=0; k<numComponents; k++) {
    data[c] = f[k];
    c+=1;
  }
}

void copyToField(double *f, double *data, unsigned numComponents, unsigned c) {
  for (unsigned k=0; k<numComponents; k++) {
    f[k] = data[c];
    c+=1;
  }
}
