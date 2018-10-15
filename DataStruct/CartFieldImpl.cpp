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

void gkylCopyFromField(double *data, double *f, unsigned numComponents, unsigned offset) {
  for (unsigned k=0; k<numComponents; k++) {
    data[k+offset] = f[k];
  }
}

void gkylCopyToField(double *f, double *data, unsigned numComponents, unsigned offset) {
  for (unsigned k=0; k<numComponents; k++) {
    f[k] = data[k+offset];
  }
}

void gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out) {
  for (unsigned n=s; n<(s+nv); ++n) out[n] = val;
}
