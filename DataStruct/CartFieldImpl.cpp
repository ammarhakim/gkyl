// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for CartField
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <cmath>
#include <stdio.h>
#include <CartFieldImpl.h>

void gkylCartFieldAccumulate(unsigned s, unsigned nv, double fact, const double *inp, double *out) {
  for (unsigned n=s; n<(s+nv); ++n)
    out[n] += fact*inp[n];
}

void gkylCartFieldAccumulateOffset(unsigned sInp, unsigned sOut, unsigned nCells, unsigned compStart, unsigned nCompInp, unsigned nCompOut, double fact, const double *inp, double *out) {
  if (nCompInp < nCompOut) {
    for (unsigned i=0; i<nCells; ++i) {
      for (unsigned c=0; c<nCompInp; ++c) {
        out[sOut + i*nCompOut + compStart + c] += fact*inp[sInp + i*nCompInp + c];
      }
    }
  }
  else {
    for (unsigned i=0; i<nCells; ++i) {
      for (unsigned c=0; c<nCompOut; ++c) {
        out[sOut + i*nCompOut + c] += fact*inp[sInp + i*nCompInp + compStart + c];
      }
    }
  }
}

void gkylCartFieldAssign(unsigned s, unsigned nv, double fact, const double *inp, double *out) {
  for (unsigned n=s; n<(s+nv); ++n)
    out[n] = fact*inp[n];
}

void gkylCartFieldAssignOffset(unsigned sInp, unsigned sOut, unsigned nCells, unsigned compStart, unsigned nCompInp, unsigned nCompOut, double fact, const double *inp, double *out) {
  if (nCompInp < nCompOut) {
    for (unsigned i=0; i<nCells; ++i) {
      for (unsigned c=0; c<nCompInp; ++c) {
        out[sOut + i*nCompOut + compStart + c] = fact*inp[sInp + i*nCompInp + c];
      }
    }
  }
  else {
    for (unsigned i=0; i<nCells; ++i) {
      for (unsigned c=0; c<nCompOut; ++c) {
        out[sOut + i*nCompOut + c] = fact*inp[sInp + i*nCompInp + compStart + c];
      }
    }
  }
}

void gkylCartFieldScale(unsigned s, unsigned nv, double fact, double *out) {
  for (unsigned n=s; n<(s+nv); ++n)
    out[n] *= fact;
}

void gkylCartFieldScaleByCell(unsigned s, unsigned nv, unsigned ncomp, double *fact, double *out) {
  unsigned icell = 0;
  for (unsigned n=s; n<(s+nv); n+=ncomp) {
    for (unsigned c=0; c<ncomp; c++) {
      out[n+c] *= fact[icell];
    }
    icell++;
  }
}

void gkylCartFieldAbs(unsigned s, unsigned nv, double *out) {
  for (unsigned n=s; n<(s+nv); ++n)
    out[n] = fabs(out[n]);
}

void gkylCopyFromField(double *data, double *f, unsigned numComponents, unsigned offset) {
  for (unsigned k=0; k<numComponents; ++k) {
    data[k+offset] = f[k];
  }
}

void gkylCopyToField(double *f, double *data, unsigned numComponents, unsigned offset) {
  for (unsigned k=0; k<numComponents; ++k) {
    f[k] = data[k+offset];
  }
}

void gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out) {
  for (unsigned n=s; n<(s+nv); ++n) out[n] = val;
}
