#include <CartFieldIntegratedQuantCalcImpl.h>
#include <cmath>

void
gkylCartFieldIntQuantV(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out)
{
  double vol = 1.0;
  if (nb>1)
  {
    for (unsigned d=0; d<ndim; ++d)
      vol *= dxv[d]/(std::sqrt(2.0));
  }
  else
  { // for polyOrder = 0 no normalization is applied
    for (unsigned d=0; d<ndim; ++d)
      vol *= dxv[d];
  }

  for (unsigned c=0; c<nc; ++c)
    out[c] += fIn[c*nb]*vol;
}

void
gkylCartFieldIntQuantV2(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out)
{
  double vol = 1.0;
  if (nb>1)
  {
    for (unsigned d=0; d<ndim; ++d)
      vol *= dxv[d]/2.0;
  }
  else
  { // for polyOrder = 0 no normalization is applied
    for (unsigned d=0; d<ndim; ++d)
      vol *= dxv[d];
  }

  for (unsigned c=0; c<nc; ++c)
  {
    double v2 = 0.0;
    for (unsigned b=0; b<nb; ++b)
      v2 += fIn[c*nb+b]*fIn[c*nb+b];
    out[c] += v2*vol;
  }
}
