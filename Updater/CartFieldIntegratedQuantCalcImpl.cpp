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
gkylCartFieldIntQuantAbsV(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out)
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
    out[c] += std::abs(fIn[c*nb])*vol;
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

void
gkylCartFieldIntQuantGradPerpV2(int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out)
{
  double vol = 1.0;
  // assume polyOrder >= 1
  for (unsigned d=0; d<ndim; ++d)
    vol *= dxv[d]/2.0;

  double dfac2_x = 4.0/dxv[0]/dxv[0];
  double dfac2_y = 4.0/dxv[1]/dxv[1];

  // assume 1 component
  if(ndim==2) {
    out[0] = vol*((3*fIn[3]*fIn[3]+3*fIn[2]*fIn[2])*dfac2_y+(3*fIn[3]*fIn[3]+3*fIn[1]*fIn[1])*dfac2_x);
  } else if(ndim==3) {
    out[0] = vol*((3*fIn[7]*fIn[7]+3*fIn[6]*fIn[6]+3*fIn[4]*fIn[4]+3*fIn[2]*fIn[2])*dfac2_y+(3*fIn[7]*fIn[7]+3*fIn[5]*fIn[5]+3*fIn[4]*fIn[4]+3*fIn[1]*fIn[1])*dfac2_x);
  }
}
