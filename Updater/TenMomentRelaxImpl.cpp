// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for ten-moment relaxation
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <TenMomentRelaxImpl.h>
#include <cmath>
#include <stdio.h>

void
gkylTenMomentRelaxRk3(TenMomentRelaxData_t *sd, FluidData_t *fd, double dt, double *f, double *em, double myK)
{
}

void
gkylTenMomentRelaxExplicit(TenMomentRelaxData_t *sd, FluidData_t *fd, double dt, double *f, double *em, double *staticEm, double myK)
{
  double r = f[0];
  double u = f[1]/r;
  double v = f[2]/r;
  double w = f[3]/r;
  double pxx = f[4]-r*u*u;
  double pxy = f[5]-r*u*v;
  double pxz = f[6]-r*u*w;
  double pyy = f[7]-r*v*v;
  double pyz = f[8]-r*v*w;
  double pzz = f[9]-r*w*w;

  double p = (pxx+pyy+pzz)/3.0;
  double vt = std::sqrt(p/r);
  double edt = std::exp(-vt*myK*dt);

// compute updated pressure tensor component
  f[4] = (pxx-p)*edt+p + r*u*u;
  f[5] = pxy*edt + r*u*v;
  f[6] = pxz*edt + r*u*w;
  f[7] = (pyy-p)*edt+p + r*v*v;
  f[8] = pyz*edt + r*v*w;
  f[9] = (pzz-p)*edt+p + r*w*w;
}
