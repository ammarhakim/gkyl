// Gkyl ----------------------------------------------------------------------------
//
// C++ back-end for ten-moment conservative to primitive (and vice versa) conversion
//    _______     ___
// + 6 @ |||| # P ||| +
//----------------------------------------------------------------------------------

#include <TenMomentConservToPrimImpl.h>
#include <cmath>
#include <stdio.h>

void 
gkylTenMomentConservativeToPrimitive(double *f)
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

  f[1] = u;
  f[2] = v;
  f[3] = w;
  f[4] = pxx;
  f[5] = pxy;
  f[6] = pxz;
  f[7] = pyy; 
  f[8] = pyz;
  f[9] = pzz; 
}
void 
gkylTenMomentPrimitiveToConservative(double *f)
{
  double r = f[0];
  double u = f[1];
  double v = f[2];
  double w = f[3];
  double pxx = f[4];
  double pxy = f[5];
  double pxz = f[6];
  double pyy = f[7];
  double pyz = f[8];
  double pzz = f[9];

  f[1] = r*u;
  f[2] = r*v;
  f[3] = r*w;
  f[4] = pxx+r*u*u;
  f[5] = pxy+r*u*v;
  f[6] = pxz+r*u*w;
  f[7] = pyy+r*v*v; 
  f[8] = pyz+r*v*w;
  f[9] = pzz+r*w*w;   
}