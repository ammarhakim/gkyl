#include <cmath> 
#include <CanonicalModDecl.h> 
double CanonicalVol1x1vSerP1(const double *w, const double *dxv, const double *H, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxdvInv0 = 1.0/(dxv[0]*dxv[1]); 
  double dxInv0 = 1.0/dxv[0]; 
  double dvInv0 = 1.0/dxv[1]; 
  out[1] += 6.0*(f[1]*H[3]+f[0]*H[2])*dxdvInv0; 
  out[2] += -6.0*(f[2]*H[3]+f[0]*H[1])*dxdvInv0; 
  out[3] += 6.0*(H[2]*f[2]-1.0*H[1]*f[1])*dxdvInv0; 
  double cflFreq = 0.0; 
  // calculate phase space speed alpha in each direction, and add to cfl 
  // here alpha is cell-averaged 
  cflFreq += fabs(1.732050807568877*H[2]*dvInv0)*dxInv0; 
  cflFreq += fabs(-1.732050807568877*H[1]*dxInv0)*dvInv0; 
  return cflFreq; 
} 
