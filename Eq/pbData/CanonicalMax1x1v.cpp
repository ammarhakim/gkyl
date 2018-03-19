#include <cmath> 
#include <CanonicalModDecl.h> 
double CanonicalVol1x1vMaxP1(const double *w, const double *dxv, const double *H, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double dxInv0 = 1.0/dxv[0]; 
  double dvInv0 = 1.0/dxv[1]; 
  out[1] += 1.5*f[0]*H[2]; 
  out[2] += -1.5*f[0]*H[1]; 
  double cflFreq = 0.0; 
  // calculate phase space speed alpha in each direction, and add to cfl 
  // here alpha is cell-averaged 
  double alpha = 0.0; 
  alpha = fabs(3.464101615137754*H[2]);  // alpha_x0 
  cflFreq += alpha*dxInv0; 
  alpha = fabs(-3.464101615137754*H[1]);  // alpha_v0 
  cflFreq += alpha*dvInv0; 
  return cflFreq; 
} 
double CanonicalVol1x1vMaxP2(const double *w, const double *dxv, const double *H, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double dxInv0 = 1.0/dxv[0]; 
  double dvInv0 = 1.0/dxv[1]; 
  out[1] += 1.5*(2.23606797749979*f[2]*H[5]+f[1]*H[3]+f[0]*H[2]); 
  out[2] += -1.5*(2.23606797749979*f[1]*H[4]+f[2]*H[3]+f[0]*H[1]); 
  out[3] += 1.5*(2.0*H[5]*f[5]+2.23606797749979*f[0]*H[5]-2.0*H[4]*f[4]-2.23606797749979*f[0]*H[4]+H[2]*f[2]-1.0*H[1]*f[1]); 
  out[4] += 1.5*(5.0*f[3]*H[5]+2.0*H[3]*f[4]+2.23606797749979*f[0]*H[3]+2.23606797749979*f[1]*H[2]); 
  out[5] += -1.5*(2.0*H[3]*f[5]+5.0*f[3]*H[4]+2.23606797749979*f[0]*H[3]+2.23606797749979*H[1]*f[2]); 
  double cflFreq = 0.0; 
  // calculate phase space speed alpha in each direction, and add to cfl 
  // here alpha is cell-averaged 
  double alpha = 0.0; 
  alpha = fabs(3.464101615137754*H[2]);  // alpha_x0 
  cflFreq += alpha*dxInv0; 
  alpha = fabs(-3.464101615137754*H[1]);  // alpha_v0 
  cflFreq += alpha*dvInv0; 
  return cflFreq; 
} 
