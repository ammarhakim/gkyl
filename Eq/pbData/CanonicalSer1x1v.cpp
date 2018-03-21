#include <cmath> 
#include <CanonicalModDecl.h> 
double CanonicalVol1x1vSerP1(const double *w, const double *dxv, const double *H, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double dxInv0 = 1.0/dxv[0]; 
  double dvInv0 = 1.0/dxv[1]; 
  out[1] += 1.5*(f[1]*H[3]+f[0]*H[2]); 
  out[2] += -1.5*(f[2]*H[3]+f[0]*H[1]); 
  out[3] += 1.5*(H[2]*f[2]-1.0*H[1]*f[1]); 
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
double CanonicalVol1x1vSerP2(const double *w, const double *dxv, const double *H, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxv0 = dxv[0]; 
  double dxv1 = dxv[1]; 
  double dxInv0 = 1.0/dxv[0]; 
  double dvInv0 = 1.0/dxv[1]; 
  out[1] += 0.1*(33.54101966249684*f[3]*H[7]+15.0*f[4]*H[6]+33.54101966249685*f[2]*H[5]+15.0*f[1]*H[3]+15.0*f[0]*H[2]); 
  out[2] += -0.1*(15.0*f[5]*H[7]+33.54101966249684*f[3]*H[6]+33.54101966249685*f[1]*H[4]+15.0*f[2]*H[3]+15.0*f[0]*H[1]); 
  out[3] += 0.5*(3.0*H[7]*f[7]+6.708203932499369*f[1]*H[7]-3.0*H[6]*f[6]-6.708203932499369*f[2]*H[6]+6.0*H[5]*f[5]+6.708203932499369*f[0]*H[5]-6.0*H[4]*f[4]-6.708203932499369*f[0]*H[4]+3.0*H[2]*f[2]-3.0*H[1]*f[1]); 
  out[4] += 0.1*(67.0820393249937*f[6]*H[7]+75.00000000000001*f[2]*H[7]+30.0*f[1]*H[6]+75.0*f[3]*H[5]+30.0*H[3]*f[4]+33.54101966249685*f[0]*H[3]+33.54101966249685*f[1]*H[2]); 
  out[5] += -0.1*(67.0820393249937*H[6]*f[7]+30.0*f[2]*H[7]+75.00000000000001*f[1]*H[6]+30.0*H[3]*f[5]+75.0*f[3]*H[4]+33.54101966249685*f[0]*H[3]+33.54101966249685*H[1]*f[2]); 
  out[6] += 0.1*(67.0820393249937*H[5]*f[7]+67.0820393249937*f[5]*H[7]+67.0820393249937*f[4]*H[7]+75.0*f[0]*H[7]+15.0*H[3]*f[6]+75.00000000000001*f[1]*H[5]-15.0*H[1]*f[4]-30.0*f[1]*H[4]+33.54101966249684*H[2]*f[3]+33.54101966249684*f[2]*H[3]); 
  out[7] += -0.1*(15.0*H[3]*f[7]+67.0820393249937*H[4]*f[6]+67.0820393249937*f[5]*H[6]+67.0820393249937*f[4]*H[6]+75.0*f[0]*H[6]-15.0*H[2]*f[5]-30.0*f[2]*H[5]+75.00000000000001*f[2]*H[4]+33.54101966249684*H[1]*f[3]+33.54101966249684*f[1]*H[3]); 
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
