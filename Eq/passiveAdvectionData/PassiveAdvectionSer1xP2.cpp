#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionVol1xSerP2(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing.
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  const double *v1 = &f[3]; 
  double cflRate = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.3535533905932737*(2.23606797749979*v1[2]-1.732050807568877*v1[1]+v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.3535533905932737*(2.23606797749979*v1[2]+1.732050807568877*v1[1]+v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.5*(1.58113883008419*v1[2]-1.224744871391589*v1[1]+0.7071067811865475*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.5*(1.58113883008419*v1[2]-1.224744871391589*v1[1]+0.7071067811865475*v1[0])*dfac1; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.5*(1.58113883008419*v1[2]+1.224744871391589*v1[1]+0.7071067811865475*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.5*(1.58113883008419*v1[2]+1.224744871391589*v1[1]+0.7071067811865475*v1[0])*dfac1; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 1.224744871391589*(f[2]*v1[2]+f[1]*v1[1]+f[0]*v1[0])*dfac1; 
  out[2] += 0.7071067811865475*(3.464101615137754*(f[1]*v1[2]+v1[1]*f[2])+3.872983346207417*(f[0]*v1[1]+v1[0]*f[1]))*dfac1; 
  return cflRate; 
} 
