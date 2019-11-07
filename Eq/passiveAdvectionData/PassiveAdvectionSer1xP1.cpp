#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionVol1xSerP1(const double *w, const double *dxv, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing.
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  const double *v1 = &f[2]; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.25*(2.449489742783178*v1[1]-1.414213562373095*v1[0])*dfac1; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(2.449489742783178*v1[1]+1.414213562373095*v1[0])*dfac1; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.5*(0.7071067811865475*v1[0]-1.224744871391589*v1[1])*dfac1; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.5*(1.224744871391589*v1[1]+0.7071067811865475*v1[0])*dfac1; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 1.224744871391589*(f[1]*v1[1]+f[0]*v1[0])*dfac1; 
  return cflFreq; 
} 
