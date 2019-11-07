#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionVol2xSerP1(const double *w, const double *dxv, double *cflFreqCtrl, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing.
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  double dfac2 = 2.0/dxv[1]; 
  double w2 = w[1]; 
  const double *v1 = &f[4]; 
  const double *v2 = &f[8]; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.25*(1.732050807568877*v1[1]-1.0*v1[0])*dfac1; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(1.732050807568877*v1[1]+v1[0])*dfac1; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.8660254037844386*v1[3]-0.5*v1[2]-0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*((-0.8660254037844386*v1[3])+0.5*v1[2]-0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*((-0.8660254037844386*v1[3])-0.5*v1[2]+0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.8660254037844386*v1[3]+0.5*v1[2]+0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.25*(1.732050807568877*v2[2]-1.0*v2[0])*dfac2; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(1.732050807568877*v2[2]+v2[0])*dfac2; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.8660254037844386*v2[3]-0.8660254037844386*v2[2]-0.5*v2[1]+0.5*v2[0])*dfac2; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.5*(v2[1]+v2[0])-0.8660254037844386*(v2[3]+v2[2]))*dfac2; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*((-0.8660254037844386*v2[3])+0.8660254037844386*v2[2]-0.5*v2[1]+0.5*v2[0])*dfac2; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.8660254037844386*(v2[3]+v2[2])+0.5*(v2[1]+v2[0]))*dfac2; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphaCtrl;
  alphaCtrl = (0.5*v1[3]-0.2886751345948129*v1[2]-0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflFreqCtrl[0] = -0.5*(alphaCtrl-std::abs(alphaCtrl)); 
  alphaCtrl = (0.5*v2[3]-0.8660254037844386*v2[2]-0.2886751345948129*v2[1]+0.5*v2[0])*dfac2; 
  cflFreqCtrl[0] += -0.5*(alphaCtrl-std::abs(alphaCtrl)); 
  alphaCtrl = ((-0.5*v1[3])-0.2886751345948129*v1[2]+0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflFreqCtrl[1] = 0.5*(alphaCtrl+std::abs(alphaCtrl)); 
  alphaCtrl = ((-0.5*v2[3])-0.8660254037844386*v2[2]+0.2886751345948129*v2[1]+0.5*v2[0])*dfac2; 
  cflFreqCtrl[1] += -0.5*(alphaCtrl-std::abs(alphaCtrl)); 
  alphaCtrl = ((-0.5*v1[3])+0.2886751345948129*v1[2]-0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflFreqCtrl[2] = -0.5*(alphaCtrl-std::abs(alphaCtrl)); 
  alphaCtrl = ((-0.5*v2[3])+0.8660254037844386*v2[2]-0.2886751345948129*v2[1]+0.5*v2[0])*dfac2; 
  cflFreqCtrl[2] += 0.5*(alphaCtrl+std::abs(alphaCtrl)); 
  alphaCtrl = (0.5*v1[3]+0.2886751345948129*v1[2]+0.8660254037844386*v1[1]+0.5*v1[0])*dfac1; 
  cflFreqCtrl[3] = 0.5*(alphaCtrl+std::abs(alphaCtrl)); 
  alphaCtrl = (0.5*v2[3]+0.8660254037844386*v2[2]+0.2886751345948129*v2[1]+0.5*v2[0])*dfac2; 
  cflFreqCtrl[3] += 0.5*(alphaCtrl+std::abs(alphaCtrl)); 
  out[1] += 0.8660254037844386*(f[3]*v1[3]+f[2]*v1[2]+f[1]*v1[1]+f[0]*v1[0])*dfac1; 
  out[2] += 0.8660254037844386*(f[3]*v2[3]+f[2]*v2[2]+f[1]*v2[1]+f[0]*v2[0])*dfac2; 
  out[3] += 0.8660254037844386*((f[2]*v2[3]+v2[2]*f[3]+f[0]*v2[1]+v2[0]*f[1])*dfac2+(f[1]*v1[3]+v1[1]*f[3]+f[0]*v1[2]+v1[0]*f[2])*dfac1); 
  return cflFreq; 
} 
