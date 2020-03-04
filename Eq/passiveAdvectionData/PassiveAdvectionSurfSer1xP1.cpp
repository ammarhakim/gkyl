#include <PassiveAdvectionModDecl.h> 
double PassiveAdvectionSurf1xSer_X1_P1(const double *cflFreqCtrlL, const double *cflFreqCtrlR, const double *w, const double *dxv, const double dt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. 
  double dfac1 = 2.0/dxv[0]; 
  double w1 = w[0]; 
  const double *v1 = &fr[2]; 
  double incr[2]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.3535533905932737*(1.732050807568877*v1[1]-1.0*v1[0]); 

  double alpha[1]; 
  alpha[0] = -0.5*(2.449489742783178*v1[1]-1.414213562373095*v1[0]); 
  if (alpha0>0) { 
  incr[0] = 0.5*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac1; 
  incr[1] = -0.5*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac1; 
  } else { 
  incr[0] = -0.5*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac1; 
  incr[1] = 0.5*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac1; 
  }
  outr[0] += incr[0]; 
  outr[1] += incr[1]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  return std::abs(alpha0); 
} 
