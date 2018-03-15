#include <MaxwellModDecl.h> 
double MaxwellVol1xSerP1(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[2]; 
  const double *ez = &q[4]; 
  const double *bx = &q[6]; 
  const double *by = &q[8]; 
  const double *bz = &q[10]; 
  const double *ph = &q[12]; 
  const double *ps = &q[14]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[2]; 
  double *outEz = &out[4]; 
  double *outBx = &out[6]; 
  double *outBy = &out[8]; 
  double *outBz = &out[10]; 
  double *outPh = &out[12]; 
  double *outPs = &out[14]; 
 
  double dx0 = 2.0/dx[0]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  return cflFreq; 
} 
double MaxwellVol1xSerP2(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[3]; 
  const double *ez = &q[6]; 
  const double *bx = &q[9]; 
  const double *by = &q[12]; 
  const double *bz = &q[15]; 
  const double *ph = &q[18]; 
  const double *ps = &q[21]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[3]; 
  double *outEz = &out[6]; 
  double *outBx = &out[9]; 
  double *outBy = &out[12]; 
  double *outBz = &out[15]; 
  double *outPh = &out[18]; 
  double *outPs = &out[21]; 
 
  double dx0 = 2.0/dx[0]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += 3.872983346207417*ph[1]*c2chi*dx0; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 3.872983346207417*bz[1]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += -3.872983346207417*by[1]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 3.872983346207417*ps[1]*dx0*gamma; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += -3.872983346207417*ez[1]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += 3.872983346207417*ey[1]*dx0; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 3.872983346207417*ex[1]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 3.872983346207417*bx[1]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  return cflFreq; 
} 
double MaxwellVol1xSerP3(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[4]; 
  const double *ez = &q[8]; 
  const double *bx = &q[12]; 
  const double *by = &q[16]; 
  const double *bz = &q[20]; 
  const double *ph = &q[24]; 
  const double *ps = &q[28]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[4]; 
  double *outEz = &out[8]; 
  double *outBx = &out[12]; 
  double *outBy = &out[16]; 
  double *outBz = &out[20]; 
  double *outPh = &out[24]; 
  double *outPs = &out[28]; 
 
  double dx0 = 2.0/dx[0]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += 3.872983346207417*ph[1]*c2chi*dx0; 
  outEx[3] += 5.916079783099617*ph[2]*c2chi*dx0+2.645751311064591*ph[0]*c2chi*dx0; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 3.872983346207417*bz[1]*c2*dx0; 
  outEy[3] += 5.916079783099617*bz[2]*c2*dx0+2.645751311064591*bz[0]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += -3.872983346207417*by[1]*c2*dx0; 
  outEz[3] += (-5.916079783099617*by[2]*c2*dx0)-2.645751311064591*by[0]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 3.872983346207417*ps[1]*dx0*gamma; 
  outBx[3] += 5.916079783099617*ps[2]*dx0*gamma+2.645751311064591*ps[0]*dx0*gamma; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += -3.872983346207417*ez[1]*dx0; 
  outBy[3] += (-5.916079783099617*ez[2]*dx0)-2.645751311064591*ez[0]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += 3.872983346207417*ey[1]*dx0; 
  outBz[3] += 5.916079783099617*ey[2]*dx0+2.645751311064591*ey[0]*dx0; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 3.872983346207417*ex[1]*chi*dx0; 
  outPh[3] += 5.916079783099617*ex[2]*chi*dx0+2.645751311064591*ex[0]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 3.872983346207417*bx[1]*c2gamma*dx0; 
  outPs[3] += 5.916079783099617*bx[2]*c2gamma*dx0+2.645751311064591*bx[0]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  return cflFreq; 
} 
double MaxwellVol1xSerP4(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[5]; 
  const double *ez = &q[10]; 
  const double *bx = &q[15]; 
  const double *by = &q[20]; 
  const double *bz = &q[25]; 
  const double *ph = &q[30]; 
  const double *ps = &q[35]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[5]; 
  double *outEz = &out[10]; 
  double *outBx = &out[15]; 
  double *outBy = &out[20]; 
  double *outBz = &out[25]; 
  double *outPh = &out[30]; 
  double *outPs = &out[35]; 
 
  double dx0 = 2.0/dx[0]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += 3.872983346207417*ph[1]*c2chi*dx0; 
  outEx[3] += 5.916079783099617*ph[2]*c2chi*dx0+2.645751311064591*ph[0]*c2chi*dx0; 
  outEx[4] += 7.937253933193772*ph[3]*c2chi*dx0+5.196152422706631*ph[1]*c2chi*dx0; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 3.872983346207417*bz[1]*c2*dx0; 
  outEy[3] += 5.916079783099617*bz[2]*c2*dx0+2.645751311064591*bz[0]*c2*dx0; 
  outEy[4] += 7.937253933193772*bz[3]*c2*dx0+5.196152422706631*bz[1]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += -3.872983346207417*by[1]*c2*dx0; 
  outEz[3] += (-5.916079783099617*by[2]*c2*dx0)-2.645751311064591*by[0]*c2*dx0; 
  outEz[4] += (-7.937253933193772*by[3]*c2*dx0)-5.196152422706631*by[1]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 3.872983346207417*ps[1]*dx0*gamma; 
  outBx[3] += 5.916079783099617*ps[2]*dx0*gamma+2.645751311064591*ps[0]*dx0*gamma; 
  outBx[4] += 7.937253933193772*ps[3]*dx0*gamma+5.196152422706631*ps[1]*dx0*gamma; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += -3.872983346207417*ez[1]*dx0; 
  outBy[3] += (-5.916079783099617*ez[2]*dx0)-2.645751311064591*ez[0]*dx0; 
  outBy[4] += (-7.937253933193772*ez[3]*dx0)-5.196152422706631*ez[1]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += 3.872983346207417*ey[1]*dx0; 
  outBz[3] += 5.916079783099617*ey[2]*dx0+2.645751311064591*ey[0]*dx0; 
  outBz[4] += 7.937253933193772*ey[3]*dx0+5.196152422706631*ey[1]*dx0; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 3.872983346207417*ex[1]*chi*dx0; 
  outPh[3] += 5.916079783099617*ex[2]*chi*dx0+2.645751311064591*ex[0]*chi*dx0; 
  outPh[4] += 7.937253933193772*ex[3]*chi*dx0+5.196152422706631*ex[1]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 3.872983346207417*bx[1]*c2gamma*dx0; 
  outPs[3] += 5.916079783099617*bx[2]*c2gamma*dx0+2.645751311064591*bx[0]*c2gamma*dx0; 
  outPs[4] += 7.937253933193772*bx[3]*c2gamma*dx0+5.196152422706631*bx[1]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  return cflFreq; 
} 
