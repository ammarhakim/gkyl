#include <MaxwellModDecl.h> 
double MaxwellVol3xMaxP1(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
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
  double dx1 = 2.0/dx[1]; 
  double dx2 = 2.0/dx[2]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*by[0]*c2*dx2; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += -1.732050807568877*bx[0]*c2*dx2; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*ph[0]*c2chi*dx2; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += -1.732050807568877*ey[0]*dx2; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ex[0]*dx2; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ps[0]*dx2*gamma; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ez[0]*chi*dx2; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*bz[0]*c2gamma*dx2; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  cflFreq += meq->c/dx[1]; 
  cflFreq += meq->c/dx[2]; 
  return cflFreq; 
} 
double MaxwellVol3xMaxP2(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[10]; 
  const double *ez = &q[20]; 
  const double *bx = &q[30]; 
  const double *by = &q[40]; 
  const double *bz = &q[50]; 
  const double *ph = &q[60]; 
  const double *ps = &q[70]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[10]; 
  double *outEz = &out[20]; 
  double *outBx = &out[30]; 
  double *outBy = &out[40]; 
  double *outBz = &out[50]; 
  double *outPh = &out[60]; 
  double *outPs = &out[70]; 
 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 
  double dx2 = 2.0/dx[2]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*by[0]*c2*dx2; 
  outEx[4] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 
  outEx[5] += 1.732050807568877*by[1]*c2*dx2+1.732050807568877*ph[3]*c2chi*dx0; 
  outEx[6] += 1.732050807568877*by[2]*c2*dx2-1.732050807568877*bz[3]*c2*dx1; 
  outEx[7] += 3.872983346207417*ph[1]*c2chi*dx0; 
  outEx[8] += -3.872983346207417*bz[2]*c2*dx1; 
  outEx[9] += 3.872983346207417*by[3]*c2*dx2; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += -1.732050807568877*bx[0]*c2*dx2; 
  outEy[4] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 
  outEy[5] += 1.732050807568877*bz[3]*c2*dx0-1.732050807568877*bx[1]*c2*dx2; 
  outEy[6] += 1.732050807568877*ph[3]*c2chi*dx1-1.732050807568877*bx[2]*c2*dx2; 
  outEy[7] += 3.872983346207417*bz[1]*c2*dx0; 
  outEy[8] += 3.872983346207417*ph[2]*c2chi*dx1; 
  outEy[9] += -3.872983346207417*bx[3]*c2*dx2; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*ph[0]*c2chi*dx2; 
  outEz[4] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 
  outEz[5] += 1.732050807568877*ph[1]*c2chi*dx2-1.732050807568877*by[3]*c2*dx0; 
  outEz[6] += 1.732050807568877*ph[2]*c2chi*dx2+1.732050807568877*bx[3]*c2*dx1; 
  outEz[7] += -3.872983346207417*by[1]*c2*dx0; 
  outEz[8] += 3.872983346207417*bx[2]*c2*dx1; 
  outEz[9] += 3.872983346207417*ph[3]*c2chi*dx2; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += -1.732050807568877*ey[0]*dx2; 
  outBx[4] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 
  outBx[5] += 1.732050807568877*ps[3]*dx0*gamma-1.732050807568877*ey[1]*dx2; 
  outBx[6] += 1.732050807568877*ez[3]*dx1-1.732050807568877*ey[2]*dx2; 
  outBx[7] += 3.872983346207417*ps[1]*dx0*gamma; 
  outBx[8] += 3.872983346207417*ez[2]*dx1; 
  outBx[9] += -3.872983346207417*ey[3]*dx2; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ex[0]*dx2; 
  outBy[4] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 
  outBy[5] += 1.732050807568877*ex[1]*dx2-1.732050807568877*ez[3]*dx0; 
  outBy[6] += 1.732050807568877*ps[3]*dx1*gamma+1.732050807568877*ex[2]*dx2; 
  outBy[7] += -3.872983346207417*ez[1]*dx0; 
  outBy[8] += 3.872983346207417*ps[2]*dx1*gamma; 
  outBy[9] += 3.872983346207417*ex[3]*dx2; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ps[0]*dx2*gamma; 
  outBz[4] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 
  outBz[5] += 1.732050807568877*ps[1]*dx2*gamma+1.732050807568877*ey[3]*dx0; 
  outBz[6] += 1.732050807568877*ps[2]*dx2*gamma-1.732050807568877*ex[3]*dx1; 
  outBz[7] += 3.872983346207417*ey[1]*dx0; 
  outBz[8] += -3.872983346207417*ex[2]*dx1; 
  outBz[9] += 3.872983346207417*ps[3]*dx2*gamma; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ez[0]*chi*dx2; 
  outPh[4] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 
  outPh[5] += 1.732050807568877*ez[1]*chi*dx2+1.732050807568877*ex[3]*chi*dx0; 
  outPh[6] += 1.732050807568877*ez[2]*chi*dx2+1.732050807568877*ey[3]*chi*dx1; 
  outPh[7] += 3.872983346207417*ex[1]*chi*dx0; 
  outPh[8] += 3.872983346207417*ey[2]*chi*dx1; 
  outPh[9] += 3.872983346207417*ez[3]*chi*dx2; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*bz[0]*c2gamma*dx2; 
  outPs[4] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 
  outPs[5] += 1.732050807568877*bz[1]*c2gamma*dx2+1.732050807568877*bx[3]*c2gamma*dx0; 
  outPs[6] += 1.732050807568877*bz[2]*c2gamma*dx2+1.732050807568877*by[3]*c2gamma*dx1; 
  outPs[7] += 3.872983346207417*bx[1]*c2gamma*dx0; 
  outPs[8] += 3.872983346207417*by[2]*c2gamma*dx1; 
  outPs[9] += 3.872983346207417*bz[3]*c2gamma*dx2; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  cflFreq += meq->c/dx[1]; 
  cflFreq += meq->c/dx[2]; 
  return cflFreq; 
} 
double MaxwellVol3xMaxP3(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[20]; 
  const double *ez = &q[40]; 
  const double *bx = &q[60]; 
  const double *by = &q[80]; 
  const double *bz = &q[100]; 
  const double *ph = &q[120]; 
  const double *ps = &q[140]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[20]; 
  double *outEz = &out[40]; 
  double *outBx = &out[60]; 
  double *outBy = &out[80]; 
  double *outBz = &out[100]; 
  double *outPh = &out[120]; 
  double *outPs = &out[140]; 
 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 
  double dx2 = 2.0/dx[2]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*by[0]*c2*dx2; 
  outEx[4] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 
  outEx[5] += 1.732050807568877*by[1]*c2*dx2+1.732050807568877*ph[3]*c2chi*dx0; 
  outEx[6] += 1.732050807568877*by[2]*c2*dx2-1.732050807568877*bz[3]*c2*dx1; 
  outEx[7] += 3.872983346207417*ph[1]*c2chi*dx0; 
  outEx[8] += -3.872983346207417*bz[2]*c2*dx1; 
  outEx[9] += 3.872983346207417*by[3]*c2*dx2; 
  outEx[10] += 1.732050807568877*by[4]*c2*dx2-1.732050807568877*bz[5]*c2*dx1+1.732050807568877*ph[6]*c2chi*dx0; 
  outEx[11] += 3.872983346207417*ph[4]*c2chi*dx0-1.732050807568877*bz[7]*c2*dx1; 
  outEx[12] += 1.732050807568877*ph[8]*c2chi*dx0-3.872983346207417*bz[4]*c2*dx1; 
  outEx[13] += 1.732050807568877*by[7]*c2*dx2+3.872983346207417*ph[5]*c2chi*dx0; 
  outEx[14] += 1.732050807568877*by[8]*c2*dx2-3.872983346207417*bz[6]*c2*dx1; 
  outEx[15] += 3.872983346207417*by[5]*c2*dx2+1.732050807568877*ph[9]*c2chi*dx0; 
  outEx[16] += 3.872983346207417*by[6]*c2*dx2-1.732050807568877*bz[9]*c2*dx1; 
  outEx[17] += 5.916079783099617*ph[7]*c2chi*dx0+2.645751311064591*ph[0]*c2chi*dx0; 
  outEx[18] += (-5.916079783099617*bz[8]*c2*dx1)-2.645751311064591*bz[0]*c2*dx1; 
  outEx[19] += 5.916079783099617*by[9]*c2*dx2+2.645751311064591*by[0]*c2*dx2; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += -1.732050807568877*bx[0]*c2*dx2; 
  outEy[4] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 
  outEy[5] += 1.732050807568877*bz[3]*c2*dx0-1.732050807568877*bx[1]*c2*dx2; 
  outEy[6] += 1.732050807568877*ph[3]*c2chi*dx1-1.732050807568877*bx[2]*c2*dx2; 
  outEy[7] += 3.872983346207417*bz[1]*c2*dx0; 
  outEy[8] += 3.872983346207417*ph[2]*c2chi*dx1; 
  outEy[9] += -3.872983346207417*bx[3]*c2*dx2; 
  outEy[10] += (-1.732050807568877*bx[4]*c2*dx2)+1.732050807568877*ph[5]*c2chi*dx1+1.732050807568877*bz[6]*c2*dx0; 
  outEy[11] += 1.732050807568877*ph[7]*c2chi*dx1+3.872983346207417*bz[4]*c2*dx0; 
  outEy[12] += 3.872983346207417*ph[4]*c2chi*dx1+1.732050807568877*bz[8]*c2*dx0; 
  outEy[13] += 3.872983346207417*bz[5]*c2*dx0-1.732050807568877*bx[7]*c2*dx2; 
  outEy[14] += 3.872983346207417*ph[6]*c2chi*dx1-1.732050807568877*bx[8]*c2*dx2; 
  outEy[15] += 1.732050807568877*bz[9]*c2*dx0-3.872983346207417*bx[5]*c2*dx2; 
  outEy[16] += 1.732050807568877*ph[9]*c2chi*dx1-3.872983346207417*bx[6]*c2*dx2; 
  outEy[17] += 5.916079783099617*bz[7]*c2*dx0+2.645751311064591*bz[0]*c2*dx0; 
  outEy[18] += 5.916079783099617*ph[8]*c2chi*dx1+2.645751311064591*ph[0]*c2chi*dx1; 
  outEy[19] += (-5.916079783099617*bx[9]*c2*dx2)-2.645751311064591*bx[0]*c2*dx2; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*ph[0]*c2chi*dx2; 
  outEz[4] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 
  outEz[5] += 1.732050807568877*ph[1]*c2chi*dx2-1.732050807568877*by[3]*c2*dx0; 
  outEz[6] += 1.732050807568877*ph[2]*c2chi*dx2+1.732050807568877*bx[3]*c2*dx1; 
  outEz[7] += -3.872983346207417*by[1]*c2*dx0; 
  outEz[8] += 3.872983346207417*bx[2]*c2*dx1; 
  outEz[9] += 3.872983346207417*ph[3]*c2chi*dx2; 
  outEz[10] += 1.732050807568877*ph[4]*c2chi*dx2+1.732050807568877*bx[5]*c2*dx1-1.732050807568877*by[6]*c2*dx0; 
  outEz[11] += 1.732050807568877*bx[7]*c2*dx1-3.872983346207417*by[4]*c2*dx0; 
  outEz[12] += 3.872983346207417*bx[4]*c2*dx1-1.732050807568877*by[8]*c2*dx0; 
  outEz[13] += 1.732050807568877*ph[7]*c2chi*dx2-3.872983346207417*by[5]*c2*dx0; 
  outEz[14] += 1.732050807568877*ph[8]*c2chi*dx2+3.872983346207417*bx[6]*c2*dx1; 
  outEz[15] += 3.872983346207417*ph[5]*c2chi*dx2-1.732050807568877*by[9]*c2*dx0; 
  outEz[16] += 3.872983346207417*ph[6]*c2chi*dx2+1.732050807568877*bx[9]*c2*dx1; 
  outEz[17] += (-5.916079783099617*by[7]*c2*dx0)-2.645751311064591*by[0]*c2*dx0; 
  outEz[18] += 5.916079783099617*bx[8]*c2*dx1+2.645751311064591*bx[0]*c2*dx1; 
  outEz[19] += 5.916079783099617*ph[9]*c2chi*dx2+2.645751311064591*ph[0]*c2chi*dx2; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += -1.732050807568877*ey[0]*dx2; 
  outBx[4] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 
  outBx[5] += 1.732050807568877*ps[3]*dx0*gamma-1.732050807568877*ey[1]*dx2; 
  outBx[6] += 1.732050807568877*ez[3]*dx1-1.732050807568877*ey[2]*dx2; 
  outBx[7] += 3.872983346207417*ps[1]*dx0*gamma; 
  outBx[8] += 3.872983346207417*ez[2]*dx1; 
  outBx[9] += -3.872983346207417*ey[3]*dx2; 
  outBx[10] += 1.732050807568877*ps[6]*dx0*gamma-1.732050807568877*ey[4]*dx2+1.732050807568877*ez[5]*dx1; 
  outBx[11] += 3.872983346207417*ps[4]*dx0*gamma+1.732050807568877*ez[7]*dx1; 
  outBx[12] += 1.732050807568877*ps[8]*dx0*gamma+3.872983346207417*ez[4]*dx1; 
  outBx[13] += 3.872983346207417*ps[5]*dx0*gamma-1.732050807568877*ey[7]*dx2; 
  outBx[14] += 3.872983346207417*ez[6]*dx1-1.732050807568877*ey[8]*dx2; 
  outBx[15] += 1.732050807568877*ps[9]*dx0*gamma-3.872983346207417*ey[5]*dx2; 
  outBx[16] += 1.732050807568877*ez[9]*dx1-3.872983346207417*ey[6]*dx2; 
  outBx[17] += 5.916079783099617*ps[7]*dx0*gamma+2.645751311064591*ps[0]*dx0*gamma; 
  outBx[18] += 5.916079783099617*ez[8]*dx1+2.645751311064591*ez[0]*dx1; 
  outBx[19] += (-5.916079783099617*ey[9]*dx2)-2.645751311064591*ey[0]*dx2; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ex[0]*dx2; 
  outBy[4] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 
  outBy[5] += 1.732050807568877*ex[1]*dx2-1.732050807568877*ez[3]*dx0; 
  outBy[6] += 1.732050807568877*ps[3]*dx1*gamma+1.732050807568877*ex[2]*dx2; 
  outBy[7] += -3.872983346207417*ez[1]*dx0; 
  outBy[8] += 3.872983346207417*ps[2]*dx1*gamma; 
  outBy[9] += 3.872983346207417*ex[3]*dx2; 
  outBy[10] += 1.732050807568877*ps[5]*dx1*gamma+1.732050807568877*ex[4]*dx2-1.732050807568877*ez[6]*dx0; 
  outBy[11] += 1.732050807568877*ps[7]*dx1*gamma-3.872983346207417*ez[4]*dx0; 
  outBy[12] += 3.872983346207417*ps[4]*dx1*gamma-1.732050807568877*ez[8]*dx0; 
  outBy[13] += 1.732050807568877*ex[7]*dx2-3.872983346207417*ez[5]*dx0; 
  outBy[14] += 3.872983346207417*ps[6]*dx1*gamma+1.732050807568877*ex[8]*dx2; 
  outBy[15] += 3.872983346207417*ex[5]*dx2-1.732050807568877*ez[9]*dx0; 
  outBy[16] += 1.732050807568877*ps[9]*dx1*gamma+3.872983346207417*ex[6]*dx2; 
  outBy[17] += (-5.916079783099617*ez[7]*dx0)-2.645751311064591*ez[0]*dx0; 
  outBy[18] += 5.916079783099617*ps[8]*dx1*gamma+2.645751311064591*ps[0]*dx1*gamma; 
  outBy[19] += 5.916079783099617*ex[9]*dx2+2.645751311064591*ex[0]*dx2; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ps[0]*dx2*gamma; 
  outBz[4] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 
  outBz[5] += 1.732050807568877*ps[1]*dx2*gamma+1.732050807568877*ey[3]*dx0; 
  outBz[6] += 1.732050807568877*ps[2]*dx2*gamma-1.732050807568877*ex[3]*dx1; 
  outBz[7] += 3.872983346207417*ey[1]*dx0; 
  outBz[8] += -3.872983346207417*ex[2]*dx1; 
  outBz[9] += 3.872983346207417*ps[3]*dx2*gamma; 
  outBz[10] += 1.732050807568877*ps[4]*dx2*gamma-1.732050807568877*ex[5]*dx1+1.732050807568877*ey[6]*dx0; 
  outBz[11] += 3.872983346207417*ey[4]*dx0-1.732050807568877*ex[7]*dx1; 
  outBz[12] += 1.732050807568877*ey[8]*dx0-3.872983346207417*ex[4]*dx1; 
  outBz[13] += 1.732050807568877*ps[7]*dx2*gamma+3.872983346207417*ey[5]*dx0; 
  outBz[14] += 1.732050807568877*ps[8]*dx2*gamma-3.872983346207417*ex[6]*dx1; 
  outBz[15] += 3.872983346207417*ps[5]*dx2*gamma+1.732050807568877*ey[9]*dx0; 
  outBz[16] += 3.872983346207417*ps[6]*dx2*gamma-1.732050807568877*ex[9]*dx1; 
  outBz[17] += 5.916079783099617*ey[7]*dx0+2.645751311064591*ey[0]*dx0; 
  outBz[18] += (-5.916079783099617*ex[8]*dx1)-2.645751311064591*ex[0]*dx1; 
  outBz[19] += 5.916079783099617*ps[9]*dx2*gamma+2.645751311064591*ps[0]*dx2*gamma; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ez[0]*chi*dx2; 
  outPh[4] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 
  outPh[5] += 1.732050807568877*ez[1]*chi*dx2+1.732050807568877*ex[3]*chi*dx0; 
  outPh[6] += 1.732050807568877*ez[2]*chi*dx2+1.732050807568877*ey[3]*chi*dx1; 
  outPh[7] += 3.872983346207417*ex[1]*chi*dx0; 
  outPh[8] += 3.872983346207417*ey[2]*chi*dx1; 
  outPh[9] += 3.872983346207417*ez[3]*chi*dx2; 
  outPh[10] += 1.732050807568877*ez[4]*chi*dx2+1.732050807568877*ey[5]*chi*dx1+1.732050807568877*ex[6]*chi*dx0; 
  outPh[11] += 1.732050807568877*ey[7]*chi*dx1+3.872983346207417*ex[4]*chi*dx0; 
  outPh[12] += 3.872983346207417*ey[4]*chi*dx1+1.732050807568877*ex[8]*chi*dx0; 
  outPh[13] += 1.732050807568877*ez[7]*chi*dx2+3.872983346207417*ex[5]*chi*dx0; 
  outPh[14] += 1.732050807568877*ez[8]*chi*dx2+3.872983346207417*ey[6]*chi*dx1; 
  outPh[15] += 3.872983346207417*ez[5]*chi*dx2+1.732050807568877*ex[9]*chi*dx0; 
  outPh[16] += 3.872983346207417*ez[6]*chi*dx2+1.732050807568877*ey[9]*chi*dx1; 
  outPh[17] += 5.916079783099617*ex[7]*chi*dx0+2.645751311064591*ex[0]*chi*dx0; 
  outPh[18] += 5.916079783099617*ey[8]*chi*dx1+2.645751311064591*ey[0]*chi*dx1; 
  outPh[19] += 5.916079783099617*ez[9]*chi*dx2+2.645751311064591*ez[0]*chi*dx2; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*bz[0]*c2gamma*dx2; 
  outPs[4] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 
  outPs[5] += 1.732050807568877*bz[1]*c2gamma*dx2+1.732050807568877*bx[3]*c2gamma*dx0; 
  outPs[6] += 1.732050807568877*bz[2]*c2gamma*dx2+1.732050807568877*by[3]*c2gamma*dx1; 
  outPs[7] += 3.872983346207417*bx[1]*c2gamma*dx0; 
  outPs[8] += 3.872983346207417*by[2]*c2gamma*dx1; 
  outPs[9] += 3.872983346207417*bz[3]*c2gamma*dx2; 
  outPs[10] += 1.732050807568877*bz[4]*c2gamma*dx2+1.732050807568877*by[5]*c2gamma*dx1+1.732050807568877*bx[6]*c2gamma*dx0; 
  outPs[11] += 1.732050807568877*by[7]*c2gamma*dx1+3.872983346207417*bx[4]*c2gamma*dx0; 
  outPs[12] += 3.872983346207417*by[4]*c2gamma*dx1+1.732050807568877*bx[8]*c2gamma*dx0; 
  outPs[13] += 1.732050807568877*bz[7]*c2gamma*dx2+3.872983346207417*bx[5]*c2gamma*dx0; 
  outPs[14] += 1.732050807568877*bz[8]*c2gamma*dx2+3.872983346207417*by[6]*c2gamma*dx1; 
  outPs[15] += 3.872983346207417*bz[5]*c2gamma*dx2+1.732050807568877*bx[9]*c2gamma*dx0; 
  outPs[16] += 3.872983346207417*bz[6]*c2gamma*dx2+1.732050807568877*by[9]*c2gamma*dx1; 
  outPs[17] += 5.916079783099617*bx[7]*c2gamma*dx0+2.645751311064591*bx[0]*c2gamma*dx0; 
  outPs[18] += 5.916079783099617*by[8]*c2gamma*dx1+2.645751311064591*by[0]*c2gamma*dx1; 
  outPs[19] += 5.916079783099617*bz[9]*c2gamma*dx2+2.645751311064591*bz[0]*c2gamma*dx2; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  cflFreq += meq->c/dx[1]; 
  cflFreq += meq->c/dx[2]; 
  return cflFreq; 
} 
