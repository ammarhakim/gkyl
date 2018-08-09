#include <MaxwellModDecl.h> 
double MaxwellVol2xTensorP1(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
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

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  cflFreq += meq->c/dx[1]; 
  return cflFreq; 
} 
double MaxwellVol2xTensorP2(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[9]; 
  const double *ez = &q[18]; 
  const double *bx = &q[27]; 
  const double *by = &q[36]; 
  const double *bz = &q[45]; 
  const double *ph = &q[54]; 
  const double *ps = &q[63]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[9]; 
  double *outEz = &out[18]; 
  double *outBx = &out[27]; 
  double *outBy = &out[36]; 
  double *outBz = &out[45]; 
  double *outPh = &out[54]; 
  double *outPs = &out[63]; 
 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 
  outEx[4] += 3.872983346207417*ph[1]*c2chi*dx0; 
  outEx[5] += -3.872983346207417*bz[2]*c2*dx1; 
  outEx[6] += 3.872983346207417*ph[3]*c2chi*dx0-1.732050807568877*bz[4]*c2*dx1; 
  outEx[7] += 1.732050807568877*ph[5]*c2chi*dx0-3.872983346207417*bz[3]*c2*dx1; 
  outEx[8] += 3.872983346207417*ph[7]*c2chi*dx0-3.872983346207417*bz[6]*c2*dx1; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 
  outEy[4] += 3.872983346207417*bz[1]*c2*dx0; 
  outEy[5] += 3.872983346207417*ph[2]*c2chi*dx1; 
  outEy[6] += 1.732050807568877*ph[4]*c2chi*dx1+3.872983346207417*bz[3]*c2*dx0; 
  outEy[7] += 3.872983346207417*ph[3]*c2chi*dx1+1.732050807568877*bz[5]*c2*dx0; 
  outEy[8] += 3.872983346207417*ph[6]*c2chi*dx1+3.872983346207417*bz[7]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 
  outEz[4] += -3.872983346207417*by[1]*c2*dx0; 
  outEz[5] += 3.872983346207417*bx[2]*c2*dx1; 
  outEz[6] += 1.732050807568877*bx[4]*c2*dx1-3.872983346207417*by[3]*c2*dx0; 
  outEz[7] += 3.872983346207417*bx[3]*c2*dx1-1.732050807568877*by[5]*c2*dx0; 
  outEz[8] += 3.872983346207417*bx[6]*c2*dx1-3.872983346207417*by[7]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 
  outBx[4] += 3.872983346207417*ps[1]*dx0*gamma; 
  outBx[5] += 3.872983346207417*ez[2]*dx1; 
  outBx[6] += 3.872983346207417*ps[3]*dx0*gamma+1.732050807568877*ez[4]*dx1; 
  outBx[7] += 1.732050807568877*ps[5]*dx0*gamma+3.872983346207417*ez[3]*dx1; 
  outBx[8] += 3.872983346207417*ps[7]*dx0*gamma+3.872983346207417*ez[6]*dx1; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 
  outBy[4] += -3.872983346207417*ez[1]*dx0; 
  outBy[5] += 3.872983346207417*ps[2]*dx1*gamma; 
  outBy[6] += 1.732050807568877*ps[4]*dx1*gamma-3.872983346207417*ez[3]*dx0; 
  outBy[7] += 3.872983346207417*ps[3]*dx1*gamma-1.732050807568877*ez[5]*dx0; 
  outBy[8] += 3.872983346207417*ps[6]*dx1*gamma-3.872983346207417*ez[7]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 
  outBz[4] += 3.872983346207417*ey[1]*dx0; 
  outBz[5] += -3.872983346207417*ex[2]*dx1; 
  outBz[6] += 3.872983346207417*ey[3]*dx0-1.732050807568877*ex[4]*dx1; 
  outBz[7] += 1.732050807568877*ey[5]*dx0-3.872983346207417*ex[3]*dx1; 
  outBz[8] += 3.872983346207417*ey[7]*dx0-3.872983346207417*ex[6]*dx1; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 
  outPh[4] += 3.872983346207417*ex[1]*chi*dx0; 
  outPh[5] += 3.872983346207417*ey[2]*chi*dx1; 
  outPh[6] += 1.732050807568877*ey[4]*chi*dx1+3.872983346207417*ex[3]*chi*dx0; 
  outPh[7] += 3.872983346207417*ey[3]*chi*dx1+1.732050807568877*ex[5]*chi*dx0; 
  outPh[8] += 3.872983346207417*ey[6]*chi*dx1+3.872983346207417*ex[7]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 
  outPs[4] += 3.872983346207417*bx[1]*c2gamma*dx0; 
  outPs[5] += 3.872983346207417*by[2]*c2gamma*dx1; 
  outPs[6] += 1.732050807568877*by[4]*c2gamma*dx1+3.872983346207417*bx[3]*c2gamma*dx0; 
  outPs[7] += 3.872983346207417*by[3]*c2gamma*dx1+1.732050807568877*bx[5]*c2gamma*dx0; 
  outPs[8] += 3.872983346207417*by[6]*c2gamma*dx1+3.872983346207417*bx[7]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  cflFreq += meq->c/dx[1]; 
  return cflFreq; 
} 
double MaxwellVol2xTensorP3(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[16]; 
  const double *ez = &q[32]; 
  const double *bx = &q[48]; 
  const double *by = &q[64]; 
  const double *bz = &q[80]; 
  const double *ph = &q[96]; 
  const double *ps = &q[112]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[16]; 
  double *outEz = &out[32]; 
  double *outBx = &out[48]; 
  double *outBy = &out[64]; 
  double *outBz = &out[80]; 
  double *outPh = &out[96]; 
  double *outPs = &out[112]; 
 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 
  outEx[4] += 3.872983346207417*ph[1]*c2chi*dx0; 
  outEx[5] += -3.872983346207417*bz[2]*c2*dx1; 
  outEx[6] += 3.872983346207417*ph[3]*c2chi*dx0-1.732050807568877*bz[4]*c2*dx1; 
  outEx[7] += 1.732050807568877*ph[5]*c2chi*dx0-3.872983346207417*bz[3]*c2*dx1; 
  outEx[8] += 5.916079783099617*ph[4]*c2chi*dx0+2.645751311064591*ph[0]*c2chi*dx0; 
  outEx[9] += (-5.916079783099617*bz[5]*c2*dx1)-2.645751311064591*bz[0]*c2*dx1; 
  outEx[10] += 3.872983346207417*ph[7]*c2chi*dx0-3.872983346207417*bz[6]*c2*dx1; 
  outEx[11] += (-1.732050807568877*bz[8]*c2*dx1)+5.916079783099615*ph[6]*c2chi*dx0+2.645751311064591*ph[2]*c2chi*dx0; 
  outEx[12] += (-5.916079783099615*bz[7]*c2*dx1)-2.645751311064591*bz[1]*c2*dx1+1.732050807568877*ph[9]*c2chi*dx0; 
  outEx[13] += (-3.872983346207417*bz[11]*c2*dx1)+5.916079783099616*ph[10]*c2chi*dx0+2.645751311064591*ph[5]*c2chi*dx0; 
  outEx[14] += (-5.916079783099616*bz[10]*c2*dx1)-2.645751311064591*bz[4]*c2*dx1+3.872983346207417*ph[12]*c2chi*dx0; 
  outEx[15] += (-5.916079783099616*bz[13]*c2*dx1)-2.645751311064591*bz[8]*c2*dx1+5.916079783099616*ph[14]*c2chi*dx0+2.645751311064591*ph[9]*c2chi*dx0; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 
  outEy[4] += 3.872983346207417*bz[1]*c2*dx0; 
  outEy[5] += 3.872983346207417*ph[2]*c2chi*dx1; 
  outEy[6] += 1.732050807568877*ph[4]*c2chi*dx1+3.872983346207417*bz[3]*c2*dx0; 
  outEy[7] += 3.872983346207417*ph[3]*c2chi*dx1+1.732050807568877*bz[5]*c2*dx0; 
  outEy[8] += 5.916079783099617*bz[4]*c2*dx0+2.645751311064591*bz[0]*c2*dx0; 
  outEy[9] += 5.916079783099617*ph[5]*c2chi*dx1+2.645751311064591*ph[0]*c2chi*dx1; 
  outEy[10] += 3.872983346207417*ph[6]*c2chi*dx1+3.872983346207417*bz[7]*c2*dx0; 
  outEy[11] += 1.732050807568877*ph[8]*c2chi*dx1+5.916079783099615*bz[6]*c2*dx0+2.645751311064591*bz[2]*c2*dx0; 
  outEy[12] += 5.916079783099615*ph[7]*c2chi*dx1+2.645751311064591*ph[1]*c2chi*dx1+1.732050807568877*bz[9]*c2*dx0; 
  outEy[13] += 3.872983346207417*ph[11]*c2chi*dx1+5.916079783099616*bz[10]*c2*dx0+2.645751311064591*bz[5]*c2*dx0; 
  outEy[14] += 5.916079783099616*ph[10]*c2chi*dx1+2.645751311064591*ph[4]*c2chi*dx1+3.872983346207417*bz[12]*c2*dx0; 
  outEy[15] += 5.916079783099616*ph[13]*c2chi*dx1+2.645751311064591*ph[8]*c2chi*dx1+5.916079783099616*bz[14]*c2*dx0+2.645751311064591*bz[9]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 
  outEz[4] += -3.872983346207417*by[1]*c2*dx0; 
  outEz[5] += 3.872983346207417*bx[2]*c2*dx1; 
  outEz[6] += 1.732050807568877*bx[4]*c2*dx1-3.872983346207417*by[3]*c2*dx0; 
  outEz[7] += 3.872983346207417*bx[3]*c2*dx1-1.732050807568877*by[5]*c2*dx0; 
  outEz[8] += (-5.916079783099617*by[4]*c2*dx0)-2.645751311064591*by[0]*c2*dx0; 
  outEz[9] += 5.916079783099617*bx[5]*c2*dx1+2.645751311064591*bx[0]*c2*dx1; 
  outEz[10] += 3.872983346207417*bx[6]*c2*dx1-3.872983346207417*by[7]*c2*dx0; 
  outEz[11] += 1.732050807568877*bx[8]*c2*dx1-5.916079783099615*by[6]*c2*dx0-2.645751311064591*by[2]*c2*dx0; 
  outEz[12] += 5.916079783099615*bx[7]*c2*dx1+2.645751311064591*bx[1]*c2*dx1-1.732050807568877*by[9]*c2*dx0; 
  outEz[13] += 3.872983346207417*bx[11]*c2*dx1-5.916079783099616*by[10]*c2*dx0-2.645751311064591*by[5]*c2*dx0; 
  outEz[14] += 5.916079783099616*bx[10]*c2*dx1+2.645751311064591*bx[4]*c2*dx1-3.872983346207417*by[12]*c2*dx0; 
  outEz[15] += 5.916079783099616*bx[13]*c2*dx1+2.645751311064591*bx[8]*c2*dx1-5.916079783099616*by[14]*c2*dx0-2.645751311064591*by[9]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 
  outBx[4] += 3.872983346207417*ps[1]*dx0*gamma; 
  outBx[5] += 3.872983346207417*ez[2]*dx1; 
  outBx[6] += 3.872983346207417*ps[3]*dx0*gamma+1.732050807568877*ez[4]*dx1; 
  outBx[7] += 1.732050807568877*ps[5]*dx0*gamma+3.872983346207417*ez[3]*dx1; 
  outBx[8] += 5.916079783099617*ps[4]*dx0*gamma+2.645751311064591*ps[0]*dx0*gamma; 
  outBx[9] += 5.916079783099617*ez[5]*dx1+2.645751311064591*ez[0]*dx1; 
  outBx[10] += 3.872983346207417*ps[7]*dx0*gamma+3.872983346207417*ez[6]*dx1; 
  outBx[11] += 5.916079783099615*ps[6]*dx0*gamma+2.645751311064591*ps[2]*dx0*gamma+1.732050807568877*ez[8]*dx1; 
  outBx[12] += 1.732050807568877*ps[9]*dx0*gamma+5.916079783099615*ez[7]*dx1+2.645751311064591*ez[1]*dx1; 
  outBx[13] += 5.916079783099616*ps[10]*dx0*gamma+2.645751311064591*ps[5]*dx0*gamma+3.872983346207417*ez[11]*dx1; 
  outBx[14] += 3.872983346207417*ps[12]*dx0*gamma+5.916079783099616*ez[10]*dx1+2.645751311064591*ez[4]*dx1; 
  outBx[15] += 5.916079783099616*ps[14]*dx0*gamma+2.645751311064591*ps[9]*dx0*gamma+5.916079783099616*ez[13]*dx1+2.645751311064591*ez[8]*dx1; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 
  outBy[4] += -3.872983346207417*ez[1]*dx0; 
  outBy[5] += 3.872983346207417*ps[2]*dx1*gamma; 
  outBy[6] += 1.732050807568877*ps[4]*dx1*gamma-3.872983346207417*ez[3]*dx0; 
  outBy[7] += 3.872983346207417*ps[3]*dx1*gamma-1.732050807568877*ez[5]*dx0; 
  outBy[8] += (-5.916079783099617*ez[4]*dx0)-2.645751311064591*ez[0]*dx0; 
  outBy[9] += 5.916079783099617*ps[5]*dx1*gamma+2.645751311064591*ps[0]*dx1*gamma; 
  outBy[10] += 3.872983346207417*ps[6]*dx1*gamma-3.872983346207417*ez[7]*dx0; 
  outBy[11] += 1.732050807568877*ps[8]*dx1*gamma-5.916079783099615*ez[6]*dx0-2.645751311064591*ez[2]*dx0; 
  outBy[12] += 5.916079783099615*ps[7]*dx1*gamma+2.645751311064591*ps[1]*dx1*gamma-1.732050807568877*ez[9]*dx0; 
  outBy[13] += 3.872983346207417*ps[11]*dx1*gamma-5.916079783099616*ez[10]*dx0-2.645751311064591*ez[5]*dx0; 
  outBy[14] += 5.916079783099616*ps[10]*dx1*gamma+2.645751311064591*ps[4]*dx1*gamma-3.872983346207417*ez[12]*dx0; 
  outBy[15] += 5.916079783099616*ps[13]*dx1*gamma+2.645751311064591*ps[8]*dx1*gamma-5.916079783099616*ez[14]*dx0-2.645751311064591*ez[9]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 
  outBz[4] += 3.872983346207417*ey[1]*dx0; 
  outBz[5] += -3.872983346207417*ex[2]*dx1; 
  outBz[6] += 3.872983346207417*ey[3]*dx0-1.732050807568877*ex[4]*dx1; 
  outBz[7] += 1.732050807568877*ey[5]*dx0-3.872983346207417*ex[3]*dx1; 
  outBz[8] += 5.916079783099617*ey[4]*dx0+2.645751311064591*ey[0]*dx0; 
  outBz[9] += (-5.916079783099617*ex[5]*dx1)-2.645751311064591*ex[0]*dx1; 
  outBz[10] += 3.872983346207417*ey[7]*dx0-3.872983346207417*ex[6]*dx1; 
  outBz[11] += (-1.732050807568877*ex[8]*dx1)+5.916079783099615*ey[6]*dx0+2.645751311064591*ey[2]*dx0; 
  outBz[12] += (-5.916079783099615*ex[7]*dx1)-2.645751311064591*ex[1]*dx1+1.732050807568877*ey[9]*dx0; 
  outBz[13] += (-3.872983346207417*ex[11]*dx1)+5.916079783099616*ey[10]*dx0+2.645751311064591*ey[5]*dx0; 
  outBz[14] += (-5.916079783099616*ex[10]*dx1)-2.645751311064591*ex[4]*dx1+3.872983346207417*ey[12]*dx0; 
  outBz[15] += (-5.916079783099616*ex[13]*dx1)-2.645751311064591*ex[8]*dx1+5.916079783099616*ey[14]*dx0+2.645751311064591*ey[9]*dx0; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 
  outPh[4] += 3.872983346207417*ex[1]*chi*dx0; 
  outPh[5] += 3.872983346207417*ey[2]*chi*dx1; 
  outPh[6] += 1.732050807568877*ey[4]*chi*dx1+3.872983346207417*ex[3]*chi*dx0; 
  outPh[7] += 3.872983346207417*ey[3]*chi*dx1+1.732050807568877*ex[5]*chi*dx0; 
  outPh[8] += 5.916079783099617*ex[4]*chi*dx0+2.645751311064591*ex[0]*chi*dx0; 
  outPh[9] += 5.916079783099617*ey[5]*chi*dx1+2.645751311064591*ey[0]*chi*dx1; 
  outPh[10] += 3.872983346207417*ey[6]*chi*dx1+3.872983346207417*ex[7]*chi*dx0; 
  outPh[11] += 1.732050807568877*ey[8]*chi*dx1+5.916079783099615*ex[6]*chi*dx0+2.645751311064591*ex[2]*chi*dx0; 
  outPh[12] += 5.916079783099615*ey[7]*chi*dx1+2.645751311064591*ey[1]*chi*dx1+1.732050807568877*ex[9]*chi*dx0; 
  outPh[13] += 3.872983346207417*ey[11]*chi*dx1+5.916079783099616*ex[10]*chi*dx0+2.645751311064591*ex[5]*chi*dx0; 
  outPh[14] += 5.916079783099616*ey[10]*chi*dx1+2.645751311064591*ey[4]*chi*dx1+3.872983346207417*ex[12]*chi*dx0; 
  outPh[15] += 5.916079783099616*ey[13]*chi*dx1+2.645751311064591*ey[8]*chi*dx1+5.916079783099616*ex[14]*chi*dx0+2.645751311064591*ex[9]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 
  outPs[4] += 3.872983346207417*bx[1]*c2gamma*dx0; 
  outPs[5] += 3.872983346207417*by[2]*c2gamma*dx1; 
  outPs[6] += 1.732050807568877*by[4]*c2gamma*dx1+3.872983346207417*bx[3]*c2gamma*dx0; 
  outPs[7] += 3.872983346207417*by[3]*c2gamma*dx1+1.732050807568877*bx[5]*c2gamma*dx0; 
  outPs[8] += 5.916079783099617*bx[4]*c2gamma*dx0+2.645751311064591*bx[0]*c2gamma*dx0; 
  outPs[9] += 5.916079783099617*by[5]*c2gamma*dx1+2.645751311064591*by[0]*c2gamma*dx1; 
  outPs[10] += 3.872983346207417*by[6]*c2gamma*dx1+3.872983346207417*bx[7]*c2gamma*dx0; 
  outPs[11] += 1.732050807568877*by[8]*c2gamma*dx1+5.916079783099615*bx[6]*c2gamma*dx0+2.645751311064591*bx[2]*c2gamma*dx0; 
  outPs[12] += 5.916079783099615*by[7]*c2gamma*dx1+2.645751311064591*by[1]*c2gamma*dx1+1.732050807568877*bx[9]*c2gamma*dx0; 
  outPs[13] += 3.872983346207417*by[11]*c2gamma*dx1+5.916079783099616*bx[10]*c2gamma*dx0+2.645751311064591*bx[5]*c2gamma*dx0; 
  outPs[14] += 5.916079783099616*by[10]*c2gamma*dx1+2.645751311064591*by[4]*c2gamma*dx1+3.872983346207417*bx[12]*c2gamma*dx0; 
  outPs[15] += 5.916079783099616*by[13]*c2gamma*dx1+2.645751311064591*by[8]*c2gamma*dx1+5.916079783099616*bx[14]*c2gamma*dx0+2.645751311064591*bx[9]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  cflFreq += meq->c/dx[1]; 
  return cflFreq; 
} 
