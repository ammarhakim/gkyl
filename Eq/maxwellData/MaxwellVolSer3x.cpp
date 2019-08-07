#include <MaxwellModDecl.h> 
double MaxwellVol3xSerP1(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[8]; 
  const double *ez = &q[16]; 
  const double *bx = &q[24]; 
  const double *by = &q[32]; 
  const double *bz = &q[40]; 
  const double *ph = &q[48]; 
  const double *ps = &q[56]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[8]; 
  double *outEz = &out[16]; 
  double *outBx = &out[24]; 
  double *outBy = &out[32]; 
  double *outBz = &out[40]; 
  double *outPh = &out[48]; 
  double *outPs = &out[56]; 
 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 
  double dx2 = 2.0/dx[2]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*by[0]*c2*dx2; 
  outEx[4] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 
  outEx[5] += 1.732050807568877*by[1]*c2*dx2+1.732050807568877*ph[3]*c2chi*dx0; 
  outEx[6] += 1.732050807568877*by[2]*c2*dx2-1.732050807568877*bz[3]*c2*dx1; 
  outEx[7] += 1.732050807568877*by[4]*c2*dx2-1.732050807568877*bz[5]*c2*dx1+1.732050807568877*ph[6]*c2chi*dx0; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += -1.732050807568877*bx[0]*c2*dx2; 
  outEy[4] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 
  outEy[5] += 1.732050807568877*bz[3]*c2*dx0-1.732050807568877*bx[1]*c2*dx2; 
  outEy[6] += 1.732050807568877*ph[3]*c2chi*dx1-1.732050807568877*bx[2]*c2*dx2; 
  outEy[7] += (-1.732050807568877*bx[4]*c2*dx2)+1.732050807568877*ph[5]*c2chi*dx1+1.732050807568877*bz[6]*c2*dx0; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*ph[0]*c2chi*dx2; 
  outEz[4] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 
  outEz[5] += 1.732050807568877*ph[1]*c2chi*dx2-1.732050807568877*by[3]*c2*dx0; 
  outEz[6] += 1.732050807568877*ph[2]*c2chi*dx2+1.732050807568877*bx[3]*c2*dx1; 
  outEz[7] += 1.732050807568877*ph[4]*c2chi*dx2+1.732050807568877*bx[5]*c2*dx1-1.732050807568877*by[6]*c2*dx0; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += -1.732050807568877*ey[0]*dx2; 
  outBx[4] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 
  outBx[5] += 1.732050807568877*ps[3]*dx0*gamma-1.732050807568877*ey[1]*dx2; 
  outBx[6] += 1.732050807568877*ez[3]*dx1-1.732050807568877*ey[2]*dx2; 
  outBx[7] += 1.732050807568877*ps[6]*dx0*gamma-1.732050807568877*ey[4]*dx2+1.732050807568877*ez[5]*dx1; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ex[0]*dx2; 
  outBy[4] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 
  outBy[5] += 1.732050807568877*ex[1]*dx2-1.732050807568877*ez[3]*dx0; 
  outBy[6] += 1.732050807568877*ps[3]*dx1*gamma+1.732050807568877*ex[2]*dx2; 
  outBy[7] += 1.732050807568877*ps[5]*dx1*gamma+1.732050807568877*ex[4]*dx2-1.732050807568877*ez[6]*dx0; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ps[0]*dx2*gamma; 
  outBz[4] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 
  outBz[5] += 1.732050807568877*ps[1]*dx2*gamma+1.732050807568877*ey[3]*dx0; 
  outBz[6] += 1.732050807568877*ps[2]*dx2*gamma-1.732050807568877*ex[3]*dx1; 
  outBz[7] += 1.732050807568877*ps[4]*dx2*gamma-1.732050807568877*ex[5]*dx1+1.732050807568877*ey[6]*dx0; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ez[0]*chi*dx2; 
  outPh[4] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 
  outPh[5] += 1.732050807568877*ez[1]*chi*dx2+1.732050807568877*ex[3]*chi*dx0; 
  outPh[6] += 1.732050807568877*ez[2]*chi*dx2+1.732050807568877*ey[3]*chi*dx1; 
  outPh[7] += 1.732050807568877*ez[4]*chi*dx2+1.732050807568877*ey[5]*chi*dx1+1.732050807568877*ex[6]*chi*dx0; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*bz[0]*c2gamma*dx2; 
  outPs[4] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 
  outPs[5] += 1.732050807568877*bz[1]*c2gamma*dx2+1.732050807568877*bx[3]*c2gamma*dx0; 
  outPs[6] += 1.732050807568877*bz[2]*c2gamma*dx2+1.732050807568877*by[3]*c2gamma*dx1; 
  outPs[7] += 1.732050807568877*bz[4]*c2gamma*dx2+1.732050807568877*by[5]*c2gamma*dx1+1.732050807568877*bx[6]*c2gamma*dx0; 

  double cflFreq = 0.0; 
  cflFreq += meq->c/dx[0]; 
  cflFreq += meq->c/dx[1]; 
  cflFreq += meq->c/dx[2]; 
  return cflFreq; 
} 
