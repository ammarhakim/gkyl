#include <MaxwellModDecl.h> 
void MaxwellVol3xMaxP1(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*meq->chi, c2gamma = c2*meq->gamma; 
 
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
  outBx[2] += 1.732050807568877*ez[0]*c2*dx1; 
  outBx[3] += -1.732050807568877*ey[0]*c2*dx2; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*c2gamma*dx1; 
  outBy[3] += 4.898979485566357*c2*dx2; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*c2*dx1; 
  outBz[3] += 1.732050807568877*ex[0]*dx2; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ps[0]*c2gamma*dx2; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*ez[0]*chi*dx2; 

} 
void MaxwellVol3xMaxP2(const MaxwellEq_t *meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*meq->chi, c2gamma = c2*meq->gamma; 
 
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
  outBx[2] += 1.732050807568877*ez[0]*c2*dx1; 
  outBx[3] += -1.732050807568877*ey[0]*c2*dx2; 
  outBx[4] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*c2*dx1; 
  outBx[5] += 1.732050807568877*ps[3]*dx0*gamma-1.732050807568877*ey[1]*c2*dx2; 
  outBx[6] += 1.732050807568877*ez[3]*c2*dx1-1.732050807568877*ey[2]*c2*dx2; 
  outBx[7] += 3.872983346207417*ps[1]*dx0*gamma; 
  outBx[8] += 3.872983346207417*ez[2]*c2*dx1; 
  outBx[9] += -3.872983346207417*ey[3]*c2*dx2; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*c2gamma*dx1; 
  outBy[3] += 4.898979485566357*c2*dx2; 
  outBy[4] += 1.732050807568877*ps[1]*c2gamma*dx1-1.732050807568877*ez[2]*dx0; 
  outBy[5] += -1.732050807568877*ez[3]*dx0; 
  outBy[6] += 1.732050807568877*ps[3]*c2gamma*dx1; 
  outBy[7] += -3.872983346207417*ez[1]*dx0; 
  outBy[8] += 3.872983346207417*ps[2]*c2gamma*dx1; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*c2*dx1; 
  outBz[3] += 1.732050807568877*ex[0]*dx2; 
  outBz[4] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*c2*dx1; 
  outBz[5] += 1.732050807568877*ex[1]*dx2+1.732050807568877*ey[3]*dx0; 
  outBz[6] += 1.732050807568877*ex[2]*dx2-1.732050807568877*ex[3]*c2*dx1; 
  outBz[7] += 3.872983346207417*ey[1]*dx0; 
  outBz[8] += -3.872983346207417*ex[2]*c2*dx1; 
  outBz[9] += 3.872983346207417*ex[3]*dx2; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ps[0]*c2gamma*dx2; 
  outPh[4] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 
  outPh[5] += 1.732050807568877*ps[1]*c2gamma*dx2+1.732050807568877*ex[3]*chi*dx0; 
  outPh[6] += 1.732050807568877*ps[2]*c2gamma*dx2+1.732050807568877*ey[3]*chi*dx1; 
  outPh[7] += 3.872983346207417*ex[1]*chi*dx0; 
  outPh[8] += 3.872983346207417*ey[2]*chi*dx1; 
  outPh[9] += 3.872983346207417*ps[3]*c2gamma*dx2; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*ez[0]*chi*dx2; 
  outPs[4] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 
  outPs[5] += 1.732050807568877*ez[1]*chi*dx2+1.732050807568877*bx[3]*c2gamma*dx0; 
  outPs[6] += 1.732050807568877*ez[2]*chi*dx2+1.732050807568877*by[3]*c2gamma*dx1; 
  outPs[7] += 3.872983346207417*bx[1]*c2gamma*dx0; 
  outPs[8] += 3.872983346207417*by[2]*c2gamma*dx1; 
  outPs[9] += 3.872983346207417*ez[3]*chi*dx2; 

} 
