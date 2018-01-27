#include <MaxwellModDecl.h> 
void MaxwellVol2xMaxP0(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[1]; 
  const double *ez = &q[2]; 
  const double *bx = &q[3]; 
  const double *by = &q[4]; 
  const double *bz = &q[5]; 
  const double *ph = &q[6]; 
  const double *ps = &q[7]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[1]; 
  double *outEz = &out[2]; 
  double *outBx = &out[3]; 
  double *outBy = &out[4]; 
  double *outBz = &out[5]; 
  double *outPh = &out[6]; 
  double *outPs = &out[7]; 
 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 

 








} 
void MaxwellVol2xMaxP1(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
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
  double dx1 = 2.0/dx[1]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 

} 
void MaxwellVol2xMaxP2(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[6]; 
  const double *ez = &q[12]; 
  const double *bx = &q[18]; 
  const double *by = &q[24]; 
  const double *bz = &q[30]; 
  const double *ph = &q[36]; 
  const double *ps = &q[42]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[6]; 
  double *outEz = &out[12]; 
  double *outBx = &out[18]; 
  double *outBy = &out[24]; 
  double *outBz = &out[30]; 
  double *outPh = &out[36]; 
  double *outPs = &out[42]; 
 
  double dx0 = 2.0/dx[0]; 
  double dx1 = 2.0/dx[1]; 

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 
  outEx[4] += 3.872983346207417*ph[1]*c2chi*dx0; 
  outEx[5] += -3.872983346207417*bz[2]*c2*dx1; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 
  outEy[4] += 3.872983346207417*bz[1]*c2*dx0; 
  outEy[5] += 3.872983346207417*ph[2]*c2chi*dx1; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 
  outEz[4] += -3.872983346207417*by[1]*c2*dx0; 
  outEz[5] += 3.872983346207417*bx[2]*c2*dx1; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 
  outBx[4] += 3.872983346207417*ps[1]*dx0*gamma; 
  outBx[5] += 3.872983346207417*ez[2]*dx1; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 
  outBy[4] += -3.872983346207417*ez[1]*dx0; 
  outBy[5] += 3.872983346207417*ps[2]*dx1*gamma; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 
  outBz[4] += 3.872983346207417*ey[1]*dx0; 
  outBz[5] += -3.872983346207417*ex[2]*dx1; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 
  outPh[4] += 3.872983346207417*ex[1]*chi*dx0; 
  outPh[5] += 3.872983346207417*ey[2]*chi*dx1; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 
  outPs[4] += 3.872983346207417*bx[1]*c2gamma*dx0; 
  outPs[5] += 3.872983346207417*by[2]*c2gamma*dx1; 

} 
void MaxwellVol2xMaxP3(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
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

 
  outEx[1] += 1.732050807568877*ph[0]*c2chi*dx0; 
  outEx[2] += -1.732050807568877*bz[0]*c2*dx1; 
  outEx[3] += 1.732050807568877*ph[2]*c2chi*dx0-1.732050807568877*bz[1]*c2*dx1; 
  outEx[4] += 3.872983346207417*ph[1]*c2chi*dx0; 
  outEx[5] += -3.872983346207417*bz[2]*c2*dx1; 
  outEx[6] += 3.872983346207417*ph[3]*c2chi*dx0-1.732050807568877*bz[4]*c2*dx1; 
  outEx[7] += 1.732050807568877*ph[5]*c2chi*dx0-3.872983346207417*bz[3]*c2*dx1; 
  outEx[8] += 5.916079783099617*ph[4]*c2chi*dx0+2.645751311064591*ph[0]*c2chi*dx0; 
  outEx[9] += (-5.916079783099617*bz[5]*c2*dx1)-2.645751311064591*bz[0]*c2*dx1; 

  outEy[1] += 1.732050807568877*bz[0]*c2*dx0; 
  outEy[2] += 1.732050807568877*ph[0]*c2chi*dx1; 
  outEy[3] += 1.732050807568877*ph[1]*c2chi*dx1+1.732050807568877*bz[2]*c2*dx0; 
  outEy[4] += 3.872983346207417*bz[1]*c2*dx0; 
  outEy[5] += 3.872983346207417*ph[2]*c2chi*dx1; 
  outEy[6] += 1.732050807568877*ph[4]*c2chi*dx1+3.872983346207417*bz[3]*c2*dx0; 
  outEy[7] += 3.872983346207417*ph[3]*c2chi*dx1+1.732050807568877*bz[5]*c2*dx0; 
  outEy[8] += 5.916079783099617*bz[4]*c2*dx0+2.645751311064591*bz[0]*c2*dx0; 
  outEy[9] += 5.916079783099617*ph[5]*c2chi*dx1+2.645751311064591*ph[0]*c2chi*dx1; 

  outEz[1] += -1.732050807568877*by[0]*c2*dx0; 
  outEz[2] += 1.732050807568877*bx[0]*c2*dx1; 
  outEz[3] += 1.732050807568877*bx[1]*c2*dx1-1.732050807568877*by[2]*c2*dx0; 
  outEz[4] += -3.872983346207417*by[1]*c2*dx0; 
  outEz[5] += 3.872983346207417*bx[2]*c2*dx1; 
  outEz[6] += 1.732050807568877*bx[4]*c2*dx1-3.872983346207417*by[3]*c2*dx0; 
  outEz[7] += 3.872983346207417*bx[3]*c2*dx1-1.732050807568877*by[5]*c2*dx0; 
  outEz[8] += (-5.916079783099617*by[4]*c2*dx0)-2.645751311064591*by[0]*c2*dx0; 
  outEz[9] += 5.916079783099617*bx[5]*c2*dx1+2.645751311064591*bx[0]*c2*dx1; 

  outBx[1] += 1.732050807568877*ps[0]*dx0*gamma; 
  outBx[2] += 1.732050807568877*ez[0]*dx1; 
  outBx[3] += 1.732050807568877*ps[2]*dx0*gamma+1.732050807568877*ez[1]*dx1; 
  outBx[4] += 3.872983346207417*ps[1]*dx0*gamma; 
  outBx[5] += 3.872983346207417*ez[2]*dx1; 
  outBx[6] += 3.872983346207417*ps[3]*dx0*gamma+1.732050807568877*ez[4]*dx1; 
  outBx[7] += 1.732050807568877*ps[5]*dx0*gamma+3.872983346207417*ez[3]*dx1; 
  outBx[8] += 5.916079783099617*ps[4]*dx0*gamma+2.645751311064591*ps[0]*dx0*gamma; 
  outBx[9] += 5.916079783099617*ez[5]*dx1+2.645751311064591*ez[0]*dx1; 

  outBy[1] += -1.732050807568877*ez[0]*dx0; 
  outBy[2] += 1.732050807568877*ps[0]*dx1*gamma; 
  outBy[3] += 1.732050807568877*ps[1]*dx1*gamma-1.732050807568877*ez[2]*dx0; 
  outBy[4] += -3.872983346207417*ez[1]*dx0; 
  outBy[5] += 3.872983346207417*ps[2]*dx1*gamma; 
  outBy[6] += 1.732050807568877*ps[4]*dx1*gamma-3.872983346207417*ez[3]*dx0; 
  outBy[7] += 3.872983346207417*ps[3]*dx1*gamma-1.732050807568877*ez[5]*dx0; 
  outBy[8] += (-5.916079783099617*ez[4]*dx0)-2.645751311064591*ez[0]*dx0; 
  outBy[9] += 5.916079783099617*ps[5]*dx1*gamma+2.645751311064591*ps[0]*dx1*gamma; 

  outBz[1] += 1.732050807568877*ey[0]*dx0; 
  outBz[2] += -1.732050807568877*ex[0]*dx1; 
  outBz[3] += 1.732050807568877*ey[2]*dx0-1.732050807568877*ex[1]*dx1; 
  outBz[4] += 3.872983346207417*ey[1]*dx0; 
  outBz[5] += -3.872983346207417*ex[2]*dx1; 
  outBz[6] += 3.872983346207417*ey[3]*dx0-1.732050807568877*ex[4]*dx1; 
  outBz[7] += 1.732050807568877*ey[5]*dx0-3.872983346207417*ex[3]*dx1; 
  outBz[8] += 5.916079783099617*ey[4]*dx0+2.645751311064591*ey[0]*dx0; 
  outBz[9] += (-5.916079783099617*ex[5]*dx1)-2.645751311064591*ex[0]*dx1; 

  outPh[1] += 1.732050807568877*ex[0]*chi*dx0; 
  outPh[2] += 1.732050807568877*ey[0]*chi*dx1; 
  outPh[3] += 1.732050807568877*ey[1]*chi*dx1+1.732050807568877*ex[2]*chi*dx0; 
  outPh[4] += 3.872983346207417*ex[1]*chi*dx0; 
  outPh[5] += 3.872983346207417*ey[2]*chi*dx1; 
  outPh[6] += 1.732050807568877*ey[4]*chi*dx1+3.872983346207417*ex[3]*chi*dx0; 
  outPh[7] += 3.872983346207417*ey[3]*chi*dx1+1.732050807568877*ex[5]*chi*dx0; 
  outPh[8] += 5.916079783099617*ex[4]*chi*dx0+2.645751311064591*ex[0]*chi*dx0; 
  outPh[9] += 5.916079783099617*ey[5]*chi*dx1+2.645751311064591*ey[0]*chi*dx1; 

  outPs[1] += 1.732050807568877*bx[0]*c2gamma*dx0; 
  outPs[2] += 1.732050807568877*by[0]*c2gamma*dx1; 
  outPs[3] += 1.732050807568877*by[1]*c2gamma*dx1+1.732050807568877*bx[2]*c2gamma*dx0; 
  outPs[4] += 3.872983346207417*bx[1]*c2gamma*dx0; 
  outPs[5] += 3.872983346207417*by[2]*c2gamma*dx1; 
  outPs[6] += 1.732050807568877*by[4]*c2gamma*dx1+3.872983346207417*bx[3]*c2gamma*dx0; 
  outPs[7] += 3.872983346207417*by[3]*c2gamma*dx1+1.732050807568877*bx[5]*c2gamma*dx0; 
  outPs[8] += 5.916079783099617*bx[4]*c2gamma*dx0+2.645751311064591*bx[0]*c2gamma*dx0; 
  outPs[9] += 5.916079783099617*by[5]*c2gamma*dx1+2.645751311064591*by[0]*c2gamma*dx1; 

} 
void MaxwellVol2xMaxP4(const MaxwellEq_t * const meq, const double *w, const double *dx, const double *q, double *out) 
{ 
  const double c2 = meq->c*meq->c, chi = meq->chi, gamma = meq->gamma; 
  const double c2chi = c2*chi, c2gamma = c2*gamma; 
 
  const double *ex = &q[0]; 
  const double *ey = &q[15]; 
  const double *ez = &q[30]; 
  const double *bx = &q[45]; 
  const double *by = &q[60]; 
  const double *bz = &q[75]; 
  const double *ph = &q[90]; 
  const double *ps = &q[105]; 
 
  double *outEx = &out[0]; 
  double *outEy = &out[15]; 
  double *outEz = &out[30]; 
  double *outBx = &out[45]; 
  double *outBy = &out[60]; 
  double *outBz = &out[75]; 
  double *outPh = &out[90]; 
  double *outPs = &out[105]; 
 
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
  outEx[13] += 7.937253933193772*ph[8]*c2chi*dx0+5.196152422706631*ph[1]*c2chi*dx0; 
  outEx[14] += (-7.937253933193772*bz[9]*c2*dx1)-5.196152422706631*bz[2]*c2*dx1; 

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
  outEy[13] += 7.937253933193772*bz[8]*c2*dx0+5.196152422706631*bz[1]*c2*dx0; 
  outEy[14] += 7.937253933193772*ph[9]*c2chi*dx1+5.196152422706631*ph[2]*c2chi*dx1; 

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
  outEz[13] += (-7.937253933193772*by[8]*c2*dx0)-5.196152422706631*by[1]*c2*dx0; 
  outEz[14] += 7.937253933193772*bx[9]*c2*dx1+5.196152422706631*bx[2]*c2*dx1; 

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
  outBx[13] += 7.937253933193772*ps[8]*dx0*gamma+5.196152422706631*ps[1]*dx0*gamma; 
  outBx[14] += 7.937253933193772*ez[9]*dx1+5.196152422706631*ez[2]*dx1; 

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
  outBy[13] += (-7.937253933193772*ez[8]*dx0)-5.196152422706631*ez[1]*dx0; 
  outBy[14] += 7.937253933193772*ps[9]*dx1*gamma+5.196152422706631*ps[2]*dx1*gamma; 

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
  outBz[13] += 7.937253933193772*ey[8]*dx0+5.196152422706631*ey[1]*dx0; 
  outBz[14] += (-7.937253933193772*ex[9]*dx1)-5.196152422706631*ex[2]*dx1; 

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
  outPh[13] += 7.937253933193772*ex[8]*chi*dx0+5.196152422706631*ex[1]*chi*dx0; 
  outPh[14] += 7.937253933193772*ey[9]*chi*dx1+5.196152422706631*ey[2]*chi*dx1; 

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
  outPs[13] += 7.937253933193772*bx[8]*c2gamma*dx0+5.196152422706631*bx[1]*c2gamma*dx0; 
  outPs[14] += 7.937253933193772*by[9]*c2gamma*dx1+5.196152422706631*by[2]*c2gamma*dx1; 

} 
