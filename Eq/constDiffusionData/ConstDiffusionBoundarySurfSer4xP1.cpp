#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf4xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  if (idxr[0] == 1) {

  incr2[1] = 0.8660254037844386*fr[0]-1.5*fr[1]; 
  incr2[5] = 0.8660254037844386*fr[2]-1.5*fr[5]; 
  incr2[6] = 0.8660254037844386*fr[3]-1.5*fr[6]; 
  incr2[8] = 0.8660254037844386*fr[4]-1.5*fr[8]; 
  incr2[11] = 0.8660254037844386*fr[7]-1.5*fr[11]; 
  incr2[12] = 0.8660254037844386*fr[9]-1.5*fr[12]; 
  incr2[13] = 0.8660254037844386*fr[10]-1.5*fr[13]; 
  incr2[15] = 0.8660254037844386*fr[14]-1.5*fr[15]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 

  } else {

  incr2[1] = 1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[5] = 1.5*fl[5]+0.8660254037844386*fl[2]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844386*fl[3]; 
  incr2[8] = 1.5*fl[8]+0.8660254037844386*fl[4]; 
  incr2[11] = 1.5*fl[11]+0.8660254037844386*fl[7]; 
  incr2[12] = 1.5*fl[12]+0.8660254037844386*fl[9]; 
  incr2[13] = 1.5*fl[13]+0.8660254037844386*fl[10]; 
  incr2[15] = 1.5*fl[15]+0.8660254037844386*fl[14]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf4xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  if (idxr[1] == 1) {

  incr2[2] = 0.8660254037844386*fr[0]-1.5*fr[2]; 
  incr2[5] = 0.8660254037844386*fr[1]-1.5*fr[5]; 
  incr2[7] = 0.8660254037844386*fr[3]-1.5*fr[7]; 
  incr2[9] = 0.8660254037844386*fr[4]-1.5*fr[9]; 
  incr2[11] = 0.8660254037844386*fr[6]-1.5*fr[11]; 
  incr2[12] = 0.8660254037844386*fr[8]-1.5*fr[12]; 
  incr2[14] = 0.8660254037844386*fr[10]-1.5*fr[14]; 
  incr2[15] = 0.8660254037844386*fr[13]-1.5*fr[15]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 

  } else {

  incr2[2] = 1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[5] = 1.5*fl[5]+0.8660254037844386*fl[1]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844386*fl[3]; 
  incr2[9] = 1.5*fl[9]+0.8660254037844386*fl[4]; 
  incr2[11] = 1.5*fl[11]+0.8660254037844386*fl[6]; 
  incr2[12] = 1.5*fl[12]+0.8660254037844386*fl[8]; 
  incr2[14] = 1.5*fl[14]+0.8660254037844386*fl[10]; 
  incr2[15] = 1.5*fl[15]+0.8660254037844386*fl[13]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf4xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  if (idxr[2] == 1) {

  incr2[3] = 0.8660254037844386*fr[0]-1.5*fr[3]; 
  incr2[6] = 0.8660254037844386*fr[1]-1.5*fr[6]; 
  incr2[7] = 0.8660254037844386*fr[2]-1.5*fr[7]; 
  incr2[10] = 0.8660254037844386*fr[4]-1.5*fr[10]; 
  incr2[11] = 0.8660254037844386*fr[5]-1.5*fr[11]; 
  incr2[13] = 0.8660254037844386*fr[8]-1.5*fr[13]; 
  incr2[14] = 0.8660254037844386*fr[9]-1.5*fr[14]; 
  incr2[15] = 0.8660254037844386*fr[12]-1.5*fr[15]; 

  outr[3] += incr2[3]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 

  } else {

  incr2[3] = 1.5*fl[3]+0.8660254037844386*fl[0]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844386*fl[1]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844386*fl[2]; 
  incr2[10] = 1.5*fl[10]+0.8660254037844386*fl[4]; 
  incr2[11] = 1.5*fl[11]+0.8660254037844386*fl[5]; 
  incr2[13] = 1.5*fl[13]+0.8660254037844386*fl[8]; 
  incr2[14] = 1.5*fl[14]+0.8660254037844386*fl[9]; 
  incr2[15] = 1.5*fl[15]+0.8660254037844386*fl[12]; 

  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf4xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxFnur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  if (idxr[3] == 1) {

  incr2[4] = 0.8660254037844386*fr[0]-1.5*fr[4]; 
  incr2[8] = 0.8660254037844386*fr[1]-1.5*fr[8]; 
  incr2[9] = 0.8660254037844386*fr[2]-1.5*fr[9]; 
  incr2[10] = 0.8660254037844386*fr[3]-1.5*fr[10]; 
  incr2[12] = 0.8660254037844386*fr[5]-1.5*fr[12]; 
  incr2[13] = 0.8660254037844386*fr[6]-1.5*fr[13]; 
  incr2[14] = 0.8660254037844386*fr[7]-1.5*fr[14]; 
  incr2[15] = 0.8660254037844386*fr[11]-1.5*fr[15]; 

  outr[4] += incr2[4]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 

  } else {

  incr2[4] = 1.5*fl[4]+0.8660254037844386*fl[0]; 
  incr2[8] = 1.5*fl[8]+0.8660254037844386*fl[1]; 
  incr2[9] = 1.5*fl[9]+0.8660254037844386*fl[2]; 
  incr2[10] = 1.5*fl[10]+0.8660254037844386*fl[3]; 
  incr2[12] = 1.5*fl[12]+0.8660254037844386*fl[5]; 
  incr2[13] = 1.5*fl[13]+0.8660254037844386*fl[6]; 
  incr2[14] = 1.5*fl[14]+0.8660254037844386*fl[7]; 
  incr2[15] = 1.5*fl[15]+0.8660254037844386*fl[11]; 

  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 

  }

} 
