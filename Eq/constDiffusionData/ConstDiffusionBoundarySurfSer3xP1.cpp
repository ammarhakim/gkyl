#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf3xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 

  if (idxr[0] == 1) {

  incr2[1] = 0.8660254037844386*fr[0]-1.5*fr[1]; 
  incr2[4] = 0.8660254037844386*fr[2]-1.5*fr[4]; 
  incr2[5] = 0.8660254037844386*fr[3]-1.5*fr[5]; 
  incr2[7] = 0.8660254037844386*fr[6]-1.5*fr[7]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 

  } else {

  incr2[1] = 1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[4] = 1.5*fl[4]+0.8660254037844386*fl[2]; 
  incr2[5] = 1.5*fl[5]+0.8660254037844386*fl[3]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844386*fl[6]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf3xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 

  if (idxr[1] == 1) {

  incr2[2] = 0.8660254037844386*fr[0]-1.5*fr[2]; 
  incr2[4] = 0.8660254037844386*fr[1]-1.5*fr[4]; 
  incr2[6] = 0.8660254037844386*fr[3]-1.5*fr[6]; 
  incr2[7] = 0.8660254037844386*fr[5]-1.5*fr[7]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 

  } else {

  incr2[2] = 1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[4] = 1.5*fl[4]+0.8660254037844386*fl[1]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844386*fl[3]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844386*fl[5]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf3xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[8]; 
  double incr2[8]; 

  if (idxr[2] == 1) {

  incr2[3] = 0.8660254037844386*fr[0]-1.5*fr[3]; 
  incr2[5] = 0.8660254037844386*fr[1]-1.5*fr[5]; 
  incr2[6] = 0.8660254037844386*fr[2]-1.5*fr[6]; 
  incr2[7] = 0.8660254037844386*fr[4]-1.5*fr[7]; 

  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 

  } else {

  incr2[3] = 1.5*fl[3]+0.8660254037844386*fl[0]; 
  incr2[5] = 1.5*fl[5]+0.8660254037844386*fl[1]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844386*fl[2]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844386*fl[4]; 

  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 

  }

} 
