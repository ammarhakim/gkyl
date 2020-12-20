#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf2xMaxP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[3]; 
  double incr2[3]; 

  if (idxr[0] == 1) {

  incr2[1] = 0.8660254037844386*fr[0]-1.5*fr[1]; 

  outr[1] += incr2[1]*rdxFnur; 

  } else {

  incr2[1] = 1.5*fl[1]+0.8660254037844386*fl[0]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf2xMaxP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[3]; 
  double incr2[3]; 

  if (idxr[1] == 1) {

  incr2[2] = 0.8660254037844386*fr[0]-1.5*fr[2]; 

  outr[2] += incr2[2]*rdxFnur; 

  } else {

  incr2[2] = 1.5*fl[2]+0.8660254037844386*fl[0]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 

  }

} 
