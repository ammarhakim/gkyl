#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf2xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[6]; 
  double incr2[6]; 

  if (idxr[0] == 1) {

  incr2[1] = 1.936491673103708*fr[4]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[3] = 0.8660254037844386*fr[2]-1.5*fr[3]; 
  incr2[4] = (-7.5*fr[4])+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 

  } else {

  incr2[1] = 1.936491673103708*fl[4]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[3] = 1.5*fl[3]+0.8660254037844386*fl[2]; 
  incr2[4] = (-7.5*fl[4])-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf2xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[6]; 
  double incr2[6]; 

  if (idxr[1] == 1) {

  incr2[2] = 1.936491673103708*fr[5]-1.5*fr[2]+0.8660254037844386*fr[0]; 
  incr2[3] = 0.8660254037844386*fr[1]-1.5*fr[3]; 
  incr2[5] = (-7.5*fr[5])+5.809475019311125*fr[2]-3.354101966249685*fr[0]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 

  } else {

  incr2[2] = 1.936491673103708*fl[5]+1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[3] = 1.5*fl[3]+0.8660254037844386*fl[1]; 
  incr2[5] = (-7.5*fl[5])-5.809475019311125*fl[2]-3.354101966249685*fl[0]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[6]; 
  double incr2[6]; 
  double incr3[6]; 
  double incr4[6]; 

  if (idxr[0] == 1) {

  incr2[1] = 5.809475019311125*fr[4]; 
  incr2[4] = -22.5*fr[4]; 

  incr3[4] = 22.5*fr[4]-5.809475019311125*fr[1]; 


  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[4] += (-1.0*incr3[4]*rdxFnur)-1.0*incr2[4]*rdxFnur; 

  } else {

  incr2[1] = 5.809475019311125*fl[4]; 
  incr2[4] = -22.5*fl[4]; 

  incr3[4] = (-22.5*fl[4])-5.809475019311125*fl[1]; 


  outl[1] += incr2[1]*rdxFnul; 
  outl[4] += incr3[4]*rdxFnul-1.0*incr2[4]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[6]; 
  double incr2[6]; 
  double incr3[6]; 
  double incr4[6]; 

  if (idxr[1] == 1) {

  incr2[2] = 5.809475019311125*fr[5]; 
  incr2[5] = -22.5*fr[5]; 

  incr3[5] = 22.5*fr[5]-5.809475019311125*fr[2]; 


  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[5] += (-1.0*incr3[5]*rdxFnur)-1.0*incr2[5]*rdxFnur; 

  } else {

  incr2[2] = 5.809475019311125*fl[5]; 
  incr2[5] = -22.5*fl[5]; 

  incr3[5] = (-22.5*fl[5])-5.809475019311125*fl[2]; 


  outl[2] += incr2[2]*rdxFnul; 
  outl[5] += incr3[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf2xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[6]; 
  double incr2[6]; 
  double incr3[6]; 
  double incr4[6]; 
  double incr5[6]; 
  double incr6[6]; 

  if (idxr[0] == 1) {







  } else {







  }

} 
void ConstHyperDiffusion6BoundarySurf2xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[6]; 
  double incr2[6]; 
  double incr3[6]; 
  double incr4[6]; 
  double incr5[6]; 
  double incr6[6]; 

  if (idxr[1] == 1) {







  } else {







  }

} 
