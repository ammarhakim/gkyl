#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf3xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[10]; 
  double incr2[10]; 

  if (idxr[0] == 1) {

  incr2[1] = 1.936491673103708*fr[7]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[4] = 0.8660254037844386*fr[2]-1.5*fr[4]; 
  incr2[5] = 0.8660254037844386*fr[3]-1.5*fr[5]; 
  incr2[7] = (-7.5*fr[7])+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 

  } else {

  incr2[1] = 1.936491673103708*fl[7]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[4] = 1.5*fl[4]+0.8660254037844386*fl[2]; 
  incr2[5] = 1.5*fl[5]+0.8660254037844386*fl[3]; 
  incr2[7] = (-7.5*fl[7])-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf3xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[10]; 
  double incr2[10]; 

  if (idxr[1] == 1) {

  incr2[2] = 1.936491673103708*fr[8]-1.5*fr[2]+0.8660254037844386*fr[0]; 
  incr2[4] = 0.8660254037844386*fr[1]-1.5*fr[4]; 
  incr2[6] = 0.8660254037844386*fr[3]-1.5*fr[6]; 
  incr2[8] = (-7.5*fr[8])+5.809475019311125*fr[2]-3.354101966249685*fr[0]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 

  } else {

  incr2[2] = 1.936491673103708*fl[8]+1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[4] = 1.5*fl[4]+0.8660254037844386*fl[1]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844386*fl[3]; 
  incr2[8] = (-7.5*fl[8])-5.809475019311125*fl[2]-3.354101966249685*fl[0]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf3xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[10]; 
  double incr2[10]; 

  if (idxr[2] == 1) {

  incr2[3] = 1.936491673103708*fr[9]-1.5*fr[3]+0.8660254037844386*fr[0]; 
  incr2[5] = 0.8660254037844386*fr[1]-1.5*fr[5]; 
  incr2[6] = 0.8660254037844386*fr[2]-1.5*fr[6]; 
  incr2[9] = (-7.5*fr[9])+5.809475019311125*fr[3]-3.354101966249685*fr[0]; 

  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 

  } else {

  incr2[3] = 1.936491673103708*fl[9]+1.5*fl[3]+0.8660254037844386*fl[0]; 
  incr2[5] = 1.5*fl[5]+0.8660254037844386*fl[1]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844386*fl[2]; 
  incr2[9] = (-7.5*fl[9])-5.809475019311125*fl[3]-3.354101966249685*fl[0]; 

  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 

  if (idxr[0] == 1) {

  incr2[1] = 5.809475019311125*fr[7]; 
  incr2[7] = -22.5*fr[7]; 

  incr3[7] = 22.5*fr[7]-5.809475019311125*fr[1]; 


  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[7] += (-1.0*incr3[7]*rdxFnur)-1.0*incr2[7]*rdxFnur; 

  } else {

  incr2[1] = 5.809475019311125*fl[7]; 
  incr2[7] = -22.5*fl[7]; 

  incr3[7] = (-22.5*fl[7])-5.809475019311125*fl[1]; 


  outl[1] += incr2[1]*rdxFnul; 
  outl[7] += incr3[7]*rdxFnul-1.0*incr2[7]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 

  if (idxr[1] == 1) {

  incr2[2] = 5.809475019311125*fr[8]; 
  incr2[8] = -22.5*fr[8]; 

  incr3[8] = 22.5*fr[8]-5.809475019311125*fr[2]; 


  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[8] += (-1.0*incr3[8]*rdxFnur)-1.0*incr2[8]*rdxFnur; 

  } else {

  incr2[2] = 5.809475019311125*fl[8]; 
  incr2[8] = -22.5*fl[8]; 

  incr3[8] = (-22.5*fl[8])-5.809475019311125*fl[2]; 


  outl[2] += incr2[2]*rdxFnul; 
  outl[8] += incr3[8]*rdxFnul-1.0*incr2[8]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 

  if (idxr[2] == 1) {

  incr2[3] = 5.809475019311125*fr[9]; 
  incr2[9] = -22.5*fr[9]; 

  incr3[9] = 22.5*fr[9]-5.809475019311125*fr[3]; 


  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[9] += (-1.0*incr3[9]*rdxFnur)-1.0*incr2[9]*rdxFnur; 

  } else {

  incr2[3] = 5.809475019311125*fl[9]; 
  incr2[9] = -22.5*fl[9]; 

  incr3[9] = (-22.5*fl[9])-5.809475019311125*fl[3]; 


  outl[3] += incr2[3]*rdxFnul; 
  outl[9] += incr3[9]*rdxFnul-1.0*incr2[9]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf3xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 
  double incr5[10]; 
  double incr6[10]; 

  if (idxr[0] == 1) {







  } else {







  }

} 
void ConstHyperDiffusion6BoundarySurf3xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 
  double incr5[10]; 
  double incr6[10]; 

  if (idxr[1] == 1) {







  } else {







  }

} 
void ConstHyperDiffusion6BoundarySurf3xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 64.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[10]; 
  double incr2[10]; 
  double incr3[10]; 
  double incr4[10]; 
  double incr5[10]; 
  double incr6[10]; 

  if (idxr[2] == 1) {







  } else {







  }

} 
