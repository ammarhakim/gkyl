#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf4xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[15]; 
  double incr2[15]; 

  if (idxr[0] == 1) {

  incr2[1] = 1.936491673103708*fr[11]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[5] = 0.8660254037844386*fr[2]-1.5*fr[5]; 
  incr2[6] = 0.8660254037844386*fr[3]-1.5*fr[6]; 
  incr2[8] = 0.8660254037844386*fr[4]-1.5*fr[8]; 
  incr2[11] = (-7.5*fr[11])+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 

  } else {

  incr2[1] = 1.936491673103708*fl[11]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[5] = 1.5*fl[5]+0.8660254037844386*fl[2]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844386*fl[3]; 
  incr2[8] = 1.5*fl[8]+0.8660254037844386*fl[4]; 
  incr2[11] = (-7.5*fl[11])-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf4xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[15]; 
  double incr2[15]; 

  if (idxr[1] == 1) {

  incr2[2] = 1.936491673103708*fr[12]-1.5*fr[2]+0.8660254037844386*fr[0]; 
  incr2[5] = 0.8660254037844386*fr[1]-1.5*fr[5]; 
  incr2[7] = 0.8660254037844386*fr[3]-1.5*fr[7]; 
  incr2[9] = 0.8660254037844386*fr[4]-1.5*fr[9]; 
  incr2[12] = (-7.5*fr[12])+5.809475019311125*fr[2]-3.354101966249685*fr[0]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 

  } else {

  incr2[2] = 1.936491673103708*fl[12]+1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[5] = 1.5*fl[5]+0.8660254037844386*fl[1]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844386*fl[3]; 
  incr2[9] = 1.5*fl[9]+0.8660254037844386*fl[4]; 
  incr2[12] = (-7.5*fl[12])-5.809475019311125*fl[2]-3.354101966249685*fl[0]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf4xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[15]; 
  double incr2[15]; 

  if (idxr[2] == 1) {

  incr2[3] = 1.936491673103708*fr[13]-1.5*fr[3]+0.8660254037844386*fr[0]; 
  incr2[6] = 0.8660254037844386*fr[1]-1.5*fr[6]; 
  incr2[7] = 0.8660254037844386*fr[2]-1.5*fr[7]; 
  incr2[10] = 0.8660254037844386*fr[4]-1.5*fr[10]; 
  incr2[13] = (-7.5*fr[13])+5.809475019311125*fr[3]-3.354101966249685*fr[0]; 

  outr[3] += incr2[3]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 

  } else {

  incr2[3] = 1.936491673103708*fl[13]+1.5*fl[3]+0.8660254037844386*fl[0]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844386*fl[1]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844386*fl[2]; 
  incr2[10] = 1.5*fl[10]+0.8660254037844386*fl[4]; 
  incr2[13] = (-7.5*fl[13])-5.809475019311125*fl[3]-3.354101966249685*fl[0]; 

  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf4xMaxP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxFnur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

  double incr1[15]; 
  double incr2[15]; 

  if (idxr[3] == 1) {

  incr2[4] = 1.936491673103708*fr[14]-1.5*fr[4]+0.8660254037844386*fr[0]; 
  incr2[8] = 0.8660254037844386*fr[1]-1.5*fr[8]; 
  incr2[9] = 0.8660254037844386*fr[2]-1.5*fr[9]; 
  incr2[10] = 0.8660254037844386*fr[3]-1.5*fr[10]; 
  incr2[14] = (-7.5*fr[14])+5.809475019311125*fr[4]-3.354101966249685*fr[0]; 

  outr[4] += incr2[4]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 

  } else {

  incr2[4] = 1.936491673103708*fl[14]+1.5*fl[4]+0.8660254037844386*fl[0]; 
  incr2[8] = 1.5*fl[8]+0.8660254037844386*fl[1]; 
  incr2[9] = 1.5*fl[9]+0.8660254037844386*fl[2]; 
  incr2[10] = 1.5*fl[10]+0.8660254037844386*fl[3]; 
  incr2[14] = (-7.5*fl[14])-5.809475019311125*fl[4]-3.354101966249685*fl[0]; 

  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf4xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 

  if (idxr[0] == 1) {

  incr2[1] = 5.809475019311125*fr[11]; 
  incr2[11] = -22.5*fr[11]; 

  incr3[11] = 22.5*fr[11]-5.809475019311125*fr[1]; 


  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[11] += (-1.0*incr3[11]*rdxFnur)-1.0*incr2[11]*rdxFnur; 

  } else {

  incr2[1] = 5.809475019311125*fl[11]; 
  incr2[11] = -22.5*fl[11]; 

  incr3[11] = (-22.5*fl[11])-5.809475019311125*fl[1]; 


  outl[1] += incr2[1]*rdxFnul; 
  outl[11] += incr3[11]*rdxFnul-1.0*incr2[11]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf4xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 

  if (idxr[1] == 1) {

  incr2[2] = 5.809475019311125*fr[12]; 
  incr2[12] = -22.5*fr[12]; 

  incr3[12] = 22.5*fr[12]-5.809475019311125*fr[2]; 


  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[12] += (-1.0*incr3[12]*rdxFnur)-1.0*incr2[12]*rdxFnur; 

  } else {

  incr2[2] = 5.809475019311125*fl[12]; 
  incr2[12] = -22.5*fl[12]; 

  incr3[12] = (-22.5*fl[12])-5.809475019311125*fl[2]; 


  outl[2] += incr2[2]*rdxFnul; 
  outl[12] += incr3[12]*rdxFnul-1.0*incr2[12]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf4xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 

  if (idxr[2] == 1) {

  incr2[3] = 5.809475019311125*fr[13]; 
  incr2[13] = -22.5*fr[13]; 

  incr3[13] = 22.5*fr[13]-5.809475019311125*fr[3]; 


  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[13] += (-1.0*incr3[13]*rdxFnur)-1.0*incr2[13]*rdxFnur; 

  } else {

  incr2[3] = 5.809475019311125*fl[13]; 
  incr2[13] = -22.5*fl[13]; 

  incr3[13] = (-22.5*fl[13])-5.809475019311125*fl[3]; 


  outl[3] += incr2[3]*rdxFnul; 
  outl[13] += incr3[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf4xMaxP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 16.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 

  if (idxr[3] == 1) {

  incr2[4] = 5.809475019311125*fr[14]; 
  incr2[14] = -22.5*fr[14]; 

  incr3[14] = 22.5*fr[14]-5.809475019311125*fr[4]; 


  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[14] += (-1.0*incr3[14]*rdxFnur)-1.0*incr2[14]*rdxFnur; 

  } else {

  incr2[4] = 5.809475019311125*fl[14]; 
  incr2[14] = -22.5*fl[14]; 

  incr3[14] = (-22.5*fl[14])-5.809475019311125*fl[4]; 


  outl[4] += incr2[4]*rdxFnul; 
  outl[14] += incr3[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf4xMaxP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 
  double incr5[15]; 
  double incr6[15]; 

  if (idxr[0] == 1) {







  } else {







  }

} 
void ConstHyperDiffusion6BoundarySurf4xMaxP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 
  double incr5[15]; 
  double incr6[15]; 

  if (idxr[1] == 1) {







  } else {







  }

} 
void ConstHyperDiffusion6BoundarySurf4xMaxP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 64.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 
  double incr5[15]; 
  double incr6[15]; 

  if (idxr[2] == 1) {







  } else {







  }

} 
void ConstHyperDiffusion6BoundarySurf4xMaxP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 64.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[15]; 
  double incr2[15]; 
  double incr3[15]; 
  double incr4[15]; 
  double incr5[15]; 
  double incr6[15]; 

  if (idxr[3] == 1) {







  } else {







  }

} 
