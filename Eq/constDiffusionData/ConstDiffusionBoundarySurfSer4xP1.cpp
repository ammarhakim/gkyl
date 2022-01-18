#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf4xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  if (edge < 0) {

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
void ConstDiffusionBoundarySurf4xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  if (edge < 0) {

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
void ConstDiffusionBoundarySurf4xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  if (edge < 0) {

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
void ConstDiffusionBoundarySurf4xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxFnur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  if (edge < 0) {

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
void ConstHyperDiffusion4BoundarySurf4xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 

  if (edge < 0) {


  incr2[1] = 1.5*fr[1]; 
  incr2[5] = 1.5*fr[5]; 
  incr2[6] = 1.5*fr[6]; 
  incr2[8] = 1.5*fr[8]; 
  incr2[11] = 1.5*fr[11]; 
  incr2[12] = 1.5*fr[12]; 
  incr2[13] = 1.5*fr[13]; 
  incr2[15] = 1.5*fr[15]; 



  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[8] += -1.0*incr2[8]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 

  } else {


  incr2[1] = 1.5*fl[1]; 
  incr2[5] = 1.5*fl[5]; 
  incr2[6] = 1.5*fl[6]; 
  incr2[8] = 1.5*fl[8]; 
  incr2[11] = 1.5*fl[11]; 
  incr2[12] = 1.5*fl[12]; 
  incr2[13] = 1.5*fl[13]; 
  incr2[15] = 1.5*fl[15]; 



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
void ConstHyperDiffusion4BoundarySurf4xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 

  if (edge < 0) {


  incr2[2] = 1.5*fr[2]; 
  incr2[5] = 1.5*fr[5]; 
  incr2[7] = 1.5*fr[7]; 
  incr2[9] = 1.5*fr[9]; 
  incr2[11] = 1.5*fr[11]; 
  incr2[12] = 1.5*fr[12]; 
  incr2[14] = 1.5*fr[14]; 
  incr2[15] = 1.5*fr[15]; 



  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[9] += -1.0*incr2[9]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 

  } else {


  incr2[2] = 1.5*fl[2]; 
  incr2[5] = 1.5*fl[5]; 
  incr2[7] = 1.5*fl[7]; 
  incr2[9] = 1.5*fl[9]; 
  incr2[11] = 1.5*fl[11]; 
  incr2[12] = 1.5*fl[12]; 
  incr2[14] = 1.5*fl[14]; 
  incr2[15] = 1.5*fl[15]; 



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
void ConstHyperDiffusion4BoundarySurf4xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 

  if (edge < 0) {


  incr2[3] = 1.5*fr[3]; 
  incr2[6] = 1.5*fr[6]; 
  incr2[7] = 1.5*fr[7]; 
  incr2[10] = 1.5*fr[10]; 
  incr2[11] = 1.5*fr[11]; 
  incr2[13] = 1.5*fr[13]; 
  incr2[14] = 1.5*fr[14]; 
  incr2[15] = 1.5*fr[15]; 



  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 

  } else {


  incr2[3] = 1.5*fl[3]; 
  incr2[6] = 1.5*fl[6]; 
  incr2[7] = 1.5*fl[7]; 
  incr2[10] = 1.5*fl[10]; 
  incr2[11] = 1.5*fl[11]; 
  incr2[13] = 1.5*fl[13]; 
  incr2[14] = 1.5*fl[14]; 
  incr2[15] = 1.5*fl[15]; 



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
void ConstHyperDiffusion4BoundarySurf4xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 16.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 

  if (edge < 0) {


  incr2[4] = 1.5*fr[4]; 
  incr2[8] = 1.5*fr[8]; 
  incr2[9] = 1.5*fr[9]; 
  incr2[10] = 1.5*fr[10]; 
  incr2[12] = 1.5*fr[12]; 
  incr2[13] = 1.5*fr[13]; 
  incr2[14] = 1.5*fr[14]; 
  incr2[15] = 1.5*fr[15]; 



  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[8] += -1.0*incr2[8]*rdxFnur; 
  outr[9] += -1.0*incr2[9]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 

  } else {


  incr2[4] = 1.5*fl[4]; 
  incr2[8] = 1.5*fl[8]; 
  incr2[9] = 1.5*fl[9]; 
  incr2[10] = 1.5*fl[10]; 
  incr2[12] = 1.5*fl[12]; 
  incr2[13] = 1.5*fl[13]; 
  incr2[14] = 1.5*fl[14]; 
  incr2[15] = 1.5*fl[15]; 



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
void ConstDiffusionVarCoeffBoundarySurf4xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[64]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  if (edge < 0) {

  incr2[1] = (-0.6495190528383289*fr[1]*nul[1])+0.375*fr[0]*nul[1]-0.375*nul[0]*fr[1]+0.2165063509461096*fr[0]*nul[0]; 
  incr2[5] = (-0.6495190528383289*nul[1]*fr[5])-0.375*nul[0]*fr[5]+0.375*nul[1]*fr[2]+0.2165063509461096*nul[0]*fr[2]; 
  incr2[6] = (-0.6495190528383289*nul[1]*fr[6])-0.375*nul[0]*fr[6]+0.375*nul[1]*fr[3]+0.2165063509461096*nul[0]*fr[3]; 
  incr2[8] = (-0.6495190528383289*nul[1]*fr[8])-0.375*nul[0]*fr[8]+0.375*nul[1]*fr[4]+0.2165063509461096*nul[0]*fr[4]; 
  incr2[11] = (-0.6495190528383289*nul[1]*fr[11])-0.375*nul[0]*fr[11]+0.375*nul[1]*fr[7]+0.2165063509461096*nul[0]*fr[7]; 
  incr2[12] = (-0.6495190528383289*nul[1]*fr[12])-0.375*nul[0]*fr[12]+0.375*nul[1]*fr[9]+0.2165063509461096*nul[0]*fr[9]; 
  incr2[13] = (-0.6495190528383289*nul[1]*fr[13])-0.375*nul[0]*fr[13]+0.375*nul[1]*fr[10]+0.2165063509461096*nul[0]*fr[10]; 
  incr2[15] = (-0.6495190528383289*nul[1]*fr[15])-0.375*nul[0]*fr[15]+0.375*nul[1]*fr[14]+0.2165063509461096*nul[0]*fr[14]; 

  outr[1] += incr2[1]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[8] += incr2[8]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 
  outr[12] += incr2[12]*rdxFr; 
  outr[13] += incr2[13]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 

  } else {

  incr2[1] = 0.6495190528383289*fl[1]*nul[1]+0.375*fl[0]*nul[1]+0.375*nul[0]*fl[1]+0.2165063509461096*fl[0]*nul[0]; 
  incr2[5] = 0.6495190528383289*nul[1]*fl[5]+0.375*nul[0]*fl[5]+0.375*nul[1]*fl[2]+0.2165063509461096*nul[0]*fl[2]; 
  incr2[6] = 0.6495190528383289*nul[1]*fl[6]+0.375*nul[0]*fl[6]+0.375*nul[1]*fl[3]+0.2165063509461096*nul[0]*fl[3]; 
  incr2[8] = 0.6495190528383289*nul[1]*fl[8]+0.375*nul[0]*fl[8]+0.375*nul[1]*fl[4]+0.2165063509461096*nul[0]*fl[4]; 
  incr2[11] = 0.6495190528383289*nul[1]*fl[11]+0.375*nul[0]*fl[11]+0.375*nul[1]*fl[7]+0.2165063509461096*nul[0]*fl[7]; 
  incr2[12] = 0.6495190528383289*nul[1]*fl[12]+0.375*nul[0]*fl[12]+0.375*nul[1]*fl[9]+0.2165063509461096*nul[0]*fl[9]; 
  incr2[13] = 0.6495190528383289*nul[1]*fl[13]+0.375*nul[0]*fl[13]+0.375*nul[1]*fl[10]+0.2165063509461096*nul[0]*fl[10]; 
  incr2[15] = 0.6495190528383289*nul[1]*fl[15]+0.375*nul[0]*fl[15]+0.375*nul[1]*fl[14]+0.2165063509461096*nul[0]*fl[14]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 
  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[8] += -1.0*incr2[8]*rdxFl; 
  outl[11] += -1.0*incr2[11]*rdxFl; 
  outl[12] += -1.0*incr2[12]*rdxFl; 
  outl[13] += -1.0*incr2[13]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf4xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[64]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[1]*dxl[1]); 
  double rdxFr = 4.0/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  if (edge < 0) {

  incr2[2] = (-0.6495190528383289*fr[2]*nul[18])+0.375*fr[0]*nul[18]-0.375*fr[2]*nul[16]+0.2165063509461096*fr[0]*nul[16]; 
  incr2[5] = (-0.6495190528383289*fr[5]*nul[18])+0.375*fr[1]*nul[18]-0.375*fr[5]*nul[16]+0.2165063509461096*fr[1]*nul[16]; 
  incr2[7] = (-0.6495190528383289*fr[7]*nul[18])+0.375*fr[3]*nul[18]-0.375*fr[7]*nul[16]+0.2165063509461096*fr[3]*nul[16]; 
  incr2[9] = (-0.6495190528383289*fr[9]*nul[18])+0.375*fr[4]*nul[18]-0.375*fr[9]*nul[16]+0.2165063509461096*fr[4]*nul[16]; 
  incr2[11] = (-0.6495190528383289*fr[11]*nul[18])+0.375*fr[6]*nul[18]-0.375*fr[11]*nul[16]+0.2165063509461096*fr[6]*nul[16]; 
  incr2[12] = (-0.6495190528383289*fr[12]*nul[18])+0.375*fr[8]*nul[18]-0.375*fr[12]*nul[16]+0.2165063509461096*fr[8]*nul[16]; 
  incr2[14] = (-0.6495190528383289*fr[14]*nul[18])+0.375*fr[10]*nul[18]-0.375*fr[14]*nul[16]+0.2165063509461096*fr[10]*nul[16]; 
  incr2[15] = (-0.6495190528383289*fr[15]*nul[18])+0.375*fr[13]*nul[18]-0.375*fr[15]*nul[16]+0.2165063509461096*fr[13]*nul[16]; 

  outr[2] += incr2[2]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 
  outr[9] += incr2[9]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 
  outr[12] += incr2[12]*rdxFr; 
  outr[14] += incr2[14]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 

  } else {

  incr2[2] = 0.6495190528383289*fl[2]*nul[18]+0.375*fl[0]*nul[18]+0.375*fl[2]*nul[16]+0.2165063509461096*fl[0]*nul[16]; 
  incr2[5] = 0.6495190528383289*fl[5]*nul[18]+0.375*fl[1]*nul[18]+0.375*fl[5]*nul[16]+0.2165063509461096*fl[1]*nul[16]; 
  incr2[7] = 0.6495190528383289*fl[7]*nul[18]+0.375*fl[3]*nul[18]+0.375*fl[7]*nul[16]+0.2165063509461096*fl[3]*nul[16]; 
  incr2[9] = 0.6495190528383289*fl[9]*nul[18]+0.375*fl[4]*nul[18]+0.375*fl[9]*nul[16]+0.2165063509461096*fl[4]*nul[16]; 
  incr2[11] = 0.6495190528383289*fl[11]*nul[18]+0.375*fl[6]*nul[18]+0.375*fl[11]*nul[16]+0.2165063509461096*fl[6]*nul[16]; 
  incr2[12] = 0.6495190528383289*fl[12]*nul[18]+0.375*fl[8]*nul[18]+0.375*fl[12]*nul[16]+0.2165063509461096*fl[8]*nul[16]; 
  incr2[14] = 0.6495190528383289*fl[14]*nul[18]+0.375*fl[10]*nul[18]+0.375*fl[14]*nul[16]+0.2165063509461096*fl[10]*nul[16]; 
  incr2[15] = 0.6495190528383289*fl[15]*nul[18]+0.375*fl[13]*nul[18]+0.375*fl[15]*nul[16]+0.2165063509461096*fl[13]*nul[16]; 

  outl[2] += -1.0*incr2[2]*rdxFl; 
  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 
  outl[9] += -1.0*incr2[9]*rdxFl; 
  outl[11] += -1.0*incr2[11]*rdxFl; 
  outl[12] += -1.0*incr2[12]*rdxFl; 
  outl[14] += -1.0*incr2[14]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf4xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[64]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[2]*dxl[2]); 
  double rdxFr = 4.0/(dxr[2]*dxr[2]); 

  double incr1[16]; 
  double incr2[16]; 

  if (edge < 0) {

  incr2[3] = (-0.6495190528383289*fr[3]*nul[35])+0.375*fr[0]*nul[35]-0.375*fr[3]*nul[32]+0.2165063509461096*fr[0]*nul[32]; 
  incr2[6] = (-0.6495190528383289*fr[6]*nul[35])+0.375*fr[1]*nul[35]-0.375*fr[6]*nul[32]+0.2165063509461096*fr[1]*nul[32]; 
  incr2[7] = (-0.6495190528383289*fr[7]*nul[35])+0.375*fr[2]*nul[35]-0.375*fr[7]*nul[32]+0.2165063509461096*fr[2]*nul[32]; 
  incr2[10] = (-0.6495190528383289*fr[10]*nul[35])+0.375*fr[4]*nul[35]-0.375*fr[10]*nul[32]+0.2165063509461096*fr[4]*nul[32]; 
  incr2[11] = (-0.6495190528383289*fr[11]*nul[35])+0.375*fr[5]*nul[35]-0.375*fr[11]*nul[32]+0.2165063509461096*fr[5]*nul[32]; 
  incr2[13] = (-0.6495190528383289*fr[13]*nul[35])+0.375*fr[8]*nul[35]-0.375*fr[13]*nul[32]+0.2165063509461096*fr[8]*nul[32]; 
  incr2[14] = (-0.6495190528383289*fr[14]*nul[35])+0.375*fr[9]*nul[35]-0.375*fr[14]*nul[32]+0.2165063509461096*fr[9]*nul[32]; 
  incr2[15] = (-0.6495190528383289*fr[15]*nul[35])+0.375*fr[12]*nul[35]-0.375*fr[15]*nul[32]+0.2165063509461096*fr[12]*nul[32]; 

  outr[3] += incr2[3]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 
  outr[13] += incr2[13]*rdxFr; 
  outr[14] += incr2[14]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 

  } else {

  incr2[3] = 0.6495190528383289*fl[3]*nul[35]+0.375*fl[0]*nul[35]+0.375*fl[3]*nul[32]+0.2165063509461096*fl[0]*nul[32]; 
  incr2[6] = 0.6495190528383289*fl[6]*nul[35]+0.375*fl[1]*nul[35]+0.375*fl[6]*nul[32]+0.2165063509461096*fl[1]*nul[32]; 
  incr2[7] = 0.6495190528383289*fl[7]*nul[35]+0.375*fl[2]*nul[35]+0.375*fl[7]*nul[32]+0.2165063509461096*fl[2]*nul[32]; 
  incr2[10] = 0.6495190528383289*fl[10]*nul[35]+0.375*fl[4]*nul[35]+0.375*fl[10]*nul[32]+0.2165063509461096*fl[4]*nul[32]; 
  incr2[11] = 0.6495190528383289*fl[11]*nul[35]+0.375*fl[5]*nul[35]+0.375*fl[11]*nul[32]+0.2165063509461096*fl[5]*nul[32]; 
  incr2[13] = 0.6495190528383289*fl[13]*nul[35]+0.375*fl[8]*nul[35]+0.375*fl[13]*nul[32]+0.2165063509461096*fl[8]*nul[32]; 
  incr2[14] = 0.6495190528383289*fl[14]*nul[35]+0.375*fl[9]*nul[35]+0.375*fl[14]*nul[32]+0.2165063509461096*fl[9]*nul[32]; 
  incr2[15] = 0.6495190528383289*fl[15]*nul[35]+0.375*fl[12]*nul[35]+0.375*fl[15]*nul[32]+0.2165063509461096*fl[12]*nul[32]; 

  outl[3] += -1.0*incr2[3]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[11] += -1.0*incr2[11]*rdxFl; 
  outl[13] += -1.0*incr2[13]*rdxFl; 
  outl[14] += -1.0*incr2[14]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf4xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[64]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[3]*dxl[3]); 
  double rdxFr = 4.0/(dxr[3]*dxr[3]); 

  double incr1[16]; 
  double incr2[16]; 

  if (edge < 0) {

  incr2[4] = (-0.6495190528383289*fr[4]*nul[52])+0.375*fr[0]*nul[52]-0.375*fr[4]*nul[48]+0.2165063509461096*fr[0]*nul[48]; 
  incr2[8] = (-0.6495190528383289*fr[8]*nul[52])+0.375*fr[1]*nul[52]-0.375*fr[8]*nul[48]+0.2165063509461096*fr[1]*nul[48]; 
  incr2[9] = (-0.6495190528383289*fr[9]*nul[52])+0.375*fr[2]*nul[52]-0.375*fr[9]*nul[48]+0.2165063509461096*fr[2]*nul[48]; 
  incr2[10] = (-0.6495190528383289*fr[10]*nul[52])+0.375*fr[3]*nul[52]-0.375*fr[10]*nul[48]+0.2165063509461096*fr[3]*nul[48]; 
  incr2[12] = (-0.6495190528383289*fr[12]*nul[52])+0.375*fr[5]*nul[52]-0.375*fr[12]*nul[48]+0.2165063509461096*fr[5]*nul[48]; 
  incr2[13] = (-0.6495190528383289*fr[13]*nul[52])+0.375*fr[6]*nul[52]-0.375*fr[13]*nul[48]+0.2165063509461096*fr[6]*nul[48]; 
  incr2[14] = (-0.6495190528383289*fr[14]*nul[52])+0.375*fr[7]*nul[52]-0.375*fr[14]*nul[48]+0.2165063509461096*fr[7]*nul[48]; 
  incr2[15] = (-0.6495190528383289*fr[15]*nul[52])+0.375*fr[11]*nul[52]-0.375*fr[15]*nul[48]+0.2165063509461096*fr[11]*nul[48]; 

  outr[4] += incr2[4]*rdxFr; 
  outr[8] += incr2[8]*rdxFr; 
  outr[9] += incr2[9]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[12] += incr2[12]*rdxFr; 
  outr[13] += incr2[13]*rdxFr; 
  outr[14] += incr2[14]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 

  } else {

  incr2[4] = 0.6495190528383289*fl[4]*nul[52]+0.375*fl[0]*nul[52]+0.375*fl[4]*nul[48]+0.2165063509461096*fl[0]*nul[48]; 
  incr2[8] = 0.6495190528383289*fl[8]*nul[52]+0.375*fl[1]*nul[52]+0.375*fl[8]*nul[48]+0.2165063509461096*fl[1]*nul[48]; 
  incr2[9] = 0.6495190528383289*fl[9]*nul[52]+0.375*fl[2]*nul[52]+0.375*fl[9]*nul[48]+0.2165063509461096*fl[2]*nul[48]; 
  incr2[10] = 0.6495190528383289*fl[10]*nul[52]+0.375*fl[3]*nul[52]+0.375*fl[10]*nul[48]+0.2165063509461096*fl[3]*nul[48]; 
  incr2[12] = 0.6495190528383289*fl[12]*nul[52]+0.375*fl[5]*nul[52]+0.375*fl[12]*nul[48]+0.2165063509461096*fl[5]*nul[48]; 
  incr2[13] = 0.6495190528383289*fl[13]*nul[52]+0.375*fl[6]*nul[52]+0.375*fl[13]*nul[48]+0.2165063509461096*fl[6]*nul[48]; 
  incr2[14] = 0.6495190528383289*fl[14]*nul[52]+0.375*fl[7]*nul[52]+0.375*fl[14]*nul[48]+0.2165063509461096*fl[7]*nul[48]; 
  incr2[15] = 0.6495190528383289*fl[15]*nul[52]+0.375*fl[11]*nul[52]+0.375*fl[15]*nul[48]+0.2165063509461096*fl[11]*nul[48]; 

  outl[4] += -1.0*incr2[4]*rdxFl; 
  outl[8] += -1.0*incr2[8]*rdxFl; 
  outl[9] += -1.0*incr2[9]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[12] += -1.0*incr2[12]*rdxFl; 
  outl[13] += -1.0*incr2[13]*rdxFl; 
  outl[14] += -1.0*incr2[14]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 

  }

} 
