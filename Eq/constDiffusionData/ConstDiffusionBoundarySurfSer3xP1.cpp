#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf3xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 

  if (edge < 0) {

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
void ConstDiffusionBoundarySurf3xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 

  if (edge < 0) {

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
void ConstDiffusionBoundarySurf3xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[8]; 
  double incr2[8]; 

  if (edge < 0) {

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
void ConstHyperDiffusion4BoundarySurf3xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 
  double incr3[8]; 
  double incr4[8]; 

  if (edge < 0) {


  incr2[1] = 1.5*fr[1]; 
  incr2[4] = 1.5*fr[4]; 
  incr2[5] = 1.5*fr[5]; 
  incr2[7] = 1.5*fr[7]; 



  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 

  } else {


  incr2[1] = 1.5*fl[1]; 
  incr2[4] = 1.5*fl[4]; 
  incr2[5] = 1.5*fl[5]; 
  incr2[7] = 1.5*fl[7]; 



  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 
  double incr3[8]; 
  double incr4[8]; 

  if (edge < 0) {


  incr2[2] = 1.5*fr[2]; 
  incr2[4] = 1.5*fr[4]; 
  incr2[6] = 1.5*fr[6]; 
  incr2[7] = 1.5*fr[7]; 



  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 

  } else {


  incr2[2] = 1.5*fl[2]; 
  incr2[4] = 1.5*fl[4]; 
  incr2[6] = 1.5*fl[6]; 
  incr2[7] = 1.5*fl[7]; 



  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[8]; 
  double incr2[8]; 
  double incr3[8]; 
  double incr4[8]; 

  if (edge < 0) {


  incr2[3] = 1.5*fr[3]; 
  incr2[5] = 1.5*fr[5]; 
  incr2[6] = 1.5*fr[6]; 
  incr2[7] = 1.5*fr[7]; 



  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 

  } else {


  incr2[3] = 1.5*fl[3]; 
  incr2[5] = 1.5*fl[5]; 
  incr2[6] = 1.5*fl[6]; 
  incr2[7] = 1.5*fl[7]; 



  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf3xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[24]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[8]; 
  double incr2[8]; 

  if (edge < 0) {

  incr2[1] = (-0.9185586535436913*fr[1]*nul[1])+0.5303300858899105*fr[0]*nul[1]-0.5303300858899105*nul[0]*fr[1]+0.3061862178478971*fr[0]*nul[0]; 
  incr2[4] = (-0.9185586535436913*nul[1]*fr[4])-0.5303300858899105*nul[0]*fr[4]+0.5303300858899105*nul[1]*fr[2]+0.3061862178478971*nul[0]*fr[2]; 
  incr2[5] = (-0.9185586535436913*nul[1]*fr[5])-0.5303300858899105*nul[0]*fr[5]+0.5303300858899105*nul[1]*fr[3]+0.3061862178478971*nul[0]*fr[3]; 
  incr2[7] = (-0.9185586535436913*nul[1]*fr[7])-0.5303300858899105*nul[0]*fr[7]+0.5303300858899105*nul[1]*fr[6]+0.3061862178478971*nul[0]*fr[6]; 

  outr[1] += incr2[1]*rdxFr; 
  outr[4] += incr2[4]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 

  } else {

  incr2[1] = 0.9185586535436913*fl[1]*nul[1]+0.5303300858899105*fl[0]*nul[1]+0.5303300858899105*nul[0]*fl[1]+0.3061862178478971*fl[0]*nul[0]; 
  incr2[4] = 0.9185586535436913*nul[1]*fl[4]+0.5303300858899105*nul[0]*fl[4]+0.5303300858899105*nul[1]*fl[2]+0.3061862178478971*nul[0]*fl[2]; 
  incr2[5] = 0.9185586535436913*nul[1]*fl[5]+0.5303300858899105*nul[0]*fl[5]+0.5303300858899105*nul[1]*fl[3]+0.3061862178478971*nul[0]*fl[3]; 
  incr2[7] = 0.9185586535436913*nul[1]*fl[7]+0.5303300858899105*nul[0]*fl[7]+0.5303300858899105*nul[1]*fl[6]+0.3061862178478971*nul[0]*fl[6]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 
  outl[4] += -1.0*incr2[4]*rdxFl; 
  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf3xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[24]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[1]*dxl[1]); 
  double rdxFr = 4.0/(dxr[1]*dxr[1]); 

  double incr1[8]; 
  double incr2[8]; 

  if (edge < 0) {

  incr2[2] = (-0.9185586535436913*fr[2]*nul[10])+0.5303300858899105*fr[0]*nul[10]-0.5303300858899105*fr[2]*nul[8]+0.3061862178478971*fr[0]*nul[8]; 
  incr2[4] = (-0.9185586535436913*fr[4]*nul[10])+0.5303300858899105*fr[1]*nul[10]-0.5303300858899105*fr[4]*nul[8]+0.3061862178478971*fr[1]*nul[8]; 
  incr2[6] = (-0.9185586535436913*fr[6]*nul[10])+0.5303300858899105*fr[3]*nul[10]-0.5303300858899105*fr[6]*nul[8]+0.3061862178478971*fr[3]*nul[8]; 
  incr2[7] = (-0.9185586535436913*fr[7]*nul[10])+0.5303300858899105*fr[5]*nul[10]-0.5303300858899105*fr[7]*nul[8]+0.3061862178478971*fr[5]*nul[8]; 

  outr[2] += incr2[2]*rdxFr; 
  outr[4] += incr2[4]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 

  } else {

  incr2[2] = 0.9185586535436913*fl[2]*nul[10]+0.5303300858899105*fl[0]*nul[10]+0.5303300858899105*fl[2]*nul[8]+0.3061862178478971*fl[0]*nul[8]; 
  incr2[4] = 0.9185586535436913*fl[4]*nul[10]+0.5303300858899105*fl[1]*nul[10]+0.5303300858899105*fl[4]*nul[8]+0.3061862178478971*fl[1]*nul[8]; 
  incr2[6] = 0.9185586535436913*fl[6]*nul[10]+0.5303300858899105*fl[3]*nul[10]+0.5303300858899105*fl[6]*nul[8]+0.3061862178478971*fl[3]*nul[8]; 
  incr2[7] = 0.9185586535436913*fl[7]*nul[10]+0.5303300858899105*fl[5]*nul[10]+0.5303300858899105*fl[7]*nul[8]+0.3061862178478971*fl[5]*nul[8]; 

  outl[2] += -1.0*incr2[2]*rdxFl; 
  outl[4] += -1.0*incr2[4]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf3xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[24]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[2]*dxl[2]); 
  double rdxFr = 4.0/(dxr[2]*dxr[2]); 

  double incr1[8]; 
  double incr2[8]; 

  if (edge < 0) {

  incr2[3] = (-0.9185586535436913*fr[3]*nul[19])+0.5303300858899105*fr[0]*nul[19]-0.5303300858899105*fr[3]*nul[16]+0.3061862178478971*fr[0]*nul[16]; 
  incr2[5] = (-0.9185586535436913*fr[5]*nul[19])+0.5303300858899105*fr[1]*nul[19]-0.5303300858899105*fr[5]*nul[16]+0.3061862178478971*fr[1]*nul[16]; 
  incr2[6] = (-0.9185586535436913*fr[6]*nul[19])+0.5303300858899105*fr[2]*nul[19]-0.5303300858899105*fr[6]*nul[16]+0.3061862178478971*fr[2]*nul[16]; 
  incr2[7] = (-0.9185586535436913*fr[7]*nul[19])+0.5303300858899105*fr[4]*nul[19]-0.5303300858899105*fr[7]*nul[16]+0.3061862178478971*fr[4]*nul[16]; 

  outr[3] += incr2[3]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 

  } else {

  incr2[3] = 0.9185586535436913*fl[3]*nul[19]+0.5303300858899105*fl[0]*nul[19]+0.5303300858899105*fl[3]*nul[16]+0.3061862178478971*fl[0]*nul[16]; 
  incr2[5] = 0.9185586535436913*fl[5]*nul[19]+0.5303300858899105*fl[1]*nul[19]+0.5303300858899105*fl[5]*nul[16]+0.3061862178478971*fl[1]*nul[16]; 
  incr2[6] = 0.9185586535436913*fl[6]*nul[19]+0.5303300858899105*fl[2]*nul[19]+0.5303300858899105*fl[6]*nul[16]+0.3061862178478971*fl[2]*nul[16]; 
  incr2[7] = 0.9185586535436913*fl[7]*nul[19]+0.5303300858899105*fl[4]*nul[19]+0.5303300858899105*fl[7]*nul[16]+0.3061862178478971*fl[4]*nul[16]; 

  outl[3] += -1.0*incr2[3]*rdxFl; 
  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 

  }

} 
