#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf2xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 

  if (idxr[0] == 1) {

  incr2[1] = 0.8660254037844386*fr[0]-1.5*fr[1]; 
  incr2[3] = 0.8660254037844386*fr[2]-1.5*fr[3]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 

  } else {

  incr2[1] = 1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[3] = 1.5*fl[3]+0.8660254037844386*fl[2]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf2xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[4]; 
  double incr2[4]; 

  if (idxr[1] == 1) {

  incr2[2] = 0.8660254037844386*fr[0]-1.5*fr[2]; 
  incr2[3] = 0.8660254037844386*fr[1]-1.5*fr[3]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 

  } else {

  incr2[2] = 1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[3] = 1.5*fl[3]+0.8660254037844386*fl[1]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 
  double incr3[4]; 
  double incr4[4]; 

  if (idxr[0] == 1) {


  incr2[1] = 1.5*fr[1]; 
  incr2[3] = 1.5*fr[3]; 



  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 

  } else {


  incr2[1] = 1.5*fl[1]; 
  incr2[3] = 1.5*fl[3]; 



  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[4]; 
  double incr2[4]; 
  double incr3[4]; 
  double incr4[4]; 

  if (idxr[1] == 1) {


  incr2[2] = 1.5*fr[2]; 
  incr2[3] = 1.5*fr[3]; 



  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 

  } else {


  incr2[2] = 1.5*fl[2]; 
  incr2[3] = 1.5*fl[3]; 



  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf2xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[8]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 

  if (idxr[0] == 1) {

  incr2[1] = (-1.299038105676658*fr[3]*nul[3])+0.75*fr[2]*nul[3]-0.75*nul[2]*fr[3]+0.4330127018922193*fr[2]*nul[2]-1.299038105676658*fr[1]*nul[1]+0.75*fr[0]*nul[1]-0.75*nul[0]*fr[1]+0.4330127018922193*fr[0]*nul[0]; 
  incr2[3] = (-1.299038105676658*fr[1]*nul[3])+0.75*fr[0]*nul[3]-1.299038105676658*nul[1]*fr[3]-0.75*nul[0]*fr[3]-0.75*fr[1]*nul[2]+0.4330127018922193*fr[0]*nul[2]+0.75*nul[1]*fr[2]+0.4330127018922193*nul[0]*fr[2]; 

  outr[1] += incr2[1]*rdxFr; 
  outr[3] += incr2[3]*rdxFr; 

  } else {

  incr2[1] = 1.299038105676658*fl[3]*nul[3]+0.75*fl[2]*nul[3]+0.75*nul[2]*fl[3]+0.4330127018922193*fl[2]*nul[2]+1.299038105676658*fl[1]*nul[1]+0.75*fl[0]*nul[1]+0.75*nul[0]*fl[1]+0.4330127018922193*fl[0]*nul[0]; 
  incr2[3] = 1.299038105676658*fl[1]*nul[3]+0.75*fl[0]*nul[3]+1.299038105676658*nul[1]*fl[3]+0.75*nul[0]*fl[3]+0.75*fl[1]*nul[2]+0.4330127018922193*fl[0]*nul[2]+0.75*nul[1]*fl[2]+0.4330127018922193*nul[0]*fl[2]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 
  outl[3] += -1.0*incr2[3]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf2xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[8]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[1]*dxl[1]); 
  double rdxFr = 4.0/(dxr[1]*dxr[1]); 

  double incr1[4]; 
  double incr2[4]; 

  if (idxr[1] == 1) {

  incr2[2] = (-1.299038105676658*fr[3]*nul[7])+0.75*fr[1]*nul[7]-1.299038105676658*fr[2]*nul[6]+0.75*fr[0]*nul[6]-0.75*fr[3]*nul[5]+0.4330127018922193*fr[1]*nul[5]-0.75*fr[2]*nul[4]+0.4330127018922193*fr[0]*nul[4]; 
  incr2[3] = (-1.299038105676658*fr[2]*nul[7])+0.75*fr[0]*nul[7]-1.299038105676658*fr[3]*nul[6]+0.75*fr[1]*nul[6]-0.75*fr[2]*nul[5]+0.4330127018922193*fr[0]*nul[5]-0.75*fr[3]*nul[4]+0.4330127018922193*fr[1]*nul[4]; 

  outr[2] += incr2[2]*rdxFr; 
  outr[3] += incr2[3]*rdxFr; 

  } else {

  incr2[2] = 1.299038105676658*fl[3]*nul[7]+0.75*fl[1]*nul[7]+1.299038105676658*fl[2]*nul[6]+0.75*fl[0]*nul[6]+0.75*fl[3]*nul[5]+0.4330127018922193*fl[1]*nul[5]+0.75*fl[2]*nul[4]+0.4330127018922193*fl[0]*nul[4]; 
  incr2[3] = 1.299038105676658*fl[2]*nul[7]+0.75*fl[0]*nul[7]+1.299038105676658*fl[3]*nul[6]+0.75*fl[1]*nul[6]+0.75*fl[2]*nul[5]+0.4330127018922193*fl[0]*nul[5]+0.75*fl[3]*nul[4]+0.4330127018922193*fl[1]*nul[4]; 

  outl[2] += -1.0*incr2[2]*rdxFl; 
  outl[3] += -1.0*incr2[3]*rdxFl; 

  }

} 
