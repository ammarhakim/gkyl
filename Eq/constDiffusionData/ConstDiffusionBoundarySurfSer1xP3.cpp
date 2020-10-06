#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf1xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[1]:    current grid index.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 

  if (idxr[0] == 1) {

  incr2[1] = (-2.29128784747792*fr[3])+1.936491673103708*fr[2]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[2] = 8.874119674649426*fr[3]-7.5*fr[2]+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 
  incr2[3] = (-21.0*fr[3])+17.74823934929885*fr[2]-13.74772708486752*fr[1]+7.937253933193772*fr[0]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 

  } else {

  incr2[1] = 2.29128784747792*fl[3]+1.936491673103708*fl[2]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[2] = (-8.874119674649426*fl[3])-7.5*fl[2]-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 
  incr2[3] = 21.0*fl[3]+17.74823934929885*fl[2]+13.74772708486752*fl[1]+7.937253933193772*fl[0]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf1xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[1]:    current grid index.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 
  double incr3[4]; 
  double incr4[4]; 

  if (idxr[0] == 1) {

  incr2[1] = 5.809475019311125*fr[2]-34.3693177121688*fr[3]; 
  incr2[2] = 133.1117951197414*fr[3]-22.5*fr[2]; 
  incr2[3] = 53.24471804789655*fr[2]-315.0*fr[3]; 

  incr3[2] = (-53.24471804789656*fr[3])+22.5*fr[2]-5.809475019311125*fr[1]; 
  incr3[3] = 315.0*fr[3]-133.1117951197414*fr[2]+34.3693177121688*fr[1]; 

  incr4[3] = (-52.5*fr[3])+44.37059837324713*fr[2]-34.3693177121688*fr[1]+19.84313483298443*fr[0]; 

  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[2] += (-1.0*incr3[2]*rdxFnur)-1.0*incr2[2]*rdxFnur; 
  outr[3] += (-1.0*incr4[3]*rdxFnur)-1.0*incr3[3]*rdxFnur-1.0*incr2[3]*rdxFnur; 

  } else {

  incr2[1] = 34.3693177121688*fl[3]+5.809475019311125*fl[2]; 
  incr2[2] = (-133.1117951197414*fl[3])-22.5*fl[2]; 
  incr2[3] = 315.0*fl[3]+53.24471804789655*fl[2]; 

  incr3[2] = (-53.24471804789656*fl[3])-22.5*fl[2]-5.809475019311125*fl[1]; 
  incr3[3] = 315.0*fl[3]+133.1117951197414*fl[2]+34.3693177121688*fl[1]; 

  incr4[3] = 52.5*fl[3]+44.37059837324713*fl[2]+34.3693177121688*fl[1]+19.84313483298443*fl[0]; 

  outl[1] += incr2[1]*rdxFnul; 
  outl[2] += incr3[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
  outl[3] += incr4[3]*rdxFnul-1.0*incr3[3]*rdxFnul+incr2[3]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf1xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // idx[1]:    current grid index.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 
  double incr3[4]; 
  double incr4[4]; 
  double incr5[4]; 
  double incr6[4]; 

  if (idxr[0] == 1) {


  incr3[2] = -133.1117951197414*fr[3]; 
  incr3[3] = 787.5*fr[3]; 

  incr4[3] = 133.1117951197414*fr[2]-787.5*fr[3]; 



  outr[2] += incr3[2]*rdxFnur; 
  outr[3] += incr4[3]*rdxFnur+incr3[3]*rdxFnur; 

  } else {


  incr3[2] = -133.1117951197414*fl[3]; 
  incr3[3] = 787.5*fl[3]; 

  incr4[3] = 787.5*fl[3]+133.1117951197414*fl[2]; 



  outl[2] += -1.0*incr3[2]*rdxFnul; 
  outl[3] += incr3[3]*rdxFnul-1.0*incr4[3]*rdxFnul; 

  }

} 
