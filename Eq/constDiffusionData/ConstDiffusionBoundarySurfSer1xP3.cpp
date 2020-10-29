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


  incr2[1] = 11.78376607274358*fr[3]-11.61895003862225*fr[2]+4.5*fr[1]; 
  incr2[2] = (-45.6383297553399*fr[3])+45.0*fr[2]-17.42842505793337*fr[1]; 
  incr2[3] = 108.0*fr[3]-106.4894360957931*fr[2]+41.24318125460255*fr[1]; 


  incr4[3] = (-16.2*fr[3])+28.39718295887816*fr[2]-30.24499958670854*fr[1]+19.84313483298443*fr[0]; 

  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[3] += (-1.0*incr4[3]*rdxFnur)-1.0*incr2[3]*rdxFnur; 

  } else {


  incr2[1] = 11.78376607274358*fl[3]+11.61895003862225*fl[2]+4.5*fl[1]; 
  incr2[2] = 45.6383297553399*fl[3]+45.0*fl[2]+17.42842505793337*fl[1]; 
  incr2[3] = 108.0*fl[3]+106.4894360957931*fl[2]+41.24318125460255*fl[1]; 


  incr4[3] = (-16.2*fl[3])-28.39718295887816*fl[2]-30.24499958670854*fl[1]-19.84313483298443*fl[0]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += (-1.0*incr4[3]*rdxFnul)-1.0*incr2[3]*rdxFnul; 

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


  incr2[1] = (-103.1079531365064*fr[3])+76.2493596284585*fr[2]-19.6875*fr[1]; 
  incr2[2] = 399.3353853592242*fr[3]-295.3125*fr[2]+76.2493596284585*fr[1]; 
  incr2[3] = (-945.0*fr[3])+698.8369243786424*fr[2]-180.4389179888861*fr[1]; 


  incr4[3] = 225.0*fr[3]-241.2651286545313*fr[2]+96.66370606547471*fr[1]; 



  outr[1] += incr2[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr4[3]*rdxFnur+incr2[3]*rdxFnur; 

  } else {


  incr2[1] = (-103.1079531365064*fl[3])-76.2493596284585*fl[2]-19.6875*fl[1]; 
  incr2[2] = (-399.3353853592242*fl[3])-295.3125*fl[2]-76.2493596284585*fl[1]; 
  incr2[3] = (-945.0*fl[3])-698.8369243786424*fl[2]-180.4389179888861*fl[1]; 


  incr4[3] = 225.0*fl[3]+241.2651286545313*fl[2]+96.66370606547471*fl[1]; 



  outl[1] += incr2[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul; 
  outl[3] += incr4[3]*rdxFnul+incr2[3]*rdxFnul; 

  }

} 
