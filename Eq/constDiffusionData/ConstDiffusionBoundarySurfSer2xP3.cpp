#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf2xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[12]; 
  double incr2[12]; 

  if (idxr[0] == 1) {

  incr2[1] = (-2.29128784747792*fr[8])+1.936491673103708*fr[4]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[3] = (-2.29128784747792*fr[10])+1.936491673103709*fr[6]-1.5*fr[3]+0.8660254037844386*fr[2]; 
  incr2[4] = 8.874119674649426*fr[8]-7.5*fr[4]+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 
  incr2[6] = 8.874119674649425*fr[10]-7.5*fr[6]+5.809475019311126*fr[3]-3.354101966249684*fr[2]; 
  incr2[7] = 0.8660254037844387*fr[5]-1.5*fr[7]; 
  incr2[8] = (-21.0*fr[8])+17.74823934929885*fr[4]-13.74772708486752*fr[1]+7.937253933193772*fr[0]; 
  incr2[10] = (-21.0*fr[10])+17.74823934929885*fr[6]-13.74772708486752*fr[3]+7.937253933193771*fr[2]; 
  incr2[11] = 0.8660254037844386*fr[9]-1.5*fr[11]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 

  } else {

  incr2[1] = 2.29128784747792*fl[8]+1.936491673103708*fl[4]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[3] = 2.29128784747792*fl[10]+1.936491673103709*fl[6]+1.5*fl[3]+0.8660254037844386*fl[2]; 
  incr2[4] = (-8.874119674649426*fl[8])-7.5*fl[4]-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 
  incr2[6] = (-8.874119674649425*fl[10])-7.5*fl[6]-5.809475019311126*fl[3]-3.354101966249684*fl[2]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844387*fl[5]; 
  incr2[8] = 21.0*fl[8]+17.74823934929885*fl[4]+13.74772708486752*fl[1]+7.937253933193772*fl[0]; 
  incr2[10] = 21.0*fl[10]+17.74823934929885*fl[6]+13.74772708486752*fl[3]+7.937253933193771*fl[2]; 
  incr2[11] = 1.5*fl[11]+0.8660254037844386*fl[9]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf2xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[12]; 
  double incr2[12]; 

  if (idxr[1] == 1) {

  incr2[2] = (-2.29128784747792*fr[9])+1.936491673103708*fr[5]-1.5*fr[2]+0.8660254037844386*fr[0]; 
  incr2[3] = (-2.29128784747792*fr[11])+1.936491673103709*fr[7]-1.5*fr[3]+0.8660254037844386*fr[1]; 
  incr2[5] = 8.874119674649426*fr[9]-7.5*fr[5]+5.809475019311125*fr[2]-3.354101966249685*fr[0]; 
  incr2[6] = 0.8660254037844387*fr[4]-1.5*fr[6]; 
  incr2[7] = 8.874119674649425*fr[11]-7.5*fr[7]+5.809475019311126*fr[3]-3.354101966249684*fr[1]; 
  incr2[9] = (-21.0*fr[9])+17.74823934929885*fr[5]-13.74772708486752*fr[2]+7.937253933193772*fr[0]; 
  incr2[10] = 0.8660254037844386*fr[8]-1.5*fr[10]; 
  incr2[11] = (-21.0*fr[11])+17.74823934929885*fr[7]-13.74772708486752*fr[3]+7.937253933193771*fr[1]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 

  } else {

  incr2[2] = 2.29128784747792*fl[9]+1.936491673103708*fl[5]+1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[3] = 2.29128784747792*fl[11]+1.936491673103709*fl[7]+1.5*fl[3]+0.8660254037844386*fl[1]; 
  incr2[5] = (-8.874119674649426*fl[9])-7.5*fl[5]-5.809475019311125*fl[2]-3.354101966249685*fl[0]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844387*fl[4]; 
  incr2[7] = (-8.874119674649425*fl[11])-7.5*fl[7]-5.809475019311126*fl[3]-3.354101966249684*fl[1]; 
  incr2[9] = 21.0*fl[9]+17.74823934929885*fl[5]+13.74772708486752*fl[2]+7.937253933193772*fl[0]; 
  incr2[10] = 1.5*fl[10]+0.8660254037844386*fl[8]; 
  incr2[11] = 21.0*fl[11]+17.74823934929885*fl[7]+13.74772708486752*fl[3]+7.937253933193771*fl[1]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[12]; 
  double incr2[12]; 
  double incr3[12]; 
  double incr4[12]; 

  if (idxr[0] == 1) {

  incr2[1] = 5.809475019311125*fr[4]-34.3693177121688*fr[8]; 
  incr2[3] = 5.809475019311126*fr[6]-34.3693177121688*fr[10]; 
  incr2[4] = 133.1117951197414*fr[8]-22.5*fr[4]; 
  incr2[6] = 133.1117951197414*fr[10]-22.5*fr[6]; 
  incr2[8] = 53.24471804789655*fr[4]-315.0*fr[8]; 
  incr2[10] = 53.24471804789655*fr[6]-315.0*fr[10]; 

  incr3[4] = (-53.24471804789656*fr[8])+22.5*fr[4]-5.809475019311125*fr[1]; 
  incr3[6] = (-53.24471804789655*fr[10])+22.5*fr[6]-5.809475019311126*fr[3]; 
  incr3[8] = 315.0*fr[8]-133.1117951197414*fr[4]+34.3693177121688*fr[1]; 
  incr3[10] = 315.0*fr[10]-133.1117951197414*fr[6]+34.3693177121688*fr[3]; 

  incr4[8] = (-52.5*fr[8])+44.37059837324713*fr[4]-34.3693177121688*fr[1]+19.84313483298443*fr[0]; 
  incr4[10] = (-52.5*fr[10])+44.37059837324713*fr[6]-34.3693177121688*fr[3]+19.84313483298443*fr[2]; 

  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[4] += (-1.0*incr3[4]*rdxFnur)-1.0*incr2[4]*rdxFnur; 
  outr[6] += (-1.0*incr3[6]*rdxFnur)-1.0*incr2[6]*rdxFnur; 
  outr[8] += (-1.0*incr4[8]*rdxFnur)-1.0*incr3[8]*rdxFnur-1.0*incr2[8]*rdxFnur; 
  outr[10] += (-1.0*incr4[10]*rdxFnur)-1.0*incr3[10]*rdxFnur-1.0*incr2[10]*rdxFnur; 

  } else {

  incr2[1] = 34.3693177121688*fl[8]+5.809475019311125*fl[4]; 
  incr2[3] = 34.3693177121688*fl[10]+5.809475019311126*fl[6]; 
  incr2[4] = (-133.1117951197414*fl[8])-22.5*fl[4]; 
  incr2[6] = (-133.1117951197414*fl[10])-22.5*fl[6]; 
  incr2[8] = 315.0*fl[8]+53.24471804789655*fl[4]; 
  incr2[10] = 315.0*fl[10]+53.24471804789655*fl[6]; 

  incr3[4] = (-53.24471804789656*fl[8])-22.5*fl[4]-5.809475019311125*fl[1]; 
  incr3[6] = (-53.24471804789655*fl[10])-22.5*fl[6]-5.809475019311126*fl[3]; 
  incr3[8] = 315.0*fl[8]+133.1117951197414*fl[4]+34.3693177121688*fl[1]; 
  incr3[10] = 315.0*fl[10]+133.1117951197414*fl[6]+34.3693177121688*fl[3]; 

  incr4[8] = 52.5*fl[8]+44.37059837324713*fl[4]+34.3693177121688*fl[1]+19.84313483298443*fl[0]; 
  incr4[10] = 52.5*fl[10]+44.37059837324713*fl[6]+34.3693177121688*fl[3]+19.84313483298443*fl[2]; 

  outl[1] += incr2[1]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul; 
  outl[4] += incr3[4]*rdxFnul-1.0*incr2[4]*rdxFnul; 
  outl[6] += incr3[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[8] += incr4[8]*rdxFnul-1.0*incr3[8]*rdxFnul+incr2[8]*rdxFnul; 
  outl[10] += incr4[10]*rdxFnul-1.0*incr3[10]*rdxFnul+incr2[10]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[12]; 
  double incr2[12]; 
  double incr3[12]; 
  double incr4[12]; 

  if (idxr[1] == 1) {

  incr2[2] = 5.809475019311125*fr[5]-34.3693177121688*fr[9]; 
  incr2[3] = 5.809475019311126*fr[7]-34.3693177121688*fr[11]; 
  incr2[5] = 133.1117951197414*fr[9]-22.5*fr[5]; 
  incr2[7] = 133.1117951197414*fr[11]-22.5*fr[7]; 
  incr2[9] = 53.24471804789655*fr[5]-315.0*fr[9]; 
  incr2[11] = 53.24471804789655*fr[7]-315.0*fr[11]; 

  incr3[5] = (-53.24471804789656*fr[9])+22.5*fr[5]-5.809475019311125*fr[2]; 
  incr3[7] = (-53.24471804789655*fr[11])+22.5*fr[7]-5.809475019311126*fr[3]; 
  incr3[9] = 315.0*fr[9]-133.1117951197414*fr[5]+34.3693177121688*fr[2]; 
  incr3[11] = 315.0*fr[11]-133.1117951197414*fr[7]+34.3693177121688*fr[3]; 

  incr4[9] = (-52.5*fr[9])+44.37059837324713*fr[5]-34.3693177121688*fr[2]+19.84313483298443*fr[0]; 
  incr4[11] = (-52.5*fr[11])+44.37059837324713*fr[7]-34.3693177121688*fr[3]+19.84313483298443*fr[1]; 

  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[5] += (-1.0*incr3[5]*rdxFnur)-1.0*incr2[5]*rdxFnur; 
  outr[7] += (-1.0*incr3[7]*rdxFnur)-1.0*incr2[7]*rdxFnur; 
  outr[9] += (-1.0*incr4[9]*rdxFnur)-1.0*incr3[9]*rdxFnur-1.0*incr2[9]*rdxFnur; 
  outr[11] += (-1.0*incr4[11]*rdxFnur)-1.0*incr3[11]*rdxFnur-1.0*incr2[11]*rdxFnur; 

  } else {

  incr2[2] = 34.3693177121688*fl[9]+5.809475019311125*fl[5]; 
  incr2[3] = 34.3693177121688*fl[11]+5.809475019311126*fl[7]; 
  incr2[5] = (-133.1117951197414*fl[9])-22.5*fl[5]; 
  incr2[7] = (-133.1117951197414*fl[11])-22.5*fl[7]; 
  incr2[9] = 315.0*fl[9]+53.24471804789655*fl[5]; 
  incr2[11] = 315.0*fl[11]+53.24471804789655*fl[7]; 

  incr3[5] = (-53.24471804789656*fl[9])-22.5*fl[5]-5.809475019311125*fl[2]; 
  incr3[7] = (-53.24471804789655*fl[11])-22.5*fl[7]-5.809475019311126*fl[3]; 
  incr3[9] = 315.0*fl[9]+133.1117951197414*fl[5]+34.3693177121688*fl[2]; 
  incr3[11] = 315.0*fl[11]+133.1117951197414*fl[7]+34.3693177121688*fl[3]; 

  incr4[9] = 52.5*fl[9]+44.37059837324713*fl[5]+34.3693177121688*fl[2]+19.84313483298443*fl[0]; 
  incr4[11] = 52.5*fl[11]+44.37059837324713*fl[7]+34.3693177121688*fl[3]+19.84313483298443*fl[1]; 

  outl[2] += incr2[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul; 
  outl[5] += incr3[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[7] += incr3[7]*rdxFnul-1.0*incr2[7]*rdxFnul; 
  outl[9] += incr4[9]*rdxFnul-1.0*incr3[9]*rdxFnul+incr2[9]*rdxFnul; 
  outl[11] += incr4[11]*rdxFnul-1.0*incr3[11]*rdxFnul+incr2[11]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf2xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[12]; 
  double incr2[12]; 
  double incr3[12]; 
  double incr4[12]; 
  double incr5[12]; 
  double incr6[12]; 

  if (idxr[0] == 1) {


  incr3[4] = -133.1117951197414*fr[8]; 
  incr3[6] = -133.1117951197414*fr[10]; 
  incr3[8] = 787.5*fr[8]; 
  incr3[10] = 787.5*fr[10]; 

  incr4[8] = 133.1117951197414*fr[4]-787.5*fr[8]; 
  incr4[10] = 133.1117951197414*fr[6]-787.5*fr[10]; 



  outr[4] += incr3[4]*rdxFnur; 
  outr[6] += incr3[6]*rdxFnur; 
  outr[8] += incr4[8]*rdxFnur+incr3[8]*rdxFnur; 
  outr[10] += incr4[10]*rdxFnur+incr3[10]*rdxFnur; 

  } else {


  incr3[4] = -133.1117951197414*fl[8]; 
  incr3[6] = -133.1117951197414*fl[10]; 
  incr3[8] = 787.5*fl[8]; 
  incr3[10] = 787.5*fl[10]; 

  incr4[8] = 787.5*fl[8]+133.1117951197414*fl[4]; 
  incr4[10] = 787.5*fl[10]+133.1117951197414*fl[6]; 



  outl[4] += -1.0*incr3[4]*rdxFnul; 
  outl[6] += -1.0*incr3[6]*rdxFnul; 
  outl[8] += incr3[8]*rdxFnul-1.0*incr4[8]*rdxFnul; 
  outl[10] += incr3[10]*rdxFnul-1.0*incr4[10]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf2xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[12]; 
  double incr2[12]; 
  double incr3[12]; 
  double incr4[12]; 
  double incr5[12]; 
  double incr6[12]; 

  if (idxr[1] == 1) {


  incr3[5] = -133.1117951197414*fr[9]; 
  incr3[7] = -133.1117951197414*fr[11]; 
  incr3[9] = 787.5*fr[9]; 
  incr3[11] = 787.5*fr[11]; 

  incr4[9] = 133.1117951197414*fr[5]-787.5*fr[9]; 
  incr4[11] = 133.1117951197414*fr[7]-787.5*fr[11]; 



  outr[5] += incr3[5]*rdxFnur; 
  outr[7] += incr3[7]*rdxFnur; 
  outr[9] += incr4[9]*rdxFnur+incr3[9]*rdxFnur; 
  outr[11] += incr4[11]*rdxFnur+incr3[11]*rdxFnur; 

  } else {


  incr3[5] = -133.1117951197414*fl[9]; 
  incr3[7] = -133.1117951197414*fl[11]; 
  incr3[9] = 787.5*fl[9]; 
  incr3[11] = 787.5*fl[11]; 

  incr4[9] = 787.5*fl[9]+133.1117951197414*fl[5]; 
  incr4[11] = 787.5*fl[11]+133.1117951197414*fl[7]; 



  outl[5] += -1.0*incr3[5]*rdxFnul; 
  outl[7] += -1.0*incr3[7]*rdxFnul; 
  outl[9] += incr3[9]*rdxFnul-1.0*incr4[9]*rdxFnul; 
  outl[11] += incr3[11]*rdxFnul-1.0*incr4[11]*rdxFnul; 

  }

} 
