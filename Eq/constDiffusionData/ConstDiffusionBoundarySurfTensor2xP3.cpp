#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf2xTensorP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 

  if (idxr[0] == 1) {

  incr2[1] = (-2.29128784747792*fr[8])+1.936491673103708*fr[4]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[3] = (-2.29128784747792*fr[11])+1.936491673103709*fr[6]-1.5*fr[3]+0.8660254037844386*fr[2]; 
  incr2[4] = 8.874119674649426*fr[8]-7.5*fr[4]+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 
  incr2[6] = 8.874119674649425*fr[11]-7.5*fr[6]+5.809475019311126*fr[3]-3.354101966249684*fr[2]; 
  incr2[7] = (-2.29128784747792*fr[13])+1.936491673103709*fr[10]-1.5*fr[7]+0.8660254037844387*fr[5]; 
  incr2[8] = (-21.0*fr[8])+17.74823934929885*fr[4]-13.74772708486752*fr[1]+7.937253933193772*fr[0]; 
  incr2[10] = 8.874119674649425*fr[13]-7.5*fr[10]+5.809475019311126*fr[7]-3.354101966249685*fr[5]; 
  incr2[11] = (-21.0*fr[11])+17.74823934929885*fr[6]-13.74772708486752*fr[3]+7.937253933193771*fr[2]; 
  incr2[12] = (-2.29128784747792*fr[15])+1.936491673103708*fr[14]-1.5*fr[12]+0.8660254037844386*fr[9]; 
  incr2[13] = (-21.0*fr[13])+17.74823934929885*fr[10]-13.74772708486752*fr[7]+7.937253933193772*fr[5]; 
  incr2[14] = 8.874119674649423*fr[15]-7.5*fr[14]+5.809475019311125*fr[12]-3.354101966249684*fr[9]; 
  incr2[15] = (-21.0*fr[15])+17.74823934929885*fr[14]-13.74772708486752*fr[12]+7.937253933193772*fr[9]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 

  } else {

  incr2[1] = 2.29128784747792*fl[8]+1.936491673103708*fl[4]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[3] = 2.29128784747792*fl[11]+1.936491673103709*fl[6]+1.5*fl[3]+0.8660254037844386*fl[2]; 
  incr2[4] = (-8.874119674649426*fl[8])-7.5*fl[4]-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 
  incr2[6] = (-8.874119674649425*fl[11])-7.5*fl[6]-5.809475019311126*fl[3]-3.354101966249684*fl[2]; 
  incr2[7] = 2.29128784747792*fl[13]+1.936491673103709*fl[10]+1.5*fl[7]+0.8660254037844387*fl[5]; 
  incr2[8] = 21.0*fl[8]+17.74823934929885*fl[4]+13.74772708486752*fl[1]+7.937253933193772*fl[0]; 
  incr2[10] = (-8.874119674649425*fl[13])-7.5*fl[10]-5.809475019311126*fl[7]-3.354101966249685*fl[5]; 
  incr2[11] = 21.0*fl[11]+17.74823934929885*fl[6]+13.74772708486752*fl[3]+7.937253933193771*fl[2]; 
  incr2[12] = 2.29128784747792*fl[15]+1.936491673103708*fl[14]+1.5*fl[12]+0.8660254037844386*fl[9]; 
  incr2[13] = 21.0*fl[13]+17.74823934929885*fl[10]+13.74772708486752*fl[7]+7.937253933193772*fl[5]; 
  incr2[14] = (-8.874119674649423*fl[15])-7.5*fl[14]-5.809475019311125*fl[12]-3.354101966249684*fl[9]; 
  incr2[15] = 21.0*fl[15]+17.74823934929885*fl[14]+13.74772708486752*fl[12]+7.937253933193772*fl[9]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf2xTensorP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 

  if (idxr[1] == 1) {

  incr2[2] = (-2.29128784747792*fr[9])+1.936491673103708*fr[5]-1.5*fr[2]+0.8660254037844386*fr[0]; 
  incr2[3] = (-2.29128784747792*fr[12])+1.936491673103709*fr[7]-1.5*fr[3]+0.8660254037844386*fr[1]; 
  incr2[5] = 8.874119674649426*fr[9]-7.5*fr[5]+5.809475019311125*fr[2]-3.354101966249685*fr[0]; 
  incr2[6] = (-2.29128784747792*fr[14])+1.936491673103709*fr[10]-1.5*fr[6]+0.8660254037844387*fr[4]; 
  incr2[7] = 8.874119674649425*fr[12]-7.5*fr[7]+5.809475019311126*fr[3]-3.354101966249684*fr[1]; 
  incr2[9] = (-21.0*fr[9])+17.74823934929885*fr[5]-13.74772708486752*fr[2]+7.937253933193772*fr[0]; 
  incr2[10] = 8.874119674649425*fr[14]-7.5*fr[10]+5.809475019311126*fr[6]-3.354101966249685*fr[4]; 
  incr2[11] = (-2.29128784747792*fr[15])+1.936491673103708*fr[13]-1.5*fr[11]+0.8660254037844386*fr[8]; 
  incr2[12] = (-21.0*fr[12])+17.74823934929885*fr[7]-13.74772708486752*fr[3]+7.937253933193771*fr[1]; 
  incr2[13] = 8.874119674649423*fr[15]-7.5*fr[13]+5.809475019311125*fr[11]-3.354101966249684*fr[8]; 
  incr2[14] = (-21.0*fr[14])+17.74823934929885*fr[10]-13.74772708486752*fr[6]+7.937253933193772*fr[4]; 
  incr2[15] = (-21.0*fr[15])+17.74823934929885*fr[13]-13.74772708486752*fr[11]+7.937253933193772*fr[8]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 

  } else {

  incr2[2] = 2.29128784747792*fl[9]+1.936491673103708*fl[5]+1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[3] = 2.29128784747792*fl[12]+1.936491673103709*fl[7]+1.5*fl[3]+0.8660254037844386*fl[1]; 
  incr2[5] = (-8.874119674649426*fl[9])-7.5*fl[5]-5.809475019311125*fl[2]-3.354101966249685*fl[0]; 
  incr2[6] = 2.29128784747792*fl[14]+1.936491673103709*fl[10]+1.5*fl[6]+0.8660254037844387*fl[4]; 
  incr2[7] = (-8.874119674649425*fl[12])-7.5*fl[7]-5.809475019311126*fl[3]-3.354101966249684*fl[1]; 
  incr2[9] = 21.0*fl[9]+17.74823934929885*fl[5]+13.74772708486752*fl[2]+7.937253933193772*fl[0]; 
  incr2[10] = (-8.874119674649425*fl[14])-7.5*fl[10]-5.809475019311126*fl[6]-3.354101966249685*fl[4]; 
  incr2[11] = 2.29128784747792*fl[15]+1.936491673103708*fl[13]+1.5*fl[11]+0.8660254037844386*fl[8]; 
  incr2[12] = 21.0*fl[12]+17.74823934929885*fl[7]+13.74772708486752*fl[3]+7.937253933193771*fl[1]; 
  incr2[13] = (-8.874119674649423*fl[15])-7.5*fl[13]-5.809475019311125*fl[11]-3.354101966249684*fl[8]; 
  incr2[14] = 21.0*fl[14]+17.74823934929885*fl[10]+13.74772708486752*fl[6]+7.937253933193772*fl[4]; 
  incr2[15] = 21.0*fl[15]+17.74823934929885*fl[13]+13.74772708486752*fl[11]+7.937253933193772*fl[8]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xTensorP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 

  if (idxr[0] == 1) {

  incr2[1] = 5.809475019311125*fr[4]-34.3693177121688*fr[8]; 
  incr2[3] = 5.809475019311126*fr[6]-34.3693177121688*fr[11]; 
  incr2[4] = 133.1117951197414*fr[8]-22.5*fr[4]; 
  incr2[6] = 133.1117951197414*fr[11]-22.5*fr[6]; 
  incr2[7] = 5.809475019311126*fr[10]-34.3693177121688*fr[13]; 
  incr2[8] = 53.24471804789655*fr[4]-315.0*fr[8]; 
  incr2[10] = 133.1117951197414*fr[13]-22.5*fr[10]; 
  incr2[11] = 53.24471804789655*fr[6]-315.0*fr[11]; 
  incr2[12] = 5.809475019311125*fr[14]-34.36931771216879*fr[15]; 
  incr2[13] = 53.24471804789655*fr[10]-315.0*fr[13]; 
  incr2[14] = 133.1117951197414*fr[15]-22.5*fr[14]; 
  incr2[15] = 53.24471804789655*fr[14]-315.0*fr[15]; 

  incr3[4] = (-53.24471804789656*fr[8])+22.5*fr[4]-5.809475019311125*fr[1]; 
  incr3[6] = (-53.24471804789655*fr[11])+22.5*fr[6]-5.809475019311126*fr[3]; 
  incr3[8] = 315.0*fr[8]-133.1117951197414*fr[4]+34.3693177121688*fr[1]; 
  incr3[10] = (-53.24471804789655*fr[13])+22.5*fr[10]-5.809475019311126*fr[7]; 
  incr3[11] = 315.0*fr[11]-133.1117951197414*fr[6]+34.3693177121688*fr[3]; 
  incr3[13] = 315.0*fr[13]-133.1117951197414*fr[10]+34.3693177121688*fr[7]; 
  incr3[14] = (-53.24471804789655*fr[15])+22.5*fr[14]-5.809475019311125*fr[12]; 
  incr3[15] = 315.0*fr[15]-133.1117951197414*fr[14]+34.3693177121688*fr[12]; 

  incr4[8] = (-52.5*fr[8])+44.37059837324713*fr[4]-34.3693177121688*fr[1]+19.84313483298443*fr[0]; 
  incr4[11] = (-52.5*fr[11])+44.37059837324713*fr[6]-34.3693177121688*fr[3]+19.84313483298443*fr[2]; 
  incr4[13] = (-52.5*fr[13])+44.37059837324712*fr[10]-34.3693177121688*fr[7]+19.84313483298443*fr[5]; 
  incr4[15] = (-52.5*fr[15])+44.37059837324712*fr[14]-34.3693177121688*fr[12]+19.84313483298443*fr[9]; 

  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[4] += (-1.0*incr3[4]*rdxFnur)-1.0*incr2[4]*rdxFnur; 
  outr[6] += (-1.0*incr3[6]*rdxFnur)-1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[8] += (-1.0*incr4[8]*rdxFnur)-1.0*incr3[8]*rdxFnur-1.0*incr2[8]*rdxFnur; 
  outr[10] += (-1.0*incr3[10]*rdxFnur)-1.0*incr2[10]*rdxFnur; 
  outr[11] += (-1.0*incr4[11]*rdxFnur)-1.0*incr3[11]*rdxFnur-1.0*incr2[11]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[13] += (-1.0*incr4[13]*rdxFnur)-1.0*incr3[13]*rdxFnur-1.0*incr2[13]*rdxFnur; 
  outr[14] += (-1.0*incr3[14]*rdxFnur)-1.0*incr2[14]*rdxFnur; 
  outr[15] += (-1.0*incr4[15]*rdxFnur)-1.0*incr3[15]*rdxFnur-1.0*incr2[15]*rdxFnur; 

  } else {

  incr2[1] = 34.3693177121688*fl[8]+5.809475019311125*fl[4]; 
  incr2[3] = 34.3693177121688*fl[11]+5.809475019311126*fl[6]; 
  incr2[4] = (-133.1117951197414*fl[8])-22.5*fl[4]; 
  incr2[6] = (-133.1117951197414*fl[11])-22.5*fl[6]; 
  incr2[7] = 34.3693177121688*fl[13]+5.809475019311126*fl[10]; 
  incr2[8] = 315.0*fl[8]+53.24471804789655*fl[4]; 
  incr2[10] = (-133.1117951197414*fl[13])-22.5*fl[10]; 
  incr2[11] = 315.0*fl[11]+53.24471804789655*fl[6]; 
  incr2[12] = 34.36931771216879*fl[15]+5.809475019311125*fl[14]; 
  incr2[13] = 315.0*fl[13]+53.24471804789655*fl[10]; 
  incr2[14] = (-133.1117951197414*fl[15])-22.5*fl[14]; 
  incr2[15] = 315.0*fl[15]+53.24471804789655*fl[14]; 

  incr3[4] = (-53.24471804789656*fl[8])-22.5*fl[4]-5.809475019311125*fl[1]; 
  incr3[6] = (-53.24471804789655*fl[11])-22.5*fl[6]-5.809475019311126*fl[3]; 
  incr3[8] = 315.0*fl[8]+133.1117951197414*fl[4]+34.3693177121688*fl[1]; 
  incr3[10] = (-53.24471804789655*fl[13])-22.5*fl[10]-5.809475019311126*fl[7]; 
  incr3[11] = 315.0*fl[11]+133.1117951197414*fl[6]+34.3693177121688*fl[3]; 
  incr3[13] = 315.0*fl[13]+133.1117951197414*fl[10]+34.3693177121688*fl[7]; 
  incr3[14] = (-53.24471804789655*fl[15])-22.5*fl[14]-5.809475019311125*fl[12]; 
  incr3[15] = 315.0*fl[15]+133.1117951197414*fl[14]+34.3693177121688*fl[12]; 

  incr4[8] = 52.5*fl[8]+44.37059837324713*fl[4]+34.3693177121688*fl[1]+19.84313483298443*fl[0]; 
  incr4[11] = 52.5*fl[11]+44.37059837324713*fl[6]+34.3693177121688*fl[3]+19.84313483298443*fl[2]; 
  incr4[13] = 52.5*fl[13]+44.37059837324712*fl[10]+34.3693177121688*fl[7]+19.84313483298443*fl[5]; 
  incr4[15] = 52.5*fl[15]+44.37059837324712*fl[14]+34.3693177121688*fl[12]+19.84313483298443*fl[9]; 

  outl[1] += incr2[1]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul; 
  outl[4] += incr3[4]*rdxFnul-1.0*incr2[4]*rdxFnul; 
  outl[6] += incr3[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[8] += incr4[8]*rdxFnul-1.0*incr3[8]*rdxFnul+incr2[8]*rdxFnul; 
  outl[10] += incr3[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += incr4[11]*rdxFnul-1.0*incr3[11]*rdxFnul+incr2[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 
  outl[13] += incr4[13]*rdxFnul-1.0*incr3[13]*rdxFnul+incr2[13]*rdxFnul; 
  outl[14] += incr3[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 
  outl[15] += incr4[15]*rdxFnul-1.0*incr3[15]*rdxFnul+incr2[15]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xTensorP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 

  if (idxr[1] == 1) {

  incr2[2] = 5.809475019311125*fr[5]-34.3693177121688*fr[9]; 
  incr2[3] = 5.809475019311126*fr[7]-34.3693177121688*fr[12]; 
  incr2[5] = 133.1117951197414*fr[9]-22.5*fr[5]; 
  incr2[6] = 5.809475019311126*fr[10]-34.3693177121688*fr[14]; 
  incr2[7] = 133.1117951197414*fr[12]-22.5*fr[7]; 
  incr2[9] = 53.24471804789655*fr[5]-315.0*fr[9]; 
  incr2[10] = 133.1117951197414*fr[14]-22.5*fr[10]; 
  incr2[11] = 5.809475019311125*fr[13]-34.36931771216879*fr[15]; 
  incr2[12] = 53.24471804789655*fr[7]-315.0*fr[12]; 
  incr2[13] = 133.1117951197414*fr[15]-22.5*fr[13]; 
  incr2[14] = 53.24471804789655*fr[10]-315.0*fr[14]; 
  incr2[15] = 53.24471804789655*fr[13]-315.0*fr[15]; 

  incr3[5] = (-53.24471804789656*fr[9])+22.5*fr[5]-5.809475019311125*fr[2]; 
  incr3[7] = (-53.24471804789655*fr[12])+22.5*fr[7]-5.809475019311126*fr[3]; 
  incr3[9] = 315.0*fr[9]-133.1117951197414*fr[5]+34.3693177121688*fr[2]; 
  incr3[10] = (-53.24471804789655*fr[14])+22.5*fr[10]-5.809475019311126*fr[6]; 
  incr3[12] = 315.0*fr[12]-133.1117951197414*fr[7]+34.3693177121688*fr[3]; 
  incr3[13] = (-53.24471804789655*fr[15])+22.5*fr[13]-5.809475019311125*fr[11]; 
  incr3[14] = 315.0*fr[14]-133.1117951197414*fr[10]+34.3693177121688*fr[6]; 
  incr3[15] = 315.0*fr[15]-133.1117951197414*fr[13]+34.3693177121688*fr[11]; 

  incr4[9] = (-52.5*fr[9])+44.37059837324713*fr[5]-34.3693177121688*fr[2]+19.84313483298443*fr[0]; 
  incr4[12] = (-52.5*fr[12])+44.37059837324713*fr[7]-34.3693177121688*fr[3]+19.84313483298443*fr[1]; 
  incr4[14] = (-52.5*fr[14])+44.37059837324712*fr[10]-34.3693177121688*fr[6]+19.84313483298443*fr[4]; 
  incr4[15] = (-52.5*fr[15])+44.37059837324712*fr[13]-34.3693177121688*fr[11]+19.84313483298443*fr[8]; 

  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[5] += (-1.0*incr3[5]*rdxFnur)-1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += (-1.0*incr3[7]*rdxFnur)-1.0*incr2[7]*rdxFnur; 
  outr[9] += (-1.0*incr4[9]*rdxFnur)-1.0*incr3[9]*rdxFnur-1.0*incr2[9]*rdxFnur; 
  outr[10] += (-1.0*incr3[10]*rdxFnur)-1.0*incr2[10]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[12] += (-1.0*incr4[12]*rdxFnur)-1.0*incr3[12]*rdxFnur-1.0*incr2[12]*rdxFnur; 
  outr[13] += (-1.0*incr3[13]*rdxFnur)-1.0*incr2[13]*rdxFnur; 
  outr[14] += (-1.0*incr4[14]*rdxFnur)-1.0*incr3[14]*rdxFnur-1.0*incr2[14]*rdxFnur; 
  outr[15] += (-1.0*incr4[15]*rdxFnur)-1.0*incr3[15]*rdxFnur-1.0*incr2[15]*rdxFnur; 

  } else {

  incr2[2] = 34.3693177121688*fl[9]+5.809475019311125*fl[5]; 
  incr2[3] = 34.3693177121688*fl[12]+5.809475019311126*fl[7]; 
  incr2[5] = (-133.1117951197414*fl[9])-22.5*fl[5]; 
  incr2[6] = 34.3693177121688*fl[14]+5.809475019311126*fl[10]; 
  incr2[7] = (-133.1117951197414*fl[12])-22.5*fl[7]; 
  incr2[9] = 315.0*fl[9]+53.24471804789655*fl[5]; 
  incr2[10] = (-133.1117951197414*fl[14])-22.5*fl[10]; 
  incr2[11] = 34.36931771216879*fl[15]+5.809475019311125*fl[13]; 
  incr2[12] = 315.0*fl[12]+53.24471804789655*fl[7]; 
  incr2[13] = (-133.1117951197414*fl[15])-22.5*fl[13]; 
  incr2[14] = 315.0*fl[14]+53.24471804789655*fl[10]; 
  incr2[15] = 315.0*fl[15]+53.24471804789655*fl[13]; 

  incr3[5] = (-53.24471804789656*fl[9])-22.5*fl[5]-5.809475019311125*fl[2]; 
  incr3[7] = (-53.24471804789655*fl[12])-22.5*fl[7]-5.809475019311126*fl[3]; 
  incr3[9] = 315.0*fl[9]+133.1117951197414*fl[5]+34.3693177121688*fl[2]; 
  incr3[10] = (-53.24471804789655*fl[14])-22.5*fl[10]-5.809475019311126*fl[6]; 
  incr3[12] = 315.0*fl[12]+133.1117951197414*fl[7]+34.3693177121688*fl[3]; 
  incr3[13] = (-53.24471804789655*fl[15])-22.5*fl[13]-5.809475019311125*fl[11]; 
  incr3[14] = 315.0*fl[14]+133.1117951197414*fl[10]+34.3693177121688*fl[6]; 
  incr3[15] = 315.0*fl[15]+133.1117951197414*fl[13]+34.3693177121688*fl[11]; 

  incr4[9] = 52.5*fl[9]+44.37059837324713*fl[5]+34.3693177121688*fl[2]+19.84313483298443*fl[0]; 
  incr4[12] = 52.5*fl[12]+44.37059837324713*fl[7]+34.3693177121688*fl[3]+19.84313483298443*fl[1]; 
  incr4[14] = 52.5*fl[14]+44.37059837324712*fl[10]+34.3693177121688*fl[6]+19.84313483298443*fl[4]; 
  incr4[15] = 52.5*fl[15]+44.37059837324712*fl[13]+34.3693177121688*fl[11]+19.84313483298443*fl[8]; 

  outl[2] += incr2[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul; 
  outl[5] += incr3[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += incr3[7]*rdxFnul-1.0*incr2[7]*rdxFnul; 
  outl[9] += incr4[9]*rdxFnul-1.0*incr3[9]*rdxFnul+incr2[9]*rdxFnul; 
  outl[10] += incr3[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 
  outl[12] += incr4[12]*rdxFnul-1.0*incr3[12]*rdxFnul+incr2[12]*rdxFnul; 
  outl[13] += incr3[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 
  outl[14] += incr4[14]*rdxFnul-1.0*incr3[14]*rdxFnul+incr2[14]*rdxFnul; 
  outl[15] += incr4[15]*rdxFnul-1.0*incr3[15]*rdxFnul+incr2[15]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf2xTensorP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 
  double incr5[16]; 
  double incr6[16]; 

  if (idxr[0] == 1) {


  incr3[4] = -133.1117951197414*fr[8]; 
  incr3[6] = -133.1117951197414*fr[11]; 
  incr3[8] = 787.5*fr[8]; 
  incr3[10] = -133.1117951197414*fr[13]; 
  incr3[11] = 787.5*fr[11]; 
  incr3[13] = 787.5*fr[13]; 
  incr3[14] = -133.1117951197414*fr[15]; 
  incr3[15] = 787.5*fr[15]; 

  incr4[8] = 133.1117951197414*fr[4]-787.5*fr[8]; 
  incr4[11] = 133.1117951197414*fr[6]-787.5*fr[11]; 
  incr4[13] = 133.1117951197414*fr[10]-787.5*fr[13]; 
  incr4[15] = 133.1117951197414*fr[14]-787.5*fr[15]; 



  outr[4] += incr3[4]*rdxFnur; 
  outr[6] += incr3[6]*rdxFnur; 
  outr[8] += incr4[8]*rdxFnur+incr3[8]*rdxFnur; 
  outr[10] += incr3[10]*rdxFnur; 
  outr[11] += incr4[11]*rdxFnur+incr3[11]*rdxFnur; 
  outr[13] += incr4[13]*rdxFnur+incr3[13]*rdxFnur; 
  outr[14] += incr3[14]*rdxFnur; 
  outr[15] += incr4[15]*rdxFnur+incr3[15]*rdxFnur; 

  } else {


  incr3[4] = -133.1117951197414*fl[8]; 
  incr3[6] = -133.1117951197414*fl[11]; 
  incr3[8] = 787.5*fl[8]; 
  incr3[10] = -133.1117951197414*fl[13]; 
  incr3[11] = 787.5*fl[11]; 
  incr3[13] = 787.5*fl[13]; 
  incr3[14] = -133.1117951197414*fl[15]; 
  incr3[15] = 787.5*fl[15]; 

  incr4[8] = 787.5*fl[8]+133.1117951197414*fl[4]; 
  incr4[11] = 787.5*fl[11]+133.1117951197414*fl[6]; 
  incr4[13] = 787.5*fl[13]+133.1117951197414*fl[10]; 
  incr4[15] = 787.5*fl[15]+133.1117951197414*fl[14]; 



  outl[4] += -1.0*incr3[4]*rdxFnul; 
  outl[6] += -1.0*incr3[6]*rdxFnul; 
  outl[8] += incr3[8]*rdxFnul-1.0*incr4[8]*rdxFnul; 
  outl[10] += -1.0*incr3[10]*rdxFnul; 
  outl[11] += incr3[11]*rdxFnul-1.0*incr4[11]*rdxFnul; 
  outl[13] += incr3[13]*rdxFnul-1.0*incr4[13]*rdxFnul; 
  outl[14] += -1.0*incr3[14]*rdxFnul; 
  outl[15] += incr3[15]*rdxFnul-1.0*incr4[15]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf2xTensorP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[2]:    current grid index.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[16]; 
  double incr2[16]; 
  double incr3[16]; 
  double incr4[16]; 
  double incr5[16]; 
  double incr6[16]; 

  if (idxr[1] == 1) {


  incr3[5] = -133.1117951197414*fr[9]; 
  incr3[7] = -133.1117951197414*fr[12]; 
  incr3[9] = 787.5*fr[9]; 
  incr3[10] = -133.1117951197414*fr[14]; 
  incr3[12] = 787.5*fr[12]; 
  incr3[13] = -133.1117951197414*fr[15]; 
  incr3[14] = 787.5*fr[14]; 
  incr3[15] = 787.5*fr[15]; 

  incr4[9] = 133.1117951197414*fr[5]-787.5*fr[9]; 
  incr4[12] = 133.1117951197414*fr[7]-787.5*fr[12]; 
  incr4[14] = 133.1117951197414*fr[10]-787.5*fr[14]; 
  incr4[15] = 133.1117951197414*fr[13]-787.5*fr[15]; 



  outr[5] += incr3[5]*rdxFnur; 
  outr[7] += incr3[7]*rdxFnur; 
  outr[9] += incr4[9]*rdxFnur+incr3[9]*rdxFnur; 
  outr[10] += incr3[10]*rdxFnur; 
  outr[12] += incr4[12]*rdxFnur+incr3[12]*rdxFnur; 
  outr[13] += incr3[13]*rdxFnur; 
  outr[14] += incr4[14]*rdxFnur+incr3[14]*rdxFnur; 
  outr[15] += incr4[15]*rdxFnur+incr3[15]*rdxFnur; 

  } else {


  incr3[5] = -133.1117951197414*fl[9]; 
  incr3[7] = -133.1117951197414*fl[12]; 
  incr3[9] = 787.5*fl[9]; 
  incr3[10] = -133.1117951197414*fl[14]; 
  incr3[12] = 787.5*fl[12]; 
  incr3[13] = -133.1117951197414*fl[15]; 
  incr3[14] = 787.5*fl[14]; 
  incr3[15] = 787.5*fl[15]; 

  incr4[9] = 787.5*fl[9]+133.1117951197414*fl[5]; 
  incr4[12] = 787.5*fl[12]+133.1117951197414*fl[7]; 
  incr4[14] = 787.5*fl[14]+133.1117951197414*fl[10]; 
  incr4[15] = 787.5*fl[15]+133.1117951197414*fl[13]; 



  outl[5] += -1.0*incr3[5]*rdxFnul; 
  outl[7] += -1.0*incr3[7]*rdxFnul; 
  outl[9] += incr3[9]*rdxFnul-1.0*incr4[9]*rdxFnul; 
  outl[10] += -1.0*incr3[10]*rdxFnul; 
  outl[12] += incr3[12]*rdxFnul-1.0*incr4[12]*rdxFnul; 
  outl[13] += -1.0*incr3[13]*rdxFnul; 
  outl[14] += incr3[14]*rdxFnul-1.0*incr4[14]*rdxFnul; 
  outl[15] += incr3[15]*rdxFnul-1.0*incr4[15]*rdxFnul; 

  }

} 
