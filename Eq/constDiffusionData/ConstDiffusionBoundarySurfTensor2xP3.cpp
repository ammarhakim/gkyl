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


  incr2[1] = 11.78376607274358*fr[8]-11.61895003862225*fr[4]+4.5*fr[1]; 
  incr2[3] = 11.78376607274359*fr[11]-11.61895003862225*fr[6]+4.5*fr[3]; 
  incr2[4] = (-45.6383297553399*fr[8])+45.0*fr[4]-17.42842505793337*fr[1]; 
  incr2[6] = (-45.6383297553399*fr[11])+45.0*fr[6]-17.42842505793338*fr[3]; 
  incr2[7] = 11.78376607274359*fr[13]-11.61895003862225*fr[10]+4.5*fr[7]; 
  incr2[8] = 108.0*fr[8]-106.4894360957931*fr[4]+41.24318125460255*fr[1]; 
  incr2[10] = (-45.6383297553399*fr[13])+45.0*fr[10]-17.42842505793338*fr[7]; 
  incr2[11] = 108.0*fr[11]-106.4894360957931*fr[6]+41.24318125460256*fr[3]; 
  incr2[12] = 11.78376607274359*fr[15]-11.61895003862225*fr[14]+4.5*fr[12]; 
  incr2[13] = 108.0*fr[13]-106.4894360957931*fr[10]+41.24318125460256*fr[7]; 
  incr2[14] = (-45.6383297553399*fr[15])+45.0*fr[14]-17.42842505793338*fr[12]; 
  incr2[15] = 108.0*fr[15]-106.4894360957931*fr[14]+41.24318125460256*fr[12]; 


  incr4[8] = (-16.2*fr[8])+28.39718295887816*fr[4]-30.24499958670854*fr[1]+19.84313483298443*fr[0]; 
  incr4[11] = (-16.2*fr[11])+28.39718295887816*fr[6]-30.24499958670854*fr[3]+19.84313483298443*fr[2]; 
  incr4[13] = (-16.2*fr[13])+28.39718295887816*fr[10]-30.24499958670855*fr[7]+19.84313483298443*fr[5]; 
  incr4[15] = (-16.2*fr[15])+28.39718295887816*fr[14]-30.24499958670854*fr[12]+19.84313483298443*fr[9]; 

  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[8] += (-1.0*incr4[8]*rdxFnur)-1.0*incr2[8]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[11] += (-1.0*incr4[11]*rdxFnur)-1.0*incr2[11]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[13] += (-1.0*incr4[13]*rdxFnur)-1.0*incr2[13]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[15] += (-1.0*incr4[15]*rdxFnur)-1.0*incr2[15]*rdxFnur; 

  } else {


  incr2[1] = 11.78376607274358*fl[8]+11.61895003862225*fl[4]+4.5*fl[1]; 
  incr2[3] = 11.78376607274359*fl[11]+11.61895003862225*fl[6]+4.5*fl[3]; 
  incr2[4] = 45.6383297553399*fl[8]+45.0*fl[4]+17.42842505793337*fl[1]; 
  incr2[6] = 45.6383297553399*fl[11]+45.0*fl[6]+17.42842505793338*fl[3]; 
  incr2[7] = 11.78376607274359*fl[13]+11.61895003862225*fl[10]+4.5*fl[7]; 
  incr2[8] = 108.0*fl[8]+106.4894360957931*fl[4]+41.24318125460255*fl[1]; 
  incr2[10] = 45.6383297553399*fl[13]+45.0*fl[10]+17.42842505793338*fl[7]; 
  incr2[11] = 108.0*fl[11]+106.4894360957931*fl[6]+41.24318125460256*fl[3]; 
  incr2[12] = 11.78376607274359*fl[15]+11.61895003862225*fl[14]+4.5*fl[12]; 
  incr2[13] = 108.0*fl[13]+106.4894360957931*fl[10]+41.24318125460256*fl[7]; 
  incr2[14] = 45.6383297553399*fl[15]+45.0*fl[14]+17.42842505793338*fl[12]; 
  incr2[15] = 108.0*fl[15]+106.4894360957931*fl[14]+41.24318125460256*fl[12]; 


  incr4[8] = (-16.2*fl[8])-28.39718295887816*fl[4]-30.24499958670854*fl[1]-19.84313483298443*fl[0]; 
  incr4[11] = (-16.2*fl[11])-28.39718295887816*fl[6]-30.24499958670854*fl[3]-19.84313483298443*fl[2]; 
  incr4[13] = (-16.2*fl[13])-28.39718295887816*fl[10]-30.24499958670855*fl[7]-19.84313483298443*fl[5]; 
  incr4[15] = (-16.2*fl[15])-28.39718295887816*fl[14]-30.24499958670854*fl[12]-19.84313483298443*fl[9]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[8] += (-1.0*incr4[8]*rdxFnul)-1.0*incr2[8]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += (-1.0*incr4[11]*rdxFnul)-1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += (-1.0*incr4[13]*rdxFnul)-1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += (-1.0*incr4[15]*rdxFnul)-1.0*incr2[15]*rdxFnul; 

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


  incr2[2] = 11.78376607274358*fr[9]-11.61895003862225*fr[5]+4.5*fr[2]; 
  incr2[3] = 11.78376607274359*fr[12]-11.61895003862225*fr[7]+4.5*fr[3]; 
  incr2[5] = (-45.6383297553399*fr[9])+45.0*fr[5]-17.42842505793337*fr[2]; 
  incr2[6] = 11.78376607274359*fr[14]-11.61895003862225*fr[10]+4.5*fr[6]; 
  incr2[7] = (-45.6383297553399*fr[12])+45.0*fr[7]-17.42842505793338*fr[3]; 
  incr2[9] = 108.0*fr[9]-106.4894360957931*fr[5]+41.24318125460255*fr[2]; 
  incr2[10] = (-45.6383297553399*fr[14])+45.0*fr[10]-17.42842505793338*fr[6]; 
  incr2[11] = 11.78376607274359*fr[15]-11.61895003862225*fr[13]+4.5*fr[11]; 
  incr2[12] = 108.0*fr[12]-106.4894360957931*fr[7]+41.24318125460256*fr[3]; 
  incr2[13] = (-45.6383297553399*fr[15])+45.0*fr[13]-17.42842505793338*fr[11]; 
  incr2[14] = 108.0*fr[14]-106.4894360957931*fr[10]+41.24318125460256*fr[6]; 
  incr2[15] = 108.0*fr[15]-106.4894360957931*fr[13]+41.24318125460256*fr[11]; 


  incr4[9] = (-16.2*fr[9])+28.39718295887816*fr[5]-30.24499958670854*fr[2]+19.84313483298443*fr[0]; 
  incr4[12] = (-16.2*fr[12])+28.39718295887816*fr[7]-30.24499958670854*fr[3]+19.84313483298443*fr[1]; 
  incr4[14] = (-16.2*fr[14])+28.39718295887816*fr[10]-30.24499958670855*fr[6]+19.84313483298443*fr[4]; 
  incr4[15] = (-16.2*fr[15])+28.39718295887816*fr[13]-30.24499958670854*fr[11]+19.84313483298443*fr[8]; 

  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[9] += (-1.0*incr4[9]*rdxFnur)-1.0*incr2[9]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[12] += (-1.0*incr4[12]*rdxFnur)-1.0*incr2[12]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[14] += (-1.0*incr4[14]*rdxFnur)-1.0*incr2[14]*rdxFnur; 
  outr[15] += (-1.0*incr4[15]*rdxFnur)-1.0*incr2[15]*rdxFnur; 

  } else {


  incr2[2] = 11.78376607274358*fl[9]+11.61895003862225*fl[5]+4.5*fl[2]; 
  incr2[3] = 11.78376607274359*fl[12]+11.61895003862225*fl[7]+4.5*fl[3]; 
  incr2[5] = 45.6383297553399*fl[9]+45.0*fl[5]+17.42842505793337*fl[2]; 
  incr2[6] = 11.78376607274359*fl[14]+11.61895003862225*fl[10]+4.5*fl[6]; 
  incr2[7] = 45.6383297553399*fl[12]+45.0*fl[7]+17.42842505793338*fl[3]; 
  incr2[9] = 108.0*fl[9]+106.4894360957931*fl[5]+41.24318125460255*fl[2]; 
  incr2[10] = 45.6383297553399*fl[14]+45.0*fl[10]+17.42842505793338*fl[6]; 
  incr2[11] = 11.78376607274359*fl[15]+11.61895003862225*fl[13]+4.5*fl[11]; 
  incr2[12] = 108.0*fl[12]+106.4894360957931*fl[7]+41.24318125460256*fl[3]; 
  incr2[13] = 45.6383297553399*fl[15]+45.0*fl[13]+17.42842505793338*fl[11]; 
  incr2[14] = 108.0*fl[14]+106.4894360957931*fl[10]+41.24318125460256*fl[6]; 
  incr2[15] = 108.0*fl[15]+106.4894360957931*fl[13]+41.24318125460256*fl[11]; 


  incr4[9] = (-16.2*fl[9])-28.39718295887816*fl[5]-30.24499958670854*fl[2]-19.84313483298443*fl[0]; 
  incr4[12] = (-16.2*fl[12])-28.39718295887816*fl[7]-30.24499958670854*fl[3]-19.84313483298443*fl[1]; 
  incr4[14] = (-16.2*fl[14])-28.39718295887816*fl[10]-30.24499958670855*fl[6]-19.84313483298443*fl[4]; 
  incr4[15] = (-16.2*fl[15])-28.39718295887816*fl[13]-30.24499958670854*fl[11]-19.84313483298443*fl[8]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[9] += (-1.0*incr4[9]*rdxFnul)-1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += (-1.0*incr4[12]*rdxFnul)-1.0*incr2[12]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[14] += (-1.0*incr4[14]*rdxFnul)-1.0*incr2[14]*rdxFnul; 
  outl[15] += (-1.0*incr4[15]*rdxFnul)-1.0*incr2[15]*rdxFnul; 

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


  incr2[1] = (-103.1079531365064*fr[8])+76.2493596284585*fr[4]-19.6875*fr[1]; 
  incr2[3] = (-103.1079531365064*fr[11])+76.24935962845854*fr[6]-19.6875*fr[3]; 
  incr2[4] = 399.3353853592242*fr[8]-295.3125*fr[4]+76.2493596284585*fr[1]; 
  incr2[6] = 399.3353853592241*fr[11]-295.3125*fr[6]+76.24935962845854*fr[3]; 
  incr2[7] = (-103.1079531365064*fr[13])+76.24935962845854*fr[10]-19.6875*fr[7]; 
  incr2[8] = (-945.0*fr[8])+698.8369243786424*fr[4]-180.4389179888861*fr[1]; 
  incr2[10] = 399.3353853592241*fr[13]-295.3125*fr[10]+76.24935962845854*fr[7]; 
  incr2[11] = (-945.0*fr[11])+698.8369243786422*fr[6]-180.4389179888862*fr[3]; 
  incr2[12] = (-103.1079531365064*fr[15])+76.24935962845852*fr[14]-19.6875*fr[12]; 
  incr2[13] = (-945.0*fr[13])+698.8369243786422*fr[10]-180.4389179888862*fr[7]; 
  incr2[14] = 399.3353853592241*fr[15]-295.3125*fr[14]+76.24935962845852*fr[12]; 
  incr2[15] = (-945.0*fr[15])+698.8369243786422*fr[14]-180.4389179888862*fr[12]; 


  incr4[8] = 225.0*fr[8]-241.2651286545313*fr[4]+96.66370606547471*fr[1]; 
  incr4[11] = 225.0*fr[11]-241.2651286545313*fr[6]+96.66370606547474*fr[3]; 
  incr4[13] = 225.0*fr[13]-241.2651286545312*fr[10]+96.66370606547476*fr[7]; 
  incr4[15] = 225.0*fr[15]-241.2651286545312*fr[14]+96.66370606547474*fr[12]; 



  outr[1] += incr2[1]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[8] += incr4[8]*rdxFnur+incr2[8]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr4[11]*rdxFnur+incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[13] += incr4[13]*rdxFnur+incr2[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr4[15]*rdxFnur+incr2[15]*rdxFnur; 

  } else {


  incr2[1] = (-103.1079531365064*fl[8])-76.2493596284585*fl[4]-19.6875*fl[1]; 
  incr2[3] = (-103.1079531365064*fl[11])-76.24935962845854*fl[6]-19.6875*fl[3]; 
  incr2[4] = (-399.3353853592242*fl[8])-295.3125*fl[4]-76.2493596284585*fl[1]; 
  incr2[6] = (-399.3353853592241*fl[11])-295.3125*fl[6]-76.24935962845854*fl[3]; 
  incr2[7] = (-103.1079531365064*fl[13])-76.24935962845854*fl[10]-19.6875*fl[7]; 
  incr2[8] = (-945.0*fl[8])-698.8369243786424*fl[4]-180.4389179888861*fl[1]; 
  incr2[10] = (-399.3353853592241*fl[13])-295.3125*fl[10]-76.24935962845854*fl[7]; 
  incr2[11] = (-945.0*fl[11])-698.8369243786422*fl[6]-180.4389179888862*fl[3]; 
  incr2[12] = (-103.1079531365064*fl[15])-76.24935962845852*fl[14]-19.6875*fl[12]; 
  incr2[13] = (-945.0*fl[13])-698.8369243786422*fl[10]-180.4389179888862*fl[7]; 
  incr2[14] = (-399.3353853592241*fl[15])-295.3125*fl[14]-76.24935962845852*fl[12]; 
  incr2[15] = (-945.0*fl[15])-698.8369243786422*fl[14]-180.4389179888862*fl[12]; 


  incr4[8] = 225.0*fl[8]+241.2651286545313*fl[4]+96.66370606547471*fl[1]; 
  incr4[11] = 225.0*fl[11]+241.2651286545313*fl[6]+96.66370606547474*fl[3]; 
  incr4[13] = 225.0*fl[13]+241.2651286545312*fl[10]+96.66370606547476*fl[7]; 
  incr4[15] = 225.0*fl[15]+241.2651286545312*fl[14]+96.66370606547474*fl[12]; 



  outl[1] += incr2[1]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[8] += incr4[8]*rdxFnul+incr2[8]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[11] += incr4[11]*rdxFnul+incr2[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 
  outl[13] += incr4[13]*rdxFnul+incr2[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 
  outl[15] += incr4[15]*rdxFnul+incr2[15]*rdxFnul; 

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


  incr2[2] = (-103.1079531365064*fr[9])+76.2493596284585*fr[5]-19.6875*fr[2]; 
  incr2[3] = (-103.1079531365064*fr[12])+76.24935962845854*fr[7]-19.6875*fr[3]; 
  incr2[5] = 399.3353853592242*fr[9]-295.3125*fr[5]+76.2493596284585*fr[2]; 
  incr2[6] = (-103.1079531365064*fr[14])+76.24935962845854*fr[10]-19.6875*fr[6]; 
  incr2[7] = 399.3353853592241*fr[12]-295.3125*fr[7]+76.24935962845854*fr[3]; 
  incr2[9] = (-945.0*fr[9])+698.8369243786424*fr[5]-180.4389179888861*fr[2]; 
  incr2[10] = 399.3353853592241*fr[14]-295.3125*fr[10]+76.24935962845854*fr[6]; 
  incr2[11] = (-103.1079531365064*fr[15])+76.24935962845852*fr[13]-19.6875*fr[11]; 
  incr2[12] = (-945.0*fr[12])+698.8369243786422*fr[7]-180.4389179888862*fr[3]; 
  incr2[13] = 399.3353853592241*fr[15]-295.3125*fr[13]+76.24935962845852*fr[11]; 
  incr2[14] = (-945.0*fr[14])+698.8369243786422*fr[10]-180.4389179888862*fr[6]; 
  incr2[15] = (-945.0*fr[15])+698.8369243786422*fr[13]-180.4389179888862*fr[11]; 


  incr4[9] = 225.0*fr[9]-241.2651286545313*fr[5]+96.66370606547471*fr[2]; 
  incr4[12] = 225.0*fr[12]-241.2651286545313*fr[7]+96.66370606547474*fr[3]; 
  incr4[14] = 225.0*fr[14]-241.2651286545312*fr[10]+96.66370606547476*fr[6]; 
  incr4[15] = 225.0*fr[15]-241.2651286545312*fr[13]+96.66370606547474*fr[11]; 



  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[9] += incr4[9]*rdxFnur+incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr4[12]*rdxFnur+incr2[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[14] += incr4[14]*rdxFnur+incr2[14]*rdxFnur; 
  outr[15] += incr4[15]*rdxFnur+incr2[15]*rdxFnur; 

  } else {


  incr2[2] = (-103.1079531365064*fl[9])-76.2493596284585*fl[5]-19.6875*fl[2]; 
  incr2[3] = (-103.1079531365064*fl[12])-76.24935962845854*fl[7]-19.6875*fl[3]; 
  incr2[5] = (-399.3353853592242*fl[9])-295.3125*fl[5]-76.2493596284585*fl[2]; 
  incr2[6] = (-103.1079531365064*fl[14])-76.24935962845854*fl[10]-19.6875*fl[6]; 
  incr2[7] = (-399.3353853592241*fl[12])-295.3125*fl[7]-76.24935962845854*fl[3]; 
  incr2[9] = (-945.0*fl[9])-698.8369243786424*fl[5]-180.4389179888861*fl[2]; 
  incr2[10] = (-399.3353853592241*fl[14])-295.3125*fl[10]-76.24935962845854*fl[6]; 
  incr2[11] = (-103.1079531365064*fl[15])-76.24935962845852*fl[13]-19.6875*fl[11]; 
  incr2[12] = (-945.0*fl[12])-698.8369243786422*fl[7]-180.4389179888862*fl[3]; 
  incr2[13] = (-399.3353853592241*fl[15])-295.3125*fl[13]-76.24935962845852*fl[11]; 
  incr2[14] = (-945.0*fl[14])-698.8369243786422*fl[10]-180.4389179888862*fl[6]; 
  incr2[15] = (-945.0*fl[15])-698.8369243786422*fl[13]-180.4389179888862*fl[11]; 


  incr4[9] = 225.0*fl[9]+241.2651286545313*fl[5]+96.66370606547471*fl[2]; 
  incr4[12] = 225.0*fl[12]+241.2651286545313*fl[7]+96.66370606547474*fl[3]; 
  incr4[14] = 225.0*fl[14]+241.2651286545312*fl[10]+96.66370606547476*fl[6]; 
  incr4[15] = 225.0*fl[15]+241.2651286545312*fl[13]+96.66370606547474*fl[11]; 



  outl[2] += incr2[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[9] += incr4[9]*rdxFnul+incr2[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 
  outl[12] += incr4[12]*rdxFnul+incr2[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 
  outl[14] += incr4[14]*rdxFnul+incr2[14]*rdxFnul; 
  outl[15] += incr4[15]*rdxFnul+incr2[15]*rdxFnul; 

  }

} 
