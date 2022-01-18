#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf2xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[12]; 
  double incr2[12]; 

  if (edge < 0) {

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
void ConstDiffusionBoundarySurf2xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[12]; 
  double incr2[12]; 

  if (edge < 0) {

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
void ConstHyperDiffusion4BoundarySurf2xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[12]; 
  double incr2[12]; 
  double incr3[12]; 
  double incr4[12]; 

  if (edge < 0) {


  incr2[1] = 11.78376607274358*fr[8]-11.61895003862225*fr[4]+4.5*fr[1]; 
  incr2[3] = 11.78376607274359*fr[10]-11.61895003862225*fr[6]+4.5*fr[3]; 
  incr2[4] = (-45.6383297553399*fr[8])+45.0*fr[4]-17.42842505793337*fr[1]; 
  incr2[6] = (-45.6383297553399*fr[10])+45.0*fr[6]-17.42842505793338*fr[3]; 
  incr2[7] = 4.5*fr[7]; 
  incr2[8] = 108.0*fr[8]-106.4894360957931*fr[4]+41.24318125460255*fr[1]; 
  incr2[10] = 108.0*fr[10]-106.4894360957931*fr[6]+41.24318125460256*fr[3]; 
  incr2[11] = 4.5*fr[11]; 


  incr4[8] = (-16.2*fr[8])+28.39718295887816*fr[4]-30.24499958670854*fr[1]+19.84313483298443*fr[0]; 
  incr4[10] = (-16.2*fr[10])+28.39718295887816*fr[6]-30.24499958670854*fr[3]+19.84313483298443*fr[2]; 

  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[8] += (-1.0*incr4[8]*rdxFnur)-1.0*incr2[8]*rdxFnur; 
  outr[10] += (-1.0*incr4[10]*rdxFnur)-1.0*incr2[10]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 

  } else {


  incr2[1] = 11.78376607274358*fl[8]+11.61895003862225*fl[4]+4.5*fl[1]; 
  incr2[3] = 11.78376607274359*fl[10]+11.61895003862225*fl[6]+4.5*fl[3]; 
  incr2[4] = 45.6383297553399*fl[8]+45.0*fl[4]+17.42842505793337*fl[1]; 
  incr2[6] = 45.6383297553399*fl[10]+45.0*fl[6]+17.42842505793338*fl[3]; 
  incr2[7] = 4.5*fl[7]; 
  incr2[8] = 108.0*fl[8]+106.4894360957931*fl[4]+41.24318125460255*fl[1]; 
  incr2[10] = 108.0*fl[10]+106.4894360957931*fl[6]+41.24318125460256*fl[3]; 
  incr2[11] = 4.5*fl[11]; 


  incr4[8] = (-16.2*fl[8])-28.39718295887816*fl[4]-30.24499958670854*fl[1]-19.84313483298443*fl[0]; 
  incr4[10] = (-16.2*fl[10])-28.39718295887816*fl[6]-30.24499958670854*fl[3]-19.84313483298443*fl[2]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[8] += (-1.0*incr4[8]*rdxFnul)-1.0*incr2[8]*rdxFnul; 
  outl[10] += (-1.0*incr4[10]*rdxFnul)-1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf2xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[12]; 
  double incr2[12]; 
  double incr3[12]; 
  double incr4[12]; 

  if (edge < 0) {


  incr2[2] = 11.78376607274358*fr[9]-11.61895003862225*fr[5]+4.5*fr[2]; 
  incr2[3] = 11.78376607274359*fr[11]-11.61895003862225*fr[7]+4.5*fr[3]; 
  incr2[5] = (-45.6383297553399*fr[9])+45.0*fr[5]-17.42842505793337*fr[2]; 
  incr2[6] = 4.5*fr[6]; 
  incr2[7] = (-45.6383297553399*fr[11])+45.0*fr[7]-17.42842505793338*fr[3]; 
  incr2[9] = 108.0*fr[9]-106.4894360957931*fr[5]+41.24318125460255*fr[2]; 
  incr2[10] = 4.5*fr[10]; 
  incr2[11] = 108.0*fr[11]-106.4894360957931*fr[7]+41.24318125460256*fr[3]; 


  incr4[9] = (-16.2*fr[9])+28.39718295887816*fr[5]-30.24499958670854*fr[2]+19.84313483298443*fr[0]; 
  incr4[11] = (-16.2*fr[11])+28.39718295887816*fr[7]-30.24499958670854*fr[3]+19.84313483298443*fr[1]; 

  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[9] += (-1.0*incr4[9]*rdxFnur)-1.0*incr2[9]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[11] += (-1.0*incr4[11]*rdxFnur)-1.0*incr2[11]*rdxFnur; 

  } else {


  incr2[2] = 11.78376607274358*fl[9]+11.61895003862225*fl[5]+4.5*fl[2]; 
  incr2[3] = 11.78376607274359*fl[11]+11.61895003862225*fl[7]+4.5*fl[3]; 
  incr2[5] = 45.6383297553399*fl[9]+45.0*fl[5]+17.42842505793337*fl[2]; 
  incr2[6] = 4.5*fl[6]; 
  incr2[7] = 45.6383297553399*fl[11]+45.0*fl[7]+17.42842505793338*fl[3]; 
  incr2[9] = 108.0*fl[9]+106.4894360957931*fl[5]+41.24318125460255*fl[2]; 
  incr2[10] = 4.5*fl[10]; 
  incr2[11] = 108.0*fl[11]+106.4894360957931*fl[7]+41.24318125460256*fl[3]; 


  incr4[9] = (-16.2*fl[9])-28.39718295887816*fl[5]-30.24499958670854*fl[2]-19.84313483298443*fl[0]; 
  incr4[11] = (-16.2*fl[11])-28.39718295887816*fl[7]-30.24499958670854*fl[3]-19.84313483298443*fl[1]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[9] += (-1.0*incr4[9]*rdxFnul)-1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += (-1.0*incr4[11]*rdxFnul)-1.0*incr2[11]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf2xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
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

  if (edge < 0) {


  incr2[1] = (-103.1079531365064*fr[8])+76.2493596284585*fr[4]-19.6875*fr[1]; 
  incr2[3] = (-103.1079531365064*fr[10])+76.24935962845854*fr[6]-19.6875*fr[3]; 
  incr2[4] = 399.3353853592242*fr[8]-295.3125*fr[4]+76.2493596284585*fr[1]; 
  incr2[6] = 399.3353853592241*fr[10]-295.3125*fr[6]+76.24935962845854*fr[3]; 
  incr2[7] = -19.6875*fr[7]; 
  incr2[8] = (-945.0*fr[8])+698.8369243786424*fr[4]-180.4389179888861*fr[1]; 
  incr2[10] = (-945.0*fr[10])+698.8369243786422*fr[6]-180.4389179888862*fr[3]; 
  incr2[11] = -19.6875*fr[11]; 


  incr4[8] = 225.0*fr[8]-241.2651286545313*fr[4]+96.66370606547471*fr[1]; 
  incr4[10] = 225.0*fr[10]-241.2651286545313*fr[6]+96.66370606547474*fr[3]; 



  outr[1] += incr2[1]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[8] += incr4[8]*rdxFnur+incr2[8]*rdxFnur; 
  outr[10] += incr4[10]*rdxFnur+incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 

  } else {


  incr2[1] = (-103.1079531365064*fl[8])-76.2493596284585*fl[4]-19.6875*fl[1]; 
  incr2[3] = (-103.1079531365064*fl[10])-76.24935962845854*fl[6]-19.6875*fl[3]; 
  incr2[4] = (-399.3353853592242*fl[8])-295.3125*fl[4]-76.2493596284585*fl[1]; 
  incr2[6] = (-399.3353853592241*fl[10])-295.3125*fl[6]-76.24935962845854*fl[3]; 
  incr2[7] = -19.6875*fl[7]; 
  incr2[8] = (-945.0*fl[8])-698.8369243786424*fl[4]-180.4389179888861*fl[1]; 
  incr2[10] = (-945.0*fl[10])-698.8369243786422*fl[6]-180.4389179888862*fl[3]; 
  incr2[11] = -19.6875*fl[11]; 


  incr4[8] = 225.0*fl[8]+241.2651286545313*fl[4]+96.66370606547471*fl[1]; 
  incr4[10] = 225.0*fl[10]+241.2651286545313*fl[6]+96.66370606547474*fl[3]; 



  outl[1] += incr2[1]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[8] += incr4[8]*rdxFnul+incr2[8]*rdxFnul; 
  outl[10] += incr4[10]*rdxFnul+incr2[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf2xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
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

  if (edge < 0) {


  incr2[2] = (-103.1079531365064*fr[9])+76.2493596284585*fr[5]-19.6875*fr[2]; 
  incr2[3] = (-103.1079531365064*fr[11])+76.24935962845854*fr[7]-19.6875*fr[3]; 
  incr2[5] = 399.3353853592242*fr[9]-295.3125*fr[5]+76.2493596284585*fr[2]; 
  incr2[6] = -19.6875*fr[6]; 
  incr2[7] = 399.3353853592241*fr[11]-295.3125*fr[7]+76.24935962845854*fr[3]; 
  incr2[9] = (-945.0*fr[9])+698.8369243786424*fr[5]-180.4389179888861*fr[2]; 
  incr2[10] = -19.6875*fr[10]; 
  incr2[11] = (-945.0*fr[11])+698.8369243786422*fr[7]-180.4389179888862*fr[3]; 


  incr4[9] = 225.0*fr[9]-241.2651286545313*fr[5]+96.66370606547471*fr[2]; 
  incr4[11] = 225.0*fr[11]-241.2651286545313*fr[7]+96.66370606547474*fr[3]; 



  outr[2] += incr2[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[9] += incr4[9]*rdxFnur+incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr4[11]*rdxFnur+incr2[11]*rdxFnur; 

  } else {


  incr2[2] = (-103.1079531365064*fl[9])-76.2493596284585*fl[5]-19.6875*fl[2]; 
  incr2[3] = (-103.1079531365064*fl[11])-76.24935962845854*fl[7]-19.6875*fl[3]; 
  incr2[5] = (-399.3353853592242*fl[9])-295.3125*fl[5]-76.2493596284585*fl[2]; 
  incr2[6] = -19.6875*fl[6]; 
  incr2[7] = (-399.3353853592241*fl[11])-295.3125*fl[7]-76.24935962845854*fl[3]; 
  incr2[9] = (-945.0*fl[9])-698.8369243786424*fl[5]-180.4389179888861*fl[2]; 
  incr2[10] = -19.6875*fl[10]; 
  incr2[11] = (-945.0*fl[11])-698.8369243786422*fl[7]-180.4389179888862*fl[3]; 


  incr4[9] = 225.0*fl[9]+241.2651286545313*fl[5]+96.66370606547471*fl[2]; 
  incr4[11] = 225.0*fl[11]+241.2651286545313*fl[7]+96.66370606547474*fl[3]; 



  outl[2] += incr2[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[9] += incr4[9]*rdxFnul+incr2[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[11] += incr4[11]*rdxFnul+incr2[11]*rdxFnul; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf2xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[24]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[12]; 
  double incr2[12]; 

  if (edge < 0) {

  incr2[1] = (-3.031088913245536*fr[8]*nul[8])+2.561737691489899*fr[4]*nul[8]-1.984313483298443*fr[1]*nul[8]+1.14564392373896*fr[0]*nul[8]-2.561737691489901*nul[4]*fr[8]-1.984313483298443*nul[1]*fr[8]-1.14564392373896*nul[0]*fr[8]+2.165063509461096*fr[4]*nul[4]-1.677050983124842*fr[1]*nul[4]+0.9682458365518543*fr[0]*nul[4]+1.677050983124842*nul[1]*fr[4]+0.9682458365518541*nul[0]*fr[4]-1.299038105676658*fr[1]*nul[1]+0.75*fr[0]*nul[1]-0.75*nul[0]*fr[1]+0.4330127018922193*fr[0]*nul[0]; 
  incr2[3] = (-3.031088913245535*nul[8]*fr[10])-2.561737691489901*nul[4]*fr[10]-1.984313483298442*nul[1]*fr[10]-1.14564392373896*nul[0]*fr[10]+2.5617376914899*fr[6]*nul[8]-1.984313483298443*fr[3]*nul[8]+1.14564392373896*fr[2]*nul[8]+2.165063509461097*nul[4]*fr[6]+1.677050983124842*nul[1]*fr[6]+0.9682458365518543*nul[0]*fr[6]-1.677050983124842*fr[3]*nul[4]+0.9682458365518543*fr[2]*nul[4]-1.299038105676658*nul[1]*fr[3]-0.75*nul[0]*fr[3]+0.75*nul[1]*fr[2]+0.4330127018922193*nul[0]*fr[2]; 
  incr2[4] = 11.7393568818739*fr[8]*nul[8]-9.921567416492215*fr[4]*nul[8]+7.685213074469699*fr[1]*nul[8]-4.437059837324712*fr[0]*nul[8]+9.921567416492215*nul[4]*fr[8]+7.685213074469706*nul[1]*fr[8]+4.437059837324713*nul[0]*fr[8]-8.385254915624213*fr[4]*nul[4]+6.495190528383289*fr[1]*nul[4]-3.75*fr[0]*nul[4]-6.495190528383286*nul[1]*fr[4]-3.75*nul[0]*fr[4]+5.031152949374527*fr[1]*nul[1]-2.904737509655563*fr[0]*nul[1]+2.904737509655563*nul[0]*fr[1]-1.677050983124842*fr[0]*nul[0]; 
  incr2[6] = 11.7393568818739*nul[8]*fr[10]+9.921567416492222*nul[4]*fr[10]+7.685213074469699*nul[1]*fr[10]+4.437059837324712*nul[0]*fr[10]-9.921567416492215*fr[6]*nul[8]+7.6852130744697*fr[3]*nul[8]-4.437059837324712*fr[2]*nul[8]-8.385254915624213*nul[4]*fr[6]-6.495190528383286*nul[1]*fr[6]-3.75*nul[0]*fr[6]+6.49519052838329*fr[3]*nul[4]-3.75*fr[2]*nul[4]+5.031152949374526*nul[1]*fr[3]+2.904737509655563*nul[0]*fr[3]-2.904737509655563*nul[1]*fr[2]-1.677050983124842*nul[0]*fr[2]; 
  incr2[7] = (-1.984313483298443*fr[7]*nul[8])+1.14564392373896*fr[5]*nul[8]-1.677050983124842*nul[4]*fr[7]-1.299038105676658*nul[1]*fr[7]-0.75*nul[0]*fr[7]+0.9682458365518543*nul[4]*fr[5]+0.75*nul[1]*fr[5]+0.4330127018922194*nul[0]*fr[5]; 
  incr2[8] = (-27.78038876617821*fr[8]*nul[8])+23.47871376374779*fr[4]*nul[8]-18.18653347947321*fr[1]*nul[8]+10.5*fr[0]*nul[8]-23.47871376374781*nul[4]*fr[8]-18.18653347947322*nul[1]*fr[8]-10.5*nul[0]*fr[8]+19.84313483298443*fr[4]*nul[4]-15.3704261489394*fr[1]*nul[4]+8.874119674649425*fr[0]*nul[4]+15.37042614893939*nul[1]*fr[4]+8.874119674649425*nul[0]*fr[4]-11.90588089979066*fr[1]*nul[1]+6.873863542433759*fr[0]*nul[1]-6.873863542433759*nul[0]*fr[1]+3.968626966596886*fr[0]*nul[0]; 
  incr2[10] = (-27.78038876617821*nul[8]*fr[10])-23.47871376374781*nul[4]*fr[10]-18.18653347947322*nul[1]*fr[10]-10.5*nul[0]*fr[10]+23.47871376374779*fr[6]*nul[8]-18.18653347947321*fr[3]*nul[8]+10.5*fr[2]*nul[8]+19.84313483298443*nul[4]*fr[6]+15.37042614893939*nul[1]*fr[6]+8.874119674649425*nul[0]*fr[6]-15.3704261489394*fr[3]*nul[4]+8.874119674649425*fr[2]*nul[4]-11.90588089979066*nul[1]*fr[3]-6.87386354243376*nul[0]*fr[3]+6.87386354243376*nul[1]*fr[2]+3.968626966596886*nul[0]*fr[2]; 
  incr2[11] = (-1.984313483298443*nul[8]*fr[11])-1.677050983124842*nul[4]*fr[11]-1.299038105676658*nul[1]*fr[11]-0.75*nul[0]*fr[11]+1.14564392373896*nul[8]*fr[9]+0.9682458365518541*nul[4]*fr[9]+0.7499999999999999*nul[1]*fr[9]+0.4330127018922193*nul[0]*fr[9]; 

  outr[1] += incr2[1]*rdxFr; 
  outr[3] += incr2[3]*rdxFr; 
  outr[4] += incr2[4]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 
  outr[8] += incr2[8]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 

  } else {

  incr2[1] = 3.031088913245536*fl[8]*nul[8]+2.561737691489899*fl[4]*nul[8]+1.984313483298443*fl[1]*nul[8]+1.14564392373896*fl[0]*nul[8]+2.561737691489901*nul[4]*fl[8]+1.984313483298443*nul[1]*fl[8]+1.14564392373896*nul[0]*fl[8]+2.165063509461096*fl[4]*nul[4]+1.677050983124842*fl[1]*nul[4]+0.9682458365518543*fl[0]*nul[4]+1.677050983124842*nul[1]*fl[4]+0.9682458365518541*nul[0]*fl[4]+1.299038105676658*fl[1]*nul[1]+0.75*fl[0]*nul[1]+0.75*nul[0]*fl[1]+0.4330127018922193*fl[0]*nul[0]; 
  incr2[3] = 3.031088913245535*nul[8]*fl[10]+2.561737691489901*nul[4]*fl[10]+1.984313483298442*nul[1]*fl[10]+1.14564392373896*nul[0]*fl[10]+2.5617376914899*fl[6]*nul[8]+1.984313483298443*fl[3]*nul[8]+1.14564392373896*fl[2]*nul[8]+2.165063509461097*nul[4]*fl[6]+1.677050983124842*nul[1]*fl[6]+0.9682458365518543*nul[0]*fl[6]+1.677050983124842*fl[3]*nul[4]+0.9682458365518543*fl[2]*nul[4]+1.299038105676658*nul[1]*fl[3]+0.75*nul[0]*fl[3]+0.75*nul[1]*fl[2]+0.4330127018922193*nul[0]*fl[2]; 
  incr2[4] = (-11.7393568818739*fl[8]*nul[8])-9.921567416492215*fl[4]*nul[8]-7.685213074469699*fl[1]*nul[8]-4.437059837324712*fl[0]*nul[8]-9.921567416492215*nul[4]*fl[8]-7.685213074469706*nul[1]*fl[8]-4.437059837324713*nul[0]*fl[8]-8.385254915624213*fl[4]*nul[4]-6.495190528383289*fl[1]*nul[4]-3.75*fl[0]*nul[4]-6.495190528383286*nul[1]*fl[4]-3.75*nul[0]*fl[4]-5.031152949374527*fl[1]*nul[1]-2.904737509655563*fl[0]*nul[1]-2.904737509655563*nul[0]*fl[1]-1.677050983124842*fl[0]*nul[0]; 
  incr2[6] = (-11.7393568818739*nul[8]*fl[10])-9.921567416492222*nul[4]*fl[10]-7.685213074469699*nul[1]*fl[10]-4.437059837324712*nul[0]*fl[10]-9.921567416492215*fl[6]*nul[8]-7.6852130744697*fl[3]*nul[8]-4.437059837324712*fl[2]*nul[8]-8.385254915624213*nul[4]*fl[6]-6.495190528383286*nul[1]*fl[6]-3.75*nul[0]*fl[6]-6.49519052838329*fl[3]*nul[4]-3.75*fl[2]*nul[4]-5.031152949374526*nul[1]*fl[3]-2.904737509655563*nul[0]*fl[3]-2.904737509655563*nul[1]*fl[2]-1.677050983124842*nul[0]*fl[2]; 
  incr2[7] = 1.984313483298443*fl[7]*nul[8]+1.14564392373896*fl[5]*nul[8]+1.677050983124842*nul[4]*fl[7]+1.299038105676658*nul[1]*fl[7]+0.75*nul[0]*fl[7]+0.9682458365518543*nul[4]*fl[5]+0.75*nul[1]*fl[5]+0.4330127018922194*nul[0]*fl[5]; 
  incr2[8] = 27.78038876617821*fl[8]*nul[8]+23.47871376374779*fl[4]*nul[8]+18.18653347947321*fl[1]*nul[8]+10.5*fl[0]*nul[8]+23.47871376374781*nul[4]*fl[8]+18.18653347947322*nul[1]*fl[8]+10.5*nul[0]*fl[8]+19.84313483298443*fl[4]*nul[4]+15.3704261489394*fl[1]*nul[4]+8.874119674649425*fl[0]*nul[4]+15.37042614893939*nul[1]*fl[4]+8.874119674649425*nul[0]*fl[4]+11.90588089979066*fl[1]*nul[1]+6.873863542433759*fl[0]*nul[1]+6.873863542433759*nul[0]*fl[1]+3.968626966596886*fl[0]*nul[0]; 
  incr2[10] = 27.78038876617821*nul[8]*fl[10]+23.47871376374781*nul[4]*fl[10]+18.18653347947322*nul[1]*fl[10]+10.5*nul[0]*fl[10]+23.47871376374779*fl[6]*nul[8]+18.18653347947321*fl[3]*nul[8]+10.5*fl[2]*nul[8]+19.84313483298443*nul[4]*fl[6]+15.37042614893939*nul[1]*fl[6]+8.874119674649425*nul[0]*fl[6]+15.3704261489394*fl[3]*nul[4]+8.874119674649425*fl[2]*nul[4]+11.90588089979066*nul[1]*fl[3]+6.87386354243376*nul[0]*fl[3]+6.87386354243376*nul[1]*fl[2]+3.968626966596886*nul[0]*fl[2]; 
  incr2[11] = 1.984313483298443*nul[8]*fl[11]+1.677050983124842*nul[4]*fl[11]+1.299038105676658*nul[1]*fl[11]+0.75*nul[0]*fl[11]+1.14564392373896*nul[8]*fl[9]+0.9682458365518541*nul[4]*fl[9]+0.7499999999999999*nul[1]*fl[9]+0.4330127018922193*nul[0]*fl[9]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 
  outl[3] += -1.0*incr2[3]*rdxFl; 
  outl[4] += incr2[4]*rdxFl; 
  outl[6] += incr2[6]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 
  outl[8] += -1.0*incr2[8]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[11] += -1.0*incr2[11]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf2xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[24]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[1]*dxl[1]); 
  double rdxFr = 4.0/(dxr[1]*dxr[1]); 

  double incr1[12]; 
  double incr2[12]; 

  if (edge < 0) {

  incr2[2] = (-3.031088913245536*fr[9]*nul[21])+2.561737691489899*fr[5]*nul[21]-1.984313483298443*fr[2]*nul[21]+1.14564392373896*fr[0]*nul[21]-2.561737691489901*fr[9]*nul[17]+2.165063509461096*fr[5]*nul[17]-1.677050983124842*fr[2]*nul[17]+0.9682458365518543*fr[0]*nul[17]-1.984313483298443*fr[9]*nul[14]+1.677050983124842*fr[5]*nul[14]-1.299038105676658*fr[2]*nul[14]+0.75*fr[0]*nul[14]-1.14564392373896*fr[9]*nul[12]+0.9682458365518541*fr[5]*nul[12]-0.75*fr[2]*nul[12]+0.4330127018922193*fr[0]*nul[12]; 
  incr2[3] = (-3.031088913245535*fr[11]*nul[21])+2.5617376914899*fr[7]*nul[21]-1.984313483298443*fr[3]*nul[21]+1.14564392373896*fr[1]*nul[21]-2.561737691489901*fr[11]*nul[17]+2.165063509461097*fr[7]*nul[17]-1.677050983124842*fr[3]*nul[17]+0.9682458365518543*fr[1]*nul[17]-1.984313483298442*fr[11]*nul[14]+1.677050983124842*fr[7]*nul[14]-1.299038105676658*fr[3]*nul[14]+0.75*fr[1]*nul[14]-1.14564392373896*fr[11]*nul[12]+0.9682458365518543*fr[7]*nul[12]-0.75*fr[3]*nul[12]+0.4330127018922193*fr[1]*nul[12]; 
  incr2[5] = 11.7393568818739*fr[9]*nul[21]-9.921567416492215*fr[5]*nul[21]+7.685213074469699*fr[2]*nul[21]-4.437059837324712*fr[0]*nul[21]+9.921567416492215*fr[9]*nul[17]-8.385254915624213*fr[5]*nul[17]+6.495190528383289*fr[2]*nul[17]-3.75*fr[0]*nul[17]+7.685213074469706*fr[9]*nul[14]-6.495190528383286*fr[5]*nul[14]+5.031152949374527*fr[2]*nul[14]-2.904737509655563*fr[0]*nul[14]+4.437059837324713*fr[9]*nul[12]-3.75*fr[5]*nul[12]+2.904737509655563*fr[2]*nul[12]-1.677050983124842*fr[0]*nul[12]; 
  incr2[6] = (-1.984313483298443*fr[6]*nul[21])+1.14564392373896*fr[4]*nul[21]-1.677050983124842*fr[6]*nul[17]+0.9682458365518543*fr[4]*nul[17]-1.299038105676658*fr[6]*nul[14]+0.75*fr[4]*nul[14]-0.75*fr[6]*nul[12]+0.4330127018922194*fr[4]*nul[12]; 
  incr2[7] = 11.7393568818739*fr[11]*nul[21]-9.921567416492215*fr[7]*nul[21]+7.6852130744697*fr[3]*nul[21]-4.437059837324712*fr[1]*nul[21]+9.921567416492222*fr[11]*nul[17]-8.385254915624213*fr[7]*nul[17]+6.49519052838329*fr[3]*nul[17]-3.75*fr[1]*nul[17]+7.685213074469699*fr[11]*nul[14]-6.495190528383286*fr[7]*nul[14]+5.031152949374526*fr[3]*nul[14]-2.904737509655563*fr[1]*nul[14]+4.437059837324712*fr[11]*nul[12]-3.75*fr[7]*nul[12]+2.904737509655563*fr[3]*nul[12]-1.677050983124842*fr[1]*nul[12]; 
  incr2[9] = (-27.78038876617821*fr[9]*nul[21])+23.47871376374779*fr[5]*nul[21]-18.18653347947321*fr[2]*nul[21]+10.5*fr[0]*nul[21]-23.47871376374781*fr[9]*nul[17]+19.84313483298443*fr[5]*nul[17]-15.3704261489394*fr[2]*nul[17]+8.874119674649425*fr[0]*nul[17]-18.18653347947322*fr[9]*nul[14]+15.37042614893939*fr[5]*nul[14]-11.90588089979066*fr[2]*nul[14]+6.873863542433759*fr[0]*nul[14]-10.5*fr[9]*nul[12]+8.874119674649425*fr[5]*nul[12]-6.873863542433759*fr[2]*nul[12]+3.968626966596886*fr[0]*nul[12]; 
  incr2[10] = (-1.984313483298443*fr[10]*nul[21])+1.14564392373896*fr[8]*nul[21]-1.677050983124842*fr[10]*nul[17]+0.9682458365518541*fr[8]*nul[17]-1.299038105676658*fr[10]*nul[14]+0.7499999999999999*fr[8]*nul[14]-0.75*fr[10]*nul[12]+0.4330127018922193*fr[8]*nul[12]; 
  incr2[11] = (-27.78038876617821*fr[11]*nul[21])+23.47871376374779*fr[7]*nul[21]-18.18653347947321*fr[3]*nul[21]+10.5*fr[1]*nul[21]-23.47871376374781*fr[11]*nul[17]+19.84313483298443*fr[7]*nul[17]-15.3704261489394*fr[3]*nul[17]+8.874119674649425*fr[1]*nul[17]-18.18653347947322*fr[11]*nul[14]+15.37042614893939*fr[7]*nul[14]-11.90588089979066*fr[3]*nul[14]+6.87386354243376*fr[1]*nul[14]-10.5*fr[11]*nul[12]+8.874119674649425*fr[7]*nul[12]-6.87386354243376*fr[3]*nul[12]+3.968626966596886*fr[1]*nul[12]; 

  outr[2] += incr2[2]*rdxFr; 
  outr[3] += incr2[3]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 
  outr[9] += incr2[9]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 

  } else {

  incr2[2] = 3.031088913245536*fl[9]*nul[21]+2.561737691489899*fl[5]*nul[21]+1.984313483298443*fl[2]*nul[21]+1.14564392373896*fl[0]*nul[21]+2.561737691489901*fl[9]*nul[17]+2.165063509461096*fl[5]*nul[17]+1.677050983124842*fl[2]*nul[17]+0.9682458365518543*fl[0]*nul[17]+1.984313483298443*fl[9]*nul[14]+1.677050983124842*fl[5]*nul[14]+1.299038105676658*fl[2]*nul[14]+0.75*fl[0]*nul[14]+1.14564392373896*fl[9]*nul[12]+0.9682458365518541*fl[5]*nul[12]+0.75*fl[2]*nul[12]+0.4330127018922193*fl[0]*nul[12]; 
  incr2[3] = 3.031088913245535*fl[11]*nul[21]+2.5617376914899*fl[7]*nul[21]+1.984313483298443*fl[3]*nul[21]+1.14564392373896*fl[1]*nul[21]+2.561737691489901*fl[11]*nul[17]+2.165063509461097*fl[7]*nul[17]+1.677050983124842*fl[3]*nul[17]+0.9682458365518543*fl[1]*nul[17]+1.984313483298442*fl[11]*nul[14]+1.677050983124842*fl[7]*nul[14]+1.299038105676658*fl[3]*nul[14]+0.75*fl[1]*nul[14]+1.14564392373896*fl[11]*nul[12]+0.9682458365518543*fl[7]*nul[12]+0.75*fl[3]*nul[12]+0.4330127018922193*fl[1]*nul[12]; 
  incr2[5] = (-11.7393568818739*fl[9]*nul[21])-9.921567416492215*fl[5]*nul[21]-7.685213074469699*fl[2]*nul[21]-4.437059837324712*fl[0]*nul[21]-9.921567416492215*fl[9]*nul[17]-8.385254915624213*fl[5]*nul[17]-6.495190528383289*fl[2]*nul[17]-3.75*fl[0]*nul[17]-7.685213074469706*fl[9]*nul[14]-6.495190528383286*fl[5]*nul[14]-5.031152949374527*fl[2]*nul[14]-2.904737509655563*fl[0]*nul[14]-4.437059837324713*fl[9]*nul[12]-3.75*fl[5]*nul[12]-2.904737509655563*fl[2]*nul[12]-1.677050983124842*fl[0]*nul[12]; 
  incr2[6] = 1.984313483298443*fl[6]*nul[21]+1.14564392373896*fl[4]*nul[21]+1.677050983124842*fl[6]*nul[17]+0.9682458365518543*fl[4]*nul[17]+1.299038105676658*fl[6]*nul[14]+0.75*fl[4]*nul[14]+0.75*fl[6]*nul[12]+0.4330127018922194*fl[4]*nul[12]; 
  incr2[7] = (-11.7393568818739*fl[11]*nul[21])-9.921567416492215*fl[7]*nul[21]-7.6852130744697*fl[3]*nul[21]-4.437059837324712*fl[1]*nul[21]-9.921567416492222*fl[11]*nul[17]-8.385254915624213*fl[7]*nul[17]-6.49519052838329*fl[3]*nul[17]-3.75*fl[1]*nul[17]-7.685213074469699*fl[11]*nul[14]-6.495190528383286*fl[7]*nul[14]-5.031152949374526*fl[3]*nul[14]-2.904737509655563*fl[1]*nul[14]-4.437059837324712*fl[11]*nul[12]-3.75*fl[7]*nul[12]-2.904737509655563*fl[3]*nul[12]-1.677050983124842*fl[1]*nul[12]; 
  incr2[9] = 27.78038876617821*fl[9]*nul[21]+23.47871376374779*fl[5]*nul[21]+18.18653347947321*fl[2]*nul[21]+10.5*fl[0]*nul[21]+23.47871376374781*fl[9]*nul[17]+19.84313483298443*fl[5]*nul[17]+15.3704261489394*fl[2]*nul[17]+8.874119674649425*fl[0]*nul[17]+18.18653347947322*fl[9]*nul[14]+15.37042614893939*fl[5]*nul[14]+11.90588089979066*fl[2]*nul[14]+6.873863542433759*fl[0]*nul[14]+10.5*fl[9]*nul[12]+8.874119674649425*fl[5]*nul[12]+6.873863542433759*fl[2]*nul[12]+3.968626966596886*fl[0]*nul[12]; 
  incr2[10] = 1.984313483298443*fl[10]*nul[21]+1.14564392373896*fl[8]*nul[21]+1.677050983124842*fl[10]*nul[17]+0.9682458365518541*fl[8]*nul[17]+1.299038105676658*fl[10]*nul[14]+0.7499999999999999*fl[8]*nul[14]+0.75*fl[10]*nul[12]+0.4330127018922193*fl[8]*nul[12]; 
  incr2[11] = 27.78038876617821*fl[11]*nul[21]+23.47871376374779*fl[7]*nul[21]+18.18653347947321*fl[3]*nul[21]+10.5*fl[1]*nul[21]+23.47871376374781*fl[11]*nul[17]+19.84313483298443*fl[7]*nul[17]+15.3704261489394*fl[3]*nul[17]+8.874119674649425*fl[1]*nul[17]+18.18653347947322*fl[11]*nul[14]+15.37042614893939*fl[7]*nul[14]+11.90588089979066*fl[3]*nul[14]+6.87386354243376*fl[1]*nul[14]+10.5*fl[11]*nul[12]+8.874119674649425*fl[7]*nul[12]+6.87386354243376*fl[3]*nul[12]+3.968626966596886*fl[1]*nul[12]; 

  outl[2] += -1.0*incr2[2]*rdxFl; 
  outl[3] += -1.0*incr2[3]*rdxFl; 
  outl[5] += incr2[5]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[7] += incr2[7]*rdxFl; 
  outl[9] += -1.0*incr2[9]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[11] += -1.0*incr2[11]*rdxFl; 

  }

} 
