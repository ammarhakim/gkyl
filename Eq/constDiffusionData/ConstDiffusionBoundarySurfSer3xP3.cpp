#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf3xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[32]; 
  double incr2[32]; 

  if (edge < 0) {

  incr2[1] = (-2.29128784747792*fr[17])+1.936491673103708*fr[7]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[4] = (-2.29128784747792*fr[23])+1.936491673103709*fr[11]-1.5*fr[4]+0.8660254037844386*fr[2]; 
  incr2[5] = (-2.29128784747792*fr[25])+1.936491673103709*fr[13]-1.5*fr[5]+0.8660254037844386*fr[3]; 
  incr2[7] = 8.874119674649426*fr[17]-7.5*fr[7]+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 
  incr2[10] = (-2.29128784747792*fr[29])+1.936491673103708*fr[20]-1.5*fr[10]+0.8660254037844386*fr[6]; 
  incr2[11] = 8.874119674649425*fr[23]-7.5*fr[11]+5.809475019311126*fr[4]-3.354101966249684*fr[2]; 
  incr2[12] = 0.8660254037844387*fr[8]-1.5*fr[12]; 
  incr2[13] = 8.874119674649425*fr[25]-7.5*fr[13]+5.809475019311126*fr[5]-3.354101966249684*fr[3]; 
  incr2[15] = 0.8660254037844387*fr[9]-1.5*fr[15]; 
  incr2[17] = (-21.0*fr[17])+17.74823934929885*fr[7]-13.74772708486752*fr[1]+7.937253933193772*fr[0]; 
  incr2[20] = 8.874119674649426*fr[29]-7.5*fr[20]+5.809475019311125*fr[10]-3.354101966249685*fr[6]; 
  incr2[21] = 0.8660254037844387*fr[14]-1.5*fr[21]; 
  incr2[22] = 0.8660254037844387*fr[16]-1.5*fr[22]; 
  incr2[23] = (-21.0*fr[23])+17.74823934929885*fr[11]-13.74772708486752*fr[4]+7.937253933193771*fr[2]; 
  incr2[24] = 0.8660254037844386*fr[18]-1.5*fr[24]; 
  incr2[25] = (-21.0*fr[25])+17.74823934929885*fr[13]-13.74772708486752*fr[5]+7.937253933193771*fr[3]; 
  incr2[27] = 0.8660254037844386*fr[19]-1.5*fr[27]; 
  incr2[29] = (-21.0*fr[29])+17.74823934929885*fr[20]-13.74772708486752*fr[10]+7.937253933193772*fr[6]; 
  incr2[30] = 0.8660254037844386*fr[26]-1.5*fr[30]; 
  incr2[31] = 0.8660254037844386*fr[28]-1.5*fr[31]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 

  } else {

  incr2[1] = 2.29128784747792*fl[17]+1.936491673103708*fl[7]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[4] = 2.29128784747792*fl[23]+1.936491673103709*fl[11]+1.5*fl[4]+0.8660254037844386*fl[2]; 
  incr2[5] = 2.29128784747792*fl[25]+1.936491673103709*fl[13]+1.5*fl[5]+0.8660254037844386*fl[3]; 
  incr2[7] = (-8.874119674649426*fl[17])-7.5*fl[7]-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 
  incr2[10] = 2.29128784747792*fl[29]+1.936491673103708*fl[20]+1.5*fl[10]+0.8660254037844386*fl[6]; 
  incr2[11] = (-8.874119674649425*fl[23])-7.5*fl[11]-5.809475019311126*fl[4]-3.354101966249684*fl[2]; 
  incr2[12] = 1.5*fl[12]+0.8660254037844387*fl[8]; 
  incr2[13] = (-8.874119674649425*fl[25])-7.5*fl[13]-5.809475019311126*fl[5]-3.354101966249684*fl[3]; 
  incr2[15] = 1.5*fl[15]+0.8660254037844387*fl[9]; 
  incr2[17] = 21.0*fl[17]+17.74823934929885*fl[7]+13.74772708486752*fl[1]+7.937253933193772*fl[0]; 
  incr2[20] = (-8.874119674649426*fl[29])-7.5*fl[20]-5.809475019311125*fl[10]-3.354101966249685*fl[6]; 
  incr2[21] = 1.5*fl[21]+0.8660254037844387*fl[14]; 
  incr2[22] = 1.5*fl[22]+0.8660254037844387*fl[16]; 
  incr2[23] = 21.0*fl[23]+17.74823934929885*fl[11]+13.74772708486752*fl[4]+7.937253933193771*fl[2]; 
  incr2[24] = 1.5*fl[24]+0.8660254037844386*fl[18]; 
  incr2[25] = 21.0*fl[25]+17.74823934929885*fl[13]+13.74772708486752*fl[5]+7.937253933193771*fl[3]; 
  incr2[27] = 1.5*fl[27]+0.8660254037844386*fl[19]; 
  incr2[29] = 21.0*fl[29]+17.74823934929885*fl[20]+13.74772708486752*fl[10]+7.937253933193772*fl[6]; 
  incr2[30] = 1.5*fl[30]+0.8660254037844386*fl[26]; 
  incr2[31] = 1.5*fl[31]+0.8660254037844386*fl[28]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf3xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[32]; 
  double incr2[32]; 

  if (edge < 0) {

  incr2[2] = (-2.29128784747792*fr[18])+1.936491673103708*fr[8]-1.5*fr[2]+0.8660254037844386*fr[0]; 
  incr2[4] = (-2.29128784747792*fr[24])+1.936491673103709*fr[12]-1.5*fr[4]+0.8660254037844386*fr[1]; 
  incr2[6] = (-2.29128784747792*fr[26])+1.936491673103709*fr[14]-1.5*fr[6]+0.8660254037844386*fr[3]; 
  incr2[8] = 8.874119674649426*fr[18]-7.5*fr[8]+5.809475019311125*fr[2]-3.354101966249685*fr[0]; 
  incr2[10] = (-2.29128784747792*fr[30])+1.936491673103708*fr[21]-1.5*fr[10]+0.8660254037844386*fr[5]; 
  incr2[11] = 0.8660254037844387*fr[7]-1.5*fr[11]; 
  incr2[12] = 8.874119674649425*fr[24]-7.5*fr[12]+5.809475019311126*fr[4]-3.354101966249684*fr[1]; 
  incr2[14] = 8.874119674649425*fr[26]-7.5*fr[14]+5.809475019311126*fr[6]-3.354101966249684*fr[3]; 
  incr2[16] = 0.8660254037844387*fr[9]-1.5*fr[16]; 
  incr2[18] = (-21.0*fr[18])+17.74823934929885*fr[8]-13.74772708486752*fr[2]+7.937253933193772*fr[0]; 
  incr2[20] = 0.8660254037844387*fr[13]-1.5*fr[20]; 
  incr2[21] = 8.874119674649426*fr[30]-7.5*fr[21]+5.809475019311125*fr[10]-3.354101966249685*fr[5]; 
  incr2[22] = 0.8660254037844387*fr[15]-1.5*fr[22]; 
  incr2[23] = 0.8660254037844386*fr[17]-1.5*fr[23]; 
  incr2[24] = (-21.0*fr[24])+17.74823934929885*fr[12]-13.74772708486752*fr[4]+7.937253933193771*fr[1]; 
  incr2[26] = (-21.0*fr[26])+17.74823934929885*fr[14]-13.74772708486752*fr[6]+7.937253933193771*fr[3]; 
  incr2[28] = 0.8660254037844386*fr[19]-1.5*fr[28]; 
  incr2[29] = 0.8660254037844386*fr[25]-1.5*fr[29]; 
  incr2[30] = (-21.0*fr[30])+17.74823934929885*fr[21]-13.74772708486752*fr[10]+7.937253933193772*fr[5]; 
  incr2[31] = 0.8660254037844386*fr[27]-1.5*fr[31]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 

  } else {

  incr2[2] = 2.29128784747792*fl[18]+1.936491673103708*fl[8]+1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[4] = 2.29128784747792*fl[24]+1.936491673103709*fl[12]+1.5*fl[4]+0.8660254037844386*fl[1]; 
  incr2[6] = 2.29128784747792*fl[26]+1.936491673103709*fl[14]+1.5*fl[6]+0.8660254037844386*fl[3]; 
  incr2[8] = (-8.874119674649426*fl[18])-7.5*fl[8]-5.809475019311125*fl[2]-3.354101966249685*fl[0]; 
  incr2[10] = 2.29128784747792*fl[30]+1.936491673103708*fl[21]+1.5*fl[10]+0.8660254037844386*fl[5]; 
  incr2[11] = 1.5*fl[11]+0.8660254037844387*fl[7]; 
  incr2[12] = (-8.874119674649425*fl[24])-7.5*fl[12]-5.809475019311126*fl[4]-3.354101966249684*fl[1]; 
  incr2[14] = (-8.874119674649425*fl[26])-7.5*fl[14]-5.809475019311126*fl[6]-3.354101966249684*fl[3]; 
  incr2[16] = 1.5*fl[16]+0.8660254037844387*fl[9]; 
  incr2[18] = 21.0*fl[18]+17.74823934929885*fl[8]+13.74772708486752*fl[2]+7.937253933193772*fl[0]; 
  incr2[20] = 1.5*fl[20]+0.8660254037844387*fl[13]; 
  incr2[21] = (-8.874119674649426*fl[30])-7.5*fl[21]-5.809475019311125*fl[10]-3.354101966249685*fl[5]; 
  incr2[22] = 1.5*fl[22]+0.8660254037844387*fl[15]; 
  incr2[23] = 1.5*fl[23]+0.8660254037844386*fl[17]; 
  incr2[24] = 21.0*fl[24]+17.74823934929885*fl[12]+13.74772708486752*fl[4]+7.937253933193771*fl[1]; 
  incr2[26] = 21.0*fl[26]+17.74823934929885*fl[14]+13.74772708486752*fl[6]+7.937253933193771*fl[3]; 
  incr2[28] = 1.5*fl[28]+0.8660254037844386*fl[19]; 
  incr2[29] = 1.5*fl[29]+0.8660254037844386*fl[25]; 
  incr2[30] = 21.0*fl[30]+17.74823934929885*fl[21]+13.74772708486752*fl[10]+7.937253933193772*fl[5]; 
  incr2[31] = 1.5*fl[31]+0.8660254037844386*fl[27]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf3xSerP3_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[32]; 
  double incr2[32]; 

  if (edge < 0) {

  incr2[3] = (-2.29128784747792*fr[19])+1.936491673103708*fr[9]-1.5*fr[3]+0.8660254037844386*fr[0]; 
  incr2[5] = (-2.29128784747792*fr[27])+1.936491673103709*fr[15]-1.5*fr[5]+0.8660254037844386*fr[1]; 
  incr2[6] = (-2.29128784747792*fr[28])+1.936491673103709*fr[16]-1.5*fr[6]+0.8660254037844386*fr[2]; 
  incr2[9] = 8.874119674649426*fr[19]-7.5*fr[9]+5.809475019311125*fr[3]-3.354101966249685*fr[0]; 
  incr2[10] = (-2.29128784747792*fr[31])+1.936491673103708*fr[22]-1.5*fr[10]+0.8660254037844386*fr[4]; 
  incr2[13] = 0.8660254037844387*fr[7]-1.5*fr[13]; 
  incr2[14] = 0.8660254037844387*fr[8]-1.5*fr[14]; 
  incr2[15] = 8.874119674649425*fr[27]-7.5*fr[15]+5.809475019311126*fr[5]-3.354101966249684*fr[1]; 
  incr2[16] = 8.874119674649425*fr[28]-7.5*fr[16]+5.809475019311126*fr[6]-3.354101966249684*fr[2]; 
  incr2[19] = (-21.0*fr[19])+17.74823934929885*fr[9]-13.74772708486752*fr[3]+7.937253933193772*fr[0]; 
  incr2[20] = 0.8660254037844387*fr[11]-1.5*fr[20]; 
  incr2[21] = 0.8660254037844387*fr[12]-1.5*fr[21]; 
  incr2[22] = 8.874119674649426*fr[31]-7.5*fr[22]+5.809475019311125*fr[10]-3.354101966249685*fr[4]; 
  incr2[25] = 0.8660254037844386*fr[17]-1.5*fr[25]; 
  incr2[26] = 0.8660254037844386*fr[18]-1.5*fr[26]; 
  incr2[27] = (-21.0*fr[27])+17.74823934929885*fr[15]-13.74772708486752*fr[5]+7.937253933193771*fr[1]; 
  incr2[28] = (-21.0*fr[28])+17.74823934929885*fr[16]-13.74772708486752*fr[6]+7.937253933193771*fr[2]; 
  incr2[29] = 0.8660254037844386*fr[23]-1.5*fr[29]; 
  incr2[30] = 0.8660254037844386*fr[24]-1.5*fr[30]; 
  incr2[31] = (-21.0*fr[31])+17.74823934929885*fr[22]-13.74772708486752*fr[10]+7.937253933193772*fr[4]; 

  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 

  } else {

  incr2[3] = 2.29128784747792*fl[19]+1.936491673103708*fl[9]+1.5*fl[3]+0.8660254037844386*fl[0]; 
  incr2[5] = 2.29128784747792*fl[27]+1.936491673103709*fl[15]+1.5*fl[5]+0.8660254037844386*fl[1]; 
  incr2[6] = 2.29128784747792*fl[28]+1.936491673103709*fl[16]+1.5*fl[6]+0.8660254037844386*fl[2]; 
  incr2[9] = (-8.874119674649426*fl[19])-7.5*fl[9]-5.809475019311125*fl[3]-3.354101966249685*fl[0]; 
  incr2[10] = 2.29128784747792*fl[31]+1.936491673103708*fl[22]+1.5*fl[10]+0.8660254037844386*fl[4]; 
  incr2[13] = 1.5*fl[13]+0.8660254037844387*fl[7]; 
  incr2[14] = 1.5*fl[14]+0.8660254037844387*fl[8]; 
  incr2[15] = (-8.874119674649425*fl[27])-7.5*fl[15]-5.809475019311126*fl[5]-3.354101966249684*fl[1]; 
  incr2[16] = (-8.874119674649425*fl[28])-7.5*fl[16]-5.809475019311126*fl[6]-3.354101966249684*fl[2]; 
  incr2[19] = 21.0*fl[19]+17.74823934929885*fl[9]+13.74772708486752*fl[3]+7.937253933193772*fl[0]; 
  incr2[20] = 1.5*fl[20]+0.8660254037844387*fl[11]; 
  incr2[21] = 1.5*fl[21]+0.8660254037844387*fl[12]; 
  incr2[22] = (-8.874119674649426*fl[31])-7.5*fl[22]-5.809475019311125*fl[10]-3.354101966249685*fl[4]; 
  incr2[25] = 1.5*fl[25]+0.8660254037844386*fl[17]; 
  incr2[26] = 1.5*fl[26]+0.8660254037844386*fl[18]; 
  incr2[27] = 21.0*fl[27]+17.74823934929885*fl[15]+13.74772708486752*fl[5]+7.937253933193771*fl[1]; 
  incr2[28] = 21.0*fl[28]+17.74823934929885*fl[16]+13.74772708486752*fl[6]+7.937253933193771*fl[2]; 
  incr2[29] = 1.5*fl[29]+0.8660254037844386*fl[23]; 
  incr2[30] = 1.5*fl[30]+0.8660254037844386*fl[24]; 
  incr2[31] = 21.0*fl[31]+17.74823934929885*fl[22]+13.74772708486752*fl[10]+7.937253933193772*fl[4]; 

  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (edge < 0) {


  incr2[1] = 11.78376607274358*fr[17]-11.61895003862225*fr[7]+4.5*fr[1]; 
  incr2[4] = 11.78376607274359*fr[23]-11.61895003862225*fr[11]+4.5*fr[4]; 
  incr2[5] = 11.78376607274359*fr[25]-11.61895003862225*fr[13]+4.5*fr[5]; 
  incr2[7] = (-45.6383297553399*fr[17])+45.0*fr[7]-17.42842505793337*fr[1]; 
  incr2[10] = 11.78376607274358*fr[29]-11.61895003862225*fr[20]+4.5*fr[10]; 
  incr2[11] = (-45.6383297553399*fr[23])+45.0*fr[11]-17.42842505793338*fr[4]; 
  incr2[12] = 4.5*fr[12]; 
  incr2[13] = (-45.6383297553399*fr[25])+45.0*fr[13]-17.42842505793338*fr[5]; 
  incr2[15] = 4.5*fr[15]; 
  incr2[17] = 108.0*fr[17]-106.4894360957931*fr[7]+41.24318125460255*fr[1]; 
  incr2[20] = (-45.6383297553399*fr[29])+45.0*fr[20]-17.42842505793337*fr[10]; 
  incr2[21] = 4.5*fr[21]; 
  incr2[22] = 4.5*fr[22]; 
  incr2[23] = 108.0*fr[23]-106.4894360957931*fr[11]+41.24318125460256*fr[4]; 
  incr2[24] = 4.5*fr[24]; 
  incr2[25] = 108.0*fr[25]-106.4894360957931*fr[13]+41.24318125460256*fr[5]; 
  incr2[27] = 4.5*fr[27]; 
  incr2[29] = 108.0*fr[29]-106.4894360957931*fr[20]+41.24318125460255*fr[10]; 
  incr2[30] = 4.5*fr[30]; 
  incr2[31] = 4.5*fr[31]; 


  incr4[17] = (-16.2*fr[17])+28.39718295887816*fr[7]-30.24499958670854*fr[1]+19.84313483298443*fr[0]; 
  incr4[23] = (-16.2*fr[23])+28.39718295887816*fr[11]-30.24499958670854*fr[4]+19.84313483298443*fr[2]; 
  incr4[25] = (-16.2*fr[25])+28.39718295887816*fr[13]-30.24499958670854*fr[5]+19.84313483298443*fr[3]; 
  incr4[29] = (-16.2*fr[29])+28.39718295887816*fr[20]-30.24499958670854*fr[10]+19.84313483298443*fr[6]; 

  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[17] += (-1.0*incr4[17]*rdxFnur)-1.0*incr2[17]*rdxFnur; 
  outr[20] += -1.0*incr2[20]*rdxFnur; 
  outr[21] += -1.0*incr2[21]*rdxFnur; 
  outr[22] += -1.0*incr2[22]*rdxFnur; 
  outr[23] += (-1.0*incr4[23]*rdxFnur)-1.0*incr2[23]*rdxFnur; 
  outr[24] += -1.0*incr2[24]*rdxFnur; 
  outr[25] += (-1.0*incr4[25]*rdxFnur)-1.0*incr2[25]*rdxFnur; 
  outr[27] += -1.0*incr2[27]*rdxFnur; 
  outr[29] += (-1.0*incr4[29]*rdxFnur)-1.0*incr2[29]*rdxFnur; 
  outr[30] += -1.0*incr2[30]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 

  } else {


  incr2[1] = 11.78376607274358*fl[17]+11.61895003862225*fl[7]+4.5*fl[1]; 
  incr2[4] = 11.78376607274359*fl[23]+11.61895003862225*fl[11]+4.5*fl[4]; 
  incr2[5] = 11.78376607274359*fl[25]+11.61895003862225*fl[13]+4.5*fl[5]; 
  incr2[7] = 45.6383297553399*fl[17]+45.0*fl[7]+17.42842505793337*fl[1]; 
  incr2[10] = 11.78376607274358*fl[29]+11.61895003862225*fl[20]+4.5*fl[10]; 
  incr2[11] = 45.6383297553399*fl[23]+45.0*fl[11]+17.42842505793338*fl[4]; 
  incr2[12] = 4.5*fl[12]; 
  incr2[13] = 45.6383297553399*fl[25]+45.0*fl[13]+17.42842505793338*fl[5]; 
  incr2[15] = 4.5*fl[15]; 
  incr2[17] = 108.0*fl[17]+106.4894360957931*fl[7]+41.24318125460255*fl[1]; 
  incr2[20] = 45.6383297553399*fl[29]+45.0*fl[20]+17.42842505793337*fl[10]; 
  incr2[21] = 4.5*fl[21]; 
  incr2[22] = 4.5*fl[22]; 
  incr2[23] = 108.0*fl[23]+106.4894360957931*fl[11]+41.24318125460256*fl[4]; 
  incr2[24] = 4.5*fl[24]; 
  incr2[25] = 108.0*fl[25]+106.4894360957931*fl[13]+41.24318125460256*fl[5]; 
  incr2[27] = 4.5*fl[27]; 
  incr2[29] = 108.0*fl[29]+106.4894360957931*fl[20]+41.24318125460255*fl[10]; 
  incr2[30] = 4.5*fl[30]; 
  incr2[31] = 4.5*fl[31]; 


  incr4[17] = (-16.2*fl[17])-28.39718295887816*fl[7]-30.24499958670854*fl[1]-19.84313483298443*fl[0]; 
  incr4[23] = (-16.2*fl[23])-28.39718295887816*fl[11]-30.24499958670854*fl[4]-19.84313483298443*fl[2]; 
  incr4[25] = (-16.2*fl[25])-28.39718295887816*fl[13]-30.24499958670854*fl[5]-19.84313483298443*fl[3]; 
  incr4[29] = (-16.2*fl[29])-28.39718295887816*fl[20]-30.24499958670854*fl[10]-19.84313483298443*fl[6]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[17] += (-1.0*incr4[17]*rdxFnul)-1.0*incr2[17]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[23] += (-1.0*incr4[23]*rdxFnul)-1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[25] += (-1.0*incr4[25]*rdxFnul)-1.0*incr2[25]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[29] += (-1.0*incr4[29]*rdxFnul)-1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (edge < 0) {


  incr2[2] = 11.78376607274358*fr[18]-11.61895003862225*fr[8]+4.5*fr[2]; 
  incr2[4] = 11.78376607274359*fr[24]-11.61895003862225*fr[12]+4.5*fr[4]; 
  incr2[6] = 11.78376607274359*fr[26]-11.61895003862225*fr[14]+4.5*fr[6]; 
  incr2[8] = (-45.6383297553399*fr[18])+45.0*fr[8]-17.42842505793337*fr[2]; 
  incr2[10] = 11.78376607274358*fr[30]-11.61895003862225*fr[21]+4.5*fr[10]; 
  incr2[11] = 4.5*fr[11]; 
  incr2[12] = (-45.6383297553399*fr[24])+45.0*fr[12]-17.42842505793338*fr[4]; 
  incr2[14] = (-45.6383297553399*fr[26])+45.0*fr[14]-17.42842505793338*fr[6]; 
  incr2[16] = 4.5*fr[16]; 
  incr2[18] = 108.0*fr[18]-106.4894360957931*fr[8]+41.24318125460255*fr[2]; 
  incr2[20] = 4.5*fr[20]; 
  incr2[21] = (-45.6383297553399*fr[30])+45.0*fr[21]-17.42842505793337*fr[10]; 
  incr2[22] = 4.5*fr[22]; 
  incr2[23] = 4.5*fr[23]; 
  incr2[24] = 108.0*fr[24]-106.4894360957931*fr[12]+41.24318125460256*fr[4]; 
  incr2[26] = 108.0*fr[26]-106.4894360957931*fr[14]+41.24318125460256*fr[6]; 
  incr2[28] = 4.5*fr[28]; 
  incr2[29] = 4.5*fr[29]; 
  incr2[30] = 108.0*fr[30]-106.4894360957931*fr[21]+41.24318125460255*fr[10]; 
  incr2[31] = 4.5*fr[31]; 


  incr4[18] = (-16.2*fr[18])+28.39718295887816*fr[8]-30.24499958670854*fr[2]+19.84313483298443*fr[0]; 
  incr4[24] = (-16.2*fr[24])+28.39718295887816*fr[12]-30.24499958670854*fr[4]+19.84313483298443*fr[1]; 
  incr4[26] = (-16.2*fr[26])+28.39718295887816*fr[14]-30.24499958670854*fr[6]+19.84313483298443*fr[3]; 
  incr4[30] = (-16.2*fr[30])+28.39718295887816*fr[21]-30.24499958670854*fr[10]+19.84313483298443*fr[5]; 

  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[8] += -1.0*incr2[8]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[18] += (-1.0*incr4[18]*rdxFnur)-1.0*incr2[18]*rdxFnur; 
  outr[20] += -1.0*incr2[20]*rdxFnur; 
  outr[21] += -1.0*incr2[21]*rdxFnur; 
  outr[22] += -1.0*incr2[22]*rdxFnur; 
  outr[23] += -1.0*incr2[23]*rdxFnur; 
  outr[24] += (-1.0*incr4[24]*rdxFnur)-1.0*incr2[24]*rdxFnur; 
  outr[26] += (-1.0*incr4[26]*rdxFnur)-1.0*incr2[26]*rdxFnur; 
  outr[28] += -1.0*incr2[28]*rdxFnur; 
  outr[29] += -1.0*incr2[29]*rdxFnur; 
  outr[30] += (-1.0*incr4[30]*rdxFnur)-1.0*incr2[30]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 

  } else {


  incr2[2] = 11.78376607274358*fl[18]+11.61895003862225*fl[8]+4.5*fl[2]; 
  incr2[4] = 11.78376607274359*fl[24]+11.61895003862225*fl[12]+4.5*fl[4]; 
  incr2[6] = 11.78376607274359*fl[26]+11.61895003862225*fl[14]+4.5*fl[6]; 
  incr2[8] = 45.6383297553399*fl[18]+45.0*fl[8]+17.42842505793337*fl[2]; 
  incr2[10] = 11.78376607274358*fl[30]+11.61895003862225*fl[21]+4.5*fl[10]; 
  incr2[11] = 4.5*fl[11]; 
  incr2[12] = 45.6383297553399*fl[24]+45.0*fl[12]+17.42842505793338*fl[4]; 
  incr2[14] = 45.6383297553399*fl[26]+45.0*fl[14]+17.42842505793338*fl[6]; 
  incr2[16] = 4.5*fl[16]; 
  incr2[18] = 108.0*fl[18]+106.4894360957931*fl[8]+41.24318125460255*fl[2]; 
  incr2[20] = 4.5*fl[20]; 
  incr2[21] = 45.6383297553399*fl[30]+45.0*fl[21]+17.42842505793337*fl[10]; 
  incr2[22] = 4.5*fl[22]; 
  incr2[23] = 4.5*fl[23]; 
  incr2[24] = 108.0*fl[24]+106.4894360957931*fl[12]+41.24318125460256*fl[4]; 
  incr2[26] = 108.0*fl[26]+106.4894360957931*fl[14]+41.24318125460256*fl[6]; 
  incr2[28] = 4.5*fl[28]; 
  incr2[29] = 4.5*fl[29]; 
  incr2[30] = 108.0*fl[30]+106.4894360957931*fl[21]+41.24318125460255*fl[10]; 
  incr2[31] = 4.5*fl[31]; 


  incr4[18] = (-16.2*fl[18])-28.39718295887816*fl[8]-30.24499958670854*fl[2]-19.84313483298443*fl[0]; 
  incr4[24] = (-16.2*fl[24])-28.39718295887816*fl[12]-30.24499958670854*fl[4]-19.84313483298443*fl[1]; 
  incr4[26] = (-16.2*fl[26])-28.39718295887816*fl[14]-30.24499958670854*fl[6]-19.84313483298443*fl[3]; 
  incr4[30] = (-16.2*fl[30])-28.39718295887816*fl[21]-30.24499958670854*fl[10]-19.84313483298443*fl[5]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[18] += (-1.0*incr4[18]*rdxFnul)-1.0*incr2[18]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += (-1.0*incr4[24]*rdxFnul)-1.0*incr2[24]*rdxFnul; 
  outl[26] += (-1.0*incr4[26]*rdxFnul)-1.0*incr2[26]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += (-1.0*incr4[30]*rdxFnul)-1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xSerP3_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (edge < 0) {


  incr2[3] = 11.78376607274358*fr[19]-11.61895003862225*fr[9]+4.5*fr[3]; 
  incr2[5] = 11.78376607274359*fr[27]-11.61895003862225*fr[15]+4.5*fr[5]; 
  incr2[6] = 11.78376607274359*fr[28]-11.61895003862225*fr[16]+4.5*fr[6]; 
  incr2[9] = (-45.6383297553399*fr[19])+45.0*fr[9]-17.42842505793337*fr[3]; 
  incr2[10] = 11.78376607274358*fr[31]-11.61895003862225*fr[22]+4.5*fr[10]; 
  incr2[13] = 4.5*fr[13]; 
  incr2[14] = 4.5*fr[14]; 
  incr2[15] = (-45.6383297553399*fr[27])+45.0*fr[15]-17.42842505793338*fr[5]; 
  incr2[16] = (-45.6383297553399*fr[28])+45.0*fr[16]-17.42842505793338*fr[6]; 
  incr2[19] = 108.0*fr[19]-106.4894360957931*fr[9]+41.24318125460255*fr[3]; 
  incr2[20] = 4.5*fr[20]; 
  incr2[21] = 4.5*fr[21]; 
  incr2[22] = (-45.6383297553399*fr[31])+45.0*fr[22]-17.42842505793337*fr[10]; 
  incr2[25] = 4.5*fr[25]; 
  incr2[26] = 4.5*fr[26]; 
  incr2[27] = 108.0*fr[27]-106.4894360957931*fr[15]+41.24318125460256*fr[5]; 
  incr2[28] = 108.0*fr[28]-106.4894360957931*fr[16]+41.24318125460256*fr[6]; 
  incr2[29] = 4.5*fr[29]; 
  incr2[30] = 4.5*fr[30]; 
  incr2[31] = 108.0*fr[31]-106.4894360957931*fr[22]+41.24318125460255*fr[10]; 


  incr4[19] = (-16.2*fr[19])+28.39718295887816*fr[9]-30.24499958670854*fr[3]+19.84313483298443*fr[0]; 
  incr4[27] = (-16.2*fr[27])+28.39718295887816*fr[15]-30.24499958670854*fr[5]+19.84313483298443*fr[1]; 
  incr4[28] = (-16.2*fr[28])+28.39718295887816*fr[16]-30.24499958670854*fr[6]+19.84313483298443*fr[2]; 
  incr4[31] = (-16.2*fr[31])+28.39718295887816*fr[22]-30.24499958670854*fr[10]+19.84313483298443*fr[4]; 

  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[9] += -1.0*incr2[9]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[19] += (-1.0*incr4[19]*rdxFnur)-1.0*incr2[19]*rdxFnur; 
  outr[20] += -1.0*incr2[20]*rdxFnur; 
  outr[21] += -1.0*incr2[21]*rdxFnur; 
  outr[22] += -1.0*incr2[22]*rdxFnur; 
  outr[25] += -1.0*incr2[25]*rdxFnur; 
  outr[26] += -1.0*incr2[26]*rdxFnur; 
  outr[27] += (-1.0*incr4[27]*rdxFnur)-1.0*incr2[27]*rdxFnur; 
  outr[28] += (-1.0*incr4[28]*rdxFnur)-1.0*incr2[28]*rdxFnur; 
  outr[29] += -1.0*incr2[29]*rdxFnur; 
  outr[30] += -1.0*incr2[30]*rdxFnur; 
  outr[31] += (-1.0*incr4[31]*rdxFnur)-1.0*incr2[31]*rdxFnur; 

  } else {


  incr2[3] = 11.78376607274358*fl[19]+11.61895003862225*fl[9]+4.5*fl[3]; 
  incr2[5] = 11.78376607274359*fl[27]+11.61895003862225*fl[15]+4.5*fl[5]; 
  incr2[6] = 11.78376607274359*fl[28]+11.61895003862225*fl[16]+4.5*fl[6]; 
  incr2[9] = 45.6383297553399*fl[19]+45.0*fl[9]+17.42842505793337*fl[3]; 
  incr2[10] = 11.78376607274358*fl[31]+11.61895003862225*fl[22]+4.5*fl[10]; 
  incr2[13] = 4.5*fl[13]; 
  incr2[14] = 4.5*fl[14]; 
  incr2[15] = 45.6383297553399*fl[27]+45.0*fl[15]+17.42842505793338*fl[5]; 
  incr2[16] = 45.6383297553399*fl[28]+45.0*fl[16]+17.42842505793338*fl[6]; 
  incr2[19] = 108.0*fl[19]+106.4894360957931*fl[9]+41.24318125460255*fl[3]; 
  incr2[20] = 4.5*fl[20]; 
  incr2[21] = 4.5*fl[21]; 
  incr2[22] = 45.6383297553399*fl[31]+45.0*fl[22]+17.42842505793337*fl[10]; 
  incr2[25] = 4.5*fl[25]; 
  incr2[26] = 4.5*fl[26]; 
  incr2[27] = 108.0*fl[27]+106.4894360957931*fl[15]+41.24318125460256*fl[5]; 
  incr2[28] = 108.0*fl[28]+106.4894360957931*fl[16]+41.24318125460256*fl[6]; 
  incr2[29] = 4.5*fl[29]; 
  incr2[30] = 4.5*fl[30]; 
  incr2[31] = 108.0*fl[31]+106.4894360957931*fl[22]+41.24318125460255*fl[10]; 


  incr4[19] = (-16.2*fl[19])-28.39718295887816*fl[9]-30.24499958670854*fl[3]-19.84313483298443*fl[0]; 
  incr4[27] = (-16.2*fl[27])-28.39718295887816*fl[15]-30.24499958670854*fl[5]-19.84313483298443*fl[1]; 
  incr4[28] = (-16.2*fl[28])-28.39718295887816*fl[16]-30.24499958670854*fl[6]-19.84313483298443*fl[2]; 
  incr4[31] = (-16.2*fl[31])-28.39718295887816*fl[22]-30.24499958670854*fl[10]-19.84313483298443*fl[4]; 

  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[19] += (-1.0*incr4[19]*rdxFnul)-1.0*incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[27] += (-1.0*incr4[27]*rdxFnul)-1.0*incr2[27]*rdxFnul; 
  outl[28] += (-1.0*incr4[28]*rdxFnul)-1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += (-1.0*incr4[31]*rdxFnul)-1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf3xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 
  double incr5[32]; 
  double incr6[32]; 

  if (edge < 0) {


  incr2[1] = (-103.1079531365064*fr[17])+76.2493596284585*fr[7]-19.6875*fr[1]; 
  incr2[4] = (-103.1079531365064*fr[23])+76.24935962845854*fr[11]-19.6875*fr[4]; 
  incr2[5] = (-103.1079531365064*fr[25])+76.24935962845854*fr[13]-19.6875*fr[5]; 
  incr2[7] = 399.3353853592242*fr[17]-295.3125*fr[7]+76.2493596284585*fr[1]; 
  incr2[10] = (-103.1079531365064*fr[29])+76.2493596284585*fr[20]-19.6875*fr[10]; 
  incr2[11] = 399.3353853592241*fr[23]-295.3125*fr[11]+76.24935962845854*fr[4]; 
  incr2[12] = -19.6875*fr[12]; 
  incr2[13] = 399.3353853592241*fr[25]-295.3125*fr[13]+76.24935962845854*fr[5]; 
  incr2[15] = -19.6875*fr[15]; 
  incr2[17] = (-945.0*fr[17])+698.8369243786424*fr[7]-180.4389179888861*fr[1]; 
  incr2[20] = 399.3353853592242*fr[29]-295.3125*fr[20]+76.2493596284585*fr[10]; 
  incr2[21] = -19.6875*fr[21]; 
  incr2[22] = -19.6875*fr[22]; 
  incr2[23] = (-945.0*fr[23])+698.8369243786422*fr[11]-180.4389179888862*fr[4]; 
  incr2[24] = -19.6875*fr[24]; 
  incr2[25] = (-945.0*fr[25])+698.8369243786422*fr[13]-180.4389179888862*fr[5]; 
  incr2[27] = -19.6875*fr[27]; 
  incr2[29] = (-945.0*fr[29])+698.8369243786424*fr[20]-180.4389179888861*fr[10]; 
  incr2[30] = -19.6875*fr[30]; 
  incr2[31] = -19.6875*fr[31]; 


  incr4[17] = 225.0*fr[17]-241.2651286545313*fr[7]+96.66370606547471*fr[1]; 
  incr4[23] = 225.0*fr[23]-241.2651286545313*fr[11]+96.66370606547474*fr[4]; 
  incr4[25] = 225.0*fr[25]-241.2651286545313*fr[13]+96.66370606547474*fr[5]; 
  incr4[29] = 225.0*fr[29]-241.2651286545313*fr[20]+96.66370606547471*fr[10]; 



  outr[1] += incr2[1]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[17] += incr4[17]*rdxFnur+incr2[17]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr4[23]*rdxFnur+incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[25] += incr4[25]*rdxFnur+incr2[25]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[29] += incr4[29]*rdxFnur+incr2[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 

  } else {


  incr2[1] = (-103.1079531365064*fl[17])-76.2493596284585*fl[7]-19.6875*fl[1]; 
  incr2[4] = (-103.1079531365064*fl[23])-76.24935962845854*fl[11]-19.6875*fl[4]; 
  incr2[5] = (-103.1079531365064*fl[25])-76.24935962845854*fl[13]-19.6875*fl[5]; 
  incr2[7] = (-399.3353853592242*fl[17])-295.3125*fl[7]-76.2493596284585*fl[1]; 
  incr2[10] = (-103.1079531365064*fl[29])-76.2493596284585*fl[20]-19.6875*fl[10]; 
  incr2[11] = (-399.3353853592241*fl[23])-295.3125*fl[11]-76.24935962845854*fl[4]; 
  incr2[12] = -19.6875*fl[12]; 
  incr2[13] = (-399.3353853592241*fl[25])-295.3125*fl[13]-76.24935962845854*fl[5]; 
  incr2[15] = -19.6875*fl[15]; 
  incr2[17] = (-945.0*fl[17])-698.8369243786424*fl[7]-180.4389179888861*fl[1]; 
  incr2[20] = (-399.3353853592242*fl[29])-295.3125*fl[20]-76.2493596284585*fl[10]; 
  incr2[21] = -19.6875*fl[21]; 
  incr2[22] = -19.6875*fl[22]; 
  incr2[23] = (-945.0*fl[23])-698.8369243786422*fl[11]-180.4389179888862*fl[4]; 
  incr2[24] = -19.6875*fl[24]; 
  incr2[25] = (-945.0*fl[25])-698.8369243786422*fl[13]-180.4389179888862*fl[5]; 
  incr2[27] = -19.6875*fl[27]; 
  incr2[29] = (-945.0*fl[29])-698.8369243786424*fl[20]-180.4389179888861*fl[10]; 
  incr2[30] = -19.6875*fl[30]; 
  incr2[31] = -19.6875*fl[31]; 


  incr4[17] = 225.0*fl[17]+241.2651286545313*fl[7]+96.66370606547471*fl[1]; 
  incr4[23] = 225.0*fl[23]+241.2651286545313*fl[11]+96.66370606547474*fl[4]; 
  incr4[25] = 225.0*fl[25]+241.2651286545313*fl[13]+96.66370606547474*fl[5]; 
  incr4[29] = 225.0*fl[29]+241.2651286545313*fl[20]+96.66370606547471*fl[10]; 



  outl[1] += incr2[1]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[17] += incr4[17]*rdxFnul+incr2[17]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[23] += incr4[23]*rdxFnul+incr2[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul; 
  outl[25] += incr4[25]*rdxFnul+incr2[25]*rdxFnul; 
  outl[27] += incr2[27]*rdxFnul; 
  outl[29] += incr4[29]*rdxFnul+incr2[29]*rdxFnul; 
  outl[30] += incr2[30]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf3xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 
  double incr5[32]; 
  double incr6[32]; 

  if (edge < 0) {


  incr2[2] = (-103.1079531365064*fr[18])+76.2493596284585*fr[8]-19.6875*fr[2]; 
  incr2[4] = (-103.1079531365064*fr[24])+76.24935962845854*fr[12]-19.6875*fr[4]; 
  incr2[6] = (-103.1079531365064*fr[26])+76.24935962845854*fr[14]-19.6875*fr[6]; 
  incr2[8] = 399.3353853592242*fr[18]-295.3125*fr[8]+76.2493596284585*fr[2]; 
  incr2[10] = (-103.1079531365064*fr[30])+76.2493596284585*fr[21]-19.6875*fr[10]; 
  incr2[11] = -19.6875*fr[11]; 
  incr2[12] = 399.3353853592241*fr[24]-295.3125*fr[12]+76.24935962845854*fr[4]; 
  incr2[14] = 399.3353853592241*fr[26]-295.3125*fr[14]+76.24935962845854*fr[6]; 
  incr2[16] = -19.6875*fr[16]; 
  incr2[18] = (-945.0*fr[18])+698.8369243786424*fr[8]-180.4389179888861*fr[2]; 
  incr2[20] = -19.6875*fr[20]; 
  incr2[21] = 399.3353853592242*fr[30]-295.3125*fr[21]+76.2493596284585*fr[10]; 
  incr2[22] = -19.6875*fr[22]; 
  incr2[23] = -19.6875*fr[23]; 
  incr2[24] = (-945.0*fr[24])+698.8369243786422*fr[12]-180.4389179888862*fr[4]; 
  incr2[26] = (-945.0*fr[26])+698.8369243786422*fr[14]-180.4389179888862*fr[6]; 
  incr2[28] = -19.6875*fr[28]; 
  incr2[29] = -19.6875*fr[29]; 
  incr2[30] = (-945.0*fr[30])+698.8369243786424*fr[21]-180.4389179888861*fr[10]; 
  incr2[31] = -19.6875*fr[31]; 


  incr4[18] = 225.0*fr[18]-241.2651286545313*fr[8]+96.66370606547471*fr[2]; 
  incr4[24] = 225.0*fr[24]-241.2651286545313*fr[12]+96.66370606547474*fr[4]; 
  incr4[26] = 225.0*fr[26]-241.2651286545313*fr[14]+96.66370606547474*fr[6]; 
  incr4[30] = 225.0*fr[30]-241.2651286545313*fr[21]+96.66370606547471*fr[10]; 



  outr[2] += incr2[2]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[18] += incr4[18]*rdxFnur+incr2[18]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr4[24]*rdxFnur+incr2[24]*rdxFnur; 
  outr[26] += incr4[26]*rdxFnur+incr2[26]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[30] += incr4[30]*rdxFnur+incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 

  } else {


  incr2[2] = (-103.1079531365064*fl[18])-76.2493596284585*fl[8]-19.6875*fl[2]; 
  incr2[4] = (-103.1079531365064*fl[24])-76.24935962845854*fl[12]-19.6875*fl[4]; 
  incr2[6] = (-103.1079531365064*fl[26])-76.24935962845854*fl[14]-19.6875*fl[6]; 
  incr2[8] = (-399.3353853592242*fl[18])-295.3125*fl[8]-76.2493596284585*fl[2]; 
  incr2[10] = (-103.1079531365064*fl[30])-76.2493596284585*fl[21]-19.6875*fl[10]; 
  incr2[11] = -19.6875*fl[11]; 
  incr2[12] = (-399.3353853592241*fl[24])-295.3125*fl[12]-76.24935962845854*fl[4]; 
  incr2[14] = (-399.3353853592241*fl[26])-295.3125*fl[14]-76.24935962845854*fl[6]; 
  incr2[16] = -19.6875*fl[16]; 
  incr2[18] = (-945.0*fl[18])-698.8369243786424*fl[8]-180.4389179888861*fl[2]; 
  incr2[20] = -19.6875*fl[20]; 
  incr2[21] = (-399.3353853592242*fl[30])-295.3125*fl[21]-76.2493596284585*fl[10]; 
  incr2[22] = -19.6875*fl[22]; 
  incr2[23] = -19.6875*fl[23]; 
  incr2[24] = (-945.0*fl[24])-698.8369243786422*fl[12]-180.4389179888862*fl[4]; 
  incr2[26] = (-945.0*fl[26])-698.8369243786422*fl[14]-180.4389179888862*fl[6]; 
  incr2[28] = -19.6875*fl[28]; 
  incr2[29] = -19.6875*fl[29]; 
  incr2[30] = (-945.0*fl[30])-698.8369243786424*fl[21]-180.4389179888861*fl[10]; 
  incr2[31] = -19.6875*fl[31]; 


  incr4[18] = 225.0*fl[18]+241.2651286545313*fl[8]+96.66370606547471*fl[2]; 
  incr4[24] = 225.0*fl[24]+241.2651286545313*fl[12]+96.66370606547474*fl[4]; 
  incr4[26] = 225.0*fl[26]+241.2651286545313*fl[14]+96.66370606547474*fl[6]; 
  incr4[30] = 225.0*fl[30]+241.2651286545313*fl[21]+96.66370606547471*fl[10]; 



  outl[2] += incr2[2]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[18] += incr4[18]*rdxFnul+incr2[18]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul; 
  outl[24] += incr4[24]*rdxFnul+incr2[24]*rdxFnul; 
  outl[26] += incr4[26]*rdxFnul+incr2[26]*rdxFnul; 
  outl[28] += incr2[28]*rdxFnul; 
  outl[29] += incr2[29]*rdxFnul; 
  outl[30] += incr4[30]*rdxFnul+incr2[30]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf3xSerP3_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 64.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 
  double incr5[32]; 
  double incr6[32]; 

  if (edge < 0) {


  incr2[3] = (-103.1079531365064*fr[19])+76.2493596284585*fr[9]-19.6875*fr[3]; 
  incr2[5] = (-103.1079531365064*fr[27])+76.24935962845854*fr[15]-19.6875*fr[5]; 
  incr2[6] = (-103.1079531365064*fr[28])+76.24935962845854*fr[16]-19.6875*fr[6]; 
  incr2[9] = 399.3353853592242*fr[19]-295.3125*fr[9]+76.2493596284585*fr[3]; 
  incr2[10] = (-103.1079531365064*fr[31])+76.2493596284585*fr[22]-19.6875*fr[10]; 
  incr2[13] = -19.6875*fr[13]; 
  incr2[14] = -19.6875*fr[14]; 
  incr2[15] = 399.3353853592241*fr[27]-295.3125*fr[15]+76.24935962845854*fr[5]; 
  incr2[16] = 399.3353853592241*fr[28]-295.3125*fr[16]+76.24935962845854*fr[6]; 
  incr2[19] = (-945.0*fr[19])+698.8369243786424*fr[9]-180.4389179888861*fr[3]; 
  incr2[20] = -19.6875*fr[20]; 
  incr2[21] = -19.6875*fr[21]; 
  incr2[22] = 399.3353853592242*fr[31]-295.3125*fr[22]+76.2493596284585*fr[10]; 
  incr2[25] = -19.6875*fr[25]; 
  incr2[26] = -19.6875*fr[26]; 
  incr2[27] = (-945.0*fr[27])+698.8369243786422*fr[15]-180.4389179888862*fr[5]; 
  incr2[28] = (-945.0*fr[28])+698.8369243786422*fr[16]-180.4389179888862*fr[6]; 
  incr2[29] = -19.6875*fr[29]; 
  incr2[30] = -19.6875*fr[30]; 
  incr2[31] = (-945.0*fr[31])+698.8369243786424*fr[22]-180.4389179888861*fr[10]; 


  incr4[19] = 225.0*fr[19]-241.2651286545313*fr[9]+96.66370606547471*fr[3]; 
  incr4[27] = 225.0*fr[27]-241.2651286545313*fr[15]+96.66370606547474*fr[5]; 
  incr4[28] = 225.0*fr[28]-241.2651286545313*fr[16]+96.66370606547474*fr[6]; 
  incr4[31] = 225.0*fr[31]-241.2651286545313*fr[22]+96.66370606547471*fr[10]; 



  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[19] += incr4[19]*rdxFnur+incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[27] += incr4[27]*rdxFnur+incr2[27]*rdxFnur; 
  outr[28] += incr4[28]*rdxFnur+incr2[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr4[31]*rdxFnur+incr2[31]*rdxFnur; 

  } else {


  incr2[3] = (-103.1079531365064*fl[19])-76.2493596284585*fl[9]-19.6875*fl[3]; 
  incr2[5] = (-103.1079531365064*fl[27])-76.24935962845854*fl[15]-19.6875*fl[5]; 
  incr2[6] = (-103.1079531365064*fl[28])-76.24935962845854*fl[16]-19.6875*fl[6]; 
  incr2[9] = (-399.3353853592242*fl[19])-295.3125*fl[9]-76.2493596284585*fl[3]; 
  incr2[10] = (-103.1079531365064*fl[31])-76.2493596284585*fl[22]-19.6875*fl[10]; 
  incr2[13] = -19.6875*fl[13]; 
  incr2[14] = -19.6875*fl[14]; 
  incr2[15] = (-399.3353853592241*fl[27])-295.3125*fl[15]-76.24935962845854*fl[5]; 
  incr2[16] = (-399.3353853592241*fl[28])-295.3125*fl[16]-76.24935962845854*fl[6]; 
  incr2[19] = (-945.0*fl[19])-698.8369243786424*fl[9]-180.4389179888861*fl[3]; 
  incr2[20] = -19.6875*fl[20]; 
  incr2[21] = -19.6875*fl[21]; 
  incr2[22] = (-399.3353853592242*fl[31])-295.3125*fl[22]-76.2493596284585*fl[10]; 
  incr2[25] = -19.6875*fl[25]; 
  incr2[26] = -19.6875*fl[26]; 
  incr2[27] = (-945.0*fl[27])-698.8369243786422*fl[15]-180.4389179888862*fl[5]; 
  incr2[28] = (-945.0*fl[28])-698.8369243786422*fl[16]-180.4389179888862*fl[6]; 
  incr2[29] = -19.6875*fl[29]; 
  incr2[30] = -19.6875*fl[30]; 
  incr2[31] = (-945.0*fl[31])-698.8369243786424*fl[22]-180.4389179888861*fl[10]; 


  incr4[19] = 225.0*fl[19]+241.2651286545313*fl[9]+96.66370606547471*fl[3]; 
  incr4[27] = 225.0*fl[27]+241.2651286545313*fl[15]+96.66370606547474*fl[5]; 
  incr4[28] = 225.0*fl[28]+241.2651286545313*fl[16]+96.66370606547474*fl[6]; 
  incr4[31] = 225.0*fl[31]+241.2651286545313*fl[22]+96.66370606547471*fl[10]; 



  outl[3] += incr2[3]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[19] += incr4[19]*rdxFnul+incr2[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul; 
  outl[27] += incr4[27]*rdxFnul+incr2[27]*rdxFnul; 
  outl[28] += incr4[28]*rdxFnul+incr2[28]*rdxFnul; 
  outl[29] += incr2[29]*rdxFnul; 
  outl[30] += incr2[30]*rdxFnul; 
  outl[31] += incr4[31]*rdxFnul+incr2[31]*rdxFnul; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf3xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[96]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[32]; 
  double incr2[32]; 

  if (edge < 0) {

  incr2[1] = (-2.14330352493528*fr[17]*nul[17])+1.811422093273679*fr[7]*nul[17]-1.403121520040228*fr[1]*nul[17]+0.8100925873009822*fr[0]*nul[17]-1.81142209327368*nul[7]*fr[17]-1.403121520040228*nul[1]*fr[17]-0.8100925873009825*nul[0]*fr[17]+1.530931089239485*fr[7]*nul[7]-1.185854122563142*fr[1]*nul[7]+0.6846531968814573*fr[0]*nul[7]+1.185854122563142*nul[1]*fr[7]+0.6846531968814573*nul[0]*fr[7]-0.9185586535436913*fr[1]*nul[1]+0.5303300858899105*fr[0]*nul[1]-0.5303300858899105*nul[0]*fr[1]+0.3061862178478971*fr[0]*nul[0]; 
  incr2[4] = (-2.14330352493528*nul[17]*fr[23])-1.811422093273679*nul[7]*fr[23]-1.403121520040228*nul[1]*fr[23]-0.8100925873009822*nul[0]*fr[23]+1.811422093273679*fr[11]*nul[17]-1.403121520040228*fr[4]*nul[17]+0.8100925873009822*fr[2]*nul[17]+1.530931089239486*nul[7]*fr[11]+1.185854122563142*nul[1]*fr[11]+0.6846531968814574*nul[0]*fr[11]-1.185854122563142*fr[4]*nul[7]+0.6846531968814573*fr[2]*nul[7]-0.9185586535436913*nul[1]*fr[4]-0.5303300858899105*nul[0]*fr[4]+0.5303300858899105*nul[1]*fr[2]+0.3061862178478971*nul[0]*fr[2]; 
  incr2[5] = (-2.14330352493528*nul[17]*fr[25])-1.811422093273679*nul[7]*fr[25]-1.403121520040228*nul[1]*fr[25]-0.8100925873009822*nul[0]*fr[25]+1.811422093273679*fr[13]*nul[17]-1.403121520040228*fr[5]*nul[17]+0.8100925873009822*fr[3]*nul[17]+1.530931089239486*nul[7]*fr[13]+1.185854122563142*nul[1]*fr[13]+0.6846531968814574*nul[0]*fr[13]-1.185854122563142*fr[5]*nul[7]+0.6846531968814573*fr[3]*nul[7]-0.9185586535436913*nul[1]*fr[5]-0.5303300858899105*nul[0]*fr[5]+0.5303300858899105*nul[1]*fr[3]+0.3061862178478971*nul[0]*fr[3]; 
  incr2[7] = 8.300978857941997*fr[17]*nul[17]-7.015607600201137*fr[7]*nul[17]+5.434266279821037*fr[1]*nul[17]-3.137475099502783*fr[0]*nul[17]+7.015607600201137*nul[7]*fr[17]+5.43426627982104*nul[1]*fr[17]+3.137475099502784*nul[0]*fr[17]-5.929270612815711*fr[7]*nul[7]+4.592793267718456*fr[1]*nul[7]-2.651650429449552*fr[0]*nul[7]-4.592793267718455*nul[1]*fr[7]-2.651650429449552*nul[0]*fr[7]+3.557562367689425*fr[1]*nul[1]-2.053959590644372*fr[0]*nul[1]+2.053959590644372*nul[0]*fr[1]-1.185854122563142*fr[0]*nul[0]; 
  incr2[10] = (-2.14330352493528*nul[17]*fr[29])-1.81142209327368*nul[7]*fr[29]-1.403121520040228*nul[1]*fr[29]-0.8100925873009825*nul[0]*fr[29]+1.811422093273679*nul[17]*fr[20]+1.530931089239485*nul[7]*fr[20]+1.185854122563142*nul[1]*fr[20]+0.6846531968814573*nul[0]*fr[20]-1.403121520040228*fr[10]*nul[17]+0.8100925873009822*fr[6]*nul[17]-1.185854122563142*nul[7]*fr[10]-0.9185586535436913*nul[1]*fr[10]-0.5303300858899105*nul[0]*fr[10]+0.6846531968814573*fr[6]*nul[7]+0.5303300858899105*nul[1]*fr[6]+0.3061862178478971*nul[0]*fr[6]; 
  incr2[11] = 8.300978857941994*nul[17]*fr[23]+7.015607600201138*nul[7]*fr[23]+5.434266279821037*nul[1]*fr[23]+3.137475099502782*nul[0]*fr[23]-7.015607600201137*fr[11]*nul[17]+5.434266279821038*fr[4]*nul[17]-3.137475099502782*fr[2]*nul[17]-5.929270612815711*nul[7]*fr[11]-4.592793267718455*nul[1]*fr[11]-2.651650429449552*nul[0]*fr[11]+4.592793267718458*fr[4]*nul[7]-2.651650429449552*fr[2]*nul[7]+3.557562367689425*nul[1]*fr[4]+2.053959590644372*nul[0]*fr[4]-2.053959590644372*nul[1]*fr[2]-1.185854122563142*nul[0]*fr[2]; 
  incr2[12] = (-1.403121520040228*fr[12]*nul[17])+0.8100925873009821*fr[8]*nul[17]-1.185854122563142*nul[7]*fr[12]-0.9185586535436913*nul[1]*fr[12]-0.5303300858899105*nul[0]*fr[12]+0.6846531968814574*nul[7]*fr[8]+0.5303300858899104*nul[1]*fr[8]+0.3061862178478971*nul[0]*fr[8]; 
  incr2[13] = 8.300978857941994*nul[17]*fr[25]+7.015607600201138*nul[7]*fr[25]+5.434266279821037*nul[1]*fr[25]+3.137475099502782*nul[0]*fr[25]-7.015607600201137*fr[13]*nul[17]+5.434266279821038*fr[5]*nul[17]-3.137475099502782*fr[3]*nul[17]-5.929270612815711*nul[7]*fr[13]-4.592793267718455*nul[1]*fr[13]-2.651650429449552*nul[0]*fr[13]+4.592793267718458*fr[5]*nul[7]-2.651650429449552*fr[3]*nul[7]+3.557562367689425*nul[1]*fr[5]+2.053959590644372*nul[0]*fr[5]-2.053959590644372*nul[1]*fr[3]-1.185854122563142*nul[0]*fr[3]; 
  incr2[15] = (-1.403121520040228*fr[15]*nul[17])+0.8100925873009821*fr[9]*nul[17]-1.185854122563142*nul[7]*fr[15]-0.9185586535436913*nul[1]*fr[15]-0.5303300858899105*nul[0]*fr[15]+0.6846531968814574*nul[7]*fr[9]+0.5303300858899104*nul[1]*fr[9]+0.3061862178478971*nul[0]*fr[9]; 
  incr2[17] = (-19.64370128056319*fr[17]*nul[17])+16.60195771588399*fr[7]*nul[17]-12.85982114961168*fr[1]*nul[17]+7.424621202458747*fr[0]*nul[17]-16.60195771588399*nul[7]*fr[17]-12.85982114961168*nul[1]*fr[17]-7.424621202458747*nul[0]*fr[17]+14.03121520040228*fr[7]*nul[7]-10.86853255964208*fr[1]*nul[7]+6.274950199005565*fr[0]*nul[7]+10.86853255964207*nul[1]*fr[7]+6.274950199005565*nul[0]*fr[7]-8.418729120241366*fr[1]*nul[1]+4.860555523805894*fr[0]*nul[1]-4.860555523805894*nul[0]*fr[1]+2.806243040080455*fr[0]*nul[0]; 
  incr2[20] = 8.300978857941997*nul[17]*fr[29]+7.015607600201137*nul[7]*fr[29]+5.43426627982104*nul[1]*fr[29]+3.137475099502784*nul[0]*fr[29]-7.015607600201137*nul[17]*fr[20]-5.929270612815711*nul[7]*fr[20]-4.592793267718455*nul[1]*fr[20]-2.651650429449552*nul[0]*fr[20]+5.434266279821037*fr[10]*nul[17]-3.137475099502783*fr[6]*nul[17]+4.592793267718456*nul[7]*fr[10]+3.557562367689425*nul[1]*fr[10]+2.053959590644372*nul[0]*fr[10]-2.651650429449552*fr[6]*nul[7]-2.053959590644372*nul[1]*fr[6]-1.185854122563142*nul[0]*fr[6]; 
  incr2[21] = (-1.403121520040228*nul[17]*fr[21])-1.185854122563142*nul[7]*fr[21]-0.9185586535436913*nul[1]*fr[21]-0.5303300858899105*nul[0]*fr[21]+0.8100925873009821*fr[14]*nul[17]+0.6846531968814574*nul[7]*fr[14]+0.5303300858899104*nul[1]*fr[14]+0.3061862178478971*nul[0]*fr[14]; 
  incr2[22] = (-1.403121520040228*nul[17]*fr[22])-1.185854122563142*nul[7]*fr[22]-0.9185586535436913*nul[1]*fr[22]-0.5303300858899105*nul[0]*fr[22]+0.8100925873009821*fr[16]*nul[17]+0.6846531968814574*nul[7]*fr[16]+0.5303300858899104*nul[1]*fr[16]+0.3061862178478971*nul[0]*fr[16]; 
  incr2[23] = (-19.64370128056319*nul[17]*fr[23])-16.60195771588399*nul[7]*fr[23]-12.85982114961168*nul[1]*fr[23]-7.424621202458747*nul[0]*fr[23]+16.60195771588399*fr[11]*nul[17]-12.85982114961168*fr[4]*nul[17]+7.424621202458747*fr[2]*nul[17]+14.03121520040228*nul[7]*fr[11]+10.86853255964207*nul[1]*fr[11]+6.274950199005565*nul[0]*fr[11]-10.86853255964208*fr[4]*nul[7]+6.274950199005565*fr[2]*nul[7]-8.418729120241364*nul[1]*fr[4]-4.860555523805894*nul[0]*fr[4]+4.860555523805894*nul[1]*fr[2]+2.806243040080455*nul[0]*fr[2]; 
  incr2[24] = (-1.403121520040228*nul[17]*fr[24])-1.185854122563142*nul[7]*fr[24]-0.9185586535436913*nul[1]*fr[24]-0.5303300858899105*nul[0]*fr[24]+0.8100925873009822*nul[17]*fr[18]+0.6846531968814573*nul[7]*fr[18]+0.5303300858899104*nul[1]*fr[18]+0.3061862178478971*nul[0]*fr[18]; 
  incr2[25] = (-19.64370128056319*nul[17]*fr[25])-16.60195771588399*nul[7]*fr[25]-12.85982114961168*nul[1]*fr[25]-7.424621202458747*nul[0]*fr[25]+16.60195771588399*fr[13]*nul[17]-12.85982114961168*fr[5]*nul[17]+7.424621202458747*fr[3]*nul[17]+14.03121520040228*nul[7]*fr[13]+10.86853255964207*nul[1]*fr[13]+6.274950199005565*nul[0]*fr[13]-10.86853255964208*fr[5]*nul[7]+6.274950199005565*fr[3]*nul[7]-8.418729120241364*nul[1]*fr[5]-4.860555523805894*nul[0]*fr[5]+4.860555523805894*nul[1]*fr[3]+2.806243040080455*nul[0]*fr[3]; 
  incr2[27] = (-1.403121520040228*nul[17]*fr[27])-1.185854122563142*nul[7]*fr[27]-0.9185586535436913*nul[1]*fr[27]-0.5303300858899105*nul[0]*fr[27]+0.8100925873009822*nul[17]*fr[19]+0.6846531968814573*nul[7]*fr[19]+0.5303300858899104*nul[1]*fr[19]+0.3061862178478971*nul[0]*fr[19]; 
  incr2[29] = (-19.64370128056319*nul[17]*fr[29])-16.60195771588399*nul[7]*fr[29]-12.85982114961168*nul[1]*fr[29]-7.424621202458747*nul[0]*fr[29]+16.60195771588399*nul[17]*fr[20]+14.03121520040228*nul[7]*fr[20]+10.86853255964207*nul[1]*fr[20]+6.274950199005565*nul[0]*fr[20]-12.85982114961168*fr[10]*nul[17]+7.424621202458747*fr[6]*nul[17]-10.86853255964208*nul[7]*fr[10]-8.418729120241366*nul[1]*fr[10]-4.860555523805894*nul[0]*fr[10]+6.274950199005565*fr[6]*nul[7]+4.860555523805894*nul[1]*fr[6]+2.806243040080455*nul[0]*fr[6]; 
  incr2[30] = (-1.403121520040228*nul[17]*fr[30])-1.185854122563142*nul[7]*fr[30]-0.9185586535436913*nul[1]*fr[30]-0.5303300858899105*nul[0]*fr[30]+0.8100925873009822*nul[17]*fr[26]+0.6846531968814573*nul[7]*fr[26]+0.5303300858899104*nul[1]*fr[26]+0.3061862178478971*nul[0]*fr[26]; 
  incr2[31] = (-1.403121520040228*nul[17]*fr[31])-1.185854122563142*nul[7]*fr[31]-0.9185586535436913*nul[1]*fr[31]-0.5303300858899105*nul[0]*fr[31]+0.8100925873009822*nul[17]*fr[28]+0.6846531968814573*nul[7]*fr[28]+0.5303300858899104*nul[1]*fr[28]+0.3061862178478971*nul[0]*fr[28]; 

  outr[1] += incr2[1]*rdxFr; 
  outr[4] += incr2[4]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 
  outr[12] += incr2[12]*rdxFr; 
  outr[13] += incr2[13]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 
  outr[17] += incr2[17]*rdxFr; 
  outr[20] += incr2[20]*rdxFr; 
  outr[21] += incr2[21]*rdxFr; 
  outr[22] += incr2[22]*rdxFr; 
  outr[23] += incr2[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr; 
  outr[25] += incr2[25]*rdxFr; 
  outr[27] += incr2[27]*rdxFr; 
  outr[29] += incr2[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 

  } else {

  incr2[1] = 2.14330352493528*fl[17]*nul[17]+1.811422093273679*fl[7]*nul[17]+1.403121520040228*fl[1]*nul[17]+0.8100925873009822*fl[0]*nul[17]+1.81142209327368*nul[7]*fl[17]+1.403121520040228*nul[1]*fl[17]+0.8100925873009825*nul[0]*fl[17]+1.530931089239485*fl[7]*nul[7]+1.185854122563142*fl[1]*nul[7]+0.6846531968814573*fl[0]*nul[7]+1.185854122563142*nul[1]*fl[7]+0.6846531968814573*nul[0]*fl[7]+0.9185586535436913*fl[1]*nul[1]+0.5303300858899105*fl[0]*nul[1]+0.5303300858899105*nul[0]*fl[1]+0.3061862178478971*fl[0]*nul[0]; 
  incr2[4] = 2.14330352493528*nul[17]*fl[23]+1.811422093273679*nul[7]*fl[23]+1.403121520040228*nul[1]*fl[23]+0.8100925873009822*nul[0]*fl[23]+1.811422093273679*fl[11]*nul[17]+1.403121520040228*fl[4]*nul[17]+0.8100925873009822*fl[2]*nul[17]+1.530931089239486*nul[7]*fl[11]+1.185854122563142*nul[1]*fl[11]+0.6846531968814574*nul[0]*fl[11]+1.185854122563142*fl[4]*nul[7]+0.6846531968814573*fl[2]*nul[7]+0.9185586535436913*nul[1]*fl[4]+0.5303300858899105*nul[0]*fl[4]+0.5303300858899105*nul[1]*fl[2]+0.3061862178478971*nul[0]*fl[2]; 
  incr2[5] = 2.14330352493528*nul[17]*fl[25]+1.811422093273679*nul[7]*fl[25]+1.403121520040228*nul[1]*fl[25]+0.8100925873009822*nul[0]*fl[25]+1.811422093273679*fl[13]*nul[17]+1.403121520040228*fl[5]*nul[17]+0.8100925873009822*fl[3]*nul[17]+1.530931089239486*nul[7]*fl[13]+1.185854122563142*nul[1]*fl[13]+0.6846531968814574*nul[0]*fl[13]+1.185854122563142*fl[5]*nul[7]+0.6846531968814573*fl[3]*nul[7]+0.9185586535436913*nul[1]*fl[5]+0.5303300858899105*nul[0]*fl[5]+0.5303300858899105*nul[1]*fl[3]+0.3061862178478971*nul[0]*fl[3]; 
  incr2[7] = (-8.300978857941997*fl[17]*nul[17])-7.015607600201137*fl[7]*nul[17]-5.434266279821037*fl[1]*nul[17]-3.137475099502783*fl[0]*nul[17]-7.015607600201137*nul[7]*fl[17]-5.43426627982104*nul[1]*fl[17]-3.137475099502784*nul[0]*fl[17]-5.929270612815711*fl[7]*nul[7]-4.592793267718456*fl[1]*nul[7]-2.651650429449552*fl[0]*nul[7]-4.592793267718455*nul[1]*fl[7]-2.651650429449552*nul[0]*fl[7]-3.557562367689425*fl[1]*nul[1]-2.053959590644372*fl[0]*nul[1]-2.053959590644372*nul[0]*fl[1]-1.185854122563142*fl[0]*nul[0]; 
  incr2[10] = 2.14330352493528*nul[17]*fl[29]+1.81142209327368*nul[7]*fl[29]+1.403121520040228*nul[1]*fl[29]+0.8100925873009825*nul[0]*fl[29]+1.811422093273679*nul[17]*fl[20]+1.530931089239485*nul[7]*fl[20]+1.185854122563142*nul[1]*fl[20]+0.6846531968814573*nul[0]*fl[20]+1.403121520040228*fl[10]*nul[17]+0.8100925873009822*fl[6]*nul[17]+1.185854122563142*nul[7]*fl[10]+0.9185586535436913*nul[1]*fl[10]+0.5303300858899105*nul[0]*fl[10]+0.6846531968814573*fl[6]*nul[7]+0.5303300858899105*nul[1]*fl[6]+0.3061862178478971*nul[0]*fl[6]; 
  incr2[11] = (-8.300978857941994*nul[17]*fl[23])-7.015607600201138*nul[7]*fl[23]-5.434266279821037*nul[1]*fl[23]-3.137475099502782*nul[0]*fl[23]-7.015607600201137*fl[11]*nul[17]-5.434266279821038*fl[4]*nul[17]-3.137475099502782*fl[2]*nul[17]-5.929270612815711*nul[7]*fl[11]-4.592793267718455*nul[1]*fl[11]-2.651650429449552*nul[0]*fl[11]-4.592793267718458*fl[4]*nul[7]-2.651650429449552*fl[2]*nul[7]-3.557562367689425*nul[1]*fl[4]-2.053959590644372*nul[0]*fl[4]-2.053959590644372*nul[1]*fl[2]-1.185854122563142*nul[0]*fl[2]; 
  incr2[12] = 1.403121520040228*fl[12]*nul[17]+0.8100925873009821*fl[8]*nul[17]+1.185854122563142*nul[7]*fl[12]+0.9185586535436913*nul[1]*fl[12]+0.5303300858899105*nul[0]*fl[12]+0.6846531968814574*nul[7]*fl[8]+0.5303300858899104*nul[1]*fl[8]+0.3061862178478971*nul[0]*fl[8]; 
  incr2[13] = (-8.300978857941994*nul[17]*fl[25])-7.015607600201138*nul[7]*fl[25]-5.434266279821037*nul[1]*fl[25]-3.137475099502782*nul[0]*fl[25]-7.015607600201137*fl[13]*nul[17]-5.434266279821038*fl[5]*nul[17]-3.137475099502782*fl[3]*nul[17]-5.929270612815711*nul[7]*fl[13]-4.592793267718455*nul[1]*fl[13]-2.651650429449552*nul[0]*fl[13]-4.592793267718458*fl[5]*nul[7]-2.651650429449552*fl[3]*nul[7]-3.557562367689425*nul[1]*fl[5]-2.053959590644372*nul[0]*fl[5]-2.053959590644372*nul[1]*fl[3]-1.185854122563142*nul[0]*fl[3]; 
  incr2[15] = 1.403121520040228*fl[15]*nul[17]+0.8100925873009821*fl[9]*nul[17]+1.185854122563142*nul[7]*fl[15]+0.9185586535436913*nul[1]*fl[15]+0.5303300858899105*nul[0]*fl[15]+0.6846531968814574*nul[7]*fl[9]+0.5303300858899104*nul[1]*fl[9]+0.3061862178478971*nul[0]*fl[9]; 
  incr2[17] = 19.64370128056319*fl[17]*nul[17]+16.60195771588399*fl[7]*nul[17]+12.85982114961168*fl[1]*nul[17]+7.424621202458747*fl[0]*nul[17]+16.60195771588399*nul[7]*fl[17]+12.85982114961168*nul[1]*fl[17]+7.424621202458747*nul[0]*fl[17]+14.03121520040228*fl[7]*nul[7]+10.86853255964208*fl[1]*nul[7]+6.274950199005565*fl[0]*nul[7]+10.86853255964207*nul[1]*fl[7]+6.274950199005565*nul[0]*fl[7]+8.418729120241366*fl[1]*nul[1]+4.860555523805894*fl[0]*nul[1]+4.860555523805894*nul[0]*fl[1]+2.806243040080455*fl[0]*nul[0]; 
  incr2[20] = (-8.300978857941997*nul[17]*fl[29])-7.015607600201137*nul[7]*fl[29]-5.43426627982104*nul[1]*fl[29]-3.137475099502784*nul[0]*fl[29]-7.015607600201137*nul[17]*fl[20]-5.929270612815711*nul[7]*fl[20]-4.592793267718455*nul[1]*fl[20]-2.651650429449552*nul[0]*fl[20]-5.434266279821037*fl[10]*nul[17]-3.137475099502783*fl[6]*nul[17]-4.592793267718456*nul[7]*fl[10]-3.557562367689425*nul[1]*fl[10]-2.053959590644372*nul[0]*fl[10]-2.651650429449552*fl[6]*nul[7]-2.053959590644372*nul[1]*fl[6]-1.185854122563142*nul[0]*fl[6]; 
  incr2[21] = 1.403121520040228*nul[17]*fl[21]+1.185854122563142*nul[7]*fl[21]+0.9185586535436913*nul[1]*fl[21]+0.5303300858899105*nul[0]*fl[21]+0.8100925873009821*fl[14]*nul[17]+0.6846531968814574*nul[7]*fl[14]+0.5303300858899104*nul[1]*fl[14]+0.3061862178478971*nul[0]*fl[14]; 
  incr2[22] = 1.403121520040228*nul[17]*fl[22]+1.185854122563142*nul[7]*fl[22]+0.9185586535436913*nul[1]*fl[22]+0.5303300858899105*nul[0]*fl[22]+0.8100925873009821*fl[16]*nul[17]+0.6846531968814574*nul[7]*fl[16]+0.5303300858899104*nul[1]*fl[16]+0.3061862178478971*nul[0]*fl[16]; 
  incr2[23] = 19.64370128056319*nul[17]*fl[23]+16.60195771588399*nul[7]*fl[23]+12.85982114961168*nul[1]*fl[23]+7.424621202458747*nul[0]*fl[23]+16.60195771588399*fl[11]*nul[17]+12.85982114961168*fl[4]*nul[17]+7.424621202458747*fl[2]*nul[17]+14.03121520040228*nul[7]*fl[11]+10.86853255964207*nul[1]*fl[11]+6.274950199005565*nul[0]*fl[11]+10.86853255964208*fl[4]*nul[7]+6.274950199005565*fl[2]*nul[7]+8.418729120241364*nul[1]*fl[4]+4.860555523805894*nul[0]*fl[4]+4.860555523805894*nul[1]*fl[2]+2.806243040080455*nul[0]*fl[2]; 
  incr2[24] = 1.403121520040228*nul[17]*fl[24]+1.185854122563142*nul[7]*fl[24]+0.9185586535436913*nul[1]*fl[24]+0.5303300858899105*nul[0]*fl[24]+0.8100925873009822*nul[17]*fl[18]+0.6846531968814573*nul[7]*fl[18]+0.5303300858899104*nul[1]*fl[18]+0.3061862178478971*nul[0]*fl[18]; 
  incr2[25] = 19.64370128056319*nul[17]*fl[25]+16.60195771588399*nul[7]*fl[25]+12.85982114961168*nul[1]*fl[25]+7.424621202458747*nul[0]*fl[25]+16.60195771588399*fl[13]*nul[17]+12.85982114961168*fl[5]*nul[17]+7.424621202458747*fl[3]*nul[17]+14.03121520040228*nul[7]*fl[13]+10.86853255964207*nul[1]*fl[13]+6.274950199005565*nul[0]*fl[13]+10.86853255964208*fl[5]*nul[7]+6.274950199005565*fl[3]*nul[7]+8.418729120241364*nul[1]*fl[5]+4.860555523805894*nul[0]*fl[5]+4.860555523805894*nul[1]*fl[3]+2.806243040080455*nul[0]*fl[3]; 
  incr2[27] = 1.403121520040228*nul[17]*fl[27]+1.185854122563142*nul[7]*fl[27]+0.9185586535436913*nul[1]*fl[27]+0.5303300858899105*nul[0]*fl[27]+0.8100925873009822*nul[17]*fl[19]+0.6846531968814573*nul[7]*fl[19]+0.5303300858899104*nul[1]*fl[19]+0.3061862178478971*nul[0]*fl[19]; 
  incr2[29] = 19.64370128056319*nul[17]*fl[29]+16.60195771588399*nul[7]*fl[29]+12.85982114961168*nul[1]*fl[29]+7.424621202458747*nul[0]*fl[29]+16.60195771588399*nul[17]*fl[20]+14.03121520040228*nul[7]*fl[20]+10.86853255964207*nul[1]*fl[20]+6.274950199005565*nul[0]*fl[20]+12.85982114961168*fl[10]*nul[17]+7.424621202458747*fl[6]*nul[17]+10.86853255964208*nul[7]*fl[10]+8.418729120241366*nul[1]*fl[10]+4.860555523805894*nul[0]*fl[10]+6.274950199005565*fl[6]*nul[7]+4.860555523805894*nul[1]*fl[6]+2.806243040080455*nul[0]*fl[6]; 
  incr2[30] = 1.403121520040228*nul[17]*fl[30]+1.185854122563142*nul[7]*fl[30]+0.9185586535436913*nul[1]*fl[30]+0.5303300858899105*nul[0]*fl[30]+0.8100925873009822*nul[17]*fl[26]+0.6846531968814573*nul[7]*fl[26]+0.5303300858899104*nul[1]*fl[26]+0.3061862178478971*nul[0]*fl[26]; 
  incr2[31] = 1.403121520040228*nul[17]*fl[31]+1.185854122563142*nul[7]*fl[31]+0.9185586535436913*nul[1]*fl[31]+0.5303300858899105*nul[0]*fl[31]+0.8100925873009822*nul[17]*fl[28]+0.6846531968814573*nul[7]*fl[28]+0.5303300858899104*nul[1]*fl[28]+0.3061862178478971*nul[0]*fl[28]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 
  outl[4] += -1.0*incr2[4]*rdxFl; 
  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[7] += incr2[7]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[11] += incr2[11]*rdxFl; 
  outl[12] += -1.0*incr2[12]*rdxFl; 
  outl[13] += incr2[13]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 
  outl[17] += -1.0*incr2[17]*rdxFl; 
  outl[20] += incr2[20]*rdxFl; 
  outl[21] += -1.0*incr2[21]*rdxFl; 
  outl[22] += -1.0*incr2[22]*rdxFl; 
  outl[23] += -1.0*incr2[23]*rdxFl; 
  outl[24] += -1.0*incr2[24]*rdxFl; 
  outl[25] += -1.0*incr2[25]*rdxFl; 
  outl[27] += -1.0*incr2[27]*rdxFl; 
  outl[29] += -1.0*incr2[29]*rdxFl; 
  outl[30] += -1.0*incr2[30]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf3xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[96]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[1]*dxl[1]); 
  double rdxFr = 4.0/(dxr[1]*dxr[1]); 

  double incr1[32]; 
  double incr2[32]; 

  if (edge < 0) {

  incr2[2] = (-2.14330352493528*fr[18]*nul[50])+1.811422093273679*fr[8]*nul[50]-1.403121520040228*fr[2]*nul[50]+0.8100925873009822*fr[0]*nul[50]-1.81142209327368*fr[18]*nul[40]+1.530931089239485*fr[8]*nul[40]-1.185854122563142*fr[2]*nul[40]+0.6846531968814573*fr[0]*nul[40]-1.403121520040228*fr[18]*nul[34]+1.185854122563142*fr[8]*nul[34]-0.9185586535436913*fr[2]*nul[34]+0.5303300858899105*fr[0]*nul[34]-0.8100925873009825*fr[18]*nul[32]+0.6846531968814573*fr[8]*nul[32]-0.5303300858899105*fr[2]*nul[32]+0.3061862178478971*fr[0]*nul[32]; 
  incr2[4] = (-2.14330352493528*fr[24]*nul[50])+1.811422093273679*fr[12]*nul[50]-1.403121520040228*fr[4]*nul[50]+0.8100925873009822*fr[1]*nul[50]-1.811422093273679*fr[24]*nul[40]+1.530931089239486*fr[12]*nul[40]-1.185854122563142*fr[4]*nul[40]+0.6846531968814573*fr[1]*nul[40]-1.403121520040228*fr[24]*nul[34]+1.185854122563142*fr[12]*nul[34]-0.9185586535436913*fr[4]*nul[34]+0.5303300858899105*fr[1]*nul[34]-0.8100925873009822*fr[24]*nul[32]+0.6846531968814574*fr[12]*nul[32]-0.5303300858899105*fr[4]*nul[32]+0.3061862178478971*fr[1]*nul[32]; 
  incr2[6] = (-2.14330352493528*fr[26]*nul[50])+1.811422093273679*fr[14]*nul[50]-1.403121520040228*fr[6]*nul[50]+0.8100925873009822*fr[3]*nul[50]-1.811422093273679*fr[26]*nul[40]+1.530931089239486*fr[14]*nul[40]-1.185854122563142*fr[6]*nul[40]+0.6846531968814573*fr[3]*nul[40]-1.403121520040228*fr[26]*nul[34]+1.185854122563142*fr[14]*nul[34]-0.9185586535436913*fr[6]*nul[34]+0.5303300858899105*fr[3]*nul[34]-0.8100925873009822*fr[26]*nul[32]+0.6846531968814574*fr[14]*nul[32]-0.5303300858899105*fr[6]*nul[32]+0.3061862178478971*fr[3]*nul[32]; 
  incr2[8] = 8.300978857941997*fr[18]*nul[50]-7.015607600201137*fr[8]*nul[50]+5.434266279821037*fr[2]*nul[50]-3.137475099502783*fr[0]*nul[50]+7.015607600201137*fr[18]*nul[40]-5.929270612815711*fr[8]*nul[40]+4.592793267718456*fr[2]*nul[40]-2.651650429449552*fr[0]*nul[40]+5.43426627982104*fr[18]*nul[34]-4.592793267718455*fr[8]*nul[34]+3.557562367689425*fr[2]*nul[34]-2.053959590644372*fr[0]*nul[34]+3.137475099502784*fr[18]*nul[32]-2.651650429449552*fr[8]*nul[32]+2.053959590644372*fr[2]*nul[32]-1.185854122563142*fr[0]*nul[32]; 
  incr2[10] = (-2.14330352493528*fr[30]*nul[50])+1.811422093273679*fr[21]*nul[50]-1.403121520040228*fr[10]*nul[50]+0.8100925873009822*fr[5]*nul[50]-1.81142209327368*fr[30]*nul[40]+1.530931089239485*fr[21]*nul[40]-1.185854122563142*fr[10]*nul[40]+0.6846531968814573*fr[5]*nul[40]-1.403121520040228*fr[30]*nul[34]+1.185854122563142*fr[21]*nul[34]-0.9185586535436913*fr[10]*nul[34]+0.5303300858899105*fr[5]*nul[34]-0.8100925873009825*fr[30]*nul[32]+0.6846531968814573*fr[21]*nul[32]-0.5303300858899105*fr[10]*nul[32]+0.3061862178478971*fr[5]*nul[32]; 
  incr2[11] = (-1.403121520040228*fr[11]*nul[50])+0.8100925873009821*fr[7]*nul[50]-1.185854122563142*fr[11]*nul[40]+0.6846531968814574*fr[7]*nul[40]-0.9185586535436913*fr[11]*nul[34]+0.5303300858899104*fr[7]*nul[34]-0.5303300858899105*fr[11]*nul[32]+0.3061862178478971*fr[7]*nul[32]; 
  incr2[12] = 8.300978857941994*fr[24]*nul[50]-7.015607600201137*fr[12]*nul[50]+5.434266279821038*fr[4]*nul[50]-3.137475099502782*fr[1]*nul[50]+7.015607600201138*fr[24]*nul[40]-5.929270612815711*fr[12]*nul[40]+4.592793267718458*fr[4]*nul[40]-2.651650429449552*fr[1]*nul[40]+5.434266279821037*fr[24]*nul[34]-4.592793267718455*fr[12]*nul[34]+3.557562367689425*fr[4]*nul[34]-2.053959590644372*fr[1]*nul[34]+3.137475099502782*fr[24]*nul[32]-2.651650429449552*fr[12]*nul[32]+2.053959590644372*fr[4]*nul[32]-1.185854122563142*fr[1]*nul[32]; 
  incr2[14] = 8.300978857941994*fr[26]*nul[50]-7.015607600201137*fr[14]*nul[50]+5.434266279821038*fr[6]*nul[50]-3.137475099502782*fr[3]*nul[50]+7.015607600201138*fr[26]*nul[40]-5.929270612815711*fr[14]*nul[40]+4.592793267718458*fr[6]*nul[40]-2.651650429449552*fr[3]*nul[40]+5.434266279821037*fr[26]*nul[34]-4.592793267718455*fr[14]*nul[34]+3.557562367689425*fr[6]*nul[34]-2.053959590644372*fr[3]*nul[34]+3.137475099502782*fr[26]*nul[32]-2.651650429449552*fr[14]*nul[32]+2.053959590644372*fr[6]*nul[32]-1.185854122563142*fr[3]*nul[32]; 
  incr2[16] = (-1.403121520040228*fr[16]*nul[50])+0.8100925873009821*fr[9]*nul[50]-1.185854122563142*fr[16]*nul[40]+0.6846531968814574*fr[9]*nul[40]-0.9185586535436913*fr[16]*nul[34]+0.5303300858899104*fr[9]*nul[34]-0.5303300858899105*fr[16]*nul[32]+0.3061862178478971*fr[9]*nul[32]; 
  incr2[18] = (-19.64370128056319*fr[18]*nul[50])+16.60195771588399*fr[8]*nul[50]-12.85982114961168*fr[2]*nul[50]+7.424621202458747*fr[0]*nul[50]-16.60195771588399*fr[18]*nul[40]+14.03121520040228*fr[8]*nul[40]-10.86853255964208*fr[2]*nul[40]+6.274950199005565*fr[0]*nul[40]-12.85982114961168*fr[18]*nul[34]+10.86853255964207*fr[8]*nul[34]-8.418729120241366*fr[2]*nul[34]+4.860555523805894*fr[0]*nul[34]-7.424621202458747*fr[18]*nul[32]+6.274950199005565*fr[8]*nul[32]-4.860555523805894*fr[2]*nul[32]+2.806243040080455*fr[0]*nul[32]; 
  incr2[20] = (-1.403121520040228*fr[20]*nul[50])+0.8100925873009821*fr[13]*nul[50]-1.185854122563142*fr[20]*nul[40]+0.6846531968814574*fr[13]*nul[40]-0.9185586535436913*fr[20]*nul[34]+0.5303300858899104*fr[13]*nul[34]-0.5303300858899105*fr[20]*nul[32]+0.3061862178478971*fr[13]*nul[32]; 
  incr2[21] = 8.300978857941997*fr[30]*nul[50]-7.015607600201137*fr[21]*nul[50]+5.434266279821037*fr[10]*nul[50]-3.137475099502783*fr[5]*nul[50]+7.015607600201137*fr[30]*nul[40]-5.929270612815711*fr[21]*nul[40]+4.592793267718456*fr[10]*nul[40]-2.651650429449552*fr[5]*nul[40]+5.43426627982104*fr[30]*nul[34]-4.592793267718455*fr[21]*nul[34]+3.557562367689425*fr[10]*nul[34]-2.053959590644372*fr[5]*nul[34]+3.137475099502784*fr[30]*nul[32]-2.651650429449552*fr[21]*nul[32]+2.053959590644372*fr[10]*nul[32]-1.185854122563142*fr[5]*nul[32]; 
  incr2[22] = (-1.403121520040228*fr[22]*nul[50])+0.8100925873009821*fr[15]*nul[50]-1.185854122563142*fr[22]*nul[40]+0.6846531968814574*fr[15]*nul[40]-0.9185586535436913*fr[22]*nul[34]+0.5303300858899104*fr[15]*nul[34]-0.5303300858899105*fr[22]*nul[32]+0.3061862178478971*fr[15]*nul[32]; 
  incr2[23] = (-1.403121520040228*fr[23]*nul[50])+0.8100925873009822*fr[17]*nul[50]-1.185854122563142*fr[23]*nul[40]+0.6846531968814573*fr[17]*nul[40]-0.9185586535436913*fr[23]*nul[34]+0.5303300858899104*fr[17]*nul[34]-0.5303300858899105*fr[23]*nul[32]+0.3061862178478971*fr[17]*nul[32]; 
  incr2[24] = (-19.64370128056319*fr[24]*nul[50])+16.60195771588399*fr[12]*nul[50]-12.85982114961168*fr[4]*nul[50]+7.424621202458747*fr[1]*nul[50]-16.60195771588399*fr[24]*nul[40]+14.03121520040228*fr[12]*nul[40]-10.86853255964208*fr[4]*nul[40]+6.274950199005565*fr[1]*nul[40]-12.85982114961168*fr[24]*nul[34]+10.86853255964207*fr[12]*nul[34]-8.418729120241364*fr[4]*nul[34]+4.860555523805894*fr[1]*nul[34]-7.424621202458747*fr[24]*nul[32]+6.274950199005565*fr[12]*nul[32]-4.860555523805894*fr[4]*nul[32]+2.806243040080455*fr[1]*nul[32]; 
  incr2[26] = (-19.64370128056319*fr[26]*nul[50])+16.60195771588399*fr[14]*nul[50]-12.85982114961168*fr[6]*nul[50]+7.424621202458747*fr[3]*nul[50]-16.60195771588399*fr[26]*nul[40]+14.03121520040228*fr[14]*nul[40]-10.86853255964208*fr[6]*nul[40]+6.274950199005565*fr[3]*nul[40]-12.85982114961168*fr[26]*nul[34]+10.86853255964207*fr[14]*nul[34]-8.418729120241364*fr[6]*nul[34]+4.860555523805894*fr[3]*nul[34]-7.424621202458747*fr[26]*nul[32]+6.274950199005565*fr[14]*nul[32]-4.860555523805894*fr[6]*nul[32]+2.806243040080455*fr[3]*nul[32]; 
  incr2[28] = (-1.403121520040228*fr[28]*nul[50])+0.8100925873009822*fr[19]*nul[50]-1.185854122563142*fr[28]*nul[40]+0.6846531968814573*fr[19]*nul[40]-0.9185586535436913*fr[28]*nul[34]+0.5303300858899104*fr[19]*nul[34]-0.5303300858899105*fr[28]*nul[32]+0.3061862178478971*fr[19]*nul[32]; 
  incr2[29] = (-1.403121520040228*fr[29]*nul[50])+0.8100925873009822*fr[25]*nul[50]-1.185854122563142*fr[29]*nul[40]+0.6846531968814573*fr[25]*nul[40]-0.9185586535436913*fr[29]*nul[34]+0.5303300858899104*fr[25]*nul[34]-0.5303300858899105*fr[29]*nul[32]+0.3061862178478971*fr[25]*nul[32]; 
  incr2[30] = (-19.64370128056319*fr[30]*nul[50])+16.60195771588399*fr[21]*nul[50]-12.85982114961168*fr[10]*nul[50]+7.424621202458747*fr[5]*nul[50]-16.60195771588399*fr[30]*nul[40]+14.03121520040228*fr[21]*nul[40]-10.86853255964208*fr[10]*nul[40]+6.274950199005565*fr[5]*nul[40]-12.85982114961168*fr[30]*nul[34]+10.86853255964207*fr[21]*nul[34]-8.418729120241366*fr[10]*nul[34]+4.860555523805894*fr[5]*nul[34]-7.424621202458747*fr[30]*nul[32]+6.274950199005565*fr[21]*nul[32]-4.860555523805894*fr[10]*nul[32]+2.806243040080455*fr[5]*nul[32]; 
  incr2[31] = (-1.403121520040228*fr[31]*nul[50])+0.8100925873009822*fr[27]*nul[50]-1.185854122563142*fr[31]*nul[40]+0.6846531968814573*fr[27]*nul[40]-0.9185586535436913*fr[31]*nul[34]+0.5303300858899104*fr[27]*nul[34]-0.5303300858899105*fr[31]*nul[32]+0.3061862178478971*fr[27]*nul[32]; 

  outr[2] += incr2[2]*rdxFr; 
  outr[4] += incr2[4]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[8] += incr2[8]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 
  outr[12] += incr2[12]*rdxFr; 
  outr[14] += incr2[14]*rdxFr; 
  outr[16] += incr2[16]*rdxFr; 
  outr[18] += incr2[18]*rdxFr; 
  outr[20] += incr2[20]*rdxFr; 
  outr[21] += incr2[21]*rdxFr; 
  outr[22] += incr2[22]*rdxFr; 
  outr[23] += incr2[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 
  outr[28] += incr2[28]*rdxFr; 
  outr[29] += incr2[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 

  } else {

  incr2[2] = 2.14330352493528*fl[18]*nul[50]+1.811422093273679*fl[8]*nul[50]+1.403121520040228*fl[2]*nul[50]+0.8100925873009822*fl[0]*nul[50]+1.81142209327368*fl[18]*nul[40]+1.530931089239485*fl[8]*nul[40]+1.185854122563142*fl[2]*nul[40]+0.6846531968814573*fl[0]*nul[40]+1.403121520040228*fl[18]*nul[34]+1.185854122563142*fl[8]*nul[34]+0.9185586535436913*fl[2]*nul[34]+0.5303300858899105*fl[0]*nul[34]+0.8100925873009825*fl[18]*nul[32]+0.6846531968814573*fl[8]*nul[32]+0.5303300858899105*fl[2]*nul[32]+0.3061862178478971*fl[0]*nul[32]; 
  incr2[4] = 2.14330352493528*fl[24]*nul[50]+1.811422093273679*fl[12]*nul[50]+1.403121520040228*fl[4]*nul[50]+0.8100925873009822*fl[1]*nul[50]+1.811422093273679*fl[24]*nul[40]+1.530931089239486*fl[12]*nul[40]+1.185854122563142*fl[4]*nul[40]+0.6846531968814573*fl[1]*nul[40]+1.403121520040228*fl[24]*nul[34]+1.185854122563142*fl[12]*nul[34]+0.9185586535436913*fl[4]*nul[34]+0.5303300858899105*fl[1]*nul[34]+0.8100925873009822*fl[24]*nul[32]+0.6846531968814574*fl[12]*nul[32]+0.5303300858899105*fl[4]*nul[32]+0.3061862178478971*fl[1]*nul[32]; 
  incr2[6] = 2.14330352493528*fl[26]*nul[50]+1.811422093273679*fl[14]*nul[50]+1.403121520040228*fl[6]*nul[50]+0.8100925873009822*fl[3]*nul[50]+1.811422093273679*fl[26]*nul[40]+1.530931089239486*fl[14]*nul[40]+1.185854122563142*fl[6]*nul[40]+0.6846531968814573*fl[3]*nul[40]+1.403121520040228*fl[26]*nul[34]+1.185854122563142*fl[14]*nul[34]+0.9185586535436913*fl[6]*nul[34]+0.5303300858899105*fl[3]*nul[34]+0.8100925873009822*fl[26]*nul[32]+0.6846531968814574*fl[14]*nul[32]+0.5303300858899105*fl[6]*nul[32]+0.3061862178478971*fl[3]*nul[32]; 
  incr2[8] = (-8.300978857941997*fl[18]*nul[50])-7.015607600201137*fl[8]*nul[50]-5.434266279821037*fl[2]*nul[50]-3.137475099502783*fl[0]*nul[50]-7.015607600201137*fl[18]*nul[40]-5.929270612815711*fl[8]*nul[40]-4.592793267718456*fl[2]*nul[40]-2.651650429449552*fl[0]*nul[40]-5.43426627982104*fl[18]*nul[34]-4.592793267718455*fl[8]*nul[34]-3.557562367689425*fl[2]*nul[34]-2.053959590644372*fl[0]*nul[34]-3.137475099502784*fl[18]*nul[32]-2.651650429449552*fl[8]*nul[32]-2.053959590644372*fl[2]*nul[32]-1.185854122563142*fl[0]*nul[32]; 
  incr2[10] = 2.14330352493528*fl[30]*nul[50]+1.811422093273679*fl[21]*nul[50]+1.403121520040228*fl[10]*nul[50]+0.8100925873009822*fl[5]*nul[50]+1.81142209327368*fl[30]*nul[40]+1.530931089239485*fl[21]*nul[40]+1.185854122563142*fl[10]*nul[40]+0.6846531968814573*fl[5]*nul[40]+1.403121520040228*fl[30]*nul[34]+1.185854122563142*fl[21]*nul[34]+0.9185586535436913*fl[10]*nul[34]+0.5303300858899105*fl[5]*nul[34]+0.8100925873009825*fl[30]*nul[32]+0.6846531968814573*fl[21]*nul[32]+0.5303300858899105*fl[10]*nul[32]+0.3061862178478971*fl[5]*nul[32]; 
  incr2[11] = 1.403121520040228*fl[11]*nul[50]+0.8100925873009821*fl[7]*nul[50]+1.185854122563142*fl[11]*nul[40]+0.6846531968814574*fl[7]*nul[40]+0.9185586535436913*fl[11]*nul[34]+0.5303300858899104*fl[7]*nul[34]+0.5303300858899105*fl[11]*nul[32]+0.3061862178478971*fl[7]*nul[32]; 
  incr2[12] = (-8.300978857941994*fl[24]*nul[50])-7.015607600201137*fl[12]*nul[50]-5.434266279821038*fl[4]*nul[50]-3.137475099502782*fl[1]*nul[50]-7.015607600201138*fl[24]*nul[40]-5.929270612815711*fl[12]*nul[40]-4.592793267718458*fl[4]*nul[40]-2.651650429449552*fl[1]*nul[40]-5.434266279821037*fl[24]*nul[34]-4.592793267718455*fl[12]*nul[34]-3.557562367689425*fl[4]*nul[34]-2.053959590644372*fl[1]*nul[34]-3.137475099502782*fl[24]*nul[32]-2.651650429449552*fl[12]*nul[32]-2.053959590644372*fl[4]*nul[32]-1.185854122563142*fl[1]*nul[32]; 
  incr2[14] = (-8.300978857941994*fl[26]*nul[50])-7.015607600201137*fl[14]*nul[50]-5.434266279821038*fl[6]*nul[50]-3.137475099502782*fl[3]*nul[50]-7.015607600201138*fl[26]*nul[40]-5.929270612815711*fl[14]*nul[40]-4.592793267718458*fl[6]*nul[40]-2.651650429449552*fl[3]*nul[40]-5.434266279821037*fl[26]*nul[34]-4.592793267718455*fl[14]*nul[34]-3.557562367689425*fl[6]*nul[34]-2.053959590644372*fl[3]*nul[34]-3.137475099502782*fl[26]*nul[32]-2.651650429449552*fl[14]*nul[32]-2.053959590644372*fl[6]*nul[32]-1.185854122563142*fl[3]*nul[32]; 
  incr2[16] = 1.403121520040228*fl[16]*nul[50]+0.8100925873009821*fl[9]*nul[50]+1.185854122563142*fl[16]*nul[40]+0.6846531968814574*fl[9]*nul[40]+0.9185586535436913*fl[16]*nul[34]+0.5303300858899104*fl[9]*nul[34]+0.5303300858899105*fl[16]*nul[32]+0.3061862178478971*fl[9]*nul[32]; 
  incr2[18] = 19.64370128056319*fl[18]*nul[50]+16.60195771588399*fl[8]*nul[50]+12.85982114961168*fl[2]*nul[50]+7.424621202458747*fl[0]*nul[50]+16.60195771588399*fl[18]*nul[40]+14.03121520040228*fl[8]*nul[40]+10.86853255964208*fl[2]*nul[40]+6.274950199005565*fl[0]*nul[40]+12.85982114961168*fl[18]*nul[34]+10.86853255964207*fl[8]*nul[34]+8.418729120241366*fl[2]*nul[34]+4.860555523805894*fl[0]*nul[34]+7.424621202458747*fl[18]*nul[32]+6.274950199005565*fl[8]*nul[32]+4.860555523805894*fl[2]*nul[32]+2.806243040080455*fl[0]*nul[32]; 
  incr2[20] = 1.403121520040228*fl[20]*nul[50]+0.8100925873009821*fl[13]*nul[50]+1.185854122563142*fl[20]*nul[40]+0.6846531968814574*fl[13]*nul[40]+0.9185586535436913*fl[20]*nul[34]+0.5303300858899104*fl[13]*nul[34]+0.5303300858899105*fl[20]*nul[32]+0.3061862178478971*fl[13]*nul[32]; 
  incr2[21] = (-8.300978857941997*fl[30]*nul[50])-7.015607600201137*fl[21]*nul[50]-5.434266279821037*fl[10]*nul[50]-3.137475099502783*fl[5]*nul[50]-7.015607600201137*fl[30]*nul[40]-5.929270612815711*fl[21]*nul[40]-4.592793267718456*fl[10]*nul[40]-2.651650429449552*fl[5]*nul[40]-5.43426627982104*fl[30]*nul[34]-4.592793267718455*fl[21]*nul[34]-3.557562367689425*fl[10]*nul[34]-2.053959590644372*fl[5]*nul[34]-3.137475099502784*fl[30]*nul[32]-2.651650429449552*fl[21]*nul[32]-2.053959590644372*fl[10]*nul[32]-1.185854122563142*fl[5]*nul[32]; 
  incr2[22] = 1.403121520040228*fl[22]*nul[50]+0.8100925873009821*fl[15]*nul[50]+1.185854122563142*fl[22]*nul[40]+0.6846531968814574*fl[15]*nul[40]+0.9185586535436913*fl[22]*nul[34]+0.5303300858899104*fl[15]*nul[34]+0.5303300858899105*fl[22]*nul[32]+0.3061862178478971*fl[15]*nul[32]; 
  incr2[23] = 1.403121520040228*fl[23]*nul[50]+0.8100925873009822*fl[17]*nul[50]+1.185854122563142*fl[23]*nul[40]+0.6846531968814573*fl[17]*nul[40]+0.9185586535436913*fl[23]*nul[34]+0.5303300858899104*fl[17]*nul[34]+0.5303300858899105*fl[23]*nul[32]+0.3061862178478971*fl[17]*nul[32]; 
  incr2[24] = 19.64370128056319*fl[24]*nul[50]+16.60195771588399*fl[12]*nul[50]+12.85982114961168*fl[4]*nul[50]+7.424621202458747*fl[1]*nul[50]+16.60195771588399*fl[24]*nul[40]+14.03121520040228*fl[12]*nul[40]+10.86853255964208*fl[4]*nul[40]+6.274950199005565*fl[1]*nul[40]+12.85982114961168*fl[24]*nul[34]+10.86853255964207*fl[12]*nul[34]+8.418729120241364*fl[4]*nul[34]+4.860555523805894*fl[1]*nul[34]+7.424621202458747*fl[24]*nul[32]+6.274950199005565*fl[12]*nul[32]+4.860555523805894*fl[4]*nul[32]+2.806243040080455*fl[1]*nul[32]; 
  incr2[26] = 19.64370128056319*fl[26]*nul[50]+16.60195771588399*fl[14]*nul[50]+12.85982114961168*fl[6]*nul[50]+7.424621202458747*fl[3]*nul[50]+16.60195771588399*fl[26]*nul[40]+14.03121520040228*fl[14]*nul[40]+10.86853255964208*fl[6]*nul[40]+6.274950199005565*fl[3]*nul[40]+12.85982114961168*fl[26]*nul[34]+10.86853255964207*fl[14]*nul[34]+8.418729120241364*fl[6]*nul[34]+4.860555523805894*fl[3]*nul[34]+7.424621202458747*fl[26]*nul[32]+6.274950199005565*fl[14]*nul[32]+4.860555523805894*fl[6]*nul[32]+2.806243040080455*fl[3]*nul[32]; 
  incr2[28] = 1.403121520040228*fl[28]*nul[50]+0.8100925873009822*fl[19]*nul[50]+1.185854122563142*fl[28]*nul[40]+0.6846531968814573*fl[19]*nul[40]+0.9185586535436913*fl[28]*nul[34]+0.5303300858899104*fl[19]*nul[34]+0.5303300858899105*fl[28]*nul[32]+0.3061862178478971*fl[19]*nul[32]; 
  incr2[29] = 1.403121520040228*fl[29]*nul[50]+0.8100925873009822*fl[25]*nul[50]+1.185854122563142*fl[29]*nul[40]+0.6846531968814573*fl[25]*nul[40]+0.9185586535436913*fl[29]*nul[34]+0.5303300858899104*fl[25]*nul[34]+0.5303300858899105*fl[29]*nul[32]+0.3061862178478971*fl[25]*nul[32]; 
  incr2[30] = 19.64370128056319*fl[30]*nul[50]+16.60195771588399*fl[21]*nul[50]+12.85982114961168*fl[10]*nul[50]+7.424621202458747*fl[5]*nul[50]+16.60195771588399*fl[30]*nul[40]+14.03121520040228*fl[21]*nul[40]+10.86853255964208*fl[10]*nul[40]+6.274950199005565*fl[5]*nul[40]+12.85982114961168*fl[30]*nul[34]+10.86853255964207*fl[21]*nul[34]+8.418729120241366*fl[10]*nul[34]+4.860555523805894*fl[5]*nul[34]+7.424621202458747*fl[30]*nul[32]+6.274950199005565*fl[21]*nul[32]+4.860555523805894*fl[10]*nul[32]+2.806243040080455*fl[5]*nul[32]; 
  incr2[31] = 1.403121520040228*fl[31]*nul[50]+0.8100925873009822*fl[27]*nul[50]+1.185854122563142*fl[31]*nul[40]+0.6846531968814573*fl[27]*nul[40]+0.9185586535436913*fl[31]*nul[34]+0.5303300858899104*fl[27]*nul[34]+0.5303300858899105*fl[31]*nul[32]+0.3061862178478971*fl[27]*nul[32]; 

  outl[2] += -1.0*incr2[2]*rdxFl; 
  outl[4] += -1.0*incr2[4]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[8] += incr2[8]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[11] += -1.0*incr2[11]*rdxFl; 
  outl[12] += incr2[12]*rdxFl; 
  outl[14] += incr2[14]*rdxFl; 
  outl[16] += -1.0*incr2[16]*rdxFl; 
  outl[18] += -1.0*incr2[18]*rdxFl; 
  outl[20] += -1.0*incr2[20]*rdxFl; 
  outl[21] += incr2[21]*rdxFl; 
  outl[22] += -1.0*incr2[22]*rdxFl; 
  outl[23] += -1.0*incr2[23]*rdxFl; 
  outl[24] += -1.0*incr2[24]*rdxFl; 
  outl[26] += -1.0*incr2[26]*rdxFl; 
  outl[28] += -1.0*incr2[28]*rdxFl; 
  outl[29] += -1.0*incr2[29]*rdxFl; 
  outl[30] += -1.0*incr2[30]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf3xSerP3_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[96]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[2]*dxl[2]); 
  double rdxFr = 4.0/(dxr[2]*dxr[2]); 

  double incr1[32]; 
  double incr2[32]; 

  if (edge < 0) {

  incr2[3] = (-2.14330352493528*fr[19]*nul[83])+1.811422093273679*fr[9]*nul[83]-1.403121520040228*fr[3]*nul[83]+0.8100925873009822*fr[0]*nul[83]-1.81142209327368*fr[19]*nul[73]+1.530931089239485*fr[9]*nul[73]-1.185854122563142*fr[3]*nul[73]+0.6846531968814573*fr[0]*nul[73]-1.403121520040228*fr[19]*nul[67]+1.185854122563142*fr[9]*nul[67]-0.9185586535436913*fr[3]*nul[67]+0.5303300858899105*fr[0]*nul[67]-0.8100925873009825*fr[19]*nul[64]+0.6846531968814573*fr[9]*nul[64]-0.5303300858899105*fr[3]*nul[64]+0.3061862178478971*fr[0]*nul[64]; 
  incr2[5] = (-2.14330352493528*fr[27]*nul[83])+1.811422093273679*fr[15]*nul[83]-1.403121520040228*fr[5]*nul[83]+0.8100925873009822*fr[1]*nul[83]-1.811422093273679*fr[27]*nul[73]+1.530931089239486*fr[15]*nul[73]-1.185854122563142*fr[5]*nul[73]+0.6846531968814573*fr[1]*nul[73]-1.403121520040228*fr[27]*nul[67]+1.185854122563142*fr[15]*nul[67]-0.9185586535436913*fr[5]*nul[67]+0.5303300858899105*fr[1]*nul[67]-0.8100925873009822*fr[27]*nul[64]+0.6846531968814574*fr[15]*nul[64]-0.5303300858899105*fr[5]*nul[64]+0.3061862178478971*fr[1]*nul[64]; 
  incr2[6] = (-2.14330352493528*fr[28]*nul[83])+1.811422093273679*fr[16]*nul[83]-1.403121520040228*fr[6]*nul[83]+0.8100925873009822*fr[2]*nul[83]-1.811422093273679*fr[28]*nul[73]+1.530931089239486*fr[16]*nul[73]-1.185854122563142*fr[6]*nul[73]+0.6846531968814573*fr[2]*nul[73]-1.403121520040228*fr[28]*nul[67]+1.185854122563142*fr[16]*nul[67]-0.9185586535436913*fr[6]*nul[67]+0.5303300858899105*fr[2]*nul[67]-0.8100925873009822*fr[28]*nul[64]+0.6846531968814574*fr[16]*nul[64]-0.5303300858899105*fr[6]*nul[64]+0.3061862178478971*fr[2]*nul[64]; 
  incr2[9] = 8.300978857941997*fr[19]*nul[83]-7.015607600201137*fr[9]*nul[83]+5.434266279821037*fr[3]*nul[83]-3.137475099502783*fr[0]*nul[83]+7.015607600201137*fr[19]*nul[73]-5.929270612815711*fr[9]*nul[73]+4.592793267718456*fr[3]*nul[73]-2.651650429449552*fr[0]*nul[73]+5.43426627982104*fr[19]*nul[67]-4.592793267718455*fr[9]*nul[67]+3.557562367689425*fr[3]*nul[67]-2.053959590644372*fr[0]*nul[67]+3.137475099502784*fr[19]*nul[64]-2.651650429449552*fr[9]*nul[64]+2.053959590644372*fr[3]*nul[64]-1.185854122563142*fr[0]*nul[64]; 
  incr2[10] = (-2.14330352493528*fr[31]*nul[83])+1.811422093273679*fr[22]*nul[83]-1.403121520040228*fr[10]*nul[83]+0.8100925873009822*fr[4]*nul[83]-1.81142209327368*fr[31]*nul[73]+1.530931089239485*fr[22]*nul[73]-1.185854122563142*fr[10]*nul[73]+0.6846531968814573*fr[4]*nul[73]-1.403121520040228*fr[31]*nul[67]+1.185854122563142*fr[22]*nul[67]-0.9185586535436913*fr[10]*nul[67]+0.5303300858899105*fr[4]*nul[67]-0.8100925873009825*fr[31]*nul[64]+0.6846531968814573*fr[22]*nul[64]-0.5303300858899105*fr[10]*nul[64]+0.3061862178478971*fr[4]*nul[64]; 
  incr2[13] = (-1.403121520040228*fr[13]*nul[83])+0.8100925873009821*fr[7]*nul[83]-1.185854122563142*fr[13]*nul[73]+0.6846531968814574*fr[7]*nul[73]-0.9185586535436913*fr[13]*nul[67]+0.5303300858899104*fr[7]*nul[67]-0.5303300858899105*fr[13]*nul[64]+0.3061862178478971*fr[7]*nul[64]; 
  incr2[14] = (-1.403121520040228*fr[14]*nul[83])+0.8100925873009821*fr[8]*nul[83]-1.185854122563142*fr[14]*nul[73]+0.6846531968814574*fr[8]*nul[73]-0.9185586535436913*fr[14]*nul[67]+0.5303300858899104*fr[8]*nul[67]-0.5303300858899105*fr[14]*nul[64]+0.3061862178478971*fr[8]*nul[64]; 
  incr2[15] = 8.300978857941994*fr[27]*nul[83]-7.015607600201137*fr[15]*nul[83]+5.434266279821038*fr[5]*nul[83]-3.137475099502782*fr[1]*nul[83]+7.015607600201138*fr[27]*nul[73]-5.929270612815711*fr[15]*nul[73]+4.592793267718458*fr[5]*nul[73]-2.651650429449552*fr[1]*nul[73]+5.434266279821037*fr[27]*nul[67]-4.592793267718455*fr[15]*nul[67]+3.557562367689425*fr[5]*nul[67]-2.053959590644372*fr[1]*nul[67]+3.137475099502782*fr[27]*nul[64]-2.651650429449552*fr[15]*nul[64]+2.053959590644372*fr[5]*nul[64]-1.185854122563142*fr[1]*nul[64]; 
  incr2[16] = 8.300978857941994*fr[28]*nul[83]-7.015607600201137*fr[16]*nul[83]+5.434266279821038*fr[6]*nul[83]-3.137475099502782*fr[2]*nul[83]+7.015607600201138*fr[28]*nul[73]-5.929270612815711*fr[16]*nul[73]+4.592793267718458*fr[6]*nul[73]-2.651650429449552*fr[2]*nul[73]+5.434266279821037*fr[28]*nul[67]-4.592793267718455*fr[16]*nul[67]+3.557562367689425*fr[6]*nul[67]-2.053959590644372*fr[2]*nul[67]+3.137475099502782*fr[28]*nul[64]-2.651650429449552*fr[16]*nul[64]+2.053959590644372*fr[6]*nul[64]-1.185854122563142*fr[2]*nul[64]; 
  incr2[19] = (-19.64370128056319*fr[19]*nul[83])+16.60195771588399*fr[9]*nul[83]-12.85982114961168*fr[3]*nul[83]+7.424621202458747*fr[0]*nul[83]-16.60195771588399*fr[19]*nul[73]+14.03121520040228*fr[9]*nul[73]-10.86853255964208*fr[3]*nul[73]+6.274950199005565*fr[0]*nul[73]-12.85982114961168*fr[19]*nul[67]+10.86853255964207*fr[9]*nul[67]-8.418729120241366*fr[3]*nul[67]+4.860555523805894*fr[0]*nul[67]-7.424621202458747*fr[19]*nul[64]+6.274950199005565*fr[9]*nul[64]-4.860555523805894*fr[3]*nul[64]+2.806243040080455*fr[0]*nul[64]; 
  incr2[20] = (-1.403121520040228*fr[20]*nul[83])+0.8100925873009821*fr[11]*nul[83]-1.185854122563142*fr[20]*nul[73]+0.6846531968814574*fr[11]*nul[73]-0.9185586535436913*fr[20]*nul[67]+0.5303300858899104*fr[11]*nul[67]-0.5303300858899105*fr[20]*nul[64]+0.3061862178478971*fr[11]*nul[64]; 
  incr2[21] = (-1.403121520040228*fr[21]*nul[83])+0.8100925873009821*fr[12]*nul[83]-1.185854122563142*fr[21]*nul[73]+0.6846531968814574*fr[12]*nul[73]-0.9185586535436913*fr[21]*nul[67]+0.5303300858899104*fr[12]*nul[67]-0.5303300858899105*fr[21]*nul[64]+0.3061862178478971*fr[12]*nul[64]; 
  incr2[22] = 8.300978857941997*fr[31]*nul[83]-7.015607600201137*fr[22]*nul[83]+5.434266279821037*fr[10]*nul[83]-3.137475099502783*fr[4]*nul[83]+7.015607600201137*fr[31]*nul[73]-5.929270612815711*fr[22]*nul[73]+4.592793267718456*fr[10]*nul[73]-2.651650429449552*fr[4]*nul[73]+5.43426627982104*fr[31]*nul[67]-4.592793267718455*fr[22]*nul[67]+3.557562367689425*fr[10]*nul[67]-2.053959590644372*fr[4]*nul[67]+3.137475099502784*fr[31]*nul[64]-2.651650429449552*fr[22]*nul[64]+2.053959590644372*fr[10]*nul[64]-1.185854122563142*fr[4]*nul[64]; 
  incr2[25] = (-1.403121520040228*fr[25]*nul[83])+0.8100925873009822*fr[17]*nul[83]-1.185854122563142*fr[25]*nul[73]+0.6846531968814573*fr[17]*nul[73]-0.9185586535436913*fr[25]*nul[67]+0.5303300858899104*fr[17]*nul[67]-0.5303300858899105*fr[25]*nul[64]+0.3061862178478971*fr[17]*nul[64]; 
  incr2[26] = (-1.403121520040228*fr[26]*nul[83])+0.8100925873009822*fr[18]*nul[83]-1.185854122563142*fr[26]*nul[73]+0.6846531968814573*fr[18]*nul[73]-0.9185586535436913*fr[26]*nul[67]+0.5303300858899104*fr[18]*nul[67]-0.5303300858899105*fr[26]*nul[64]+0.3061862178478971*fr[18]*nul[64]; 
  incr2[27] = (-19.64370128056319*fr[27]*nul[83])+16.60195771588399*fr[15]*nul[83]-12.85982114961168*fr[5]*nul[83]+7.424621202458747*fr[1]*nul[83]-16.60195771588399*fr[27]*nul[73]+14.03121520040228*fr[15]*nul[73]-10.86853255964208*fr[5]*nul[73]+6.274950199005565*fr[1]*nul[73]-12.85982114961168*fr[27]*nul[67]+10.86853255964207*fr[15]*nul[67]-8.418729120241364*fr[5]*nul[67]+4.860555523805894*fr[1]*nul[67]-7.424621202458747*fr[27]*nul[64]+6.274950199005565*fr[15]*nul[64]-4.860555523805894*fr[5]*nul[64]+2.806243040080455*fr[1]*nul[64]; 
  incr2[28] = (-19.64370128056319*fr[28]*nul[83])+16.60195771588399*fr[16]*nul[83]-12.85982114961168*fr[6]*nul[83]+7.424621202458747*fr[2]*nul[83]-16.60195771588399*fr[28]*nul[73]+14.03121520040228*fr[16]*nul[73]-10.86853255964208*fr[6]*nul[73]+6.274950199005565*fr[2]*nul[73]-12.85982114961168*fr[28]*nul[67]+10.86853255964207*fr[16]*nul[67]-8.418729120241364*fr[6]*nul[67]+4.860555523805894*fr[2]*nul[67]-7.424621202458747*fr[28]*nul[64]+6.274950199005565*fr[16]*nul[64]-4.860555523805894*fr[6]*nul[64]+2.806243040080455*fr[2]*nul[64]; 
  incr2[29] = (-1.403121520040228*fr[29]*nul[83])+0.8100925873009822*fr[23]*nul[83]-1.185854122563142*fr[29]*nul[73]+0.6846531968814573*fr[23]*nul[73]-0.9185586535436913*fr[29]*nul[67]+0.5303300858899104*fr[23]*nul[67]-0.5303300858899105*fr[29]*nul[64]+0.3061862178478971*fr[23]*nul[64]; 
  incr2[30] = (-1.403121520040228*fr[30]*nul[83])+0.8100925873009822*fr[24]*nul[83]-1.185854122563142*fr[30]*nul[73]+0.6846531968814573*fr[24]*nul[73]-0.9185586535436913*fr[30]*nul[67]+0.5303300858899104*fr[24]*nul[67]-0.5303300858899105*fr[30]*nul[64]+0.3061862178478971*fr[24]*nul[64]; 
  incr2[31] = (-19.64370128056319*fr[31]*nul[83])+16.60195771588399*fr[22]*nul[83]-12.85982114961168*fr[10]*nul[83]+7.424621202458747*fr[4]*nul[83]-16.60195771588399*fr[31]*nul[73]+14.03121520040228*fr[22]*nul[73]-10.86853255964208*fr[10]*nul[73]+6.274950199005565*fr[4]*nul[73]-12.85982114961168*fr[31]*nul[67]+10.86853255964207*fr[22]*nul[67]-8.418729120241366*fr[10]*nul[67]+4.860555523805894*fr[4]*nul[67]-7.424621202458747*fr[31]*nul[64]+6.274950199005565*fr[22]*nul[64]-4.860555523805894*fr[10]*nul[64]+2.806243040080455*fr[4]*nul[64]; 

  outr[3] += incr2[3]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[9] += incr2[9]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[13] += incr2[13]*rdxFr; 
  outr[14] += incr2[14]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 
  outr[16] += incr2[16]*rdxFr; 
  outr[19] += incr2[19]*rdxFr; 
  outr[20] += incr2[20]*rdxFr; 
  outr[21] += incr2[21]*rdxFr; 
  outr[22] += incr2[22]*rdxFr; 
  outr[25] += incr2[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 
  outr[27] += incr2[27]*rdxFr; 
  outr[28] += incr2[28]*rdxFr; 
  outr[29] += incr2[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 

  } else {

  incr2[3] = 2.14330352493528*fl[19]*nul[83]+1.811422093273679*fl[9]*nul[83]+1.403121520040228*fl[3]*nul[83]+0.8100925873009822*fl[0]*nul[83]+1.81142209327368*fl[19]*nul[73]+1.530931089239485*fl[9]*nul[73]+1.185854122563142*fl[3]*nul[73]+0.6846531968814573*fl[0]*nul[73]+1.403121520040228*fl[19]*nul[67]+1.185854122563142*fl[9]*nul[67]+0.9185586535436913*fl[3]*nul[67]+0.5303300858899105*fl[0]*nul[67]+0.8100925873009825*fl[19]*nul[64]+0.6846531968814573*fl[9]*nul[64]+0.5303300858899105*fl[3]*nul[64]+0.3061862178478971*fl[0]*nul[64]; 
  incr2[5] = 2.14330352493528*fl[27]*nul[83]+1.811422093273679*fl[15]*nul[83]+1.403121520040228*fl[5]*nul[83]+0.8100925873009822*fl[1]*nul[83]+1.811422093273679*fl[27]*nul[73]+1.530931089239486*fl[15]*nul[73]+1.185854122563142*fl[5]*nul[73]+0.6846531968814573*fl[1]*nul[73]+1.403121520040228*fl[27]*nul[67]+1.185854122563142*fl[15]*nul[67]+0.9185586535436913*fl[5]*nul[67]+0.5303300858899105*fl[1]*nul[67]+0.8100925873009822*fl[27]*nul[64]+0.6846531968814574*fl[15]*nul[64]+0.5303300858899105*fl[5]*nul[64]+0.3061862178478971*fl[1]*nul[64]; 
  incr2[6] = 2.14330352493528*fl[28]*nul[83]+1.811422093273679*fl[16]*nul[83]+1.403121520040228*fl[6]*nul[83]+0.8100925873009822*fl[2]*nul[83]+1.811422093273679*fl[28]*nul[73]+1.530931089239486*fl[16]*nul[73]+1.185854122563142*fl[6]*nul[73]+0.6846531968814573*fl[2]*nul[73]+1.403121520040228*fl[28]*nul[67]+1.185854122563142*fl[16]*nul[67]+0.9185586535436913*fl[6]*nul[67]+0.5303300858899105*fl[2]*nul[67]+0.8100925873009822*fl[28]*nul[64]+0.6846531968814574*fl[16]*nul[64]+0.5303300858899105*fl[6]*nul[64]+0.3061862178478971*fl[2]*nul[64]; 
  incr2[9] = (-8.300978857941997*fl[19]*nul[83])-7.015607600201137*fl[9]*nul[83]-5.434266279821037*fl[3]*nul[83]-3.137475099502783*fl[0]*nul[83]-7.015607600201137*fl[19]*nul[73]-5.929270612815711*fl[9]*nul[73]-4.592793267718456*fl[3]*nul[73]-2.651650429449552*fl[0]*nul[73]-5.43426627982104*fl[19]*nul[67]-4.592793267718455*fl[9]*nul[67]-3.557562367689425*fl[3]*nul[67]-2.053959590644372*fl[0]*nul[67]-3.137475099502784*fl[19]*nul[64]-2.651650429449552*fl[9]*nul[64]-2.053959590644372*fl[3]*nul[64]-1.185854122563142*fl[0]*nul[64]; 
  incr2[10] = 2.14330352493528*fl[31]*nul[83]+1.811422093273679*fl[22]*nul[83]+1.403121520040228*fl[10]*nul[83]+0.8100925873009822*fl[4]*nul[83]+1.81142209327368*fl[31]*nul[73]+1.530931089239485*fl[22]*nul[73]+1.185854122563142*fl[10]*nul[73]+0.6846531968814573*fl[4]*nul[73]+1.403121520040228*fl[31]*nul[67]+1.185854122563142*fl[22]*nul[67]+0.9185586535436913*fl[10]*nul[67]+0.5303300858899105*fl[4]*nul[67]+0.8100925873009825*fl[31]*nul[64]+0.6846531968814573*fl[22]*nul[64]+0.5303300858899105*fl[10]*nul[64]+0.3061862178478971*fl[4]*nul[64]; 
  incr2[13] = 1.403121520040228*fl[13]*nul[83]+0.8100925873009821*fl[7]*nul[83]+1.185854122563142*fl[13]*nul[73]+0.6846531968814574*fl[7]*nul[73]+0.9185586535436913*fl[13]*nul[67]+0.5303300858899104*fl[7]*nul[67]+0.5303300858899105*fl[13]*nul[64]+0.3061862178478971*fl[7]*nul[64]; 
  incr2[14] = 1.403121520040228*fl[14]*nul[83]+0.8100925873009821*fl[8]*nul[83]+1.185854122563142*fl[14]*nul[73]+0.6846531968814574*fl[8]*nul[73]+0.9185586535436913*fl[14]*nul[67]+0.5303300858899104*fl[8]*nul[67]+0.5303300858899105*fl[14]*nul[64]+0.3061862178478971*fl[8]*nul[64]; 
  incr2[15] = (-8.300978857941994*fl[27]*nul[83])-7.015607600201137*fl[15]*nul[83]-5.434266279821038*fl[5]*nul[83]-3.137475099502782*fl[1]*nul[83]-7.015607600201138*fl[27]*nul[73]-5.929270612815711*fl[15]*nul[73]-4.592793267718458*fl[5]*nul[73]-2.651650429449552*fl[1]*nul[73]-5.434266279821037*fl[27]*nul[67]-4.592793267718455*fl[15]*nul[67]-3.557562367689425*fl[5]*nul[67]-2.053959590644372*fl[1]*nul[67]-3.137475099502782*fl[27]*nul[64]-2.651650429449552*fl[15]*nul[64]-2.053959590644372*fl[5]*nul[64]-1.185854122563142*fl[1]*nul[64]; 
  incr2[16] = (-8.300978857941994*fl[28]*nul[83])-7.015607600201137*fl[16]*nul[83]-5.434266279821038*fl[6]*nul[83]-3.137475099502782*fl[2]*nul[83]-7.015607600201138*fl[28]*nul[73]-5.929270612815711*fl[16]*nul[73]-4.592793267718458*fl[6]*nul[73]-2.651650429449552*fl[2]*nul[73]-5.434266279821037*fl[28]*nul[67]-4.592793267718455*fl[16]*nul[67]-3.557562367689425*fl[6]*nul[67]-2.053959590644372*fl[2]*nul[67]-3.137475099502782*fl[28]*nul[64]-2.651650429449552*fl[16]*nul[64]-2.053959590644372*fl[6]*nul[64]-1.185854122563142*fl[2]*nul[64]; 
  incr2[19] = 19.64370128056319*fl[19]*nul[83]+16.60195771588399*fl[9]*nul[83]+12.85982114961168*fl[3]*nul[83]+7.424621202458747*fl[0]*nul[83]+16.60195771588399*fl[19]*nul[73]+14.03121520040228*fl[9]*nul[73]+10.86853255964208*fl[3]*nul[73]+6.274950199005565*fl[0]*nul[73]+12.85982114961168*fl[19]*nul[67]+10.86853255964207*fl[9]*nul[67]+8.418729120241366*fl[3]*nul[67]+4.860555523805894*fl[0]*nul[67]+7.424621202458747*fl[19]*nul[64]+6.274950199005565*fl[9]*nul[64]+4.860555523805894*fl[3]*nul[64]+2.806243040080455*fl[0]*nul[64]; 
  incr2[20] = 1.403121520040228*fl[20]*nul[83]+0.8100925873009821*fl[11]*nul[83]+1.185854122563142*fl[20]*nul[73]+0.6846531968814574*fl[11]*nul[73]+0.9185586535436913*fl[20]*nul[67]+0.5303300858899104*fl[11]*nul[67]+0.5303300858899105*fl[20]*nul[64]+0.3061862178478971*fl[11]*nul[64]; 
  incr2[21] = 1.403121520040228*fl[21]*nul[83]+0.8100925873009821*fl[12]*nul[83]+1.185854122563142*fl[21]*nul[73]+0.6846531968814574*fl[12]*nul[73]+0.9185586535436913*fl[21]*nul[67]+0.5303300858899104*fl[12]*nul[67]+0.5303300858899105*fl[21]*nul[64]+0.3061862178478971*fl[12]*nul[64]; 
  incr2[22] = (-8.300978857941997*fl[31]*nul[83])-7.015607600201137*fl[22]*nul[83]-5.434266279821037*fl[10]*nul[83]-3.137475099502783*fl[4]*nul[83]-7.015607600201137*fl[31]*nul[73]-5.929270612815711*fl[22]*nul[73]-4.592793267718456*fl[10]*nul[73]-2.651650429449552*fl[4]*nul[73]-5.43426627982104*fl[31]*nul[67]-4.592793267718455*fl[22]*nul[67]-3.557562367689425*fl[10]*nul[67]-2.053959590644372*fl[4]*nul[67]-3.137475099502784*fl[31]*nul[64]-2.651650429449552*fl[22]*nul[64]-2.053959590644372*fl[10]*nul[64]-1.185854122563142*fl[4]*nul[64]; 
  incr2[25] = 1.403121520040228*fl[25]*nul[83]+0.8100925873009822*fl[17]*nul[83]+1.185854122563142*fl[25]*nul[73]+0.6846531968814573*fl[17]*nul[73]+0.9185586535436913*fl[25]*nul[67]+0.5303300858899104*fl[17]*nul[67]+0.5303300858899105*fl[25]*nul[64]+0.3061862178478971*fl[17]*nul[64]; 
  incr2[26] = 1.403121520040228*fl[26]*nul[83]+0.8100925873009822*fl[18]*nul[83]+1.185854122563142*fl[26]*nul[73]+0.6846531968814573*fl[18]*nul[73]+0.9185586535436913*fl[26]*nul[67]+0.5303300858899104*fl[18]*nul[67]+0.5303300858899105*fl[26]*nul[64]+0.3061862178478971*fl[18]*nul[64]; 
  incr2[27] = 19.64370128056319*fl[27]*nul[83]+16.60195771588399*fl[15]*nul[83]+12.85982114961168*fl[5]*nul[83]+7.424621202458747*fl[1]*nul[83]+16.60195771588399*fl[27]*nul[73]+14.03121520040228*fl[15]*nul[73]+10.86853255964208*fl[5]*nul[73]+6.274950199005565*fl[1]*nul[73]+12.85982114961168*fl[27]*nul[67]+10.86853255964207*fl[15]*nul[67]+8.418729120241364*fl[5]*nul[67]+4.860555523805894*fl[1]*nul[67]+7.424621202458747*fl[27]*nul[64]+6.274950199005565*fl[15]*nul[64]+4.860555523805894*fl[5]*nul[64]+2.806243040080455*fl[1]*nul[64]; 
  incr2[28] = 19.64370128056319*fl[28]*nul[83]+16.60195771588399*fl[16]*nul[83]+12.85982114961168*fl[6]*nul[83]+7.424621202458747*fl[2]*nul[83]+16.60195771588399*fl[28]*nul[73]+14.03121520040228*fl[16]*nul[73]+10.86853255964208*fl[6]*nul[73]+6.274950199005565*fl[2]*nul[73]+12.85982114961168*fl[28]*nul[67]+10.86853255964207*fl[16]*nul[67]+8.418729120241364*fl[6]*nul[67]+4.860555523805894*fl[2]*nul[67]+7.424621202458747*fl[28]*nul[64]+6.274950199005565*fl[16]*nul[64]+4.860555523805894*fl[6]*nul[64]+2.806243040080455*fl[2]*nul[64]; 
  incr2[29] = 1.403121520040228*fl[29]*nul[83]+0.8100925873009822*fl[23]*nul[83]+1.185854122563142*fl[29]*nul[73]+0.6846531968814573*fl[23]*nul[73]+0.9185586535436913*fl[29]*nul[67]+0.5303300858899104*fl[23]*nul[67]+0.5303300858899105*fl[29]*nul[64]+0.3061862178478971*fl[23]*nul[64]; 
  incr2[30] = 1.403121520040228*fl[30]*nul[83]+0.8100925873009822*fl[24]*nul[83]+1.185854122563142*fl[30]*nul[73]+0.6846531968814573*fl[24]*nul[73]+0.9185586535436913*fl[30]*nul[67]+0.5303300858899104*fl[24]*nul[67]+0.5303300858899105*fl[30]*nul[64]+0.3061862178478971*fl[24]*nul[64]; 
  incr2[31] = 19.64370128056319*fl[31]*nul[83]+16.60195771588399*fl[22]*nul[83]+12.85982114961168*fl[10]*nul[83]+7.424621202458747*fl[4]*nul[83]+16.60195771588399*fl[31]*nul[73]+14.03121520040228*fl[22]*nul[73]+10.86853255964208*fl[10]*nul[73]+6.274950199005565*fl[4]*nul[73]+12.85982114961168*fl[31]*nul[67]+10.86853255964207*fl[22]*nul[67]+8.418729120241366*fl[10]*nul[67]+4.860555523805894*fl[4]*nul[67]+7.424621202458747*fl[31]*nul[64]+6.274950199005565*fl[22]*nul[64]+4.860555523805894*fl[10]*nul[64]+2.806243040080455*fl[4]*nul[64]; 

  outl[3] += -1.0*incr2[3]*rdxFl; 
  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[9] += incr2[9]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[13] += -1.0*incr2[13]*rdxFl; 
  outl[14] += -1.0*incr2[14]*rdxFl; 
  outl[15] += incr2[15]*rdxFl; 
  outl[16] += incr2[16]*rdxFl; 
  outl[19] += -1.0*incr2[19]*rdxFl; 
  outl[20] += -1.0*incr2[20]*rdxFl; 
  outl[21] += -1.0*incr2[21]*rdxFl; 
  outl[22] += incr2[22]*rdxFl; 
  outl[25] += -1.0*incr2[25]*rdxFl; 
  outl[26] += -1.0*incr2[26]*rdxFl; 
  outl[27] += -1.0*incr2[27]*rdxFl; 
  outl[28] += -1.0*incr2[28]*rdxFl; 
  outl[29] += -1.0*incr2[29]*rdxFl; 
  outl[30] += -1.0*incr2[30]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 

  }

} 
