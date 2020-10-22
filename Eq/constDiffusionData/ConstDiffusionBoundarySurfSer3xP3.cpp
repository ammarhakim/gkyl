#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf3xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[0] == 1) {

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
void ConstDiffusionBoundarySurf3xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[1] == 1) {

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
void ConstDiffusionBoundarySurf3xSerP3_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[2] == 1) {

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
void ConstHyperDiffusion4BoundarySurf3xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (idxr[0] == 1) {

  incr2[1] = 5.809475019311125*fr[7]-34.3693177121688*fr[17]; 
  incr2[4] = 5.809475019311126*fr[11]-34.3693177121688*fr[23]; 
  incr2[5] = 5.809475019311126*fr[13]-34.3693177121688*fr[25]; 
  incr2[7] = 133.1117951197414*fr[17]-22.5*fr[7]; 
  incr2[10] = 5.809475019311125*fr[20]-34.3693177121688*fr[29]; 
  incr2[11] = 133.1117951197414*fr[23]-22.5*fr[11]; 
  incr2[13] = 133.1117951197414*fr[25]-22.5*fr[13]; 
  incr2[17] = 53.24471804789655*fr[7]-315.0*fr[17]; 
  incr2[20] = 133.1117951197414*fr[29]-22.5*fr[20]; 
  incr2[23] = 53.24471804789655*fr[11]-315.0*fr[23]; 
  incr2[25] = 53.24471804789655*fr[13]-315.0*fr[25]; 
  incr2[29] = 53.24471804789655*fr[20]-315.0*fr[29]; 

  incr3[7] = (-53.24471804789656*fr[17])+22.5*fr[7]-5.809475019311125*fr[1]; 
  incr3[11] = (-53.24471804789655*fr[23])+22.5*fr[11]-5.809475019311126*fr[4]; 
  incr3[13] = (-53.24471804789655*fr[25])+22.5*fr[13]-5.809475019311126*fr[5]; 
  incr3[17] = 315.0*fr[17]-133.1117951197414*fr[7]+34.3693177121688*fr[1]; 
  incr3[20] = (-53.24471804789656*fr[29])+22.5*fr[20]-5.809475019311125*fr[10]; 
  incr3[23] = 315.0*fr[23]-133.1117951197414*fr[11]+34.3693177121688*fr[4]; 
  incr3[25] = 315.0*fr[25]-133.1117951197414*fr[13]+34.3693177121688*fr[5]; 
  incr3[29] = 315.0*fr[29]-133.1117951197414*fr[20]+34.3693177121688*fr[10]; 

  incr4[17] = (-52.5*fr[17])+44.37059837324713*fr[7]-34.3693177121688*fr[1]+19.84313483298443*fr[0]; 
  incr4[23] = (-52.5*fr[23])+44.37059837324713*fr[11]-34.3693177121688*fr[4]+19.84313483298443*fr[2]; 
  incr4[25] = (-52.5*fr[25])+44.37059837324713*fr[13]-34.3693177121688*fr[5]+19.84313483298443*fr[3]; 
  incr4[29] = (-52.5*fr[29])+44.37059837324713*fr[20]-34.3693177121688*fr[10]+19.84313483298443*fr[6]; 

  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[7] += (-1.0*incr3[7]*rdxFnur)-1.0*incr2[7]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[11] += (-1.0*incr3[11]*rdxFnur)-1.0*incr2[11]*rdxFnur; 
  outr[13] += (-1.0*incr3[13]*rdxFnur)-1.0*incr2[13]*rdxFnur; 
  outr[17] += (-1.0*incr4[17]*rdxFnur)-1.0*incr3[17]*rdxFnur-1.0*incr2[17]*rdxFnur; 
  outr[20] += (-1.0*incr3[20]*rdxFnur)-1.0*incr2[20]*rdxFnur; 
  outr[23] += (-1.0*incr4[23]*rdxFnur)-1.0*incr3[23]*rdxFnur-1.0*incr2[23]*rdxFnur; 
  outr[25] += (-1.0*incr4[25]*rdxFnur)-1.0*incr3[25]*rdxFnur-1.0*incr2[25]*rdxFnur; 
  outr[29] += (-1.0*incr4[29]*rdxFnur)-1.0*incr3[29]*rdxFnur-1.0*incr2[29]*rdxFnur; 

  } else {

  incr2[1] = 34.3693177121688*fl[17]+5.809475019311125*fl[7]; 
  incr2[4] = 34.3693177121688*fl[23]+5.809475019311126*fl[11]; 
  incr2[5] = 34.3693177121688*fl[25]+5.809475019311126*fl[13]; 
  incr2[7] = (-133.1117951197414*fl[17])-22.5*fl[7]; 
  incr2[10] = 34.3693177121688*fl[29]+5.809475019311125*fl[20]; 
  incr2[11] = (-133.1117951197414*fl[23])-22.5*fl[11]; 
  incr2[13] = (-133.1117951197414*fl[25])-22.5*fl[13]; 
  incr2[17] = 315.0*fl[17]+53.24471804789655*fl[7]; 
  incr2[20] = (-133.1117951197414*fl[29])-22.5*fl[20]; 
  incr2[23] = 315.0*fl[23]+53.24471804789655*fl[11]; 
  incr2[25] = 315.0*fl[25]+53.24471804789655*fl[13]; 
  incr2[29] = 315.0*fl[29]+53.24471804789655*fl[20]; 

  incr3[7] = (-53.24471804789656*fl[17])-22.5*fl[7]-5.809475019311125*fl[1]; 
  incr3[11] = (-53.24471804789655*fl[23])-22.5*fl[11]-5.809475019311126*fl[4]; 
  incr3[13] = (-53.24471804789655*fl[25])-22.5*fl[13]-5.809475019311126*fl[5]; 
  incr3[17] = 315.0*fl[17]+133.1117951197414*fl[7]+34.3693177121688*fl[1]; 
  incr3[20] = (-53.24471804789656*fl[29])-22.5*fl[20]-5.809475019311125*fl[10]; 
  incr3[23] = 315.0*fl[23]+133.1117951197414*fl[11]+34.3693177121688*fl[4]; 
  incr3[25] = 315.0*fl[25]+133.1117951197414*fl[13]+34.3693177121688*fl[5]; 
  incr3[29] = 315.0*fl[29]+133.1117951197414*fl[20]+34.3693177121688*fl[10]; 

  incr4[17] = 52.5*fl[17]+44.37059837324713*fl[7]+34.3693177121688*fl[1]+19.84313483298443*fl[0]; 
  incr4[23] = 52.5*fl[23]+44.37059837324713*fl[11]+34.3693177121688*fl[4]+19.84313483298443*fl[2]; 
  incr4[25] = 52.5*fl[25]+44.37059837324713*fl[13]+34.3693177121688*fl[5]+19.84313483298443*fl[3]; 
  incr4[29] = 52.5*fl[29]+44.37059837324713*fl[20]+34.3693177121688*fl[10]+19.84313483298443*fl[6]; 

  outl[1] += incr2[1]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[7] += incr3[7]*rdxFnul-1.0*incr2[7]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[11] += incr3[11]*rdxFnul-1.0*incr2[11]*rdxFnul; 
  outl[13] += incr3[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 
  outl[17] += incr4[17]*rdxFnul-1.0*incr3[17]*rdxFnul+incr2[17]*rdxFnul; 
  outl[20] += incr3[20]*rdxFnul-1.0*incr2[20]*rdxFnul; 
  outl[23] += incr4[23]*rdxFnul-1.0*incr3[23]*rdxFnul+incr2[23]*rdxFnul; 
  outl[25] += incr4[25]*rdxFnul-1.0*incr3[25]*rdxFnul+incr2[25]*rdxFnul; 
  outl[29] += incr4[29]*rdxFnul-1.0*incr3[29]*rdxFnul+incr2[29]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (idxr[1] == 1) {

  incr2[2] = 5.809475019311125*fr[8]-34.3693177121688*fr[18]; 
  incr2[4] = 5.809475019311126*fr[12]-34.3693177121688*fr[24]; 
  incr2[6] = 5.809475019311126*fr[14]-34.3693177121688*fr[26]; 
  incr2[8] = 133.1117951197414*fr[18]-22.5*fr[8]; 
  incr2[10] = 5.809475019311125*fr[21]-34.3693177121688*fr[30]; 
  incr2[12] = 133.1117951197414*fr[24]-22.5*fr[12]; 
  incr2[14] = 133.1117951197414*fr[26]-22.5*fr[14]; 
  incr2[18] = 53.24471804789655*fr[8]-315.0*fr[18]; 
  incr2[21] = 133.1117951197414*fr[30]-22.5*fr[21]; 
  incr2[24] = 53.24471804789655*fr[12]-315.0*fr[24]; 
  incr2[26] = 53.24471804789655*fr[14]-315.0*fr[26]; 
  incr2[30] = 53.24471804789655*fr[21]-315.0*fr[30]; 

  incr3[8] = (-53.24471804789656*fr[18])+22.5*fr[8]-5.809475019311125*fr[2]; 
  incr3[12] = (-53.24471804789655*fr[24])+22.5*fr[12]-5.809475019311126*fr[4]; 
  incr3[14] = (-53.24471804789655*fr[26])+22.5*fr[14]-5.809475019311126*fr[6]; 
  incr3[18] = 315.0*fr[18]-133.1117951197414*fr[8]+34.3693177121688*fr[2]; 
  incr3[21] = (-53.24471804789656*fr[30])+22.5*fr[21]-5.809475019311125*fr[10]; 
  incr3[24] = 315.0*fr[24]-133.1117951197414*fr[12]+34.3693177121688*fr[4]; 
  incr3[26] = 315.0*fr[26]-133.1117951197414*fr[14]+34.3693177121688*fr[6]; 
  incr3[30] = 315.0*fr[30]-133.1117951197414*fr[21]+34.3693177121688*fr[10]; 

  incr4[18] = (-52.5*fr[18])+44.37059837324713*fr[8]-34.3693177121688*fr[2]+19.84313483298443*fr[0]; 
  incr4[24] = (-52.5*fr[24])+44.37059837324713*fr[12]-34.3693177121688*fr[4]+19.84313483298443*fr[1]; 
  incr4[26] = (-52.5*fr[26])+44.37059837324713*fr[14]-34.3693177121688*fr[6]+19.84313483298443*fr[3]; 
  incr4[30] = (-52.5*fr[30])+44.37059837324713*fr[21]-34.3693177121688*fr[10]+19.84313483298443*fr[5]; 

  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[8] += (-1.0*incr3[8]*rdxFnur)-1.0*incr2[8]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[12] += (-1.0*incr3[12]*rdxFnur)-1.0*incr2[12]*rdxFnur; 
  outr[14] += (-1.0*incr3[14]*rdxFnur)-1.0*incr2[14]*rdxFnur; 
  outr[18] += (-1.0*incr4[18]*rdxFnur)-1.0*incr3[18]*rdxFnur-1.0*incr2[18]*rdxFnur; 
  outr[21] += (-1.0*incr3[21]*rdxFnur)-1.0*incr2[21]*rdxFnur; 
  outr[24] += (-1.0*incr4[24]*rdxFnur)-1.0*incr3[24]*rdxFnur-1.0*incr2[24]*rdxFnur; 
  outr[26] += (-1.0*incr4[26]*rdxFnur)-1.0*incr3[26]*rdxFnur-1.0*incr2[26]*rdxFnur; 
  outr[30] += (-1.0*incr4[30]*rdxFnur)-1.0*incr3[30]*rdxFnur-1.0*incr2[30]*rdxFnur; 

  } else {

  incr2[2] = 34.3693177121688*fl[18]+5.809475019311125*fl[8]; 
  incr2[4] = 34.3693177121688*fl[24]+5.809475019311126*fl[12]; 
  incr2[6] = 34.3693177121688*fl[26]+5.809475019311126*fl[14]; 
  incr2[8] = (-133.1117951197414*fl[18])-22.5*fl[8]; 
  incr2[10] = 34.3693177121688*fl[30]+5.809475019311125*fl[21]; 
  incr2[12] = (-133.1117951197414*fl[24])-22.5*fl[12]; 
  incr2[14] = (-133.1117951197414*fl[26])-22.5*fl[14]; 
  incr2[18] = 315.0*fl[18]+53.24471804789655*fl[8]; 
  incr2[21] = (-133.1117951197414*fl[30])-22.5*fl[21]; 
  incr2[24] = 315.0*fl[24]+53.24471804789655*fl[12]; 
  incr2[26] = 315.0*fl[26]+53.24471804789655*fl[14]; 
  incr2[30] = 315.0*fl[30]+53.24471804789655*fl[21]; 

  incr3[8] = (-53.24471804789656*fl[18])-22.5*fl[8]-5.809475019311125*fl[2]; 
  incr3[12] = (-53.24471804789655*fl[24])-22.5*fl[12]-5.809475019311126*fl[4]; 
  incr3[14] = (-53.24471804789655*fl[26])-22.5*fl[14]-5.809475019311126*fl[6]; 
  incr3[18] = 315.0*fl[18]+133.1117951197414*fl[8]+34.3693177121688*fl[2]; 
  incr3[21] = (-53.24471804789656*fl[30])-22.5*fl[21]-5.809475019311125*fl[10]; 
  incr3[24] = 315.0*fl[24]+133.1117951197414*fl[12]+34.3693177121688*fl[4]; 
  incr3[26] = 315.0*fl[26]+133.1117951197414*fl[14]+34.3693177121688*fl[6]; 
  incr3[30] = 315.0*fl[30]+133.1117951197414*fl[21]+34.3693177121688*fl[10]; 

  incr4[18] = 52.5*fl[18]+44.37059837324713*fl[8]+34.3693177121688*fl[2]+19.84313483298443*fl[0]; 
  incr4[24] = 52.5*fl[24]+44.37059837324713*fl[12]+34.3693177121688*fl[4]+19.84313483298443*fl[1]; 
  incr4[26] = 52.5*fl[26]+44.37059837324713*fl[14]+34.3693177121688*fl[6]+19.84313483298443*fl[3]; 
  incr4[30] = 52.5*fl[30]+44.37059837324713*fl[21]+34.3693177121688*fl[10]+19.84313483298443*fl[5]; 

  outl[2] += incr2[2]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[8] += incr3[8]*rdxFnul-1.0*incr2[8]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[12] += incr3[12]*rdxFnul-1.0*incr2[12]*rdxFnul; 
  outl[14] += incr3[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 
  outl[18] += incr4[18]*rdxFnul-1.0*incr3[18]*rdxFnul+incr2[18]*rdxFnul; 
  outl[21] += incr3[21]*rdxFnul-1.0*incr2[21]*rdxFnul; 
  outl[24] += incr4[24]*rdxFnul-1.0*incr3[24]*rdxFnul+incr2[24]*rdxFnul; 
  outl[26] += incr4[26]*rdxFnul-1.0*incr3[26]*rdxFnul+incr2[26]*rdxFnul; 
  outl[30] += incr4[30]*rdxFnul-1.0*incr3[30]*rdxFnul+incr2[30]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xSerP3_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
  // nu[3]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (idxr[2] == 1) {

  incr2[3] = 5.809475019311125*fr[9]-34.3693177121688*fr[19]; 
  incr2[5] = 5.809475019311126*fr[15]-34.3693177121688*fr[27]; 
  incr2[6] = 5.809475019311126*fr[16]-34.3693177121688*fr[28]; 
  incr2[9] = 133.1117951197414*fr[19]-22.5*fr[9]; 
  incr2[10] = 5.809475019311125*fr[22]-34.3693177121688*fr[31]; 
  incr2[15] = 133.1117951197414*fr[27]-22.5*fr[15]; 
  incr2[16] = 133.1117951197414*fr[28]-22.5*fr[16]; 
  incr2[19] = 53.24471804789655*fr[9]-315.0*fr[19]; 
  incr2[22] = 133.1117951197414*fr[31]-22.5*fr[22]; 
  incr2[27] = 53.24471804789655*fr[15]-315.0*fr[27]; 
  incr2[28] = 53.24471804789655*fr[16]-315.0*fr[28]; 
  incr2[31] = 53.24471804789655*fr[22]-315.0*fr[31]; 

  incr3[9] = (-53.24471804789656*fr[19])+22.5*fr[9]-5.809475019311125*fr[3]; 
  incr3[15] = (-53.24471804789655*fr[27])+22.5*fr[15]-5.809475019311126*fr[5]; 
  incr3[16] = (-53.24471804789655*fr[28])+22.5*fr[16]-5.809475019311126*fr[6]; 
  incr3[19] = 315.0*fr[19]-133.1117951197414*fr[9]+34.3693177121688*fr[3]; 
  incr3[22] = (-53.24471804789656*fr[31])+22.5*fr[22]-5.809475019311125*fr[10]; 
  incr3[27] = 315.0*fr[27]-133.1117951197414*fr[15]+34.3693177121688*fr[5]; 
  incr3[28] = 315.0*fr[28]-133.1117951197414*fr[16]+34.3693177121688*fr[6]; 
  incr3[31] = 315.0*fr[31]-133.1117951197414*fr[22]+34.3693177121688*fr[10]; 

  incr4[19] = (-52.5*fr[19])+44.37059837324713*fr[9]-34.3693177121688*fr[3]+19.84313483298443*fr[0]; 
  incr4[27] = (-52.5*fr[27])+44.37059837324713*fr[15]-34.3693177121688*fr[5]+19.84313483298443*fr[1]; 
  incr4[28] = (-52.5*fr[28])+44.37059837324713*fr[16]-34.3693177121688*fr[6]+19.84313483298443*fr[2]; 
  incr4[31] = (-52.5*fr[31])+44.37059837324713*fr[22]-34.3693177121688*fr[10]+19.84313483298443*fr[4]; 

  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[9] += (-1.0*incr3[9]*rdxFnur)-1.0*incr2[9]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[15] += (-1.0*incr3[15]*rdxFnur)-1.0*incr2[15]*rdxFnur; 
  outr[16] += (-1.0*incr3[16]*rdxFnur)-1.0*incr2[16]*rdxFnur; 
  outr[19] += (-1.0*incr4[19]*rdxFnur)-1.0*incr3[19]*rdxFnur-1.0*incr2[19]*rdxFnur; 
  outr[22] += (-1.0*incr3[22]*rdxFnur)-1.0*incr2[22]*rdxFnur; 
  outr[27] += (-1.0*incr4[27]*rdxFnur)-1.0*incr3[27]*rdxFnur-1.0*incr2[27]*rdxFnur; 
  outr[28] += (-1.0*incr4[28]*rdxFnur)-1.0*incr3[28]*rdxFnur-1.0*incr2[28]*rdxFnur; 
  outr[31] += (-1.0*incr4[31]*rdxFnur)-1.0*incr3[31]*rdxFnur-1.0*incr2[31]*rdxFnur; 

  } else {

  incr2[3] = 34.3693177121688*fl[19]+5.809475019311125*fl[9]; 
  incr2[5] = 34.3693177121688*fl[27]+5.809475019311126*fl[15]; 
  incr2[6] = 34.3693177121688*fl[28]+5.809475019311126*fl[16]; 
  incr2[9] = (-133.1117951197414*fl[19])-22.5*fl[9]; 
  incr2[10] = 34.3693177121688*fl[31]+5.809475019311125*fl[22]; 
  incr2[15] = (-133.1117951197414*fl[27])-22.5*fl[15]; 
  incr2[16] = (-133.1117951197414*fl[28])-22.5*fl[16]; 
  incr2[19] = 315.0*fl[19]+53.24471804789655*fl[9]; 
  incr2[22] = (-133.1117951197414*fl[31])-22.5*fl[22]; 
  incr2[27] = 315.0*fl[27]+53.24471804789655*fl[15]; 
  incr2[28] = 315.0*fl[28]+53.24471804789655*fl[16]; 
  incr2[31] = 315.0*fl[31]+53.24471804789655*fl[22]; 

  incr3[9] = (-53.24471804789656*fl[19])-22.5*fl[9]-5.809475019311125*fl[3]; 
  incr3[15] = (-53.24471804789655*fl[27])-22.5*fl[15]-5.809475019311126*fl[5]; 
  incr3[16] = (-53.24471804789655*fl[28])-22.5*fl[16]-5.809475019311126*fl[6]; 
  incr3[19] = 315.0*fl[19]+133.1117951197414*fl[9]+34.3693177121688*fl[3]; 
  incr3[22] = (-53.24471804789656*fl[31])-22.5*fl[22]-5.809475019311125*fl[10]; 
  incr3[27] = 315.0*fl[27]+133.1117951197414*fl[15]+34.3693177121688*fl[5]; 
  incr3[28] = 315.0*fl[28]+133.1117951197414*fl[16]+34.3693177121688*fl[6]; 
  incr3[31] = 315.0*fl[31]+133.1117951197414*fl[22]+34.3693177121688*fl[10]; 

  incr4[19] = 52.5*fl[19]+44.37059837324713*fl[9]+34.3693177121688*fl[3]+19.84313483298443*fl[0]; 
  incr4[27] = 52.5*fl[27]+44.37059837324713*fl[15]+34.3693177121688*fl[5]+19.84313483298443*fl[1]; 
  incr4[28] = 52.5*fl[28]+44.37059837324713*fl[16]+34.3693177121688*fl[6]+19.84313483298443*fl[2]; 
  incr4[31] = 52.5*fl[31]+44.37059837324713*fl[22]+34.3693177121688*fl[10]+19.84313483298443*fl[4]; 

  outl[3] += incr2[3]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[9] += incr3[9]*rdxFnul-1.0*incr2[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[15] += incr3[15]*rdxFnul-1.0*incr2[15]*rdxFnul; 
  outl[16] += incr3[16]*rdxFnul-1.0*incr2[16]*rdxFnul; 
  outl[19] += incr4[19]*rdxFnul-1.0*incr3[19]*rdxFnul+incr2[19]*rdxFnul; 
  outl[22] += incr3[22]*rdxFnul-1.0*incr2[22]*rdxFnul; 
  outl[27] += incr4[27]*rdxFnul-1.0*incr3[27]*rdxFnul+incr2[27]*rdxFnul; 
  outl[28] += incr4[28]*rdxFnul-1.0*incr3[28]*rdxFnul+incr2[28]*rdxFnul; 
  outl[31] += incr4[31]*rdxFnul-1.0*incr3[31]*rdxFnul+incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf3xSerP3_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
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

  if (idxr[0] == 1) {


  incr3[7] = -133.1117951197414*fr[17]; 
  incr3[11] = -133.1117951197414*fr[23]; 
  incr3[13] = -133.1117951197414*fr[25]; 
  incr3[17] = 787.5*fr[17]; 
  incr3[20] = -133.1117951197414*fr[29]; 
  incr3[23] = 787.5*fr[23]; 
  incr3[25] = 787.5*fr[25]; 
  incr3[29] = 787.5*fr[29]; 

  incr4[17] = 133.1117951197414*fr[7]-787.5*fr[17]; 
  incr4[23] = 133.1117951197414*fr[11]-787.5*fr[23]; 
  incr4[25] = 133.1117951197414*fr[13]-787.5*fr[25]; 
  incr4[29] = 133.1117951197414*fr[20]-787.5*fr[29]; 



  outr[7] += incr3[7]*rdxFnur; 
  outr[11] += incr3[11]*rdxFnur; 
  outr[13] += incr3[13]*rdxFnur; 
  outr[17] += incr4[17]*rdxFnur+incr3[17]*rdxFnur; 
  outr[20] += incr3[20]*rdxFnur; 
  outr[23] += incr4[23]*rdxFnur+incr3[23]*rdxFnur; 
  outr[25] += incr4[25]*rdxFnur+incr3[25]*rdxFnur; 
  outr[29] += incr4[29]*rdxFnur+incr3[29]*rdxFnur; 

  } else {


  incr3[7] = -133.1117951197414*fl[17]; 
  incr3[11] = -133.1117951197414*fl[23]; 
  incr3[13] = -133.1117951197414*fl[25]; 
  incr3[17] = 787.5*fl[17]; 
  incr3[20] = -133.1117951197414*fl[29]; 
  incr3[23] = 787.5*fl[23]; 
  incr3[25] = 787.5*fl[25]; 
  incr3[29] = 787.5*fl[29]; 

  incr4[17] = 787.5*fl[17]+133.1117951197414*fl[7]; 
  incr4[23] = 787.5*fl[23]+133.1117951197414*fl[11]; 
  incr4[25] = 787.5*fl[25]+133.1117951197414*fl[13]; 
  incr4[29] = 787.5*fl[29]+133.1117951197414*fl[20]; 



  outl[7] += -1.0*incr3[7]*rdxFnul; 
  outl[11] += -1.0*incr3[11]*rdxFnul; 
  outl[13] += -1.0*incr3[13]*rdxFnul; 
  outl[17] += incr3[17]*rdxFnul-1.0*incr4[17]*rdxFnul; 
  outl[20] += -1.0*incr3[20]*rdxFnul; 
  outl[23] += incr3[23]*rdxFnul-1.0*incr4[23]*rdxFnul; 
  outl[25] += incr3[25]*rdxFnul-1.0*incr4[25]*rdxFnul; 
  outl[29] += incr3[29]*rdxFnul-1.0*incr4[29]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf3xSerP3_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
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

  if (idxr[1] == 1) {


  incr3[8] = -133.1117951197414*fr[18]; 
  incr3[12] = -133.1117951197414*fr[24]; 
  incr3[14] = -133.1117951197414*fr[26]; 
  incr3[18] = 787.5*fr[18]; 
  incr3[21] = -133.1117951197414*fr[30]; 
  incr3[24] = 787.5*fr[24]; 
  incr3[26] = 787.5*fr[26]; 
  incr3[30] = 787.5*fr[30]; 

  incr4[18] = 133.1117951197414*fr[8]-787.5*fr[18]; 
  incr4[24] = 133.1117951197414*fr[12]-787.5*fr[24]; 
  incr4[26] = 133.1117951197414*fr[14]-787.5*fr[26]; 
  incr4[30] = 133.1117951197414*fr[21]-787.5*fr[30]; 



  outr[8] += incr3[8]*rdxFnur; 
  outr[12] += incr3[12]*rdxFnur; 
  outr[14] += incr3[14]*rdxFnur; 
  outr[18] += incr4[18]*rdxFnur+incr3[18]*rdxFnur; 
  outr[21] += incr3[21]*rdxFnur; 
  outr[24] += incr4[24]*rdxFnur+incr3[24]*rdxFnur; 
  outr[26] += incr4[26]*rdxFnur+incr3[26]*rdxFnur; 
  outr[30] += incr4[30]*rdxFnur+incr3[30]*rdxFnur; 

  } else {


  incr3[8] = -133.1117951197414*fl[18]; 
  incr3[12] = -133.1117951197414*fl[24]; 
  incr3[14] = -133.1117951197414*fl[26]; 
  incr3[18] = 787.5*fl[18]; 
  incr3[21] = -133.1117951197414*fl[30]; 
  incr3[24] = 787.5*fl[24]; 
  incr3[26] = 787.5*fl[26]; 
  incr3[30] = 787.5*fl[30]; 

  incr4[18] = 787.5*fl[18]+133.1117951197414*fl[8]; 
  incr4[24] = 787.5*fl[24]+133.1117951197414*fl[12]; 
  incr4[26] = 787.5*fl[26]+133.1117951197414*fl[14]; 
  incr4[30] = 787.5*fl[30]+133.1117951197414*fl[21]; 



  outl[8] += -1.0*incr3[8]*rdxFnul; 
  outl[12] += -1.0*incr3[12]*rdxFnul; 
  outl[14] += -1.0*incr3[14]*rdxFnul; 
  outl[18] += incr3[18]*rdxFnul-1.0*incr4[18]*rdxFnul; 
  outl[21] += -1.0*incr3[21]*rdxFnul; 
  outl[24] += incr3[24]*rdxFnul-1.0*incr4[24]*rdxFnul; 
  outl[26] += incr3[26]*rdxFnul-1.0*incr4[26]*rdxFnul; 
  outl[30] += incr3[30]*rdxFnul-1.0*incr4[30]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf3xSerP3_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[3]:    current grid index.
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

  if (idxr[2] == 1) {


  incr3[9] = -133.1117951197414*fr[19]; 
  incr3[15] = -133.1117951197414*fr[27]; 
  incr3[16] = -133.1117951197414*fr[28]; 
  incr3[19] = 787.5*fr[19]; 
  incr3[22] = -133.1117951197414*fr[31]; 
  incr3[27] = 787.5*fr[27]; 
  incr3[28] = 787.5*fr[28]; 
  incr3[31] = 787.5*fr[31]; 

  incr4[19] = 133.1117951197414*fr[9]-787.5*fr[19]; 
  incr4[27] = 133.1117951197414*fr[15]-787.5*fr[27]; 
  incr4[28] = 133.1117951197414*fr[16]-787.5*fr[28]; 
  incr4[31] = 133.1117951197414*fr[22]-787.5*fr[31]; 



  outr[9] += incr3[9]*rdxFnur; 
  outr[15] += incr3[15]*rdxFnur; 
  outr[16] += incr3[16]*rdxFnur; 
  outr[19] += incr4[19]*rdxFnur+incr3[19]*rdxFnur; 
  outr[22] += incr3[22]*rdxFnur; 
  outr[27] += incr4[27]*rdxFnur+incr3[27]*rdxFnur; 
  outr[28] += incr4[28]*rdxFnur+incr3[28]*rdxFnur; 
  outr[31] += incr4[31]*rdxFnur+incr3[31]*rdxFnur; 

  } else {


  incr3[9] = -133.1117951197414*fl[19]; 
  incr3[15] = -133.1117951197414*fl[27]; 
  incr3[16] = -133.1117951197414*fl[28]; 
  incr3[19] = 787.5*fl[19]; 
  incr3[22] = -133.1117951197414*fl[31]; 
  incr3[27] = 787.5*fl[27]; 
  incr3[28] = 787.5*fl[28]; 
  incr3[31] = 787.5*fl[31]; 

  incr4[19] = 787.5*fl[19]+133.1117951197414*fl[9]; 
  incr4[27] = 787.5*fl[27]+133.1117951197414*fl[15]; 
  incr4[28] = 787.5*fl[28]+133.1117951197414*fl[16]; 
  incr4[31] = 787.5*fl[31]+133.1117951197414*fl[22]; 



  outl[9] += -1.0*incr3[9]*rdxFnul; 
  outl[15] += -1.0*incr3[15]*rdxFnul; 
  outl[16] += -1.0*incr3[16]*rdxFnul; 
  outl[19] += incr3[19]*rdxFnul-1.0*incr4[19]*rdxFnul; 
  outl[22] += -1.0*incr3[22]*rdxFnul; 
  outl[27] += incr3[27]*rdxFnul-1.0*incr4[27]*rdxFnul; 
  outl[28] += incr3[28]*rdxFnul-1.0*incr4[28]*rdxFnul; 
  outl[31] += incr3[31]*rdxFnul-1.0*incr4[31]*rdxFnul; 

  }

} 
