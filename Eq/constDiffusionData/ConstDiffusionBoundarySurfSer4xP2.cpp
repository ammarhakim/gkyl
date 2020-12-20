#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf4xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[48]; 
  double incr2[48]; 

  if (idxr[0] == 1) {

  incr2[1] = 1.936491673103708*fr[11]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[5] = 1.936491673103709*fr[19]-1.5*fr[5]+0.8660254037844386*fr[2]; 
  incr2[6] = 1.936491673103709*fr[21]-1.5*fr[6]+0.8660254037844386*fr[3]; 
  incr2[8] = 1.936491673103709*fr[25]-1.5*fr[8]+0.8660254037844386*fr[4]; 
  incr2[11] = (-7.5*fr[11])+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 
  incr2[15] = 1.936491673103708*fr[32]-1.5*fr[15]+0.8660254037844386*fr[7]; 
  incr2[16] = 1.936491673103708*fr[35]-1.5*fr[16]+0.8660254037844386*fr[9]; 
  incr2[17] = 1.936491673103708*fr[37]-1.5*fr[17]+0.8660254037844386*fr[10]; 
  incr2[19] = (-7.5*fr[19])+5.809475019311126*fr[5]-3.354101966249684*fr[2]; 
  incr2[20] = 0.8660254037844387*fr[12]-1.5*fr[20]; 
  incr2[21] = (-7.5*fr[21])+5.809475019311126*fr[6]-3.354101966249684*fr[3]; 
  incr2[23] = 0.8660254037844387*fr[13]-1.5*fr[23]; 
  incr2[25] = (-7.5*fr[25])+5.809475019311126*fr[8]-3.354101966249684*fr[4]; 
  incr2[28] = 0.8660254037844387*fr[14]-1.5*fr[28]; 
  incr2[31] = 1.936491673103709*fr[44]-1.5*fr[31]+0.8660254037844386*fr[18]; 
  incr2[32] = (-7.5*fr[32])+5.809475019311125*fr[15]-3.354101966249685*fr[7]; 
  incr2[33] = 0.8660254037844387*fr[22]-1.5*fr[33]; 
  incr2[34] = 0.8660254037844387*fr[24]-1.5*fr[34]; 
  incr2[35] = (-7.5*fr[35])+5.809475019311125*fr[16]-3.354101966249685*fr[9]; 
  incr2[36] = 0.8660254037844387*fr[26]-1.5*fr[36]; 
  incr2[37] = (-7.5*fr[37])+5.809475019311125*fr[17]-3.354101966249685*fr[10]; 
  incr2[39] = 0.8660254037844387*fr[27]-1.5*fr[39]; 
  incr2[41] = 0.8660254037844387*fr[29]-1.5*fr[41]; 
  incr2[42] = 0.8660254037844387*fr[30]-1.5*fr[42]; 
  incr2[44] = (-7.5*fr[44])+5.809475019311126*fr[31]-3.354101966249684*fr[18]; 
  incr2[45] = 0.8660254037844387*fr[38]-1.5*fr[45]; 
  incr2[46] = 0.8660254037844387*fr[40]-1.5*fr[46]; 
  incr2[47] = 0.8660254037844387*fr[43]-1.5*fr[47]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 
  outr[32] += incr2[32]*rdxFnur; 
  outr[33] += incr2[33]*rdxFnur; 
  outr[34] += incr2[34]*rdxFnur; 
  outr[35] += incr2[35]*rdxFnur; 
  outr[36] += incr2[36]*rdxFnur; 
  outr[37] += incr2[37]*rdxFnur; 
  outr[39] += incr2[39]*rdxFnur; 
  outr[41] += incr2[41]*rdxFnur; 
  outr[42] += incr2[42]*rdxFnur; 
  outr[44] += incr2[44]*rdxFnur; 
  outr[45] += incr2[45]*rdxFnur; 
  outr[46] += incr2[46]*rdxFnur; 
  outr[47] += incr2[47]*rdxFnur; 

  } else {

  incr2[1] = 1.936491673103708*fl[11]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[5] = 1.936491673103709*fl[19]+1.5*fl[5]+0.8660254037844386*fl[2]; 
  incr2[6] = 1.936491673103709*fl[21]+1.5*fl[6]+0.8660254037844386*fl[3]; 
  incr2[8] = 1.936491673103709*fl[25]+1.5*fl[8]+0.8660254037844386*fl[4]; 
  incr2[11] = (-7.5*fl[11])-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 
  incr2[15] = 1.936491673103708*fl[32]+1.5*fl[15]+0.8660254037844386*fl[7]; 
  incr2[16] = 1.936491673103708*fl[35]+1.5*fl[16]+0.8660254037844386*fl[9]; 
  incr2[17] = 1.936491673103708*fl[37]+1.5*fl[17]+0.8660254037844386*fl[10]; 
  incr2[19] = (-7.5*fl[19])-5.809475019311126*fl[5]-3.354101966249684*fl[2]; 
  incr2[20] = 1.5*fl[20]+0.8660254037844387*fl[12]; 
  incr2[21] = (-7.5*fl[21])-5.809475019311126*fl[6]-3.354101966249684*fl[3]; 
  incr2[23] = 1.5*fl[23]+0.8660254037844387*fl[13]; 
  incr2[25] = (-7.5*fl[25])-5.809475019311126*fl[8]-3.354101966249684*fl[4]; 
  incr2[28] = 1.5*fl[28]+0.8660254037844387*fl[14]; 
  incr2[31] = 1.936491673103709*fl[44]+1.5*fl[31]+0.8660254037844386*fl[18]; 
  incr2[32] = (-7.5*fl[32])-5.809475019311125*fl[15]-3.354101966249685*fl[7]; 
  incr2[33] = 1.5*fl[33]+0.8660254037844387*fl[22]; 
  incr2[34] = 1.5*fl[34]+0.8660254037844387*fl[24]; 
  incr2[35] = (-7.5*fl[35])-5.809475019311125*fl[16]-3.354101966249685*fl[9]; 
  incr2[36] = 1.5*fl[36]+0.8660254037844387*fl[26]; 
  incr2[37] = (-7.5*fl[37])-5.809475019311125*fl[17]-3.354101966249685*fl[10]; 
  incr2[39] = 1.5*fl[39]+0.8660254037844387*fl[27]; 
  incr2[41] = 1.5*fl[41]+0.8660254037844387*fl[29]; 
  incr2[42] = 1.5*fl[42]+0.8660254037844387*fl[30]; 
  incr2[44] = (-7.5*fl[44])-5.809475019311126*fl[31]-3.354101966249684*fl[18]; 
  incr2[45] = 1.5*fl[45]+0.8660254037844387*fl[38]; 
  incr2[46] = 1.5*fl[46]+0.8660254037844387*fl[40]; 
  incr2[47] = 1.5*fl[47]+0.8660254037844387*fl[43]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 
  outl[32] += incr2[32]*rdxFnul; 
  outl[33] += -1.0*incr2[33]*rdxFnul; 
  outl[34] += -1.0*incr2[34]*rdxFnul; 
  outl[35] += incr2[35]*rdxFnul; 
  outl[36] += -1.0*incr2[36]*rdxFnul; 
  outl[37] += incr2[37]*rdxFnul; 
  outl[39] += -1.0*incr2[39]*rdxFnul; 
  outl[41] += -1.0*incr2[41]*rdxFnul; 
  outl[42] += -1.0*incr2[42]*rdxFnul; 
  outl[44] += incr2[44]*rdxFnul; 
  outl[45] += -1.0*incr2[45]*rdxFnul; 
  outl[46] += -1.0*incr2[46]*rdxFnul; 
  outl[47] += -1.0*incr2[47]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf4xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[48]; 
  double incr2[48]; 

  if (idxr[1] == 1) {

  incr2[2] = 1.936491673103708*fr[12]-1.5*fr[2]+0.8660254037844386*fr[0]; 
  incr2[5] = 1.936491673103709*fr[20]-1.5*fr[5]+0.8660254037844386*fr[1]; 
  incr2[7] = 1.936491673103709*fr[22]-1.5*fr[7]+0.8660254037844386*fr[3]; 
  incr2[9] = 1.936491673103709*fr[26]-1.5*fr[9]+0.8660254037844386*fr[4]; 
  incr2[12] = (-7.5*fr[12])+5.809475019311125*fr[2]-3.354101966249685*fr[0]; 
  incr2[15] = 1.936491673103708*fr[33]-1.5*fr[15]+0.8660254037844386*fr[6]; 
  incr2[16] = 1.936491673103708*fr[36]-1.5*fr[16]+0.8660254037844386*fr[8]; 
  incr2[18] = 1.936491673103708*fr[38]-1.5*fr[18]+0.8660254037844386*fr[10]; 
  incr2[19] = 0.8660254037844387*fr[11]-1.5*fr[19]; 
  incr2[20] = (-7.5*fr[20])+5.809475019311126*fr[5]-3.354101966249684*fr[1]; 
  incr2[22] = (-7.5*fr[22])+5.809475019311126*fr[7]-3.354101966249684*fr[3]; 
  incr2[24] = 0.8660254037844387*fr[13]-1.5*fr[24]; 
  incr2[26] = (-7.5*fr[26])+5.809475019311126*fr[9]-3.354101966249684*fr[4]; 
  incr2[29] = 0.8660254037844387*fr[14]-1.5*fr[29]; 
  incr2[31] = 1.936491673103709*fr[45]-1.5*fr[31]+0.8660254037844386*fr[17]; 
  incr2[32] = 0.8660254037844387*fr[21]-1.5*fr[32]; 
  incr2[33] = (-7.5*fr[33])+5.809475019311125*fr[15]-3.354101966249685*fr[6]; 
  incr2[34] = 0.8660254037844387*fr[23]-1.5*fr[34]; 
  incr2[35] = 0.8660254037844387*fr[25]-1.5*fr[35]; 
  incr2[36] = (-7.5*fr[36])+5.809475019311125*fr[16]-3.354101966249685*fr[8]; 
  incr2[38] = (-7.5*fr[38])+5.809475019311125*fr[18]-3.354101966249685*fr[10]; 
  incr2[40] = 0.8660254037844387*fr[27]-1.5*fr[40]; 
  incr2[41] = 0.8660254037844387*fr[28]-1.5*fr[41]; 
  incr2[43] = 0.8660254037844387*fr[30]-1.5*fr[43]; 
  incr2[44] = 0.8660254037844387*fr[37]-1.5*fr[44]; 
  incr2[45] = (-7.5*fr[45])+5.809475019311126*fr[31]-3.354101966249684*fr[17]; 
  incr2[46] = 0.8660254037844387*fr[39]-1.5*fr[46]; 
  incr2[47] = 0.8660254037844387*fr[42]-1.5*fr[47]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 
  outr[32] += incr2[32]*rdxFnur; 
  outr[33] += incr2[33]*rdxFnur; 
  outr[34] += incr2[34]*rdxFnur; 
  outr[35] += incr2[35]*rdxFnur; 
  outr[36] += incr2[36]*rdxFnur; 
  outr[38] += incr2[38]*rdxFnur; 
  outr[40] += incr2[40]*rdxFnur; 
  outr[41] += incr2[41]*rdxFnur; 
  outr[43] += incr2[43]*rdxFnur; 
  outr[44] += incr2[44]*rdxFnur; 
  outr[45] += incr2[45]*rdxFnur; 
  outr[46] += incr2[46]*rdxFnur; 
  outr[47] += incr2[47]*rdxFnur; 

  } else {

  incr2[2] = 1.936491673103708*fl[12]+1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[5] = 1.936491673103709*fl[20]+1.5*fl[5]+0.8660254037844386*fl[1]; 
  incr2[7] = 1.936491673103709*fl[22]+1.5*fl[7]+0.8660254037844386*fl[3]; 
  incr2[9] = 1.936491673103709*fl[26]+1.5*fl[9]+0.8660254037844386*fl[4]; 
  incr2[12] = (-7.5*fl[12])-5.809475019311125*fl[2]-3.354101966249685*fl[0]; 
  incr2[15] = 1.936491673103708*fl[33]+1.5*fl[15]+0.8660254037844386*fl[6]; 
  incr2[16] = 1.936491673103708*fl[36]+1.5*fl[16]+0.8660254037844386*fl[8]; 
  incr2[18] = 1.936491673103708*fl[38]+1.5*fl[18]+0.8660254037844386*fl[10]; 
  incr2[19] = 1.5*fl[19]+0.8660254037844387*fl[11]; 
  incr2[20] = (-7.5*fl[20])-5.809475019311126*fl[5]-3.354101966249684*fl[1]; 
  incr2[22] = (-7.5*fl[22])-5.809475019311126*fl[7]-3.354101966249684*fl[3]; 
  incr2[24] = 1.5*fl[24]+0.8660254037844387*fl[13]; 
  incr2[26] = (-7.5*fl[26])-5.809475019311126*fl[9]-3.354101966249684*fl[4]; 
  incr2[29] = 1.5*fl[29]+0.8660254037844387*fl[14]; 
  incr2[31] = 1.936491673103709*fl[45]+1.5*fl[31]+0.8660254037844386*fl[17]; 
  incr2[32] = 1.5*fl[32]+0.8660254037844387*fl[21]; 
  incr2[33] = (-7.5*fl[33])-5.809475019311125*fl[15]-3.354101966249685*fl[6]; 
  incr2[34] = 1.5*fl[34]+0.8660254037844387*fl[23]; 
  incr2[35] = 1.5*fl[35]+0.8660254037844387*fl[25]; 
  incr2[36] = (-7.5*fl[36])-5.809475019311125*fl[16]-3.354101966249685*fl[8]; 
  incr2[38] = (-7.5*fl[38])-5.809475019311125*fl[18]-3.354101966249685*fl[10]; 
  incr2[40] = 1.5*fl[40]+0.8660254037844387*fl[27]; 
  incr2[41] = 1.5*fl[41]+0.8660254037844387*fl[28]; 
  incr2[43] = 1.5*fl[43]+0.8660254037844387*fl[30]; 
  incr2[44] = 1.5*fl[44]+0.8660254037844387*fl[37]; 
  incr2[45] = (-7.5*fl[45])-5.809475019311126*fl[31]-3.354101966249684*fl[17]; 
  incr2[46] = 1.5*fl[46]+0.8660254037844387*fl[39]; 
  incr2[47] = 1.5*fl[47]+0.8660254037844387*fl[42]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 
  outl[32] += -1.0*incr2[32]*rdxFnul; 
  outl[33] += incr2[33]*rdxFnul; 
  outl[34] += -1.0*incr2[34]*rdxFnul; 
  outl[35] += -1.0*incr2[35]*rdxFnul; 
  outl[36] += incr2[36]*rdxFnul; 
  outl[38] += incr2[38]*rdxFnul; 
  outl[40] += -1.0*incr2[40]*rdxFnul; 
  outl[41] += -1.0*incr2[41]*rdxFnul; 
  outl[43] += -1.0*incr2[43]*rdxFnul; 
  outl[44] += -1.0*incr2[44]*rdxFnul; 
  outl[45] += incr2[45]*rdxFnul; 
  outl[46] += -1.0*incr2[46]*rdxFnul; 
  outl[47] += -1.0*incr2[47]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf4xSerP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[48]; 
  double incr2[48]; 

  if (idxr[2] == 1) {

  incr2[3] = 1.936491673103708*fr[13]-1.5*fr[3]+0.8660254037844386*fr[0]; 
  incr2[6] = 1.936491673103709*fr[23]-1.5*fr[6]+0.8660254037844386*fr[1]; 
  incr2[7] = 1.936491673103709*fr[24]-1.5*fr[7]+0.8660254037844386*fr[2]; 
  incr2[10] = 1.936491673103709*fr[27]-1.5*fr[10]+0.8660254037844386*fr[4]; 
  incr2[13] = (-7.5*fr[13])+5.809475019311125*fr[3]-3.354101966249685*fr[0]; 
  incr2[15] = 1.936491673103708*fr[34]-1.5*fr[15]+0.8660254037844386*fr[5]; 
  incr2[17] = 1.936491673103708*fr[39]-1.5*fr[17]+0.8660254037844386*fr[8]; 
  incr2[18] = 1.936491673103708*fr[40]-1.5*fr[18]+0.8660254037844386*fr[9]; 
  incr2[21] = 0.8660254037844387*fr[11]-1.5*fr[21]; 
  incr2[22] = 0.8660254037844387*fr[12]-1.5*fr[22]; 
  incr2[23] = (-7.5*fr[23])+5.809475019311126*fr[6]-3.354101966249684*fr[1]; 
  incr2[24] = (-7.5*fr[24])+5.809475019311126*fr[7]-3.354101966249684*fr[2]; 
  incr2[27] = (-7.5*fr[27])+5.809475019311126*fr[10]-3.354101966249684*fr[4]; 
  incr2[30] = 0.8660254037844387*fr[14]-1.5*fr[30]; 
  incr2[31] = 1.936491673103709*fr[46]-1.5*fr[31]+0.8660254037844386*fr[16]; 
  incr2[32] = 0.8660254037844387*fr[19]-1.5*fr[32]; 
  incr2[33] = 0.8660254037844387*fr[20]-1.5*fr[33]; 
  incr2[34] = (-7.5*fr[34])+5.809475019311125*fr[15]-3.354101966249685*fr[5]; 
  incr2[37] = 0.8660254037844387*fr[25]-1.5*fr[37]; 
  incr2[38] = 0.8660254037844387*fr[26]-1.5*fr[38]; 
  incr2[39] = (-7.5*fr[39])+5.809475019311125*fr[17]-3.354101966249685*fr[8]; 
  incr2[40] = (-7.5*fr[40])+5.809475019311125*fr[18]-3.354101966249685*fr[9]; 
  incr2[42] = 0.8660254037844387*fr[28]-1.5*fr[42]; 
  incr2[43] = 0.8660254037844387*fr[29]-1.5*fr[43]; 
  incr2[44] = 0.8660254037844387*fr[35]-1.5*fr[44]; 
  incr2[45] = 0.8660254037844387*fr[36]-1.5*fr[45]; 
  incr2[46] = (-7.5*fr[46])+5.809475019311126*fr[31]-3.354101966249684*fr[16]; 
  incr2[47] = 0.8660254037844387*fr[41]-1.5*fr[47]; 

  outr[3] += incr2[3]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 
  outr[32] += incr2[32]*rdxFnur; 
  outr[33] += incr2[33]*rdxFnur; 
  outr[34] += incr2[34]*rdxFnur; 
  outr[37] += incr2[37]*rdxFnur; 
  outr[38] += incr2[38]*rdxFnur; 
  outr[39] += incr2[39]*rdxFnur; 
  outr[40] += incr2[40]*rdxFnur; 
  outr[42] += incr2[42]*rdxFnur; 
  outr[43] += incr2[43]*rdxFnur; 
  outr[44] += incr2[44]*rdxFnur; 
  outr[45] += incr2[45]*rdxFnur; 
  outr[46] += incr2[46]*rdxFnur; 
  outr[47] += incr2[47]*rdxFnur; 

  } else {

  incr2[3] = 1.936491673103708*fl[13]+1.5*fl[3]+0.8660254037844386*fl[0]; 
  incr2[6] = 1.936491673103709*fl[23]+1.5*fl[6]+0.8660254037844386*fl[1]; 
  incr2[7] = 1.936491673103709*fl[24]+1.5*fl[7]+0.8660254037844386*fl[2]; 
  incr2[10] = 1.936491673103709*fl[27]+1.5*fl[10]+0.8660254037844386*fl[4]; 
  incr2[13] = (-7.5*fl[13])-5.809475019311125*fl[3]-3.354101966249685*fl[0]; 
  incr2[15] = 1.936491673103708*fl[34]+1.5*fl[15]+0.8660254037844386*fl[5]; 
  incr2[17] = 1.936491673103708*fl[39]+1.5*fl[17]+0.8660254037844386*fl[8]; 
  incr2[18] = 1.936491673103708*fl[40]+1.5*fl[18]+0.8660254037844386*fl[9]; 
  incr2[21] = 1.5*fl[21]+0.8660254037844387*fl[11]; 
  incr2[22] = 1.5*fl[22]+0.8660254037844387*fl[12]; 
  incr2[23] = (-7.5*fl[23])-5.809475019311126*fl[6]-3.354101966249684*fl[1]; 
  incr2[24] = (-7.5*fl[24])-5.809475019311126*fl[7]-3.354101966249684*fl[2]; 
  incr2[27] = (-7.5*fl[27])-5.809475019311126*fl[10]-3.354101966249684*fl[4]; 
  incr2[30] = 1.5*fl[30]+0.8660254037844387*fl[14]; 
  incr2[31] = 1.936491673103709*fl[46]+1.5*fl[31]+0.8660254037844386*fl[16]; 
  incr2[32] = 1.5*fl[32]+0.8660254037844387*fl[19]; 
  incr2[33] = 1.5*fl[33]+0.8660254037844387*fl[20]; 
  incr2[34] = (-7.5*fl[34])-5.809475019311125*fl[15]-3.354101966249685*fl[5]; 
  incr2[37] = 1.5*fl[37]+0.8660254037844387*fl[25]; 
  incr2[38] = 1.5*fl[38]+0.8660254037844387*fl[26]; 
  incr2[39] = (-7.5*fl[39])-5.809475019311125*fl[17]-3.354101966249685*fl[8]; 
  incr2[40] = (-7.5*fl[40])-5.809475019311125*fl[18]-3.354101966249685*fl[9]; 
  incr2[42] = 1.5*fl[42]+0.8660254037844387*fl[28]; 
  incr2[43] = 1.5*fl[43]+0.8660254037844387*fl[29]; 
  incr2[44] = 1.5*fl[44]+0.8660254037844387*fl[35]; 
  incr2[45] = 1.5*fl[45]+0.8660254037844387*fl[36]; 
  incr2[46] = (-7.5*fl[46])-5.809475019311126*fl[31]-3.354101966249684*fl[16]; 
  incr2[47] = 1.5*fl[47]+0.8660254037844387*fl[41]; 

  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul; 
  outl[27] += incr2[27]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 
  outl[32] += -1.0*incr2[32]*rdxFnul; 
  outl[33] += -1.0*incr2[33]*rdxFnul; 
  outl[34] += incr2[34]*rdxFnul; 
  outl[37] += -1.0*incr2[37]*rdxFnul; 
  outl[38] += -1.0*incr2[38]*rdxFnul; 
  outl[39] += incr2[39]*rdxFnul; 
  outl[40] += incr2[40]*rdxFnul; 
  outl[42] += -1.0*incr2[42]*rdxFnul; 
  outl[43] += -1.0*incr2[43]*rdxFnul; 
  outl[44] += -1.0*incr2[44]*rdxFnul; 
  outl[45] += -1.0*incr2[45]*rdxFnul; 
  outl[46] += incr2[46]*rdxFnul; 
  outl[47] += -1.0*incr2[47]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf4xSerP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxFnur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

  double incr1[48]; 
  double incr2[48]; 

  if (idxr[3] == 1) {

  incr2[4] = 1.936491673103708*fr[14]-1.5*fr[4]+0.8660254037844386*fr[0]; 
  incr2[8] = 1.936491673103709*fr[28]-1.5*fr[8]+0.8660254037844386*fr[1]; 
  incr2[9] = 1.936491673103709*fr[29]-1.5*fr[9]+0.8660254037844386*fr[2]; 
  incr2[10] = 1.936491673103709*fr[30]-1.5*fr[10]+0.8660254037844386*fr[3]; 
  incr2[14] = (-7.5*fr[14])+5.809475019311125*fr[4]-3.354101966249685*fr[0]; 
  incr2[16] = 1.936491673103708*fr[41]-1.5*fr[16]+0.8660254037844386*fr[5]; 
  incr2[17] = 1.936491673103708*fr[42]-1.5*fr[17]+0.8660254037844386*fr[6]; 
  incr2[18] = 1.936491673103708*fr[43]-1.5*fr[18]+0.8660254037844386*fr[7]; 
  incr2[25] = 0.8660254037844387*fr[11]-1.5*fr[25]; 
  incr2[26] = 0.8660254037844387*fr[12]-1.5*fr[26]; 
  incr2[27] = 0.8660254037844387*fr[13]-1.5*fr[27]; 
  incr2[28] = (-7.5*fr[28])+5.809475019311126*fr[8]-3.354101966249684*fr[1]; 
  incr2[29] = (-7.5*fr[29])+5.809475019311126*fr[9]-3.354101966249684*fr[2]; 
  incr2[30] = (-7.5*fr[30])+5.809475019311126*fr[10]-3.354101966249684*fr[3]; 
  incr2[31] = 1.936491673103709*fr[47]-1.5*fr[31]+0.8660254037844386*fr[15]; 
  incr2[35] = 0.8660254037844387*fr[19]-1.5*fr[35]; 
  incr2[36] = 0.8660254037844387*fr[20]-1.5*fr[36]; 
  incr2[37] = 0.8660254037844387*fr[21]-1.5*fr[37]; 
  incr2[38] = 0.8660254037844387*fr[22]-1.5*fr[38]; 
  incr2[39] = 0.8660254037844387*fr[23]-1.5*fr[39]; 
  incr2[40] = 0.8660254037844387*fr[24]-1.5*fr[40]; 
  incr2[41] = (-7.5*fr[41])+5.809475019311125*fr[16]-3.354101966249685*fr[5]; 
  incr2[42] = (-7.5*fr[42])+5.809475019311125*fr[17]-3.354101966249685*fr[6]; 
  incr2[43] = (-7.5*fr[43])+5.809475019311125*fr[18]-3.354101966249685*fr[7]; 
  incr2[44] = 0.8660254037844387*fr[32]-1.5*fr[44]; 
  incr2[45] = 0.8660254037844387*fr[33]-1.5*fr[45]; 
  incr2[46] = 0.8660254037844387*fr[34]-1.5*fr[46]; 
  incr2[47] = (-7.5*fr[47])+5.809475019311126*fr[31]-3.354101966249684*fr[15]; 

  outr[4] += incr2[4]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 
  outr[35] += incr2[35]*rdxFnur; 
  outr[36] += incr2[36]*rdxFnur; 
  outr[37] += incr2[37]*rdxFnur; 
  outr[38] += incr2[38]*rdxFnur; 
  outr[39] += incr2[39]*rdxFnur; 
  outr[40] += incr2[40]*rdxFnur; 
  outr[41] += incr2[41]*rdxFnur; 
  outr[42] += incr2[42]*rdxFnur; 
  outr[43] += incr2[43]*rdxFnur; 
  outr[44] += incr2[44]*rdxFnur; 
  outr[45] += incr2[45]*rdxFnur; 
  outr[46] += incr2[46]*rdxFnur; 
  outr[47] += incr2[47]*rdxFnur; 

  } else {

  incr2[4] = 1.936491673103708*fl[14]+1.5*fl[4]+0.8660254037844386*fl[0]; 
  incr2[8] = 1.936491673103709*fl[28]+1.5*fl[8]+0.8660254037844386*fl[1]; 
  incr2[9] = 1.936491673103709*fl[29]+1.5*fl[9]+0.8660254037844386*fl[2]; 
  incr2[10] = 1.936491673103709*fl[30]+1.5*fl[10]+0.8660254037844386*fl[3]; 
  incr2[14] = (-7.5*fl[14])-5.809475019311125*fl[4]-3.354101966249685*fl[0]; 
  incr2[16] = 1.936491673103708*fl[41]+1.5*fl[16]+0.8660254037844386*fl[5]; 
  incr2[17] = 1.936491673103708*fl[42]+1.5*fl[17]+0.8660254037844386*fl[6]; 
  incr2[18] = 1.936491673103708*fl[43]+1.5*fl[18]+0.8660254037844386*fl[7]; 
  incr2[25] = 1.5*fl[25]+0.8660254037844387*fl[11]; 
  incr2[26] = 1.5*fl[26]+0.8660254037844387*fl[12]; 
  incr2[27] = 1.5*fl[27]+0.8660254037844387*fl[13]; 
  incr2[28] = (-7.5*fl[28])-5.809475019311126*fl[8]-3.354101966249684*fl[1]; 
  incr2[29] = (-7.5*fl[29])-5.809475019311126*fl[9]-3.354101966249684*fl[2]; 
  incr2[30] = (-7.5*fl[30])-5.809475019311126*fl[10]-3.354101966249684*fl[3]; 
  incr2[31] = 1.936491673103709*fl[47]+1.5*fl[31]+0.8660254037844386*fl[15]; 
  incr2[35] = 1.5*fl[35]+0.8660254037844387*fl[19]; 
  incr2[36] = 1.5*fl[36]+0.8660254037844387*fl[20]; 
  incr2[37] = 1.5*fl[37]+0.8660254037844387*fl[21]; 
  incr2[38] = 1.5*fl[38]+0.8660254037844387*fl[22]; 
  incr2[39] = 1.5*fl[39]+0.8660254037844387*fl[23]; 
  incr2[40] = 1.5*fl[40]+0.8660254037844387*fl[24]; 
  incr2[41] = (-7.5*fl[41])-5.809475019311125*fl[16]-3.354101966249685*fl[5]; 
  incr2[42] = (-7.5*fl[42])-5.809475019311125*fl[17]-3.354101966249685*fl[6]; 
  incr2[43] = (-7.5*fl[43])-5.809475019311125*fl[18]-3.354101966249685*fl[7]; 
  incr2[44] = 1.5*fl[44]+0.8660254037844387*fl[32]; 
  incr2[45] = 1.5*fl[45]+0.8660254037844387*fl[33]; 
  incr2[46] = 1.5*fl[46]+0.8660254037844387*fl[34]; 
  incr2[47] = (-7.5*fl[47])-5.809475019311126*fl[31]-3.354101966249684*fl[15]; 

  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[28] += incr2[28]*rdxFnul; 
  outl[29] += incr2[29]*rdxFnul; 
  outl[30] += incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 
  outl[35] += -1.0*incr2[35]*rdxFnul; 
  outl[36] += -1.0*incr2[36]*rdxFnul; 
  outl[37] += -1.0*incr2[37]*rdxFnul; 
  outl[38] += -1.0*incr2[38]*rdxFnul; 
  outl[39] += -1.0*incr2[39]*rdxFnul; 
  outl[40] += -1.0*incr2[40]*rdxFnul; 
  outl[41] += incr2[41]*rdxFnul; 
  outl[42] += incr2[42]*rdxFnul; 
  outl[43] += incr2[43]*rdxFnul; 
  outl[44] += -1.0*incr2[44]*rdxFnul; 
  outl[45] += -1.0*incr2[45]*rdxFnul; 
  outl[46] += -1.0*incr2[46]*rdxFnul; 
  outl[47] += incr2[47]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf4xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 

  if (idxr[0] == 1) {


  incr2[1] = 2.8125*fr[1]-5.083290641897234*fr[11]; 
  incr2[5] = 2.8125*fr[5]-5.083290641897235*fr[19]; 
  incr2[6] = 2.8125*fr[6]-5.083290641897235*fr[21]; 
  incr2[8] = 2.8125*fr[8]-5.083290641897235*fr[25]; 
  incr2[11] = 19.6875*fr[11]-10.89276566120836*fr[1]; 
  incr2[15] = 2.8125*fr[15]-5.083290641897234*fr[32]; 
  incr2[16] = 2.8125*fr[16]-5.083290641897234*fr[35]; 
  incr2[17] = 2.8125*fr[17]-5.083290641897234*fr[37]; 
  incr2[19] = 19.6875*fr[19]-10.89276566120836*fr[5]; 
  incr2[20] = 2.8125*fr[20]; 
  incr2[21] = 19.6875*fr[21]-10.89276566120836*fr[6]; 
  incr2[23] = 2.8125*fr[23]; 
  incr2[25] = 19.6875*fr[25]-10.89276566120836*fr[8]; 
  incr2[28] = 2.8125*fr[28]; 
  incr2[31] = 2.8125*fr[31]-5.083290641897235*fr[44]; 
  incr2[32] = 19.6875*fr[32]-10.89276566120836*fr[15]; 
  incr2[33] = 2.8125*fr[33]; 
  incr2[34] = 2.8125*fr[34]; 
  incr2[35] = 19.6875*fr[35]-10.89276566120836*fr[16]; 
  incr2[36] = 2.8125*fr[36]; 
  incr2[37] = 19.6875*fr[37]-10.89276566120836*fr[17]; 
  incr2[39] = 2.8125*fr[39]; 
  incr2[41] = 2.8125*fr[41]; 
  incr2[42] = 2.8125*fr[42]; 
  incr2[44] = 19.6875*fr[44]-10.89276566120836*fr[31]; 
  incr2[45] = 2.8125*fr[45]; 
  incr2[46] = 2.8125*fr[46]; 
  incr2[47] = 2.8125*fr[47]; 



  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[8] += -1.0*incr2[8]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[19] += -1.0*incr2[19]*rdxFnur; 
  outr[20] += -1.0*incr2[20]*rdxFnur; 
  outr[21] += -1.0*incr2[21]*rdxFnur; 
  outr[23] += -1.0*incr2[23]*rdxFnur; 
  outr[25] += -1.0*incr2[25]*rdxFnur; 
  outr[28] += -1.0*incr2[28]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 
  outr[32] += -1.0*incr2[32]*rdxFnur; 
  outr[33] += -1.0*incr2[33]*rdxFnur; 
  outr[34] += -1.0*incr2[34]*rdxFnur; 
  outr[35] += -1.0*incr2[35]*rdxFnur; 
  outr[36] += -1.0*incr2[36]*rdxFnur; 
  outr[37] += -1.0*incr2[37]*rdxFnur; 
  outr[39] += -1.0*incr2[39]*rdxFnur; 
  outr[41] += -1.0*incr2[41]*rdxFnur; 
  outr[42] += -1.0*incr2[42]*rdxFnur; 
  outr[44] += -1.0*incr2[44]*rdxFnur; 
  outr[45] += -1.0*incr2[45]*rdxFnur; 
  outr[46] += -1.0*incr2[46]*rdxFnur; 
  outr[47] += -1.0*incr2[47]*rdxFnur; 

  } else {


  incr2[1] = 5.083290641897234*fl[11]+2.8125*fl[1]; 
  incr2[5] = 5.083290641897235*fl[19]+2.8125*fl[5]; 
  incr2[6] = 5.083290641897235*fl[21]+2.8125*fl[6]; 
  incr2[8] = 5.083290641897235*fl[25]+2.8125*fl[8]; 
  incr2[11] = 19.6875*fl[11]+10.89276566120836*fl[1]; 
  incr2[15] = 5.083290641897234*fl[32]+2.8125*fl[15]; 
  incr2[16] = 5.083290641897234*fl[35]+2.8125*fl[16]; 
  incr2[17] = 5.083290641897234*fl[37]+2.8125*fl[17]; 
  incr2[19] = 19.6875*fl[19]+10.89276566120836*fl[5]; 
  incr2[20] = 2.8125*fl[20]; 
  incr2[21] = 19.6875*fl[21]+10.89276566120836*fl[6]; 
  incr2[23] = 2.8125*fl[23]; 
  incr2[25] = 19.6875*fl[25]+10.89276566120836*fl[8]; 
  incr2[28] = 2.8125*fl[28]; 
  incr2[31] = 5.083290641897235*fl[44]+2.8125*fl[31]; 
  incr2[32] = 19.6875*fl[32]+10.89276566120836*fl[15]; 
  incr2[33] = 2.8125*fl[33]; 
  incr2[34] = 2.8125*fl[34]; 
  incr2[35] = 19.6875*fl[35]+10.89276566120836*fl[16]; 
  incr2[36] = 2.8125*fl[36]; 
  incr2[37] = 19.6875*fl[37]+10.89276566120836*fl[17]; 
  incr2[39] = 2.8125*fl[39]; 
  incr2[41] = 2.8125*fl[41]; 
  incr2[42] = 2.8125*fl[42]; 
  incr2[44] = 19.6875*fl[44]+10.89276566120836*fl[31]; 
  incr2[45] = 2.8125*fl[45]; 
  incr2[46] = 2.8125*fl[46]; 
  incr2[47] = 2.8125*fl[47]; 



  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 
  outl[32] += -1.0*incr2[32]*rdxFnul; 
  outl[33] += -1.0*incr2[33]*rdxFnul; 
  outl[34] += -1.0*incr2[34]*rdxFnul; 
  outl[35] += -1.0*incr2[35]*rdxFnul; 
  outl[36] += -1.0*incr2[36]*rdxFnul; 
  outl[37] += -1.0*incr2[37]*rdxFnul; 
  outl[39] += -1.0*incr2[39]*rdxFnul; 
  outl[41] += -1.0*incr2[41]*rdxFnul; 
  outl[42] += -1.0*incr2[42]*rdxFnul; 
  outl[44] += -1.0*incr2[44]*rdxFnul; 
  outl[45] += -1.0*incr2[45]*rdxFnul; 
  outl[46] += -1.0*incr2[46]*rdxFnul; 
  outl[47] += -1.0*incr2[47]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf4xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 

  if (idxr[1] == 1) {


  incr2[2] = 2.8125*fr[2]-5.083290641897234*fr[12]; 
  incr2[5] = 2.8125*fr[5]-5.083290641897235*fr[20]; 
  incr2[7] = 2.8125*fr[7]-5.083290641897235*fr[22]; 
  incr2[9] = 2.8125*fr[9]-5.083290641897235*fr[26]; 
  incr2[12] = 19.6875*fr[12]-10.89276566120836*fr[2]; 
  incr2[15] = 2.8125*fr[15]-5.083290641897234*fr[33]; 
  incr2[16] = 2.8125*fr[16]-5.083290641897234*fr[36]; 
  incr2[18] = 2.8125*fr[18]-5.083290641897234*fr[38]; 
  incr2[19] = 2.8125*fr[19]; 
  incr2[20] = 19.6875*fr[20]-10.89276566120836*fr[5]; 
  incr2[22] = 19.6875*fr[22]-10.89276566120836*fr[7]; 
  incr2[24] = 2.8125*fr[24]; 
  incr2[26] = 19.6875*fr[26]-10.89276566120836*fr[9]; 
  incr2[29] = 2.8125*fr[29]; 
  incr2[31] = 2.8125*fr[31]-5.083290641897235*fr[45]; 
  incr2[32] = 2.8125*fr[32]; 
  incr2[33] = 19.6875*fr[33]-10.89276566120836*fr[15]; 
  incr2[34] = 2.8125*fr[34]; 
  incr2[35] = 2.8125*fr[35]; 
  incr2[36] = 19.6875*fr[36]-10.89276566120836*fr[16]; 
  incr2[38] = 19.6875*fr[38]-10.89276566120836*fr[18]; 
  incr2[40] = 2.8125*fr[40]; 
  incr2[41] = 2.8125*fr[41]; 
  incr2[43] = 2.8125*fr[43]; 
  incr2[44] = 2.8125*fr[44]; 
  incr2[45] = 19.6875*fr[45]-10.89276566120836*fr[31]; 
  incr2[46] = 2.8125*fr[46]; 
  incr2[47] = 2.8125*fr[47]; 



  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[9] += -1.0*incr2[9]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[19] += -1.0*incr2[19]*rdxFnur; 
  outr[20] += -1.0*incr2[20]*rdxFnur; 
  outr[22] += -1.0*incr2[22]*rdxFnur; 
  outr[24] += -1.0*incr2[24]*rdxFnur; 
  outr[26] += -1.0*incr2[26]*rdxFnur; 
  outr[29] += -1.0*incr2[29]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 
  outr[32] += -1.0*incr2[32]*rdxFnur; 
  outr[33] += -1.0*incr2[33]*rdxFnur; 
  outr[34] += -1.0*incr2[34]*rdxFnur; 
  outr[35] += -1.0*incr2[35]*rdxFnur; 
  outr[36] += -1.0*incr2[36]*rdxFnur; 
  outr[38] += -1.0*incr2[38]*rdxFnur; 
  outr[40] += -1.0*incr2[40]*rdxFnur; 
  outr[41] += -1.0*incr2[41]*rdxFnur; 
  outr[43] += -1.0*incr2[43]*rdxFnur; 
  outr[44] += -1.0*incr2[44]*rdxFnur; 
  outr[45] += -1.0*incr2[45]*rdxFnur; 
  outr[46] += -1.0*incr2[46]*rdxFnur; 
  outr[47] += -1.0*incr2[47]*rdxFnur; 

  } else {


  incr2[2] = 5.083290641897234*fl[12]+2.8125*fl[2]; 
  incr2[5] = 5.083290641897235*fl[20]+2.8125*fl[5]; 
  incr2[7] = 5.083290641897235*fl[22]+2.8125*fl[7]; 
  incr2[9] = 5.083290641897235*fl[26]+2.8125*fl[9]; 
  incr2[12] = 19.6875*fl[12]+10.89276566120836*fl[2]; 
  incr2[15] = 5.083290641897234*fl[33]+2.8125*fl[15]; 
  incr2[16] = 5.083290641897234*fl[36]+2.8125*fl[16]; 
  incr2[18] = 5.083290641897234*fl[38]+2.8125*fl[18]; 
  incr2[19] = 2.8125*fl[19]; 
  incr2[20] = 19.6875*fl[20]+10.89276566120836*fl[5]; 
  incr2[22] = 19.6875*fl[22]+10.89276566120836*fl[7]; 
  incr2[24] = 2.8125*fl[24]; 
  incr2[26] = 19.6875*fl[26]+10.89276566120836*fl[9]; 
  incr2[29] = 2.8125*fl[29]; 
  incr2[31] = 5.083290641897235*fl[45]+2.8125*fl[31]; 
  incr2[32] = 2.8125*fl[32]; 
  incr2[33] = 19.6875*fl[33]+10.89276566120836*fl[15]; 
  incr2[34] = 2.8125*fl[34]; 
  incr2[35] = 2.8125*fl[35]; 
  incr2[36] = 19.6875*fl[36]+10.89276566120836*fl[16]; 
  incr2[38] = 19.6875*fl[38]+10.89276566120836*fl[18]; 
  incr2[40] = 2.8125*fl[40]; 
  incr2[41] = 2.8125*fl[41]; 
  incr2[43] = 2.8125*fl[43]; 
  incr2[44] = 2.8125*fl[44]; 
  incr2[45] = 19.6875*fl[45]+10.89276566120836*fl[31]; 
  incr2[46] = 2.8125*fl[46]; 
  incr2[47] = 2.8125*fl[47]; 



  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 
  outl[32] += -1.0*incr2[32]*rdxFnul; 
  outl[33] += -1.0*incr2[33]*rdxFnul; 
  outl[34] += -1.0*incr2[34]*rdxFnul; 
  outl[35] += -1.0*incr2[35]*rdxFnul; 
  outl[36] += -1.0*incr2[36]*rdxFnul; 
  outl[38] += -1.0*incr2[38]*rdxFnul; 
  outl[40] += -1.0*incr2[40]*rdxFnul; 
  outl[41] += -1.0*incr2[41]*rdxFnul; 
  outl[43] += -1.0*incr2[43]*rdxFnul; 
  outl[44] += -1.0*incr2[44]*rdxFnul; 
  outl[45] += -1.0*incr2[45]*rdxFnul; 
  outl[46] += -1.0*incr2[46]*rdxFnul; 
  outl[47] += -1.0*incr2[47]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf4xSerP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 

  if (idxr[2] == 1) {


  incr2[3] = 2.8125*fr[3]-5.083290641897234*fr[13]; 
  incr2[6] = 2.8125*fr[6]-5.083290641897235*fr[23]; 
  incr2[7] = 2.8125*fr[7]-5.083290641897235*fr[24]; 
  incr2[10] = 2.8125*fr[10]-5.083290641897235*fr[27]; 
  incr2[13] = 19.6875*fr[13]-10.89276566120836*fr[3]; 
  incr2[15] = 2.8125*fr[15]-5.083290641897234*fr[34]; 
  incr2[17] = 2.8125*fr[17]-5.083290641897234*fr[39]; 
  incr2[18] = 2.8125*fr[18]-5.083290641897234*fr[40]; 
  incr2[21] = 2.8125*fr[21]; 
  incr2[22] = 2.8125*fr[22]; 
  incr2[23] = 19.6875*fr[23]-10.89276566120836*fr[6]; 
  incr2[24] = 19.6875*fr[24]-10.89276566120836*fr[7]; 
  incr2[27] = 19.6875*fr[27]-10.89276566120836*fr[10]; 
  incr2[30] = 2.8125*fr[30]; 
  incr2[31] = 2.8125*fr[31]-5.083290641897235*fr[46]; 
  incr2[32] = 2.8125*fr[32]; 
  incr2[33] = 2.8125*fr[33]; 
  incr2[34] = 19.6875*fr[34]-10.89276566120836*fr[15]; 
  incr2[37] = 2.8125*fr[37]; 
  incr2[38] = 2.8125*fr[38]; 
  incr2[39] = 19.6875*fr[39]-10.89276566120836*fr[17]; 
  incr2[40] = 19.6875*fr[40]-10.89276566120836*fr[18]; 
  incr2[42] = 2.8125*fr[42]; 
  incr2[43] = 2.8125*fr[43]; 
  incr2[44] = 2.8125*fr[44]; 
  incr2[45] = 2.8125*fr[45]; 
  incr2[46] = 19.6875*fr[46]-10.89276566120836*fr[31]; 
  incr2[47] = 2.8125*fr[47]; 



  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[21] += -1.0*incr2[21]*rdxFnur; 
  outr[22] += -1.0*incr2[22]*rdxFnur; 
  outr[23] += -1.0*incr2[23]*rdxFnur; 
  outr[24] += -1.0*incr2[24]*rdxFnur; 
  outr[27] += -1.0*incr2[27]*rdxFnur; 
  outr[30] += -1.0*incr2[30]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 
  outr[32] += -1.0*incr2[32]*rdxFnur; 
  outr[33] += -1.0*incr2[33]*rdxFnur; 
  outr[34] += -1.0*incr2[34]*rdxFnur; 
  outr[37] += -1.0*incr2[37]*rdxFnur; 
  outr[38] += -1.0*incr2[38]*rdxFnur; 
  outr[39] += -1.0*incr2[39]*rdxFnur; 
  outr[40] += -1.0*incr2[40]*rdxFnur; 
  outr[42] += -1.0*incr2[42]*rdxFnur; 
  outr[43] += -1.0*incr2[43]*rdxFnur; 
  outr[44] += -1.0*incr2[44]*rdxFnur; 
  outr[45] += -1.0*incr2[45]*rdxFnur; 
  outr[46] += -1.0*incr2[46]*rdxFnur; 
  outr[47] += -1.0*incr2[47]*rdxFnur; 

  } else {


  incr2[3] = 5.083290641897234*fl[13]+2.8125*fl[3]; 
  incr2[6] = 5.083290641897235*fl[23]+2.8125*fl[6]; 
  incr2[7] = 5.083290641897235*fl[24]+2.8125*fl[7]; 
  incr2[10] = 5.083290641897235*fl[27]+2.8125*fl[10]; 
  incr2[13] = 19.6875*fl[13]+10.89276566120836*fl[3]; 
  incr2[15] = 5.083290641897234*fl[34]+2.8125*fl[15]; 
  incr2[17] = 5.083290641897234*fl[39]+2.8125*fl[17]; 
  incr2[18] = 5.083290641897234*fl[40]+2.8125*fl[18]; 
  incr2[21] = 2.8125*fl[21]; 
  incr2[22] = 2.8125*fl[22]; 
  incr2[23] = 19.6875*fl[23]+10.89276566120836*fl[6]; 
  incr2[24] = 19.6875*fl[24]+10.89276566120836*fl[7]; 
  incr2[27] = 19.6875*fl[27]+10.89276566120836*fl[10]; 
  incr2[30] = 2.8125*fl[30]; 
  incr2[31] = 5.083290641897235*fl[46]+2.8125*fl[31]; 
  incr2[32] = 2.8125*fl[32]; 
  incr2[33] = 2.8125*fl[33]; 
  incr2[34] = 19.6875*fl[34]+10.89276566120836*fl[15]; 
  incr2[37] = 2.8125*fl[37]; 
  incr2[38] = 2.8125*fl[38]; 
  incr2[39] = 19.6875*fl[39]+10.89276566120836*fl[17]; 
  incr2[40] = 19.6875*fl[40]+10.89276566120836*fl[18]; 
  incr2[42] = 2.8125*fl[42]; 
  incr2[43] = 2.8125*fl[43]; 
  incr2[44] = 2.8125*fl[44]; 
  incr2[45] = 2.8125*fl[45]; 
  incr2[46] = 19.6875*fl[46]+10.89276566120836*fl[31]; 
  incr2[47] = 2.8125*fl[47]; 



  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 
  outl[32] += -1.0*incr2[32]*rdxFnul; 
  outl[33] += -1.0*incr2[33]*rdxFnul; 
  outl[34] += -1.0*incr2[34]*rdxFnul; 
  outl[37] += -1.0*incr2[37]*rdxFnul; 
  outl[38] += -1.0*incr2[38]*rdxFnul; 
  outl[39] += -1.0*incr2[39]*rdxFnul; 
  outl[40] += -1.0*incr2[40]*rdxFnul; 
  outl[42] += -1.0*incr2[42]*rdxFnul; 
  outl[43] += -1.0*incr2[43]*rdxFnul; 
  outl[44] += -1.0*incr2[44]*rdxFnul; 
  outl[45] += -1.0*incr2[45]*rdxFnul; 
  outl[46] += -1.0*incr2[46]*rdxFnul; 
  outl[47] += -1.0*incr2[47]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf4xSerP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 16.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 

  if (idxr[3] == 1) {


  incr2[4] = 2.8125*fr[4]-5.083290641897234*fr[14]; 
  incr2[8] = 2.8125*fr[8]-5.083290641897235*fr[28]; 
  incr2[9] = 2.8125*fr[9]-5.083290641897235*fr[29]; 
  incr2[10] = 2.8125*fr[10]-5.083290641897235*fr[30]; 
  incr2[14] = 19.6875*fr[14]-10.89276566120836*fr[4]; 
  incr2[16] = 2.8125*fr[16]-5.083290641897234*fr[41]; 
  incr2[17] = 2.8125*fr[17]-5.083290641897234*fr[42]; 
  incr2[18] = 2.8125*fr[18]-5.083290641897234*fr[43]; 
  incr2[25] = 2.8125*fr[25]; 
  incr2[26] = 2.8125*fr[26]; 
  incr2[27] = 2.8125*fr[27]; 
  incr2[28] = 19.6875*fr[28]-10.89276566120836*fr[8]; 
  incr2[29] = 19.6875*fr[29]-10.89276566120836*fr[9]; 
  incr2[30] = 19.6875*fr[30]-10.89276566120836*fr[10]; 
  incr2[31] = 2.8125*fr[31]-5.083290641897235*fr[47]; 
  incr2[35] = 2.8125*fr[35]; 
  incr2[36] = 2.8125*fr[36]; 
  incr2[37] = 2.8125*fr[37]; 
  incr2[38] = 2.8125*fr[38]; 
  incr2[39] = 2.8125*fr[39]; 
  incr2[40] = 2.8125*fr[40]; 
  incr2[41] = 19.6875*fr[41]-10.89276566120836*fr[16]; 
  incr2[42] = 19.6875*fr[42]-10.89276566120836*fr[17]; 
  incr2[43] = 19.6875*fr[43]-10.89276566120836*fr[18]; 
  incr2[44] = 2.8125*fr[44]; 
  incr2[45] = 2.8125*fr[45]; 
  incr2[46] = 2.8125*fr[46]; 
  incr2[47] = 19.6875*fr[47]-10.89276566120836*fr[31]; 



  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[8] += -1.0*incr2[8]*rdxFnur; 
  outr[9] += -1.0*incr2[9]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[25] += -1.0*incr2[25]*rdxFnur; 
  outr[26] += -1.0*incr2[26]*rdxFnur; 
  outr[27] += -1.0*incr2[27]*rdxFnur; 
  outr[28] += -1.0*incr2[28]*rdxFnur; 
  outr[29] += -1.0*incr2[29]*rdxFnur; 
  outr[30] += -1.0*incr2[30]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 
  outr[35] += -1.0*incr2[35]*rdxFnur; 
  outr[36] += -1.0*incr2[36]*rdxFnur; 
  outr[37] += -1.0*incr2[37]*rdxFnur; 
  outr[38] += -1.0*incr2[38]*rdxFnur; 
  outr[39] += -1.0*incr2[39]*rdxFnur; 
  outr[40] += -1.0*incr2[40]*rdxFnur; 
  outr[41] += -1.0*incr2[41]*rdxFnur; 
  outr[42] += -1.0*incr2[42]*rdxFnur; 
  outr[43] += -1.0*incr2[43]*rdxFnur; 
  outr[44] += -1.0*incr2[44]*rdxFnur; 
  outr[45] += -1.0*incr2[45]*rdxFnur; 
  outr[46] += -1.0*incr2[46]*rdxFnur; 
  outr[47] += -1.0*incr2[47]*rdxFnur; 

  } else {


  incr2[4] = 5.083290641897234*fl[14]+2.8125*fl[4]; 
  incr2[8] = 5.083290641897235*fl[28]+2.8125*fl[8]; 
  incr2[9] = 5.083290641897235*fl[29]+2.8125*fl[9]; 
  incr2[10] = 5.083290641897235*fl[30]+2.8125*fl[10]; 
  incr2[14] = 19.6875*fl[14]+10.89276566120836*fl[4]; 
  incr2[16] = 5.083290641897234*fl[41]+2.8125*fl[16]; 
  incr2[17] = 5.083290641897234*fl[42]+2.8125*fl[17]; 
  incr2[18] = 5.083290641897234*fl[43]+2.8125*fl[18]; 
  incr2[25] = 2.8125*fl[25]; 
  incr2[26] = 2.8125*fl[26]; 
  incr2[27] = 2.8125*fl[27]; 
  incr2[28] = 19.6875*fl[28]+10.89276566120836*fl[8]; 
  incr2[29] = 19.6875*fl[29]+10.89276566120836*fl[9]; 
  incr2[30] = 19.6875*fl[30]+10.89276566120836*fl[10]; 
  incr2[31] = 5.083290641897235*fl[47]+2.8125*fl[31]; 
  incr2[35] = 2.8125*fl[35]; 
  incr2[36] = 2.8125*fl[36]; 
  incr2[37] = 2.8125*fl[37]; 
  incr2[38] = 2.8125*fl[38]; 
  incr2[39] = 2.8125*fl[39]; 
  incr2[40] = 2.8125*fl[40]; 
  incr2[41] = 19.6875*fl[41]+10.89276566120836*fl[16]; 
  incr2[42] = 19.6875*fl[42]+10.89276566120836*fl[17]; 
  incr2[43] = 19.6875*fl[43]+10.89276566120836*fl[18]; 
  incr2[44] = 2.8125*fl[44]; 
  incr2[45] = 2.8125*fl[45]; 
  incr2[46] = 2.8125*fl[46]; 
  incr2[47] = 19.6875*fl[47]+10.89276566120836*fl[31]; 



  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 
  outl[35] += -1.0*incr2[35]*rdxFnul; 
  outl[36] += -1.0*incr2[36]*rdxFnul; 
  outl[37] += -1.0*incr2[37]*rdxFnul; 
  outl[38] += -1.0*incr2[38]*rdxFnul; 
  outl[39] += -1.0*incr2[39]*rdxFnul; 
  outl[40] += -1.0*incr2[40]*rdxFnul; 
  outl[41] += -1.0*incr2[41]*rdxFnul; 
  outl[42] += -1.0*incr2[42]*rdxFnul; 
  outl[43] += -1.0*incr2[43]*rdxFnul; 
  outl[44] += -1.0*incr2[44]*rdxFnul; 
  outl[45] += -1.0*incr2[45]*rdxFnul; 
  outl[46] += -1.0*incr2[46]*rdxFnul; 
  outl[47] += -1.0*incr2[47]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf4xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 64.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 
  double incr5[48]; 
  double incr6[48]; 

  if (idxr[0] == 1) {


  incr2[1] = 19.06233990711463*fr[11]-4.921875*fr[1]; 
  incr2[5] = 19.06233990711463*fr[19]-4.921875*fr[5]; 
  incr2[6] = 19.06233990711463*fr[21]-4.921875*fr[6]; 
  incr2[8] = 19.06233990711463*fr[25]-4.921875*fr[8]; 
  incr2[11] = 19.06233990711463*fr[1]-73.828125*fr[11]; 
  incr2[15] = 19.06233990711463*fr[32]-4.921875*fr[15]; 
  incr2[16] = 19.06233990711463*fr[35]-4.921875*fr[16]; 
  incr2[17] = 19.06233990711463*fr[37]-4.921875*fr[17]; 
  incr2[19] = 19.06233990711463*fr[5]-73.828125*fr[19]; 
  incr2[20] = -4.921875*fr[20]; 
  incr2[21] = 19.06233990711463*fr[6]-73.828125*fr[21]; 
  incr2[23] = -4.921875*fr[23]; 
  incr2[25] = 19.06233990711463*fr[8]-73.828125*fr[25]; 
  incr2[28] = -4.921875*fr[28]; 
  incr2[31] = 19.06233990711463*fr[44]-4.921875*fr[31]; 
  incr2[32] = 19.06233990711463*fr[15]-73.828125*fr[32]; 
  incr2[33] = -4.921875*fr[33]; 
  incr2[34] = -4.921875*fr[34]; 
  incr2[35] = 19.06233990711463*fr[16]-73.828125*fr[35]; 
  incr2[36] = -4.921875*fr[36]; 
  incr2[37] = 19.06233990711463*fr[17]-73.828125*fr[37]; 
  incr2[39] = -4.921875*fr[39]; 
  incr2[41] = -4.921875*fr[41]; 
  incr2[42] = -4.921875*fr[42]; 
  incr2[44] = 19.06233990711463*fr[31]-73.828125*fr[44]; 
  incr2[45] = -4.921875*fr[45]; 
  incr2[46] = -4.921875*fr[46]; 
  incr2[47] = -4.921875*fr[47]; 





  outr[1] += incr2[1]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 
  outr[32] += incr2[32]*rdxFnur; 
  outr[33] += incr2[33]*rdxFnur; 
  outr[34] += incr2[34]*rdxFnur; 
  outr[35] += incr2[35]*rdxFnur; 
  outr[36] += incr2[36]*rdxFnur; 
  outr[37] += incr2[37]*rdxFnur; 
  outr[39] += incr2[39]*rdxFnur; 
  outr[41] += incr2[41]*rdxFnur; 
  outr[42] += incr2[42]*rdxFnur; 
  outr[44] += incr2[44]*rdxFnur; 
  outr[45] += incr2[45]*rdxFnur; 
  outr[46] += incr2[46]*rdxFnur; 
  outr[47] += incr2[47]*rdxFnur; 

  } else {


  incr2[1] = (-19.06233990711463*fl[11])-4.921875*fl[1]; 
  incr2[5] = (-19.06233990711463*fl[19])-4.921875*fl[5]; 
  incr2[6] = (-19.06233990711463*fl[21])-4.921875*fl[6]; 
  incr2[8] = (-19.06233990711463*fl[25])-4.921875*fl[8]; 
  incr2[11] = (-73.828125*fl[11])-19.06233990711463*fl[1]; 
  incr2[15] = (-19.06233990711463*fl[32])-4.921875*fl[15]; 
  incr2[16] = (-19.06233990711463*fl[35])-4.921875*fl[16]; 
  incr2[17] = (-19.06233990711463*fl[37])-4.921875*fl[17]; 
  incr2[19] = (-73.828125*fl[19])-19.06233990711463*fl[5]; 
  incr2[20] = -4.921875*fl[20]; 
  incr2[21] = (-73.828125*fl[21])-19.06233990711463*fl[6]; 
  incr2[23] = -4.921875*fl[23]; 
  incr2[25] = (-73.828125*fl[25])-19.06233990711463*fl[8]; 
  incr2[28] = -4.921875*fl[28]; 
  incr2[31] = (-19.06233990711463*fl[44])-4.921875*fl[31]; 
  incr2[32] = (-73.828125*fl[32])-19.06233990711463*fl[15]; 
  incr2[33] = -4.921875*fl[33]; 
  incr2[34] = -4.921875*fl[34]; 
  incr2[35] = (-73.828125*fl[35])-19.06233990711463*fl[16]; 
  incr2[36] = -4.921875*fl[36]; 
  incr2[37] = (-73.828125*fl[37])-19.06233990711463*fl[17]; 
  incr2[39] = -4.921875*fl[39]; 
  incr2[41] = -4.921875*fl[41]; 
  incr2[42] = -4.921875*fl[42]; 
  incr2[44] = (-73.828125*fl[44])-19.06233990711463*fl[31]; 
  incr2[45] = -4.921875*fl[45]; 
  incr2[46] = -4.921875*fl[46]; 
  incr2[47] = -4.921875*fl[47]; 





  outl[1] += incr2[1]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul; 
  outl[28] += incr2[28]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul; 
  outl[32] += incr2[32]*rdxFnul; 
  outl[33] += incr2[33]*rdxFnul; 
  outl[34] += incr2[34]*rdxFnul; 
  outl[35] += incr2[35]*rdxFnul; 
  outl[36] += incr2[36]*rdxFnul; 
  outl[37] += incr2[37]*rdxFnul; 
  outl[39] += incr2[39]*rdxFnul; 
  outl[41] += incr2[41]*rdxFnul; 
  outl[42] += incr2[42]*rdxFnul; 
  outl[44] += incr2[44]*rdxFnul; 
  outl[45] += incr2[45]*rdxFnul; 
  outl[46] += incr2[46]*rdxFnul; 
  outl[47] += incr2[47]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf4xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 64.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 
  double incr5[48]; 
  double incr6[48]; 

  if (idxr[1] == 1) {


  incr2[2] = 19.06233990711463*fr[12]-4.921875*fr[2]; 
  incr2[5] = 19.06233990711463*fr[20]-4.921875*fr[5]; 
  incr2[7] = 19.06233990711463*fr[22]-4.921875*fr[7]; 
  incr2[9] = 19.06233990711463*fr[26]-4.921875*fr[9]; 
  incr2[12] = 19.06233990711463*fr[2]-73.828125*fr[12]; 
  incr2[15] = 19.06233990711463*fr[33]-4.921875*fr[15]; 
  incr2[16] = 19.06233990711463*fr[36]-4.921875*fr[16]; 
  incr2[18] = 19.06233990711463*fr[38]-4.921875*fr[18]; 
  incr2[19] = -4.921875*fr[19]; 
  incr2[20] = 19.06233990711463*fr[5]-73.828125*fr[20]; 
  incr2[22] = 19.06233990711463*fr[7]-73.828125*fr[22]; 
  incr2[24] = -4.921875*fr[24]; 
  incr2[26] = 19.06233990711463*fr[9]-73.828125*fr[26]; 
  incr2[29] = -4.921875*fr[29]; 
  incr2[31] = 19.06233990711463*fr[45]-4.921875*fr[31]; 
  incr2[32] = -4.921875*fr[32]; 
  incr2[33] = 19.06233990711463*fr[15]-73.828125*fr[33]; 
  incr2[34] = -4.921875*fr[34]; 
  incr2[35] = -4.921875*fr[35]; 
  incr2[36] = 19.06233990711463*fr[16]-73.828125*fr[36]; 
  incr2[38] = 19.06233990711463*fr[18]-73.828125*fr[38]; 
  incr2[40] = -4.921875*fr[40]; 
  incr2[41] = -4.921875*fr[41]; 
  incr2[43] = -4.921875*fr[43]; 
  incr2[44] = -4.921875*fr[44]; 
  incr2[45] = 19.06233990711463*fr[31]-73.828125*fr[45]; 
  incr2[46] = -4.921875*fr[46]; 
  incr2[47] = -4.921875*fr[47]; 





  outr[2] += incr2[2]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 
  outr[32] += incr2[32]*rdxFnur; 
  outr[33] += incr2[33]*rdxFnur; 
  outr[34] += incr2[34]*rdxFnur; 
  outr[35] += incr2[35]*rdxFnur; 
  outr[36] += incr2[36]*rdxFnur; 
  outr[38] += incr2[38]*rdxFnur; 
  outr[40] += incr2[40]*rdxFnur; 
  outr[41] += incr2[41]*rdxFnur; 
  outr[43] += incr2[43]*rdxFnur; 
  outr[44] += incr2[44]*rdxFnur; 
  outr[45] += incr2[45]*rdxFnur; 
  outr[46] += incr2[46]*rdxFnur; 
  outr[47] += incr2[47]*rdxFnur; 

  } else {


  incr2[2] = (-19.06233990711463*fl[12])-4.921875*fl[2]; 
  incr2[5] = (-19.06233990711463*fl[20])-4.921875*fl[5]; 
  incr2[7] = (-19.06233990711463*fl[22])-4.921875*fl[7]; 
  incr2[9] = (-19.06233990711463*fl[26])-4.921875*fl[9]; 
  incr2[12] = (-73.828125*fl[12])-19.06233990711463*fl[2]; 
  incr2[15] = (-19.06233990711463*fl[33])-4.921875*fl[15]; 
  incr2[16] = (-19.06233990711463*fl[36])-4.921875*fl[16]; 
  incr2[18] = (-19.06233990711463*fl[38])-4.921875*fl[18]; 
  incr2[19] = -4.921875*fl[19]; 
  incr2[20] = (-73.828125*fl[20])-19.06233990711463*fl[5]; 
  incr2[22] = (-73.828125*fl[22])-19.06233990711463*fl[7]; 
  incr2[24] = -4.921875*fl[24]; 
  incr2[26] = (-73.828125*fl[26])-19.06233990711463*fl[9]; 
  incr2[29] = -4.921875*fl[29]; 
  incr2[31] = (-19.06233990711463*fl[45])-4.921875*fl[31]; 
  incr2[32] = -4.921875*fl[32]; 
  incr2[33] = (-73.828125*fl[33])-19.06233990711463*fl[15]; 
  incr2[34] = -4.921875*fl[34]; 
  incr2[35] = -4.921875*fl[35]; 
  incr2[36] = (-73.828125*fl[36])-19.06233990711463*fl[16]; 
  incr2[38] = (-73.828125*fl[38])-19.06233990711463*fl[18]; 
  incr2[40] = -4.921875*fl[40]; 
  incr2[41] = -4.921875*fl[41]; 
  incr2[43] = -4.921875*fl[43]; 
  incr2[44] = -4.921875*fl[44]; 
  incr2[45] = (-73.828125*fl[45])-19.06233990711463*fl[31]; 
  incr2[46] = -4.921875*fl[46]; 
  incr2[47] = -4.921875*fl[47]; 





  outl[2] += incr2[2]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul; 
  outl[29] += incr2[29]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul; 
  outl[32] += incr2[32]*rdxFnul; 
  outl[33] += incr2[33]*rdxFnul; 
  outl[34] += incr2[34]*rdxFnul; 
  outl[35] += incr2[35]*rdxFnul; 
  outl[36] += incr2[36]*rdxFnul; 
  outl[38] += incr2[38]*rdxFnul; 
  outl[40] += incr2[40]*rdxFnul; 
  outl[41] += incr2[41]*rdxFnul; 
  outl[43] += incr2[43]*rdxFnul; 
  outl[44] += incr2[44]*rdxFnul; 
  outl[45] += incr2[45]*rdxFnul; 
  outl[46] += incr2[46]*rdxFnul; 
  outl[47] += incr2[47]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf4xSerP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 64.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 
  double incr5[48]; 
  double incr6[48]; 

  if (idxr[2] == 1) {


  incr2[3] = 19.06233990711463*fr[13]-4.921875*fr[3]; 
  incr2[6] = 19.06233990711463*fr[23]-4.921875*fr[6]; 
  incr2[7] = 19.06233990711463*fr[24]-4.921875*fr[7]; 
  incr2[10] = 19.06233990711463*fr[27]-4.921875*fr[10]; 
  incr2[13] = 19.06233990711463*fr[3]-73.828125*fr[13]; 
  incr2[15] = 19.06233990711463*fr[34]-4.921875*fr[15]; 
  incr2[17] = 19.06233990711463*fr[39]-4.921875*fr[17]; 
  incr2[18] = 19.06233990711463*fr[40]-4.921875*fr[18]; 
  incr2[21] = -4.921875*fr[21]; 
  incr2[22] = -4.921875*fr[22]; 
  incr2[23] = 19.06233990711463*fr[6]-73.828125*fr[23]; 
  incr2[24] = 19.06233990711463*fr[7]-73.828125*fr[24]; 
  incr2[27] = 19.06233990711463*fr[10]-73.828125*fr[27]; 
  incr2[30] = -4.921875*fr[30]; 
  incr2[31] = 19.06233990711463*fr[46]-4.921875*fr[31]; 
  incr2[32] = -4.921875*fr[32]; 
  incr2[33] = -4.921875*fr[33]; 
  incr2[34] = 19.06233990711463*fr[15]-73.828125*fr[34]; 
  incr2[37] = -4.921875*fr[37]; 
  incr2[38] = -4.921875*fr[38]; 
  incr2[39] = 19.06233990711463*fr[17]-73.828125*fr[39]; 
  incr2[40] = 19.06233990711463*fr[18]-73.828125*fr[40]; 
  incr2[42] = -4.921875*fr[42]; 
  incr2[43] = -4.921875*fr[43]; 
  incr2[44] = -4.921875*fr[44]; 
  incr2[45] = -4.921875*fr[45]; 
  incr2[46] = 19.06233990711463*fr[31]-73.828125*fr[46]; 
  incr2[47] = -4.921875*fr[47]; 





  outr[3] += incr2[3]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 
  outr[32] += incr2[32]*rdxFnur; 
  outr[33] += incr2[33]*rdxFnur; 
  outr[34] += incr2[34]*rdxFnur; 
  outr[37] += incr2[37]*rdxFnur; 
  outr[38] += incr2[38]*rdxFnur; 
  outr[39] += incr2[39]*rdxFnur; 
  outr[40] += incr2[40]*rdxFnur; 
  outr[42] += incr2[42]*rdxFnur; 
  outr[43] += incr2[43]*rdxFnur; 
  outr[44] += incr2[44]*rdxFnur; 
  outr[45] += incr2[45]*rdxFnur; 
  outr[46] += incr2[46]*rdxFnur; 
  outr[47] += incr2[47]*rdxFnur; 

  } else {


  incr2[3] = (-19.06233990711463*fl[13])-4.921875*fl[3]; 
  incr2[6] = (-19.06233990711463*fl[23])-4.921875*fl[6]; 
  incr2[7] = (-19.06233990711463*fl[24])-4.921875*fl[7]; 
  incr2[10] = (-19.06233990711463*fl[27])-4.921875*fl[10]; 
  incr2[13] = (-73.828125*fl[13])-19.06233990711463*fl[3]; 
  incr2[15] = (-19.06233990711463*fl[34])-4.921875*fl[15]; 
  incr2[17] = (-19.06233990711463*fl[39])-4.921875*fl[17]; 
  incr2[18] = (-19.06233990711463*fl[40])-4.921875*fl[18]; 
  incr2[21] = -4.921875*fl[21]; 
  incr2[22] = -4.921875*fl[22]; 
  incr2[23] = (-73.828125*fl[23])-19.06233990711463*fl[6]; 
  incr2[24] = (-73.828125*fl[24])-19.06233990711463*fl[7]; 
  incr2[27] = (-73.828125*fl[27])-19.06233990711463*fl[10]; 
  incr2[30] = -4.921875*fl[30]; 
  incr2[31] = (-19.06233990711463*fl[46])-4.921875*fl[31]; 
  incr2[32] = -4.921875*fl[32]; 
  incr2[33] = -4.921875*fl[33]; 
  incr2[34] = (-73.828125*fl[34])-19.06233990711463*fl[15]; 
  incr2[37] = -4.921875*fl[37]; 
  incr2[38] = -4.921875*fl[38]; 
  incr2[39] = (-73.828125*fl[39])-19.06233990711463*fl[17]; 
  incr2[40] = (-73.828125*fl[40])-19.06233990711463*fl[18]; 
  incr2[42] = -4.921875*fl[42]; 
  incr2[43] = -4.921875*fl[43]; 
  incr2[44] = -4.921875*fl[44]; 
  incr2[45] = -4.921875*fl[45]; 
  incr2[46] = (-73.828125*fl[46])-19.06233990711463*fl[31]; 
  incr2[47] = -4.921875*fl[47]; 





  outl[3] += incr2[3]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul; 
  outl[27] += incr2[27]*rdxFnul; 
  outl[30] += incr2[30]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul; 
  outl[32] += incr2[32]*rdxFnul; 
  outl[33] += incr2[33]*rdxFnul; 
  outl[34] += incr2[34]*rdxFnul; 
  outl[37] += incr2[37]*rdxFnul; 
  outl[38] += incr2[38]*rdxFnul; 
  outl[39] += incr2[39]*rdxFnul; 
  outl[40] += incr2[40]*rdxFnul; 
  outl[42] += incr2[42]*rdxFnul; 
  outl[43] += incr2[43]*rdxFnul; 
  outl[44] += incr2[44]*rdxFnul; 
  outl[45] += incr2[45]*rdxFnul; 
  outl[46] += incr2[46]*rdxFnul; 
  outl[47] += incr2[47]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf4xSerP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[4]:    current grid index.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 64.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 64.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 
  double incr5[48]; 
  double incr6[48]; 

  if (idxr[3] == 1) {


  incr2[4] = 19.06233990711463*fr[14]-4.921875*fr[4]; 
  incr2[8] = 19.06233990711463*fr[28]-4.921875*fr[8]; 
  incr2[9] = 19.06233990711463*fr[29]-4.921875*fr[9]; 
  incr2[10] = 19.06233990711463*fr[30]-4.921875*fr[10]; 
  incr2[14] = 19.06233990711463*fr[4]-73.828125*fr[14]; 
  incr2[16] = 19.06233990711463*fr[41]-4.921875*fr[16]; 
  incr2[17] = 19.06233990711463*fr[42]-4.921875*fr[17]; 
  incr2[18] = 19.06233990711463*fr[43]-4.921875*fr[18]; 
  incr2[25] = -4.921875*fr[25]; 
  incr2[26] = -4.921875*fr[26]; 
  incr2[27] = -4.921875*fr[27]; 
  incr2[28] = 19.06233990711463*fr[8]-73.828125*fr[28]; 
  incr2[29] = 19.06233990711463*fr[9]-73.828125*fr[29]; 
  incr2[30] = 19.06233990711463*fr[10]-73.828125*fr[30]; 
  incr2[31] = 19.06233990711463*fr[47]-4.921875*fr[31]; 
  incr2[35] = -4.921875*fr[35]; 
  incr2[36] = -4.921875*fr[36]; 
  incr2[37] = -4.921875*fr[37]; 
  incr2[38] = -4.921875*fr[38]; 
  incr2[39] = -4.921875*fr[39]; 
  incr2[40] = -4.921875*fr[40]; 
  incr2[41] = 19.06233990711463*fr[16]-73.828125*fr[41]; 
  incr2[42] = 19.06233990711463*fr[17]-73.828125*fr[42]; 
  incr2[43] = 19.06233990711463*fr[18]-73.828125*fr[43]; 
  incr2[44] = -4.921875*fr[44]; 
  incr2[45] = -4.921875*fr[45]; 
  incr2[46] = -4.921875*fr[46]; 
  incr2[47] = 19.06233990711463*fr[31]-73.828125*fr[47]; 





  outr[4] += incr2[4]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 
  outr[35] += incr2[35]*rdxFnur; 
  outr[36] += incr2[36]*rdxFnur; 
  outr[37] += incr2[37]*rdxFnur; 
  outr[38] += incr2[38]*rdxFnur; 
  outr[39] += incr2[39]*rdxFnur; 
  outr[40] += incr2[40]*rdxFnur; 
  outr[41] += incr2[41]*rdxFnur; 
  outr[42] += incr2[42]*rdxFnur; 
  outr[43] += incr2[43]*rdxFnur; 
  outr[44] += incr2[44]*rdxFnur; 
  outr[45] += incr2[45]*rdxFnur; 
  outr[46] += incr2[46]*rdxFnur; 
  outr[47] += incr2[47]*rdxFnur; 

  } else {


  incr2[4] = (-19.06233990711463*fl[14])-4.921875*fl[4]; 
  incr2[8] = (-19.06233990711463*fl[28])-4.921875*fl[8]; 
  incr2[9] = (-19.06233990711463*fl[29])-4.921875*fl[9]; 
  incr2[10] = (-19.06233990711463*fl[30])-4.921875*fl[10]; 
  incr2[14] = (-73.828125*fl[14])-19.06233990711463*fl[4]; 
  incr2[16] = (-19.06233990711463*fl[41])-4.921875*fl[16]; 
  incr2[17] = (-19.06233990711463*fl[42])-4.921875*fl[17]; 
  incr2[18] = (-19.06233990711463*fl[43])-4.921875*fl[18]; 
  incr2[25] = -4.921875*fl[25]; 
  incr2[26] = -4.921875*fl[26]; 
  incr2[27] = -4.921875*fl[27]; 
  incr2[28] = (-73.828125*fl[28])-19.06233990711463*fl[8]; 
  incr2[29] = (-73.828125*fl[29])-19.06233990711463*fl[9]; 
  incr2[30] = (-73.828125*fl[30])-19.06233990711463*fl[10]; 
  incr2[31] = (-19.06233990711463*fl[47])-4.921875*fl[31]; 
  incr2[35] = -4.921875*fl[35]; 
  incr2[36] = -4.921875*fl[36]; 
  incr2[37] = -4.921875*fl[37]; 
  incr2[38] = -4.921875*fl[38]; 
  incr2[39] = -4.921875*fl[39]; 
  incr2[40] = -4.921875*fl[40]; 
  incr2[41] = (-73.828125*fl[41])-19.06233990711463*fl[16]; 
  incr2[42] = (-73.828125*fl[42])-19.06233990711463*fl[17]; 
  incr2[43] = (-73.828125*fl[43])-19.06233990711463*fl[18]; 
  incr2[44] = -4.921875*fl[44]; 
  incr2[45] = -4.921875*fl[45]; 
  incr2[46] = -4.921875*fl[46]; 
  incr2[47] = (-73.828125*fl[47])-19.06233990711463*fl[31]; 





  outl[4] += incr2[4]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul; 
  outl[27] += incr2[27]*rdxFnul; 
  outl[28] += incr2[28]*rdxFnul; 
  outl[29] += incr2[29]*rdxFnul; 
  outl[30] += incr2[30]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul; 
  outl[35] += incr2[35]*rdxFnul; 
  outl[36] += incr2[36]*rdxFnul; 
  outl[37] += incr2[37]*rdxFnul; 
  outl[38] += incr2[38]*rdxFnul; 
  outl[39] += incr2[39]*rdxFnul; 
  outl[40] += incr2[40]*rdxFnul; 
  outl[41] += incr2[41]*rdxFnul; 
  outl[42] += incr2[42]*rdxFnul; 
  outl[43] += incr2[43]*rdxFnul; 
  outl[44] += incr2[44]*rdxFnul; 
  outl[45] += incr2[45]*rdxFnul; 
  outl[46] += incr2[46]*rdxFnul; 
  outl[47] += incr2[47]*rdxFnul; 

  }

} 
