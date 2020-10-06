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

  incr2[1] = 5.809475019311125*fr[11]; 
  incr2[5] = 5.809475019311126*fr[19]; 
  incr2[6] = 5.809475019311126*fr[21]; 
  incr2[8] = 5.809475019311126*fr[25]; 
  incr2[11] = -22.5*fr[11]; 
  incr2[15] = 5.809475019311125*fr[32]; 
  incr2[16] = 5.809475019311125*fr[35]; 
  incr2[17] = 5.809475019311125*fr[37]; 
  incr2[19] = -22.5*fr[19]; 
  incr2[21] = -22.5*fr[21]; 
  incr2[25] = -22.5*fr[25]; 
  incr2[31] = 5.809475019311126*fr[44]; 
  incr2[32] = -22.5*fr[32]; 
  incr2[35] = -22.5*fr[35]; 
  incr2[37] = -22.5*fr[37]; 
  incr2[44] = -22.5*fr[44]; 

  incr3[11] = 22.5*fr[11]-5.809475019311125*fr[1]; 
  incr3[19] = 22.5*fr[19]-5.809475019311126*fr[5]; 
  incr3[21] = 22.5*fr[21]-5.809475019311126*fr[6]; 
  incr3[25] = 22.5*fr[25]-5.809475019311126*fr[8]; 
  incr3[32] = 22.5*fr[32]-5.809475019311125*fr[15]; 
  incr3[35] = 22.5*fr[35]-5.809475019311125*fr[16]; 
  incr3[37] = 22.5*fr[37]-5.809475019311125*fr[17]; 
  incr3[44] = 22.5*fr[44]-5.809475019311126*fr[31]; 


  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[8] += -1.0*incr2[8]*rdxFnur; 
  outr[11] += (-1.0*incr3[11]*rdxFnur)-1.0*incr2[11]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[19] += (-1.0*incr3[19]*rdxFnur)-1.0*incr2[19]*rdxFnur; 
  outr[21] += (-1.0*incr3[21]*rdxFnur)-1.0*incr2[21]*rdxFnur; 
  outr[25] += (-1.0*incr3[25]*rdxFnur)-1.0*incr2[25]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 
  outr[32] += (-1.0*incr3[32]*rdxFnur)-1.0*incr2[32]*rdxFnur; 
  outr[35] += (-1.0*incr3[35]*rdxFnur)-1.0*incr2[35]*rdxFnur; 
  outr[37] += (-1.0*incr3[37]*rdxFnur)-1.0*incr2[37]*rdxFnur; 
  outr[44] += (-1.0*incr3[44]*rdxFnur)-1.0*incr2[44]*rdxFnur; 

  } else {

  incr2[1] = 5.809475019311125*fl[11]; 
  incr2[5] = 5.809475019311126*fl[19]; 
  incr2[6] = 5.809475019311126*fl[21]; 
  incr2[8] = 5.809475019311126*fl[25]; 
  incr2[11] = -22.5*fl[11]; 
  incr2[15] = 5.809475019311125*fl[32]; 
  incr2[16] = 5.809475019311125*fl[35]; 
  incr2[17] = 5.809475019311125*fl[37]; 
  incr2[19] = -22.5*fl[19]; 
  incr2[21] = -22.5*fl[21]; 
  incr2[25] = -22.5*fl[25]; 
  incr2[31] = 5.809475019311126*fl[44]; 
  incr2[32] = -22.5*fl[32]; 
  incr2[35] = -22.5*fl[35]; 
  incr2[37] = -22.5*fl[37]; 
  incr2[44] = -22.5*fl[44]; 

  incr3[11] = (-22.5*fl[11])-5.809475019311125*fl[1]; 
  incr3[19] = (-22.5*fl[19])-5.809475019311126*fl[5]; 
  incr3[21] = (-22.5*fl[21])-5.809475019311126*fl[6]; 
  incr3[25] = (-22.5*fl[25])-5.809475019311126*fl[8]; 
  incr3[32] = (-22.5*fl[32])-5.809475019311125*fl[15]; 
  incr3[35] = (-22.5*fl[35])-5.809475019311125*fl[16]; 
  incr3[37] = (-22.5*fl[37])-5.809475019311125*fl[17]; 
  incr3[44] = (-22.5*fl[44])-5.809475019311126*fl[31]; 


  outl[1] += incr2[1]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul; 
  outl[11] += incr3[11]*rdxFnul-1.0*incr2[11]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul; 
  outl[19] += incr3[19]*rdxFnul-1.0*incr2[19]*rdxFnul; 
  outl[21] += incr3[21]*rdxFnul-1.0*incr2[21]*rdxFnul; 
  outl[25] += incr3[25]*rdxFnul-1.0*incr2[25]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul; 
  outl[32] += incr3[32]*rdxFnul-1.0*incr2[32]*rdxFnul; 
  outl[35] += incr3[35]*rdxFnul-1.0*incr2[35]*rdxFnul; 
  outl[37] += incr3[37]*rdxFnul-1.0*incr2[37]*rdxFnul; 
  outl[44] += incr3[44]*rdxFnul-1.0*incr2[44]*rdxFnul; 

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

  incr2[2] = 5.809475019311125*fr[12]; 
  incr2[5] = 5.809475019311126*fr[20]; 
  incr2[7] = 5.809475019311126*fr[22]; 
  incr2[9] = 5.809475019311126*fr[26]; 
  incr2[12] = -22.5*fr[12]; 
  incr2[15] = 5.809475019311125*fr[33]; 
  incr2[16] = 5.809475019311125*fr[36]; 
  incr2[18] = 5.809475019311125*fr[38]; 
  incr2[20] = -22.5*fr[20]; 
  incr2[22] = -22.5*fr[22]; 
  incr2[26] = -22.5*fr[26]; 
  incr2[31] = 5.809475019311126*fr[45]; 
  incr2[33] = -22.5*fr[33]; 
  incr2[36] = -22.5*fr[36]; 
  incr2[38] = -22.5*fr[38]; 
  incr2[45] = -22.5*fr[45]; 

  incr3[12] = 22.5*fr[12]-5.809475019311125*fr[2]; 
  incr3[20] = 22.5*fr[20]-5.809475019311126*fr[5]; 
  incr3[22] = 22.5*fr[22]-5.809475019311126*fr[7]; 
  incr3[26] = 22.5*fr[26]-5.809475019311126*fr[9]; 
  incr3[33] = 22.5*fr[33]-5.809475019311125*fr[15]; 
  incr3[36] = 22.5*fr[36]-5.809475019311125*fr[16]; 
  incr3[38] = 22.5*fr[38]-5.809475019311125*fr[18]; 
  incr3[45] = 22.5*fr[45]-5.809475019311126*fr[31]; 


  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[9] += -1.0*incr2[9]*rdxFnur; 
  outr[12] += (-1.0*incr3[12]*rdxFnur)-1.0*incr2[12]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[20] += (-1.0*incr3[20]*rdxFnur)-1.0*incr2[20]*rdxFnur; 
  outr[22] += (-1.0*incr3[22]*rdxFnur)-1.0*incr2[22]*rdxFnur; 
  outr[26] += (-1.0*incr3[26]*rdxFnur)-1.0*incr2[26]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 
  outr[33] += (-1.0*incr3[33]*rdxFnur)-1.0*incr2[33]*rdxFnur; 
  outr[36] += (-1.0*incr3[36]*rdxFnur)-1.0*incr2[36]*rdxFnur; 
  outr[38] += (-1.0*incr3[38]*rdxFnur)-1.0*incr2[38]*rdxFnur; 
  outr[45] += (-1.0*incr3[45]*rdxFnur)-1.0*incr2[45]*rdxFnur; 

  } else {

  incr2[2] = 5.809475019311125*fl[12]; 
  incr2[5] = 5.809475019311126*fl[20]; 
  incr2[7] = 5.809475019311126*fl[22]; 
  incr2[9] = 5.809475019311126*fl[26]; 
  incr2[12] = -22.5*fl[12]; 
  incr2[15] = 5.809475019311125*fl[33]; 
  incr2[16] = 5.809475019311125*fl[36]; 
  incr2[18] = 5.809475019311125*fl[38]; 
  incr2[20] = -22.5*fl[20]; 
  incr2[22] = -22.5*fl[22]; 
  incr2[26] = -22.5*fl[26]; 
  incr2[31] = 5.809475019311126*fl[45]; 
  incr2[33] = -22.5*fl[33]; 
  incr2[36] = -22.5*fl[36]; 
  incr2[38] = -22.5*fl[38]; 
  incr2[45] = -22.5*fl[45]; 

  incr3[12] = (-22.5*fl[12])-5.809475019311125*fl[2]; 
  incr3[20] = (-22.5*fl[20])-5.809475019311126*fl[5]; 
  incr3[22] = (-22.5*fl[22])-5.809475019311126*fl[7]; 
  incr3[26] = (-22.5*fl[26])-5.809475019311126*fl[9]; 
  incr3[33] = (-22.5*fl[33])-5.809475019311125*fl[15]; 
  incr3[36] = (-22.5*fl[36])-5.809475019311125*fl[16]; 
  incr3[38] = (-22.5*fl[38])-5.809475019311125*fl[18]; 
  incr3[45] = (-22.5*fl[45])-5.809475019311126*fl[31]; 


  outl[2] += incr2[2]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul; 
  outl[12] += incr3[12]*rdxFnul-1.0*incr2[12]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul; 
  outl[20] += incr3[20]*rdxFnul-1.0*incr2[20]*rdxFnul; 
  outl[22] += incr3[22]*rdxFnul-1.0*incr2[22]*rdxFnul; 
  outl[26] += incr3[26]*rdxFnul-1.0*incr2[26]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul; 
  outl[33] += incr3[33]*rdxFnul-1.0*incr2[33]*rdxFnul; 
  outl[36] += incr3[36]*rdxFnul-1.0*incr2[36]*rdxFnul; 
  outl[38] += incr3[38]*rdxFnul-1.0*incr2[38]*rdxFnul; 
  outl[45] += incr3[45]*rdxFnul-1.0*incr2[45]*rdxFnul; 

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

  incr2[3] = 5.809475019311125*fr[13]; 
  incr2[6] = 5.809475019311126*fr[23]; 
  incr2[7] = 5.809475019311126*fr[24]; 
  incr2[10] = 5.809475019311126*fr[27]; 
  incr2[13] = -22.5*fr[13]; 
  incr2[15] = 5.809475019311125*fr[34]; 
  incr2[17] = 5.809475019311125*fr[39]; 
  incr2[18] = 5.809475019311125*fr[40]; 
  incr2[23] = -22.5*fr[23]; 
  incr2[24] = -22.5*fr[24]; 
  incr2[27] = -22.5*fr[27]; 
  incr2[31] = 5.809475019311126*fr[46]; 
  incr2[34] = -22.5*fr[34]; 
  incr2[39] = -22.5*fr[39]; 
  incr2[40] = -22.5*fr[40]; 
  incr2[46] = -22.5*fr[46]; 

  incr3[13] = 22.5*fr[13]-5.809475019311125*fr[3]; 
  incr3[23] = 22.5*fr[23]-5.809475019311126*fr[6]; 
  incr3[24] = 22.5*fr[24]-5.809475019311126*fr[7]; 
  incr3[27] = 22.5*fr[27]-5.809475019311126*fr[10]; 
  incr3[34] = 22.5*fr[34]-5.809475019311125*fr[15]; 
  incr3[39] = 22.5*fr[39]-5.809475019311125*fr[17]; 
  incr3[40] = 22.5*fr[40]-5.809475019311125*fr[18]; 
  incr3[46] = 22.5*fr[46]-5.809475019311126*fr[31]; 


  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[13] += (-1.0*incr3[13]*rdxFnur)-1.0*incr2[13]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[23] += (-1.0*incr3[23]*rdxFnur)-1.0*incr2[23]*rdxFnur; 
  outr[24] += (-1.0*incr3[24]*rdxFnur)-1.0*incr2[24]*rdxFnur; 
  outr[27] += (-1.0*incr3[27]*rdxFnur)-1.0*incr2[27]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 
  outr[34] += (-1.0*incr3[34]*rdxFnur)-1.0*incr2[34]*rdxFnur; 
  outr[39] += (-1.0*incr3[39]*rdxFnur)-1.0*incr2[39]*rdxFnur; 
  outr[40] += (-1.0*incr3[40]*rdxFnur)-1.0*incr2[40]*rdxFnur; 
  outr[46] += (-1.0*incr3[46]*rdxFnur)-1.0*incr2[46]*rdxFnur; 

  } else {

  incr2[3] = 5.809475019311125*fl[13]; 
  incr2[6] = 5.809475019311126*fl[23]; 
  incr2[7] = 5.809475019311126*fl[24]; 
  incr2[10] = 5.809475019311126*fl[27]; 
  incr2[13] = -22.5*fl[13]; 
  incr2[15] = 5.809475019311125*fl[34]; 
  incr2[17] = 5.809475019311125*fl[39]; 
  incr2[18] = 5.809475019311125*fl[40]; 
  incr2[23] = -22.5*fl[23]; 
  incr2[24] = -22.5*fl[24]; 
  incr2[27] = -22.5*fl[27]; 
  incr2[31] = 5.809475019311126*fl[46]; 
  incr2[34] = -22.5*fl[34]; 
  incr2[39] = -22.5*fl[39]; 
  incr2[40] = -22.5*fl[40]; 
  incr2[46] = -22.5*fl[46]; 

  incr3[13] = (-22.5*fl[13])-5.809475019311125*fl[3]; 
  incr3[23] = (-22.5*fl[23])-5.809475019311126*fl[6]; 
  incr3[24] = (-22.5*fl[24])-5.809475019311126*fl[7]; 
  incr3[27] = (-22.5*fl[27])-5.809475019311126*fl[10]; 
  incr3[34] = (-22.5*fl[34])-5.809475019311125*fl[15]; 
  incr3[39] = (-22.5*fl[39])-5.809475019311125*fl[17]; 
  incr3[40] = (-22.5*fl[40])-5.809475019311125*fl[18]; 
  incr3[46] = (-22.5*fl[46])-5.809475019311126*fl[31]; 


  outl[3] += incr2[3]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[13] += incr3[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul; 
  outl[23] += incr3[23]*rdxFnul-1.0*incr2[23]*rdxFnul; 
  outl[24] += incr3[24]*rdxFnul-1.0*incr2[24]*rdxFnul; 
  outl[27] += incr3[27]*rdxFnul-1.0*incr2[27]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul; 
  outl[34] += incr3[34]*rdxFnul-1.0*incr2[34]*rdxFnul; 
  outl[39] += incr3[39]*rdxFnul-1.0*incr2[39]*rdxFnul; 
  outl[40] += incr3[40]*rdxFnul-1.0*incr2[40]*rdxFnul; 
  outl[46] += incr3[46]*rdxFnul-1.0*incr2[46]*rdxFnul; 

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

  incr2[4] = 5.809475019311125*fr[14]; 
  incr2[8] = 5.809475019311126*fr[28]; 
  incr2[9] = 5.809475019311126*fr[29]; 
  incr2[10] = 5.809475019311126*fr[30]; 
  incr2[14] = -22.5*fr[14]; 
  incr2[16] = 5.809475019311125*fr[41]; 
  incr2[17] = 5.809475019311125*fr[42]; 
  incr2[18] = 5.809475019311125*fr[43]; 
  incr2[28] = -22.5*fr[28]; 
  incr2[29] = -22.5*fr[29]; 
  incr2[30] = -22.5*fr[30]; 
  incr2[31] = 5.809475019311126*fr[47]; 
  incr2[41] = -22.5*fr[41]; 
  incr2[42] = -22.5*fr[42]; 
  incr2[43] = -22.5*fr[43]; 
  incr2[47] = -22.5*fr[47]; 

  incr3[14] = 22.5*fr[14]-5.809475019311125*fr[4]; 
  incr3[28] = 22.5*fr[28]-5.809475019311126*fr[8]; 
  incr3[29] = 22.5*fr[29]-5.809475019311126*fr[9]; 
  incr3[30] = 22.5*fr[30]-5.809475019311126*fr[10]; 
  incr3[41] = 22.5*fr[41]-5.809475019311125*fr[16]; 
  incr3[42] = 22.5*fr[42]-5.809475019311125*fr[17]; 
  incr3[43] = 22.5*fr[43]-5.809475019311125*fr[18]; 
  incr3[47] = 22.5*fr[47]-5.809475019311126*fr[31]; 


  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[8] += -1.0*incr2[8]*rdxFnur; 
  outr[9] += -1.0*incr2[9]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[14] += (-1.0*incr3[14]*rdxFnur)-1.0*incr2[14]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[28] += (-1.0*incr3[28]*rdxFnur)-1.0*incr2[28]*rdxFnur; 
  outr[29] += (-1.0*incr3[29]*rdxFnur)-1.0*incr2[29]*rdxFnur; 
  outr[30] += (-1.0*incr3[30]*rdxFnur)-1.0*incr2[30]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 
  outr[41] += (-1.0*incr3[41]*rdxFnur)-1.0*incr2[41]*rdxFnur; 
  outr[42] += (-1.0*incr3[42]*rdxFnur)-1.0*incr2[42]*rdxFnur; 
  outr[43] += (-1.0*incr3[43]*rdxFnur)-1.0*incr2[43]*rdxFnur; 
  outr[47] += (-1.0*incr3[47]*rdxFnur)-1.0*incr2[47]*rdxFnur; 

  } else {

  incr2[4] = 5.809475019311125*fl[14]; 
  incr2[8] = 5.809475019311126*fl[28]; 
  incr2[9] = 5.809475019311126*fl[29]; 
  incr2[10] = 5.809475019311126*fl[30]; 
  incr2[14] = -22.5*fl[14]; 
  incr2[16] = 5.809475019311125*fl[41]; 
  incr2[17] = 5.809475019311125*fl[42]; 
  incr2[18] = 5.809475019311125*fl[43]; 
  incr2[28] = -22.5*fl[28]; 
  incr2[29] = -22.5*fl[29]; 
  incr2[30] = -22.5*fl[30]; 
  incr2[31] = 5.809475019311126*fl[47]; 
  incr2[41] = -22.5*fl[41]; 
  incr2[42] = -22.5*fl[42]; 
  incr2[43] = -22.5*fl[43]; 
  incr2[47] = -22.5*fl[47]; 

  incr3[14] = (-22.5*fl[14])-5.809475019311125*fl[4]; 
  incr3[28] = (-22.5*fl[28])-5.809475019311126*fl[8]; 
  incr3[29] = (-22.5*fl[29])-5.809475019311126*fl[9]; 
  incr3[30] = (-22.5*fl[30])-5.809475019311126*fl[10]; 
  incr3[41] = (-22.5*fl[41])-5.809475019311125*fl[16]; 
  incr3[42] = (-22.5*fl[42])-5.809475019311125*fl[17]; 
  incr3[43] = (-22.5*fl[43])-5.809475019311125*fl[18]; 
  incr3[47] = (-22.5*fl[47])-5.809475019311126*fl[31]; 


  outl[4] += incr2[4]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[14] += incr3[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul; 
  outl[28] += incr3[28]*rdxFnul-1.0*incr2[28]*rdxFnul; 
  outl[29] += incr3[29]*rdxFnul-1.0*incr2[29]*rdxFnul; 
  outl[30] += incr3[30]*rdxFnul-1.0*incr2[30]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul; 
  outl[41] += incr3[41]*rdxFnul-1.0*incr2[41]*rdxFnul; 
  outl[42] += incr3[42]*rdxFnul-1.0*incr2[42]*rdxFnul; 
  outl[43] += incr3[43]*rdxFnul-1.0*incr2[43]*rdxFnul; 
  outl[47] += incr3[47]*rdxFnul-1.0*incr2[47]*rdxFnul; 

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







  } else {







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







  } else {







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







  } else {







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







  } else {







  }

} 
