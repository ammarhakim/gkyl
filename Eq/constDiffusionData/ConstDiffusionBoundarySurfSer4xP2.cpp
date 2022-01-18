#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf4xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[48]; 
  double incr2[48]; 

  if (edge < 0) {

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
void ConstDiffusionBoundarySurf4xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[48]; 
  double incr2[48]; 

  if (edge < 0) {

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
void ConstDiffusionBoundarySurf4xSerP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[48]; 
  double incr2[48]; 

  if (edge < 0) {

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
void ConstDiffusionBoundarySurf4xSerP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxFnur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

  double incr1[48]; 
  double incr2[48]; 

  if (edge < 0) {

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
void ConstHyperDiffusion4BoundarySurf4xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 

  if (edge < 0) {


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
void ConstHyperDiffusion4BoundarySurf4xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 

  if (edge < 0) {


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
void ConstHyperDiffusion4BoundarySurf4xSerP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 

  if (edge < 0) {


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
void ConstHyperDiffusion4BoundarySurf4xSerP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[4]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 16.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[48]; 
  double incr2[48]; 
  double incr3[48]; 
  double incr4[48]; 

  if (edge < 0) {


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
void ConstHyperDiffusion6BoundarySurf4xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
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

  if (edge < 0) {


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
void ConstHyperDiffusion6BoundarySurf4xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
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

  if (edge < 0) {


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
void ConstHyperDiffusion6BoundarySurf4xSerP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
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

  if (edge < 0) {


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
void ConstHyperDiffusion6BoundarySurf4xSerP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
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

  if (edge < 0) {


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
void ConstDiffusionVarCoeffBoundarySurf4xSerP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[192]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[48]; 
  double incr2[48]; 

  if (edge < 0) {

  incr2[1] = 1.082531754730548*fr[11]*nul[11]-0.8385254915624212*fr[1]*nul[11]+0.4841229182759271*fr[0]*nul[11]+0.8385254915624212*nul[1]*fr[11]+0.4841229182759271*nul[0]*fr[11]-0.6495190528383289*fr[1]*nul[1]+0.375*fr[0]*nul[1]-0.375*nul[0]*fr[1]+0.2165063509461096*fr[0]*nul[0]; 
  incr2[5] = 1.082531754730548*nul[11]*fr[19]+0.838525491562421*nul[1]*fr[19]+0.4841229182759271*nul[0]*fr[19]-0.8385254915624212*fr[5]*nul[11]+0.4841229182759271*fr[2]*nul[11]-0.6495190528383289*nul[1]*fr[5]-0.375*nul[0]*fr[5]+0.375*nul[1]*fr[2]+0.2165063509461096*nul[0]*fr[2]; 
  incr2[6] = 1.082531754730548*nul[11]*fr[21]+0.838525491562421*nul[1]*fr[21]+0.4841229182759271*nul[0]*fr[21]-0.8385254915624212*fr[6]*nul[11]+0.4841229182759271*fr[3]*nul[11]-0.6495190528383289*nul[1]*fr[6]-0.375*nul[0]*fr[6]+0.375*nul[1]*fr[3]+0.2165063509461096*nul[0]*fr[3]; 
  incr2[8] = 1.082531754730548*nul[11]*fr[25]+0.838525491562421*nul[1]*fr[25]+0.4841229182759271*nul[0]*fr[25]-0.8385254915624212*fr[8]*nul[11]+0.4841229182759271*fr[4]*nul[11]-0.6495190528383289*nul[1]*fr[8]-0.375*nul[0]*fr[8]+0.375*nul[1]*fr[4]+0.2165063509461096*nul[0]*fr[4]; 
  incr2[11] = (-4.192627457812106*fr[11]*nul[11])+3.247595264191645*fr[1]*nul[11]-1.875*fr[0]*nul[11]-3.247595264191643*nul[1]*fr[11]-1.875*nul[0]*fr[11]+2.515576474687264*fr[1]*nul[1]-1.452368754827781*fr[0]*nul[1]+1.452368754827781*nul[0]*fr[1]-0.8385254915624212*fr[0]*nul[0]; 
  incr2[15] = 1.082531754730548*nul[11]*fr[32]+0.8385254915624212*nul[1]*fr[32]+0.4841229182759271*nul[0]*fr[32]-0.8385254915624212*nul[11]*fr[15]-0.6495190528383289*nul[1]*fr[15]-0.375*nul[0]*fr[15]+0.4841229182759271*fr[7]*nul[11]+0.375*nul[1]*fr[7]+0.2165063509461096*nul[0]*fr[7]; 
  incr2[16] = 1.082531754730548*nul[11]*fr[35]+0.8385254915624212*nul[1]*fr[35]+0.4841229182759271*nul[0]*fr[35]-0.8385254915624212*nul[11]*fr[16]-0.6495190528383289*nul[1]*fr[16]-0.375*nul[0]*fr[16]+0.4841229182759271*fr[9]*nul[11]+0.375*nul[1]*fr[9]+0.2165063509461096*nul[0]*fr[9]; 
  incr2[17] = 1.082531754730548*nul[11]*fr[37]+0.8385254915624212*nul[1]*fr[37]+0.4841229182759271*nul[0]*fr[37]-0.8385254915624212*nul[11]*fr[17]-0.6495190528383289*nul[1]*fr[17]-0.375*nul[0]*fr[17]+0.4841229182759271*fr[10]*nul[11]+0.375*nul[1]*fr[10]+0.2165063509461096*nul[0]*fr[10]; 
  incr2[19] = (-4.192627457812106*nul[11]*fr[19])-3.247595264191643*nul[1]*fr[19]-1.875*nul[0]*fr[19]+3.247595264191645*fr[5]*nul[11]-1.875*fr[2]*nul[11]+2.515576474687263*nul[1]*fr[5]+1.452368754827781*nul[0]*fr[5]-1.452368754827781*nul[1]*fr[2]-0.8385254915624211*nul[0]*fr[2]; 
  incr2[20] = (-0.8385254915624212*nul[11]*fr[20])-0.6495190528383289*nul[1]*fr[20]-0.375*nul[0]*fr[20]+0.4841229182759271*nul[11]*fr[12]+0.375*nul[1]*fr[12]+0.2165063509461097*nul[0]*fr[12]; 
  incr2[21] = (-4.192627457812106*nul[11]*fr[21])-3.247595264191643*nul[1]*fr[21]-1.875*nul[0]*fr[21]+3.247595264191645*fr[6]*nul[11]-1.875*fr[3]*nul[11]+2.515576474687263*nul[1]*fr[6]+1.452368754827781*nul[0]*fr[6]-1.452368754827781*nul[1]*fr[3]-0.8385254915624211*nul[0]*fr[3]; 
  incr2[23] = (-0.8385254915624212*nul[11]*fr[23])-0.6495190528383289*nul[1]*fr[23]-0.375*nul[0]*fr[23]+0.4841229182759271*nul[11]*fr[13]+0.375*nul[1]*fr[13]+0.2165063509461097*nul[0]*fr[13]; 
  incr2[25] = (-4.192627457812106*nul[11]*fr[25])-3.247595264191643*nul[1]*fr[25]-1.875*nul[0]*fr[25]+3.247595264191645*fr[8]*nul[11]-1.875*fr[4]*nul[11]+2.515576474687263*nul[1]*fr[8]+1.452368754827781*nul[0]*fr[8]-1.452368754827781*nul[1]*fr[4]-0.8385254915624211*nul[0]*fr[4]; 
  incr2[28] = (-0.8385254915624212*nul[11]*fr[28])-0.6495190528383289*nul[1]*fr[28]-0.375*nul[0]*fr[28]+0.4841229182759271*nul[11]*fr[14]+0.375*nul[1]*fr[14]+0.2165063509461097*nul[0]*fr[14]; 
  incr2[31] = 1.082531754730548*nul[11]*fr[44]+0.838525491562421*nul[1]*fr[44]+0.4841229182759271*nul[0]*fr[44]-0.8385254915624212*nul[11]*fr[31]-0.6495190528383289*nul[1]*fr[31]-0.375*nul[0]*fr[31]+0.4841229182759271*nul[11]*fr[18]+0.375*nul[1]*fr[18]+0.2165063509461096*nul[0]*fr[18]; 
  incr2[32] = (-4.192627457812106*nul[11]*fr[32])-3.247595264191643*nul[1]*fr[32]-1.875*nul[0]*fr[32]+3.247595264191645*nul[11]*fr[15]+2.515576474687264*nul[1]*fr[15]+1.452368754827781*nul[0]*fr[15]-1.875*fr[7]*nul[11]-1.452368754827781*nul[1]*fr[7]-0.8385254915624212*nul[0]*fr[7]; 
  incr2[33] = (-0.8385254915624212*nul[11]*fr[33])-0.6495190528383289*nul[1]*fr[33]-0.375*nul[0]*fr[33]+0.4841229182759271*nul[11]*fr[22]+0.375*nul[1]*fr[22]+0.2165063509461097*nul[0]*fr[22]; 
  incr2[34] = (-0.8385254915624212*nul[11]*fr[34])-0.6495190528383289*nul[1]*fr[34]-0.375*nul[0]*fr[34]+0.4841229182759271*nul[11]*fr[24]+0.375*nul[1]*fr[24]+0.2165063509461097*nul[0]*fr[24]; 
  incr2[35] = (-4.192627457812106*nul[11]*fr[35])-3.247595264191643*nul[1]*fr[35]-1.875*nul[0]*fr[35]+3.247595264191645*nul[11]*fr[16]+2.515576474687264*nul[1]*fr[16]+1.452368754827781*nul[0]*fr[16]-1.875*fr[9]*nul[11]-1.452368754827781*nul[1]*fr[9]-0.8385254915624212*nul[0]*fr[9]; 
  incr2[36] = (-0.8385254915624212*nul[11]*fr[36])-0.6495190528383289*nul[1]*fr[36]-0.375*nul[0]*fr[36]+0.4841229182759271*nul[11]*fr[26]+0.375*nul[1]*fr[26]+0.2165063509461097*nul[0]*fr[26]; 
  incr2[37] = (-4.192627457812106*nul[11]*fr[37])-3.247595264191643*nul[1]*fr[37]-1.875*nul[0]*fr[37]+3.247595264191645*nul[11]*fr[17]+2.515576474687264*nul[1]*fr[17]+1.452368754827781*nul[0]*fr[17]-1.875*fr[10]*nul[11]-1.452368754827781*nul[1]*fr[10]-0.8385254915624212*nul[0]*fr[10]; 
  incr2[39] = (-0.8385254915624212*nul[11]*fr[39])-0.6495190528383289*nul[1]*fr[39]-0.375*nul[0]*fr[39]+0.4841229182759271*nul[11]*fr[27]+0.375*nul[1]*fr[27]+0.2165063509461097*nul[0]*fr[27]; 
  incr2[41] = (-0.8385254915624212*nul[11]*fr[41])-0.6495190528383289*nul[1]*fr[41]-0.375*nul[0]*fr[41]+0.4841229182759271*nul[11]*fr[29]+0.375*nul[1]*fr[29]+0.2165063509461097*nul[0]*fr[29]; 
  incr2[42] = (-0.8385254915624212*nul[11]*fr[42])-0.6495190528383289*nul[1]*fr[42]-0.375*nul[0]*fr[42]+0.4841229182759271*nul[11]*fr[30]+0.375*nul[1]*fr[30]+0.2165063509461097*nul[0]*fr[30]; 
  incr2[44] = (-4.192627457812106*nul[11]*fr[44])-3.247595264191643*nul[1]*fr[44]-1.875*nul[0]*fr[44]+3.247595264191645*nul[11]*fr[31]+2.515576474687263*nul[1]*fr[31]+1.452368754827781*nul[0]*fr[31]-1.875*nul[11]*fr[18]-1.452368754827781*nul[1]*fr[18]-0.8385254915624211*nul[0]*fr[18]; 
  incr2[45] = (-0.8385254915624212*nul[11]*fr[45])-0.6495190528383289*nul[1]*fr[45]-0.375*nul[0]*fr[45]+0.4841229182759271*nul[11]*fr[38]+0.375*nul[1]*fr[38]+0.2165063509461097*nul[0]*fr[38]; 
  incr2[46] = (-0.8385254915624212*nul[11]*fr[46])-0.6495190528383289*nul[1]*fr[46]-0.375*nul[0]*fr[46]+0.4841229182759271*nul[11]*fr[40]+0.375*nul[1]*fr[40]+0.2165063509461097*nul[0]*fr[40]; 
  incr2[47] = (-0.8385254915624212*nul[11]*fr[47])-0.6495190528383289*nul[1]*fr[47]-0.375*nul[0]*fr[47]+0.4841229182759271*nul[11]*fr[43]+0.375*nul[1]*fr[43]+0.2165063509461097*nul[0]*fr[43]; 

  outr[1] += incr2[1]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[8] += incr2[8]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 
  outr[16] += incr2[16]*rdxFr; 
  outr[17] += incr2[17]*rdxFr; 
  outr[19] += incr2[19]*rdxFr; 
  outr[20] += incr2[20]*rdxFr; 
  outr[21] += incr2[21]*rdxFr; 
  outr[23] += incr2[23]*rdxFr; 
  outr[25] += incr2[25]*rdxFr; 
  outr[28] += incr2[28]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 
  outr[32] += incr2[32]*rdxFr; 
  outr[33] += incr2[33]*rdxFr; 
  outr[34] += incr2[34]*rdxFr; 
  outr[35] += incr2[35]*rdxFr; 
  outr[36] += incr2[36]*rdxFr; 
  outr[37] += incr2[37]*rdxFr; 
  outr[39] += incr2[39]*rdxFr; 
  outr[41] += incr2[41]*rdxFr; 
  outr[42] += incr2[42]*rdxFr; 
  outr[44] += incr2[44]*rdxFr; 
  outr[45] += incr2[45]*rdxFr; 
  outr[46] += incr2[46]*rdxFr; 
  outr[47] += incr2[47]*rdxFr; 

  } else {

  incr2[1] = 1.082531754730548*fl[11]*nul[11]+0.8385254915624212*fl[1]*nul[11]+0.4841229182759271*fl[0]*nul[11]+0.8385254915624212*nul[1]*fl[11]+0.4841229182759271*nul[0]*fl[11]+0.6495190528383289*fl[1]*nul[1]+0.375*fl[0]*nul[1]+0.375*nul[0]*fl[1]+0.2165063509461096*fl[0]*nul[0]; 
  incr2[5] = 1.082531754730548*nul[11]*fl[19]+0.838525491562421*nul[1]*fl[19]+0.4841229182759271*nul[0]*fl[19]+0.8385254915624212*fl[5]*nul[11]+0.4841229182759271*fl[2]*nul[11]+0.6495190528383289*nul[1]*fl[5]+0.375*nul[0]*fl[5]+0.375*nul[1]*fl[2]+0.2165063509461096*nul[0]*fl[2]; 
  incr2[6] = 1.082531754730548*nul[11]*fl[21]+0.838525491562421*nul[1]*fl[21]+0.4841229182759271*nul[0]*fl[21]+0.8385254915624212*fl[6]*nul[11]+0.4841229182759271*fl[3]*nul[11]+0.6495190528383289*nul[1]*fl[6]+0.375*nul[0]*fl[6]+0.375*nul[1]*fl[3]+0.2165063509461096*nul[0]*fl[3]; 
  incr2[8] = 1.082531754730548*nul[11]*fl[25]+0.838525491562421*nul[1]*fl[25]+0.4841229182759271*nul[0]*fl[25]+0.8385254915624212*fl[8]*nul[11]+0.4841229182759271*fl[4]*nul[11]+0.6495190528383289*nul[1]*fl[8]+0.375*nul[0]*fl[8]+0.375*nul[1]*fl[4]+0.2165063509461096*nul[0]*fl[4]; 
  incr2[11] = (-4.192627457812106*fl[11]*nul[11])-3.247595264191645*fl[1]*nul[11]-1.875*fl[0]*nul[11]-3.247595264191643*nul[1]*fl[11]-1.875*nul[0]*fl[11]-2.515576474687264*fl[1]*nul[1]-1.452368754827781*fl[0]*nul[1]-1.452368754827781*nul[0]*fl[1]-0.8385254915624212*fl[0]*nul[0]; 
  incr2[15] = 1.082531754730548*nul[11]*fl[32]+0.8385254915624212*nul[1]*fl[32]+0.4841229182759271*nul[0]*fl[32]+0.8385254915624212*nul[11]*fl[15]+0.6495190528383289*nul[1]*fl[15]+0.375*nul[0]*fl[15]+0.4841229182759271*fl[7]*nul[11]+0.375*nul[1]*fl[7]+0.2165063509461096*nul[0]*fl[7]; 
  incr2[16] = 1.082531754730548*nul[11]*fl[35]+0.8385254915624212*nul[1]*fl[35]+0.4841229182759271*nul[0]*fl[35]+0.8385254915624212*nul[11]*fl[16]+0.6495190528383289*nul[1]*fl[16]+0.375*nul[0]*fl[16]+0.4841229182759271*fl[9]*nul[11]+0.375*nul[1]*fl[9]+0.2165063509461096*nul[0]*fl[9]; 
  incr2[17] = 1.082531754730548*nul[11]*fl[37]+0.8385254915624212*nul[1]*fl[37]+0.4841229182759271*nul[0]*fl[37]+0.8385254915624212*nul[11]*fl[17]+0.6495190528383289*nul[1]*fl[17]+0.375*nul[0]*fl[17]+0.4841229182759271*fl[10]*nul[11]+0.375*nul[1]*fl[10]+0.2165063509461096*nul[0]*fl[10]; 
  incr2[19] = (-4.192627457812106*nul[11]*fl[19])-3.247595264191643*nul[1]*fl[19]-1.875*nul[0]*fl[19]-3.247595264191645*fl[5]*nul[11]-1.875*fl[2]*nul[11]-2.515576474687263*nul[1]*fl[5]-1.452368754827781*nul[0]*fl[5]-1.452368754827781*nul[1]*fl[2]-0.8385254915624211*nul[0]*fl[2]; 
  incr2[20] = 0.8385254915624212*nul[11]*fl[20]+0.6495190528383289*nul[1]*fl[20]+0.375*nul[0]*fl[20]+0.4841229182759271*nul[11]*fl[12]+0.375*nul[1]*fl[12]+0.2165063509461097*nul[0]*fl[12]; 
  incr2[21] = (-4.192627457812106*nul[11]*fl[21])-3.247595264191643*nul[1]*fl[21]-1.875*nul[0]*fl[21]-3.247595264191645*fl[6]*nul[11]-1.875*fl[3]*nul[11]-2.515576474687263*nul[1]*fl[6]-1.452368754827781*nul[0]*fl[6]-1.452368754827781*nul[1]*fl[3]-0.8385254915624211*nul[0]*fl[3]; 
  incr2[23] = 0.8385254915624212*nul[11]*fl[23]+0.6495190528383289*nul[1]*fl[23]+0.375*nul[0]*fl[23]+0.4841229182759271*nul[11]*fl[13]+0.375*nul[1]*fl[13]+0.2165063509461097*nul[0]*fl[13]; 
  incr2[25] = (-4.192627457812106*nul[11]*fl[25])-3.247595264191643*nul[1]*fl[25]-1.875*nul[0]*fl[25]-3.247595264191645*fl[8]*nul[11]-1.875*fl[4]*nul[11]-2.515576474687263*nul[1]*fl[8]-1.452368754827781*nul[0]*fl[8]-1.452368754827781*nul[1]*fl[4]-0.8385254915624211*nul[0]*fl[4]; 
  incr2[28] = 0.8385254915624212*nul[11]*fl[28]+0.6495190528383289*nul[1]*fl[28]+0.375*nul[0]*fl[28]+0.4841229182759271*nul[11]*fl[14]+0.375*nul[1]*fl[14]+0.2165063509461097*nul[0]*fl[14]; 
  incr2[31] = 1.082531754730548*nul[11]*fl[44]+0.838525491562421*nul[1]*fl[44]+0.4841229182759271*nul[0]*fl[44]+0.8385254915624212*nul[11]*fl[31]+0.6495190528383289*nul[1]*fl[31]+0.375*nul[0]*fl[31]+0.4841229182759271*nul[11]*fl[18]+0.375*nul[1]*fl[18]+0.2165063509461096*nul[0]*fl[18]; 
  incr2[32] = (-4.192627457812106*nul[11]*fl[32])-3.247595264191643*nul[1]*fl[32]-1.875*nul[0]*fl[32]-3.247595264191645*nul[11]*fl[15]-2.515576474687264*nul[1]*fl[15]-1.452368754827781*nul[0]*fl[15]-1.875*fl[7]*nul[11]-1.452368754827781*nul[1]*fl[7]-0.8385254915624212*nul[0]*fl[7]; 
  incr2[33] = 0.8385254915624212*nul[11]*fl[33]+0.6495190528383289*nul[1]*fl[33]+0.375*nul[0]*fl[33]+0.4841229182759271*nul[11]*fl[22]+0.375*nul[1]*fl[22]+0.2165063509461097*nul[0]*fl[22]; 
  incr2[34] = 0.8385254915624212*nul[11]*fl[34]+0.6495190528383289*nul[1]*fl[34]+0.375*nul[0]*fl[34]+0.4841229182759271*nul[11]*fl[24]+0.375*nul[1]*fl[24]+0.2165063509461097*nul[0]*fl[24]; 
  incr2[35] = (-4.192627457812106*nul[11]*fl[35])-3.247595264191643*nul[1]*fl[35]-1.875*nul[0]*fl[35]-3.247595264191645*nul[11]*fl[16]-2.515576474687264*nul[1]*fl[16]-1.452368754827781*nul[0]*fl[16]-1.875*fl[9]*nul[11]-1.452368754827781*nul[1]*fl[9]-0.8385254915624212*nul[0]*fl[9]; 
  incr2[36] = 0.8385254915624212*nul[11]*fl[36]+0.6495190528383289*nul[1]*fl[36]+0.375*nul[0]*fl[36]+0.4841229182759271*nul[11]*fl[26]+0.375*nul[1]*fl[26]+0.2165063509461097*nul[0]*fl[26]; 
  incr2[37] = (-4.192627457812106*nul[11]*fl[37])-3.247595264191643*nul[1]*fl[37]-1.875*nul[0]*fl[37]-3.247595264191645*nul[11]*fl[17]-2.515576474687264*nul[1]*fl[17]-1.452368754827781*nul[0]*fl[17]-1.875*fl[10]*nul[11]-1.452368754827781*nul[1]*fl[10]-0.8385254915624212*nul[0]*fl[10]; 
  incr2[39] = 0.8385254915624212*nul[11]*fl[39]+0.6495190528383289*nul[1]*fl[39]+0.375*nul[0]*fl[39]+0.4841229182759271*nul[11]*fl[27]+0.375*nul[1]*fl[27]+0.2165063509461097*nul[0]*fl[27]; 
  incr2[41] = 0.8385254915624212*nul[11]*fl[41]+0.6495190528383289*nul[1]*fl[41]+0.375*nul[0]*fl[41]+0.4841229182759271*nul[11]*fl[29]+0.375*nul[1]*fl[29]+0.2165063509461097*nul[0]*fl[29]; 
  incr2[42] = 0.8385254915624212*nul[11]*fl[42]+0.6495190528383289*nul[1]*fl[42]+0.375*nul[0]*fl[42]+0.4841229182759271*nul[11]*fl[30]+0.375*nul[1]*fl[30]+0.2165063509461097*nul[0]*fl[30]; 
  incr2[44] = (-4.192627457812106*nul[11]*fl[44])-3.247595264191643*nul[1]*fl[44]-1.875*nul[0]*fl[44]-3.247595264191645*nul[11]*fl[31]-2.515576474687263*nul[1]*fl[31]-1.452368754827781*nul[0]*fl[31]-1.875*nul[11]*fl[18]-1.452368754827781*nul[1]*fl[18]-0.8385254915624211*nul[0]*fl[18]; 
  incr2[45] = 0.8385254915624212*nul[11]*fl[45]+0.6495190528383289*nul[1]*fl[45]+0.375*nul[0]*fl[45]+0.4841229182759271*nul[11]*fl[38]+0.375*nul[1]*fl[38]+0.2165063509461097*nul[0]*fl[38]; 
  incr2[46] = 0.8385254915624212*nul[11]*fl[46]+0.6495190528383289*nul[1]*fl[46]+0.375*nul[0]*fl[46]+0.4841229182759271*nul[11]*fl[40]+0.375*nul[1]*fl[40]+0.2165063509461097*nul[0]*fl[40]; 
  incr2[47] = 0.8385254915624212*nul[11]*fl[47]+0.6495190528383289*nul[1]*fl[47]+0.375*nul[0]*fl[47]+0.4841229182759271*nul[11]*fl[43]+0.375*nul[1]*fl[43]+0.2165063509461097*nul[0]*fl[43]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 
  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[8] += -1.0*incr2[8]*rdxFl; 
  outl[11] += incr2[11]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 
  outl[16] += -1.0*incr2[16]*rdxFl; 
  outl[17] += -1.0*incr2[17]*rdxFl; 
  outl[19] += incr2[19]*rdxFl; 
  outl[20] += -1.0*incr2[20]*rdxFl; 
  outl[21] += incr2[21]*rdxFl; 
  outl[23] += -1.0*incr2[23]*rdxFl; 
  outl[25] += incr2[25]*rdxFl; 
  outl[28] += -1.0*incr2[28]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 
  outl[32] += incr2[32]*rdxFl; 
  outl[33] += -1.0*incr2[33]*rdxFl; 
  outl[34] += -1.0*incr2[34]*rdxFl; 
  outl[35] += incr2[35]*rdxFl; 
  outl[36] += -1.0*incr2[36]*rdxFl; 
  outl[37] += incr2[37]*rdxFl; 
  outl[39] += -1.0*incr2[39]*rdxFl; 
  outl[41] += -1.0*incr2[41]*rdxFl; 
  outl[42] += -1.0*incr2[42]*rdxFl; 
  outl[44] += incr2[44]*rdxFl; 
  outl[45] += -1.0*incr2[45]*rdxFl; 
  outl[46] += -1.0*incr2[46]*rdxFl; 
  outl[47] += -1.0*incr2[47]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf4xSerP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[192]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[1]*dxl[1]); 
  double rdxFr = 4.0/(dxr[1]*dxr[1]); 

  double incr1[48]; 
  double incr2[48]; 

  if (edge < 0) {

  incr2[2] = 1.082531754730548*fr[12]*nul[60]-0.8385254915624212*fr[2]*nul[60]+0.4841229182759271*fr[0]*nul[60]+0.8385254915624212*fr[12]*nul[50]-0.6495190528383289*fr[2]*nul[50]+0.375*fr[0]*nul[50]+0.4841229182759271*fr[12]*nul[48]-0.375*fr[2]*nul[48]+0.2165063509461096*fr[0]*nul[48]; 
  incr2[5] = 1.082531754730548*fr[20]*nul[60]-0.8385254915624212*fr[5]*nul[60]+0.4841229182759271*fr[1]*nul[60]+0.838525491562421*fr[20]*nul[50]-0.6495190528383289*fr[5]*nul[50]+0.375*fr[1]*nul[50]+0.4841229182759271*fr[20]*nul[48]-0.375*fr[5]*nul[48]+0.2165063509461096*fr[1]*nul[48]; 
  incr2[7] = 1.082531754730548*fr[22]*nul[60]-0.8385254915624212*fr[7]*nul[60]+0.4841229182759271*fr[3]*nul[60]+0.838525491562421*fr[22]*nul[50]-0.6495190528383289*fr[7]*nul[50]+0.375*fr[3]*nul[50]+0.4841229182759271*fr[22]*nul[48]-0.375*fr[7]*nul[48]+0.2165063509461096*fr[3]*nul[48]; 
  incr2[9] = 1.082531754730548*fr[26]*nul[60]-0.8385254915624212*fr[9]*nul[60]+0.4841229182759271*fr[4]*nul[60]+0.838525491562421*fr[26]*nul[50]-0.6495190528383289*fr[9]*nul[50]+0.375*fr[4]*nul[50]+0.4841229182759271*fr[26]*nul[48]-0.375*fr[9]*nul[48]+0.2165063509461096*fr[4]*nul[48]; 
  incr2[12] = (-4.192627457812106*fr[12]*nul[60])+3.247595264191645*fr[2]*nul[60]-1.875*fr[0]*nul[60]-3.247595264191643*fr[12]*nul[50]+2.515576474687264*fr[2]*nul[50]-1.452368754827781*fr[0]*nul[50]-1.875*fr[12]*nul[48]+1.452368754827781*fr[2]*nul[48]-0.8385254915624212*fr[0]*nul[48]; 
  incr2[15] = 1.082531754730548*fr[33]*nul[60]-0.8385254915624212*fr[15]*nul[60]+0.4841229182759271*fr[6]*nul[60]+0.8385254915624212*fr[33]*nul[50]-0.6495190528383289*fr[15]*nul[50]+0.375*fr[6]*nul[50]+0.4841229182759271*fr[33]*nul[48]-0.375*fr[15]*nul[48]+0.2165063509461096*fr[6]*nul[48]; 
  incr2[16] = 1.082531754730548*fr[36]*nul[60]-0.8385254915624212*fr[16]*nul[60]+0.4841229182759271*fr[8]*nul[60]+0.8385254915624212*fr[36]*nul[50]-0.6495190528383289*fr[16]*nul[50]+0.375*fr[8]*nul[50]+0.4841229182759271*fr[36]*nul[48]-0.375*fr[16]*nul[48]+0.2165063509461096*fr[8]*nul[48]; 
  incr2[18] = 1.082531754730548*fr[38]*nul[60]-0.8385254915624212*fr[18]*nul[60]+0.4841229182759271*fr[10]*nul[60]+0.8385254915624212*fr[38]*nul[50]-0.6495190528383289*fr[18]*nul[50]+0.375*fr[10]*nul[50]+0.4841229182759271*fr[38]*nul[48]-0.375*fr[18]*nul[48]+0.2165063509461096*fr[10]*nul[48]; 
  incr2[19] = (-0.8385254915624212*fr[19]*nul[60])+0.4841229182759271*fr[11]*nul[60]-0.6495190528383289*fr[19]*nul[50]+0.375*fr[11]*nul[50]-0.375*fr[19]*nul[48]+0.2165063509461097*fr[11]*nul[48]; 
  incr2[20] = (-4.192627457812106*fr[20]*nul[60])+3.247595264191645*fr[5]*nul[60]-1.875*fr[1]*nul[60]-3.247595264191643*fr[20]*nul[50]+2.515576474687263*fr[5]*nul[50]-1.452368754827781*fr[1]*nul[50]-1.875*fr[20]*nul[48]+1.452368754827781*fr[5]*nul[48]-0.8385254915624211*fr[1]*nul[48]; 
  incr2[22] = (-4.192627457812106*fr[22]*nul[60])+3.247595264191645*fr[7]*nul[60]-1.875*fr[3]*nul[60]-3.247595264191643*fr[22]*nul[50]+2.515576474687263*fr[7]*nul[50]-1.452368754827781*fr[3]*nul[50]-1.875*fr[22]*nul[48]+1.452368754827781*fr[7]*nul[48]-0.8385254915624211*fr[3]*nul[48]; 
  incr2[24] = (-0.8385254915624212*fr[24]*nul[60])+0.4841229182759271*fr[13]*nul[60]-0.6495190528383289*fr[24]*nul[50]+0.375*fr[13]*nul[50]-0.375*fr[24]*nul[48]+0.2165063509461097*fr[13]*nul[48]; 
  incr2[26] = (-4.192627457812106*fr[26]*nul[60])+3.247595264191645*fr[9]*nul[60]-1.875*fr[4]*nul[60]-3.247595264191643*fr[26]*nul[50]+2.515576474687263*fr[9]*nul[50]-1.452368754827781*fr[4]*nul[50]-1.875*fr[26]*nul[48]+1.452368754827781*fr[9]*nul[48]-0.8385254915624211*fr[4]*nul[48]; 
  incr2[29] = (-0.8385254915624212*fr[29]*nul[60])+0.4841229182759271*fr[14]*nul[60]-0.6495190528383289*fr[29]*nul[50]+0.375*fr[14]*nul[50]-0.375*fr[29]*nul[48]+0.2165063509461097*fr[14]*nul[48]; 
  incr2[31] = 1.082531754730548*fr[45]*nul[60]-0.8385254915624212*fr[31]*nul[60]+0.4841229182759271*fr[17]*nul[60]+0.838525491562421*fr[45]*nul[50]-0.6495190528383289*fr[31]*nul[50]+0.375*fr[17]*nul[50]+0.4841229182759271*fr[45]*nul[48]-0.375*fr[31]*nul[48]+0.2165063509461096*fr[17]*nul[48]; 
  incr2[32] = (-0.8385254915624212*fr[32]*nul[60])+0.4841229182759271*fr[21]*nul[60]-0.6495190528383289*fr[32]*nul[50]+0.375*fr[21]*nul[50]-0.375*fr[32]*nul[48]+0.2165063509461097*fr[21]*nul[48]; 
  incr2[33] = (-4.192627457812106*fr[33]*nul[60])+3.247595264191645*fr[15]*nul[60]-1.875*fr[6]*nul[60]-3.247595264191643*fr[33]*nul[50]+2.515576474687264*fr[15]*nul[50]-1.452368754827781*fr[6]*nul[50]-1.875*fr[33]*nul[48]+1.452368754827781*fr[15]*nul[48]-0.8385254915624212*fr[6]*nul[48]; 
  incr2[34] = (-0.8385254915624212*fr[34]*nul[60])+0.4841229182759271*fr[23]*nul[60]-0.6495190528383289*fr[34]*nul[50]+0.375*fr[23]*nul[50]-0.375*fr[34]*nul[48]+0.2165063509461097*fr[23]*nul[48]; 
  incr2[35] = (-0.8385254915624212*fr[35]*nul[60])+0.4841229182759271*fr[25]*nul[60]-0.6495190528383289*fr[35]*nul[50]+0.375*fr[25]*nul[50]-0.375*fr[35]*nul[48]+0.2165063509461097*fr[25]*nul[48]; 
  incr2[36] = (-4.192627457812106*fr[36]*nul[60])+3.247595264191645*fr[16]*nul[60]-1.875*fr[8]*nul[60]-3.247595264191643*fr[36]*nul[50]+2.515576474687264*fr[16]*nul[50]-1.452368754827781*fr[8]*nul[50]-1.875*fr[36]*nul[48]+1.452368754827781*fr[16]*nul[48]-0.8385254915624212*fr[8]*nul[48]; 
  incr2[38] = (-4.192627457812106*fr[38]*nul[60])+3.247595264191645*fr[18]*nul[60]-1.875*fr[10]*nul[60]-3.247595264191643*fr[38]*nul[50]+2.515576474687264*fr[18]*nul[50]-1.452368754827781*fr[10]*nul[50]-1.875*fr[38]*nul[48]+1.452368754827781*fr[18]*nul[48]-0.8385254915624212*fr[10]*nul[48]; 
  incr2[40] = (-0.8385254915624212*fr[40]*nul[60])+0.4841229182759271*fr[27]*nul[60]-0.6495190528383289*fr[40]*nul[50]+0.375*fr[27]*nul[50]-0.375*fr[40]*nul[48]+0.2165063509461097*fr[27]*nul[48]; 
  incr2[41] = (-0.8385254915624212*fr[41]*nul[60])+0.4841229182759271*fr[28]*nul[60]-0.6495190528383289*fr[41]*nul[50]+0.375*fr[28]*nul[50]-0.375*fr[41]*nul[48]+0.2165063509461097*fr[28]*nul[48]; 
  incr2[43] = (-0.8385254915624212*fr[43]*nul[60])+0.4841229182759271*fr[30]*nul[60]-0.6495190528383289*fr[43]*nul[50]+0.375*fr[30]*nul[50]-0.375*fr[43]*nul[48]+0.2165063509461097*fr[30]*nul[48]; 
  incr2[44] = (-0.8385254915624212*fr[44]*nul[60])+0.4841229182759271*fr[37]*nul[60]-0.6495190528383289*fr[44]*nul[50]+0.375*fr[37]*nul[50]-0.375*fr[44]*nul[48]+0.2165063509461097*fr[37]*nul[48]; 
  incr2[45] = (-4.192627457812106*fr[45]*nul[60])+3.247595264191645*fr[31]*nul[60]-1.875*fr[17]*nul[60]-3.247595264191643*fr[45]*nul[50]+2.515576474687263*fr[31]*nul[50]-1.452368754827781*fr[17]*nul[50]-1.875*fr[45]*nul[48]+1.452368754827781*fr[31]*nul[48]-0.8385254915624211*fr[17]*nul[48]; 
  incr2[46] = (-0.8385254915624212*fr[46]*nul[60])+0.4841229182759271*fr[39]*nul[60]-0.6495190528383289*fr[46]*nul[50]+0.375*fr[39]*nul[50]-0.375*fr[46]*nul[48]+0.2165063509461097*fr[39]*nul[48]; 
  incr2[47] = (-0.8385254915624212*fr[47]*nul[60])+0.4841229182759271*fr[42]*nul[60]-0.6495190528383289*fr[47]*nul[50]+0.375*fr[42]*nul[50]-0.375*fr[47]*nul[48]+0.2165063509461097*fr[42]*nul[48]; 

  outr[2] += incr2[2]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 
  outr[9] += incr2[9]*rdxFr; 
  outr[12] += incr2[12]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 
  outr[16] += incr2[16]*rdxFr; 
  outr[18] += incr2[18]*rdxFr; 
  outr[19] += incr2[19]*rdxFr; 
  outr[20] += incr2[20]*rdxFr; 
  outr[22] += incr2[22]*rdxFr; 
  outr[24] += incr2[24]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 
  outr[29] += incr2[29]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 
  outr[32] += incr2[32]*rdxFr; 
  outr[33] += incr2[33]*rdxFr; 
  outr[34] += incr2[34]*rdxFr; 
  outr[35] += incr2[35]*rdxFr; 
  outr[36] += incr2[36]*rdxFr; 
  outr[38] += incr2[38]*rdxFr; 
  outr[40] += incr2[40]*rdxFr; 
  outr[41] += incr2[41]*rdxFr; 
  outr[43] += incr2[43]*rdxFr; 
  outr[44] += incr2[44]*rdxFr; 
  outr[45] += incr2[45]*rdxFr; 
  outr[46] += incr2[46]*rdxFr; 
  outr[47] += incr2[47]*rdxFr; 

  } else {

  incr2[2] = 1.082531754730548*fl[12]*nul[60]+0.8385254915624212*fl[2]*nul[60]+0.4841229182759271*fl[0]*nul[60]+0.8385254915624212*fl[12]*nul[50]+0.6495190528383289*fl[2]*nul[50]+0.375*fl[0]*nul[50]+0.4841229182759271*fl[12]*nul[48]+0.375*fl[2]*nul[48]+0.2165063509461096*fl[0]*nul[48]; 
  incr2[5] = 1.082531754730548*fl[20]*nul[60]+0.8385254915624212*fl[5]*nul[60]+0.4841229182759271*fl[1]*nul[60]+0.838525491562421*fl[20]*nul[50]+0.6495190528383289*fl[5]*nul[50]+0.375*fl[1]*nul[50]+0.4841229182759271*fl[20]*nul[48]+0.375*fl[5]*nul[48]+0.2165063509461096*fl[1]*nul[48]; 
  incr2[7] = 1.082531754730548*fl[22]*nul[60]+0.8385254915624212*fl[7]*nul[60]+0.4841229182759271*fl[3]*nul[60]+0.838525491562421*fl[22]*nul[50]+0.6495190528383289*fl[7]*nul[50]+0.375*fl[3]*nul[50]+0.4841229182759271*fl[22]*nul[48]+0.375*fl[7]*nul[48]+0.2165063509461096*fl[3]*nul[48]; 
  incr2[9] = 1.082531754730548*fl[26]*nul[60]+0.8385254915624212*fl[9]*nul[60]+0.4841229182759271*fl[4]*nul[60]+0.838525491562421*fl[26]*nul[50]+0.6495190528383289*fl[9]*nul[50]+0.375*fl[4]*nul[50]+0.4841229182759271*fl[26]*nul[48]+0.375*fl[9]*nul[48]+0.2165063509461096*fl[4]*nul[48]; 
  incr2[12] = (-4.192627457812106*fl[12]*nul[60])-3.247595264191645*fl[2]*nul[60]-1.875*fl[0]*nul[60]-3.247595264191643*fl[12]*nul[50]-2.515576474687264*fl[2]*nul[50]-1.452368754827781*fl[0]*nul[50]-1.875*fl[12]*nul[48]-1.452368754827781*fl[2]*nul[48]-0.8385254915624212*fl[0]*nul[48]; 
  incr2[15] = 1.082531754730548*fl[33]*nul[60]+0.8385254915624212*fl[15]*nul[60]+0.4841229182759271*fl[6]*nul[60]+0.8385254915624212*fl[33]*nul[50]+0.6495190528383289*fl[15]*nul[50]+0.375*fl[6]*nul[50]+0.4841229182759271*fl[33]*nul[48]+0.375*fl[15]*nul[48]+0.2165063509461096*fl[6]*nul[48]; 
  incr2[16] = 1.082531754730548*fl[36]*nul[60]+0.8385254915624212*fl[16]*nul[60]+0.4841229182759271*fl[8]*nul[60]+0.8385254915624212*fl[36]*nul[50]+0.6495190528383289*fl[16]*nul[50]+0.375*fl[8]*nul[50]+0.4841229182759271*fl[36]*nul[48]+0.375*fl[16]*nul[48]+0.2165063509461096*fl[8]*nul[48]; 
  incr2[18] = 1.082531754730548*fl[38]*nul[60]+0.8385254915624212*fl[18]*nul[60]+0.4841229182759271*fl[10]*nul[60]+0.8385254915624212*fl[38]*nul[50]+0.6495190528383289*fl[18]*nul[50]+0.375*fl[10]*nul[50]+0.4841229182759271*fl[38]*nul[48]+0.375*fl[18]*nul[48]+0.2165063509461096*fl[10]*nul[48]; 
  incr2[19] = 0.8385254915624212*fl[19]*nul[60]+0.4841229182759271*fl[11]*nul[60]+0.6495190528383289*fl[19]*nul[50]+0.375*fl[11]*nul[50]+0.375*fl[19]*nul[48]+0.2165063509461097*fl[11]*nul[48]; 
  incr2[20] = (-4.192627457812106*fl[20]*nul[60])-3.247595264191645*fl[5]*nul[60]-1.875*fl[1]*nul[60]-3.247595264191643*fl[20]*nul[50]-2.515576474687263*fl[5]*nul[50]-1.452368754827781*fl[1]*nul[50]-1.875*fl[20]*nul[48]-1.452368754827781*fl[5]*nul[48]-0.8385254915624211*fl[1]*nul[48]; 
  incr2[22] = (-4.192627457812106*fl[22]*nul[60])-3.247595264191645*fl[7]*nul[60]-1.875*fl[3]*nul[60]-3.247595264191643*fl[22]*nul[50]-2.515576474687263*fl[7]*nul[50]-1.452368754827781*fl[3]*nul[50]-1.875*fl[22]*nul[48]-1.452368754827781*fl[7]*nul[48]-0.8385254915624211*fl[3]*nul[48]; 
  incr2[24] = 0.8385254915624212*fl[24]*nul[60]+0.4841229182759271*fl[13]*nul[60]+0.6495190528383289*fl[24]*nul[50]+0.375*fl[13]*nul[50]+0.375*fl[24]*nul[48]+0.2165063509461097*fl[13]*nul[48]; 
  incr2[26] = (-4.192627457812106*fl[26]*nul[60])-3.247595264191645*fl[9]*nul[60]-1.875*fl[4]*nul[60]-3.247595264191643*fl[26]*nul[50]-2.515576474687263*fl[9]*nul[50]-1.452368754827781*fl[4]*nul[50]-1.875*fl[26]*nul[48]-1.452368754827781*fl[9]*nul[48]-0.8385254915624211*fl[4]*nul[48]; 
  incr2[29] = 0.8385254915624212*fl[29]*nul[60]+0.4841229182759271*fl[14]*nul[60]+0.6495190528383289*fl[29]*nul[50]+0.375*fl[14]*nul[50]+0.375*fl[29]*nul[48]+0.2165063509461097*fl[14]*nul[48]; 
  incr2[31] = 1.082531754730548*fl[45]*nul[60]+0.8385254915624212*fl[31]*nul[60]+0.4841229182759271*fl[17]*nul[60]+0.838525491562421*fl[45]*nul[50]+0.6495190528383289*fl[31]*nul[50]+0.375*fl[17]*nul[50]+0.4841229182759271*fl[45]*nul[48]+0.375*fl[31]*nul[48]+0.2165063509461096*fl[17]*nul[48]; 
  incr2[32] = 0.8385254915624212*fl[32]*nul[60]+0.4841229182759271*fl[21]*nul[60]+0.6495190528383289*fl[32]*nul[50]+0.375*fl[21]*nul[50]+0.375*fl[32]*nul[48]+0.2165063509461097*fl[21]*nul[48]; 
  incr2[33] = (-4.192627457812106*fl[33]*nul[60])-3.247595264191645*fl[15]*nul[60]-1.875*fl[6]*nul[60]-3.247595264191643*fl[33]*nul[50]-2.515576474687264*fl[15]*nul[50]-1.452368754827781*fl[6]*nul[50]-1.875*fl[33]*nul[48]-1.452368754827781*fl[15]*nul[48]-0.8385254915624212*fl[6]*nul[48]; 
  incr2[34] = 0.8385254915624212*fl[34]*nul[60]+0.4841229182759271*fl[23]*nul[60]+0.6495190528383289*fl[34]*nul[50]+0.375*fl[23]*nul[50]+0.375*fl[34]*nul[48]+0.2165063509461097*fl[23]*nul[48]; 
  incr2[35] = 0.8385254915624212*fl[35]*nul[60]+0.4841229182759271*fl[25]*nul[60]+0.6495190528383289*fl[35]*nul[50]+0.375*fl[25]*nul[50]+0.375*fl[35]*nul[48]+0.2165063509461097*fl[25]*nul[48]; 
  incr2[36] = (-4.192627457812106*fl[36]*nul[60])-3.247595264191645*fl[16]*nul[60]-1.875*fl[8]*nul[60]-3.247595264191643*fl[36]*nul[50]-2.515576474687264*fl[16]*nul[50]-1.452368754827781*fl[8]*nul[50]-1.875*fl[36]*nul[48]-1.452368754827781*fl[16]*nul[48]-0.8385254915624212*fl[8]*nul[48]; 
  incr2[38] = (-4.192627457812106*fl[38]*nul[60])-3.247595264191645*fl[18]*nul[60]-1.875*fl[10]*nul[60]-3.247595264191643*fl[38]*nul[50]-2.515576474687264*fl[18]*nul[50]-1.452368754827781*fl[10]*nul[50]-1.875*fl[38]*nul[48]-1.452368754827781*fl[18]*nul[48]-0.8385254915624212*fl[10]*nul[48]; 
  incr2[40] = 0.8385254915624212*fl[40]*nul[60]+0.4841229182759271*fl[27]*nul[60]+0.6495190528383289*fl[40]*nul[50]+0.375*fl[27]*nul[50]+0.375*fl[40]*nul[48]+0.2165063509461097*fl[27]*nul[48]; 
  incr2[41] = 0.8385254915624212*fl[41]*nul[60]+0.4841229182759271*fl[28]*nul[60]+0.6495190528383289*fl[41]*nul[50]+0.375*fl[28]*nul[50]+0.375*fl[41]*nul[48]+0.2165063509461097*fl[28]*nul[48]; 
  incr2[43] = 0.8385254915624212*fl[43]*nul[60]+0.4841229182759271*fl[30]*nul[60]+0.6495190528383289*fl[43]*nul[50]+0.375*fl[30]*nul[50]+0.375*fl[43]*nul[48]+0.2165063509461097*fl[30]*nul[48]; 
  incr2[44] = 0.8385254915624212*fl[44]*nul[60]+0.4841229182759271*fl[37]*nul[60]+0.6495190528383289*fl[44]*nul[50]+0.375*fl[37]*nul[50]+0.375*fl[44]*nul[48]+0.2165063509461097*fl[37]*nul[48]; 
  incr2[45] = (-4.192627457812106*fl[45]*nul[60])-3.247595264191645*fl[31]*nul[60]-1.875*fl[17]*nul[60]-3.247595264191643*fl[45]*nul[50]-2.515576474687263*fl[31]*nul[50]-1.452368754827781*fl[17]*nul[50]-1.875*fl[45]*nul[48]-1.452368754827781*fl[31]*nul[48]-0.8385254915624211*fl[17]*nul[48]; 
  incr2[46] = 0.8385254915624212*fl[46]*nul[60]+0.4841229182759271*fl[39]*nul[60]+0.6495190528383289*fl[46]*nul[50]+0.375*fl[39]*nul[50]+0.375*fl[46]*nul[48]+0.2165063509461097*fl[39]*nul[48]; 
  incr2[47] = 0.8385254915624212*fl[47]*nul[60]+0.4841229182759271*fl[42]*nul[60]+0.6495190528383289*fl[47]*nul[50]+0.375*fl[42]*nul[50]+0.375*fl[47]*nul[48]+0.2165063509461097*fl[42]*nul[48]; 

  outl[2] += -1.0*incr2[2]*rdxFl; 
  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 
  outl[9] += -1.0*incr2[9]*rdxFl; 
  outl[12] += incr2[12]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 
  outl[16] += -1.0*incr2[16]*rdxFl; 
  outl[18] += -1.0*incr2[18]*rdxFl; 
  outl[19] += -1.0*incr2[19]*rdxFl; 
  outl[20] += incr2[20]*rdxFl; 
  outl[22] += incr2[22]*rdxFl; 
  outl[24] += -1.0*incr2[24]*rdxFl; 
  outl[26] += incr2[26]*rdxFl; 
  outl[29] += -1.0*incr2[29]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 
  outl[32] += -1.0*incr2[32]*rdxFl; 
  outl[33] += incr2[33]*rdxFl; 
  outl[34] += -1.0*incr2[34]*rdxFl; 
  outl[35] += -1.0*incr2[35]*rdxFl; 
  outl[36] += incr2[36]*rdxFl; 
  outl[38] += incr2[38]*rdxFl; 
  outl[40] += -1.0*incr2[40]*rdxFl; 
  outl[41] += -1.0*incr2[41]*rdxFl; 
  outl[43] += -1.0*incr2[43]*rdxFl; 
  outl[44] += -1.0*incr2[44]*rdxFl; 
  outl[45] += incr2[45]*rdxFl; 
  outl[46] += -1.0*incr2[46]*rdxFl; 
  outl[47] += -1.0*incr2[47]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf4xSerP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[192]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[2]*dxl[2]); 
  double rdxFr = 4.0/(dxr[2]*dxr[2]); 

  double incr1[48]; 
  double incr2[48]; 

  if (edge < 0) {

  incr2[3] = 1.082531754730548*fr[13]*nul[109]-0.8385254915624212*fr[3]*nul[109]+0.4841229182759271*fr[0]*nul[109]+0.8385254915624212*fr[13]*nul[99]-0.6495190528383289*fr[3]*nul[99]+0.375*fr[0]*nul[99]+0.4841229182759271*fr[13]*nul[96]-0.375*fr[3]*nul[96]+0.2165063509461096*fr[0]*nul[96]; 
  incr2[6] = 1.082531754730548*fr[23]*nul[109]-0.8385254915624212*fr[6]*nul[109]+0.4841229182759271*fr[1]*nul[109]+0.838525491562421*fr[23]*nul[99]-0.6495190528383289*fr[6]*nul[99]+0.375*fr[1]*nul[99]+0.4841229182759271*fr[23]*nul[96]-0.375*fr[6]*nul[96]+0.2165063509461096*fr[1]*nul[96]; 
  incr2[7] = 1.082531754730548*fr[24]*nul[109]-0.8385254915624212*fr[7]*nul[109]+0.4841229182759271*fr[2]*nul[109]+0.838525491562421*fr[24]*nul[99]-0.6495190528383289*fr[7]*nul[99]+0.375*fr[2]*nul[99]+0.4841229182759271*fr[24]*nul[96]-0.375*fr[7]*nul[96]+0.2165063509461096*fr[2]*nul[96]; 
  incr2[10] = 1.082531754730548*fr[27]*nul[109]-0.8385254915624212*fr[10]*nul[109]+0.4841229182759271*fr[4]*nul[109]+0.838525491562421*fr[27]*nul[99]-0.6495190528383289*fr[10]*nul[99]+0.375*fr[4]*nul[99]+0.4841229182759271*fr[27]*nul[96]-0.375*fr[10]*nul[96]+0.2165063509461096*fr[4]*nul[96]; 
  incr2[13] = (-4.192627457812106*fr[13]*nul[109])+3.247595264191645*fr[3]*nul[109]-1.875*fr[0]*nul[109]-3.247595264191643*fr[13]*nul[99]+2.515576474687264*fr[3]*nul[99]-1.452368754827781*fr[0]*nul[99]-1.875*fr[13]*nul[96]+1.452368754827781*fr[3]*nul[96]-0.8385254915624212*fr[0]*nul[96]; 
  incr2[15] = 1.082531754730548*fr[34]*nul[109]-0.8385254915624212*fr[15]*nul[109]+0.4841229182759271*fr[5]*nul[109]+0.8385254915624212*fr[34]*nul[99]-0.6495190528383289*fr[15]*nul[99]+0.375*fr[5]*nul[99]+0.4841229182759271*fr[34]*nul[96]-0.375*fr[15]*nul[96]+0.2165063509461096*fr[5]*nul[96]; 
  incr2[17] = 1.082531754730548*fr[39]*nul[109]-0.8385254915624212*fr[17]*nul[109]+0.4841229182759271*fr[8]*nul[109]+0.8385254915624212*fr[39]*nul[99]-0.6495190528383289*fr[17]*nul[99]+0.375*fr[8]*nul[99]+0.4841229182759271*fr[39]*nul[96]-0.375*fr[17]*nul[96]+0.2165063509461096*fr[8]*nul[96]; 
  incr2[18] = 1.082531754730548*fr[40]*nul[109]-0.8385254915624212*fr[18]*nul[109]+0.4841229182759271*fr[9]*nul[109]+0.8385254915624212*fr[40]*nul[99]-0.6495190528383289*fr[18]*nul[99]+0.375*fr[9]*nul[99]+0.4841229182759271*fr[40]*nul[96]-0.375*fr[18]*nul[96]+0.2165063509461096*fr[9]*nul[96]; 
  incr2[21] = (-0.8385254915624212*fr[21]*nul[109])+0.4841229182759271*fr[11]*nul[109]-0.6495190528383289*fr[21]*nul[99]+0.375*fr[11]*nul[99]-0.375*fr[21]*nul[96]+0.2165063509461097*fr[11]*nul[96]; 
  incr2[22] = (-0.8385254915624212*fr[22]*nul[109])+0.4841229182759271*fr[12]*nul[109]-0.6495190528383289*fr[22]*nul[99]+0.375*fr[12]*nul[99]-0.375*fr[22]*nul[96]+0.2165063509461097*fr[12]*nul[96]; 
  incr2[23] = (-4.192627457812106*fr[23]*nul[109])+3.247595264191645*fr[6]*nul[109]-1.875*fr[1]*nul[109]-3.247595264191643*fr[23]*nul[99]+2.515576474687263*fr[6]*nul[99]-1.452368754827781*fr[1]*nul[99]-1.875*fr[23]*nul[96]+1.452368754827781*fr[6]*nul[96]-0.8385254915624211*fr[1]*nul[96]; 
  incr2[24] = (-4.192627457812106*fr[24]*nul[109])+3.247595264191645*fr[7]*nul[109]-1.875*fr[2]*nul[109]-3.247595264191643*fr[24]*nul[99]+2.515576474687263*fr[7]*nul[99]-1.452368754827781*fr[2]*nul[99]-1.875*fr[24]*nul[96]+1.452368754827781*fr[7]*nul[96]-0.8385254915624211*fr[2]*nul[96]; 
  incr2[27] = (-4.192627457812106*fr[27]*nul[109])+3.247595264191645*fr[10]*nul[109]-1.875*fr[4]*nul[109]-3.247595264191643*fr[27]*nul[99]+2.515576474687263*fr[10]*nul[99]-1.452368754827781*fr[4]*nul[99]-1.875*fr[27]*nul[96]+1.452368754827781*fr[10]*nul[96]-0.8385254915624211*fr[4]*nul[96]; 
  incr2[30] = (-0.8385254915624212*fr[30]*nul[109])+0.4841229182759271*fr[14]*nul[109]-0.6495190528383289*fr[30]*nul[99]+0.375*fr[14]*nul[99]-0.375*fr[30]*nul[96]+0.2165063509461097*fr[14]*nul[96]; 
  incr2[31] = 1.082531754730548*fr[46]*nul[109]-0.8385254915624212*fr[31]*nul[109]+0.4841229182759271*fr[16]*nul[109]+0.838525491562421*fr[46]*nul[99]-0.6495190528383289*fr[31]*nul[99]+0.375*fr[16]*nul[99]+0.4841229182759271*fr[46]*nul[96]-0.375*fr[31]*nul[96]+0.2165063509461096*fr[16]*nul[96]; 
  incr2[32] = (-0.8385254915624212*fr[32]*nul[109])+0.4841229182759271*fr[19]*nul[109]-0.6495190528383289*fr[32]*nul[99]+0.375*fr[19]*nul[99]-0.375*fr[32]*nul[96]+0.2165063509461097*fr[19]*nul[96]; 
  incr2[33] = (-0.8385254915624212*fr[33]*nul[109])+0.4841229182759271*fr[20]*nul[109]-0.6495190528383289*fr[33]*nul[99]+0.375*fr[20]*nul[99]-0.375*fr[33]*nul[96]+0.2165063509461097*fr[20]*nul[96]; 
  incr2[34] = (-4.192627457812106*fr[34]*nul[109])+3.247595264191645*fr[15]*nul[109]-1.875*fr[5]*nul[109]-3.247595264191643*fr[34]*nul[99]+2.515576474687264*fr[15]*nul[99]-1.452368754827781*fr[5]*nul[99]-1.875*fr[34]*nul[96]+1.452368754827781*fr[15]*nul[96]-0.8385254915624212*fr[5]*nul[96]; 
  incr2[37] = (-0.8385254915624212*fr[37]*nul[109])+0.4841229182759271*fr[25]*nul[109]-0.6495190528383289*fr[37]*nul[99]+0.375*fr[25]*nul[99]-0.375*fr[37]*nul[96]+0.2165063509461097*fr[25]*nul[96]; 
  incr2[38] = (-0.8385254915624212*fr[38]*nul[109])+0.4841229182759271*fr[26]*nul[109]-0.6495190528383289*fr[38]*nul[99]+0.375*fr[26]*nul[99]-0.375*fr[38]*nul[96]+0.2165063509461097*fr[26]*nul[96]; 
  incr2[39] = (-4.192627457812106*fr[39]*nul[109])+3.247595264191645*fr[17]*nul[109]-1.875*fr[8]*nul[109]-3.247595264191643*fr[39]*nul[99]+2.515576474687264*fr[17]*nul[99]-1.452368754827781*fr[8]*nul[99]-1.875*fr[39]*nul[96]+1.452368754827781*fr[17]*nul[96]-0.8385254915624212*fr[8]*nul[96]; 
  incr2[40] = (-4.192627457812106*fr[40]*nul[109])+3.247595264191645*fr[18]*nul[109]-1.875*fr[9]*nul[109]-3.247595264191643*fr[40]*nul[99]+2.515576474687264*fr[18]*nul[99]-1.452368754827781*fr[9]*nul[99]-1.875*fr[40]*nul[96]+1.452368754827781*fr[18]*nul[96]-0.8385254915624212*fr[9]*nul[96]; 
  incr2[42] = (-0.8385254915624212*fr[42]*nul[109])+0.4841229182759271*fr[28]*nul[109]-0.6495190528383289*fr[42]*nul[99]+0.375*fr[28]*nul[99]-0.375*fr[42]*nul[96]+0.2165063509461097*fr[28]*nul[96]; 
  incr2[43] = (-0.8385254915624212*fr[43]*nul[109])+0.4841229182759271*fr[29]*nul[109]-0.6495190528383289*fr[43]*nul[99]+0.375*fr[29]*nul[99]-0.375*fr[43]*nul[96]+0.2165063509461097*fr[29]*nul[96]; 
  incr2[44] = (-0.8385254915624212*fr[44]*nul[109])+0.4841229182759271*fr[35]*nul[109]-0.6495190528383289*fr[44]*nul[99]+0.375*fr[35]*nul[99]-0.375*fr[44]*nul[96]+0.2165063509461097*fr[35]*nul[96]; 
  incr2[45] = (-0.8385254915624212*fr[45]*nul[109])+0.4841229182759271*fr[36]*nul[109]-0.6495190528383289*fr[45]*nul[99]+0.375*fr[36]*nul[99]-0.375*fr[45]*nul[96]+0.2165063509461097*fr[36]*nul[96]; 
  incr2[46] = (-4.192627457812106*fr[46]*nul[109])+3.247595264191645*fr[31]*nul[109]-1.875*fr[16]*nul[109]-3.247595264191643*fr[46]*nul[99]+2.515576474687263*fr[31]*nul[99]-1.452368754827781*fr[16]*nul[99]-1.875*fr[46]*nul[96]+1.452368754827781*fr[31]*nul[96]-0.8385254915624211*fr[16]*nul[96]; 
  incr2[47] = (-0.8385254915624212*fr[47]*nul[109])+0.4841229182759271*fr[41]*nul[109]-0.6495190528383289*fr[47]*nul[99]+0.375*fr[41]*nul[99]-0.375*fr[47]*nul[96]+0.2165063509461097*fr[41]*nul[96]; 

  outr[3] += incr2[3]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[13] += incr2[13]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 
  outr[17] += incr2[17]*rdxFr; 
  outr[18] += incr2[18]*rdxFr; 
  outr[21] += incr2[21]*rdxFr; 
  outr[22] += incr2[22]*rdxFr; 
  outr[23] += incr2[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr; 
  outr[27] += incr2[27]*rdxFr; 
  outr[30] += incr2[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 
  outr[32] += incr2[32]*rdxFr; 
  outr[33] += incr2[33]*rdxFr; 
  outr[34] += incr2[34]*rdxFr; 
  outr[37] += incr2[37]*rdxFr; 
  outr[38] += incr2[38]*rdxFr; 
  outr[39] += incr2[39]*rdxFr; 
  outr[40] += incr2[40]*rdxFr; 
  outr[42] += incr2[42]*rdxFr; 
  outr[43] += incr2[43]*rdxFr; 
  outr[44] += incr2[44]*rdxFr; 
  outr[45] += incr2[45]*rdxFr; 
  outr[46] += incr2[46]*rdxFr; 
  outr[47] += incr2[47]*rdxFr; 

  } else {

  incr2[3] = 1.082531754730548*fl[13]*nul[109]+0.8385254915624212*fl[3]*nul[109]+0.4841229182759271*fl[0]*nul[109]+0.8385254915624212*fl[13]*nul[99]+0.6495190528383289*fl[3]*nul[99]+0.375*fl[0]*nul[99]+0.4841229182759271*fl[13]*nul[96]+0.375*fl[3]*nul[96]+0.2165063509461096*fl[0]*nul[96]; 
  incr2[6] = 1.082531754730548*fl[23]*nul[109]+0.8385254915624212*fl[6]*nul[109]+0.4841229182759271*fl[1]*nul[109]+0.838525491562421*fl[23]*nul[99]+0.6495190528383289*fl[6]*nul[99]+0.375*fl[1]*nul[99]+0.4841229182759271*fl[23]*nul[96]+0.375*fl[6]*nul[96]+0.2165063509461096*fl[1]*nul[96]; 
  incr2[7] = 1.082531754730548*fl[24]*nul[109]+0.8385254915624212*fl[7]*nul[109]+0.4841229182759271*fl[2]*nul[109]+0.838525491562421*fl[24]*nul[99]+0.6495190528383289*fl[7]*nul[99]+0.375*fl[2]*nul[99]+0.4841229182759271*fl[24]*nul[96]+0.375*fl[7]*nul[96]+0.2165063509461096*fl[2]*nul[96]; 
  incr2[10] = 1.082531754730548*fl[27]*nul[109]+0.8385254915624212*fl[10]*nul[109]+0.4841229182759271*fl[4]*nul[109]+0.838525491562421*fl[27]*nul[99]+0.6495190528383289*fl[10]*nul[99]+0.375*fl[4]*nul[99]+0.4841229182759271*fl[27]*nul[96]+0.375*fl[10]*nul[96]+0.2165063509461096*fl[4]*nul[96]; 
  incr2[13] = (-4.192627457812106*fl[13]*nul[109])-3.247595264191645*fl[3]*nul[109]-1.875*fl[0]*nul[109]-3.247595264191643*fl[13]*nul[99]-2.515576474687264*fl[3]*nul[99]-1.452368754827781*fl[0]*nul[99]-1.875*fl[13]*nul[96]-1.452368754827781*fl[3]*nul[96]-0.8385254915624212*fl[0]*nul[96]; 
  incr2[15] = 1.082531754730548*fl[34]*nul[109]+0.8385254915624212*fl[15]*nul[109]+0.4841229182759271*fl[5]*nul[109]+0.8385254915624212*fl[34]*nul[99]+0.6495190528383289*fl[15]*nul[99]+0.375*fl[5]*nul[99]+0.4841229182759271*fl[34]*nul[96]+0.375*fl[15]*nul[96]+0.2165063509461096*fl[5]*nul[96]; 
  incr2[17] = 1.082531754730548*fl[39]*nul[109]+0.8385254915624212*fl[17]*nul[109]+0.4841229182759271*fl[8]*nul[109]+0.8385254915624212*fl[39]*nul[99]+0.6495190528383289*fl[17]*nul[99]+0.375*fl[8]*nul[99]+0.4841229182759271*fl[39]*nul[96]+0.375*fl[17]*nul[96]+0.2165063509461096*fl[8]*nul[96]; 
  incr2[18] = 1.082531754730548*fl[40]*nul[109]+0.8385254915624212*fl[18]*nul[109]+0.4841229182759271*fl[9]*nul[109]+0.8385254915624212*fl[40]*nul[99]+0.6495190528383289*fl[18]*nul[99]+0.375*fl[9]*nul[99]+0.4841229182759271*fl[40]*nul[96]+0.375*fl[18]*nul[96]+0.2165063509461096*fl[9]*nul[96]; 
  incr2[21] = 0.8385254915624212*fl[21]*nul[109]+0.4841229182759271*fl[11]*nul[109]+0.6495190528383289*fl[21]*nul[99]+0.375*fl[11]*nul[99]+0.375*fl[21]*nul[96]+0.2165063509461097*fl[11]*nul[96]; 
  incr2[22] = 0.8385254915624212*fl[22]*nul[109]+0.4841229182759271*fl[12]*nul[109]+0.6495190528383289*fl[22]*nul[99]+0.375*fl[12]*nul[99]+0.375*fl[22]*nul[96]+0.2165063509461097*fl[12]*nul[96]; 
  incr2[23] = (-4.192627457812106*fl[23]*nul[109])-3.247595264191645*fl[6]*nul[109]-1.875*fl[1]*nul[109]-3.247595264191643*fl[23]*nul[99]-2.515576474687263*fl[6]*nul[99]-1.452368754827781*fl[1]*nul[99]-1.875*fl[23]*nul[96]-1.452368754827781*fl[6]*nul[96]-0.8385254915624211*fl[1]*nul[96]; 
  incr2[24] = (-4.192627457812106*fl[24]*nul[109])-3.247595264191645*fl[7]*nul[109]-1.875*fl[2]*nul[109]-3.247595264191643*fl[24]*nul[99]-2.515576474687263*fl[7]*nul[99]-1.452368754827781*fl[2]*nul[99]-1.875*fl[24]*nul[96]-1.452368754827781*fl[7]*nul[96]-0.8385254915624211*fl[2]*nul[96]; 
  incr2[27] = (-4.192627457812106*fl[27]*nul[109])-3.247595264191645*fl[10]*nul[109]-1.875*fl[4]*nul[109]-3.247595264191643*fl[27]*nul[99]-2.515576474687263*fl[10]*nul[99]-1.452368754827781*fl[4]*nul[99]-1.875*fl[27]*nul[96]-1.452368754827781*fl[10]*nul[96]-0.8385254915624211*fl[4]*nul[96]; 
  incr2[30] = 0.8385254915624212*fl[30]*nul[109]+0.4841229182759271*fl[14]*nul[109]+0.6495190528383289*fl[30]*nul[99]+0.375*fl[14]*nul[99]+0.375*fl[30]*nul[96]+0.2165063509461097*fl[14]*nul[96]; 
  incr2[31] = 1.082531754730548*fl[46]*nul[109]+0.8385254915624212*fl[31]*nul[109]+0.4841229182759271*fl[16]*nul[109]+0.838525491562421*fl[46]*nul[99]+0.6495190528383289*fl[31]*nul[99]+0.375*fl[16]*nul[99]+0.4841229182759271*fl[46]*nul[96]+0.375*fl[31]*nul[96]+0.2165063509461096*fl[16]*nul[96]; 
  incr2[32] = 0.8385254915624212*fl[32]*nul[109]+0.4841229182759271*fl[19]*nul[109]+0.6495190528383289*fl[32]*nul[99]+0.375*fl[19]*nul[99]+0.375*fl[32]*nul[96]+0.2165063509461097*fl[19]*nul[96]; 
  incr2[33] = 0.8385254915624212*fl[33]*nul[109]+0.4841229182759271*fl[20]*nul[109]+0.6495190528383289*fl[33]*nul[99]+0.375*fl[20]*nul[99]+0.375*fl[33]*nul[96]+0.2165063509461097*fl[20]*nul[96]; 
  incr2[34] = (-4.192627457812106*fl[34]*nul[109])-3.247595264191645*fl[15]*nul[109]-1.875*fl[5]*nul[109]-3.247595264191643*fl[34]*nul[99]-2.515576474687264*fl[15]*nul[99]-1.452368754827781*fl[5]*nul[99]-1.875*fl[34]*nul[96]-1.452368754827781*fl[15]*nul[96]-0.8385254915624212*fl[5]*nul[96]; 
  incr2[37] = 0.8385254915624212*fl[37]*nul[109]+0.4841229182759271*fl[25]*nul[109]+0.6495190528383289*fl[37]*nul[99]+0.375*fl[25]*nul[99]+0.375*fl[37]*nul[96]+0.2165063509461097*fl[25]*nul[96]; 
  incr2[38] = 0.8385254915624212*fl[38]*nul[109]+0.4841229182759271*fl[26]*nul[109]+0.6495190528383289*fl[38]*nul[99]+0.375*fl[26]*nul[99]+0.375*fl[38]*nul[96]+0.2165063509461097*fl[26]*nul[96]; 
  incr2[39] = (-4.192627457812106*fl[39]*nul[109])-3.247595264191645*fl[17]*nul[109]-1.875*fl[8]*nul[109]-3.247595264191643*fl[39]*nul[99]-2.515576474687264*fl[17]*nul[99]-1.452368754827781*fl[8]*nul[99]-1.875*fl[39]*nul[96]-1.452368754827781*fl[17]*nul[96]-0.8385254915624212*fl[8]*nul[96]; 
  incr2[40] = (-4.192627457812106*fl[40]*nul[109])-3.247595264191645*fl[18]*nul[109]-1.875*fl[9]*nul[109]-3.247595264191643*fl[40]*nul[99]-2.515576474687264*fl[18]*nul[99]-1.452368754827781*fl[9]*nul[99]-1.875*fl[40]*nul[96]-1.452368754827781*fl[18]*nul[96]-0.8385254915624212*fl[9]*nul[96]; 
  incr2[42] = 0.8385254915624212*fl[42]*nul[109]+0.4841229182759271*fl[28]*nul[109]+0.6495190528383289*fl[42]*nul[99]+0.375*fl[28]*nul[99]+0.375*fl[42]*nul[96]+0.2165063509461097*fl[28]*nul[96]; 
  incr2[43] = 0.8385254915624212*fl[43]*nul[109]+0.4841229182759271*fl[29]*nul[109]+0.6495190528383289*fl[43]*nul[99]+0.375*fl[29]*nul[99]+0.375*fl[43]*nul[96]+0.2165063509461097*fl[29]*nul[96]; 
  incr2[44] = 0.8385254915624212*fl[44]*nul[109]+0.4841229182759271*fl[35]*nul[109]+0.6495190528383289*fl[44]*nul[99]+0.375*fl[35]*nul[99]+0.375*fl[44]*nul[96]+0.2165063509461097*fl[35]*nul[96]; 
  incr2[45] = 0.8385254915624212*fl[45]*nul[109]+0.4841229182759271*fl[36]*nul[109]+0.6495190528383289*fl[45]*nul[99]+0.375*fl[36]*nul[99]+0.375*fl[45]*nul[96]+0.2165063509461097*fl[36]*nul[96]; 
  incr2[46] = (-4.192627457812106*fl[46]*nul[109])-3.247595264191645*fl[31]*nul[109]-1.875*fl[16]*nul[109]-3.247595264191643*fl[46]*nul[99]-2.515576474687263*fl[31]*nul[99]-1.452368754827781*fl[16]*nul[99]-1.875*fl[46]*nul[96]-1.452368754827781*fl[31]*nul[96]-0.8385254915624211*fl[16]*nul[96]; 
  incr2[47] = 0.8385254915624212*fl[47]*nul[109]+0.4841229182759271*fl[41]*nul[109]+0.6495190528383289*fl[47]*nul[99]+0.375*fl[41]*nul[99]+0.375*fl[47]*nul[96]+0.2165063509461097*fl[41]*nul[96]; 

  outl[3] += -1.0*incr2[3]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[13] += incr2[13]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 
  outl[17] += -1.0*incr2[17]*rdxFl; 
  outl[18] += -1.0*incr2[18]*rdxFl; 
  outl[21] += -1.0*incr2[21]*rdxFl; 
  outl[22] += -1.0*incr2[22]*rdxFl; 
  outl[23] += incr2[23]*rdxFl; 
  outl[24] += incr2[24]*rdxFl; 
  outl[27] += incr2[27]*rdxFl; 
  outl[30] += -1.0*incr2[30]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 
  outl[32] += -1.0*incr2[32]*rdxFl; 
  outl[33] += -1.0*incr2[33]*rdxFl; 
  outl[34] += incr2[34]*rdxFl; 
  outl[37] += -1.0*incr2[37]*rdxFl; 
  outl[38] += -1.0*incr2[38]*rdxFl; 
  outl[39] += incr2[39]*rdxFl; 
  outl[40] += incr2[40]*rdxFl; 
  outl[42] += -1.0*incr2[42]*rdxFl; 
  outl[43] += -1.0*incr2[43]*rdxFl; 
  outl[44] += -1.0*incr2[44]*rdxFl; 
  outl[45] += -1.0*incr2[45]*rdxFl; 
  outl[46] += incr2[46]*rdxFl; 
  outl[47] += -1.0*incr2[47]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf4xSerP2_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:      Cell-center coordinates.
  // dx[4]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[192]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[3]*dxl[3]); 
  double rdxFr = 4.0/(dxr[3]*dxr[3]); 

  double incr1[48]; 
  double incr2[48]; 

  if (edge < 0) {

  incr2[4] = 1.082531754730548*fr[14]*nul[158]-0.8385254915624212*fr[4]*nul[158]+0.4841229182759271*fr[0]*nul[158]+0.8385254915624212*fr[14]*nul[148]-0.6495190528383289*fr[4]*nul[148]+0.375*fr[0]*nul[148]+0.4841229182759271*fr[14]*nul[144]-0.375*fr[4]*nul[144]+0.2165063509461096*fr[0]*nul[144]; 
  incr2[8] = 1.082531754730548*fr[28]*nul[158]-0.8385254915624212*fr[8]*nul[158]+0.4841229182759271*fr[1]*nul[158]+0.838525491562421*fr[28]*nul[148]-0.6495190528383289*fr[8]*nul[148]+0.375*fr[1]*nul[148]+0.4841229182759271*fr[28]*nul[144]-0.375*fr[8]*nul[144]+0.2165063509461096*fr[1]*nul[144]; 
  incr2[9] = 1.082531754730548*fr[29]*nul[158]-0.8385254915624212*fr[9]*nul[158]+0.4841229182759271*fr[2]*nul[158]+0.838525491562421*fr[29]*nul[148]-0.6495190528383289*fr[9]*nul[148]+0.375*fr[2]*nul[148]+0.4841229182759271*fr[29]*nul[144]-0.375*fr[9]*nul[144]+0.2165063509461096*fr[2]*nul[144]; 
  incr2[10] = 1.082531754730548*fr[30]*nul[158]-0.8385254915624212*fr[10]*nul[158]+0.4841229182759271*fr[3]*nul[158]+0.838525491562421*fr[30]*nul[148]-0.6495190528383289*fr[10]*nul[148]+0.375*fr[3]*nul[148]+0.4841229182759271*fr[30]*nul[144]-0.375*fr[10]*nul[144]+0.2165063509461096*fr[3]*nul[144]; 
  incr2[14] = (-4.192627457812106*fr[14]*nul[158])+3.247595264191645*fr[4]*nul[158]-1.875*fr[0]*nul[158]-3.247595264191643*fr[14]*nul[148]+2.515576474687264*fr[4]*nul[148]-1.452368754827781*fr[0]*nul[148]-1.875*fr[14]*nul[144]+1.452368754827781*fr[4]*nul[144]-0.8385254915624212*fr[0]*nul[144]; 
  incr2[16] = 1.082531754730548*fr[41]*nul[158]-0.8385254915624212*fr[16]*nul[158]+0.4841229182759271*fr[5]*nul[158]+0.8385254915624212*fr[41]*nul[148]-0.6495190528383289*fr[16]*nul[148]+0.375*fr[5]*nul[148]+0.4841229182759271*fr[41]*nul[144]-0.375*fr[16]*nul[144]+0.2165063509461096*fr[5]*nul[144]; 
  incr2[17] = 1.082531754730548*fr[42]*nul[158]-0.8385254915624212*fr[17]*nul[158]+0.4841229182759271*fr[6]*nul[158]+0.8385254915624212*fr[42]*nul[148]-0.6495190528383289*fr[17]*nul[148]+0.375*fr[6]*nul[148]+0.4841229182759271*fr[42]*nul[144]-0.375*fr[17]*nul[144]+0.2165063509461096*fr[6]*nul[144]; 
  incr2[18] = 1.082531754730548*fr[43]*nul[158]-0.8385254915624212*fr[18]*nul[158]+0.4841229182759271*fr[7]*nul[158]+0.8385254915624212*fr[43]*nul[148]-0.6495190528383289*fr[18]*nul[148]+0.375*fr[7]*nul[148]+0.4841229182759271*fr[43]*nul[144]-0.375*fr[18]*nul[144]+0.2165063509461096*fr[7]*nul[144]; 
  incr2[25] = (-0.8385254915624212*fr[25]*nul[158])+0.4841229182759271*fr[11]*nul[158]-0.6495190528383289*fr[25]*nul[148]+0.375*fr[11]*nul[148]-0.375*fr[25]*nul[144]+0.2165063509461097*fr[11]*nul[144]; 
  incr2[26] = (-0.8385254915624212*fr[26]*nul[158])+0.4841229182759271*fr[12]*nul[158]-0.6495190528383289*fr[26]*nul[148]+0.375*fr[12]*nul[148]-0.375*fr[26]*nul[144]+0.2165063509461097*fr[12]*nul[144]; 
  incr2[27] = (-0.8385254915624212*fr[27]*nul[158])+0.4841229182759271*fr[13]*nul[158]-0.6495190528383289*fr[27]*nul[148]+0.375*fr[13]*nul[148]-0.375*fr[27]*nul[144]+0.2165063509461097*fr[13]*nul[144]; 
  incr2[28] = (-4.192627457812106*fr[28]*nul[158])+3.247595264191645*fr[8]*nul[158]-1.875*fr[1]*nul[158]-3.247595264191643*fr[28]*nul[148]+2.515576474687263*fr[8]*nul[148]-1.452368754827781*fr[1]*nul[148]-1.875*fr[28]*nul[144]+1.452368754827781*fr[8]*nul[144]-0.8385254915624211*fr[1]*nul[144]; 
  incr2[29] = (-4.192627457812106*fr[29]*nul[158])+3.247595264191645*fr[9]*nul[158]-1.875*fr[2]*nul[158]-3.247595264191643*fr[29]*nul[148]+2.515576474687263*fr[9]*nul[148]-1.452368754827781*fr[2]*nul[148]-1.875*fr[29]*nul[144]+1.452368754827781*fr[9]*nul[144]-0.8385254915624211*fr[2]*nul[144]; 
  incr2[30] = (-4.192627457812106*fr[30]*nul[158])+3.247595264191645*fr[10]*nul[158]-1.875*fr[3]*nul[158]-3.247595264191643*fr[30]*nul[148]+2.515576474687263*fr[10]*nul[148]-1.452368754827781*fr[3]*nul[148]-1.875*fr[30]*nul[144]+1.452368754827781*fr[10]*nul[144]-0.8385254915624211*fr[3]*nul[144]; 
  incr2[31] = 1.082531754730548*fr[47]*nul[158]-0.8385254915624212*fr[31]*nul[158]+0.4841229182759271*fr[15]*nul[158]+0.838525491562421*fr[47]*nul[148]-0.6495190528383289*fr[31]*nul[148]+0.375*fr[15]*nul[148]+0.4841229182759271*fr[47]*nul[144]-0.375*fr[31]*nul[144]+0.2165063509461096*fr[15]*nul[144]; 
  incr2[35] = (-0.8385254915624212*fr[35]*nul[158])+0.4841229182759271*fr[19]*nul[158]-0.6495190528383289*fr[35]*nul[148]+0.375*fr[19]*nul[148]-0.375*fr[35]*nul[144]+0.2165063509461097*fr[19]*nul[144]; 
  incr2[36] = (-0.8385254915624212*fr[36]*nul[158])+0.4841229182759271*fr[20]*nul[158]-0.6495190528383289*fr[36]*nul[148]+0.375*fr[20]*nul[148]-0.375*fr[36]*nul[144]+0.2165063509461097*fr[20]*nul[144]; 
  incr2[37] = (-0.8385254915624212*fr[37]*nul[158])+0.4841229182759271*fr[21]*nul[158]-0.6495190528383289*fr[37]*nul[148]+0.375*fr[21]*nul[148]-0.375*fr[37]*nul[144]+0.2165063509461097*fr[21]*nul[144]; 
  incr2[38] = (-0.8385254915624212*fr[38]*nul[158])+0.4841229182759271*fr[22]*nul[158]-0.6495190528383289*fr[38]*nul[148]+0.375*fr[22]*nul[148]-0.375*fr[38]*nul[144]+0.2165063509461097*fr[22]*nul[144]; 
  incr2[39] = (-0.8385254915624212*fr[39]*nul[158])+0.4841229182759271*fr[23]*nul[158]-0.6495190528383289*fr[39]*nul[148]+0.375*fr[23]*nul[148]-0.375*fr[39]*nul[144]+0.2165063509461097*fr[23]*nul[144]; 
  incr2[40] = (-0.8385254915624212*fr[40]*nul[158])+0.4841229182759271*fr[24]*nul[158]-0.6495190528383289*fr[40]*nul[148]+0.375*fr[24]*nul[148]-0.375*fr[40]*nul[144]+0.2165063509461097*fr[24]*nul[144]; 
  incr2[41] = (-4.192627457812106*fr[41]*nul[158])+3.247595264191645*fr[16]*nul[158]-1.875*fr[5]*nul[158]-3.247595264191643*fr[41]*nul[148]+2.515576474687264*fr[16]*nul[148]-1.452368754827781*fr[5]*nul[148]-1.875*fr[41]*nul[144]+1.452368754827781*fr[16]*nul[144]-0.8385254915624212*fr[5]*nul[144]; 
  incr2[42] = (-4.192627457812106*fr[42]*nul[158])+3.247595264191645*fr[17]*nul[158]-1.875*fr[6]*nul[158]-3.247595264191643*fr[42]*nul[148]+2.515576474687264*fr[17]*nul[148]-1.452368754827781*fr[6]*nul[148]-1.875*fr[42]*nul[144]+1.452368754827781*fr[17]*nul[144]-0.8385254915624212*fr[6]*nul[144]; 
  incr2[43] = (-4.192627457812106*fr[43]*nul[158])+3.247595264191645*fr[18]*nul[158]-1.875*fr[7]*nul[158]-3.247595264191643*fr[43]*nul[148]+2.515576474687264*fr[18]*nul[148]-1.452368754827781*fr[7]*nul[148]-1.875*fr[43]*nul[144]+1.452368754827781*fr[18]*nul[144]-0.8385254915624212*fr[7]*nul[144]; 
  incr2[44] = (-0.8385254915624212*fr[44]*nul[158])+0.4841229182759271*fr[32]*nul[158]-0.6495190528383289*fr[44]*nul[148]+0.375*fr[32]*nul[148]-0.375*fr[44]*nul[144]+0.2165063509461097*fr[32]*nul[144]; 
  incr2[45] = (-0.8385254915624212*fr[45]*nul[158])+0.4841229182759271*fr[33]*nul[158]-0.6495190528383289*fr[45]*nul[148]+0.375*fr[33]*nul[148]-0.375*fr[45]*nul[144]+0.2165063509461097*fr[33]*nul[144]; 
  incr2[46] = (-0.8385254915624212*fr[46]*nul[158])+0.4841229182759271*fr[34]*nul[158]-0.6495190528383289*fr[46]*nul[148]+0.375*fr[34]*nul[148]-0.375*fr[46]*nul[144]+0.2165063509461097*fr[34]*nul[144]; 
  incr2[47] = (-4.192627457812106*fr[47]*nul[158])+3.247595264191645*fr[31]*nul[158]-1.875*fr[15]*nul[158]-3.247595264191643*fr[47]*nul[148]+2.515576474687263*fr[31]*nul[148]-1.452368754827781*fr[15]*nul[148]-1.875*fr[47]*nul[144]+1.452368754827781*fr[31]*nul[144]-0.8385254915624211*fr[15]*nul[144]; 

  outr[4] += incr2[4]*rdxFr; 
  outr[8] += incr2[8]*rdxFr; 
  outr[9] += incr2[9]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[14] += incr2[14]*rdxFr; 
  outr[16] += incr2[16]*rdxFr; 
  outr[17] += incr2[17]*rdxFr; 
  outr[18] += incr2[18]*rdxFr; 
  outr[25] += incr2[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 
  outr[27] += incr2[27]*rdxFr; 
  outr[28] += incr2[28]*rdxFr; 
  outr[29] += incr2[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 
  outr[35] += incr2[35]*rdxFr; 
  outr[36] += incr2[36]*rdxFr; 
  outr[37] += incr2[37]*rdxFr; 
  outr[38] += incr2[38]*rdxFr; 
  outr[39] += incr2[39]*rdxFr; 
  outr[40] += incr2[40]*rdxFr; 
  outr[41] += incr2[41]*rdxFr; 
  outr[42] += incr2[42]*rdxFr; 
  outr[43] += incr2[43]*rdxFr; 
  outr[44] += incr2[44]*rdxFr; 
  outr[45] += incr2[45]*rdxFr; 
  outr[46] += incr2[46]*rdxFr; 
  outr[47] += incr2[47]*rdxFr; 

  } else {

  incr2[4] = 1.082531754730548*fl[14]*nul[158]+0.8385254915624212*fl[4]*nul[158]+0.4841229182759271*fl[0]*nul[158]+0.8385254915624212*fl[14]*nul[148]+0.6495190528383289*fl[4]*nul[148]+0.375*fl[0]*nul[148]+0.4841229182759271*fl[14]*nul[144]+0.375*fl[4]*nul[144]+0.2165063509461096*fl[0]*nul[144]; 
  incr2[8] = 1.082531754730548*fl[28]*nul[158]+0.8385254915624212*fl[8]*nul[158]+0.4841229182759271*fl[1]*nul[158]+0.838525491562421*fl[28]*nul[148]+0.6495190528383289*fl[8]*nul[148]+0.375*fl[1]*nul[148]+0.4841229182759271*fl[28]*nul[144]+0.375*fl[8]*nul[144]+0.2165063509461096*fl[1]*nul[144]; 
  incr2[9] = 1.082531754730548*fl[29]*nul[158]+0.8385254915624212*fl[9]*nul[158]+0.4841229182759271*fl[2]*nul[158]+0.838525491562421*fl[29]*nul[148]+0.6495190528383289*fl[9]*nul[148]+0.375*fl[2]*nul[148]+0.4841229182759271*fl[29]*nul[144]+0.375*fl[9]*nul[144]+0.2165063509461096*fl[2]*nul[144]; 
  incr2[10] = 1.082531754730548*fl[30]*nul[158]+0.8385254915624212*fl[10]*nul[158]+0.4841229182759271*fl[3]*nul[158]+0.838525491562421*fl[30]*nul[148]+0.6495190528383289*fl[10]*nul[148]+0.375*fl[3]*nul[148]+0.4841229182759271*fl[30]*nul[144]+0.375*fl[10]*nul[144]+0.2165063509461096*fl[3]*nul[144]; 
  incr2[14] = (-4.192627457812106*fl[14]*nul[158])-3.247595264191645*fl[4]*nul[158]-1.875*fl[0]*nul[158]-3.247595264191643*fl[14]*nul[148]-2.515576474687264*fl[4]*nul[148]-1.452368754827781*fl[0]*nul[148]-1.875*fl[14]*nul[144]-1.452368754827781*fl[4]*nul[144]-0.8385254915624212*fl[0]*nul[144]; 
  incr2[16] = 1.082531754730548*fl[41]*nul[158]+0.8385254915624212*fl[16]*nul[158]+0.4841229182759271*fl[5]*nul[158]+0.8385254915624212*fl[41]*nul[148]+0.6495190528383289*fl[16]*nul[148]+0.375*fl[5]*nul[148]+0.4841229182759271*fl[41]*nul[144]+0.375*fl[16]*nul[144]+0.2165063509461096*fl[5]*nul[144]; 
  incr2[17] = 1.082531754730548*fl[42]*nul[158]+0.8385254915624212*fl[17]*nul[158]+0.4841229182759271*fl[6]*nul[158]+0.8385254915624212*fl[42]*nul[148]+0.6495190528383289*fl[17]*nul[148]+0.375*fl[6]*nul[148]+0.4841229182759271*fl[42]*nul[144]+0.375*fl[17]*nul[144]+0.2165063509461096*fl[6]*nul[144]; 
  incr2[18] = 1.082531754730548*fl[43]*nul[158]+0.8385254915624212*fl[18]*nul[158]+0.4841229182759271*fl[7]*nul[158]+0.8385254915624212*fl[43]*nul[148]+0.6495190528383289*fl[18]*nul[148]+0.375*fl[7]*nul[148]+0.4841229182759271*fl[43]*nul[144]+0.375*fl[18]*nul[144]+0.2165063509461096*fl[7]*nul[144]; 
  incr2[25] = 0.8385254915624212*fl[25]*nul[158]+0.4841229182759271*fl[11]*nul[158]+0.6495190528383289*fl[25]*nul[148]+0.375*fl[11]*nul[148]+0.375*fl[25]*nul[144]+0.2165063509461097*fl[11]*nul[144]; 
  incr2[26] = 0.8385254915624212*fl[26]*nul[158]+0.4841229182759271*fl[12]*nul[158]+0.6495190528383289*fl[26]*nul[148]+0.375*fl[12]*nul[148]+0.375*fl[26]*nul[144]+0.2165063509461097*fl[12]*nul[144]; 
  incr2[27] = 0.8385254915624212*fl[27]*nul[158]+0.4841229182759271*fl[13]*nul[158]+0.6495190528383289*fl[27]*nul[148]+0.375*fl[13]*nul[148]+0.375*fl[27]*nul[144]+0.2165063509461097*fl[13]*nul[144]; 
  incr2[28] = (-4.192627457812106*fl[28]*nul[158])-3.247595264191645*fl[8]*nul[158]-1.875*fl[1]*nul[158]-3.247595264191643*fl[28]*nul[148]-2.515576474687263*fl[8]*nul[148]-1.452368754827781*fl[1]*nul[148]-1.875*fl[28]*nul[144]-1.452368754827781*fl[8]*nul[144]-0.8385254915624211*fl[1]*nul[144]; 
  incr2[29] = (-4.192627457812106*fl[29]*nul[158])-3.247595264191645*fl[9]*nul[158]-1.875*fl[2]*nul[158]-3.247595264191643*fl[29]*nul[148]-2.515576474687263*fl[9]*nul[148]-1.452368754827781*fl[2]*nul[148]-1.875*fl[29]*nul[144]-1.452368754827781*fl[9]*nul[144]-0.8385254915624211*fl[2]*nul[144]; 
  incr2[30] = (-4.192627457812106*fl[30]*nul[158])-3.247595264191645*fl[10]*nul[158]-1.875*fl[3]*nul[158]-3.247595264191643*fl[30]*nul[148]-2.515576474687263*fl[10]*nul[148]-1.452368754827781*fl[3]*nul[148]-1.875*fl[30]*nul[144]-1.452368754827781*fl[10]*nul[144]-0.8385254915624211*fl[3]*nul[144]; 
  incr2[31] = 1.082531754730548*fl[47]*nul[158]+0.8385254915624212*fl[31]*nul[158]+0.4841229182759271*fl[15]*nul[158]+0.838525491562421*fl[47]*nul[148]+0.6495190528383289*fl[31]*nul[148]+0.375*fl[15]*nul[148]+0.4841229182759271*fl[47]*nul[144]+0.375*fl[31]*nul[144]+0.2165063509461096*fl[15]*nul[144]; 
  incr2[35] = 0.8385254915624212*fl[35]*nul[158]+0.4841229182759271*fl[19]*nul[158]+0.6495190528383289*fl[35]*nul[148]+0.375*fl[19]*nul[148]+0.375*fl[35]*nul[144]+0.2165063509461097*fl[19]*nul[144]; 
  incr2[36] = 0.8385254915624212*fl[36]*nul[158]+0.4841229182759271*fl[20]*nul[158]+0.6495190528383289*fl[36]*nul[148]+0.375*fl[20]*nul[148]+0.375*fl[36]*nul[144]+0.2165063509461097*fl[20]*nul[144]; 
  incr2[37] = 0.8385254915624212*fl[37]*nul[158]+0.4841229182759271*fl[21]*nul[158]+0.6495190528383289*fl[37]*nul[148]+0.375*fl[21]*nul[148]+0.375*fl[37]*nul[144]+0.2165063509461097*fl[21]*nul[144]; 
  incr2[38] = 0.8385254915624212*fl[38]*nul[158]+0.4841229182759271*fl[22]*nul[158]+0.6495190528383289*fl[38]*nul[148]+0.375*fl[22]*nul[148]+0.375*fl[38]*nul[144]+0.2165063509461097*fl[22]*nul[144]; 
  incr2[39] = 0.8385254915624212*fl[39]*nul[158]+0.4841229182759271*fl[23]*nul[158]+0.6495190528383289*fl[39]*nul[148]+0.375*fl[23]*nul[148]+0.375*fl[39]*nul[144]+0.2165063509461097*fl[23]*nul[144]; 
  incr2[40] = 0.8385254915624212*fl[40]*nul[158]+0.4841229182759271*fl[24]*nul[158]+0.6495190528383289*fl[40]*nul[148]+0.375*fl[24]*nul[148]+0.375*fl[40]*nul[144]+0.2165063509461097*fl[24]*nul[144]; 
  incr2[41] = (-4.192627457812106*fl[41]*nul[158])-3.247595264191645*fl[16]*nul[158]-1.875*fl[5]*nul[158]-3.247595264191643*fl[41]*nul[148]-2.515576474687264*fl[16]*nul[148]-1.452368754827781*fl[5]*nul[148]-1.875*fl[41]*nul[144]-1.452368754827781*fl[16]*nul[144]-0.8385254915624212*fl[5]*nul[144]; 
  incr2[42] = (-4.192627457812106*fl[42]*nul[158])-3.247595264191645*fl[17]*nul[158]-1.875*fl[6]*nul[158]-3.247595264191643*fl[42]*nul[148]-2.515576474687264*fl[17]*nul[148]-1.452368754827781*fl[6]*nul[148]-1.875*fl[42]*nul[144]-1.452368754827781*fl[17]*nul[144]-0.8385254915624212*fl[6]*nul[144]; 
  incr2[43] = (-4.192627457812106*fl[43]*nul[158])-3.247595264191645*fl[18]*nul[158]-1.875*fl[7]*nul[158]-3.247595264191643*fl[43]*nul[148]-2.515576474687264*fl[18]*nul[148]-1.452368754827781*fl[7]*nul[148]-1.875*fl[43]*nul[144]-1.452368754827781*fl[18]*nul[144]-0.8385254915624212*fl[7]*nul[144]; 
  incr2[44] = 0.8385254915624212*fl[44]*nul[158]+0.4841229182759271*fl[32]*nul[158]+0.6495190528383289*fl[44]*nul[148]+0.375*fl[32]*nul[148]+0.375*fl[44]*nul[144]+0.2165063509461097*fl[32]*nul[144]; 
  incr2[45] = 0.8385254915624212*fl[45]*nul[158]+0.4841229182759271*fl[33]*nul[158]+0.6495190528383289*fl[45]*nul[148]+0.375*fl[33]*nul[148]+0.375*fl[45]*nul[144]+0.2165063509461097*fl[33]*nul[144]; 
  incr2[46] = 0.8385254915624212*fl[46]*nul[158]+0.4841229182759271*fl[34]*nul[158]+0.6495190528383289*fl[46]*nul[148]+0.375*fl[34]*nul[148]+0.375*fl[46]*nul[144]+0.2165063509461097*fl[34]*nul[144]; 
  incr2[47] = (-4.192627457812106*fl[47]*nul[158])-3.247595264191645*fl[31]*nul[158]-1.875*fl[15]*nul[158]-3.247595264191643*fl[47]*nul[148]-2.515576474687263*fl[31]*nul[148]-1.452368754827781*fl[15]*nul[148]-1.875*fl[47]*nul[144]-1.452368754827781*fl[31]*nul[144]-0.8385254915624211*fl[15]*nul[144]; 

  outl[4] += -1.0*incr2[4]*rdxFl; 
  outl[8] += -1.0*incr2[8]*rdxFl; 
  outl[9] += -1.0*incr2[9]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[14] += incr2[14]*rdxFl; 
  outl[16] += -1.0*incr2[16]*rdxFl; 
  outl[17] += -1.0*incr2[17]*rdxFl; 
  outl[18] += -1.0*incr2[18]*rdxFl; 
  outl[25] += -1.0*incr2[25]*rdxFl; 
  outl[26] += -1.0*incr2[26]*rdxFl; 
  outl[27] += -1.0*incr2[27]*rdxFl; 
  outl[28] += incr2[28]*rdxFl; 
  outl[29] += incr2[29]*rdxFl; 
  outl[30] += incr2[30]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 
  outl[35] += -1.0*incr2[35]*rdxFl; 
  outl[36] += -1.0*incr2[36]*rdxFl; 
  outl[37] += -1.0*incr2[37]*rdxFl; 
  outl[38] += -1.0*incr2[38]*rdxFl; 
  outl[39] += -1.0*incr2[39]*rdxFl; 
  outl[40] += -1.0*incr2[40]*rdxFl; 
  outl[41] += incr2[41]*rdxFl; 
  outl[42] += incr2[42]*rdxFl; 
  outl[43] += incr2[43]*rdxFl; 
  outl[44] += -1.0*incr2[44]*rdxFl; 
  outl[45] += -1.0*incr2[45]*rdxFl; 
  outl[46] += -1.0*incr2[46]*rdxFl; 
  outl[47] += incr2[47]*rdxFl; 

  }

} 
