#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf5xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[0] == 1) {

  incr2[1] = 0.8660254037844386*fr[0]-1.5*fr[1]; 
  incr2[6] = 0.8660254037844386*fr[2]-1.5*fr[6]; 
  incr2[7] = 0.8660254037844386*fr[3]-1.5*fr[7]; 
  incr2[9] = 0.8660254037844386*fr[4]-1.5*fr[9]; 
  incr2[12] = 0.8660254037844386*fr[5]-1.5*fr[12]; 
  incr2[16] = 0.8660254037844386*fr[8]-1.5*fr[16]; 
  incr2[17] = 0.8660254037844386*fr[10]-1.5*fr[17]; 
  incr2[18] = 0.8660254037844386*fr[11]-1.5*fr[18]; 
  incr2[20] = 0.8660254037844386*fr[13]-1.5*fr[20]; 
  incr2[21] = 0.8660254037844386*fr[14]-1.5*fr[21]; 
  incr2[23] = 0.8660254037844386*fr[15]-1.5*fr[23]; 
  incr2[26] = 0.8660254037844386*fr[19]-1.5*fr[26]; 
  incr2[27] = 0.8660254037844386*fr[22]-1.5*fr[27]; 
  incr2[28] = 0.8660254037844386*fr[24]-1.5*fr[28]; 
  incr2[29] = 0.8660254037844386*fr[25]-1.5*fr[29]; 
  incr2[31] = 0.8660254037844386*fr[30]-1.5*fr[31]; 

  outr[1] += incr2[1]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 

  } else {

  incr2[1] = 1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844386*fl[2]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844386*fl[3]; 
  incr2[9] = 1.5*fl[9]+0.8660254037844386*fl[4]; 
  incr2[12] = 1.5*fl[12]+0.8660254037844386*fl[5]; 
  incr2[16] = 1.5*fl[16]+0.8660254037844386*fl[8]; 
  incr2[17] = 1.5*fl[17]+0.8660254037844386*fl[10]; 
  incr2[18] = 1.5*fl[18]+0.8660254037844386*fl[11]; 
  incr2[20] = 1.5*fl[20]+0.8660254037844386*fl[13]; 
  incr2[21] = 1.5*fl[21]+0.8660254037844386*fl[14]; 
  incr2[23] = 1.5*fl[23]+0.8660254037844386*fl[15]; 
  incr2[26] = 1.5*fl[26]+0.8660254037844386*fl[19]; 
  incr2[27] = 1.5*fl[27]+0.8660254037844386*fl[22]; 
  incr2[28] = 1.5*fl[28]+0.8660254037844386*fl[24]; 
  incr2[29] = 1.5*fl[29]+0.8660254037844386*fl[25]; 
  incr2[31] = 1.5*fl[31]+0.8660254037844386*fl[30]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf5xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[1] == 1) {

  incr2[2] = 0.8660254037844386*fr[0]-1.5*fr[2]; 
  incr2[6] = 0.8660254037844386*fr[1]-1.5*fr[6]; 
  incr2[8] = 0.8660254037844386*fr[3]-1.5*fr[8]; 
  incr2[10] = 0.8660254037844386*fr[4]-1.5*fr[10]; 
  incr2[13] = 0.8660254037844386*fr[5]-1.5*fr[13]; 
  incr2[16] = 0.8660254037844386*fr[7]-1.5*fr[16]; 
  incr2[17] = 0.8660254037844386*fr[9]-1.5*fr[17]; 
  incr2[19] = 0.8660254037844386*fr[11]-1.5*fr[19]; 
  incr2[20] = 0.8660254037844386*fr[12]-1.5*fr[20]; 
  incr2[22] = 0.8660254037844386*fr[14]-1.5*fr[22]; 
  incr2[24] = 0.8660254037844386*fr[15]-1.5*fr[24]; 
  incr2[26] = 0.8660254037844386*fr[18]-1.5*fr[26]; 
  incr2[27] = 0.8660254037844386*fr[21]-1.5*fr[27]; 
  incr2[28] = 0.8660254037844386*fr[23]-1.5*fr[28]; 
  incr2[30] = 0.8660254037844386*fr[25]-1.5*fr[30]; 
  incr2[31] = 0.8660254037844386*fr[29]-1.5*fr[31]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 

  } else {

  incr2[2] = 1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[6] = 1.5*fl[6]+0.8660254037844386*fl[1]; 
  incr2[8] = 1.5*fl[8]+0.8660254037844386*fl[3]; 
  incr2[10] = 1.5*fl[10]+0.8660254037844386*fl[4]; 
  incr2[13] = 1.5*fl[13]+0.8660254037844386*fl[5]; 
  incr2[16] = 1.5*fl[16]+0.8660254037844386*fl[7]; 
  incr2[17] = 1.5*fl[17]+0.8660254037844386*fl[9]; 
  incr2[19] = 1.5*fl[19]+0.8660254037844386*fl[11]; 
  incr2[20] = 1.5*fl[20]+0.8660254037844386*fl[12]; 
  incr2[22] = 1.5*fl[22]+0.8660254037844386*fl[14]; 
  incr2[24] = 1.5*fl[24]+0.8660254037844386*fl[15]; 
  incr2[26] = 1.5*fl[26]+0.8660254037844386*fl[18]; 
  incr2[27] = 1.5*fl[27]+0.8660254037844386*fl[21]; 
  incr2[28] = 1.5*fl[28]+0.8660254037844386*fl[23]; 
  incr2[30] = 1.5*fl[30]+0.8660254037844386*fl[25]; 
  incr2[31] = 1.5*fl[31]+0.8660254037844386*fl[29]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf5xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[2] == 1) {

  incr2[3] = 0.8660254037844386*fr[0]-1.5*fr[3]; 
  incr2[7] = 0.8660254037844386*fr[1]-1.5*fr[7]; 
  incr2[8] = 0.8660254037844386*fr[2]-1.5*fr[8]; 
  incr2[11] = 0.8660254037844386*fr[4]-1.5*fr[11]; 
  incr2[14] = 0.8660254037844386*fr[5]-1.5*fr[14]; 
  incr2[16] = 0.8660254037844386*fr[6]-1.5*fr[16]; 
  incr2[18] = 0.8660254037844386*fr[9]-1.5*fr[18]; 
  incr2[19] = 0.8660254037844386*fr[10]-1.5*fr[19]; 
  incr2[21] = 0.8660254037844386*fr[12]-1.5*fr[21]; 
  incr2[22] = 0.8660254037844386*fr[13]-1.5*fr[22]; 
  incr2[25] = 0.8660254037844386*fr[15]-1.5*fr[25]; 
  incr2[26] = 0.8660254037844386*fr[17]-1.5*fr[26]; 
  incr2[27] = 0.8660254037844386*fr[20]-1.5*fr[27]; 
  incr2[29] = 0.8660254037844386*fr[23]-1.5*fr[29]; 
  incr2[30] = 0.8660254037844386*fr[24]-1.5*fr[30]; 
  incr2[31] = 0.8660254037844386*fr[28]-1.5*fr[31]; 

  outr[3] += incr2[3]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 

  } else {

  incr2[3] = 1.5*fl[3]+0.8660254037844386*fl[0]; 
  incr2[7] = 1.5*fl[7]+0.8660254037844386*fl[1]; 
  incr2[8] = 1.5*fl[8]+0.8660254037844386*fl[2]; 
  incr2[11] = 1.5*fl[11]+0.8660254037844386*fl[4]; 
  incr2[14] = 1.5*fl[14]+0.8660254037844386*fl[5]; 
  incr2[16] = 1.5*fl[16]+0.8660254037844386*fl[6]; 
  incr2[18] = 1.5*fl[18]+0.8660254037844386*fl[9]; 
  incr2[19] = 1.5*fl[19]+0.8660254037844386*fl[10]; 
  incr2[21] = 1.5*fl[21]+0.8660254037844386*fl[12]; 
  incr2[22] = 1.5*fl[22]+0.8660254037844386*fl[13]; 
  incr2[25] = 1.5*fl[25]+0.8660254037844386*fl[15]; 
  incr2[26] = 1.5*fl[26]+0.8660254037844386*fl[17]; 
  incr2[27] = 1.5*fl[27]+0.8660254037844386*fl[20]; 
  incr2[29] = 1.5*fl[29]+0.8660254037844386*fl[23]; 
  incr2[30] = 1.5*fl[30]+0.8660254037844386*fl[24]; 
  incr2[31] = 1.5*fl[31]+0.8660254037844386*fl[28]; 

  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf5xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxFnur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[3] == 1) {

  incr2[4] = 0.8660254037844386*fr[0]-1.5*fr[4]; 
  incr2[9] = 0.8660254037844386*fr[1]-1.5*fr[9]; 
  incr2[10] = 0.8660254037844386*fr[2]-1.5*fr[10]; 
  incr2[11] = 0.8660254037844386*fr[3]-1.5*fr[11]; 
  incr2[15] = 0.8660254037844386*fr[5]-1.5*fr[15]; 
  incr2[17] = 0.8660254037844386*fr[6]-1.5*fr[17]; 
  incr2[18] = 0.8660254037844386*fr[7]-1.5*fr[18]; 
  incr2[19] = 0.8660254037844386*fr[8]-1.5*fr[19]; 
  incr2[23] = 0.8660254037844386*fr[12]-1.5*fr[23]; 
  incr2[24] = 0.8660254037844386*fr[13]-1.5*fr[24]; 
  incr2[25] = 0.8660254037844386*fr[14]-1.5*fr[25]; 
  incr2[26] = 0.8660254037844386*fr[16]-1.5*fr[26]; 
  incr2[28] = 0.8660254037844386*fr[20]-1.5*fr[28]; 
  incr2[29] = 0.8660254037844386*fr[21]-1.5*fr[29]; 
  incr2[30] = 0.8660254037844386*fr[22]-1.5*fr[30]; 
  incr2[31] = 0.8660254037844386*fr[27]-1.5*fr[31]; 

  outr[4] += incr2[4]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 

  } else {

  incr2[4] = 1.5*fl[4]+0.8660254037844386*fl[0]; 
  incr2[9] = 1.5*fl[9]+0.8660254037844386*fl[1]; 
  incr2[10] = 1.5*fl[10]+0.8660254037844386*fl[2]; 
  incr2[11] = 1.5*fl[11]+0.8660254037844386*fl[3]; 
  incr2[15] = 1.5*fl[15]+0.8660254037844386*fl[5]; 
  incr2[17] = 1.5*fl[17]+0.8660254037844386*fl[6]; 
  incr2[18] = 1.5*fl[18]+0.8660254037844386*fl[7]; 
  incr2[19] = 1.5*fl[19]+0.8660254037844386*fl[8]; 
  incr2[23] = 1.5*fl[23]+0.8660254037844386*fl[12]; 
  incr2[24] = 1.5*fl[24]+0.8660254037844386*fl[13]; 
  incr2[25] = 1.5*fl[25]+0.8660254037844386*fl[14]; 
  incr2[26] = 1.5*fl[26]+0.8660254037844386*fl[16]; 
  incr2[28] = 1.5*fl[28]+0.8660254037844386*fl[20]; 
  incr2[29] = 1.5*fl[29]+0.8660254037844386*fl[21]; 
  incr2[30] = 1.5*fl[30]+0.8660254037844386*fl[22]; 
  incr2[31] = 1.5*fl[31]+0.8660254037844386*fl[27]; 

  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf5xSerP1_X5(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[4]/(dxl[4]*dxl[4]); 
  double rdxFnur = 4.0*nu[4]/(dxr[4]*dxr[4]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[4] == 1) {

  incr2[5] = 0.8660254037844386*fr[0]-1.5*fr[5]; 
  incr2[12] = 0.8660254037844386*fr[1]-1.5*fr[12]; 
  incr2[13] = 0.8660254037844386*fr[2]-1.5*fr[13]; 
  incr2[14] = 0.8660254037844386*fr[3]-1.5*fr[14]; 
  incr2[15] = 0.8660254037844386*fr[4]-1.5*fr[15]; 
  incr2[20] = 0.8660254037844386*fr[6]-1.5*fr[20]; 
  incr2[21] = 0.8660254037844386*fr[7]-1.5*fr[21]; 
  incr2[22] = 0.8660254037844386*fr[8]-1.5*fr[22]; 
  incr2[23] = 0.8660254037844386*fr[9]-1.5*fr[23]; 
  incr2[24] = 0.8660254037844386*fr[10]-1.5*fr[24]; 
  incr2[25] = 0.8660254037844386*fr[11]-1.5*fr[25]; 
  incr2[27] = 0.8660254037844386*fr[16]-1.5*fr[27]; 
  incr2[28] = 0.8660254037844386*fr[17]-1.5*fr[28]; 
  incr2[29] = 0.8660254037844386*fr[18]-1.5*fr[29]; 
  incr2[30] = 0.8660254037844386*fr[19]-1.5*fr[30]; 
  incr2[31] = 0.8660254037844386*fr[26]-1.5*fr[31]; 

  outr[5] += incr2[5]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur; 

  } else {

  incr2[5] = 1.5*fl[5]+0.8660254037844386*fl[0]; 
  incr2[12] = 1.5*fl[12]+0.8660254037844386*fl[1]; 
  incr2[13] = 1.5*fl[13]+0.8660254037844386*fl[2]; 
  incr2[14] = 1.5*fl[14]+0.8660254037844386*fl[3]; 
  incr2[15] = 1.5*fl[15]+0.8660254037844386*fl[4]; 
  incr2[20] = 1.5*fl[20]+0.8660254037844386*fl[6]; 
  incr2[21] = 1.5*fl[21]+0.8660254037844386*fl[7]; 
  incr2[22] = 1.5*fl[22]+0.8660254037844386*fl[8]; 
  incr2[23] = 1.5*fl[23]+0.8660254037844386*fl[9]; 
  incr2[24] = 1.5*fl[24]+0.8660254037844386*fl[10]; 
  incr2[25] = 1.5*fl[25]+0.8660254037844386*fl[11]; 
  incr2[27] = 1.5*fl[27]+0.8660254037844386*fl[16]; 
  incr2[28] = 1.5*fl[28]+0.8660254037844386*fl[17]; 
  incr2[29] = 1.5*fl[29]+0.8660254037844386*fl[18]; 
  incr2[30] = 1.5*fl[30]+0.8660254037844386*fl[19]; 
  incr2[31] = 1.5*fl[31]+0.8660254037844386*fl[26]; 

  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf5xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (idxr[0] == 1) {


  incr2[1] = 1.5*fr[1]; 
  incr2[6] = 1.5*fr[6]; 
  incr2[7] = 1.5*fr[7]; 
  incr2[9] = 1.5*fr[9]; 
  incr2[12] = 1.5*fr[12]; 
  incr2[16] = 1.5*fr[16]; 
  incr2[17] = 1.5*fr[17]; 
  incr2[18] = 1.5*fr[18]; 
  incr2[20] = 1.5*fr[20]; 
  incr2[21] = 1.5*fr[21]; 
  incr2[23] = 1.5*fr[23]; 
  incr2[26] = 1.5*fr[26]; 
  incr2[27] = 1.5*fr[27]; 
  incr2[28] = 1.5*fr[28]; 
  incr2[29] = 1.5*fr[29]; 
  incr2[31] = 1.5*fr[31]; 



  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[9] += -1.0*incr2[9]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[20] += -1.0*incr2[20]*rdxFnur; 
  outr[21] += -1.0*incr2[21]*rdxFnur; 
  outr[23] += -1.0*incr2[23]*rdxFnur; 
  outr[26] += -1.0*incr2[26]*rdxFnur; 
  outr[27] += -1.0*incr2[27]*rdxFnur; 
  outr[28] += -1.0*incr2[28]*rdxFnur; 
  outr[29] += -1.0*incr2[29]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 

  } else {


  incr2[1] = 1.5*fl[1]; 
  incr2[6] = 1.5*fl[6]; 
  incr2[7] = 1.5*fl[7]; 
  incr2[9] = 1.5*fl[9]; 
  incr2[12] = 1.5*fl[12]; 
  incr2[16] = 1.5*fl[16]; 
  incr2[17] = 1.5*fl[17]; 
  incr2[18] = 1.5*fl[18]; 
  incr2[20] = 1.5*fl[20]; 
  incr2[21] = 1.5*fl[21]; 
  incr2[23] = 1.5*fl[23]; 
  incr2[26] = 1.5*fl[26]; 
  incr2[27] = 1.5*fl[27]; 
  incr2[28] = 1.5*fl[28]; 
  incr2[29] = 1.5*fl[29]; 
  incr2[31] = 1.5*fl[31]; 



  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf5xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (idxr[1] == 1) {


  incr2[2] = 1.5*fr[2]; 
  incr2[6] = 1.5*fr[6]; 
  incr2[8] = 1.5*fr[8]; 
  incr2[10] = 1.5*fr[10]; 
  incr2[13] = 1.5*fr[13]; 
  incr2[16] = 1.5*fr[16]; 
  incr2[17] = 1.5*fr[17]; 
  incr2[19] = 1.5*fr[19]; 
  incr2[20] = 1.5*fr[20]; 
  incr2[22] = 1.5*fr[22]; 
  incr2[24] = 1.5*fr[24]; 
  incr2[26] = 1.5*fr[26]; 
  incr2[27] = 1.5*fr[27]; 
  incr2[28] = 1.5*fr[28]; 
  incr2[30] = 1.5*fr[30]; 
  incr2[31] = 1.5*fr[31]; 



  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[8] += -1.0*incr2[8]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[19] += -1.0*incr2[19]*rdxFnur; 
  outr[20] += -1.0*incr2[20]*rdxFnur; 
  outr[22] += -1.0*incr2[22]*rdxFnur; 
  outr[24] += -1.0*incr2[24]*rdxFnur; 
  outr[26] += -1.0*incr2[26]*rdxFnur; 
  outr[27] += -1.0*incr2[27]*rdxFnur; 
  outr[28] += -1.0*incr2[28]*rdxFnur; 
  outr[30] += -1.0*incr2[30]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 

  } else {


  incr2[2] = 1.5*fl[2]; 
  incr2[6] = 1.5*fl[6]; 
  incr2[8] = 1.5*fl[8]; 
  incr2[10] = 1.5*fl[10]; 
  incr2[13] = 1.5*fl[13]; 
  incr2[16] = 1.5*fl[16]; 
  incr2[17] = 1.5*fl[17]; 
  incr2[19] = 1.5*fl[19]; 
  incr2[20] = 1.5*fl[20]; 
  incr2[22] = 1.5*fl[22]; 
  incr2[24] = 1.5*fl[24]; 
  incr2[26] = 1.5*fl[26]; 
  incr2[27] = 1.5*fl[27]; 
  incr2[28] = 1.5*fl[28]; 
  incr2[30] = 1.5*fl[30]; 
  incr2[31] = 1.5*fl[31]; 



  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf5xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (idxr[2] == 1) {


  incr2[3] = 1.5*fr[3]; 
  incr2[7] = 1.5*fr[7]; 
  incr2[8] = 1.5*fr[8]; 
  incr2[11] = 1.5*fr[11]; 
  incr2[14] = 1.5*fr[14]; 
  incr2[16] = 1.5*fr[16]; 
  incr2[18] = 1.5*fr[18]; 
  incr2[19] = 1.5*fr[19]; 
  incr2[21] = 1.5*fr[21]; 
  incr2[22] = 1.5*fr[22]; 
  incr2[25] = 1.5*fr[25]; 
  incr2[26] = 1.5*fr[26]; 
  incr2[27] = 1.5*fr[27]; 
  incr2[29] = 1.5*fr[29]; 
  incr2[30] = 1.5*fr[30]; 
  incr2[31] = 1.5*fr[31]; 



  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[8] += -1.0*incr2[8]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[19] += -1.0*incr2[19]*rdxFnur; 
  outr[21] += -1.0*incr2[21]*rdxFnur; 
  outr[22] += -1.0*incr2[22]*rdxFnur; 
  outr[25] += -1.0*incr2[25]*rdxFnur; 
  outr[26] += -1.0*incr2[26]*rdxFnur; 
  outr[27] += -1.0*incr2[27]*rdxFnur; 
  outr[29] += -1.0*incr2[29]*rdxFnur; 
  outr[30] += -1.0*incr2[30]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 

  } else {


  incr2[3] = 1.5*fl[3]; 
  incr2[7] = 1.5*fl[7]; 
  incr2[8] = 1.5*fl[8]; 
  incr2[11] = 1.5*fl[11]; 
  incr2[14] = 1.5*fl[14]; 
  incr2[16] = 1.5*fl[16]; 
  incr2[18] = 1.5*fl[18]; 
  incr2[19] = 1.5*fl[19]; 
  incr2[21] = 1.5*fl[21]; 
  incr2[22] = 1.5*fl[22]; 
  incr2[25] = 1.5*fl[25]; 
  incr2[26] = 1.5*fl[26]; 
  incr2[27] = 1.5*fl[27]; 
  incr2[29] = 1.5*fl[29]; 
  incr2[30] = 1.5*fl[30]; 
  incr2[31] = 1.5*fl[31]; 



  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf5xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 16.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (idxr[3] == 1) {


  incr2[4] = 1.5*fr[4]; 
  incr2[9] = 1.5*fr[9]; 
  incr2[10] = 1.5*fr[10]; 
  incr2[11] = 1.5*fr[11]; 
  incr2[15] = 1.5*fr[15]; 
  incr2[17] = 1.5*fr[17]; 
  incr2[18] = 1.5*fr[18]; 
  incr2[19] = 1.5*fr[19]; 
  incr2[23] = 1.5*fr[23]; 
  incr2[24] = 1.5*fr[24]; 
  incr2[25] = 1.5*fr[25]; 
  incr2[26] = 1.5*fr[26]; 
  incr2[28] = 1.5*fr[28]; 
  incr2[29] = 1.5*fr[29]; 
  incr2[30] = 1.5*fr[30]; 
  incr2[31] = 1.5*fr[31]; 



  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[9] += -1.0*incr2[9]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[19] += -1.0*incr2[19]*rdxFnur; 
  outr[23] += -1.0*incr2[23]*rdxFnur; 
  outr[24] += -1.0*incr2[24]*rdxFnur; 
  outr[25] += -1.0*incr2[25]*rdxFnur; 
  outr[26] += -1.0*incr2[26]*rdxFnur; 
  outr[28] += -1.0*incr2[28]*rdxFnur; 
  outr[29] += -1.0*incr2[29]*rdxFnur; 
  outr[30] += -1.0*incr2[30]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 

  } else {


  incr2[4] = 1.5*fl[4]; 
  incr2[9] = 1.5*fl[9]; 
  incr2[10] = 1.5*fl[10]; 
  incr2[11] = 1.5*fl[11]; 
  incr2[15] = 1.5*fl[15]; 
  incr2[17] = 1.5*fl[17]; 
  incr2[18] = 1.5*fl[18]; 
  incr2[19] = 1.5*fl[19]; 
  incr2[23] = 1.5*fl[23]; 
  incr2[24] = 1.5*fl[24]; 
  incr2[25] = 1.5*fl[25]; 
  incr2[26] = 1.5*fl[26]; 
  incr2[28] = 1.5*fl[28]; 
  incr2[29] = 1.5*fl[29]; 
  incr2[30] = 1.5*fl[30]; 
  incr2[31] = 1.5*fl[31]; 



  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf5xSerP1_X5(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[4]/(dxl[4]*dxl[4]*dxl[4]*dxl[4]); 
  double rdxFnur = 16.0*nu[4]/(dxr[4]*dxr[4]*dxr[4]*dxr[4]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  if (idxr[4] == 1) {


  incr2[5] = 1.5*fr[5]; 
  incr2[12] = 1.5*fr[12]; 
  incr2[13] = 1.5*fr[13]; 
  incr2[14] = 1.5*fr[14]; 
  incr2[15] = 1.5*fr[15]; 
  incr2[20] = 1.5*fr[20]; 
  incr2[21] = 1.5*fr[21]; 
  incr2[22] = 1.5*fr[22]; 
  incr2[23] = 1.5*fr[23]; 
  incr2[24] = 1.5*fr[24]; 
  incr2[25] = 1.5*fr[25]; 
  incr2[27] = 1.5*fr[27]; 
  incr2[28] = 1.5*fr[28]; 
  incr2[29] = 1.5*fr[29]; 
  incr2[30] = 1.5*fr[30]; 
  incr2[31] = 1.5*fr[31]; 



  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[20] += -1.0*incr2[20]*rdxFnur; 
  outr[21] += -1.0*incr2[21]*rdxFnur; 
  outr[22] += -1.0*incr2[22]*rdxFnur; 
  outr[23] += -1.0*incr2[23]*rdxFnur; 
  outr[24] += -1.0*incr2[24]*rdxFnur; 
  outr[25] += -1.0*incr2[25]*rdxFnur; 
  outr[27] += -1.0*incr2[27]*rdxFnur; 
  outr[28] += -1.0*incr2[28]*rdxFnur; 
  outr[29] += -1.0*incr2[29]*rdxFnur; 
  outr[30] += -1.0*incr2[30]*rdxFnur; 
  outr[31] += -1.0*incr2[31]*rdxFnur; 

  } else {


  incr2[5] = 1.5*fl[5]; 
  incr2[12] = 1.5*fl[12]; 
  incr2[13] = 1.5*fl[13]; 
  incr2[14] = 1.5*fl[14]; 
  incr2[15] = 1.5*fl[15]; 
  incr2[20] = 1.5*fl[20]; 
  incr2[21] = 1.5*fl[21]; 
  incr2[22] = 1.5*fl[22]; 
  incr2[23] = 1.5*fl[23]; 
  incr2[24] = 1.5*fl[24]; 
  incr2[25] = 1.5*fl[25]; 
  incr2[27] = 1.5*fl[27]; 
  incr2[28] = 1.5*fl[28]; 
  incr2[29] = 1.5*fl[29]; 
  incr2[30] = 1.5*fl[30]; 
  incr2[31] = 1.5*fl[31]; 



  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[27] += -1.0*incr2[27]*rdxFnul; 
  outl[28] += -1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr2[30]*rdxFnul; 
  outl[31] += -1.0*incr2[31]*rdxFnul; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf5xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[160]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[0] == 1) {

  incr2[1] = (-0.4592793267718456*fr[7]*nul[7])+0.2651650429449552*fr[3]*nul[7]-0.2651650429449552*nul[3]*fr[7]+0.1530931089239486*fr[3]*nul[3]-0.4592793267718456*fr[1]*nul[1]+0.2651650429449552*fr[0]*nul[1]-0.2651650429449552*nul[0]*fr[1]+0.1530931089239486*fr[0]*nul[0]; 
  incr2[6] = (-0.4592793267718456*nul[7]*fr[16])-0.2651650429449552*nul[3]*fr[16]+0.2651650429449552*nul[7]*fr[8]+0.1530931089239486*nul[3]*fr[8]-0.4592793267718456*nul[1]*fr[6]-0.2651650429449552*nul[0]*fr[6]+0.2651650429449552*nul[1]*fr[2]+0.1530931089239486*nul[0]*fr[2]; 
  incr2[7] = (-0.4592793267718456*fr[1]*nul[7])+0.2651650429449552*fr[0]*nul[7]-0.4592793267718456*nul[1]*fr[7]-0.2651650429449552*nul[0]*fr[7]-0.2651650429449552*fr[1]*nul[3]+0.1530931089239486*fr[0]*nul[3]+0.2651650429449552*nul[1]*fr[3]+0.1530931089239486*nul[0]*fr[3]; 
  incr2[9] = (-0.4592793267718456*nul[7]*fr[18])-0.2651650429449552*nul[3]*fr[18]+0.2651650429449552*nul[7]*fr[11]+0.1530931089239486*nul[3]*fr[11]-0.4592793267718456*nul[1]*fr[9]-0.2651650429449552*nul[0]*fr[9]+0.2651650429449552*nul[1]*fr[4]+0.1530931089239486*nul[0]*fr[4]; 
  incr2[12] = (-0.4592793267718456*nul[7]*fr[21])-0.2651650429449552*nul[3]*fr[21]+0.2651650429449552*nul[7]*fr[14]+0.1530931089239486*nul[3]*fr[14]-0.4592793267718456*nul[1]*fr[12]-0.2651650429449552*nul[0]*fr[12]+0.2651650429449552*nul[1]*fr[5]+0.1530931089239486*nul[0]*fr[5]; 
  incr2[16] = (-0.4592793267718456*nul[1]*fr[16])-0.2651650429449552*nul[0]*fr[16]+0.2651650429449552*nul[1]*fr[8]+0.1530931089239486*nul[0]*fr[8]-0.4592793267718456*fr[6]*nul[7]+0.2651650429449552*fr[2]*nul[7]-0.2651650429449552*nul[3]*fr[6]+0.1530931089239486*fr[2]*nul[3]; 
  incr2[17] = (-0.4592793267718456*nul[7]*fr[26])-0.2651650429449552*nul[3]*fr[26]+0.2651650429449552*nul[7]*fr[19]+0.1530931089239486*nul[3]*fr[19]-0.4592793267718456*nul[1]*fr[17]-0.2651650429449552*nul[0]*fr[17]+0.2651650429449552*nul[1]*fr[10]+0.1530931089239486*nul[0]*fr[10]; 
  incr2[18] = (-0.4592793267718456*nul[1]*fr[18])-0.2651650429449552*nul[0]*fr[18]+0.2651650429449552*nul[1]*fr[11]+0.1530931089239486*nul[0]*fr[11]-0.4592793267718456*nul[7]*fr[9]-0.2651650429449552*nul[3]*fr[9]+0.2651650429449552*fr[4]*nul[7]+0.1530931089239486*nul[3]*fr[4]; 
  incr2[20] = (-0.4592793267718456*nul[7]*fr[27])-0.2651650429449552*nul[3]*fr[27]+0.2651650429449552*nul[7]*fr[22]+0.1530931089239486*nul[3]*fr[22]-0.4592793267718456*nul[1]*fr[20]-0.2651650429449552*nul[0]*fr[20]+0.2651650429449552*nul[1]*fr[13]+0.1530931089239486*nul[0]*fr[13]; 
  incr2[21] = (-0.4592793267718456*nul[1]*fr[21])-0.2651650429449552*nul[0]*fr[21]+0.2651650429449552*nul[1]*fr[14]+0.1530931089239486*nul[0]*fr[14]-0.4592793267718456*nul[7]*fr[12]-0.2651650429449552*nul[3]*fr[12]+0.2651650429449552*fr[5]*nul[7]+0.1530931089239486*nul[3]*fr[5]; 
  incr2[23] = (-0.4592793267718456*nul[7]*fr[29])-0.2651650429449552*nul[3]*fr[29]+0.2651650429449552*nul[7]*fr[25]+0.1530931089239486*nul[3]*fr[25]-0.4592793267718456*nul[1]*fr[23]-0.2651650429449552*nul[0]*fr[23]+0.2651650429449552*nul[1]*fr[15]+0.1530931089239486*nul[0]*fr[15]; 
  incr2[26] = (-0.4592793267718456*nul[1]*fr[26])-0.2651650429449552*nul[0]*fr[26]+0.2651650429449552*nul[1]*fr[19]+0.1530931089239486*nul[0]*fr[19]-0.4592793267718456*nul[7]*fr[17]-0.2651650429449552*nul[3]*fr[17]+0.2651650429449552*nul[7]*fr[10]+0.1530931089239486*nul[3]*fr[10]; 
  incr2[27] = (-0.4592793267718456*nul[1]*fr[27])-0.2651650429449552*nul[0]*fr[27]+0.2651650429449552*nul[1]*fr[22]+0.1530931089239486*nul[0]*fr[22]-0.4592793267718456*nul[7]*fr[20]-0.2651650429449552*nul[3]*fr[20]+0.2651650429449552*nul[7]*fr[13]+0.1530931089239486*nul[3]*fr[13]; 
  incr2[28] = (-0.4592793267718456*nul[7]*fr[31])-0.2651650429449552*nul[3]*fr[31]+0.2651650429449552*nul[7]*fr[30]+0.1530931089239486*nul[3]*fr[30]-0.4592793267718456*nul[1]*fr[28]-0.2651650429449552*nul[0]*fr[28]+0.2651650429449552*nul[1]*fr[24]+0.1530931089239486*nul[0]*fr[24]; 
  incr2[29] = (-0.4592793267718456*nul[1]*fr[29])-0.2651650429449552*nul[0]*fr[29]+0.2651650429449552*nul[1]*fr[25]+0.1530931089239486*nul[0]*fr[25]-0.4592793267718456*nul[7]*fr[23]-0.2651650429449552*nul[3]*fr[23]+0.2651650429449552*nul[7]*fr[15]+0.1530931089239486*nul[3]*fr[15]; 
  incr2[31] = (-0.4592793267718456*nul[1]*fr[31])-0.2651650429449552*nul[0]*fr[31]+0.2651650429449552*nul[1]*fr[30]+0.1530931089239486*nul[0]*fr[30]-0.4592793267718456*nul[7]*fr[28]-0.2651650429449552*nul[3]*fr[28]+0.2651650429449552*nul[7]*fr[24]+0.1530931089239486*nul[3]*fr[24]; 

  outr[1] += incr2[1]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 
  outr[9] += incr2[9]*rdxFr; 
  outr[12] += incr2[12]*rdxFr; 
  outr[16] += incr2[16]*rdxFr; 
  outr[17] += incr2[17]*rdxFr; 
  outr[18] += incr2[18]*rdxFr; 
  outr[20] += incr2[20]*rdxFr; 
  outr[21] += incr2[21]*rdxFr; 
  outr[23] += incr2[23]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 
  outr[27] += incr2[27]*rdxFr; 
  outr[28] += incr2[28]*rdxFr; 
  outr[29] += incr2[29]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 

  } else {

  incr2[1] = 0.4592793267718456*fl[7]*nul[7]+0.2651650429449552*fl[3]*nul[7]+0.2651650429449552*nul[3]*fl[7]+0.1530931089239486*fl[3]*nul[3]+0.4592793267718456*fl[1]*nul[1]+0.2651650429449552*fl[0]*nul[1]+0.2651650429449552*nul[0]*fl[1]+0.1530931089239486*fl[0]*nul[0]; 
  incr2[6] = 0.4592793267718456*nul[7]*fl[16]+0.2651650429449552*nul[3]*fl[16]+0.2651650429449552*nul[7]*fl[8]+0.1530931089239486*nul[3]*fl[8]+0.4592793267718456*nul[1]*fl[6]+0.2651650429449552*nul[0]*fl[6]+0.2651650429449552*nul[1]*fl[2]+0.1530931089239486*nul[0]*fl[2]; 
  incr2[7] = 0.4592793267718456*fl[1]*nul[7]+0.2651650429449552*fl[0]*nul[7]+0.4592793267718456*nul[1]*fl[7]+0.2651650429449552*nul[0]*fl[7]+0.2651650429449552*fl[1]*nul[3]+0.1530931089239486*fl[0]*nul[3]+0.2651650429449552*nul[1]*fl[3]+0.1530931089239486*nul[0]*fl[3]; 
  incr2[9] = 0.4592793267718456*nul[7]*fl[18]+0.2651650429449552*nul[3]*fl[18]+0.2651650429449552*nul[7]*fl[11]+0.1530931089239486*nul[3]*fl[11]+0.4592793267718456*nul[1]*fl[9]+0.2651650429449552*nul[0]*fl[9]+0.2651650429449552*nul[1]*fl[4]+0.1530931089239486*nul[0]*fl[4]; 
  incr2[12] = 0.4592793267718456*nul[7]*fl[21]+0.2651650429449552*nul[3]*fl[21]+0.2651650429449552*nul[7]*fl[14]+0.1530931089239486*nul[3]*fl[14]+0.4592793267718456*nul[1]*fl[12]+0.2651650429449552*nul[0]*fl[12]+0.2651650429449552*nul[1]*fl[5]+0.1530931089239486*nul[0]*fl[5]; 
  incr2[16] = 0.4592793267718456*nul[1]*fl[16]+0.2651650429449552*nul[0]*fl[16]+0.2651650429449552*nul[1]*fl[8]+0.1530931089239486*nul[0]*fl[8]+0.4592793267718456*fl[6]*nul[7]+0.2651650429449552*fl[2]*nul[7]+0.2651650429449552*nul[3]*fl[6]+0.1530931089239486*fl[2]*nul[3]; 
  incr2[17] = 0.4592793267718456*nul[7]*fl[26]+0.2651650429449552*nul[3]*fl[26]+0.2651650429449552*nul[7]*fl[19]+0.1530931089239486*nul[3]*fl[19]+0.4592793267718456*nul[1]*fl[17]+0.2651650429449552*nul[0]*fl[17]+0.2651650429449552*nul[1]*fl[10]+0.1530931089239486*nul[0]*fl[10]; 
  incr2[18] = 0.4592793267718456*nul[1]*fl[18]+0.2651650429449552*nul[0]*fl[18]+0.2651650429449552*nul[1]*fl[11]+0.1530931089239486*nul[0]*fl[11]+0.4592793267718456*nul[7]*fl[9]+0.2651650429449552*nul[3]*fl[9]+0.2651650429449552*fl[4]*nul[7]+0.1530931089239486*nul[3]*fl[4]; 
  incr2[20] = 0.4592793267718456*nul[7]*fl[27]+0.2651650429449552*nul[3]*fl[27]+0.2651650429449552*nul[7]*fl[22]+0.1530931089239486*nul[3]*fl[22]+0.4592793267718456*nul[1]*fl[20]+0.2651650429449552*nul[0]*fl[20]+0.2651650429449552*nul[1]*fl[13]+0.1530931089239486*nul[0]*fl[13]; 
  incr2[21] = 0.4592793267718456*nul[1]*fl[21]+0.2651650429449552*nul[0]*fl[21]+0.2651650429449552*nul[1]*fl[14]+0.1530931089239486*nul[0]*fl[14]+0.4592793267718456*nul[7]*fl[12]+0.2651650429449552*nul[3]*fl[12]+0.2651650429449552*fl[5]*nul[7]+0.1530931089239486*nul[3]*fl[5]; 
  incr2[23] = 0.4592793267718456*nul[7]*fl[29]+0.2651650429449552*nul[3]*fl[29]+0.2651650429449552*nul[7]*fl[25]+0.1530931089239486*nul[3]*fl[25]+0.4592793267718456*nul[1]*fl[23]+0.2651650429449552*nul[0]*fl[23]+0.2651650429449552*nul[1]*fl[15]+0.1530931089239486*nul[0]*fl[15]; 
  incr2[26] = 0.4592793267718456*nul[1]*fl[26]+0.2651650429449552*nul[0]*fl[26]+0.2651650429449552*nul[1]*fl[19]+0.1530931089239486*nul[0]*fl[19]+0.4592793267718456*nul[7]*fl[17]+0.2651650429449552*nul[3]*fl[17]+0.2651650429449552*nul[7]*fl[10]+0.1530931089239486*nul[3]*fl[10]; 
  incr2[27] = 0.4592793267718456*nul[1]*fl[27]+0.2651650429449552*nul[0]*fl[27]+0.2651650429449552*nul[1]*fl[22]+0.1530931089239486*nul[0]*fl[22]+0.4592793267718456*nul[7]*fl[20]+0.2651650429449552*nul[3]*fl[20]+0.2651650429449552*nul[7]*fl[13]+0.1530931089239486*nul[3]*fl[13]; 
  incr2[28] = 0.4592793267718456*nul[7]*fl[31]+0.2651650429449552*nul[3]*fl[31]+0.2651650429449552*nul[7]*fl[30]+0.1530931089239486*nul[3]*fl[30]+0.4592793267718456*nul[1]*fl[28]+0.2651650429449552*nul[0]*fl[28]+0.2651650429449552*nul[1]*fl[24]+0.1530931089239486*nul[0]*fl[24]; 
  incr2[29] = 0.4592793267718456*nul[1]*fl[29]+0.2651650429449552*nul[0]*fl[29]+0.2651650429449552*nul[1]*fl[25]+0.1530931089239486*nul[0]*fl[25]+0.4592793267718456*nul[7]*fl[23]+0.2651650429449552*nul[3]*fl[23]+0.2651650429449552*nul[7]*fl[15]+0.1530931089239486*nul[3]*fl[15]; 
  incr2[31] = 0.4592793267718456*nul[1]*fl[31]+0.2651650429449552*nul[0]*fl[31]+0.2651650429449552*nul[1]*fl[30]+0.1530931089239486*nul[0]*fl[30]+0.4592793267718456*nul[7]*fl[28]+0.2651650429449552*nul[3]*fl[28]+0.2651650429449552*nul[7]*fl[24]+0.1530931089239486*nul[3]*fl[24]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 
  outl[9] += -1.0*incr2[9]*rdxFl; 
  outl[12] += -1.0*incr2[12]*rdxFl; 
  outl[16] += -1.0*incr2[16]*rdxFl; 
  outl[17] += -1.0*incr2[17]*rdxFl; 
  outl[18] += -1.0*incr2[18]*rdxFl; 
  outl[20] += -1.0*incr2[20]*rdxFl; 
  outl[21] += -1.0*incr2[21]*rdxFl; 
  outl[23] += -1.0*incr2[23]*rdxFl; 
  outl[26] += -1.0*incr2[26]*rdxFl; 
  outl[27] += -1.0*incr2[27]*rdxFl; 
  outl[28] += -1.0*incr2[28]*rdxFl; 
  outl[29] += -1.0*incr2[29]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf5xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[160]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[1]*dxl[1]); 
  double rdxFr = 4.0/(dxr[1]*dxr[1]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[1] == 1) {

  incr2[2] = (-0.2651650429449552*fr[16]*nul[39])+0.1530931089239486*fr[7]*nul[39]-0.2651650429449552*fr[8]*nul[35]+0.1530931089239486*fr[3]*nul[35]-0.2651650429449552*fr[6]*nul[33]+0.1530931089239486*fr[1]*nul[33]-0.2651650429449552*fr[2]*nul[32]+0.1530931089239486*fr[0]*nul[32]; 
  incr2[6] = (-0.2651650429449552*fr[8]*nul[39])+0.1530931089239486*fr[3]*nul[39]-0.2651650429449552*fr[16]*nul[35]+0.1530931089239486*fr[7]*nul[35]-0.2651650429449552*fr[2]*nul[33]+0.1530931089239486*fr[0]*nul[33]-0.2651650429449552*fr[6]*nul[32]+0.1530931089239486*fr[1]*nul[32]; 
  incr2[8] = (-0.2651650429449552*fr[6]*nul[39])+0.1530931089239486*fr[1]*nul[39]-0.2651650429449552*fr[2]*nul[35]+0.1530931089239486*fr[0]*nul[35]-0.2651650429449552*fr[16]*nul[33]+0.1530931089239486*fr[7]*nul[33]-0.2651650429449552*fr[8]*nul[32]+0.1530931089239486*fr[3]*nul[32]; 
  incr2[10] = (-0.2651650429449552*fr[26]*nul[39])+0.1530931089239486*fr[18]*nul[39]-0.2651650429449552*fr[19]*nul[35]+0.1530931089239486*fr[11]*nul[35]-0.2651650429449552*fr[17]*nul[33]+0.1530931089239486*fr[9]*nul[33]-0.2651650429449552*fr[10]*nul[32]+0.1530931089239486*fr[4]*nul[32]; 
  incr2[13] = (-0.2651650429449552*fr[27]*nul[39])+0.1530931089239486*fr[21]*nul[39]-0.2651650429449552*fr[22]*nul[35]+0.1530931089239486*fr[14]*nul[35]-0.2651650429449552*fr[20]*nul[33]+0.1530931089239486*fr[12]*nul[33]-0.2651650429449552*fr[13]*nul[32]+0.1530931089239486*fr[5]*nul[32]; 
  incr2[16] = (-0.2651650429449552*fr[2]*nul[39])+0.1530931089239486*fr[0]*nul[39]-0.2651650429449552*fr[6]*nul[35]+0.1530931089239486*fr[1]*nul[35]-0.2651650429449552*fr[8]*nul[33]+0.1530931089239486*fr[3]*nul[33]-0.2651650429449552*fr[16]*nul[32]+0.1530931089239486*fr[7]*nul[32]; 
  incr2[17] = (-0.2651650429449552*fr[19]*nul[39])+0.1530931089239486*fr[11]*nul[39]-0.2651650429449552*fr[26]*nul[35]+0.1530931089239486*fr[18]*nul[35]-0.2651650429449552*fr[10]*nul[33]+0.1530931089239486*fr[4]*nul[33]-0.2651650429449552*fr[17]*nul[32]+0.1530931089239486*fr[9]*nul[32]; 
  incr2[19] = (-0.2651650429449552*fr[17]*nul[39])+0.1530931089239486*fr[9]*nul[39]-0.2651650429449552*fr[10]*nul[35]+0.1530931089239486*fr[4]*nul[35]-0.2651650429449552*fr[26]*nul[33]+0.1530931089239486*fr[18]*nul[33]-0.2651650429449552*fr[19]*nul[32]+0.1530931089239486*fr[11]*nul[32]; 
  incr2[20] = (-0.2651650429449552*fr[22]*nul[39])+0.1530931089239486*fr[14]*nul[39]-0.2651650429449552*fr[27]*nul[35]+0.1530931089239486*fr[21]*nul[35]-0.2651650429449552*fr[13]*nul[33]+0.1530931089239486*fr[5]*nul[33]-0.2651650429449552*fr[20]*nul[32]+0.1530931089239486*fr[12]*nul[32]; 
  incr2[22] = (-0.2651650429449552*fr[20]*nul[39])+0.1530931089239486*fr[12]*nul[39]-0.2651650429449552*fr[13]*nul[35]+0.1530931089239486*fr[5]*nul[35]-0.2651650429449552*fr[27]*nul[33]+0.1530931089239486*fr[21]*nul[33]-0.2651650429449552*fr[22]*nul[32]+0.1530931089239486*fr[14]*nul[32]; 
  incr2[24] = (-0.2651650429449552*fr[31]*nul[39])+0.1530931089239486*fr[29]*nul[39]-0.2651650429449552*fr[30]*nul[35]+0.1530931089239486*fr[25]*nul[35]-0.2651650429449552*fr[28]*nul[33]+0.1530931089239486*fr[23]*nul[33]-0.2651650429449552*fr[24]*nul[32]+0.1530931089239486*fr[15]*nul[32]; 
  incr2[26] = (-0.2651650429449552*fr[10]*nul[39])+0.1530931089239486*fr[4]*nul[39]-0.2651650429449552*fr[17]*nul[35]+0.1530931089239486*fr[9]*nul[35]-0.2651650429449552*fr[19]*nul[33]+0.1530931089239486*fr[11]*nul[33]-0.2651650429449552*fr[26]*nul[32]+0.1530931089239486*fr[18]*nul[32]; 
  incr2[27] = (-0.2651650429449552*fr[13]*nul[39])+0.1530931089239486*fr[5]*nul[39]-0.2651650429449552*fr[20]*nul[35]+0.1530931089239486*fr[12]*nul[35]-0.2651650429449552*fr[22]*nul[33]+0.1530931089239486*fr[14]*nul[33]-0.2651650429449552*fr[27]*nul[32]+0.1530931089239486*fr[21]*nul[32]; 
  incr2[28] = (-0.2651650429449552*fr[30]*nul[39])+0.1530931089239486*fr[25]*nul[39]-0.2651650429449552*fr[31]*nul[35]+0.1530931089239486*fr[29]*nul[35]-0.2651650429449552*fr[24]*nul[33]+0.1530931089239486*fr[15]*nul[33]-0.2651650429449552*fr[28]*nul[32]+0.1530931089239486*fr[23]*nul[32]; 
  incr2[30] = (-0.2651650429449552*fr[28]*nul[39])+0.1530931089239486*fr[23]*nul[39]-0.2651650429449552*fr[24]*nul[35]+0.1530931089239486*fr[15]*nul[35]-0.2651650429449552*fr[31]*nul[33]+0.1530931089239486*fr[29]*nul[33]-0.2651650429449552*fr[30]*nul[32]+0.1530931089239486*fr[25]*nul[32]; 
  incr2[31] = (-0.2651650429449552*fr[24]*nul[39])+0.1530931089239486*fr[15]*nul[39]-0.2651650429449552*fr[28]*nul[35]+0.1530931089239486*fr[23]*nul[35]-0.2651650429449552*fr[30]*nul[33]+0.1530931089239486*fr[25]*nul[33]-0.2651650429449552*fr[31]*nul[32]+0.1530931089239486*fr[29]*nul[32]; 

  outr[2] += incr2[2]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[8] += incr2[8]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[13] += incr2[13]*rdxFr; 
  outr[16] += incr2[16]*rdxFr; 
  outr[17] += incr2[17]*rdxFr; 
  outr[19] += incr2[19]*rdxFr; 
  outr[20] += incr2[20]*rdxFr; 
  outr[22] += incr2[22]*rdxFr; 
  outr[24] += incr2[24]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 
  outr[27] += incr2[27]*rdxFr; 
  outr[28] += incr2[28]*rdxFr; 
  outr[30] += incr2[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 

  } else {

  incr2[2] = 0.2651650429449552*fl[16]*nul[39]+0.1530931089239486*fl[7]*nul[39]+0.2651650429449552*fl[8]*nul[35]+0.1530931089239486*fl[3]*nul[35]+0.2651650429449552*fl[6]*nul[33]+0.1530931089239486*fl[1]*nul[33]+0.2651650429449552*fl[2]*nul[32]+0.1530931089239486*fl[0]*nul[32]; 
  incr2[6] = 0.2651650429449552*fl[8]*nul[39]+0.1530931089239486*fl[3]*nul[39]+0.2651650429449552*fl[16]*nul[35]+0.1530931089239486*fl[7]*nul[35]+0.2651650429449552*fl[2]*nul[33]+0.1530931089239486*fl[0]*nul[33]+0.2651650429449552*fl[6]*nul[32]+0.1530931089239486*fl[1]*nul[32]; 
  incr2[8] = 0.2651650429449552*fl[6]*nul[39]+0.1530931089239486*fl[1]*nul[39]+0.2651650429449552*fl[2]*nul[35]+0.1530931089239486*fl[0]*nul[35]+0.2651650429449552*fl[16]*nul[33]+0.1530931089239486*fl[7]*nul[33]+0.2651650429449552*fl[8]*nul[32]+0.1530931089239486*fl[3]*nul[32]; 
  incr2[10] = 0.2651650429449552*fl[26]*nul[39]+0.1530931089239486*fl[18]*nul[39]+0.2651650429449552*fl[19]*nul[35]+0.1530931089239486*fl[11]*nul[35]+0.2651650429449552*fl[17]*nul[33]+0.1530931089239486*fl[9]*nul[33]+0.2651650429449552*fl[10]*nul[32]+0.1530931089239486*fl[4]*nul[32]; 
  incr2[13] = 0.2651650429449552*fl[27]*nul[39]+0.1530931089239486*fl[21]*nul[39]+0.2651650429449552*fl[22]*nul[35]+0.1530931089239486*fl[14]*nul[35]+0.2651650429449552*fl[20]*nul[33]+0.1530931089239486*fl[12]*nul[33]+0.2651650429449552*fl[13]*nul[32]+0.1530931089239486*fl[5]*nul[32]; 
  incr2[16] = 0.2651650429449552*fl[2]*nul[39]+0.1530931089239486*fl[0]*nul[39]+0.2651650429449552*fl[6]*nul[35]+0.1530931089239486*fl[1]*nul[35]+0.2651650429449552*fl[8]*nul[33]+0.1530931089239486*fl[3]*nul[33]+0.2651650429449552*fl[16]*nul[32]+0.1530931089239486*fl[7]*nul[32]; 
  incr2[17] = 0.2651650429449552*fl[19]*nul[39]+0.1530931089239486*fl[11]*nul[39]+0.2651650429449552*fl[26]*nul[35]+0.1530931089239486*fl[18]*nul[35]+0.2651650429449552*fl[10]*nul[33]+0.1530931089239486*fl[4]*nul[33]+0.2651650429449552*fl[17]*nul[32]+0.1530931089239486*fl[9]*nul[32]; 
  incr2[19] = 0.2651650429449552*fl[17]*nul[39]+0.1530931089239486*fl[9]*nul[39]+0.2651650429449552*fl[10]*nul[35]+0.1530931089239486*fl[4]*nul[35]+0.2651650429449552*fl[26]*nul[33]+0.1530931089239486*fl[18]*nul[33]+0.2651650429449552*fl[19]*nul[32]+0.1530931089239486*fl[11]*nul[32]; 
  incr2[20] = 0.2651650429449552*fl[22]*nul[39]+0.1530931089239486*fl[14]*nul[39]+0.2651650429449552*fl[27]*nul[35]+0.1530931089239486*fl[21]*nul[35]+0.2651650429449552*fl[13]*nul[33]+0.1530931089239486*fl[5]*nul[33]+0.2651650429449552*fl[20]*nul[32]+0.1530931089239486*fl[12]*nul[32]; 
  incr2[22] = 0.2651650429449552*fl[20]*nul[39]+0.1530931089239486*fl[12]*nul[39]+0.2651650429449552*fl[13]*nul[35]+0.1530931089239486*fl[5]*nul[35]+0.2651650429449552*fl[27]*nul[33]+0.1530931089239486*fl[21]*nul[33]+0.2651650429449552*fl[22]*nul[32]+0.1530931089239486*fl[14]*nul[32]; 
  incr2[24] = 0.2651650429449552*fl[31]*nul[39]+0.1530931089239486*fl[29]*nul[39]+0.2651650429449552*fl[30]*nul[35]+0.1530931089239486*fl[25]*nul[35]+0.2651650429449552*fl[28]*nul[33]+0.1530931089239486*fl[23]*nul[33]+0.2651650429449552*fl[24]*nul[32]+0.1530931089239486*fl[15]*nul[32]; 
  incr2[26] = 0.2651650429449552*fl[10]*nul[39]+0.1530931089239486*fl[4]*nul[39]+0.2651650429449552*fl[17]*nul[35]+0.1530931089239486*fl[9]*nul[35]+0.2651650429449552*fl[19]*nul[33]+0.1530931089239486*fl[11]*nul[33]+0.2651650429449552*fl[26]*nul[32]+0.1530931089239486*fl[18]*nul[32]; 
  incr2[27] = 0.2651650429449552*fl[13]*nul[39]+0.1530931089239486*fl[5]*nul[39]+0.2651650429449552*fl[20]*nul[35]+0.1530931089239486*fl[12]*nul[35]+0.2651650429449552*fl[22]*nul[33]+0.1530931089239486*fl[14]*nul[33]+0.2651650429449552*fl[27]*nul[32]+0.1530931089239486*fl[21]*nul[32]; 
  incr2[28] = 0.2651650429449552*fl[30]*nul[39]+0.1530931089239486*fl[25]*nul[39]+0.2651650429449552*fl[31]*nul[35]+0.1530931089239486*fl[29]*nul[35]+0.2651650429449552*fl[24]*nul[33]+0.1530931089239486*fl[15]*nul[33]+0.2651650429449552*fl[28]*nul[32]+0.1530931089239486*fl[23]*nul[32]; 
  incr2[30] = 0.2651650429449552*fl[28]*nul[39]+0.1530931089239486*fl[23]*nul[39]+0.2651650429449552*fl[24]*nul[35]+0.1530931089239486*fl[15]*nul[35]+0.2651650429449552*fl[31]*nul[33]+0.1530931089239486*fl[29]*nul[33]+0.2651650429449552*fl[30]*nul[32]+0.1530931089239486*fl[25]*nul[32]; 
  incr2[31] = 0.2651650429449552*fl[24]*nul[39]+0.1530931089239486*fl[15]*nul[39]+0.2651650429449552*fl[28]*nul[35]+0.1530931089239486*fl[23]*nul[35]+0.2651650429449552*fl[30]*nul[33]+0.1530931089239486*fl[25]*nul[33]+0.2651650429449552*fl[31]*nul[32]+0.1530931089239486*fl[29]*nul[32]; 

  outl[2] += -1.0*incr2[2]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[8] += -1.0*incr2[8]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[13] += -1.0*incr2[13]*rdxFl; 
  outl[16] += -1.0*incr2[16]*rdxFl; 
  outl[17] += -1.0*incr2[17]*rdxFl; 
  outl[19] += -1.0*incr2[19]*rdxFl; 
  outl[20] += -1.0*incr2[20]*rdxFl; 
  outl[22] += -1.0*incr2[22]*rdxFl; 
  outl[24] += -1.0*incr2[24]*rdxFl; 
  outl[26] += -1.0*incr2[26]*rdxFl; 
  outl[27] += -1.0*incr2[27]*rdxFl; 
  outl[28] += -1.0*incr2[28]*rdxFl; 
  outl[30] += -1.0*incr2[30]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf5xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[160]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[2]*dxl[2]); 
  double rdxFr = 4.0/(dxr[2]*dxr[2]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[2] == 1) {

  incr2[3] = (-0.4592793267718456*fr[7]*nul[71])+0.2651650429449552*fr[1]*nul[71]-0.4592793267718456*fr[3]*nul[67]+0.2651650429449552*fr[0]*nul[67]-0.2651650429449552*fr[7]*nul[65]+0.1530931089239486*fr[1]*nul[65]-0.2651650429449552*fr[3]*nul[64]+0.1530931089239486*fr[0]*nul[64]; 
  incr2[7] = (-0.4592793267718456*fr[3]*nul[71])+0.2651650429449552*fr[0]*nul[71]-0.4592793267718456*fr[7]*nul[67]+0.2651650429449552*fr[1]*nul[67]-0.2651650429449552*fr[3]*nul[65]+0.1530931089239486*fr[0]*nul[65]-0.2651650429449552*fr[7]*nul[64]+0.1530931089239486*fr[1]*nul[64]; 
  incr2[8] = (-0.4592793267718456*fr[16]*nul[71])+0.2651650429449552*fr[6]*nul[71]-0.4592793267718456*fr[8]*nul[67]+0.2651650429449552*fr[2]*nul[67]-0.2651650429449552*fr[16]*nul[65]+0.1530931089239486*fr[6]*nul[65]-0.2651650429449552*fr[8]*nul[64]+0.1530931089239486*fr[2]*nul[64]; 
  incr2[11] = (-0.4592793267718456*fr[18]*nul[71])+0.2651650429449552*fr[9]*nul[71]-0.4592793267718456*fr[11]*nul[67]+0.2651650429449552*fr[4]*nul[67]-0.2651650429449552*fr[18]*nul[65]+0.1530931089239486*fr[9]*nul[65]-0.2651650429449552*fr[11]*nul[64]+0.1530931089239486*fr[4]*nul[64]; 
  incr2[14] = (-0.4592793267718456*fr[21]*nul[71])+0.2651650429449552*fr[12]*nul[71]-0.4592793267718456*fr[14]*nul[67]+0.2651650429449552*fr[5]*nul[67]-0.2651650429449552*fr[21]*nul[65]+0.1530931089239486*fr[12]*nul[65]-0.2651650429449552*fr[14]*nul[64]+0.1530931089239486*fr[5]*nul[64]; 
  incr2[16] = (-0.4592793267718456*fr[8]*nul[71])+0.2651650429449552*fr[2]*nul[71]-0.4592793267718456*fr[16]*nul[67]+0.2651650429449552*fr[6]*nul[67]-0.2651650429449552*fr[8]*nul[65]+0.1530931089239486*fr[2]*nul[65]-0.2651650429449552*fr[16]*nul[64]+0.1530931089239486*fr[6]*nul[64]; 
  incr2[18] = (-0.4592793267718456*fr[11]*nul[71])+0.2651650429449552*fr[4]*nul[71]-0.4592793267718456*fr[18]*nul[67]+0.2651650429449552*fr[9]*nul[67]-0.2651650429449552*fr[11]*nul[65]+0.1530931089239486*fr[4]*nul[65]-0.2651650429449552*fr[18]*nul[64]+0.1530931089239486*fr[9]*nul[64]; 
  incr2[19] = (-0.4592793267718456*fr[26]*nul[71])+0.2651650429449552*fr[17]*nul[71]-0.4592793267718456*fr[19]*nul[67]+0.2651650429449552*fr[10]*nul[67]-0.2651650429449552*fr[26]*nul[65]+0.1530931089239486*fr[17]*nul[65]-0.2651650429449552*fr[19]*nul[64]+0.1530931089239486*fr[10]*nul[64]; 
  incr2[21] = (-0.4592793267718456*fr[14]*nul[71])+0.2651650429449552*fr[5]*nul[71]-0.4592793267718456*fr[21]*nul[67]+0.2651650429449552*fr[12]*nul[67]-0.2651650429449552*fr[14]*nul[65]+0.1530931089239486*fr[5]*nul[65]-0.2651650429449552*fr[21]*nul[64]+0.1530931089239486*fr[12]*nul[64]; 
  incr2[22] = (-0.4592793267718456*fr[27]*nul[71])+0.2651650429449552*fr[20]*nul[71]-0.4592793267718456*fr[22]*nul[67]+0.2651650429449552*fr[13]*nul[67]-0.2651650429449552*fr[27]*nul[65]+0.1530931089239486*fr[20]*nul[65]-0.2651650429449552*fr[22]*nul[64]+0.1530931089239486*fr[13]*nul[64]; 
  incr2[25] = (-0.4592793267718456*fr[29]*nul[71])+0.2651650429449552*fr[23]*nul[71]-0.4592793267718456*fr[25]*nul[67]+0.2651650429449552*fr[15]*nul[67]-0.2651650429449552*fr[29]*nul[65]+0.1530931089239486*fr[23]*nul[65]-0.2651650429449552*fr[25]*nul[64]+0.1530931089239486*fr[15]*nul[64]; 
  incr2[26] = (-0.4592793267718456*fr[19]*nul[71])+0.2651650429449552*fr[10]*nul[71]-0.4592793267718456*fr[26]*nul[67]+0.2651650429449552*fr[17]*nul[67]-0.2651650429449552*fr[19]*nul[65]+0.1530931089239486*fr[10]*nul[65]-0.2651650429449552*fr[26]*nul[64]+0.1530931089239486*fr[17]*nul[64]; 
  incr2[27] = (-0.4592793267718456*fr[22]*nul[71])+0.2651650429449552*fr[13]*nul[71]-0.4592793267718456*fr[27]*nul[67]+0.2651650429449552*fr[20]*nul[67]-0.2651650429449552*fr[22]*nul[65]+0.1530931089239486*fr[13]*nul[65]-0.2651650429449552*fr[27]*nul[64]+0.1530931089239486*fr[20]*nul[64]; 
  incr2[29] = (-0.4592793267718456*fr[25]*nul[71])+0.2651650429449552*fr[15]*nul[71]-0.4592793267718456*fr[29]*nul[67]+0.2651650429449552*fr[23]*nul[67]-0.2651650429449552*fr[25]*nul[65]+0.1530931089239486*fr[15]*nul[65]-0.2651650429449552*fr[29]*nul[64]+0.1530931089239486*fr[23]*nul[64]; 
  incr2[30] = (-0.4592793267718456*fr[31]*nul[71])+0.2651650429449552*fr[28]*nul[71]-0.4592793267718456*fr[30]*nul[67]+0.2651650429449552*fr[24]*nul[67]-0.2651650429449552*fr[31]*nul[65]+0.1530931089239486*fr[28]*nul[65]-0.2651650429449552*fr[30]*nul[64]+0.1530931089239486*fr[24]*nul[64]; 
  incr2[31] = (-0.4592793267718456*fr[30]*nul[71])+0.2651650429449552*fr[24]*nul[71]-0.4592793267718456*fr[31]*nul[67]+0.2651650429449552*fr[28]*nul[67]-0.2651650429449552*fr[30]*nul[65]+0.1530931089239486*fr[24]*nul[65]-0.2651650429449552*fr[31]*nul[64]+0.1530931089239486*fr[28]*nul[64]; 

  outr[3] += incr2[3]*rdxFr; 
  outr[7] += incr2[7]*rdxFr; 
  outr[8] += incr2[8]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 
  outr[14] += incr2[14]*rdxFr; 
  outr[16] += incr2[16]*rdxFr; 
  outr[18] += incr2[18]*rdxFr; 
  outr[19] += incr2[19]*rdxFr; 
  outr[21] += incr2[21]*rdxFr; 
  outr[22] += incr2[22]*rdxFr; 
  outr[25] += incr2[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 
  outr[27] += incr2[27]*rdxFr; 
  outr[29] += incr2[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 

  } else {

  incr2[3] = 0.4592793267718456*fl[7]*nul[71]+0.2651650429449552*fl[1]*nul[71]+0.4592793267718456*fl[3]*nul[67]+0.2651650429449552*fl[0]*nul[67]+0.2651650429449552*fl[7]*nul[65]+0.1530931089239486*fl[1]*nul[65]+0.2651650429449552*fl[3]*nul[64]+0.1530931089239486*fl[0]*nul[64]; 
  incr2[7] = 0.4592793267718456*fl[3]*nul[71]+0.2651650429449552*fl[0]*nul[71]+0.4592793267718456*fl[7]*nul[67]+0.2651650429449552*fl[1]*nul[67]+0.2651650429449552*fl[3]*nul[65]+0.1530931089239486*fl[0]*nul[65]+0.2651650429449552*fl[7]*nul[64]+0.1530931089239486*fl[1]*nul[64]; 
  incr2[8] = 0.4592793267718456*fl[16]*nul[71]+0.2651650429449552*fl[6]*nul[71]+0.4592793267718456*fl[8]*nul[67]+0.2651650429449552*fl[2]*nul[67]+0.2651650429449552*fl[16]*nul[65]+0.1530931089239486*fl[6]*nul[65]+0.2651650429449552*fl[8]*nul[64]+0.1530931089239486*fl[2]*nul[64]; 
  incr2[11] = 0.4592793267718456*fl[18]*nul[71]+0.2651650429449552*fl[9]*nul[71]+0.4592793267718456*fl[11]*nul[67]+0.2651650429449552*fl[4]*nul[67]+0.2651650429449552*fl[18]*nul[65]+0.1530931089239486*fl[9]*nul[65]+0.2651650429449552*fl[11]*nul[64]+0.1530931089239486*fl[4]*nul[64]; 
  incr2[14] = 0.4592793267718456*fl[21]*nul[71]+0.2651650429449552*fl[12]*nul[71]+0.4592793267718456*fl[14]*nul[67]+0.2651650429449552*fl[5]*nul[67]+0.2651650429449552*fl[21]*nul[65]+0.1530931089239486*fl[12]*nul[65]+0.2651650429449552*fl[14]*nul[64]+0.1530931089239486*fl[5]*nul[64]; 
  incr2[16] = 0.4592793267718456*fl[8]*nul[71]+0.2651650429449552*fl[2]*nul[71]+0.4592793267718456*fl[16]*nul[67]+0.2651650429449552*fl[6]*nul[67]+0.2651650429449552*fl[8]*nul[65]+0.1530931089239486*fl[2]*nul[65]+0.2651650429449552*fl[16]*nul[64]+0.1530931089239486*fl[6]*nul[64]; 
  incr2[18] = 0.4592793267718456*fl[11]*nul[71]+0.2651650429449552*fl[4]*nul[71]+0.4592793267718456*fl[18]*nul[67]+0.2651650429449552*fl[9]*nul[67]+0.2651650429449552*fl[11]*nul[65]+0.1530931089239486*fl[4]*nul[65]+0.2651650429449552*fl[18]*nul[64]+0.1530931089239486*fl[9]*nul[64]; 
  incr2[19] = 0.4592793267718456*fl[26]*nul[71]+0.2651650429449552*fl[17]*nul[71]+0.4592793267718456*fl[19]*nul[67]+0.2651650429449552*fl[10]*nul[67]+0.2651650429449552*fl[26]*nul[65]+0.1530931089239486*fl[17]*nul[65]+0.2651650429449552*fl[19]*nul[64]+0.1530931089239486*fl[10]*nul[64]; 
  incr2[21] = 0.4592793267718456*fl[14]*nul[71]+0.2651650429449552*fl[5]*nul[71]+0.4592793267718456*fl[21]*nul[67]+0.2651650429449552*fl[12]*nul[67]+0.2651650429449552*fl[14]*nul[65]+0.1530931089239486*fl[5]*nul[65]+0.2651650429449552*fl[21]*nul[64]+0.1530931089239486*fl[12]*nul[64]; 
  incr2[22] = 0.4592793267718456*fl[27]*nul[71]+0.2651650429449552*fl[20]*nul[71]+0.4592793267718456*fl[22]*nul[67]+0.2651650429449552*fl[13]*nul[67]+0.2651650429449552*fl[27]*nul[65]+0.1530931089239486*fl[20]*nul[65]+0.2651650429449552*fl[22]*nul[64]+0.1530931089239486*fl[13]*nul[64]; 
  incr2[25] = 0.4592793267718456*fl[29]*nul[71]+0.2651650429449552*fl[23]*nul[71]+0.4592793267718456*fl[25]*nul[67]+0.2651650429449552*fl[15]*nul[67]+0.2651650429449552*fl[29]*nul[65]+0.1530931089239486*fl[23]*nul[65]+0.2651650429449552*fl[25]*nul[64]+0.1530931089239486*fl[15]*nul[64]; 
  incr2[26] = 0.4592793267718456*fl[19]*nul[71]+0.2651650429449552*fl[10]*nul[71]+0.4592793267718456*fl[26]*nul[67]+0.2651650429449552*fl[17]*nul[67]+0.2651650429449552*fl[19]*nul[65]+0.1530931089239486*fl[10]*nul[65]+0.2651650429449552*fl[26]*nul[64]+0.1530931089239486*fl[17]*nul[64]; 
  incr2[27] = 0.4592793267718456*fl[22]*nul[71]+0.2651650429449552*fl[13]*nul[71]+0.4592793267718456*fl[27]*nul[67]+0.2651650429449552*fl[20]*nul[67]+0.2651650429449552*fl[22]*nul[65]+0.1530931089239486*fl[13]*nul[65]+0.2651650429449552*fl[27]*nul[64]+0.1530931089239486*fl[20]*nul[64]; 
  incr2[29] = 0.4592793267718456*fl[25]*nul[71]+0.2651650429449552*fl[15]*nul[71]+0.4592793267718456*fl[29]*nul[67]+0.2651650429449552*fl[23]*nul[67]+0.2651650429449552*fl[25]*nul[65]+0.1530931089239486*fl[15]*nul[65]+0.2651650429449552*fl[29]*nul[64]+0.1530931089239486*fl[23]*nul[64]; 
  incr2[30] = 0.4592793267718456*fl[31]*nul[71]+0.2651650429449552*fl[28]*nul[71]+0.4592793267718456*fl[30]*nul[67]+0.2651650429449552*fl[24]*nul[67]+0.2651650429449552*fl[31]*nul[65]+0.1530931089239486*fl[28]*nul[65]+0.2651650429449552*fl[30]*nul[64]+0.1530931089239486*fl[24]*nul[64]; 
  incr2[31] = 0.4592793267718456*fl[30]*nul[71]+0.2651650429449552*fl[24]*nul[71]+0.4592793267718456*fl[31]*nul[67]+0.2651650429449552*fl[28]*nul[67]+0.2651650429449552*fl[30]*nul[65]+0.1530931089239486*fl[24]*nul[65]+0.2651650429449552*fl[31]*nul[64]+0.1530931089239486*fl[28]*nul[64]; 

  outl[3] += -1.0*incr2[3]*rdxFl; 
  outl[7] += -1.0*incr2[7]*rdxFl; 
  outl[8] += -1.0*incr2[8]*rdxFl; 
  outl[11] += -1.0*incr2[11]*rdxFl; 
  outl[14] += -1.0*incr2[14]*rdxFl; 
  outl[16] += -1.0*incr2[16]*rdxFl; 
  outl[18] += -1.0*incr2[18]*rdxFl; 
  outl[19] += -1.0*incr2[19]*rdxFl; 
  outl[21] += -1.0*incr2[21]*rdxFl; 
  outl[22] += -1.0*incr2[22]*rdxFl; 
  outl[25] += -1.0*incr2[25]*rdxFl; 
  outl[26] += -1.0*incr2[26]*rdxFl; 
  outl[27] += -1.0*incr2[27]*rdxFl; 
  outl[29] += -1.0*incr2[29]*rdxFl; 
  outl[30] += -1.0*incr2[30]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf5xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[160]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[3]*dxl[3]); 
  double rdxFr = 4.0/(dxr[3]*dxr[3]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[3] == 1) {

  incr2[4] = (-0.2651650429449552*fr[18]*nul[103])+0.1530931089239486*fr[7]*nul[103]-0.2651650429449552*fr[11]*nul[99]+0.1530931089239486*fr[3]*nul[99]-0.2651650429449552*fr[9]*nul[97]+0.1530931089239486*fr[1]*nul[97]-0.2651650429449552*fr[4]*nul[96]+0.1530931089239486*fr[0]*nul[96]; 
  incr2[9] = (-0.2651650429449552*fr[11]*nul[103])+0.1530931089239486*fr[3]*nul[103]-0.2651650429449552*fr[18]*nul[99]+0.1530931089239486*fr[7]*nul[99]-0.2651650429449552*fr[4]*nul[97]+0.1530931089239486*fr[0]*nul[97]-0.2651650429449552*fr[9]*nul[96]+0.1530931089239486*fr[1]*nul[96]; 
  incr2[10] = (-0.2651650429449552*fr[26]*nul[103])+0.1530931089239486*fr[16]*nul[103]-0.2651650429449552*fr[19]*nul[99]+0.1530931089239486*fr[8]*nul[99]-0.2651650429449552*fr[17]*nul[97]+0.1530931089239486*fr[6]*nul[97]-0.2651650429449552*fr[10]*nul[96]+0.1530931089239486*fr[2]*nul[96]; 
  incr2[11] = (-0.2651650429449552*fr[9]*nul[103])+0.1530931089239486*fr[1]*nul[103]-0.2651650429449552*fr[4]*nul[99]+0.1530931089239486*fr[0]*nul[99]-0.2651650429449552*fr[18]*nul[97]+0.1530931089239486*fr[7]*nul[97]-0.2651650429449552*fr[11]*nul[96]+0.1530931089239486*fr[3]*nul[96]; 
  incr2[15] = (-0.2651650429449552*fr[29]*nul[103])+0.1530931089239486*fr[21]*nul[103]-0.2651650429449552*fr[25]*nul[99]+0.1530931089239486*fr[14]*nul[99]-0.2651650429449552*fr[23]*nul[97]+0.1530931089239486*fr[12]*nul[97]-0.2651650429449552*fr[15]*nul[96]+0.1530931089239486*fr[5]*nul[96]; 
  incr2[17] = (-0.2651650429449552*fr[19]*nul[103])+0.1530931089239486*fr[8]*nul[103]-0.2651650429449552*fr[26]*nul[99]+0.1530931089239486*fr[16]*nul[99]-0.2651650429449552*fr[10]*nul[97]+0.1530931089239486*fr[2]*nul[97]-0.2651650429449552*fr[17]*nul[96]+0.1530931089239486*fr[6]*nul[96]; 
  incr2[18] = (-0.2651650429449552*fr[4]*nul[103])+0.1530931089239486*fr[0]*nul[103]-0.2651650429449552*fr[9]*nul[99]+0.1530931089239486*fr[1]*nul[99]-0.2651650429449552*fr[11]*nul[97]+0.1530931089239486*fr[3]*nul[97]-0.2651650429449552*fr[18]*nul[96]+0.1530931089239486*fr[7]*nul[96]; 
  incr2[19] = (-0.2651650429449552*fr[17]*nul[103])+0.1530931089239486*fr[6]*nul[103]-0.2651650429449552*fr[10]*nul[99]+0.1530931089239486*fr[2]*nul[99]-0.2651650429449552*fr[26]*nul[97]+0.1530931089239486*fr[16]*nul[97]-0.2651650429449552*fr[19]*nul[96]+0.1530931089239486*fr[8]*nul[96]; 
  incr2[23] = (-0.2651650429449552*fr[25]*nul[103])+0.1530931089239486*fr[14]*nul[103]-0.2651650429449552*fr[29]*nul[99]+0.1530931089239486*fr[21]*nul[99]-0.2651650429449552*fr[15]*nul[97]+0.1530931089239486*fr[5]*nul[97]-0.2651650429449552*fr[23]*nul[96]+0.1530931089239486*fr[12]*nul[96]; 
  incr2[24] = (-0.2651650429449552*fr[31]*nul[103])+0.1530931089239486*fr[27]*nul[103]-0.2651650429449552*fr[30]*nul[99]+0.1530931089239486*fr[22]*nul[99]-0.2651650429449552*fr[28]*nul[97]+0.1530931089239486*fr[20]*nul[97]-0.2651650429449552*fr[24]*nul[96]+0.1530931089239486*fr[13]*nul[96]; 
  incr2[25] = (-0.2651650429449552*fr[23]*nul[103])+0.1530931089239486*fr[12]*nul[103]-0.2651650429449552*fr[15]*nul[99]+0.1530931089239486*fr[5]*nul[99]-0.2651650429449552*fr[29]*nul[97]+0.1530931089239486*fr[21]*nul[97]-0.2651650429449552*fr[25]*nul[96]+0.1530931089239486*fr[14]*nul[96]; 
  incr2[26] = (-0.2651650429449552*fr[10]*nul[103])+0.1530931089239486*fr[2]*nul[103]-0.2651650429449552*fr[17]*nul[99]+0.1530931089239486*fr[6]*nul[99]-0.2651650429449552*fr[19]*nul[97]+0.1530931089239486*fr[8]*nul[97]-0.2651650429449552*fr[26]*nul[96]+0.1530931089239486*fr[16]*nul[96]; 
  incr2[28] = (-0.2651650429449552*fr[30]*nul[103])+0.1530931089239486*fr[22]*nul[103]-0.2651650429449552*fr[31]*nul[99]+0.1530931089239486*fr[27]*nul[99]-0.2651650429449552*fr[24]*nul[97]+0.1530931089239486*fr[13]*nul[97]-0.2651650429449552*fr[28]*nul[96]+0.1530931089239486*fr[20]*nul[96]; 
  incr2[29] = (-0.2651650429449552*fr[15]*nul[103])+0.1530931089239486*fr[5]*nul[103]-0.2651650429449552*fr[23]*nul[99]+0.1530931089239486*fr[12]*nul[99]-0.2651650429449552*fr[25]*nul[97]+0.1530931089239486*fr[14]*nul[97]-0.2651650429449552*fr[29]*nul[96]+0.1530931089239486*fr[21]*nul[96]; 
  incr2[30] = (-0.2651650429449552*fr[28]*nul[103])+0.1530931089239486*fr[20]*nul[103]-0.2651650429449552*fr[24]*nul[99]+0.1530931089239486*fr[13]*nul[99]-0.2651650429449552*fr[31]*nul[97]+0.1530931089239486*fr[27]*nul[97]-0.2651650429449552*fr[30]*nul[96]+0.1530931089239486*fr[22]*nul[96]; 
  incr2[31] = (-0.2651650429449552*fr[24]*nul[103])+0.1530931089239486*fr[13]*nul[103]-0.2651650429449552*fr[28]*nul[99]+0.1530931089239486*fr[20]*nul[99]-0.2651650429449552*fr[30]*nul[97]+0.1530931089239486*fr[22]*nul[97]-0.2651650429449552*fr[31]*nul[96]+0.1530931089239486*fr[27]*nul[96]; 

  outr[4] += incr2[4]*rdxFr; 
  outr[9] += incr2[9]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 
  outr[17] += incr2[17]*rdxFr; 
  outr[18] += incr2[18]*rdxFr; 
  outr[19] += incr2[19]*rdxFr; 
  outr[23] += incr2[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr; 
  outr[25] += incr2[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 
  outr[28] += incr2[28]*rdxFr; 
  outr[29] += incr2[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 

  } else {

  incr2[4] = 0.2651650429449552*fl[18]*nul[103]+0.1530931089239486*fl[7]*nul[103]+0.2651650429449552*fl[11]*nul[99]+0.1530931089239486*fl[3]*nul[99]+0.2651650429449552*fl[9]*nul[97]+0.1530931089239486*fl[1]*nul[97]+0.2651650429449552*fl[4]*nul[96]+0.1530931089239486*fl[0]*nul[96]; 
  incr2[9] = 0.2651650429449552*fl[11]*nul[103]+0.1530931089239486*fl[3]*nul[103]+0.2651650429449552*fl[18]*nul[99]+0.1530931089239486*fl[7]*nul[99]+0.2651650429449552*fl[4]*nul[97]+0.1530931089239486*fl[0]*nul[97]+0.2651650429449552*fl[9]*nul[96]+0.1530931089239486*fl[1]*nul[96]; 
  incr2[10] = 0.2651650429449552*fl[26]*nul[103]+0.1530931089239486*fl[16]*nul[103]+0.2651650429449552*fl[19]*nul[99]+0.1530931089239486*fl[8]*nul[99]+0.2651650429449552*fl[17]*nul[97]+0.1530931089239486*fl[6]*nul[97]+0.2651650429449552*fl[10]*nul[96]+0.1530931089239486*fl[2]*nul[96]; 
  incr2[11] = 0.2651650429449552*fl[9]*nul[103]+0.1530931089239486*fl[1]*nul[103]+0.2651650429449552*fl[4]*nul[99]+0.1530931089239486*fl[0]*nul[99]+0.2651650429449552*fl[18]*nul[97]+0.1530931089239486*fl[7]*nul[97]+0.2651650429449552*fl[11]*nul[96]+0.1530931089239486*fl[3]*nul[96]; 
  incr2[15] = 0.2651650429449552*fl[29]*nul[103]+0.1530931089239486*fl[21]*nul[103]+0.2651650429449552*fl[25]*nul[99]+0.1530931089239486*fl[14]*nul[99]+0.2651650429449552*fl[23]*nul[97]+0.1530931089239486*fl[12]*nul[97]+0.2651650429449552*fl[15]*nul[96]+0.1530931089239486*fl[5]*nul[96]; 
  incr2[17] = 0.2651650429449552*fl[19]*nul[103]+0.1530931089239486*fl[8]*nul[103]+0.2651650429449552*fl[26]*nul[99]+0.1530931089239486*fl[16]*nul[99]+0.2651650429449552*fl[10]*nul[97]+0.1530931089239486*fl[2]*nul[97]+0.2651650429449552*fl[17]*nul[96]+0.1530931089239486*fl[6]*nul[96]; 
  incr2[18] = 0.2651650429449552*fl[4]*nul[103]+0.1530931089239486*fl[0]*nul[103]+0.2651650429449552*fl[9]*nul[99]+0.1530931089239486*fl[1]*nul[99]+0.2651650429449552*fl[11]*nul[97]+0.1530931089239486*fl[3]*nul[97]+0.2651650429449552*fl[18]*nul[96]+0.1530931089239486*fl[7]*nul[96]; 
  incr2[19] = 0.2651650429449552*fl[17]*nul[103]+0.1530931089239486*fl[6]*nul[103]+0.2651650429449552*fl[10]*nul[99]+0.1530931089239486*fl[2]*nul[99]+0.2651650429449552*fl[26]*nul[97]+0.1530931089239486*fl[16]*nul[97]+0.2651650429449552*fl[19]*nul[96]+0.1530931089239486*fl[8]*nul[96]; 
  incr2[23] = 0.2651650429449552*fl[25]*nul[103]+0.1530931089239486*fl[14]*nul[103]+0.2651650429449552*fl[29]*nul[99]+0.1530931089239486*fl[21]*nul[99]+0.2651650429449552*fl[15]*nul[97]+0.1530931089239486*fl[5]*nul[97]+0.2651650429449552*fl[23]*nul[96]+0.1530931089239486*fl[12]*nul[96]; 
  incr2[24] = 0.2651650429449552*fl[31]*nul[103]+0.1530931089239486*fl[27]*nul[103]+0.2651650429449552*fl[30]*nul[99]+0.1530931089239486*fl[22]*nul[99]+0.2651650429449552*fl[28]*nul[97]+0.1530931089239486*fl[20]*nul[97]+0.2651650429449552*fl[24]*nul[96]+0.1530931089239486*fl[13]*nul[96]; 
  incr2[25] = 0.2651650429449552*fl[23]*nul[103]+0.1530931089239486*fl[12]*nul[103]+0.2651650429449552*fl[15]*nul[99]+0.1530931089239486*fl[5]*nul[99]+0.2651650429449552*fl[29]*nul[97]+0.1530931089239486*fl[21]*nul[97]+0.2651650429449552*fl[25]*nul[96]+0.1530931089239486*fl[14]*nul[96]; 
  incr2[26] = 0.2651650429449552*fl[10]*nul[103]+0.1530931089239486*fl[2]*nul[103]+0.2651650429449552*fl[17]*nul[99]+0.1530931089239486*fl[6]*nul[99]+0.2651650429449552*fl[19]*nul[97]+0.1530931089239486*fl[8]*nul[97]+0.2651650429449552*fl[26]*nul[96]+0.1530931089239486*fl[16]*nul[96]; 
  incr2[28] = 0.2651650429449552*fl[30]*nul[103]+0.1530931089239486*fl[22]*nul[103]+0.2651650429449552*fl[31]*nul[99]+0.1530931089239486*fl[27]*nul[99]+0.2651650429449552*fl[24]*nul[97]+0.1530931089239486*fl[13]*nul[97]+0.2651650429449552*fl[28]*nul[96]+0.1530931089239486*fl[20]*nul[96]; 
  incr2[29] = 0.2651650429449552*fl[15]*nul[103]+0.1530931089239486*fl[5]*nul[103]+0.2651650429449552*fl[23]*nul[99]+0.1530931089239486*fl[12]*nul[99]+0.2651650429449552*fl[25]*nul[97]+0.1530931089239486*fl[14]*nul[97]+0.2651650429449552*fl[29]*nul[96]+0.1530931089239486*fl[21]*nul[96]; 
  incr2[30] = 0.2651650429449552*fl[28]*nul[103]+0.1530931089239486*fl[20]*nul[103]+0.2651650429449552*fl[24]*nul[99]+0.1530931089239486*fl[13]*nul[99]+0.2651650429449552*fl[31]*nul[97]+0.1530931089239486*fl[27]*nul[97]+0.2651650429449552*fl[30]*nul[96]+0.1530931089239486*fl[22]*nul[96]; 
  incr2[31] = 0.2651650429449552*fl[24]*nul[103]+0.1530931089239486*fl[13]*nul[103]+0.2651650429449552*fl[28]*nul[99]+0.1530931089239486*fl[20]*nul[99]+0.2651650429449552*fl[30]*nul[97]+0.1530931089239486*fl[22]*nul[97]+0.2651650429449552*fl[31]*nul[96]+0.1530931089239486*fl[27]*nul[96]; 

  outl[4] += -1.0*incr2[4]*rdxFl; 
  outl[9] += -1.0*incr2[9]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[11] += -1.0*incr2[11]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 
  outl[17] += -1.0*incr2[17]*rdxFl; 
  outl[18] += -1.0*incr2[18]*rdxFl; 
  outl[19] += -1.0*incr2[19]*rdxFl; 
  outl[23] += -1.0*incr2[23]*rdxFl; 
  outl[24] += -1.0*incr2[24]*rdxFl; 
  outl[25] += -1.0*incr2[25]*rdxFl; 
  outl[26] += -1.0*incr2[26]*rdxFl; 
  outl[28] += -1.0*incr2[28]*rdxFl; 
  outl[29] += -1.0*incr2[29]*rdxFl; 
  outl[30] += -1.0*incr2[30]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf5xSerP1_X5(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // idx[5]:    current grid index.
  // nu[160]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[4]*dxl[4]); 
  double rdxFr = 4.0/(dxr[4]*dxr[4]); 

  double incr1[32]; 
  double incr2[32]; 

  if (idxr[4] == 1) {

  incr2[5] = (-0.2651650429449552*fr[21]*nul[135])+0.1530931089239486*fr[7]*nul[135]-0.2651650429449552*fr[14]*nul[131]+0.1530931089239486*fr[3]*nul[131]-0.2651650429449552*fr[12]*nul[129]+0.1530931089239486*fr[1]*nul[129]-0.2651650429449552*fr[5]*nul[128]+0.1530931089239486*fr[0]*nul[128]; 
  incr2[12] = (-0.2651650429449552*fr[14]*nul[135])+0.1530931089239486*fr[3]*nul[135]-0.2651650429449552*fr[21]*nul[131]+0.1530931089239486*fr[7]*nul[131]-0.2651650429449552*fr[5]*nul[129]+0.1530931089239486*fr[0]*nul[129]-0.2651650429449552*fr[12]*nul[128]+0.1530931089239486*fr[1]*nul[128]; 
  incr2[13] = (-0.2651650429449552*fr[27]*nul[135])+0.1530931089239486*fr[16]*nul[135]-0.2651650429449552*fr[22]*nul[131]+0.1530931089239486*fr[8]*nul[131]-0.2651650429449552*fr[20]*nul[129]+0.1530931089239486*fr[6]*nul[129]-0.2651650429449552*fr[13]*nul[128]+0.1530931089239486*fr[2]*nul[128]; 
  incr2[14] = (-0.2651650429449552*fr[12]*nul[135])+0.1530931089239486*fr[1]*nul[135]-0.2651650429449552*fr[5]*nul[131]+0.1530931089239486*fr[0]*nul[131]-0.2651650429449552*fr[21]*nul[129]+0.1530931089239486*fr[7]*nul[129]-0.2651650429449552*fr[14]*nul[128]+0.1530931089239486*fr[3]*nul[128]; 
  incr2[15] = (-0.2651650429449552*fr[29]*nul[135])+0.1530931089239486*fr[18]*nul[135]-0.2651650429449552*fr[25]*nul[131]+0.1530931089239486*fr[11]*nul[131]-0.2651650429449552*fr[23]*nul[129]+0.1530931089239486*fr[9]*nul[129]-0.2651650429449552*fr[15]*nul[128]+0.1530931089239486*fr[4]*nul[128]; 
  incr2[20] = (-0.2651650429449552*fr[22]*nul[135])+0.1530931089239486*fr[8]*nul[135]-0.2651650429449552*fr[27]*nul[131]+0.1530931089239486*fr[16]*nul[131]-0.2651650429449552*fr[13]*nul[129]+0.1530931089239486*fr[2]*nul[129]-0.2651650429449552*fr[20]*nul[128]+0.1530931089239486*fr[6]*nul[128]; 
  incr2[21] = (-0.2651650429449552*fr[5]*nul[135])+0.1530931089239486*fr[0]*nul[135]-0.2651650429449552*fr[12]*nul[131]+0.1530931089239486*fr[1]*nul[131]-0.2651650429449552*fr[14]*nul[129]+0.1530931089239486*fr[3]*nul[129]-0.2651650429449552*fr[21]*nul[128]+0.1530931089239486*fr[7]*nul[128]; 
  incr2[22] = (-0.2651650429449552*fr[20]*nul[135])+0.1530931089239486*fr[6]*nul[135]-0.2651650429449552*fr[13]*nul[131]+0.1530931089239486*fr[2]*nul[131]-0.2651650429449552*fr[27]*nul[129]+0.1530931089239486*fr[16]*nul[129]-0.2651650429449552*fr[22]*nul[128]+0.1530931089239486*fr[8]*nul[128]; 
  incr2[23] = (-0.2651650429449552*fr[25]*nul[135])+0.1530931089239486*fr[11]*nul[135]-0.2651650429449552*fr[29]*nul[131]+0.1530931089239486*fr[18]*nul[131]-0.2651650429449552*fr[15]*nul[129]+0.1530931089239486*fr[4]*nul[129]-0.2651650429449552*fr[23]*nul[128]+0.1530931089239486*fr[9]*nul[128]; 
  incr2[24] = (-0.2651650429449552*fr[31]*nul[135])+0.1530931089239486*fr[26]*nul[135]-0.2651650429449552*fr[30]*nul[131]+0.1530931089239486*fr[19]*nul[131]-0.2651650429449552*fr[28]*nul[129]+0.1530931089239486*fr[17]*nul[129]-0.2651650429449552*fr[24]*nul[128]+0.1530931089239486*fr[10]*nul[128]; 
  incr2[25] = (-0.2651650429449552*fr[23]*nul[135])+0.1530931089239486*fr[9]*nul[135]-0.2651650429449552*fr[15]*nul[131]+0.1530931089239486*fr[4]*nul[131]-0.2651650429449552*fr[29]*nul[129]+0.1530931089239486*fr[18]*nul[129]-0.2651650429449552*fr[25]*nul[128]+0.1530931089239486*fr[11]*nul[128]; 
  incr2[27] = (-0.2651650429449552*fr[13]*nul[135])+0.1530931089239486*fr[2]*nul[135]-0.2651650429449552*fr[20]*nul[131]+0.1530931089239486*fr[6]*nul[131]-0.2651650429449552*fr[22]*nul[129]+0.1530931089239486*fr[8]*nul[129]-0.2651650429449552*fr[27]*nul[128]+0.1530931089239486*fr[16]*nul[128]; 
  incr2[28] = (-0.2651650429449552*fr[30]*nul[135])+0.1530931089239486*fr[19]*nul[135]-0.2651650429449552*fr[31]*nul[131]+0.1530931089239486*fr[26]*nul[131]-0.2651650429449552*fr[24]*nul[129]+0.1530931089239486*fr[10]*nul[129]-0.2651650429449552*fr[28]*nul[128]+0.1530931089239486*fr[17]*nul[128]; 
  incr2[29] = (-0.2651650429449552*fr[15]*nul[135])+0.1530931089239486*fr[4]*nul[135]-0.2651650429449552*fr[23]*nul[131]+0.1530931089239486*fr[9]*nul[131]-0.2651650429449552*fr[25]*nul[129]+0.1530931089239486*fr[11]*nul[129]-0.2651650429449552*fr[29]*nul[128]+0.1530931089239486*fr[18]*nul[128]; 
  incr2[30] = (-0.2651650429449552*fr[28]*nul[135])+0.1530931089239486*fr[17]*nul[135]-0.2651650429449552*fr[24]*nul[131]+0.1530931089239486*fr[10]*nul[131]-0.2651650429449552*fr[31]*nul[129]+0.1530931089239486*fr[26]*nul[129]-0.2651650429449552*fr[30]*nul[128]+0.1530931089239486*fr[19]*nul[128]; 
  incr2[31] = (-0.2651650429449552*fr[24]*nul[135])+0.1530931089239486*fr[10]*nul[135]-0.2651650429449552*fr[28]*nul[131]+0.1530931089239486*fr[17]*nul[131]-0.2651650429449552*fr[30]*nul[129]+0.1530931089239486*fr[19]*nul[129]-0.2651650429449552*fr[31]*nul[128]+0.1530931089239486*fr[26]*nul[128]; 

  outr[5] += incr2[5]*rdxFr; 
  outr[12] += incr2[12]*rdxFr; 
  outr[13] += incr2[13]*rdxFr; 
  outr[14] += incr2[14]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 
  outr[20] += incr2[20]*rdxFr; 
  outr[21] += incr2[21]*rdxFr; 
  outr[22] += incr2[22]*rdxFr; 
  outr[23] += incr2[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr; 
  outr[25] += incr2[25]*rdxFr; 
  outr[27] += incr2[27]*rdxFr; 
  outr[28] += incr2[28]*rdxFr; 
  outr[29] += incr2[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr; 

  } else {

  incr2[5] = 0.2651650429449552*fl[21]*nul[135]+0.1530931089239486*fl[7]*nul[135]+0.2651650429449552*fl[14]*nul[131]+0.1530931089239486*fl[3]*nul[131]+0.2651650429449552*fl[12]*nul[129]+0.1530931089239486*fl[1]*nul[129]+0.2651650429449552*fl[5]*nul[128]+0.1530931089239486*fl[0]*nul[128]; 
  incr2[12] = 0.2651650429449552*fl[14]*nul[135]+0.1530931089239486*fl[3]*nul[135]+0.2651650429449552*fl[21]*nul[131]+0.1530931089239486*fl[7]*nul[131]+0.2651650429449552*fl[5]*nul[129]+0.1530931089239486*fl[0]*nul[129]+0.2651650429449552*fl[12]*nul[128]+0.1530931089239486*fl[1]*nul[128]; 
  incr2[13] = 0.2651650429449552*fl[27]*nul[135]+0.1530931089239486*fl[16]*nul[135]+0.2651650429449552*fl[22]*nul[131]+0.1530931089239486*fl[8]*nul[131]+0.2651650429449552*fl[20]*nul[129]+0.1530931089239486*fl[6]*nul[129]+0.2651650429449552*fl[13]*nul[128]+0.1530931089239486*fl[2]*nul[128]; 
  incr2[14] = 0.2651650429449552*fl[12]*nul[135]+0.1530931089239486*fl[1]*nul[135]+0.2651650429449552*fl[5]*nul[131]+0.1530931089239486*fl[0]*nul[131]+0.2651650429449552*fl[21]*nul[129]+0.1530931089239486*fl[7]*nul[129]+0.2651650429449552*fl[14]*nul[128]+0.1530931089239486*fl[3]*nul[128]; 
  incr2[15] = 0.2651650429449552*fl[29]*nul[135]+0.1530931089239486*fl[18]*nul[135]+0.2651650429449552*fl[25]*nul[131]+0.1530931089239486*fl[11]*nul[131]+0.2651650429449552*fl[23]*nul[129]+0.1530931089239486*fl[9]*nul[129]+0.2651650429449552*fl[15]*nul[128]+0.1530931089239486*fl[4]*nul[128]; 
  incr2[20] = 0.2651650429449552*fl[22]*nul[135]+0.1530931089239486*fl[8]*nul[135]+0.2651650429449552*fl[27]*nul[131]+0.1530931089239486*fl[16]*nul[131]+0.2651650429449552*fl[13]*nul[129]+0.1530931089239486*fl[2]*nul[129]+0.2651650429449552*fl[20]*nul[128]+0.1530931089239486*fl[6]*nul[128]; 
  incr2[21] = 0.2651650429449552*fl[5]*nul[135]+0.1530931089239486*fl[0]*nul[135]+0.2651650429449552*fl[12]*nul[131]+0.1530931089239486*fl[1]*nul[131]+0.2651650429449552*fl[14]*nul[129]+0.1530931089239486*fl[3]*nul[129]+0.2651650429449552*fl[21]*nul[128]+0.1530931089239486*fl[7]*nul[128]; 
  incr2[22] = 0.2651650429449552*fl[20]*nul[135]+0.1530931089239486*fl[6]*nul[135]+0.2651650429449552*fl[13]*nul[131]+0.1530931089239486*fl[2]*nul[131]+0.2651650429449552*fl[27]*nul[129]+0.1530931089239486*fl[16]*nul[129]+0.2651650429449552*fl[22]*nul[128]+0.1530931089239486*fl[8]*nul[128]; 
  incr2[23] = 0.2651650429449552*fl[25]*nul[135]+0.1530931089239486*fl[11]*nul[135]+0.2651650429449552*fl[29]*nul[131]+0.1530931089239486*fl[18]*nul[131]+0.2651650429449552*fl[15]*nul[129]+0.1530931089239486*fl[4]*nul[129]+0.2651650429449552*fl[23]*nul[128]+0.1530931089239486*fl[9]*nul[128]; 
  incr2[24] = 0.2651650429449552*fl[31]*nul[135]+0.1530931089239486*fl[26]*nul[135]+0.2651650429449552*fl[30]*nul[131]+0.1530931089239486*fl[19]*nul[131]+0.2651650429449552*fl[28]*nul[129]+0.1530931089239486*fl[17]*nul[129]+0.2651650429449552*fl[24]*nul[128]+0.1530931089239486*fl[10]*nul[128]; 
  incr2[25] = 0.2651650429449552*fl[23]*nul[135]+0.1530931089239486*fl[9]*nul[135]+0.2651650429449552*fl[15]*nul[131]+0.1530931089239486*fl[4]*nul[131]+0.2651650429449552*fl[29]*nul[129]+0.1530931089239486*fl[18]*nul[129]+0.2651650429449552*fl[25]*nul[128]+0.1530931089239486*fl[11]*nul[128]; 
  incr2[27] = 0.2651650429449552*fl[13]*nul[135]+0.1530931089239486*fl[2]*nul[135]+0.2651650429449552*fl[20]*nul[131]+0.1530931089239486*fl[6]*nul[131]+0.2651650429449552*fl[22]*nul[129]+0.1530931089239486*fl[8]*nul[129]+0.2651650429449552*fl[27]*nul[128]+0.1530931089239486*fl[16]*nul[128]; 
  incr2[28] = 0.2651650429449552*fl[30]*nul[135]+0.1530931089239486*fl[19]*nul[135]+0.2651650429449552*fl[31]*nul[131]+0.1530931089239486*fl[26]*nul[131]+0.2651650429449552*fl[24]*nul[129]+0.1530931089239486*fl[10]*nul[129]+0.2651650429449552*fl[28]*nul[128]+0.1530931089239486*fl[17]*nul[128]; 
  incr2[29] = 0.2651650429449552*fl[15]*nul[135]+0.1530931089239486*fl[4]*nul[135]+0.2651650429449552*fl[23]*nul[131]+0.1530931089239486*fl[9]*nul[131]+0.2651650429449552*fl[25]*nul[129]+0.1530931089239486*fl[11]*nul[129]+0.2651650429449552*fl[29]*nul[128]+0.1530931089239486*fl[18]*nul[128]; 
  incr2[30] = 0.2651650429449552*fl[28]*nul[135]+0.1530931089239486*fl[17]*nul[135]+0.2651650429449552*fl[24]*nul[131]+0.1530931089239486*fl[10]*nul[131]+0.2651650429449552*fl[31]*nul[129]+0.1530931089239486*fl[26]*nul[129]+0.2651650429449552*fl[30]*nul[128]+0.1530931089239486*fl[19]*nul[128]; 
  incr2[31] = 0.2651650429449552*fl[24]*nul[135]+0.1530931089239486*fl[10]*nul[135]+0.2651650429449552*fl[28]*nul[131]+0.1530931089239486*fl[17]*nul[131]+0.2651650429449552*fl[30]*nul[129]+0.1530931089239486*fl[19]*nul[129]+0.2651650429449552*fl[31]*nul[128]+0.1530931089239486*fl[26]*nul[128]; 

  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[12] += -1.0*incr2[12]*rdxFl; 
  outl[13] += -1.0*incr2[13]*rdxFl; 
  outl[14] += -1.0*incr2[14]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 
  outl[20] += -1.0*incr2[20]*rdxFl; 
  outl[21] += -1.0*incr2[21]*rdxFl; 
  outl[22] += -1.0*incr2[22]*rdxFl; 
  outl[23] += -1.0*incr2[23]*rdxFl; 
  outl[24] += -1.0*incr2[24]*rdxFl; 
  outl[25] += -1.0*incr2[25]*rdxFl; 
  outl[27] += -1.0*incr2[27]*rdxFl; 
  outl[28] += -1.0*incr2[28]*rdxFl; 
  outl[29] += -1.0*incr2[29]*rdxFl; 
  outl[30] += -1.0*incr2[30]*rdxFl; 
  outl[31] += -1.0*incr2[31]*rdxFl; 

  }

} 
