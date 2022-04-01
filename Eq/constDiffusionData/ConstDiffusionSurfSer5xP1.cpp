#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf5xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[32]; 
  double incr2[32]; 

  incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[2] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[8] = 0.5412658773652741*fr[16]+0.5412658773652741*fl[16]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[10] = 0.5412658773652741*fr[17]+0.5412658773652741*fl[17]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = 0.5412658773652741*fr[18]+0.5412658773652741*fl[18]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[13] = 0.5412658773652741*fr[20]+0.5412658773652741*fl[20]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = 0.5412658773652741*fr[21]+0.5412658773652741*fl[21]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = 0.5412658773652741*fr[23]+0.5412658773652741*fl[23]-0.5625*fr[15]+0.5625*fl[15]; 
  incr1[16] = (-0.9375*fr[16])-0.9375*fl[16]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[17] = (-0.9375*fr[17])-0.9375*fl[17]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[18] = (-0.9375*fr[18])-0.9375*fl[18]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 
  incr1[19] = 0.5412658773652741*fr[26]+0.5412658773652741*fl[26]-0.5625*fr[19]+0.5625*fl[19]; 
  incr1[20] = (-0.9375*fr[20])-0.9375*fl[20]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 
  incr1[21] = (-0.9375*fr[21])-0.9375*fl[21]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 
  incr1[22] = 0.5412658773652741*fr[27]+0.5412658773652741*fl[27]-0.5625*fr[22]+0.5625*fl[22]; 
  incr1[23] = (-0.9375*fr[23])-0.9375*fl[23]+0.9742785792574932*fr[15]-0.9742785792574932*fl[15]; 
  incr1[24] = 0.5412658773652741*fr[28]+0.5412658773652741*fl[28]-0.5625*fr[24]+0.5625*fl[24]; 
  incr1[25] = 0.5412658773652741*fr[29]+0.5412658773652741*fl[29]-0.5625*fr[25]+0.5625*fl[25]; 
  incr1[26] = (-0.9375*fr[26])-0.9375*fl[26]+0.9742785792574932*fr[19]-0.9742785792574932*fl[19]; 
  incr1[27] = (-0.9375*fr[27])-0.9375*fl[27]+0.9742785792574932*fr[22]-0.9742785792574932*fl[22]; 
  incr1[28] = (-0.9375*fr[28])-0.9375*fl[28]+0.9742785792574932*fr[24]-0.9742785792574932*fl[24]; 
  incr1[29] = (-0.9375*fr[29])-0.9375*fl[29]+0.9742785792574932*fr[25]-0.9742785792574932*fl[25]; 
  incr1[30] = 0.5412658773652741*fr[31]+0.5412658773652741*fl[31]-0.5625*fr[30]+0.5625*fl[30]; 
  incr1[31] = (-0.9375*fr[31])-0.9375*fl[31]+0.9742785792574932*fr[30]-0.9742785792574932*fl[30]; 

  incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[16] = (-0.5*fr[16])+0.5*fl[16]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[17] = (-0.5*fr[17])+0.5*fl[17]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[18] = (-0.5*fr[18])+0.5*fl[18]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 
  incr2[20] = (-0.5*fr[20])+0.5*fl[20]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 
  incr2[21] = (-0.5*fr[21])+0.5*fl[21]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 
  incr2[23] = (-0.5*fr[23])+0.5*fl[23]+0.4330127018922193*fr[15]+0.4330127018922193*fl[15]; 
  incr2[26] = (-0.5*fr[26])+0.5*fl[26]+0.4330127018922193*fr[19]+0.4330127018922193*fl[19]; 
  incr2[27] = (-0.5*fr[27])+0.5*fl[27]+0.4330127018922193*fr[22]+0.4330127018922193*fl[22]; 
  incr2[28] = (-0.5*fr[28])+0.5*fl[28]+0.4330127018922193*fr[24]+0.4330127018922193*fl[24]; 
  incr2[29] = (-0.5*fr[29])+0.5*fl[29]+0.4330127018922193*fr[25]+0.4330127018922193*fl[25]; 
  incr2[31] = (-0.5*fr[31])+0.5*fl[31]+0.4330127018922193*fr[30]+0.4330127018922193*fl[30]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr1[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur+incr1[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur+incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur+incr1[9]*rdxFnur; 
  outr[10] += incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 
  outr[15] += incr1[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur+incr1[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur+incr1[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur+incr1[18]*rdxFnur; 
  outr[19] += incr1[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur+incr1[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur+incr1[21]*rdxFnur; 
  outr[22] += incr1[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur+incr1[23]*rdxFnur; 
  outr[24] += incr1[24]*rdxFnur; 
  outr[25] += incr1[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur+incr1[26]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur+incr1[27]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur+incr1[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur+incr1[29]*rdxFnur; 
  outr[30] += incr1[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur+incr1[31]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul-1.0*incr2[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul-1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr1[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul-1.0*incr2[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 
  outl[15] += -1.0*incr1[15]*rdxFnul; 
  outl[16] += incr1[16]*rdxFnul-1.0*incr2[16]*rdxFnul; 
  outl[17] += incr1[17]*rdxFnul-1.0*incr2[17]*rdxFnul; 
  outl[18] += incr1[18]*rdxFnul-1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr1[19]*rdxFnul; 
  outl[20] += incr1[20]*rdxFnul-1.0*incr2[20]*rdxFnul; 
  outl[21] += incr1[21]*rdxFnul-1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr1[22]*rdxFnul; 
  outl[23] += incr1[23]*rdxFnul-1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr1[24]*rdxFnul; 
  outl[25] += -1.0*incr1[25]*rdxFnul; 
  outl[26] += incr1[26]*rdxFnul-1.0*incr2[26]*rdxFnul; 
  outl[27] += incr1[27]*rdxFnul-1.0*incr2[27]*rdxFnul; 
  outl[28] += incr1[28]*rdxFnul-1.0*incr2[28]*rdxFnul; 
  outl[29] += incr1[29]*rdxFnul-1.0*incr2[29]*rdxFnul; 
  outl[30] += -1.0*incr1[30]*rdxFnul; 
  outl[31] += incr1[31]*rdxFnul-1.0*incr2[31]*rdxFnul; 

} 
void ConstDiffusionSurf5xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[32]; 
  double incr2[32]; 

  incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[6]+0.5412658773652741*fl[6]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[3] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = (-0.9375*fr[6])-0.9375*fl[6]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[7] = 0.5412658773652741*fr[16]+0.5412658773652741*fl[16]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[9] = 0.5412658773652741*fr[17]+0.5412658773652741*fl[17]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[11] = 0.5412658773652741*fr[19]+0.5412658773652741*fl[19]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = 0.5412658773652741*fr[20]+0.5412658773652741*fl[20]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[14] = 0.5412658773652741*fr[22]+0.5412658773652741*fl[22]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = 0.5412658773652741*fr[24]+0.5412658773652741*fl[24]-0.5625*fr[15]+0.5625*fl[15]; 
  incr1[16] = (-0.9375*fr[16])-0.9375*fl[16]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[17] = (-0.9375*fr[17])-0.9375*fl[17]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[18] = 0.5412658773652741*fr[26]+0.5412658773652741*fl[26]-0.5625*fr[18]+0.5625*fl[18]; 
  incr1[19] = (-0.9375*fr[19])-0.9375*fl[19]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 
  incr1[20] = (-0.9375*fr[20])-0.9375*fl[20]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 
  incr1[21] = 0.5412658773652741*fr[27]+0.5412658773652741*fl[27]-0.5625*fr[21]+0.5625*fl[21]; 
  incr1[22] = (-0.9375*fr[22])-0.9375*fl[22]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 
  incr1[23] = 0.5412658773652741*fr[28]+0.5412658773652741*fl[28]-0.5625*fr[23]+0.5625*fl[23]; 
  incr1[24] = (-0.9375*fr[24])-0.9375*fl[24]+0.9742785792574932*fr[15]-0.9742785792574932*fl[15]; 
  incr1[25] = 0.5412658773652741*fr[30]+0.5412658773652741*fl[30]-0.5625*fr[25]+0.5625*fl[25]; 
  incr1[26] = (-0.9375*fr[26])-0.9375*fl[26]+0.9742785792574932*fr[18]-0.9742785792574932*fl[18]; 
  incr1[27] = (-0.9375*fr[27])-0.9375*fl[27]+0.9742785792574932*fr[21]-0.9742785792574932*fl[21]; 
  incr1[28] = (-0.9375*fr[28])-0.9375*fl[28]+0.9742785792574932*fr[23]-0.9742785792574932*fl[23]; 
  incr1[29] = 0.5412658773652741*fr[31]+0.5412658773652741*fl[31]-0.5625*fr[29]+0.5625*fl[29]; 
  incr1[30] = (-0.9375*fr[30])-0.9375*fl[30]+0.9742785792574932*fr[25]-0.9742785792574932*fl[25]; 
  incr1[31] = (-0.9375*fr[31])-0.9375*fl[31]+0.9742785792574932*fr[29]-0.9742785792574932*fl[29]; 

  incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[6] = (-0.5*fr[6])+0.5*fl[6]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[16] = (-0.5*fr[16])+0.5*fl[16]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[17] = (-0.5*fr[17])+0.5*fl[17]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[19] = (-0.5*fr[19])+0.5*fl[19]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 
  incr2[20] = (-0.5*fr[20])+0.5*fl[20]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 
  incr2[22] = (-0.5*fr[22])+0.5*fl[22]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 
  incr2[24] = (-0.5*fr[24])+0.5*fl[24]+0.4330127018922193*fr[15]+0.4330127018922193*fl[15]; 
  incr2[26] = (-0.5*fr[26])+0.5*fl[26]+0.4330127018922193*fr[18]+0.4330127018922193*fl[18]; 
  incr2[27] = (-0.5*fr[27])+0.5*fl[27]+0.4330127018922193*fr[21]+0.4330127018922193*fl[21]; 
  incr2[28] = (-0.5*fr[28])+0.5*fl[28]+0.4330127018922193*fr[23]+0.4330127018922193*fl[23]; 
  incr2[30] = (-0.5*fr[30])+0.5*fl[30]+0.4330127018922193*fr[25]+0.4330127018922193*fl[25]; 
  incr2[31] = (-0.5*fr[31])+0.5*fl[31]+0.4330127018922193*fr[29]+0.4330127018922193*fl[29]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr1[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur+incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur+incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 
  outr[15] += incr1[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur+incr1[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur+incr1[17]*rdxFnur; 
  outr[18] += incr1[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur+incr1[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur+incr1[20]*rdxFnur; 
  outr[21] += incr1[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur+incr1[22]*rdxFnur; 
  outr[23] += incr1[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur+incr1[24]*rdxFnur; 
  outr[25] += incr1[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur+incr1[26]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur+incr1[27]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur+incr1[28]*rdxFnur; 
  outr[29] += incr1[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur+incr1[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur+incr1[31]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul-1.0*incr2[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul-1.0*incr2[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 
  outl[15] += -1.0*incr1[15]*rdxFnul; 
  outl[16] += incr1[16]*rdxFnul-1.0*incr2[16]*rdxFnul; 
  outl[17] += incr1[17]*rdxFnul-1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr1[18]*rdxFnul; 
  outl[19] += incr1[19]*rdxFnul-1.0*incr2[19]*rdxFnul; 
  outl[20] += incr1[20]*rdxFnul-1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr1[21]*rdxFnul; 
  outl[22] += incr1[22]*rdxFnul-1.0*incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr1[23]*rdxFnul; 
  outl[24] += incr1[24]*rdxFnul-1.0*incr2[24]*rdxFnul; 
  outl[25] += -1.0*incr1[25]*rdxFnul; 
  outl[26] += incr1[26]*rdxFnul-1.0*incr2[26]*rdxFnul; 
  outl[27] += incr1[27]*rdxFnul-1.0*incr2[27]*rdxFnul; 
  outl[28] += incr1[28]*rdxFnul-1.0*incr2[28]*rdxFnul; 
  outl[29] += -1.0*incr1[29]*rdxFnul; 
  outl[30] += incr1[30]*rdxFnul-1.0*incr2[30]*rdxFnul; 
  outl[31] += incr1[31]*rdxFnul-1.0*incr2[31]*rdxFnul; 

} 
void ConstDiffusionSurf5xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[2]/(dxl[2]*dxl[2]); 
  double rdxFnur = 4.0*nu[2]/(dxr[2]*dxr[2]); 

  double incr1[32]; 
  double incr2[32]; 

  incr1[0] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[7]+0.5412658773652741*fl[7]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[8]+0.5412658773652741*fl[8]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[4] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = 0.5412658773652741*fr[16]+0.5412658773652741*fl[16]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = (-0.9375*fr[7])-0.9375*fl[7]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[8] = (-0.9375*fr[8])-0.9375*fl[8]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[9] = 0.5412658773652741*fr[18]+0.5412658773652741*fl[18]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = 0.5412658773652741*fr[19]+0.5412658773652741*fl[19]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[12] = 0.5412658773652741*fr[21]+0.5412658773652741*fl[21]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = 0.5412658773652741*fr[22]+0.5412658773652741*fl[22]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[15] = 0.5412658773652741*fr[25]+0.5412658773652741*fl[25]-0.5625*fr[15]+0.5625*fl[15]; 
  incr1[16] = (-0.9375*fr[16])-0.9375*fl[16]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[17] = 0.5412658773652741*fr[26]+0.5412658773652741*fl[26]-0.5625*fr[17]+0.5625*fl[17]; 
  incr1[18] = (-0.9375*fr[18])-0.9375*fl[18]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[19] = (-0.9375*fr[19])-0.9375*fl[19]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[20] = 0.5412658773652741*fr[27]+0.5412658773652741*fl[27]-0.5625*fr[20]+0.5625*fl[20]; 
  incr1[21] = (-0.9375*fr[21])-0.9375*fl[21]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 
  incr1[22] = (-0.9375*fr[22])-0.9375*fl[22]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 
  incr1[23] = 0.5412658773652741*fr[29]+0.5412658773652741*fl[29]-0.5625*fr[23]+0.5625*fl[23]; 
  incr1[24] = 0.5412658773652741*fr[30]+0.5412658773652741*fl[30]-0.5625*fr[24]+0.5625*fl[24]; 
  incr1[25] = (-0.9375*fr[25])-0.9375*fl[25]+0.9742785792574932*fr[15]-0.9742785792574932*fl[15]; 
  incr1[26] = (-0.9375*fr[26])-0.9375*fl[26]+0.9742785792574932*fr[17]-0.9742785792574932*fl[17]; 
  incr1[27] = (-0.9375*fr[27])-0.9375*fl[27]+0.9742785792574932*fr[20]-0.9742785792574932*fl[20]; 
  incr1[28] = 0.5412658773652741*fr[31]+0.5412658773652741*fl[31]-0.5625*fr[28]+0.5625*fl[28]; 
  incr1[29] = (-0.9375*fr[29])-0.9375*fl[29]+0.9742785792574932*fr[23]-0.9742785792574932*fl[23]; 
  incr1[30] = (-0.9375*fr[30])-0.9375*fl[30]+0.9742785792574932*fr[24]-0.9742785792574932*fl[24]; 
  incr1[31] = (-0.9375*fr[31])-0.9375*fl[31]+0.9742785792574932*fr[28]-0.9742785792574932*fl[28]; 

  incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[7] = (-0.5*fr[7])+0.5*fl[7]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[8] = (-0.5*fr[8])+0.5*fl[8]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[16] = (-0.5*fr[16])+0.5*fl[16]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[18] = (-0.5*fr[18])+0.5*fl[18]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[19] = (-0.5*fr[19])+0.5*fl[19]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[21] = (-0.5*fr[21])+0.5*fl[21]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 
  incr2[22] = (-0.5*fr[22])+0.5*fl[22]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 
  incr2[25] = (-0.5*fr[25])+0.5*fl[25]+0.4330127018922193*fr[15]+0.4330127018922193*fl[15]; 
  incr2[26] = (-0.5*fr[26])+0.5*fl[26]+0.4330127018922193*fr[17]+0.4330127018922193*fl[17]; 
  incr2[27] = (-0.5*fr[27])+0.5*fl[27]+0.4330127018922193*fr[20]+0.4330127018922193*fl[20]; 
  incr2[29] = (-0.5*fr[29])+0.5*fl[29]+0.4330127018922193*fr[23]+0.4330127018922193*fl[23]; 
  incr2[30] = (-0.5*fr[30])+0.5*fl[30]+0.4330127018922193*fr[24]+0.4330127018922193*fl[24]; 
  incr2[31] = (-0.5*fr[31])+0.5*fl[31]+0.4330127018922193*fr[28]+0.4330127018922193*fl[28]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr1[5]*rdxFnur; 
  outr[6] += incr1[6]*rdxFnur; 
  outr[7] += incr2[7]*rdxFnur+incr1[7]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur+incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 
  outr[10] += incr1[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur+incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur+incr1[14]*rdxFnur; 
  outr[15] += incr1[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur+incr1[16]*rdxFnur; 
  outr[17] += incr1[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur+incr1[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur+incr1[19]*rdxFnur; 
  outr[20] += incr1[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur+incr1[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur+incr1[22]*rdxFnur; 
  outr[23] += incr1[23]*rdxFnur; 
  outr[24] += incr1[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur+incr1[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur+incr1[26]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur+incr1[27]*rdxFnur; 
  outr[28] += incr1[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur+incr1[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur+incr1[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur+incr1[31]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += -1.0*incr1[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul-1.0*incr2[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul-1.0*incr2[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 
  outl[10] += -1.0*incr1[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul-1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 
  outl[15] += -1.0*incr1[15]*rdxFnul; 
  outl[16] += incr1[16]*rdxFnul-1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr1[17]*rdxFnul; 
  outl[18] += incr1[18]*rdxFnul-1.0*incr2[18]*rdxFnul; 
  outl[19] += incr1[19]*rdxFnul-1.0*incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr1[20]*rdxFnul; 
  outl[21] += incr1[21]*rdxFnul-1.0*incr2[21]*rdxFnul; 
  outl[22] += incr1[22]*rdxFnul-1.0*incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr1[23]*rdxFnul; 
  outl[24] += -1.0*incr1[24]*rdxFnul; 
  outl[25] += incr1[25]*rdxFnul-1.0*incr2[25]*rdxFnul; 
  outl[26] += incr1[26]*rdxFnul-1.0*incr2[26]*rdxFnul; 
  outl[27] += incr1[27]*rdxFnul-1.0*incr2[27]*rdxFnul; 
  outl[28] += -1.0*incr1[28]*rdxFnul; 
  outl[29] += incr1[29]*rdxFnul-1.0*incr2[29]*rdxFnul; 
  outl[30] += incr1[30]*rdxFnul-1.0*incr2[30]*rdxFnul; 
  outl[31] += incr1[31]*rdxFnul-1.0*incr2[31]*rdxFnul; 

} 
void ConstDiffusionSurf5xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[3]/(dxl[3]*dxl[3]); 
  double rdxFnur = 4.0*nu[3]/(dxr[3]*dxr[3]); 

  double incr1[32]; 
  double incr2[32]; 

  incr1[0] = 0.5412658773652741*fr[4]+0.5412658773652741*fl[4]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[9]+0.5412658773652741*fl[9]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[10]+0.5412658773652741*fl[10]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[11]+0.5412658773652741*fl[11]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = (-0.9375*fr[4])-0.9375*fl[4]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[5] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[5]+0.5625*fl[5]; 
  incr1[6] = 0.5412658773652741*fr[17]+0.5412658773652741*fl[17]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = 0.5412658773652741*fr[18]+0.5412658773652741*fl[18]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = 0.5412658773652741*fr[19]+0.5412658773652741*fl[19]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = (-0.9375*fr[9])-0.9375*fl[9]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[10] = (-0.9375*fr[10])-0.9375*fl[10]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[11] = (-0.9375*fr[11])-0.9375*fl[11]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[12] = 0.5412658773652741*fr[23]+0.5412658773652741*fl[23]-0.5625*fr[12]+0.5625*fl[12]; 
  incr1[13] = 0.5412658773652741*fr[24]+0.5412658773652741*fl[24]-0.5625*fr[13]+0.5625*fl[13]; 
  incr1[14] = 0.5412658773652741*fr[25]+0.5412658773652741*fl[25]-0.5625*fr[14]+0.5625*fl[14]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[5]-0.9742785792574932*fl[5]; 
  incr1[16] = 0.5412658773652741*fr[26]+0.5412658773652741*fl[26]-0.5625*fr[16]+0.5625*fl[16]; 
  incr1[17] = (-0.9375*fr[17])-0.9375*fl[17]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[18] = (-0.9375*fr[18])-0.9375*fl[18]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[19] = (-0.9375*fr[19])-0.9375*fl[19]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[20] = 0.5412658773652741*fr[28]+0.5412658773652741*fl[28]-0.5625*fr[20]+0.5625*fl[20]; 
  incr1[21] = 0.5412658773652741*fr[29]+0.5412658773652741*fl[29]-0.5625*fr[21]+0.5625*fl[21]; 
  incr1[22] = 0.5412658773652741*fr[30]+0.5412658773652741*fl[30]-0.5625*fr[22]+0.5625*fl[22]; 
  incr1[23] = (-0.9375*fr[23])-0.9375*fl[23]+0.9742785792574932*fr[12]-0.9742785792574932*fl[12]; 
  incr1[24] = (-0.9375*fr[24])-0.9375*fl[24]+0.9742785792574932*fr[13]-0.9742785792574932*fl[13]; 
  incr1[25] = (-0.9375*fr[25])-0.9375*fl[25]+0.9742785792574932*fr[14]-0.9742785792574932*fl[14]; 
  incr1[26] = (-0.9375*fr[26])-0.9375*fl[26]+0.9742785792574932*fr[16]-0.9742785792574932*fl[16]; 
  incr1[27] = 0.5412658773652741*fr[31]+0.5412658773652741*fl[31]-0.5625*fr[27]+0.5625*fl[27]; 
  incr1[28] = (-0.9375*fr[28])-0.9375*fl[28]+0.9742785792574932*fr[20]-0.9742785792574932*fl[20]; 
  incr1[29] = (-0.9375*fr[29])-0.9375*fl[29]+0.9742785792574932*fr[21]-0.9742785792574932*fl[21]; 
  incr1[30] = (-0.9375*fr[30])-0.9375*fl[30]+0.9742785792574932*fr[22]-0.9742785792574932*fl[22]; 
  incr1[31] = (-0.9375*fr[31])-0.9375*fl[31]+0.9742785792574932*fr[27]-0.9742785792574932*fl[27]; 

  incr2[4] = (-0.5*fr[4])+0.5*fl[4]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[9] = (-0.5*fr[9])+0.5*fl[9]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[10] = (-0.5*fr[10])+0.5*fl[10]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[11] = (-0.5*fr[11])+0.5*fl[11]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[5]+0.4330127018922193*fl[5]; 
  incr2[17] = (-0.5*fr[17])+0.5*fl[17]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[18] = (-0.5*fr[18])+0.5*fl[18]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[19] = (-0.5*fr[19])+0.5*fl[19]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[23] = (-0.5*fr[23])+0.5*fl[23]+0.4330127018922193*fr[12]+0.4330127018922193*fl[12]; 
  incr2[24] = (-0.5*fr[24])+0.5*fl[24]+0.4330127018922193*fr[13]+0.4330127018922193*fl[13]; 
  incr2[25] = (-0.5*fr[25])+0.5*fl[25]+0.4330127018922193*fr[14]+0.4330127018922193*fl[14]; 
  incr2[26] = (-0.5*fr[26])+0.5*fl[26]+0.4330127018922193*fr[16]+0.4330127018922193*fl[16]; 
  incr2[28] = (-0.5*fr[28])+0.5*fl[28]+0.4330127018922193*fr[20]+0.4330127018922193*fl[20]; 
  incr2[29] = (-0.5*fr[29])+0.5*fl[29]+0.4330127018922193*fr[21]+0.4330127018922193*fl[21]; 
  incr2[30] = (-0.5*fr[30])+0.5*fl[30]+0.4330127018922193*fr[22]+0.4330127018922193*fl[22]; 
  incr2[31] = (-0.5*fr[31])+0.5*fl[31]+0.4330127018922193*fr[27]+0.4330127018922193*fl[27]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur+incr1[4]*rdxFnur; 
  outr[5] += incr1[5]*rdxFnur; 
  outr[6] += incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur+incr1[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur+incr1[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur+incr1[11]*rdxFnur; 
  outr[12] += incr1[12]*rdxFnur; 
  outr[13] += incr1[13]*rdxFnur; 
  outr[14] += incr1[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur+incr1[15]*rdxFnur; 
  outr[16] += incr1[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur+incr1[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur+incr1[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur+incr1[19]*rdxFnur; 
  outr[20] += incr1[20]*rdxFnur; 
  outr[21] += incr1[21]*rdxFnur; 
  outr[22] += incr1[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur+incr1[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur+incr1[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur+incr1[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur+incr1[26]*rdxFnur; 
  outr[27] += incr1[27]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur+incr1[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur+incr1[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur+incr1[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur+incr1[31]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul-1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr1[5]*rdxFnul; 
  outl[6] += -1.0*incr1[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul-1.0*incr2[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul-1.0*incr2[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul-1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr1[12]*rdxFnul; 
  outl[13] += -1.0*incr1[13]*rdxFnul; 
  outl[14] += -1.0*incr1[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul-1.0*incr2[15]*rdxFnul; 
  outl[16] += -1.0*incr1[16]*rdxFnul; 
  outl[17] += incr1[17]*rdxFnul-1.0*incr2[17]*rdxFnul; 
  outl[18] += incr1[18]*rdxFnul-1.0*incr2[18]*rdxFnul; 
  outl[19] += incr1[19]*rdxFnul-1.0*incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr1[20]*rdxFnul; 
  outl[21] += -1.0*incr1[21]*rdxFnul; 
  outl[22] += -1.0*incr1[22]*rdxFnul; 
  outl[23] += incr1[23]*rdxFnul-1.0*incr2[23]*rdxFnul; 
  outl[24] += incr1[24]*rdxFnul-1.0*incr2[24]*rdxFnul; 
  outl[25] += incr1[25]*rdxFnul-1.0*incr2[25]*rdxFnul; 
  outl[26] += incr1[26]*rdxFnul-1.0*incr2[26]*rdxFnul; 
  outl[27] += -1.0*incr1[27]*rdxFnul; 
  outl[28] += incr1[28]*rdxFnul-1.0*incr2[28]*rdxFnul; 
  outl[29] += incr1[29]*rdxFnul-1.0*incr2[29]*rdxFnul; 
  outl[30] += incr1[30]*rdxFnul-1.0*incr2[30]*rdxFnul; 
  outl[31] += incr1[31]*rdxFnul-1.0*incr2[31]*rdxFnul; 

} 
void ConstDiffusionSurf5xSerP1_X5(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[4]/(dxl[4]*dxl[4]); 
  double rdxFnur = 4.0*nu[4]/(dxr[4]*dxr[4]); 

  double incr1[32]; 
  double incr2[32]; 

  incr1[0] = 0.5412658773652741*fr[5]+0.5412658773652741*fl[5]-0.5625*fr[0]+0.5625*fl[0]; 
  incr1[1] = 0.5412658773652741*fr[12]+0.5412658773652741*fl[12]-0.5625*fr[1]+0.5625*fl[1]; 
  incr1[2] = 0.5412658773652741*fr[13]+0.5412658773652741*fl[13]-0.5625*fr[2]+0.5625*fl[2]; 
  incr1[3] = 0.5412658773652741*fr[14]+0.5412658773652741*fl[14]-0.5625*fr[3]+0.5625*fl[3]; 
  incr1[4] = 0.5412658773652741*fr[15]+0.5412658773652741*fl[15]-0.5625*fr[4]+0.5625*fl[4]; 
  incr1[5] = (-0.9375*fr[5])-0.9375*fl[5]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
  incr1[6] = 0.5412658773652741*fr[20]+0.5412658773652741*fl[20]-0.5625*fr[6]+0.5625*fl[6]; 
  incr1[7] = 0.5412658773652741*fr[21]+0.5412658773652741*fl[21]-0.5625*fr[7]+0.5625*fl[7]; 
  incr1[8] = 0.5412658773652741*fr[22]+0.5412658773652741*fl[22]-0.5625*fr[8]+0.5625*fl[8]; 
  incr1[9] = 0.5412658773652741*fr[23]+0.5412658773652741*fl[23]-0.5625*fr[9]+0.5625*fl[9]; 
  incr1[10] = 0.5412658773652741*fr[24]+0.5412658773652741*fl[24]-0.5625*fr[10]+0.5625*fl[10]; 
  incr1[11] = 0.5412658773652741*fr[25]+0.5412658773652741*fl[25]-0.5625*fr[11]+0.5625*fl[11]; 
  incr1[12] = (-0.9375*fr[12])-0.9375*fl[12]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 
  incr1[13] = (-0.9375*fr[13])-0.9375*fl[13]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 
  incr1[14] = (-0.9375*fr[14])-0.9375*fl[14]+0.9742785792574932*fr[3]-0.9742785792574932*fl[3]; 
  incr1[15] = (-0.9375*fr[15])-0.9375*fl[15]+0.9742785792574932*fr[4]-0.9742785792574932*fl[4]; 
  incr1[16] = 0.5412658773652741*fr[27]+0.5412658773652741*fl[27]-0.5625*fr[16]+0.5625*fl[16]; 
  incr1[17] = 0.5412658773652741*fr[28]+0.5412658773652741*fl[28]-0.5625*fr[17]+0.5625*fl[17]; 
  incr1[18] = 0.5412658773652741*fr[29]+0.5412658773652741*fl[29]-0.5625*fr[18]+0.5625*fl[18]; 
  incr1[19] = 0.5412658773652741*fr[30]+0.5412658773652741*fl[30]-0.5625*fr[19]+0.5625*fl[19]; 
  incr1[20] = (-0.9375*fr[20])-0.9375*fl[20]+0.9742785792574932*fr[6]-0.9742785792574932*fl[6]; 
  incr1[21] = (-0.9375*fr[21])-0.9375*fl[21]+0.9742785792574932*fr[7]-0.9742785792574932*fl[7]; 
  incr1[22] = (-0.9375*fr[22])-0.9375*fl[22]+0.9742785792574932*fr[8]-0.9742785792574932*fl[8]; 
  incr1[23] = (-0.9375*fr[23])-0.9375*fl[23]+0.9742785792574932*fr[9]-0.9742785792574932*fl[9]; 
  incr1[24] = (-0.9375*fr[24])-0.9375*fl[24]+0.9742785792574932*fr[10]-0.9742785792574932*fl[10]; 
  incr1[25] = (-0.9375*fr[25])-0.9375*fl[25]+0.9742785792574932*fr[11]-0.9742785792574932*fl[11]; 
  incr1[26] = 0.5412658773652741*fr[31]+0.5412658773652741*fl[31]-0.5625*fr[26]+0.5625*fl[26]; 
  incr1[27] = (-0.9375*fr[27])-0.9375*fl[27]+0.9742785792574932*fr[16]-0.9742785792574932*fl[16]; 
  incr1[28] = (-0.9375*fr[28])-0.9375*fl[28]+0.9742785792574932*fr[17]-0.9742785792574932*fl[17]; 
  incr1[29] = (-0.9375*fr[29])-0.9375*fl[29]+0.9742785792574932*fr[18]-0.9742785792574932*fl[18]; 
  incr1[30] = (-0.9375*fr[30])-0.9375*fl[30]+0.9742785792574932*fr[19]-0.9742785792574932*fl[19]; 
  incr1[31] = (-0.9375*fr[31])-0.9375*fl[31]+0.9742785792574932*fr[26]-0.9742785792574932*fl[26]; 

  incr2[5] = (-0.5*fr[5])+0.5*fl[5]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
  incr2[12] = (-0.5*fr[12])+0.5*fl[12]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 
  incr2[13] = (-0.5*fr[13])+0.5*fl[13]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 
  incr2[14] = (-0.5*fr[14])+0.5*fl[14]+0.4330127018922193*fr[3]+0.4330127018922193*fl[3]; 
  incr2[15] = (-0.5*fr[15])+0.5*fl[15]+0.4330127018922193*fr[4]+0.4330127018922193*fl[4]; 
  incr2[20] = (-0.5*fr[20])+0.5*fl[20]+0.4330127018922193*fr[6]+0.4330127018922193*fl[6]; 
  incr2[21] = (-0.5*fr[21])+0.5*fl[21]+0.4330127018922193*fr[7]+0.4330127018922193*fl[7]; 
  incr2[22] = (-0.5*fr[22])+0.5*fl[22]+0.4330127018922193*fr[8]+0.4330127018922193*fl[8]; 
  incr2[23] = (-0.5*fr[23])+0.5*fl[23]+0.4330127018922193*fr[9]+0.4330127018922193*fl[9]; 
  incr2[24] = (-0.5*fr[24])+0.5*fl[24]+0.4330127018922193*fr[10]+0.4330127018922193*fl[10]; 
  incr2[25] = (-0.5*fr[25])+0.5*fl[25]+0.4330127018922193*fr[11]+0.4330127018922193*fl[11]; 
  incr2[27] = (-0.5*fr[27])+0.5*fl[27]+0.4330127018922193*fr[16]+0.4330127018922193*fl[16]; 
  incr2[28] = (-0.5*fr[28])+0.5*fl[28]+0.4330127018922193*fr[17]+0.4330127018922193*fl[17]; 
  incr2[29] = (-0.5*fr[29])+0.5*fl[29]+0.4330127018922193*fr[18]+0.4330127018922193*fl[18]; 
  incr2[30] = (-0.5*fr[30])+0.5*fl[30]+0.4330127018922193*fr[19]+0.4330127018922193*fl[19]; 
  incr2[31] = (-0.5*fr[31])+0.5*fl[31]+0.4330127018922193*fr[26]+0.4330127018922193*fl[26]; 

  outr[0] += incr1[0]*rdxFnur; 
  outr[1] += incr1[1]*rdxFnur; 
  outr[2] += incr1[2]*rdxFnur; 
  outr[3] += incr1[3]*rdxFnur; 
  outr[4] += incr1[4]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur+incr1[5]*rdxFnur; 
  outr[6] += incr1[6]*rdxFnur; 
  outr[7] += incr1[7]*rdxFnur; 
  outr[8] += incr1[8]*rdxFnur; 
  outr[9] += incr1[9]*rdxFnur; 
  outr[10] += incr1[10]*rdxFnur; 
  outr[11] += incr1[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur+incr1[12]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur+incr1[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur+incr1[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur+incr1[15]*rdxFnur; 
  outr[16] += incr1[16]*rdxFnur; 
  outr[17] += incr1[17]*rdxFnur; 
  outr[18] += incr1[18]*rdxFnur; 
  outr[19] += incr1[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur+incr1[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur+incr1[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur+incr1[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur+incr1[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur+incr1[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur+incr1[25]*rdxFnur; 
  outr[26] += incr1[26]*rdxFnur; 
  outr[27] += incr2[27]*rdxFnur+incr1[27]*rdxFnur; 
  outr[28] += incr2[28]*rdxFnur+incr1[28]*rdxFnur; 
  outr[29] += incr2[29]*rdxFnur+incr1[29]*rdxFnur; 
  outr[30] += incr2[30]*rdxFnur+incr1[30]*rdxFnur; 
  outr[31] += incr2[31]*rdxFnur+incr1[31]*rdxFnur; 

  outl[0] += -1.0*incr1[0]*rdxFnul; 
  outl[1] += -1.0*incr1[1]*rdxFnul; 
  outl[2] += -1.0*incr1[2]*rdxFnul; 
  outl[3] += -1.0*incr1[3]*rdxFnul; 
  outl[4] += -1.0*incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul-1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr1[6]*rdxFnul; 
  outl[7] += -1.0*incr1[7]*rdxFnul; 
  outl[8] += -1.0*incr1[8]*rdxFnul; 
  outl[9] += -1.0*incr1[9]*rdxFnul; 
  outl[10] += -1.0*incr1[10]*rdxFnul; 
  outl[11] += -1.0*incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul-1.0*incr2[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul-1.0*incr2[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul-1.0*incr2[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul-1.0*incr2[15]*rdxFnul; 
  outl[16] += -1.0*incr1[16]*rdxFnul; 
  outl[17] += -1.0*incr1[17]*rdxFnul; 
  outl[18] += -1.0*incr1[18]*rdxFnul; 
  outl[19] += -1.0*incr1[19]*rdxFnul; 
  outl[20] += incr1[20]*rdxFnul-1.0*incr2[20]*rdxFnul; 
  outl[21] += incr1[21]*rdxFnul-1.0*incr2[21]*rdxFnul; 
  outl[22] += incr1[22]*rdxFnul-1.0*incr2[22]*rdxFnul; 
  outl[23] += incr1[23]*rdxFnul-1.0*incr2[23]*rdxFnul; 
  outl[24] += incr1[24]*rdxFnul-1.0*incr2[24]*rdxFnul; 
  outl[25] += incr1[25]*rdxFnul-1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr1[26]*rdxFnul; 
  outl[27] += incr1[27]*rdxFnul-1.0*incr2[27]*rdxFnul; 
  outl[28] += incr1[28]*rdxFnul-1.0*incr2[28]*rdxFnul; 
  outl[29] += incr1[29]*rdxFnul-1.0*incr2[29]*rdxFnul; 
  outl[30] += incr1[30]*rdxFnul-1.0*incr2[30]*rdxFnul; 
  outl[31] += incr1[31]*rdxFnul-1.0*incr2[31]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf5xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[0]/(dxl[0]*dxl[0]*dxl[0]*dxl[0]); 
  double rdxFnur = 16.0*nu[0]/(dxr[0]*dxr[0]*dxr[0]*dxr[0]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  incr1[0] = (-1.623797632095822*fr[1])-1.623797632095822*fl[1]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = 2.8125*fr[1]+2.8125*fl[1]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[2] = (-1.623797632095822*fr[6])-1.623797632095822*fl[6]+0.9375*fr[2]-0.9375*fl[2]; 
  incr1[3] = (-1.623797632095822*fr[7])-1.623797632095822*fl[7]+0.9375*fr[3]-0.9375*fl[3]; 
  incr1[4] = (-1.623797632095822*fr[9])-1.623797632095822*fl[9]+0.9375*fr[4]-0.9375*fl[4]; 
  incr1[5] = (-1.623797632095822*fr[12])-1.623797632095822*fl[12]+0.9375*fr[5]-0.9375*fl[5]; 
  incr1[6] = 2.8125*fr[6]+2.8125*fl[6]-1.623797632095822*fr[2]+1.623797632095822*fl[2]; 
  incr1[7] = 2.8125*fr[7]+2.8125*fl[7]-1.623797632095822*fr[3]+1.623797632095822*fl[3]; 
  incr1[8] = (-1.623797632095822*fr[16])-1.623797632095822*fl[16]+0.9375*fr[8]-0.9375*fl[8]; 
  incr1[9] = 2.8125*fr[9]+2.8125*fl[9]-1.623797632095822*fr[4]+1.623797632095822*fl[4]; 
  incr1[10] = (-1.623797632095822*fr[17])-1.623797632095822*fl[17]+0.9375*fr[10]-0.9375*fl[10]; 
  incr1[11] = (-1.623797632095822*fr[18])-1.623797632095822*fl[18]+0.9375*fr[11]-0.9375*fl[11]; 
  incr1[12] = 2.8125*fr[12]+2.8125*fl[12]-1.623797632095822*fr[5]+1.623797632095822*fl[5]; 
  incr1[13] = (-1.623797632095822*fr[20])-1.623797632095822*fl[20]+0.9375*fr[13]-0.9375*fl[13]; 
  incr1[14] = (-1.623797632095822*fr[21])-1.623797632095822*fl[21]+0.9375*fr[14]-0.9375*fl[14]; 
  incr1[15] = (-1.623797632095822*fr[23])-1.623797632095822*fl[23]+0.9375*fr[15]-0.9375*fl[15]; 
  incr1[16] = 2.8125*fr[16]+2.8125*fl[16]-1.623797632095822*fr[8]+1.623797632095822*fl[8]; 
  incr1[17] = 2.8125*fr[17]+2.8125*fl[17]-1.623797632095822*fr[10]+1.623797632095822*fl[10]; 
  incr1[18] = 2.8125*fr[18]+2.8125*fl[18]-1.623797632095822*fr[11]+1.623797632095822*fl[11]; 
  incr1[19] = (-1.623797632095822*fr[26])-1.623797632095822*fl[26]+0.9375*fr[19]-0.9375*fl[19]; 
  incr1[20] = 2.8125*fr[20]+2.8125*fl[20]-1.623797632095822*fr[13]+1.623797632095822*fl[13]; 
  incr1[21] = 2.8125*fr[21]+2.8125*fl[21]-1.623797632095822*fr[14]+1.623797632095822*fl[14]; 
  incr1[22] = (-1.623797632095822*fr[27])-1.623797632095822*fl[27]+0.9375*fr[22]-0.9375*fl[22]; 
  incr1[23] = 2.8125*fr[23]+2.8125*fl[23]-1.623797632095822*fr[15]+1.623797632095822*fl[15]; 
  incr1[24] = (-1.623797632095822*fr[28])-1.623797632095822*fl[28]+0.9375*fr[24]-0.9375*fl[24]; 
  incr1[25] = (-1.623797632095822*fr[29])-1.623797632095822*fl[29]+0.9375*fr[25]-0.9375*fl[25]; 
  incr1[26] = 2.8125*fr[26]+2.8125*fl[26]-1.623797632095822*fr[19]+1.623797632095822*fl[19]; 
  incr1[27] = 2.8125*fr[27]+2.8125*fl[27]-1.623797632095822*fr[22]+1.623797632095822*fl[22]; 
  incr1[28] = 2.8125*fr[28]+2.8125*fl[28]-1.623797632095822*fr[24]+1.623797632095822*fl[24]; 
  incr1[29] = 2.8125*fr[29]+2.8125*fl[29]-1.623797632095822*fr[25]+1.623797632095822*fl[25]; 
  incr1[30] = (-1.623797632095822*fr[31])-1.623797632095822*fl[31]+0.9375*fr[30]-0.9375*fl[30]; 
  incr1[31] = 2.8125*fr[31]+2.8125*fl[31]-1.623797632095822*fr[30]+1.623797632095822*fl[30]; 

  incr2[1] = 0.75*fr[1]-0.75*fl[1]; 
  incr2[6] = 0.75*fr[6]-0.75*fl[6]; 
  incr2[7] = 0.75*fr[7]-0.75*fl[7]; 
  incr2[9] = 0.75*fr[9]-0.75*fl[9]; 
  incr2[12] = 0.75*fr[12]-0.75*fl[12]; 
  incr2[16] = 0.75*fr[16]-0.75*fl[16]; 
  incr2[17] = 0.75*fr[17]-0.75*fl[17]; 
  incr2[18] = 0.75*fr[18]-0.75*fl[18]; 
  incr2[20] = 0.75*fr[20]-0.75*fl[20]; 
  incr2[21] = 0.75*fr[21]-0.75*fl[21]; 
  incr2[23] = 0.75*fr[23]-0.75*fl[23]; 
  incr2[26] = 0.75*fr[26]-0.75*fl[26]; 
  incr2[27] = 0.75*fr[27]-0.75*fl[27]; 
  incr2[28] = 0.75*fr[28]-0.75*fl[28]; 
  incr2[29] = 0.75*fr[29]-0.75*fl[29]; 
  incr2[31] = 0.75*fr[31]-0.75*fl[31]; 



  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += (-1.0*incr2[1]*rdxFnur)-1.0*incr1[1]*rdxFnur; 
  outr[2] += -1.0*incr1[2]*rdxFnur; 
  outr[3] += -1.0*incr1[3]*rdxFnur; 
  outr[4] += -1.0*incr1[4]*rdxFnur; 
  outr[5] += -1.0*incr1[5]*rdxFnur; 
  outr[6] += (-1.0*incr2[6]*rdxFnur)-1.0*incr1[6]*rdxFnur; 
  outr[7] += (-1.0*incr2[7]*rdxFnur)-1.0*incr1[7]*rdxFnur; 
  outr[8] += -1.0*incr1[8]*rdxFnur; 
  outr[9] += (-1.0*incr2[9]*rdxFnur)-1.0*incr1[9]*rdxFnur; 
  outr[10] += -1.0*incr1[10]*rdxFnur; 
  outr[11] += -1.0*incr1[11]*rdxFnur; 
  outr[12] += (-1.0*incr2[12]*rdxFnur)-1.0*incr1[12]*rdxFnur; 
  outr[13] += -1.0*incr1[13]*rdxFnur; 
  outr[14] += -1.0*incr1[14]*rdxFnur; 
  outr[15] += -1.0*incr1[15]*rdxFnur; 
  outr[16] += (-1.0*incr2[16]*rdxFnur)-1.0*incr1[16]*rdxFnur; 
  outr[17] += (-1.0*incr2[17]*rdxFnur)-1.0*incr1[17]*rdxFnur; 
  outr[18] += (-1.0*incr2[18]*rdxFnur)-1.0*incr1[18]*rdxFnur; 
  outr[19] += -1.0*incr1[19]*rdxFnur; 
  outr[20] += (-1.0*incr2[20]*rdxFnur)-1.0*incr1[20]*rdxFnur; 
  outr[21] += (-1.0*incr2[21]*rdxFnur)-1.0*incr1[21]*rdxFnur; 
  outr[22] += -1.0*incr1[22]*rdxFnur; 
  outr[23] += (-1.0*incr2[23]*rdxFnur)-1.0*incr1[23]*rdxFnur; 
  outr[24] += -1.0*incr1[24]*rdxFnur; 
  outr[25] += -1.0*incr1[25]*rdxFnur; 
  outr[26] += (-1.0*incr2[26]*rdxFnur)-1.0*incr1[26]*rdxFnur; 
  outr[27] += (-1.0*incr2[27]*rdxFnur)-1.0*incr1[27]*rdxFnur; 
  outr[28] += (-1.0*incr2[28]*rdxFnur)-1.0*incr1[28]*rdxFnur; 
  outr[29] += (-1.0*incr2[29]*rdxFnur)-1.0*incr1[29]*rdxFnur; 
  outr[30] += -1.0*incr1[30]*rdxFnur; 
  outr[31] += (-1.0*incr2[31]*rdxFnur)-1.0*incr1[31]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr2[1]*rdxFnul-1.0*incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul-1.0*incr1[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul-1.0*incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul-1.0*incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul-1.0*incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul-1.0*incr1[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul-1.0*incr1[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul-1.0*incr1[18]*rdxFnul; 
  outl[19] += incr1[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul-1.0*incr1[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul-1.0*incr1[21]*rdxFnul; 
  outl[22] += incr1[22]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul-1.0*incr1[23]*rdxFnul; 
  outl[24] += incr1[24]*rdxFnul; 
  outl[25] += incr1[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul-1.0*incr1[26]*rdxFnul; 
  outl[27] += incr2[27]*rdxFnul-1.0*incr1[27]*rdxFnul; 
  outl[28] += incr2[28]*rdxFnul-1.0*incr1[28]*rdxFnul; 
  outl[29] += incr2[29]*rdxFnul-1.0*incr1[29]*rdxFnul; 
  outl[30] += incr1[30]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul-1.0*incr1[31]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf5xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[1]/(dxl[1]*dxl[1]*dxl[1]*dxl[1]); 
  double rdxFnur = 16.0*nu[1]/(dxr[1]*dxr[1]*dxr[1]*dxr[1]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  incr1[0] = (-1.623797632095822*fr[2])-1.623797632095822*fl[2]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = (-1.623797632095822*fr[6])-1.623797632095822*fl[6]+0.9375*fr[1]-0.9375*fl[1]; 
  incr1[2] = 2.8125*fr[2]+2.8125*fl[2]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[3] = (-1.623797632095822*fr[8])-1.623797632095822*fl[8]+0.9375*fr[3]-0.9375*fl[3]; 
  incr1[4] = (-1.623797632095822*fr[10])-1.623797632095822*fl[10]+0.9375*fr[4]-0.9375*fl[4]; 
  incr1[5] = (-1.623797632095822*fr[13])-1.623797632095822*fl[13]+0.9375*fr[5]-0.9375*fl[5]; 
  incr1[6] = 2.8125*fr[6]+2.8125*fl[6]-1.623797632095822*fr[1]+1.623797632095822*fl[1]; 
  incr1[7] = (-1.623797632095822*fr[16])-1.623797632095822*fl[16]+0.9375*fr[7]-0.9375*fl[7]; 
  incr1[8] = 2.8125*fr[8]+2.8125*fl[8]-1.623797632095822*fr[3]+1.623797632095822*fl[3]; 
  incr1[9] = (-1.623797632095822*fr[17])-1.623797632095822*fl[17]+0.9375*fr[9]-0.9375*fl[9]; 
  incr1[10] = 2.8125*fr[10]+2.8125*fl[10]-1.623797632095822*fr[4]+1.623797632095822*fl[4]; 
  incr1[11] = (-1.623797632095822*fr[19])-1.623797632095822*fl[19]+0.9375*fr[11]-0.9375*fl[11]; 
  incr1[12] = (-1.623797632095822*fr[20])-1.623797632095822*fl[20]+0.9375*fr[12]-0.9375*fl[12]; 
  incr1[13] = 2.8125*fr[13]+2.8125*fl[13]-1.623797632095822*fr[5]+1.623797632095822*fl[5]; 
  incr1[14] = (-1.623797632095822*fr[22])-1.623797632095822*fl[22]+0.9375*fr[14]-0.9375*fl[14]; 
  incr1[15] = (-1.623797632095822*fr[24])-1.623797632095822*fl[24]+0.9375*fr[15]-0.9375*fl[15]; 
  incr1[16] = 2.8125*fr[16]+2.8125*fl[16]-1.623797632095822*fr[7]+1.623797632095822*fl[7]; 
  incr1[17] = 2.8125*fr[17]+2.8125*fl[17]-1.623797632095822*fr[9]+1.623797632095822*fl[9]; 
  incr1[18] = (-1.623797632095822*fr[26])-1.623797632095822*fl[26]+0.9375*fr[18]-0.9375*fl[18]; 
  incr1[19] = 2.8125*fr[19]+2.8125*fl[19]-1.623797632095822*fr[11]+1.623797632095822*fl[11]; 
  incr1[20] = 2.8125*fr[20]+2.8125*fl[20]-1.623797632095822*fr[12]+1.623797632095822*fl[12]; 
  incr1[21] = (-1.623797632095822*fr[27])-1.623797632095822*fl[27]+0.9375*fr[21]-0.9375*fl[21]; 
  incr1[22] = 2.8125*fr[22]+2.8125*fl[22]-1.623797632095822*fr[14]+1.623797632095822*fl[14]; 
  incr1[23] = (-1.623797632095822*fr[28])-1.623797632095822*fl[28]+0.9375*fr[23]-0.9375*fl[23]; 
  incr1[24] = 2.8125*fr[24]+2.8125*fl[24]-1.623797632095822*fr[15]+1.623797632095822*fl[15]; 
  incr1[25] = (-1.623797632095822*fr[30])-1.623797632095822*fl[30]+0.9375*fr[25]-0.9375*fl[25]; 
  incr1[26] = 2.8125*fr[26]+2.8125*fl[26]-1.623797632095822*fr[18]+1.623797632095822*fl[18]; 
  incr1[27] = 2.8125*fr[27]+2.8125*fl[27]-1.623797632095822*fr[21]+1.623797632095822*fl[21]; 
  incr1[28] = 2.8125*fr[28]+2.8125*fl[28]-1.623797632095822*fr[23]+1.623797632095822*fl[23]; 
  incr1[29] = (-1.623797632095822*fr[31])-1.623797632095822*fl[31]+0.9375*fr[29]-0.9375*fl[29]; 
  incr1[30] = 2.8125*fr[30]+2.8125*fl[30]-1.623797632095822*fr[25]+1.623797632095822*fl[25]; 
  incr1[31] = 2.8125*fr[31]+2.8125*fl[31]-1.623797632095822*fr[29]+1.623797632095822*fl[29]; 

  incr2[2] = 0.75*fr[2]-0.75*fl[2]; 
  incr2[6] = 0.75*fr[6]-0.75*fl[6]; 
  incr2[8] = 0.75*fr[8]-0.75*fl[8]; 
  incr2[10] = 0.75*fr[10]-0.75*fl[10]; 
  incr2[13] = 0.75*fr[13]-0.75*fl[13]; 
  incr2[16] = 0.75*fr[16]-0.75*fl[16]; 
  incr2[17] = 0.75*fr[17]-0.75*fl[17]; 
  incr2[19] = 0.75*fr[19]-0.75*fl[19]; 
  incr2[20] = 0.75*fr[20]-0.75*fl[20]; 
  incr2[22] = 0.75*fr[22]-0.75*fl[22]; 
  incr2[24] = 0.75*fr[24]-0.75*fl[24]; 
  incr2[26] = 0.75*fr[26]-0.75*fl[26]; 
  incr2[27] = 0.75*fr[27]-0.75*fl[27]; 
  incr2[28] = 0.75*fr[28]-0.75*fl[28]; 
  incr2[30] = 0.75*fr[30]-0.75*fl[30]; 
  incr2[31] = 0.75*fr[31]-0.75*fl[31]; 



  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += -1.0*incr1[1]*rdxFnur; 
  outr[2] += (-1.0*incr2[2]*rdxFnur)-1.0*incr1[2]*rdxFnur; 
  outr[3] += -1.0*incr1[3]*rdxFnur; 
  outr[4] += -1.0*incr1[4]*rdxFnur; 
  outr[5] += -1.0*incr1[5]*rdxFnur; 
  outr[6] += (-1.0*incr2[6]*rdxFnur)-1.0*incr1[6]*rdxFnur; 
  outr[7] += -1.0*incr1[7]*rdxFnur; 
  outr[8] += (-1.0*incr2[8]*rdxFnur)-1.0*incr1[8]*rdxFnur; 
  outr[9] += -1.0*incr1[9]*rdxFnur; 
  outr[10] += (-1.0*incr2[10]*rdxFnur)-1.0*incr1[10]*rdxFnur; 
  outr[11] += -1.0*incr1[11]*rdxFnur; 
  outr[12] += -1.0*incr1[12]*rdxFnur; 
  outr[13] += (-1.0*incr2[13]*rdxFnur)-1.0*incr1[13]*rdxFnur; 
  outr[14] += -1.0*incr1[14]*rdxFnur; 
  outr[15] += -1.0*incr1[15]*rdxFnur; 
  outr[16] += (-1.0*incr2[16]*rdxFnur)-1.0*incr1[16]*rdxFnur; 
  outr[17] += (-1.0*incr2[17]*rdxFnur)-1.0*incr1[17]*rdxFnur; 
  outr[18] += -1.0*incr1[18]*rdxFnur; 
  outr[19] += (-1.0*incr2[19]*rdxFnur)-1.0*incr1[19]*rdxFnur; 
  outr[20] += (-1.0*incr2[20]*rdxFnur)-1.0*incr1[20]*rdxFnur; 
  outr[21] += -1.0*incr1[21]*rdxFnur; 
  outr[22] += (-1.0*incr2[22]*rdxFnur)-1.0*incr1[22]*rdxFnur; 
  outr[23] += -1.0*incr1[23]*rdxFnur; 
  outr[24] += (-1.0*incr2[24]*rdxFnur)-1.0*incr1[24]*rdxFnur; 
  outr[25] += -1.0*incr1[25]*rdxFnur; 
  outr[26] += (-1.0*incr2[26]*rdxFnur)-1.0*incr1[26]*rdxFnur; 
  outr[27] += (-1.0*incr2[27]*rdxFnur)-1.0*incr1[27]*rdxFnur; 
  outr[28] += (-1.0*incr2[28]*rdxFnur)-1.0*incr1[28]*rdxFnur; 
  outr[29] += -1.0*incr1[29]*rdxFnur; 
  outr[30] += (-1.0*incr2[30]*rdxFnur)-1.0*incr1[30]*rdxFnur; 
  outr[31] += (-1.0*incr2[31]*rdxFnur)-1.0*incr1[31]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul; 
  outl[2] += incr2[2]*rdxFnul-1.0*incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul-1.0*incr1[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul-1.0*incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul-1.0*incr1[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul-1.0*incr1[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul-1.0*incr1[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul-1.0*incr1[17]*rdxFnul; 
  outl[18] += incr1[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul-1.0*incr1[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul-1.0*incr1[20]*rdxFnul; 
  outl[21] += incr1[21]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul-1.0*incr1[22]*rdxFnul; 
  outl[23] += incr1[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul-1.0*incr1[24]*rdxFnul; 
  outl[25] += incr1[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul-1.0*incr1[26]*rdxFnul; 
  outl[27] += incr2[27]*rdxFnul-1.0*incr1[27]*rdxFnul; 
  outl[28] += incr2[28]*rdxFnul-1.0*incr1[28]*rdxFnul; 
  outl[29] += incr1[29]*rdxFnul; 
  outl[30] += incr2[30]*rdxFnul-1.0*incr1[30]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul-1.0*incr1[31]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf5xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[2]/(dxl[2]*dxl[2]*dxl[2]*dxl[2]); 
  double rdxFnur = 16.0*nu[2]/(dxr[2]*dxr[2]*dxr[2]*dxr[2]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  incr1[0] = (-1.623797632095822*fr[3])-1.623797632095822*fl[3]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = (-1.623797632095822*fr[7])-1.623797632095822*fl[7]+0.9375*fr[1]-0.9375*fl[1]; 
  incr1[2] = (-1.623797632095822*fr[8])-1.623797632095822*fl[8]+0.9375*fr[2]-0.9375*fl[2]; 
  incr1[3] = 2.8125*fr[3]+2.8125*fl[3]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[4] = (-1.623797632095822*fr[11])-1.623797632095822*fl[11]+0.9375*fr[4]-0.9375*fl[4]; 
  incr1[5] = (-1.623797632095822*fr[14])-1.623797632095822*fl[14]+0.9375*fr[5]-0.9375*fl[5]; 
  incr1[6] = (-1.623797632095822*fr[16])-1.623797632095822*fl[16]+0.9375*fr[6]-0.9375*fl[6]; 
  incr1[7] = 2.8125*fr[7]+2.8125*fl[7]-1.623797632095822*fr[1]+1.623797632095822*fl[1]; 
  incr1[8] = 2.8125*fr[8]+2.8125*fl[8]-1.623797632095822*fr[2]+1.623797632095822*fl[2]; 
  incr1[9] = (-1.623797632095822*fr[18])-1.623797632095822*fl[18]+0.9375*fr[9]-0.9375*fl[9]; 
  incr1[10] = (-1.623797632095822*fr[19])-1.623797632095822*fl[19]+0.9375*fr[10]-0.9375*fl[10]; 
  incr1[11] = 2.8125*fr[11]+2.8125*fl[11]-1.623797632095822*fr[4]+1.623797632095822*fl[4]; 
  incr1[12] = (-1.623797632095822*fr[21])-1.623797632095822*fl[21]+0.9375*fr[12]-0.9375*fl[12]; 
  incr1[13] = (-1.623797632095822*fr[22])-1.623797632095822*fl[22]+0.9375*fr[13]-0.9375*fl[13]; 
  incr1[14] = 2.8125*fr[14]+2.8125*fl[14]-1.623797632095822*fr[5]+1.623797632095822*fl[5]; 
  incr1[15] = (-1.623797632095822*fr[25])-1.623797632095822*fl[25]+0.9375*fr[15]-0.9375*fl[15]; 
  incr1[16] = 2.8125*fr[16]+2.8125*fl[16]-1.623797632095822*fr[6]+1.623797632095822*fl[6]; 
  incr1[17] = (-1.623797632095822*fr[26])-1.623797632095822*fl[26]+0.9375*fr[17]-0.9375*fl[17]; 
  incr1[18] = 2.8125*fr[18]+2.8125*fl[18]-1.623797632095822*fr[9]+1.623797632095822*fl[9]; 
  incr1[19] = 2.8125*fr[19]+2.8125*fl[19]-1.623797632095822*fr[10]+1.623797632095822*fl[10]; 
  incr1[20] = (-1.623797632095822*fr[27])-1.623797632095822*fl[27]+0.9375*fr[20]-0.9375*fl[20]; 
  incr1[21] = 2.8125*fr[21]+2.8125*fl[21]-1.623797632095822*fr[12]+1.623797632095822*fl[12]; 
  incr1[22] = 2.8125*fr[22]+2.8125*fl[22]-1.623797632095822*fr[13]+1.623797632095822*fl[13]; 
  incr1[23] = (-1.623797632095822*fr[29])-1.623797632095822*fl[29]+0.9375*fr[23]-0.9375*fl[23]; 
  incr1[24] = (-1.623797632095822*fr[30])-1.623797632095822*fl[30]+0.9375*fr[24]-0.9375*fl[24]; 
  incr1[25] = 2.8125*fr[25]+2.8125*fl[25]-1.623797632095822*fr[15]+1.623797632095822*fl[15]; 
  incr1[26] = 2.8125*fr[26]+2.8125*fl[26]-1.623797632095822*fr[17]+1.623797632095822*fl[17]; 
  incr1[27] = 2.8125*fr[27]+2.8125*fl[27]-1.623797632095822*fr[20]+1.623797632095822*fl[20]; 
  incr1[28] = (-1.623797632095822*fr[31])-1.623797632095822*fl[31]+0.9375*fr[28]-0.9375*fl[28]; 
  incr1[29] = 2.8125*fr[29]+2.8125*fl[29]-1.623797632095822*fr[23]+1.623797632095822*fl[23]; 
  incr1[30] = 2.8125*fr[30]+2.8125*fl[30]-1.623797632095822*fr[24]+1.623797632095822*fl[24]; 
  incr1[31] = 2.8125*fr[31]+2.8125*fl[31]-1.623797632095822*fr[28]+1.623797632095822*fl[28]; 

  incr2[3] = 0.75*fr[3]-0.75*fl[3]; 
  incr2[7] = 0.75*fr[7]-0.75*fl[7]; 
  incr2[8] = 0.75*fr[8]-0.75*fl[8]; 
  incr2[11] = 0.75*fr[11]-0.75*fl[11]; 
  incr2[14] = 0.75*fr[14]-0.75*fl[14]; 
  incr2[16] = 0.75*fr[16]-0.75*fl[16]; 
  incr2[18] = 0.75*fr[18]-0.75*fl[18]; 
  incr2[19] = 0.75*fr[19]-0.75*fl[19]; 
  incr2[21] = 0.75*fr[21]-0.75*fl[21]; 
  incr2[22] = 0.75*fr[22]-0.75*fl[22]; 
  incr2[25] = 0.75*fr[25]-0.75*fl[25]; 
  incr2[26] = 0.75*fr[26]-0.75*fl[26]; 
  incr2[27] = 0.75*fr[27]-0.75*fl[27]; 
  incr2[29] = 0.75*fr[29]-0.75*fl[29]; 
  incr2[30] = 0.75*fr[30]-0.75*fl[30]; 
  incr2[31] = 0.75*fr[31]-0.75*fl[31]; 



  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += -1.0*incr1[1]*rdxFnur; 
  outr[2] += -1.0*incr1[2]*rdxFnur; 
  outr[3] += (-1.0*incr2[3]*rdxFnur)-1.0*incr1[3]*rdxFnur; 
  outr[4] += -1.0*incr1[4]*rdxFnur; 
  outr[5] += -1.0*incr1[5]*rdxFnur; 
  outr[6] += -1.0*incr1[6]*rdxFnur; 
  outr[7] += (-1.0*incr2[7]*rdxFnur)-1.0*incr1[7]*rdxFnur; 
  outr[8] += (-1.0*incr2[8]*rdxFnur)-1.0*incr1[8]*rdxFnur; 
  outr[9] += -1.0*incr1[9]*rdxFnur; 
  outr[10] += -1.0*incr1[10]*rdxFnur; 
  outr[11] += (-1.0*incr2[11]*rdxFnur)-1.0*incr1[11]*rdxFnur; 
  outr[12] += -1.0*incr1[12]*rdxFnur; 
  outr[13] += -1.0*incr1[13]*rdxFnur; 
  outr[14] += (-1.0*incr2[14]*rdxFnur)-1.0*incr1[14]*rdxFnur; 
  outr[15] += -1.0*incr1[15]*rdxFnur; 
  outr[16] += (-1.0*incr2[16]*rdxFnur)-1.0*incr1[16]*rdxFnur; 
  outr[17] += -1.0*incr1[17]*rdxFnur; 
  outr[18] += (-1.0*incr2[18]*rdxFnur)-1.0*incr1[18]*rdxFnur; 
  outr[19] += (-1.0*incr2[19]*rdxFnur)-1.0*incr1[19]*rdxFnur; 
  outr[20] += -1.0*incr1[20]*rdxFnur; 
  outr[21] += (-1.0*incr2[21]*rdxFnur)-1.0*incr1[21]*rdxFnur; 
  outr[22] += (-1.0*incr2[22]*rdxFnur)-1.0*incr1[22]*rdxFnur; 
  outr[23] += -1.0*incr1[23]*rdxFnur; 
  outr[24] += -1.0*incr1[24]*rdxFnur; 
  outr[25] += (-1.0*incr2[25]*rdxFnur)-1.0*incr1[25]*rdxFnur; 
  outr[26] += (-1.0*incr2[26]*rdxFnur)-1.0*incr1[26]*rdxFnur; 
  outr[27] += (-1.0*incr2[27]*rdxFnur)-1.0*incr1[27]*rdxFnur; 
  outr[28] += -1.0*incr1[28]*rdxFnur; 
  outr[29] += (-1.0*incr2[29]*rdxFnur)-1.0*incr1[29]*rdxFnur; 
  outr[30] += (-1.0*incr2[30]*rdxFnur)-1.0*incr1[30]*rdxFnur; 
  outr[31] += (-1.0*incr2[31]*rdxFnur)-1.0*incr1[31]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul; 
  outl[3] += incr2[3]*rdxFnul-1.0*incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul-1.0*incr1[7]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul-1.0*incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul-1.0*incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul-1.0*incr1[14]*rdxFnul; 
  outl[15] += incr1[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul-1.0*incr1[16]*rdxFnul; 
  outl[17] += incr1[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul-1.0*incr1[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul-1.0*incr1[19]*rdxFnul; 
  outl[20] += incr1[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul-1.0*incr1[21]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul-1.0*incr1[22]*rdxFnul; 
  outl[23] += incr1[23]*rdxFnul; 
  outl[24] += incr1[24]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul-1.0*incr1[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul-1.0*incr1[26]*rdxFnul; 
  outl[27] += incr2[27]*rdxFnul-1.0*incr1[27]*rdxFnul; 
  outl[28] += incr1[28]*rdxFnul; 
  outl[29] += incr2[29]*rdxFnul-1.0*incr1[29]*rdxFnul; 
  outl[30] += incr2[30]*rdxFnul-1.0*incr1[30]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul-1.0*incr1[31]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf5xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[3]/(dxl[3]*dxl[3]*dxl[3]*dxl[3]); 
  double rdxFnur = 16.0*nu[3]/(dxr[3]*dxr[3]*dxr[3]*dxr[3]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  incr1[0] = (-1.623797632095822*fr[4])-1.623797632095822*fl[4]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = (-1.623797632095822*fr[9])-1.623797632095822*fl[9]+0.9375*fr[1]-0.9375*fl[1]; 
  incr1[2] = (-1.623797632095822*fr[10])-1.623797632095822*fl[10]+0.9375*fr[2]-0.9375*fl[2]; 
  incr1[3] = (-1.623797632095822*fr[11])-1.623797632095822*fl[11]+0.9375*fr[3]-0.9375*fl[3]; 
  incr1[4] = 2.8125*fr[4]+2.8125*fl[4]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[5] = (-1.623797632095822*fr[15])-1.623797632095822*fl[15]+0.9375*fr[5]-0.9375*fl[5]; 
  incr1[6] = (-1.623797632095822*fr[17])-1.623797632095822*fl[17]+0.9375*fr[6]-0.9375*fl[6]; 
  incr1[7] = (-1.623797632095822*fr[18])-1.623797632095822*fl[18]+0.9375*fr[7]-0.9375*fl[7]; 
  incr1[8] = (-1.623797632095822*fr[19])-1.623797632095822*fl[19]+0.9375*fr[8]-0.9375*fl[8]; 
  incr1[9] = 2.8125*fr[9]+2.8125*fl[9]-1.623797632095822*fr[1]+1.623797632095822*fl[1]; 
  incr1[10] = 2.8125*fr[10]+2.8125*fl[10]-1.623797632095822*fr[2]+1.623797632095822*fl[2]; 
  incr1[11] = 2.8125*fr[11]+2.8125*fl[11]-1.623797632095822*fr[3]+1.623797632095822*fl[3]; 
  incr1[12] = (-1.623797632095822*fr[23])-1.623797632095822*fl[23]+0.9375*fr[12]-0.9375*fl[12]; 
  incr1[13] = (-1.623797632095822*fr[24])-1.623797632095822*fl[24]+0.9375*fr[13]-0.9375*fl[13]; 
  incr1[14] = (-1.623797632095822*fr[25])-1.623797632095822*fl[25]+0.9375*fr[14]-0.9375*fl[14]; 
  incr1[15] = 2.8125*fr[15]+2.8125*fl[15]-1.623797632095822*fr[5]+1.623797632095822*fl[5]; 
  incr1[16] = (-1.623797632095822*fr[26])-1.623797632095822*fl[26]+0.9375*fr[16]-0.9375*fl[16]; 
  incr1[17] = 2.8125*fr[17]+2.8125*fl[17]-1.623797632095822*fr[6]+1.623797632095822*fl[6]; 
  incr1[18] = 2.8125*fr[18]+2.8125*fl[18]-1.623797632095822*fr[7]+1.623797632095822*fl[7]; 
  incr1[19] = 2.8125*fr[19]+2.8125*fl[19]-1.623797632095822*fr[8]+1.623797632095822*fl[8]; 
  incr1[20] = (-1.623797632095822*fr[28])-1.623797632095822*fl[28]+0.9375*fr[20]-0.9375*fl[20]; 
  incr1[21] = (-1.623797632095822*fr[29])-1.623797632095822*fl[29]+0.9375*fr[21]-0.9375*fl[21]; 
  incr1[22] = (-1.623797632095822*fr[30])-1.623797632095822*fl[30]+0.9375*fr[22]-0.9375*fl[22]; 
  incr1[23] = 2.8125*fr[23]+2.8125*fl[23]-1.623797632095822*fr[12]+1.623797632095822*fl[12]; 
  incr1[24] = 2.8125*fr[24]+2.8125*fl[24]-1.623797632095822*fr[13]+1.623797632095822*fl[13]; 
  incr1[25] = 2.8125*fr[25]+2.8125*fl[25]-1.623797632095822*fr[14]+1.623797632095822*fl[14]; 
  incr1[26] = 2.8125*fr[26]+2.8125*fl[26]-1.623797632095822*fr[16]+1.623797632095822*fl[16]; 
  incr1[27] = (-1.623797632095822*fr[31])-1.623797632095822*fl[31]+0.9375*fr[27]-0.9375*fl[27]; 
  incr1[28] = 2.8125*fr[28]+2.8125*fl[28]-1.623797632095822*fr[20]+1.623797632095822*fl[20]; 
  incr1[29] = 2.8125*fr[29]+2.8125*fl[29]-1.623797632095822*fr[21]+1.623797632095822*fl[21]; 
  incr1[30] = 2.8125*fr[30]+2.8125*fl[30]-1.623797632095822*fr[22]+1.623797632095822*fl[22]; 
  incr1[31] = 2.8125*fr[31]+2.8125*fl[31]-1.623797632095822*fr[27]+1.623797632095822*fl[27]; 

  incr2[4] = 0.75*fr[4]-0.75*fl[4]; 
  incr2[9] = 0.75*fr[9]-0.75*fl[9]; 
  incr2[10] = 0.75*fr[10]-0.75*fl[10]; 
  incr2[11] = 0.75*fr[11]-0.75*fl[11]; 
  incr2[15] = 0.75*fr[15]-0.75*fl[15]; 
  incr2[17] = 0.75*fr[17]-0.75*fl[17]; 
  incr2[18] = 0.75*fr[18]-0.75*fl[18]; 
  incr2[19] = 0.75*fr[19]-0.75*fl[19]; 
  incr2[23] = 0.75*fr[23]-0.75*fl[23]; 
  incr2[24] = 0.75*fr[24]-0.75*fl[24]; 
  incr2[25] = 0.75*fr[25]-0.75*fl[25]; 
  incr2[26] = 0.75*fr[26]-0.75*fl[26]; 
  incr2[28] = 0.75*fr[28]-0.75*fl[28]; 
  incr2[29] = 0.75*fr[29]-0.75*fl[29]; 
  incr2[30] = 0.75*fr[30]-0.75*fl[30]; 
  incr2[31] = 0.75*fr[31]-0.75*fl[31]; 



  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += -1.0*incr1[1]*rdxFnur; 
  outr[2] += -1.0*incr1[2]*rdxFnur; 
  outr[3] += -1.0*incr1[3]*rdxFnur; 
  outr[4] += (-1.0*incr2[4]*rdxFnur)-1.0*incr1[4]*rdxFnur; 
  outr[5] += -1.0*incr1[5]*rdxFnur; 
  outr[6] += -1.0*incr1[6]*rdxFnur; 
  outr[7] += -1.0*incr1[7]*rdxFnur; 
  outr[8] += -1.0*incr1[8]*rdxFnur; 
  outr[9] += (-1.0*incr2[9]*rdxFnur)-1.0*incr1[9]*rdxFnur; 
  outr[10] += (-1.0*incr2[10]*rdxFnur)-1.0*incr1[10]*rdxFnur; 
  outr[11] += (-1.0*incr2[11]*rdxFnur)-1.0*incr1[11]*rdxFnur; 
  outr[12] += -1.0*incr1[12]*rdxFnur; 
  outr[13] += -1.0*incr1[13]*rdxFnur; 
  outr[14] += -1.0*incr1[14]*rdxFnur; 
  outr[15] += (-1.0*incr2[15]*rdxFnur)-1.0*incr1[15]*rdxFnur; 
  outr[16] += -1.0*incr1[16]*rdxFnur; 
  outr[17] += (-1.0*incr2[17]*rdxFnur)-1.0*incr1[17]*rdxFnur; 
  outr[18] += (-1.0*incr2[18]*rdxFnur)-1.0*incr1[18]*rdxFnur; 
  outr[19] += (-1.0*incr2[19]*rdxFnur)-1.0*incr1[19]*rdxFnur; 
  outr[20] += -1.0*incr1[20]*rdxFnur; 
  outr[21] += -1.0*incr1[21]*rdxFnur; 
  outr[22] += -1.0*incr1[22]*rdxFnur; 
  outr[23] += (-1.0*incr2[23]*rdxFnur)-1.0*incr1[23]*rdxFnur; 
  outr[24] += (-1.0*incr2[24]*rdxFnur)-1.0*incr1[24]*rdxFnur; 
  outr[25] += (-1.0*incr2[25]*rdxFnur)-1.0*incr1[25]*rdxFnur; 
  outr[26] += (-1.0*incr2[26]*rdxFnur)-1.0*incr1[26]*rdxFnur; 
  outr[27] += -1.0*incr1[27]*rdxFnur; 
  outr[28] += (-1.0*incr2[28]*rdxFnur)-1.0*incr1[28]*rdxFnur; 
  outr[29] += (-1.0*incr2[29]*rdxFnur)-1.0*incr1[29]*rdxFnur; 
  outr[30] += (-1.0*incr2[30]*rdxFnur)-1.0*incr1[30]*rdxFnur; 
  outr[31] += (-1.0*incr2[31]*rdxFnur)-1.0*incr1[31]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul-1.0*incr1[4]*rdxFnul; 
  outl[5] += incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul-1.0*incr1[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul-1.0*incr1[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul-1.0*incr1[11]*rdxFnul; 
  outl[12] += incr1[12]*rdxFnul; 
  outl[13] += incr1[13]*rdxFnul; 
  outl[14] += incr1[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul-1.0*incr1[15]*rdxFnul; 
  outl[16] += incr1[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul-1.0*incr1[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul-1.0*incr1[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul-1.0*incr1[19]*rdxFnul; 
  outl[20] += incr1[20]*rdxFnul; 
  outl[21] += incr1[21]*rdxFnul; 
  outl[22] += incr1[22]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul-1.0*incr1[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul-1.0*incr1[24]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul-1.0*incr1[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul-1.0*incr1[26]*rdxFnul; 
  outl[27] += incr1[27]*rdxFnul; 
  outl[28] += incr2[28]*rdxFnul-1.0*incr1[28]*rdxFnul; 
  outl[29] += incr2[29]*rdxFnul-1.0*incr1[29]*rdxFnul; 
  outl[30] += incr2[30]*rdxFnul-1.0*incr1[30]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul-1.0*incr1[31]*rdxFnul; 

} 
void ConstHyperDiffusion4Surf5xSerP1_X5(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[5]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 16.0*nu[4]/(dxl[4]*dxl[4]*dxl[4]*dxl[4]); 
  double rdxFnur = 16.0*nu[4]/(dxr[4]*dxr[4]*dxr[4]*dxr[4]); 

  double incr1[32]; 
  double incr2[32]; 
  double incr3[32]; 
  double incr4[32]; 

  incr1[0] = (-1.623797632095822*fr[5])-1.623797632095822*fl[5]+0.9375*fr[0]-0.9375*fl[0]; 
  incr1[1] = (-1.623797632095822*fr[12])-1.623797632095822*fl[12]+0.9375*fr[1]-0.9375*fl[1]; 
  incr1[2] = (-1.623797632095822*fr[13])-1.623797632095822*fl[13]+0.9375*fr[2]-0.9375*fl[2]; 
  incr1[3] = (-1.623797632095822*fr[14])-1.623797632095822*fl[14]+0.9375*fr[3]-0.9375*fl[3]; 
  incr1[4] = (-1.623797632095822*fr[15])-1.623797632095822*fl[15]+0.9375*fr[4]-0.9375*fl[4]; 
  incr1[5] = 2.8125*fr[5]+2.8125*fl[5]-1.623797632095822*fr[0]+1.623797632095822*fl[0]; 
  incr1[6] = (-1.623797632095822*fr[20])-1.623797632095822*fl[20]+0.9375*fr[6]-0.9375*fl[6]; 
  incr1[7] = (-1.623797632095822*fr[21])-1.623797632095822*fl[21]+0.9375*fr[7]-0.9375*fl[7]; 
  incr1[8] = (-1.623797632095822*fr[22])-1.623797632095822*fl[22]+0.9375*fr[8]-0.9375*fl[8]; 
  incr1[9] = (-1.623797632095822*fr[23])-1.623797632095822*fl[23]+0.9375*fr[9]-0.9375*fl[9]; 
  incr1[10] = (-1.623797632095822*fr[24])-1.623797632095822*fl[24]+0.9375*fr[10]-0.9375*fl[10]; 
  incr1[11] = (-1.623797632095822*fr[25])-1.623797632095822*fl[25]+0.9375*fr[11]-0.9375*fl[11]; 
  incr1[12] = 2.8125*fr[12]+2.8125*fl[12]-1.623797632095822*fr[1]+1.623797632095822*fl[1]; 
  incr1[13] = 2.8125*fr[13]+2.8125*fl[13]-1.623797632095822*fr[2]+1.623797632095822*fl[2]; 
  incr1[14] = 2.8125*fr[14]+2.8125*fl[14]-1.623797632095822*fr[3]+1.623797632095822*fl[3]; 
  incr1[15] = 2.8125*fr[15]+2.8125*fl[15]-1.623797632095822*fr[4]+1.623797632095822*fl[4]; 
  incr1[16] = (-1.623797632095822*fr[27])-1.623797632095822*fl[27]+0.9375*fr[16]-0.9375*fl[16]; 
  incr1[17] = (-1.623797632095822*fr[28])-1.623797632095822*fl[28]+0.9375*fr[17]-0.9375*fl[17]; 
  incr1[18] = (-1.623797632095822*fr[29])-1.623797632095822*fl[29]+0.9375*fr[18]-0.9375*fl[18]; 
  incr1[19] = (-1.623797632095822*fr[30])-1.623797632095822*fl[30]+0.9375*fr[19]-0.9375*fl[19]; 
  incr1[20] = 2.8125*fr[20]+2.8125*fl[20]-1.623797632095822*fr[6]+1.623797632095822*fl[6]; 
  incr1[21] = 2.8125*fr[21]+2.8125*fl[21]-1.623797632095822*fr[7]+1.623797632095822*fl[7]; 
  incr1[22] = 2.8125*fr[22]+2.8125*fl[22]-1.623797632095822*fr[8]+1.623797632095822*fl[8]; 
  incr1[23] = 2.8125*fr[23]+2.8125*fl[23]-1.623797632095822*fr[9]+1.623797632095822*fl[9]; 
  incr1[24] = 2.8125*fr[24]+2.8125*fl[24]-1.623797632095822*fr[10]+1.623797632095822*fl[10]; 
  incr1[25] = 2.8125*fr[25]+2.8125*fl[25]-1.623797632095822*fr[11]+1.623797632095822*fl[11]; 
  incr1[26] = (-1.623797632095822*fr[31])-1.623797632095822*fl[31]+0.9375*fr[26]-0.9375*fl[26]; 
  incr1[27] = 2.8125*fr[27]+2.8125*fl[27]-1.623797632095822*fr[16]+1.623797632095822*fl[16]; 
  incr1[28] = 2.8125*fr[28]+2.8125*fl[28]-1.623797632095822*fr[17]+1.623797632095822*fl[17]; 
  incr1[29] = 2.8125*fr[29]+2.8125*fl[29]-1.623797632095822*fr[18]+1.623797632095822*fl[18]; 
  incr1[30] = 2.8125*fr[30]+2.8125*fl[30]-1.623797632095822*fr[19]+1.623797632095822*fl[19]; 
  incr1[31] = 2.8125*fr[31]+2.8125*fl[31]-1.623797632095822*fr[26]+1.623797632095822*fl[26]; 

  incr2[5] = 0.75*fr[5]-0.75*fl[5]; 
  incr2[12] = 0.75*fr[12]-0.75*fl[12]; 
  incr2[13] = 0.75*fr[13]-0.75*fl[13]; 
  incr2[14] = 0.75*fr[14]-0.75*fl[14]; 
  incr2[15] = 0.75*fr[15]-0.75*fl[15]; 
  incr2[20] = 0.75*fr[20]-0.75*fl[20]; 
  incr2[21] = 0.75*fr[21]-0.75*fl[21]; 
  incr2[22] = 0.75*fr[22]-0.75*fl[22]; 
  incr2[23] = 0.75*fr[23]-0.75*fl[23]; 
  incr2[24] = 0.75*fr[24]-0.75*fl[24]; 
  incr2[25] = 0.75*fr[25]-0.75*fl[25]; 
  incr2[27] = 0.75*fr[27]-0.75*fl[27]; 
  incr2[28] = 0.75*fr[28]-0.75*fl[28]; 
  incr2[29] = 0.75*fr[29]-0.75*fl[29]; 
  incr2[30] = 0.75*fr[30]-0.75*fl[30]; 
  incr2[31] = 0.75*fr[31]-0.75*fl[31]; 



  outr[0] += -1.0*incr1[0]*rdxFnur; 
  outr[1] += -1.0*incr1[1]*rdxFnur; 
  outr[2] += -1.0*incr1[2]*rdxFnur; 
  outr[3] += -1.0*incr1[3]*rdxFnur; 
  outr[4] += -1.0*incr1[4]*rdxFnur; 
  outr[5] += (-1.0*incr2[5]*rdxFnur)-1.0*incr1[5]*rdxFnur; 
  outr[6] += -1.0*incr1[6]*rdxFnur; 
  outr[7] += -1.0*incr1[7]*rdxFnur; 
  outr[8] += -1.0*incr1[8]*rdxFnur; 
  outr[9] += -1.0*incr1[9]*rdxFnur; 
  outr[10] += -1.0*incr1[10]*rdxFnur; 
  outr[11] += -1.0*incr1[11]*rdxFnur; 
  outr[12] += (-1.0*incr2[12]*rdxFnur)-1.0*incr1[12]*rdxFnur; 
  outr[13] += (-1.0*incr2[13]*rdxFnur)-1.0*incr1[13]*rdxFnur; 
  outr[14] += (-1.0*incr2[14]*rdxFnur)-1.0*incr1[14]*rdxFnur; 
  outr[15] += (-1.0*incr2[15]*rdxFnur)-1.0*incr1[15]*rdxFnur; 
  outr[16] += -1.0*incr1[16]*rdxFnur; 
  outr[17] += -1.0*incr1[17]*rdxFnur; 
  outr[18] += -1.0*incr1[18]*rdxFnur; 
  outr[19] += -1.0*incr1[19]*rdxFnur; 
  outr[20] += (-1.0*incr2[20]*rdxFnur)-1.0*incr1[20]*rdxFnur; 
  outr[21] += (-1.0*incr2[21]*rdxFnur)-1.0*incr1[21]*rdxFnur; 
  outr[22] += (-1.0*incr2[22]*rdxFnur)-1.0*incr1[22]*rdxFnur; 
  outr[23] += (-1.0*incr2[23]*rdxFnur)-1.0*incr1[23]*rdxFnur; 
  outr[24] += (-1.0*incr2[24]*rdxFnur)-1.0*incr1[24]*rdxFnur; 
  outr[25] += (-1.0*incr2[25]*rdxFnur)-1.0*incr1[25]*rdxFnur; 
  outr[26] += -1.0*incr1[26]*rdxFnur; 
  outr[27] += (-1.0*incr2[27]*rdxFnur)-1.0*incr1[27]*rdxFnur; 
  outr[28] += (-1.0*incr2[28]*rdxFnur)-1.0*incr1[28]*rdxFnur; 
  outr[29] += (-1.0*incr2[29]*rdxFnur)-1.0*incr1[29]*rdxFnur; 
  outr[30] += (-1.0*incr2[30]*rdxFnur)-1.0*incr1[30]*rdxFnur; 
  outr[31] += (-1.0*incr2[31]*rdxFnur)-1.0*incr1[31]*rdxFnur; 

  outl[0] += incr1[0]*rdxFnul; 
  outl[1] += incr1[1]*rdxFnul; 
  outl[2] += incr1[2]*rdxFnul; 
  outl[3] += incr1[3]*rdxFnul; 
  outl[4] += incr1[4]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul-1.0*incr1[5]*rdxFnul; 
  outl[6] += incr1[6]*rdxFnul; 
  outl[7] += incr1[7]*rdxFnul; 
  outl[8] += incr1[8]*rdxFnul; 
  outl[9] += incr1[9]*rdxFnul; 
  outl[10] += incr1[10]*rdxFnul; 
  outl[11] += incr1[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul-1.0*incr1[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul-1.0*incr1[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul-1.0*incr1[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul-1.0*incr1[15]*rdxFnul; 
  outl[16] += incr1[16]*rdxFnul; 
  outl[17] += incr1[17]*rdxFnul; 
  outl[18] += incr1[18]*rdxFnul; 
  outl[19] += incr1[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul-1.0*incr1[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul-1.0*incr1[21]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul-1.0*incr1[22]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul-1.0*incr1[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul-1.0*incr1[24]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul-1.0*incr1[25]*rdxFnul; 
  outl[26] += incr1[26]*rdxFnul; 
  outl[27] += incr2[27]*rdxFnul-1.0*incr1[27]*rdxFnul; 
  outl[28] += incr2[28]*rdxFnul-1.0*incr1[28]*rdxFnul; 
  outl[29] += incr2[29]*rdxFnul-1.0*incr1[29]*rdxFnul; 
  outl[30] += incr2[30]*rdxFnul-1.0*incr1[30]*rdxFnul; 
  outl[31] += incr2[31]*rdxFnul-1.0*incr1[31]*rdxFnul; 

} 
void ConstDiffusionVarCoeffSurf5xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[160]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[32]; 
  double incr2[32]; 

  incr1[0] = 0.1657281518405969*fr[7]*nul[7]+0.1657281518405969*fl[7]*nul[7]-0.172229747539442*fr[3]*nul[7]+0.172229747539442*fl[3]*nul[7]+0.09568319307746781*nul[3]*fr[7]+0.09568319307746781*nul[3]*fl[7]-0.09943689110435816*fr[3]*nul[3]+0.09943689110435816*fl[3]*nul[3]+0.1657281518405969*fr[1]*nul[1]+0.1657281518405969*fl[1]*nul[1]-0.172229747539442*fr[0]*nul[1]+0.172229747539442*fl[0]*nul[1]+0.09568319307746781*nul[0]*fr[1]+0.09568319307746781*nul[0]*fl[1]-0.09943689110435816*fr[0]*nul[0]+0.09943689110435816*fl[0]*nul[0]; 
  incr1[1] = (-0.2870495792324034*fr[7]*nul[7])-0.2870495792324034*fl[7]*nul[7]+0.2983106733130745*fr[3]*nul[7]-0.2983106733130745*fl[3]*nul[7]-0.1657281518405969*nul[3]*fr[7]-0.1657281518405969*nul[3]*fl[7]+0.172229747539442*fr[3]*nul[3]-0.172229747539442*fl[3]*nul[3]-0.2870495792324034*fr[1]*nul[1]-0.2870495792324034*fl[1]*nul[1]+0.2983106733130745*fr[0]*nul[1]-0.2983106733130745*fl[0]*nul[1]-0.1657281518405969*nul[0]*fr[1]-0.1657281518405969*nul[0]*fl[1]+0.172229747539442*fr[0]*nul[0]-0.172229747539442*fl[0]*nul[0]; 
  incr1[2] = 0.1657281518405969*nul[7]*fr[16]+0.09568319307746781*nul[3]*fr[16]+0.1657281518405969*nul[7]*fl[16]+0.09568319307746781*nul[3]*fl[16]-0.172229747539442*nul[7]*fr[8]-0.09943689110435816*nul[3]*fr[8]+0.172229747539442*nul[7]*fl[8]+0.09943689110435816*nul[3]*fl[8]+0.1657281518405969*nul[1]*fr[6]+0.09568319307746781*nul[0]*fr[6]+0.1657281518405969*nul[1]*fl[6]+0.09568319307746781*nul[0]*fl[6]-0.172229747539442*nul[1]*fr[2]-0.09943689110435816*nul[0]*fr[2]+0.172229747539442*nul[1]*fl[2]+0.09943689110435816*nul[0]*fl[2]; 
  incr1[3] = 0.1657281518405969*fr[1]*nul[7]+0.1657281518405969*fl[1]*nul[7]-0.172229747539442*fr[0]*nul[7]+0.172229747539442*fl[0]*nul[7]+0.1657281518405969*nul[1]*fr[7]+0.09568319307746781*nul[0]*fr[7]+0.1657281518405969*nul[1]*fl[7]+0.09568319307746781*nul[0]*fl[7]+0.09568319307746781*fr[1]*nul[3]+0.09568319307746781*fl[1]*nul[3]-0.09943689110435816*fr[0]*nul[3]+0.09943689110435816*fl[0]*nul[3]-0.172229747539442*nul[1]*fr[3]-0.09943689110435816*nul[0]*fr[3]+0.172229747539442*nul[1]*fl[3]+0.09943689110435816*nul[0]*fl[3]; 
  incr1[4] = 0.1657281518405969*nul[7]*fr[18]+0.09568319307746781*nul[3]*fr[18]+0.1657281518405969*nul[7]*fl[18]+0.09568319307746781*nul[3]*fl[18]-0.172229747539442*nul[7]*fr[11]-0.09943689110435816*nul[3]*fr[11]+0.172229747539442*nul[7]*fl[11]+0.09943689110435816*nul[3]*fl[11]+0.1657281518405969*nul[1]*fr[9]+0.09568319307746781*nul[0]*fr[9]+0.1657281518405969*nul[1]*fl[9]+0.09568319307746781*nul[0]*fl[9]-0.172229747539442*nul[1]*fr[4]-0.09943689110435816*nul[0]*fr[4]+0.172229747539442*nul[1]*fl[4]+0.09943689110435816*nul[0]*fl[4]; 
  incr1[5] = 0.1657281518405969*nul[7]*fr[21]+0.09568319307746781*nul[3]*fr[21]+0.1657281518405969*nul[7]*fl[21]+0.09568319307746781*nul[3]*fl[21]-0.172229747539442*nul[7]*fr[14]-0.09943689110435816*nul[3]*fr[14]+0.172229747539442*nul[7]*fl[14]+0.09943689110435816*nul[3]*fl[14]+0.1657281518405969*nul[1]*fr[12]+0.09568319307746781*nul[0]*fr[12]+0.1657281518405969*nul[1]*fl[12]+0.09568319307746781*nul[0]*fl[12]-0.172229747539442*nul[1]*fr[5]-0.09943689110435816*nul[0]*fr[5]+0.172229747539442*nul[1]*fl[5]+0.09943689110435816*nul[0]*fl[5]; 
  incr1[6] = (-0.2870495792324034*nul[7]*fr[16])-0.1657281518405969*nul[3]*fr[16]-0.2870495792324034*nul[7]*fl[16]-0.1657281518405969*nul[3]*fl[16]+0.2983106733130745*nul[7]*fr[8]+0.172229747539442*nul[3]*fr[8]-0.2983106733130745*nul[7]*fl[8]-0.172229747539442*nul[3]*fl[8]-0.2870495792324034*nul[1]*fr[6]-0.1657281518405969*nul[0]*fr[6]-0.2870495792324034*nul[1]*fl[6]-0.1657281518405969*nul[0]*fl[6]+0.2983106733130745*nul[1]*fr[2]+0.172229747539442*nul[0]*fr[2]-0.2983106733130745*nul[1]*fl[2]-0.172229747539442*nul[0]*fl[2]; 
  incr1[7] = (-0.2870495792324034*fr[1]*nul[7])-0.2870495792324034*fl[1]*nul[7]+0.2983106733130745*fr[0]*nul[7]-0.2983106733130745*fl[0]*nul[7]-0.2870495792324034*nul[1]*fr[7]-0.1657281518405969*nul[0]*fr[7]-0.2870495792324034*nul[1]*fl[7]-0.1657281518405969*nul[0]*fl[7]-0.1657281518405969*fr[1]*nul[3]-0.1657281518405969*fl[1]*nul[3]+0.172229747539442*fr[0]*nul[3]-0.172229747539442*fl[0]*nul[3]+0.2983106733130745*nul[1]*fr[3]+0.172229747539442*nul[0]*fr[3]-0.2983106733130745*nul[1]*fl[3]-0.172229747539442*nul[0]*fl[3]; 
  incr1[8] = 0.1657281518405969*nul[1]*fr[16]+0.09568319307746781*nul[0]*fr[16]+0.1657281518405969*nul[1]*fl[16]+0.09568319307746781*nul[0]*fl[16]-0.172229747539442*nul[1]*fr[8]-0.09943689110435816*nul[0]*fr[8]+0.172229747539442*nul[1]*fl[8]+0.09943689110435816*nul[0]*fl[8]+0.1657281518405969*fr[6]*nul[7]+0.1657281518405969*fl[6]*nul[7]-0.172229747539442*fr[2]*nul[7]+0.172229747539442*fl[2]*nul[7]+0.09568319307746781*nul[3]*fr[6]+0.09568319307746781*nul[3]*fl[6]-0.09943689110435816*fr[2]*nul[3]+0.09943689110435816*fl[2]*nul[3]; 
  incr1[9] = (-0.2870495792324034*nul[7]*fr[18])-0.1657281518405969*nul[3]*fr[18]-0.2870495792324034*nul[7]*fl[18]-0.1657281518405969*nul[3]*fl[18]+0.2983106733130745*nul[7]*fr[11]+0.172229747539442*nul[3]*fr[11]-0.2983106733130745*nul[7]*fl[11]-0.172229747539442*nul[3]*fl[11]-0.2870495792324034*nul[1]*fr[9]-0.1657281518405969*nul[0]*fr[9]-0.2870495792324034*nul[1]*fl[9]-0.1657281518405969*nul[0]*fl[9]+0.2983106733130745*nul[1]*fr[4]+0.172229747539442*nul[0]*fr[4]-0.2983106733130745*nul[1]*fl[4]-0.172229747539442*nul[0]*fl[4]; 
  incr1[10] = 0.1657281518405969*nul[7]*fr[26]+0.09568319307746781*nul[3]*fr[26]+0.1657281518405969*nul[7]*fl[26]+0.09568319307746781*nul[3]*fl[26]-0.172229747539442*nul[7]*fr[19]-0.09943689110435816*nul[3]*fr[19]+0.172229747539442*nul[7]*fl[19]+0.09943689110435816*nul[3]*fl[19]+0.1657281518405969*nul[1]*fr[17]+0.09568319307746781*nul[0]*fr[17]+0.1657281518405969*nul[1]*fl[17]+0.09568319307746781*nul[0]*fl[17]-0.172229747539442*nul[1]*fr[10]-0.09943689110435816*nul[0]*fr[10]+0.172229747539442*nul[1]*fl[10]+0.09943689110435816*nul[0]*fl[10]; 
  incr1[11] = 0.1657281518405969*nul[1]*fr[18]+0.09568319307746781*nul[0]*fr[18]+0.1657281518405969*nul[1]*fl[18]+0.09568319307746781*nul[0]*fl[18]-0.172229747539442*nul[1]*fr[11]-0.09943689110435816*nul[0]*fr[11]+0.172229747539442*nul[1]*fl[11]+0.09943689110435816*nul[0]*fl[11]+0.1657281518405969*nul[7]*fr[9]+0.09568319307746781*nul[3]*fr[9]+0.1657281518405969*nul[7]*fl[9]+0.09568319307746781*nul[3]*fl[9]-0.172229747539442*fr[4]*nul[7]+0.172229747539442*fl[4]*nul[7]-0.09943689110435816*nul[3]*fr[4]+0.09943689110435816*nul[3]*fl[4]; 
  incr1[12] = (-0.2870495792324034*nul[7]*fr[21])-0.1657281518405969*nul[3]*fr[21]-0.2870495792324034*nul[7]*fl[21]-0.1657281518405969*nul[3]*fl[21]+0.2983106733130745*nul[7]*fr[14]+0.172229747539442*nul[3]*fr[14]-0.2983106733130745*nul[7]*fl[14]-0.172229747539442*nul[3]*fl[14]-0.2870495792324034*nul[1]*fr[12]-0.1657281518405969*nul[0]*fr[12]-0.2870495792324034*nul[1]*fl[12]-0.1657281518405969*nul[0]*fl[12]+0.2983106733130745*nul[1]*fr[5]+0.172229747539442*nul[0]*fr[5]-0.2983106733130745*nul[1]*fl[5]-0.172229747539442*nul[0]*fl[5]; 
  incr1[13] = 0.1657281518405969*nul[7]*fr[27]+0.09568319307746781*nul[3]*fr[27]+0.1657281518405969*nul[7]*fl[27]+0.09568319307746781*nul[3]*fl[27]-0.172229747539442*nul[7]*fr[22]-0.09943689110435816*nul[3]*fr[22]+0.172229747539442*nul[7]*fl[22]+0.09943689110435816*nul[3]*fl[22]+0.1657281518405969*nul[1]*fr[20]+0.09568319307746781*nul[0]*fr[20]+0.1657281518405969*nul[1]*fl[20]+0.09568319307746781*nul[0]*fl[20]-0.172229747539442*nul[1]*fr[13]-0.09943689110435816*nul[0]*fr[13]+0.172229747539442*nul[1]*fl[13]+0.09943689110435816*nul[0]*fl[13]; 
  incr1[14] = 0.1657281518405969*nul[1]*fr[21]+0.09568319307746781*nul[0]*fr[21]+0.1657281518405969*nul[1]*fl[21]+0.09568319307746781*nul[0]*fl[21]-0.172229747539442*nul[1]*fr[14]-0.09943689110435816*nul[0]*fr[14]+0.172229747539442*nul[1]*fl[14]+0.09943689110435816*nul[0]*fl[14]+0.1657281518405969*nul[7]*fr[12]+0.09568319307746781*nul[3]*fr[12]+0.1657281518405969*nul[7]*fl[12]+0.09568319307746781*nul[3]*fl[12]-0.172229747539442*fr[5]*nul[7]+0.172229747539442*fl[5]*nul[7]-0.09943689110435816*nul[3]*fr[5]+0.09943689110435816*nul[3]*fl[5]; 
  incr1[15] = 0.1657281518405969*nul[7]*fr[29]+0.09568319307746781*nul[3]*fr[29]+0.1657281518405969*nul[7]*fl[29]+0.09568319307746781*nul[3]*fl[29]-0.172229747539442*nul[7]*fr[25]-0.09943689110435816*nul[3]*fr[25]+0.172229747539442*nul[7]*fl[25]+0.09943689110435816*nul[3]*fl[25]+0.1657281518405969*nul[1]*fr[23]+0.09568319307746781*nul[0]*fr[23]+0.1657281518405969*nul[1]*fl[23]+0.09568319307746781*nul[0]*fl[23]-0.172229747539442*nul[1]*fr[15]-0.09943689110435816*nul[0]*fr[15]+0.172229747539442*nul[1]*fl[15]+0.09943689110435816*nul[0]*fl[15]; 
  incr1[16] = (-0.2870495792324034*nul[1]*fr[16])-0.1657281518405969*nul[0]*fr[16]-0.2870495792324034*nul[1]*fl[16]-0.1657281518405969*nul[0]*fl[16]+0.2983106733130745*nul[1]*fr[8]+0.172229747539442*nul[0]*fr[8]-0.2983106733130745*nul[1]*fl[8]-0.172229747539442*nul[0]*fl[8]-0.2870495792324034*fr[6]*nul[7]-0.2870495792324034*fl[6]*nul[7]+0.2983106733130745*fr[2]*nul[7]-0.2983106733130745*fl[2]*nul[7]-0.1657281518405969*nul[3]*fr[6]-0.1657281518405969*nul[3]*fl[6]+0.172229747539442*fr[2]*nul[3]-0.172229747539442*fl[2]*nul[3]; 
  incr1[17] = (-0.2870495792324034*nul[7]*fr[26])-0.1657281518405969*nul[3]*fr[26]-0.2870495792324034*nul[7]*fl[26]-0.1657281518405969*nul[3]*fl[26]+0.2983106733130745*nul[7]*fr[19]+0.172229747539442*nul[3]*fr[19]-0.2983106733130745*nul[7]*fl[19]-0.172229747539442*nul[3]*fl[19]-0.2870495792324034*nul[1]*fr[17]-0.1657281518405969*nul[0]*fr[17]-0.2870495792324034*nul[1]*fl[17]-0.1657281518405969*nul[0]*fl[17]+0.2983106733130745*nul[1]*fr[10]+0.172229747539442*nul[0]*fr[10]-0.2983106733130745*nul[1]*fl[10]-0.172229747539442*nul[0]*fl[10]; 
  incr1[18] = (-0.2870495792324034*nul[1]*fr[18])-0.1657281518405969*nul[0]*fr[18]-0.2870495792324034*nul[1]*fl[18]-0.1657281518405969*nul[0]*fl[18]+0.2983106733130745*nul[1]*fr[11]+0.172229747539442*nul[0]*fr[11]-0.2983106733130745*nul[1]*fl[11]-0.172229747539442*nul[0]*fl[11]-0.2870495792324034*nul[7]*fr[9]-0.1657281518405969*nul[3]*fr[9]-0.2870495792324034*nul[7]*fl[9]-0.1657281518405969*nul[3]*fl[9]+0.2983106733130745*fr[4]*nul[7]-0.2983106733130745*fl[4]*nul[7]+0.172229747539442*nul[3]*fr[4]-0.172229747539442*nul[3]*fl[4]; 
  incr1[19] = 0.1657281518405969*nul[1]*fr[26]+0.09568319307746781*nul[0]*fr[26]+0.1657281518405969*nul[1]*fl[26]+0.09568319307746781*nul[0]*fl[26]-0.172229747539442*nul[1]*fr[19]-0.09943689110435816*nul[0]*fr[19]+0.172229747539442*nul[1]*fl[19]+0.09943689110435816*nul[0]*fl[19]+0.1657281518405969*nul[7]*fr[17]+0.09568319307746781*nul[3]*fr[17]+0.1657281518405969*nul[7]*fl[17]+0.09568319307746781*nul[3]*fl[17]-0.172229747539442*nul[7]*fr[10]-0.09943689110435816*nul[3]*fr[10]+0.172229747539442*nul[7]*fl[10]+0.09943689110435816*nul[3]*fl[10]; 
  incr1[20] = (-0.2870495792324034*nul[7]*fr[27])-0.1657281518405969*nul[3]*fr[27]-0.2870495792324034*nul[7]*fl[27]-0.1657281518405969*nul[3]*fl[27]+0.2983106733130745*nul[7]*fr[22]+0.172229747539442*nul[3]*fr[22]-0.2983106733130745*nul[7]*fl[22]-0.172229747539442*nul[3]*fl[22]-0.2870495792324034*nul[1]*fr[20]-0.1657281518405969*nul[0]*fr[20]-0.2870495792324034*nul[1]*fl[20]-0.1657281518405969*nul[0]*fl[20]+0.2983106733130745*nul[1]*fr[13]+0.172229747539442*nul[0]*fr[13]-0.2983106733130745*nul[1]*fl[13]-0.172229747539442*nul[0]*fl[13]; 
  incr1[21] = (-0.2870495792324034*nul[1]*fr[21])-0.1657281518405969*nul[0]*fr[21]-0.2870495792324034*nul[1]*fl[21]-0.1657281518405969*nul[0]*fl[21]+0.2983106733130745*nul[1]*fr[14]+0.172229747539442*nul[0]*fr[14]-0.2983106733130745*nul[1]*fl[14]-0.172229747539442*nul[0]*fl[14]-0.2870495792324034*nul[7]*fr[12]-0.1657281518405969*nul[3]*fr[12]-0.2870495792324034*nul[7]*fl[12]-0.1657281518405969*nul[3]*fl[12]+0.2983106733130745*fr[5]*nul[7]-0.2983106733130745*fl[5]*nul[7]+0.172229747539442*nul[3]*fr[5]-0.172229747539442*nul[3]*fl[5]; 
  incr1[22] = 0.1657281518405969*nul[1]*fr[27]+0.09568319307746781*nul[0]*fr[27]+0.1657281518405969*nul[1]*fl[27]+0.09568319307746781*nul[0]*fl[27]-0.172229747539442*nul[1]*fr[22]-0.09943689110435816*nul[0]*fr[22]+0.172229747539442*nul[1]*fl[22]+0.09943689110435816*nul[0]*fl[22]+0.1657281518405969*nul[7]*fr[20]+0.09568319307746781*nul[3]*fr[20]+0.1657281518405969*nul[7]*fl[20]+0.09568319307746781*nul[3]*fl[20]-0.172229747539442*nul[7]*fr[13]-0.09943689110435816*nul[3]*fr[13]+0.172229747539442*nul[7]*fl[13]+0.09943689110435816*nul[3]*fl[13]; 
  incr1[23] = (-0.2870495792324034*nul[7]*fr[29])-0.1657281518405969*nul[3]*fr[29]-0.2870495792324034*nul[7]*fl[29]-0.1657281518405969*nul[3]*fl[29]+0.2983106733130745*nul[7]*fr[25]+0.172229747539442*nul[3]*fr[25]-0.2983106733130745*nul[7]*fl[25]-0.172229747539442*nul[3]*fl[25]-0.2870495792324034*nul[1]*fr[23]-0.1657281518405969*nul[0]*fr[23]-0.2870495792324034*nul[1]*fl[23]-0.1657281518405969*nul[0]*fl[23]+0.2983106733130745*nul[1]*fr[15]+0.172229747539442*nul[0]*fr[15]-0.2983106733130745*nul[1]*fl[15]-0.172229747539442*nul[0]*fl[15]; 
  incr1[24] = 0.1657281518405969*nul[7]*fr[31]+0.09568319307746781*nul[3]*fr[31]+0.1657281518405969*nul[7]*fl[31]+0.09568319307746781*nul[3]*fl[31]-0.172229747539442*nul[7]*fr[30]-0.09943689110435816*nul[3]*fr[30]+0.172229747539442*nul[7]*fl[30]+0.09943689110435816*nul[3]*fl[30]+0.1657281518405969*nul[1]*fr[28]+0.09568319307746781*nul[0]*fr[28]+0.1657281518405969*nul[1]*fl[28]+0.09568319307746781*nul[0]*fl[28]-0.172229747539442*nul[1]*fr[24]-0.09943689110435816*nul[0]*fr[24]+0.172229747539442*nul[1]*fl[24]+0.09943689110435816*nul[0]*fl[24]; 
  incr1[25] = 0.1657281518405969*nul[1]*fr[29]+0.09568319307746781*nul[0]*fr[29]+0.1657281518405969*nul[1]*fl[29]+0.09568319307746781*nul[0]*fl[29]-0.172229747539442*nul[1]*fr[25]-0.09943689110435816*nul[0]*fr[25]+0.172229747539442*nul[1]*fl[25]+0.09943689110435816*nul[0]*fl[25]+0.1657281518405969*nul[7]*fr[23]+0.09568319307746781*nul[3]*fr[23]+0.1657281518405969*nul[7]*fl[23]+0.09568319307746781*nul[3]*fl[23]-0.172229747539442*nul[7]*fr[15]-0.09943689110435816*nul[3]*fr[15]+0.172229747539442*nul[7]*fl[15]+0.09943689110435816*nul[3]*fl[15]; 
  incr1[26] = (-0.2870495792324034*nul[1]*fr[26])-0.1657281518405969*nul[0]*fr[26]-0.2870495792324034*nul[1]*fl[26]-0.1657281518405969*nul[0]*fl[26]+0.2983106733130745*nul[1]*fr[19]+0.172229747539442*nul[0]*fr[19]-0.2983106733130745*nul[1]*fl[19]-0.172229747539442*nul[0]*fl[19]-0.2870495792324034*nul[7]*fr[17]-0.1657281518405969*nul[3]*fr[17]-0.2870495792324034*nul[7]*fl[17]-0.1657281518405969*nul[3]*fl[17]+0.2983106733130745*nul[7]*fr[10]+0.172229747539442*nul[3]*fr[10]-0.2983106733130745*nul[7]*fl[10]-0.172229747539442*nul[3]*fl[10]; 
  incr1[27] = (-0.2870495792324034*nul[1]*fr[27])-0.1657281518405969*nul[0]*fr[27]-0.2870495792324034*nul[1]*fl[27]-0.1657281518405969*nul[0]*fl[27]+0.2983106733130745*nul[1]*fr[22]+0.172229747539442*nul[0]*fr[22]-0.2983106733130745*nul[1]*fl[22]-0.172229747539442*nul[0]*fl[22]-0.2870495792324034*nul[7]*fr[20]-0.1657281518405969*nul[3]*fr[20]-0.2870495792324034*nul[7]*fl[20]-0.1657281518405969*nul[3]*fl[20]+0.2983106733130745*nul[7]*fr[13]+0.172229747539442*nul[3]*fr[13]-0.2983106733130745*nul[7]*fl[13]-0.172229747539442*nul[3]*fl[13]; 
  incr1[28] = (-0.2870495792324034*nul[7]*fr[31])-0.1657281518405969*nul[3]*fr[31]-0.2870495792324034*nul[7]*fl[31]-0.1657281518405969*nul[3]*fl[31]+0.2983106733130745*nul[7]*fr[30]+0.172229747539442*nul[3]*fr[30]-0.2983106733130745*nul[7]*fl[30]-0.172229747539442*nul[3]*fl[30]-0.2870495792324034*nul[1]*fr[28]-0.1657281518405969*nul[0]*fr[28]-0.2870495792324034*nul[1]*fl[28]-0.1657281518405969*nul[0]*fl[28]+0.2983106733130745*nul[1]*fr[24]+0.172229747539442*nul[0]*fr[24]-0.2983106733130745*nul[1]*fl[24]-0.172229747539442*nul[0]*fl[24]; 
  incr1[29] = (-0.2870495792324034*nul[1]*fr[29])-0.1657281518405969*nul[0]*fr[29]-0.2870495792324034*nul[1]*fl[29]-0.1657281518405969*nul[0]*fl[29]+0.2983106733130745*nul[1]*fr[25]+0.172229747539442*nul[0]*fr[25]-0.2983106733130745*nul[1]*fl[25]-0.172229747539442*nul[0]*fl[25]-0.2870495792324034*nul[7]*fr[23]-0.1657281518405969*nul[3]*fr[23]-0.2870495792324034*nul[7]*fl[23]-0.1657281518405969*nul[3]*fl[23]+0.2983106733130745*nul[7]*fr[15]+0.172229747539442*nul[3]*fr[15]-0.2983106733130745*nul[7]*fl[15]-0.172229747539442*nul[3]*fl[15]; 
  incr1[30] = 0.1657281518405969*nul[1]*fr[31]+0.09568319307746781*nul[0]*fr[31]+0.1657281518405969*nul[1]*fl[31]+0.09568319307746781*nul[0]*fl[31]-0.172229747539442*nul[1]*fr[30]-0.09943689110435816*nul[0]*fr[30]+0.172229747539442*nul[1]*fl[30]+0.09943689110435816*nul[0]*fl[30]+0.1657281518405969*nul[7]*fr[28]+0.09568319307746781*nul[3]*fr[28]+0.1657281518405969*nul[7]*fl[28]+0.09568319307746781*nul[3]*fl[28]-0.172229747539442*nul[7]*fr[24]-0.09943689110435816*nul[3]*fr[24]+0.172229747539442*nul[7]*fl[24]+0.09943689110435816*nul[3]*fl[24]; 
  incr1[31] = (-0.2870495792324034*nul[1]*fr[31])-0.1657281518405969*nul[0]*fr[31]-0.2870495792324034*nul[1]*fl[31]-0.1657281518405969*nul[0]*fl[31]+0.2983106733130745*nul[1]*fr[30]+0.172229747539442*nul[0]*fr[30]-0.2983106733130745*nul[1]*fl[30]-0.172229747539442*nul[0]*fl[30]-0.2870495792324034*nul[7]*fr[28]-0.1657281518405969*nul[3]*fr[28]-0.2870495792324034*nul[7]*fl[28]-0.1657281518405969*nul[3]*fl[28]+0.2983106733130745*nul[7]*fr[24]+0.172229747539442*nul[3]*fr[24]-0.2983106733130745*nul[7]*fl[24]-0.172229747539442*nul[3]*fl[24]; 

  incr2[1] = (-0.1530931089239486*fr[7]*nul[7])+0.1530931089239486*fl[7]*nul[7]+0.1325825214724776*fr[3]*nul[7]+0.1325825214724776*fl[3]*nul[7]-0.0883883476483184*nul[3]*fr[7]+0.0883883476483184*nul[3]*fl[7]+0.07654655446197427*fr[3]*nul[3]+0.07654655446197427*fl[3]*nul[3]-0.1530931089239486*fr[1]*nul[1]+0.1530931089239486*fl[1]*nul[1]+0.1325825214724776*fr[0]*nul[1]+0.1325825214724776*fl[0]*nul[1]-0.0883883476483184*nul[0]*fr[1]+0.0883883476483184*nul[0]*fl[1]+0.07654655446197427*fr[0]*nul[0]+0.07654655446197427*fl[0]*nul[0]; 
  incr2[6] = (-0.1530931089239486*nul[7]*fr[16])-0.0883883476483184*nul[3]*fr[16]+0.1530931089239486*nul[7]*fl[16]+0.0883883476483184*nul[3]*fl[16]+0.1325825214724776*nul[7]*fr[8]+0.07654655446197427*nul[3]*fr[8]+0.1325825214724776*nul[7]*fl[8]+0.07654655446197427*nul[3]*fl[8]-0.1530931089239486*nul[1]*fr[6]-0.0883883476483184*nul[0]*fr[6]+0.1530931089239486*nul[1]*fl[6]+0.0883883476483184*nul[0]*fl[6]+0.1325825214724776*nul[1]*fr[2]+0.07654655446197427*nul[0]*fr[2]+0.1325825214724776*nul[1]*fl[2]+0.07654655446197427*nul[0]*fl[2]; 
  incr2[7] = (-0.1530931089239486*fr[1]*nul[7])+0.1530931089239486*fl[1]*nul[7]+0.1325825214724776*fr[0]*nul[7]+0.1325825214724776*fl[0]*nul[7]-0.1530931089239486*nul[1]*fr[7]-0.0883883476483184*nul[0]*fr[7]+0.1530931089239486*nul[1]*fl[7]+0.0883883476483184*nul[0]*fl[7]-0.0883883476483184*fr[1]*nul[3]+0.0883883476483184*fl[1]*nul[3]+0.07654655446197427*fr[0]*nul[3]+0.07654655446197427*fl[0]*nul[3]+0.1325825214724776*nul[1]*fr[3]+0.07654655446197427*nul[0]*fr[3]+0.1325825214724776*nul[1]*fl[3]+0.07654655446197427*nul[0]*fl[3]; 
  incr2[9] = (-0.1530931089239486*nul[7]*fr[18])-0.0883883476483184*nul[3]*fr[18]+0.1530931089239486*nul[7]*fl[18]+0.0883883476483184*nul[3]*fl[18]+0.1325825214724776*nul[7]*fr[11]+0.07654655446197427*nul[3]*fr[11]+0.1325825214724776*nul[7]*fl[11]+0.07654655446197427*nul[3]*fl[11]-0.1530931089239486*nul[1]*fr[9]-0.0883883476483184*nul[0]*fr[9]+0.1530931089239486*nul[1]*fl[9]+0.0883883476483184*nul[0]*fl[9]+0.1325825214724776*nul[1]*fr[4]+0.07654655446197427*nul[0]*fr[4]+0.1325825214724776*nul[1]*fl[4]+0.07654655446197427*nul[0]*fl[4]; 
  incr2[12] = (-0.1530931089239486*nul[7]*fr[21])-0.0883883476483184*nul[3]*fr[21]+0.1530931089239486*nul[7]*fl[21]+0.0883883476483184*nul[3]*fl[21]+0.1325825214724776*nul[7]*fr[14]+0.07654655446197427*nul[3]*fr[14]+0.1325825214724776*nul[7]*fl[14]+0.07654655446197427*nul[3]*fl[14]-0.1530931089239486*nul[1]*fr[12]-0.0883883476483184*nul[0]*fr[12]+0.1530931089239486*nul[1]*fl[12]+0.0883883476483184*nul[0]*fl[12]+0.1325825214724776*nul[1]*fr[5]+0.07654655446197427*nul[0]*fr[5]+0.1325825214724776*nul[1]*fl[5]+0.07654655446197427*nul[0]*fl[5]; 
  incr2[16] = (-0.1530931089239486*nul[1]*fr[16])-0.0883883476483184*nul[0]*fr[16]+0.1530931089239486*nul[1]*fl[16]+0.0883883476483184*nul[0]*fl[16]+0.1325825214724776*nul[1]*fr[8]+0.07654655446197427*nul[0]*fr[8]+0.1325825214724776*nul[1]*fl[8]+0.07654655446197427*nul[0]*fl[8]-0.1530931089239486*fr[6]*nul[7]+0.1530931089239486*fl[6]*nul[7]+0.1325825214724776*fr[2]*nul[7]+0.1325825214724776*fl[2]*nul[7]-0.0883883476483184*nul[3]*fr[6]+0.0883883476483184*nul[3]*fl[6]+0.07654655446197427*fr[2]*nul[3]+0.07654655446197427*fl[2]*nul[3]; 
  incr2[17] = (-0.1530931089239486*nul[7]*fr[26])-0.0883883476483184*nul[3]*fr[26]+0.1530931089239486*nul[7]*fl[26]+0.0883883476483184*nul[3]*fl[26]+0.1325825214724776*nul[7]*fr[19]+0.07654655446197427*nul[3]*fr[19]+0.1325825214724776*nul[7]*fl[19]+0.07654655446197427*nul[3]*fl[19]-0.1530931089239486*nul[1]*fr[17]-0.0883883476483184*nul[0]*fr[17]+0.1530931089239486*nul[1]*fl[17]+0.0883883476483184*nul[0]*fl[17]+0.1325825214724776*nul[1]*fr[10]+0.07654655446197427*nul[0]*fr[10]+0.1325825214724776*nul[1]*fl[10]+0.07654655446197427*nul[0]*fl[10]; 
  incr2[18] = (-0.1530931089239486*nul[1]*fr[18])-0.0883883476483184*nul[0]*fr[18]+0.1530931089239486*nul[1]*fl[18]+0.0883883476483184*nul[0]*fl[18]+0.1325825214724776*nul[1]*fr[11]+0.07654655446197427*nul[0]*fr[11]+0.1325825214724776*nul[1]*fl[11]+0.07654655446197427*nul[0]*fl[11]-0.1530931089239486*nul[7]*fr[9]-0.0883883476483184*nul[3]*fr[9]+0.1530931089239486*nul[7]*fl[9]+0.0883883476483184*nul[3]*fl[9]+0.1325825214724776*fr[4]*nul[7]+0.1325825214724776*fl[4]*nul[7]+0.07654655446197427*nul[3]*fr[4]+0.07654655446197427*nul[3]*fl[4]; 
  incr2[20] = (-0.1530931089239486*nul[7]*fr[27])-0.0883883476483184*nul[3]*fr[27]+0.1530931089239486*nul[7]*fl[27]+0.0883883476483184*nul[3]*fl[27]+0.1325825214724776*nul[7]*fr[22]+0.07654655446197427*nul[3]*fr[22]+0.1325825214724776*nul[7]*fl[22]+0.07654655446197427*nul[3]*fl[22]-0.1530931089239486*nul[1]*fr[20]-0.0883883476483184*nul[0]*fr[20]+0.1530931089239486*nul[1]*fl[20]+0.0883883476483184*nul[0]*fl[20]+0.1325825214724776*nul[1]*fr[13]+0.07654655446197427*nul[0]*fr[13]+0.1325825214724776*nul[1]*fl[13]+0.07654655446197427*nul[0]*fl[13]; 
  incr2[21] = (-0.1530931089239486*nul[1]*fr[21])-0.0883883476483184*nul[0]*fr[21]+0.1530931089239486*nul[1]*fl[21]+0.0883883476483184*nul[0]*fl[21]+0.1325825214724776*nul[1]*fr[14]+0.07654655446197427*nul[0]*fr[14]+0.1325825214724776*nul[1]*fl[14]+0.07654655446197427*nul[0]*fl[14]-0.1530931089239486*nul[7]*fr[12]-0.0883883476483184*nul[3]*fr[12]+0.1530931089239486*nul[7]*fl[12]+0.0883883476483184*nul[3]*fl[12]+0.1325825214724776*fr[5]*nul[7]+0.1325825214724776*fl[5]*nul[7]+0.07654655446197427*nul[3]*fr[5]+0.07654655446197427*nul[3]*fl[5]; 
  incr2[23] = (-0.1530931089239486*nul[7]*fr[29])-0.0883883476483184*nul[3]*fr[29]+0.1530931089239486*nul[7]*fl[29]+0.0883883476483184*nul[3]*fl[29]+0.1325825214724776*nul[7]*fr[25]+0.07654655446197427*nul[3]*fr[25]+0.1325825214724776*nul[7]*fl[25]+0.07654655446197427*nul[3]*fl[25]-0.1530931089239486*nul[1]*fr[23]-0.0883883476483184*nul[0]*fr[23]+0.1530931089239486*nul[1]*fl[23]+0.0883883476483184*nul[0]*fl[23]+0.1325825214724776*nul[1]*fr[15]+0.07654655446197427*nul[0]*fr[15]+0.1325825214724776*nul[1]*fl[15]+0.07654655446197427*nul[0]*fl[15]; 
  incr2[26] = (-0.1530931089239486*nul[1]*fr[26])-0.0883883476483184*nul[0]*fr[26]+0.1530931089239486*nul[1]*fl[26]+0.0883883476483184*nul[0]*fl[26]+0.1325825214724776*nul[1]*fr[19]+0.07654655446197427*nul[0]*fr[19]+0.1325825214724776*nul[1]*fl[19]+0.07654655446197427*nul[0]*fl[19]-0.1530931089239486*nul[7]*fr[17]-0.0883883476483184*nul[3]*fr[17]+0.1530931089239486*nul[7]*fl[17]+0.0883883476483184*nul[3]*fl[17]+0.1325825214724776*nul[7]*fr[10]+0.07654655446197427*nul[3]*fr[10]+0.1325825214724776*nul[7]*fl[10]+0.07654655446197427*nul[3]*fl[10]; 
  incr2[27] = (-0.1530931089239486*nul[1]*fr[27])-0.0883883476483184*nul[0]*fr[27]+0.1530931089239486*nul[1]*fl[27]+0.0883883476483184*nul[0]*fl[27]+0.1325825214724776*nul[1]*fr[22]+0.07654655446197427*nul[0]*fr[22]+0.1325825214724776*nul[1]*fl[22]+0.07654655446197427*nul[0]*fl[22]-0.1530931089239486*nul[7]*fr[20]-0.0883883476483184*nul[3]*fr[20]+0.1530931089239486*nul[7]*fl[20]+0.0883883476483184*nul[3]*fl[20]+0.1325825214724776*nul[7]*fr[13]+0.07654655446197427*nul[3]*fr[13]+0.1325825214724776*nul[7]*fl[13]+0.07654655446197427*nul[3]*fl[13]; 
  incr2[28] = (-0.1530931089239486*nul[7]*fr[31])-0.0883883476483184*nul[3]*fr[31]+0.1530931089239486*nul[7]*fl[31]+0.0883883476483184*nul[3]*fl[31]+0.1325825214724776*nul[7]*fr[30]+0.07654655446197427*nul[3]*fr[30]+0.1325825214724776*nul[7]*fl[30]+0.07654655446197427*nul[3]*fl[30]-0.1530931089239486*nul[1]*fr[28]-0.0883883476483184*nul[0]*fr[28]+0.1530931089239486*nul[1]*fl[28]+0.0883883476483184*nul[0]*fl[28]+0.1325825214724776*nul[1]*fr[24]+0.07654655446197427*nul[0]*fr[24]+0.1325825214724776*nul[1]*fl[24]+0.07654655446197427*nul[0]*fl[24]; 
  incr2[29] = (-0.1530931089239486*nul[1]*fr[29])-0.0883883476483184*nul[0]*fr[29]+0.1530931089239486*nul[1]*fl[29]+0.0883883476483184*nul[0]*fl[29]+0.1325825214724776*nul[1]*fr[25]+0.07654655446197427*nul[0]*fr[25]+0.1325825214724776*nul[1]*fl[25]+0.07654655446197427*nul[0]*fl[25]-0.1530931089239486*nul[7]*fr[23]-0.0883883476483184*nul[3]*fr[23]+0.1530931089239486*nul[7]*fl[23]+0.0883883476483184*nul[3]*fl[23]+0.1325825214724776*nul[7]*fr[15]+0.07654655446197427*nul[3]*fr[15]+0.1325825214724776*nul[7]*fl[15]+0.07654655446197427*nul[3]*fl[15]; 
  incr2[31] = (-0.1530931089239486*nul[1]*fr[31])-0.0883883476483184*nul[0]*fr[31]+0.1530931089239486*nul[1]*fl[31]+0.0883883476483184*nul[0]*fl[31]+0.1325825214724776*nul[1]*fr[30]+0.07654655446197427*nul[0]*fr[30]+0.1325825214724776*nul[1]*fl[30]+0.07654655446197427*nul[0]*fl[30]-0.1530931089239486*nul[7]*fr[28]-0.0883883476483184*nul[3]*fr[28]+0.1530931089239486*nul[7]*fl[28]+0.0883883476483184*nul[3]*fl[28]+0.1325825214724776*nul[7]*fr[24]+0.07654655446197427*nul[3]*fr[24]+0.1325825214724776*nul[7]*fl[24]+0.07654655446197427*nul[3]*fl[24]; 

  outr[0] += incr1[0]*rdxFr; 
  outr[1] += incr2[1]*rdxFr+incr1[1]*rdxFr; 
  outr[2] += incr1[2]*rdxFr; 
  outr[3] += incr1[3]*rdxFr; 
  outr[4] += incr1[4]*rdxFr; 
  outr[5] += incr1[5]*rdxFr; 
  outr[6] += incr2[6]*rdxFr+incr1[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr+incr1[7]*rdxFr; 
  outr[8] += incr1[8]*rdxFr; 
  outr[9] += incr2[9]*rdxFr+incr1[9]*rdxFr; 
  outr[10] += incr1[10]*rdxFr; 
  outr[11] += incr1[11]*rdxFr; 
  outr[12] += incr2[12]*rdxFr+incr1[12]*rdxFr; 
  outr[13] += incr1[13]*rdxFr; 
  outr[14] += incr1[14]*rdxFr; 
  outr[15] += incr1[15]*rdxFr; 
  outr[16] += incr2[16]*rdxFr+incr1[16]*rdxFr; 
  outr[17] += incr2[17]*rdxFr+incr1[17]*rdxFr; 
  outr[18] += incr2[18]*rdxFr+incr1[18]*rdxFr; 
  outr[19] += incr1[19]*rdxFr; 
  outr[20] += incr2[20]*rdxFr+incr1[20]*rdxFr; 
  outr[21] += incr2[21]*rdxFr+incr1[21]*rdxFr; 
  outr[22] += incr1[22]*rdxFr; 
  outr[23] += incr2[23]*rdxFr+incr1[23]*rdxFr; 
  outr[24] += incr1[24]*rdxFr; 
  outr[25] += incr1[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr+incr1[26]*rdxFr; 
  outr[27] += incr2[27]*rdxFr+incr1[27]*rdxFr; 
  outr[28] += incr2[28]*rdxFr+incr1[28]*rdxFr; 
  outr[29] += incr2[29]*rdxFr+incr1[29]*rdxFr; 
  outr[30] += incr1[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr+incr1[31]*rdxFr; 

  outl[0] += -1.0*incr1[0]*rdxFl; 
  outl[1] += incr1[1]*rdxFl-1.0*incr2[1]*rdxFl; 
  outl[2] += -1.0*incr1[2]*rdxFl; 
  outl[3] += -1.0*incr1[3]*rdxFl; 
  outl[4] += -1.0*incr1[4]*rdxFl; 
  outl[5] += -1.0*incr1[5]*rdxFl; 
  outl[6] += incr1[6]*rdxFl-1.0*incr2[6]*rdxFl; 
  outl[7] += incr1[7]*rdxFl-1.0*incr2[7]*rdxFl; 
  outl[8] += -1.0*incr1[8]*rdxFl; 
  outl[9] += incr1[9]*rdxFl-1.0*incr2[9]*rdxFl; 
  outl[10] += -1.0*incr1[10]*rdxFl; 
  outl[11] += -1.0*incr1[11]*rdxFl; 
  outl[12] += incr1[12]*rdxFl-1.0*incr2[12]*rdxFl; 
  outl[13] += -1.0*incr1[13]*rdxFl; 
  outl[14] += -1.0*incr1[14]*rdxFl; 
  outl[15] += -1.0*incr1[15]*rdxFl; 
  outl[16] += incr1[16]*rdxFl-1.0*incr2[16]*rdxFl; 
  outl[17] += incr1[17]*rdxFl-1.0*incr2[17]*rdxFl; 
  outl[18] += incr1[18]*rdxFl-1.0*incr2[18]*rdxFl; 
  outl[19] += -1.0*incr1[19]*rdxFl; 
  outl[20] += incr1[20]*rdxFl-1.0*incr2[20]*rdxFl; 
  outl[21] += incr1[21]*rdxFl-1.0*incr2[21]*rdxFl; 
  outl[22] += -1.0*incr1[22]*rdxFl; 
  outl[23] += incr1[23]*rdxFl-1.0*incr2[23]*rdxFl; 
  outl[24] += -1.0*incr1[24]*rdxFl; 
  outl[25] += -1.0*incr1[25]*rdxFl; 
  outl[26] += incr1[26]*rdxFl-1.0*incr2[26]*rdxFl; 
  outl[27] += incr1[27]*rdxFl-1.0*incr2[27]*rdxFl; 
  outl[28] += incr1[28]*rdxFl-1.0*incr2[28]*rdxFl; 
  outl[29] += incr1[29]*rdxFl-1.0*incr2[29]*rdxFl; 
  outl[30] += -1.0*incr1[30]*rdxFl; 
  outl[31] += incr1[31]*rdxFl-1.0*incr2[31]*rdxFl; 

} 
void ConstDiffusionVarCoeffSurf5xSerP1_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[160]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[1]*dxl[1]); 
  double rdxFr = 4.0/(dxr[1]*dxr[1]); 

  double incr1[32]; 
  double incr2[32]; 

  incr1[0] = 0.09568319307746781*fr[16]*nul[39]+0.09568319307746781*fl[16]*nul[39]-0.09943689110435816*fr[7]*nul[39]+0.09943689110435816*fl[7]*nul[39]+0.09568319307746781*fr[8]*nul[35]+0.09568319307746781*fl[8]*nul[35]-0.09943689110435816*fr[3]*nul[35]+0.09943689110435816*fl[3]*nul[35]+0.09568319307746781*fr[6]*nul[33]+0.09568319307746781*fl[6]*nul[33]-0.09943689110435816*fr[1]*nul[33]+0.09943689110435816*fl[1]*nul[33]+0.09568319307746781*fr[2]*nul[32]+0.09568319307746781*fl[2]*nul[32]-0.09943689110435816*fr[0]*nul[32]+0.09943689110435816*fl[0]*nul[32]; 
  incr1[1] = 0.09568319307746781*fr[8]*nul[39]+0.09568319307746781*fl[8]*nul[39]-0.09943689110435816*fr[3]*nul[39]+0.09943689110435816*fl[3]*nul[39]+0.09568319307746781*fr[16]*nul[35]+0.09568319307746781*fl[16]*nul[35]-0.09943689110435816*fr[7]*nul[35]+0.09943689110435816*fl[7]*nul[35]+0.09568319307746781*fr[2]*nul[33]+0.09568319307746781*fl[2]*nul[33]-0.09943689110435816*fr[0]*nul[33]+0.09943689110435816*fl[0]*nul[33]+0.09568319307746781*fr[6]*nul[32]+0.09568319307746781*fl[6]*nul[32]-0.09943689110435816*fr[1]*nul[32]+0.09943689110435816*fl[1]*nul[32]; 
  incr1[2] = (-0.1657281518405969*fr[16]*nul[39])-0.1657281518405969*fl[16]*nul[39]+0.172229747539442*fr[7]*nul[39]-0.172229747539442*fl[7]*nul[39]-0.1657281518405969*fr[8]*nul[35]-0.1657281518405969*fl[8]*nul[35]+0.172229747539442*fr[3]*nul[35]-0.172229747539442*fl[3]*nul[35]-0.1657281518405969*fr[6]*nul[33]-0.1657281518405969*fl[6]*nul[33]+0.172229747539442*fr[1]*nul[33]-0.172229747539442*fl[1]*nul[33]-0.1657281518405969*fr[2]*nul[32]-0.1657281518405969*fl[2]*nul[32]+0.172229747539442*fr[0]*nul[32]-0.172229747539442*fl[0]*nul[32]; 
  incr1[3] = 0.09568319307746781*fr[6]*nul[39]+0.09568319307746781*fl[6]*nul[39]-0.09943689110435816*fr[1]*nul[39]+0.09943689110435816*fl[1]*nul[39]+0.09568319307746781*fr[2]*nul[35]+0.09568319307746781*fl[2]*nul[35]-0.09943689110435816*fr[0]*nul[35]+0.09943689110435816*fl[0]*nul[35]+0.09568319307746781*fr[16]*nul[33]+0.09568319307746781*fl[16]*nul[33]-0.09943689110435816*fr[7]*nul[33]+0.09943689110435816*fl[7]*nul[33]+0.09568319307746781*fr[8]*nul[32]+0.09568319307746781*fl[8]*nul[32]-0.09943689110435816*fr[3]*nul[32]+0.09943689110435816*fl[3]*nul[32]; 
  incr1[4] = 0.09568319307746781*fr[26]*nul[39]+0.09568319307746781*fl[26]*nul[39]-0.09943689110435816*fr[18]*nul[39]+0.09943689110435816*fl[18]*nul[39]+0.09568319307746781*fr[19]*nul[35]+0.09568319307746781*fl[19]*nul[35]-0.09943689110435816*fr[11]*nul[35]+0.09943689110435816*fl[11]*nul[35]+0.09568319307746781*fr[17]*nul[33]+0.09568319307746781*fl[17]*nul[33]-0.09943689110435816*fr[9]*nul[33]+0.09943689110435816*fl[9]*nul[33]+0.09568319307746781*fr[10]*nul[32]+0.09568319307746781*fl[10]*nul[32]-0.09943689110435816*fr[4]*nul[32]+0.09943689110435816*fl[4]*nul[32]; 
  incr1[5] = 0.09568319307746781*fr[27]*nul[39]+0.09568319307746781*fl[27]*nul[39]-0.09943689110435816*fr[21]*nul[39]+0.09943689110435816*fl[21]*nul[39]+0.09568319307746781*fr[22]*nul[35]+0.09568319307746781*fl[22]*nul[35]-0.09943689110435816*fr[14]*nul[35]+0.09943689110435816*fl[14]*nul[35]+0.09568319307746781*fr[20]*nul[33]+0.09568319307746781*fl[20]*nul[33]-0.09943689110435816*fr[12]*nul[33]+0.09943689110435816*fl[12]*nul[33]+0.09568319307746781*fr[13]*nul[32]+0.09568319307746781*fl[13]*nul[32]-0.09943689110435816*fr[5]*nul[32]+0.09943689110435816*fl[5]*nul[32]; 
  incr1[6] = (-0.1657281518405969*fr[8]*nul[39])-0.1657281518405969*fl[8]*nul[39]+0.172229747539442*fr[3]*nul[39]-0.172229747539442*fl[3]*nul[39]-0.1657281518405969*fr[16]*nul[35]-0.1657281518405969*fl[16]*nul[35]+0.172229747539442*fr[7]*nul[35]-0.172229747539442*fl[7]*nul[35]-0.1657281518405969*fr[2]*nul[33]-0.1657281518405969*fl[2]*nul[33]+0.172229747539442*fr[0]*nul[33]-0.172229747539442*fl[0]*nul[33]-0.1657281518405969*fr[6]*nul[32]-0.1657281518405969*fl[6]*nul[32]+0.172229747539442*fr[1]*nul[32]-0.172229747539442*fl[1]*nul[32]; 
  incr1[7] = 0.09568319307746781*fr[2]*nul[39]+0.09568319307746781*fl[2]*nul[39]-0.09943689110435816*fr[0]*nul[39]+0.09943689110435816*fl[0]*nul[39]+0.09568319307746781*fr[6]*nul[35]+0.09568319307746781*fl[6]*nul[35]-0.09943689110435816*fr[1]*nul[35]+0.09943689110435816*fl[1]*nul[35]+0.09568319307746781*fr[8]*nul[33]+0.09568319307746781*fl[8]*nul[33]-0.09943689110435816*fr[3]*nul[33]+0.09943689110435816*fl[3]*nul[33]+0.09568319307746781*fr[16]*nul[32]+0.09568319307746781*fl[16]*nul[32]-0.09943689110435816*fr[7]*nul[32]+0.09943689110435816*fl[7]*nul[32]; 
  incr1[8] = (-0.1657281518405969*fr[6]*nul[39])-0.1657281518405969*fl[6]*nul[39]+0.172229747539442*fr[1]*nul[39]-0.172229747539442*fl[1]*nul[39]-0.1657281518405969*fr[2]*nul[35]-0.1657281518405969*fl[2]*nul[35]+0.172229747539442*fr[0]*nul[35]-0.172229747539442*fl[0]*nul[35]-0.1657281518405969*fr[16]*nul[33]-0.1657281518405969*fl[16]*nul[33]+0.172229747539442*fr[7]*nul[33]-0.172229747539442*fl[7]*nul[33]-0.1657281518405969*fr[8]*nul[32]-0.1657281518405969*fl[8]*nul[32]+0.172229747539442*fr[3]*nul[32]-0.172229747539442*fl[3]*nul[32]; 
  incr1[9] = 0.09568319307746781*fr[19]*nul[39]+0.09568319307746781*fl[19]*nul[39]-0.09943689110435816*fr[11]*nul[39]+0.09943689110435816*fl[11]*nul[39]+0.09568319307746781*fr[26]*nul[35]+0.09568319307746781*fl[26]*nul[35]-0.09943689110435816*fr[18]*nul[35]+0.09943689110435816*fl[18]*nul[35]+0.09568319307746781*fr[10]*nul[33]+0.09568319307746781*fl[10]*nul[33]-0.09943689110435816*fr[4]*nul[33]+0.09943689110435816*fl[4]*nul[33]+0.09568319307746781*fr[17]*nul[32]+0.09568319307746781*fl[17]*nul[32]-0.09943689110435816*fr[9]*nul[32]+0.09943689110435816*fl[9]*nul[32]; 
  incr1[10] = (-0.1657281518405969*fr[26]*nul[39])-0.1657281518405969*fl[26]*nul[39]+0.172229747539442*fr[18]*nul[39]-0.172229747539442*fl[18]*nul[39]-0.1657281518405969*fr[19]*nul[35]-0.1657281518405969*fl[19]*nul[35]+0.172229747539442*fr[11]*nul[35]-0.172229747539442*fl[11]*nul[35]-0.1657281518405969*fr[17]*nul[33]-0.1657281518405969*fl[17]*nul[33]+0.172229747539442*fr[9]*nul[33]-0.172229747539442*fl[9]*nul[33]-0.1657281518405969*fr[10]*nul[32]-0.1657281518405969*fl[10]*nul[32]+0.172229747539442*fr[4]*nul[32]-0.172229747539442*fl[4]*nul[32]; 
  incr1[11] = 0.09568319307746781*fr[17]*nul[39]+0.09568319307746781*fl[17]*nul[39]-0.09943689110435816*fr[9]*nul[39]+0.09943689110435816*fl[9]*nul[39]+0.09568319307746781*fr[10]*nul[35]+0.09568319307746781*fl[10]*nul[35]-0.09943689110435816*fr[4]*nul[35]+0.09943689110435816*fl[4]*nul[35]+0.09568319307746781*fr[26]*nul[33]+0.09568319307746781*fl[26]*nul[33]-0.09943689110435816*fr[18]*nul[33]+0.09943689110435816*fl[18]*nul[33]+0.09568319307746781*fr[19]*nul[32]+0.09568319307746781*fl[19]*nul[32]-0.09943689110435816*fr[11]*nul[32]+0.09943689110435816*fl[11]*nul[32]; 
  incr1[12] = 0.09568319307746781*fr[22]*nul[39]+0.09568319307746781*fl[22]*nul[39]-0.09943689110435816*fr[14]*nul[39]+0.09943689110435816*fl[14]*nul[39]+0.09568319307746781*fr[27]*nul[35]+0.09568319307746781*fl[27]*nul[35]-0.09943689110435816*fr[21]*nul[35]+0.09943689110435816*fl[21]*nul[35]+0.09568319307746781*fr[13]*nul[33]+0.09568319307746781*fl[13]*nul[33]-0.09943689110435816*fr[5]*nul[33]+0.09943689110435816*fl[5]*nul[33]+0.09568319307746781*fr[20]*nul[32]+0.09568319307746781*fl[20]*nul[32]-0.09943689110435816*fr[12]*nul[32]+0.09943689110435816*fl[12]*nul[32]; 
  incr1[13] = (-0.1657281518405969*fr[27]*nul[39])-0.1657281518405969*fl[27]*nul[39]+0.172229747539442*fr[21]*nul[39]-0.172229747539442*fl[21]*nul[39]-0.1657281518405969*fr[22]*nul[35]-0.1657281518405969*fl[22]*nul[35]+0.172229747539442*fr[14]*nul[35]-0.172229747539442*fl[14]*nul[35]-0.1657281518405969*fr[20]*nul[33]-0.1657281518405969*fl[20]*nul[33]+0.172229747539442*fr[12]*nul[33]-0.172229747539442*fl[12]*nul[33]-0.1657281518405969*fr[13]*nul[32]-0.1657281518405969*fl[13]*nul[32]+0.172229747539442*fr[5]*nul[32]-0.172229747539442*fl[5]*nul[32]; 
  incr1[14] = 0.09568319307746781*fr[20]*nul[39]+0.09568319307746781*fl[20]*nul[39]-0.09943689110435816*fr[12]*nul[39]+0.09943689110435816*fl[12]*nul[39]+0.09568319307746781*fr[13]*nul[35]+0.09568319307746781*fl[13]*nul[35]-0.09943689110435816*fr[5]*nul[35]+0.09943689110435816*fl[5]*nul[35]+0.09568319307746781*fr[27]*nul[33]+0.09568319307746781*fl[27]*nul[33]-0.09943689110435816*fr[21]*nul[33]+0.09943689110435816*fl[21]*nul[33]+0.09568319307746781*fr[22]*nul[32]+0.09568319307746781*fl[22]*nul[32]-0.09943689110435816*fr[14]*nul[32]+0.09943689110435816*fl[14]*nul[32]; 
  incr1[15] = 0.09568319307746781*fr[31]*nul[39]+0.09568319307746781*fl[31]*nul[39]-0.09943689110435816*fr[29]*nul[39]+0.09943689110435816*fl[29]*nul[39]+0.09568319307746781*fr[30]*nul[35]+0.09568319307746781*fl[30]*nul[35]-0.09943689110435816*fr[25]*nul[35]+0.09943689110435816*fl[25]*nul[35]+0.09568319307746781*fr[28]*nul[33]+0.09568319307746781*fl[28]*nul[33]-0.09943689110435816*fr[23]*nul[33]+0.09943689110435816*fl[23]*nul[33]+0.09568319307746781*fr[24]*nul[32]+0.09568319307746781*fl[24]*nul[32]-0.09943689110435816*fr[15]*nul[32]+0.09943689110435816*fl[15]*nul[32]; 
  incr1[16] = (-0.1657281518405969*fr[2]*nul[39])-0.1657281518405969*fl[2]*nul[39]+0.172229747539442*fr[0]*nul[39]-0.172229747539442*fl[0]*nul[39]-0.1657281518405969*fr[6]*nul[35]-0.1657281518405969*fl[6]*nul[35]+0.172229747539442*fr[1]*nul[35]-0.172229747539442*fl[1]*nul[35]-0.1657281518405969*fr[8]*nul[33]-0.1657281518405969*fl[8]*nul[33]+0.172229747539442*fr[3]*nul[33]-0.172229747539442*fl[3]*nul[33]-0.1657281518405969*fr[16]*nul[32]-0.1657281518405969*fl[16]*nul[32]+0.172229747539442*fr[7]*nul[32]-0.172229747539442*fl[7]*nul[32]; 
  incr1[17] = (-0.1657281518405969*fr[19]*nul[39])-0.1657281518405969*fl[19]*nul[39]+0.172229747539442*fr[11]*nul[39]-0.172229747539442*fl[11]*nul[39]-0.1657281518405969*fr[26]*nul[35]-0.1657281518405969*fl[26]*nul[35]+0.172229747539442*fr[18]*nul[35]-0.172229747539442*fl[18]*nul[35]-0.1657281518405969*fr[10]*nul[33]-0.1657281518405969*fl[10]*nul[33]+0.172229747539442*fr[4]*nul[33]-0.172229747539442*fl[4]*nul[33]-0.1657281518405969*fr[17]*nul[32]-0.1657281518405969*fl[17]*nul[32]+0.172229747539442*fr[9]*nul[32]-0.172229747539442*fl[9]*nul[32]; 
  incr1[18] = 0.09568319307746781*fr[10]*nul[39]+0.09568319307746781*fl[10]*nul[39]-0.09943689110435816*fr[4]*nul[39]+0.09943689110435816*fl[4]*nul[39]+0.09568319307746781*fr[17]*nul[35]+0.09568319307746781*fl[17]*nul[35]-0.09943689110435816*fr[9]*nul[35]+0.09943689110435816*fl[9]*nul[35]+0.09568319307746781*fr[19]*nul[33]+0.09568319307746781*fl[19]*nul[33]-0.09943689110435816*fr[11]*nul[33]+0.09943689110435816*fl[11]*nul[33]+0.09568319307746781*fr[26]*nul[32]+0.09568319307746781*fl[26]*nul[32]-0.09943689110435816*fr[18]*nul[32]+0.09943689110435816*fl[18]*nul[32]; 
  incr1[19] = (-0.1657281518405969*fr[17]*nul[39])-0.1657281518405969*fl[17]*nul[39]+0.172229747539442*fr[9]*nul[39]-0.172229747539442*fl[9]*nul[39]-0.1657281518405969*fr[10]*nul[35]-0.1657281518405969*fl[10]*nul[35]+0.172229747539442*fr[4]*nul[35]-0.172229747539442*fl[4]*nul[35]-0.1657281518405969*fr[26]*nul[33]-0.1657281518405969*fl[26]*nul[33]+0.172229747539442*fr[18]*nul[33]-0.172229747539442*fl[18]*nul[33]-0.1657281518405969*fr[19]*nul[32]-0.1657281518405969*fl[19]*nul[32]+0.172229747539442*fr[11]*nul[32]-0.172229747539442*fl[11]*nul[32]; 
  incr1[20] = (-0.1657281518405969*fr[22]*nul[39])-0.1657281518405969*fl[22]*nul[39]+0.172229747539442*fr[14]*nul[39]-0.172229747539442*fl[14]*nul[39]-0.1657281518405969*fr[27]*nul[35]-0.1657281518405969*fl[27]*nul[35]+0.172229747539442*fr[21]*nul[35]-0.172229747539442*fl[21]*nul[35]-0.1657281518405969*fr[13]*nul[33]-0.1657281518405969*fl[13]*nul[33]+0.172229747539442*fr[5]*nul[33]-0.172229747539442*fl[5]*nul[33]-0.1657281518405969*fr[20]*nul[32]-0.1657281518405969*fl[20]*nul[32]+0.172229747539442*fr[12]*nul[32]-0.172229747539442*fl[12]*nul[32]; 
  incr1[21] = 0.09568319307746781*fr[13]*nul[39]+0.09568319307746781*fl[13]*nul[39]-0.09943689110435816*fr[5]*nul[39]+0.09943689110435816*fl[5]*nul[39]+0.09568319307746781*fr[20]*nul[35]+0.09568319307746781*fl[20]*nul[35]-0.09943689110435816*fr[12]*nul[35]+0.09943689110435816*fl[12]*nul[35]+0.09568319307746781*fr[22]*nul[33]+0.09568319307746781*fl[22]*nul[33]-0.09943689110435816*fr[14]*nul[33]+0.09943689110435816*fl[14]*nul[33]+0.09568319307746781*fr[27]*nul[32]+0.09568319307746781*fl[27]*nul[32]-0.09943689110435816*fr[21]*nul[32]+0.09943689110435816*fl[21]*nul[32]; 
  incr1[22] = (-0.1657281518405969*fr[20]*nul[39])-0.1657281518405969*fl[20]*nul[39]+0.172229747539442*fr[12]*nul[39]-0.172229747539442*fl[12]*nul[39]-0.1657281518405969*fr[13]*nul[35]-0.1657281518405969*fl[13]*nul[35]+0.172229747539442*fr[5]*nul[35]-0.172229747539442*fl[5]*nul[35]-0.1657281518405969*fr[27]*nul[33]-0.1657281518405969*fl[27]*nul[33]+0.172229747539442*fr[21]*nul[33]-0.172229747539442*fl[21]*nul[33]-0.1657281518405969*fr[22]*nul[32]-0.1657281518405969*fl[22]*nul[32]+0.172229747539442*fr[14]*nul[32]-0.172229747539442*fl[14]*nul[32]; 
  incr1[23] = 0.09568319307746781*fr[30]*nul[39]+0.09568319307746781*fl[30]*nul[39]-0.09943689110435816*fr[25]*nul[39]+0.09943689110435816*fl[25]*nul[39]+0.09568319307746781*fr[31]*nul[35]+0.09568319307746781*fl[31]*nul[35]-0.09943689110435816*fr[29]*nul[35]+0.09943689110435816*fl[29]*nul[35]+0.09568319307746781*fr[24]*nul[33]+0.09568319307746781*fl[24]*nul[33]-0.09943689110435816*fr[15]*nul[33]+0.09943689110435816*fl[15]*nul[33]+0.09568319307746781*fr[28]*nul[32]+0.09568319307746781*fl[28]*nul[32]-0.09943689110435816*fr[23]*nul[32]+0.09943689110435816*fl[23]*nul[32]; 
  incr1[24] = (-0.1657281518405969*fr[31]*nul[39])-0.1657281518405969*fl[31]*nul[39]+0.172229747539442*fr[29]*nul[39]-0.172229747539442*fl[29]*nul[39]-0.1657281518405969*fr[30]*nul[35]-0.1657281518405969*fl[30]*nul[35]+0.172229747539442*fr[25]*nul[35]-0.172229747539442*fl[25]*nul[35]-0.1657281518405969*fr[28]*nul[33]-0.1657281518405969*fl[28]*nul[33]+0.172229747539442*fr[23]*nul[33]-0.172229747539442*fl[23]*nul[33]-0.1657281518405969*fr[24]*nul[32]-0.1657281518405969*fl[24]*nul[32]+0.172229747539442*fr[15]*nul[32]-0.172229747539442*fl[15]*nul[32]; 
  incr1[25] = 0.09568319307746781*fr[28]*nul[39]+0.09568319307746781*fl[28]*nul[39]-0.09943689110435816*fr[23]*nul[39]+0.09943689110435816*fl[23]*nul[39]+0.09568319307746781*fr[24]*nul[35]+0.09568319307746781*fl[24]*nul[35]-0.09943689110435816*fr[15]*nul[35]+0.09943689110435816*fl[15]*nul[35]+0.09568319307746781*fr[31]*nul[33]+0.09568319307746781*fl[31]*nul[33]-0.09943689110435816*fr[29]*nul[33]+0.09943689110435816*fl[29]*nul[33]+0.09568319307746781*fr[30]*nul[32]+0.09568319307746781*fl[30]*nul[32]-0.09943689110435816*fr[25]*nul[32]+0.09943689110435816*fl[25]*nul[32]; 
  incr1[26] = (-0.1657281518405969*fr[10]*nul[39])-0.1657281518405969*fl[10]*nul[39]+0.172229747539442*fr[4]*nul[39]-0.172229747539442*fl[4]*nul[39]-0.1657281518405969*fr[17]*nul[35]-0.1657281518405969*fl[17]*nul[35]+0.172229747539442*fr[9]*nul[35]-0.172229747539442*fl[9]*nul[35]-0.1657281518405969*fr[19]*nul[33]-0.1657281518405969*fl[19]*nul[33]+0.172229747539442*fr[11]*nul[33]-0.172229747539442*fl[11]*nul[33]-0.1657281518405969*fr[26]*nul[32]-0.1657281518405969*fl[26]*nul[32]+0.172229747539442*fr[18]*nul[32]-0.172229747539442*fl[18]*nul[32]; 
  incr1[27] = (-0.1657281518405969*fr[13]*nul[39])-0.1657281518405969*fl[13]*nul[39]+0.172229747539442*fr[5]*nul[39]-0.172229747539442*fl[5]*nul[39]-0.1657281518405969*fr[20]*nul[35]-0.1657281518405969*fl[20]*nul[35]+0.172229747539442*fr[12]*nul[35]-0.172229747539442*fl[12]*nul[35]-0.1657281518405969*fr[22]*nul[33]-0.1657281518405969*fl[22]*nul[33]+0.172229747539442*fr[14]*nul[33]-0.172229747539442*fl[14]*nul[33]-0.1657281518405969*fr[27]*nul[32]-0.1657281518405969*fl[27]*nul[32]+0.172229747539442*fr[21]*nul[32]-0.172229747539442*fl[21]*nul[32]; 
  incr1[28] = (-0.1657281518405969*fr[30]*nul[39])-0.1657281518405969*fl[30]*nul[39]+0.172229747539442*fr[25]*nul[39]-0.172229747539442*fl[25]*nul[39]-0.1657281518405969*fr[31]*nul[35]-0.1657281518405969*fl[31]*nul[35]+0.172229747539442*fr[29]*nul[35]-0.172229747539442*fl[29]*nul[35]-0.1657281518405969*fr[24]*nul[33]-0.1657281518405969*fl[24]*nul[33]+0.172229747539442*fr[15]*nul[33]-0.172229747539442*fl[15]*nul[33]-0.1657281518405969*fr[28]*nul[32]-0.1657281518405969*fl[28]*nul[32]+0.172229747539442*fr[23]*nul[32]-0.172229747539442*fl[23]*nul[32]; 
  incr1[29] = 0.09568319307746781*fr[24]*nul[39]+0.09568319307746781*fl[24]*nul[39]-0.09943689110435816*fr[15]*nul[39]+0.09943689110435816*fl[15]*nul[39]+0.09568319307746781*fr[28]*nul[35]+0.09568319307746781*fl[28]*nul[35]-0.09943689110435816*fr[23]*nul[35]+0.09943689110435816*fl[23]*nul[35]+0.09568319307746781*fr[30]*nul[33]+0.09568319307746781*fl[30]*nul[33]-0.09943689110435816*fr[25]*nul[33]+0.09943689110435816*fl[25]*nul[33]+0.09568319307746781*fr[31]*nul[32]+0.09568319307746781*fl[31]*nul[32]-0.09943689110435816*fr[29]*nul[32]+0.09943689110435816*fl[29]*nul[32]; 
  incr1[30] = (-0.1657281518405969*fr[28]*nul[39])-0.1657281518405969*fl[28]*nul[39]+0.172229747539442*fr[23]*nul[39]-0.172229747539442*fl[23]*nul[39]-0.1657281518405969*fr[24]*nul[35]-0.1657281518405969*fl[24]*nul[35]+0.172229747539442*fr[15]*nul[35]-0.172229747539442*fl[15]*nul[35]-0.1657281518405969*fr[31]*nul[33]-0.1657281518405969*fl[31]*nul[33]+0.172229747539442*fr[29]*nul[33]-0.172229747539442*fl[29]*nul[33]-0.1657281518405969*fr[30]*nul[32]-0.1657281518405969*fl[30]*nul[32]+0.172229747539442*fr[25]*nul[32]-0.172229747539442*fl[25]*nul[32]; 
  incr1[31] = (-0.1657281518405969*fr[24]*nul[39])-0.1657281518405969*fl[24]*nul[39]+0.172229747539442*fr[15]*nul[39]-0.172229747539442*fl[15]*nul[39]-0.1657281518405969*fr[28]*nul[35]-0.1657281518405969*fl[28]*nul[35]+0.172229747539442*fr[23]*nul[35]-0.172229747539442*fl[23]*nul[35]-0.1657281518405969*fr[30]*nul[33]-0.1657281518405969*fl[30]*nul[33]+0.172229747539442*fr[25]*nul[33]-0.172229747539442*fl[25]*nul[33]-0.1657281518405969*fr[31]*nul[32]-0.1657281518405969*fl[31]*nul[32]+0.172229747539442*fr[29]*nul[32]-0.172229747539442*fl[29]*nul[32]; 

  incr2[2] = (-0.0883883476483184*fr[16]*nul[39])+0.0883883476483184*fl[16]*nul[39]+0.07654655446197427*fr[7]*nul[39]+0.07654655446197427*fl[7]*nul[39]-0.0883883476483184*fr[8]*nul[35]+0.0883883476483184*fl[8]*nul[35]+0.07654655446197427*fr[3]*nul[35]+0.07654655446197427*fl[3]*nul[35]-0.0883883476483184*fr[6]*nul[33]+0.0883883476483184*fl[6]*nul[33]+0.07654655446197427*fr[1]*nul[33]+0.07654655446197427*fl[1]*nul[33]-0.0883883476483184*fr[2]*nul[32]+0.0883883476483184*fl[2]*nul[32]+0.07654655446197427*fr[0]*nul[32]+0.07654655446197427*fl[0]*nul[32]; 
  incr2[6] = (-0.0883883476483184*fr[8]*nul[39])+0.0883883476483184*fl[8]*nul[39]+0.07654655446197427*fr[3]*nul[39]+0.07654655446197427*fl[3]*nul[39]-0.0883883476483184*fr[16]*nul[35]+0.0883883476483184*fl[16]*nul[35]+0.07654655446197427*fr[7]*nul[35]+0.07654655446197427*fl[7]*nul[35]-0.0883883476483184*fr[2]*nul[33]+0.0883883476483184*fl[2]*nul[33]+0.07654655446197427*fr[0]*nul[33]+0.07654655446197427*fl[0]*nul[33]-0.0883883476483184*fr[6]*nul[32]+0.0883883476483184*fl[6]*nul[32]+0.07654655446197427*fr[1]*nul[32]+0.07654655446197427*fl[1]*nul[32]; 
  incr2[8] = (-0.0883883476483184*fr[6]*nul[39])+0.0883883476483184*fl[6]*nul[39]+0.07654655446197427*fr[1]*nul[39]+0.07654655446197427*fl[1]*nul[39]-0.0883883476483184*fr[2]*nul[35]+0.0883883476483184*fl[2]*nul[35]+0.07654655446197427*fr[0]*nul[35]+0.07654655446197427*fl[0]*nul[35]-0.0883883476483184*fr[16]*nul[33]+0.0883883476483184*fl[16]*nul[33]+0.07654655446197427*fr[7]*nul[33]+0.07654655446197427*fl[7]*nul[33]-0.0883883476483184*fr[8]*nul[32]+0.0883883476483184*fl[8]*nul[32]+0.07654655446197427*fr[3]*nul[32]+0.07654655446197427*fl[3]*nul[32]; 
  incr2[10] = (-0.0883883476483184*fr[26]*nul[39])+0.0883883476483184*fl[26]*nul[39]+0.07654655446197427*fr[18]*nul[39]+0.07654655446197427*fl[18]*nul[39]-0.0883883476483184*fr[19]*nul[35]+0.0883883476483184*fl[19]*nul[35]+0.07654655446197427*fr[11]*nul[35]+0.07654655446197427*fl[11]*nul[35]-0.0883883476483184*fr[17]*nul[33]+0.0883883476483184*fl[17]*nul[33]+0.07654655446197427*fr[9]*nul[33]+0.07654655446197427*fl[9]*nul[33]-0.0883883476483184*fr[10]*nul[32]+0.0883883476483184*fl[10]*nul[32]+0.07654655446197427*fr[4]*nul[32]+0.07654655446197427*fl[4]*nul[32]; 
  incr2[13] = (-0.0883883476483184*fr[27]*nul[39])+0.0883883476483184*fl[27]*nul[39]+0.07654655446197427*fr[21]*nul[39]+0.07654655446197427*fl[21]*nul[39]-0.0883883476483184*fr[22]*nul[35]+0.0883883476483184*fl[22]*nul[35]+0.07654655446197427*fr[14]*nul[35]+0.07654655446197427*fl[14]*nul[35]-0.0883883476483184*fr[20]*nul[33]+0.0883883476483184*fl[20]*nul[33]+0.07654655446197427*fr[12]*nul[33]+0.07654655446197427*fl[12]*nul[33]-0.0883883476483184*fr[13]*nul[32]+0.0883883476483184*fl[13]*nul[32]+0.07654655446197427*fr[5]*nul[32]+0.07654655446197427*fl[5]*nul[32]; 
  incr2[16] = (-0.0883883476483184*fr[2]*nul[39])+0.0883883476483184*fl[2]*nul[39]+0.07654655446197427*fr[0]*nul[39]+0.07654655446197427*fl[0]*nul[39]-0.0883883476483184*fr[6]*nul[35]+0.0883883476483184*fl[6]*nul[35]+0.07654655446197427*fr[1]*nul[35]+0.07654655446197427*fl[1]*nul[35]-0.0883883476483184*fr[8]*nul[33]+0.0883883476483184*fl[8]*nul[33]+0.07654655446197427*fr[3]*nul[33]+0.07654655446197427*fl[3]*nul[33]-0.0883883476483184*fr[16]*nul[32]+0.0883883476483184*fl[16]*nul[32]+0.07654655446197427*fr[7]*nul[32]+0.07654655446197427*fl[7]*nul[32]; 
  incr2[17] = (-0.0883883476483184*fr[19]*nul[39])+0.0883883476483184*fl[19]*nul[39]+0.07654655446197427*fr[11]*nul[39]+0.07654655446197427*fl[11]*nul[39]-0.0883883476483184*fr[26]*nul[35]+0.0883883476483184*fl[26]*nul[35]+0.07654655446197427*fr[18]*nul[35]+0.07654655446197427*fl[18]*nul[35]-0.0883883476483184*fr[10]*nul[33]+0.0883883476483184*fl[10]*nul[33]+0.07654655446197427*fr[4]*nul[33]+0.07654655446197427*fl[4]*nul[33]-0.0883883476483184*fr[17]*nul[32]+0.0883883476483184*fl[17]*nul[32]+0.07654655446197427*fr[9]*nul[32]+0.07654655446197427*fl[9]*nul[32]; 
  incr2[19] = (-0.0883883476483184*fr[17]*nul[39])+0.0883883476483184*fl[17]*nul[39]+0.07654655446197427*fr[9]*nul[39]+0.07654655446197427*fl[9]*nul[39]-0.0883883476483184*fr[10]*nul[35]+0.0883883476483184*fl[10]*nul[35]+0.07654655446197427*fr[4]*nul[35]+0.07654655446197427*fl[4]*nul[35]-0.0883883476483184*fr[26]*nul[33]+0.0883883476483184*fl[26]*nul[33]+0.07654655446197427*fr[18]*nul[33]+0.07654655446197427*fl[18]*nul[33]-0.0883883476483184*fr[19]*nul[32]+0.0883883476483184*fl[19]*nul[32]+0.07654655446197427*fr[11]*nul[32]+0.07654655446197427*fl[11]*nul[32]; 
  incr2[20] = (-0.0883883476483184*fr[22]*nul[39])+0.0883883476483184*fl[22]*nul[39]+0.07654655446197427*fr[14]*nul[39]+0.07654655446197427*fl[14]*nul[39]-0.0883883476483184*fr[27]*nul[35]+0.0883883476483184*fl[27]*nul[35]+0.07654655446197427*fr[21]*nul[35]+0.07654655446197427*fl[21]*nul[35]-0.0883883476483184*fr[13]*nul[33]+0.0883883476483184*fl[13]*nul[33]+0.07654655446197427*fr[5]*nul[33]+0.07654655446197427*fl[5]*nul[33]-0.0883883476483184*fr[20]*nul[32]+0.0883883476483184*fl[20]*nul[32]+0.07654655446197427*fr[12]*nul[32]+0.07654655446197427*fl[12]*nul[32]; 
  incr2[22] = (-0.0883883476483184*fr[20]*nul[39])+0.0883883476483184*fl[20]*nul[39]+0.07654655446197427*fr[12]*nul[39]+0.07654655446197427*fl[12]*nul[39]-0.0883883476483184*fr[13]*nul[35]+0.0883883476483184*fl[13]*nul[35]+0.07654655446197427*fr[5]*nul[35]+0.07654655446197427*fl[5]*nul[35]-0.0883883476483184*fr[27]*nul[33]+0.0883883476483184*fl[27]*nul[33]+0.07654655446197427*fr[21]*nul[33]+0.07654655446197427*fl[21]*nul[33]-0.0883883476483184*fr[22]*nul[32]+0.0883883476483184*fl[22]*nul[32]+0.07654655446197427*fr[14]*nul[32]+0.07654655446197427*fl[14]*nul[32]; 
  incr2[24] = (-0.0883883476483184*fr[31]*nul[39])+0.0883883476483184*fl[31]*nul[39]+0.07654655446197427*fr[29]*nul[39]+0.07654655446197427*fl[29]*nul[39]-0.0883883476483184*fr[30]*nul[35]+0.0883883476483184*fl[30]*nul[35]+0.07654655446197427*fr[25]*nul[35]+0.07654655446197427*fl[25]*nul[35]-0.0883883476483184*fr[28]*nul[33]+0.0883883476483184*fl[28]*nul[33]+0.07654655446197427*fr[23]*nul[33]+0.07654655446197427*fl[23]*nul[33]-0.0883883476483184*fr[24]*nul[32]+0.0883883476483184*fl[24]*nul[32]+0.07654655446197427*fr[15]*nul[32]+0.07654655446197427*fl[15]*nul[32]; 
  incr2[26] = (-0.0883883476483184*fr[10]*nul[39])+0.0883883476483184*fl[10]*nul[39]+0.07654655446197427*fr[4]*nul[39]+0.07654655446197427*fl[4]*nul[39]-0.0883883476483184*fr[17]*nul[35]+0.0883883476483184*fl[17]*nul[35]+0.07654655446197427*fr[9]*nul[35]+0.07654655446197427*fl[9]*nul[35]-0.0883883476483184*fr[19]*nul[33]+0.0883883476483184*fl[19]*nul[33]+0.07654655446197427*fr[11]*nul[33]+0.07654655446197427*fl[11]*nul[33]-0.0883883476483184*fr[26]*nul[32]+0.0883883476483184*fl[26]*nul[32]+0.07654655446197427*fr[18]*nul[32]+0.07654655446197427*fl[18]*nul[32]; 
  incr2[27] = (-0.0883883476483184*fr[13]*nul[39])+0.0883883476483184*fl[13]*nul[39]+0.07654655446197427*fr[5]*nul[39]+0.07654655446197427*fl[5]*nul[39]-0.0883883476483184*fr[20]*nul[35]+0.0883883476483184*fl[20]*nul[35]+0.07654655446197427*fr[12]*nul[35]+0.07654655446197427*fl[12]*nul[35]-0.0883883476483184*fr[22]*nul[33]+0.0883883476483184*fl[22]*nul[33]+0.07654655446197427*fr[14]*nul[33]+0.07654655446197427*fl[14]*nul[33]-0.0883883476483184*fr[27]*nul[32]+0.0883883476483184*fl[27]*nul[32]+0.07654655446197427*fr[21]*nul[32]+0.07654655446197427*fl[21]*nul[32]; 
  incr2[28] = (-0.0883883476483184*fr[30]*nul[39])+0.0883883476483184*fl[30]*nul[39]+0.07654655446197427*fr[25]*nul[39]+0.07654655446197427*fl[25]*nul[39]-0.0883883476483184*fr[31]*nul[35]+0.0883883476483184*fl[31]*nul[35]+0.07654655446197427*fr[29]*nul[35]+0.07654655446197427*fl[29]*nul[35]-0.0883883476483184*fr[24]*nul[33]+0.0883883476483184*fl[24]*nul[33]+0.07654655446197427*fr[15]*nul[33]+0.07654655446197427*fl[15]*nul[33]-0.0883883476483184*fr[28]*nul[32]+0.0883883476483184*fl[28]*nul[32]+0.07654655446197427*fr[23]*nul[32]+0.07654655446197427*fl[23]*nul[32]; 
  incr2[30] = (-0.0883883476483184*fr[28]*nul[39])+0.0883883476483184*fl[28]*nul[39]+0.07654655446197427*fr[23]*nul[39]+0.07654655446197427*fl[23]*nul[39]-0.0883883476483184*fr[24]*nul[35]+0.0883883476483184*fl[24]*nul[35]+0.07654655446197427*fr[15]*nul[35]+0.07654655446197427*fl[15]*nul[35]-0.0883883476483184*fr[31]*nul[33]+0.0883883476483184*fl[31]*nul[33]+0.07654655446197427*fr[29]*nul[33]+0.07654655446197427*fl[29]*nul[33]-0.0883883476483184*fr[30]*nul[32]+0.0883883476483184*fl[30]*nul[32]+0.07654655446197427*fr[25]*nul[32]+0.07654655446197427*fl[25]*nul[32]; 
  incr2[31] = (-0.0883883476483184*fr[24]*nul[39])+0.0883883476483184*fl[24]*nul[39]+0.07654655446197427*fr[15]*nul[39]+0.07654655446197427*fl[15]*nul[39]-0.0883883476483184*fr[28]*nul[35]+0.0883883476483184*fl[28]*nul[35]+0.07654655446197427*fr[23]*nul[35]+0.07654655446197427*fl[23]*nul[35]-0.0883883476483184*fr[30]*nul[33]+0.0883883476483184*fl[30]*nul[33]+0.07654655446197427*fr[25]*nul[33]+0.07654655446197427*fl[25]*nul[33]-0.0883883476483184*fr[31]*nul[32]+0.0883883476483184*fl[31]*nul[32]+0.07654655446197427*fr[29]*nul[32]+0.07654655446197427*fl[29]*nul[32]; 

  outr[0] += incr1[0]*rdxFr; 
  outr[1] += incr1[1]*rdxFr; 
  outr[2] += incr2[2]*rdxFr+incr1[2]*rdxFr; 
  outr[3] += incr1[3]*rdxFr; 
  outr[4] += incr1[4]*rdxFr; 
  outr[5] += incr1[5]*rdxFr; 
  outr[6] += incr2[6]*rdxFr+incr1[6]*rdxFr; 
  outr[7] += incr1[7]*rdxFr; 
  outr[8] += incr2[8]*rdxFr+incr1[8]*rdxFr; 
  outr[9] += incr1[9]*rdxFr; 
  outr[10] += incr2[10]*rdxFr+incr1[10]*rdxFr; 
  outr[11] += incr1[11]*rdxFr; 
  outr[12] += incr1[12]*rdxFr; 
  outr[13] += incr2[13]*rdxFr+incr1[13]*rdxFr; 
  outr[14] += incr1[14]*rdxFr; 
  outr[15] += incr1[15]*rdxFr; 
  outr[16] += incr2[16]*rdxFr+incr1[16]*rdxFr; 
  outr[17] += incr2[17]*rdxFr+incr1[17]*rdxFr; 
  outr[18] += incr1[18]*rdxFr; 
  outr[19] += incr2[19]*rdxFr+incr1[19]*rdxFr; 
  outr[20] += incr2[20]*rdxFr+incr1[20]*rdxFr; 
  outr[21] += incr1[21]*rdxFr; 
  outr[22] += incr2[22]*rdxFr+incr1[22]*rdxFr; 
  outr[23] += incr1[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr+incr1[24]*rdxFr; 
  outr[25] += incr1[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr+incr1[26]*rdxFr; 
  outr[27] += incr2[27]*rdxFr+incr1[27]*rdxFr; 
  outr[28] += incr2[28]*rdxFr+incr1[28]*rdxFr; 
  outr[29] += incr1[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr+incr1[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr+incr1[31]*rdxFr; 

  outl[0] += -1.0*incr1[0]*rdxFl; 
  outl[1] += -1.0*incr1[1]*rdxFl; 
  outl[2] += incr1[2]*rdxFl-1.0*incr2[2]*rdxFl; 
  outl[3] += -1.0*incr1[3]*rdxFl; 
  outl[4] += -1.0*incr1[4]*rdxFl; 
  outl[5] += -1.0*incr1[5]*rdxFl; 
  outl[6] += incr1[6]*rdxFl-1.0*incr2[6]*rdxFl; 
  outl[7] += -1.0*incr1[7]*rdxFl; 
  outl[8] += incr1[8]*rdxFl-1.0*incr2[8]*rdxFl; 
  outl[9] += -1.0*incr1[9]*rdxFl; 
  outl[10] += incr1[10]*rdxFl-1.0*incr2[10]*rdxFl; 
  outl[11] += -1.0*incr1[11]*rdxFl; 
  outl[12] += -1.0*incr1[12]*rdxFl; 
  outl[13] += incr1[13]*rdxFl-1.0*incr2[13]*rdxFl; 
  outl[14] += -1.0*incr1[14]*rdxFl; 
  outl[15] += -1.0*incr1[15]*rdxFl; 
  outl[16] += incr1[16]*rdxFl-1.0*incr2[16]*rdxFl; 
  outl[17] += incr1[17]*rdxFl-1.0*incr2[17]*rdxFl; 
  outl[18] += -1.0*incr1[18]*rdxFl; 
  outl[19] += incr1[19]*rdxFl-1.0*incr2[19]*rdxFl; 
  outl[20] += incr1[20]*rdxFl-1.0*incr2[20]*rdxFl; 
  outl[21] += -1.0*incr1[21]*rdxFl; 
  outl[22] += incr1[22]*rdxFl-1.0*incr2[22]*rdxFl; 
  outl[23] += -1.0*incr1[23]*rdxFl; 
  outl[24] += incr1[24]*rdxFl-1.0*incr2[24]*rdxFl; 
  outl[25] += -1.0*incr1[25]*rdxFl; 
  outl[26] += incr1[26]*rdxFl-1.0*incr2[26]*rdxFl; 
  outl[27] += incr1[27]*rdxFl-1.0*incr2[27]*rdxFl; 
  outl[28] += incr1[28]*rdxFl-1.0*incr2[28]*rdxFl; 
  outl[29] += -1.0*incr1[29]*rdxFl; 
  outl[30] += incr1[30]*rdxFl-1.0*incr2[30]*rdxFl; 
  outl[31] += incr1[31]*rdxFl-1.0*incr2[31]*rdxFl; 

} 
void ConstDiffusionVarCoeffSurf5xSerP1_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[160]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[2]*dxl[2]); 
  double rdxFr = 4.0/(dxr[2]*dxr[2]); 

  double incr1[32]; 
  double incr2[32]; 

  incr1[0] = 0.1657281518405969*fr[7]*nul[71]+0.1657281518405969*fl[7]*nul[71]-0.172229747539442*fr[1]*nul[71]+0.172229747539442*fl[1]*nul[71]+0.1657281518405969*fr[3]*nul[67]+0.1657281518405969*fl[3]*nul[67]-0.172229747539442*fr[0]*nul[67]+0.172229747539442*fl[0]*nul[67]+0.09568319307746781*fr[7]*nul[65]+0.09568319307746781*fl[7]*nul[65]-0.09943689110435816*fr[1]*nul[65]+0.09943689110435816*fl[1]*nul[65]+0.09568319307746781*fr[3]*nul[64]+0.09568319307746781*fl[3]*nul[64]-0.09943689110435816*fr[0]*nul[64]+0.09943689110435816*fl[0]*nul[64]; 
  incr1[1] = 0.1657281518405969*fr[3]*nul[71]+0.1657281518405969*fl[3]*nul[71]-0.172229747539442*fr[0]*nul[71]+0.172229747539442*fl[0]*nul[71]+0.1657281518405969*fr[7]*nul[67]+0.1657281518405969*fl[7]*nul[67]-0.172229747539442*fr[1]*nul[67]+0.172229747539442*fl[1]*nul[67]+0.09568319307746781*fr[3]*nul[65]+0.09568319307746781*fl[3]*nul[65]-0.09943689110435816*fr[0]*nul[65]+0.09943689110435816*fl[0]*nul[65]+0.09568319307746781*fr[7]*nul[64]+0.09568319307746781*fl[7]*nul[64]-0.09943689110435816*fr[1]*nul[64]+0.09943689110435816*fl[1]*nul[64]; 
  incr1[2] = 0.1657281518405969*fr[16]*nul[71]+0.1657281518405969*fl[16]*nul[71]-0.172229747539442*fr[6]*nul[71]+0.172229747539442*fl[6]*nul[71]+0.1657281518405969*fr[8]*nul[67]+0.1657281518405969*fl[8]*nul[67]-0.172229747539442*fr[2]*nul[67]+0.172229747539442*fl[2]*nul[67]+0.09568319307746781*fr[16]*nul[65]+0.09568319307746781*fl[16]*nul[65]-0.09943689110435816*fr[6]*nul[65]+0.09943689110435816*fl[6]*nul[65]+0.09568319307746781*fr[8]*nul[64]+0.09568319307746781*fl[8]*nul[64]-0.09943689110435816*fr[2]*nul[64]+0.09943689110435816*fl[2]*nul[64]; 
  incr1[3] = (-0.2870495792324034*fr[7]*nul[71])-0.2870495792324034*fl[7]*nul[71]+0.2983106733130745*fr[1]*nul[71]-0.2983106733130745*fl[1]*nul[71]-0.2870495792324034*fr[3]*nul[67]-0.2870495792324034*fl[3]*nul[67]+0.2983106733130745*fr[0]*nul[67]-0.2983106733130745*fl[0]*nul[67]-0.1657281518405969*fr[7]*nul[65]-0.1657281518405969*fl[7]*nul[65]+0.172229747539442*fr[1]*nul[65]-0.172229747539442*fl[1]*nul[65]-0.1657281518405969*fr[3]*nul[64]-0.1657281518405969*fl[3]*nul[64]+0.172229747539442*fr[0]*nul[64]-0.172229747539442*fl[0]*nul[64]; 
  incr1[4] = 0.1657281518405969*fr[18]*nul[71]+0.1657281518405969*fl[18]*nul[71]-0.172229747539442*fr[9]*nul[71]+0.172229747539442*fl[9]*nul[71]+0.1657281518405969*fr[11]*nul[67]+0.1657281518405969*fl[11]*nul[67]-0.172229747539442*fr[4]*nul[67]+0.172229747539442*fl[4]*nul[67]+0.09568319307746781*fr[18]*nul[65]+0.09568319307746781*fl[18]*nul[65]-0.09943689110435816*fr[9]*nul[65]+0.09943689110435816*fl[9]*nul[65]+0.09568319307746781*fr[11]*nul[64]+0.09568319307746781*fl[11]*nul[64]-0.09943689110435816*fr[4]*nul[64]+0.09943689110435816*fl[4]*nul[64]; 
  incr1[5] = 0.1657281518405969*fr[21]*nul[71]+0.1657281518405969*fl[21]*nul[71]-0.172229747539442*fr[12]*nul[71]+0.172229747539442*fl[12]*nul[71]+0.1657281518405969*fr[14]*nul[67]+0.1657281518405969*fl[14]*nul[67]-0.172229747539442*fr[5]*nul[67]+0.172229747539442*fl[5]*nul[67]+0.09568319307746781*fr[21]*nul[65]+0.09568319307746781*fl[21]*nul[65]-0.09943689110435816*fr[12]*nul[65]+0.09943689110435816*fl[12]*nul[65]+0.09568319307746781*fr[14]*nul[64]+0.09568319307746781*fl[14]*nul[64]-0.09943689110435816*fr[5]*nul[64]+0.09943689110435816*fl[5]*nul[64]; 
  incr1[6] = 0.1657281518405969*fr[8]*nul[71]+0.1657281518405969*fl[8]*nul[71]-0.172229747539442*fr[2]*nul[71]+0.172229747539442*fl[2]*nul[71]+0.1657281518405969*fr[16]*nul[67]+0.1657281518405969*fl[16]*nul[67]-0.172229747539442*fr[6]*nul[67]+0.172229747539442*fl[6]*nul[67]+0.09568319307746781*fr[8]*nul[65]+0.09568319307746781*fl[8]*nul[65]-0.09943689110435816*fr[2]*nul[65]+0.09943689110435816*fl[2]*nul[65]+0.09568319307746781*fr[16]*nul[64]+0.09568319307746781*fl[16]*nul[64]-0.09943689110435816*fr[6]*nul[64]+0.09943689110435816*fl[6]*nul[64]; 
  incr1[7] = (-0.2870495792324034*fr[3]*nul[71])-0.2870495792324034*fl[3]*nul[71]+0.2983106733130745*fr[0]*nul[71]-0.2983106733130745*fl[0]*nul[71]-0.2870495792324034*fr[7]*nul[67]-0.2870495792324034*fl[7]*nul[67]+0.2983106733130745*fr[1]*nul[67]-0.2983106733130745*fl[1]*nul[67]-0.1657281518405969*fr[3]*nul[65]-0.1657281518405969*fl[3]*nul[65]+0.172229747539442*fr[0]*nul[65]-0.172229747539442*fl[0]*nul[65]-0.1657281518405969*fr[7]*nul[64]-0.1657281518405969*fl[7]*nul[64]+0.172229747539442*fr[1]*nul[64]-0.172229747539442*fl[1]*nul[64]; 
  incr1[8] = (-0.2870495792324034*fr[16]*nul[71])-0.2870495792324034*fl[16]*nul[71]+0.2983106733130745*fr[6]*nul[71]-0.2983106733130745*fl[6]*nul[71]-0.2870495792324034*fr[8]*nul[67]-0.2870495792324034*fl[8]*nul[67]+0.2983106733130745*fr[2]*nul[67]-0.2983106733130745*fl[2]*nul[67]-0.1657281518405969*fr[16]*nul[65]-0.1657281518405969*fl[16]*nul[65]+0.172229747539442*fr[6]*nul[65]-0.172229747539442*fl[6]*nul[65]-0.1657281518405969*fr[8]*nul[64]-0.1657281518405969*fl[8]*nul[64]+0.172229747539442*fr[2]*nul[64]-0.172229747539442*fl[2]*nul[64]; 
  incr1[9] = 0.1657281518405969*fr[11]*nul[71]+0.1657281518405969*fl[11]*nul[71]-0.172229747539442*fr[4]*nul[71]+0.172229747539442*fl[4]*nul[71]+0.1657281518405969*fr[18]*nul[67]+0.1657281518405969*fl[18]*nul[67]-0.172229747539442*fr[9]*nul[67]+0.172229747539442*fl[9]*nul[67]+0.09568319307746781*fr[11]*nul[65]+0.09568319307746781*fl[11]*nul[65]-0.09943689110435816*fr[4]*nul[65]+0.09943689110435816*fl[4]*nul[65]+0.09568319307746781*fr[18]*nul[64]+0.09568319307746781*fl[18]*nul[64]-0.09943689110435816*fr[9]*nul[64]+0.09943689110435816*fl[9]*nul[64]; 
  incr1[10] = 0.1657281518405969*fr[26]*nul[71]+0.1657281518405969*fl[26]*nul[71]-0.172229747539442*fr[17]*nul[71]+0.172229747539442*fl[17]*nul[71]+0.1657281518405969*fr[19]*nul[67]+0.1657281518405969*fl[19]*nul[67]-0.172229747539442*fr[10]*nul[67]+0.172229747539442*fl[10]*nul[67]+0.09568319307746781*fr[26]*nul[65]+0.09568319307746781*fl[26]*nul[65]-0.09943689110435816*fr[17]*nul[65]+0.09943689110435816*fl[17]*nul[65]+0.09568319307746781*fr[19]*nul[64]+0.09568319307746781*fl[19]*nul[64]-0.09943689110435816*fr[10]*nul[64]+0.09943689110435816*fl[10]*nul[64]; 
  incr1[11] = (-0.2870495792324034*fr[18]*nul[71])-0.2870495792324034*fl[18]*nul[71]+0.2983106733130745*fr[9]*nul[71]-0.2983106733130745*fl[9]*nul[71]-0.2870495792324034*fr[11]*nul[67]-0.2870495792324034*fl[11]*nul[67]+0.2983106733130745*fr[4]*nul[67]-0.2983106733130745*fl[4]*nul[67]-0.1657281518405969*fr[18]*nul[65]-0.1657281518405969*fl[18]*nul[65]+0.172229747539442*fr[9]*nul[65]-0.172229747539442*fl[9]*nul[65]-0.1657281518405969*fr[11]*nul[64]-0.1657281518405969*fl[11]*nul[64]+0.172229747539442*fr[4]*nul[64]-0.172229747539442*fl[4]*nul[64]; 
  incr1[12] = 0.1657281518405969*fr[14]*nul[71]+0.1657281518405969*fl[14]*nul[71]-0.172229747539442*fr[5]*nul[71]+0.172229747539442*fl[5]*nul[71]+0.1657281518405969*fr[21]*nul[67]+0.1657281518405969*fl[21]*nul[67]-0.172229747539442*fr[12]*nul[67]+0.172229747539442*fl[12]*nul[67]+0.09568319307746781*fr[14]*nul[65]+0.09568319307746781*fl[14]*nul[65]-0.09943689110435816*fr[5]*nul[65]+0.09943689110435816*fl[5]*nul[65]+0.09568319307746781*fr[21]*nul[64]+0.09568319307746781*fl[21]*nul[64]-0.09943689110435816*fr[12]*nul[64]+0.09943689110435816*fl[12]*nul[64]; 
  incr1[13] = 0.1657281518405969*fr[27]*nul[71]+0.1657281518405969*fl[27]*nul[71]-0.172229747539442*fr[20]*nul[71]+0.172229747539442*fl[20]*nul[71]+0.1657281518405969*fr[22]*nul[67]+0.1657281518405969*fl[22]*nul[67]-0.172229747539442*fr[13]*nul[67]+0.172229747539442*fl[13]*nul[67]+0.09568319307746781*fr[27]*nul[65]+0.09568319307746781*fl[27]*nul[65]-0.09943689110435816*fr[20]*nul[65]+0.09943689110435816*fl[20]*nul[65]+0.09568319307746781*fr[22]*nul[64]+0.09568319307746781*fl[22]*nul[64]-0.09943689110435816*fr[13]*nul[64]+0.09943689110435816*fl[13]*nul[64]; 
  incr1[14] = (-0.2870495792324034*fr[21]*nul[71])-0.2870495792324034*fl[21]*nul[71]+0.2983106733130745*fr[12]*nul[71]-0.2983106733130745*fl[12]*nul[71]-0.2870495792324034*fr[14]*nul[67]-0.2870495792324034*fl[14]*nul[67]+0.2983106733130745*fr[5]*nul[67]-0.2983106733130745*fl[5]*nul[67]-0.1657281518405969*fr[21]*nul[65]-0.1657281518405969*fl[21]*nul[65]+0.172229747539442*fr[12]*nul[65]-0.172229747539442*fl[12]*nul[65]-0.1657281518405969*fr[14]*nul[64]-0.1657281518405969*fl[14]*nul[64]+0.172229747539442*fr[5]*nul[64]-0.172229747539442*fl[5]*nul[64]; 
  incr1[15] = 0.1657281518405969*fr[29]*nul[71]+0.1657281518405969*fl[29]*nul[71]-0.172229747539442*fr[23]*nul[71]+0.172229747539442*fl[23]*nul[71]+0.1657281518405969*fr[25]*nul[67]+0.1657281518405969*fl[25]*nul[67]-0.172229747539442*fr[15]*nul[67]+0.172229747539442*fl[15]*nul[67]+0.09568319307746781*fr[29]*nul[65]+0.09568319307746781*fl[29]*nul[65]-0.09943689110435816*fr[23]*nul[65]+0.09943689110435816*fl[23]*nul[65]+0.09568319307746781*fr[25]*nul[64]+0.09568319307746781*fl[25]*nul[64]-0.09943689110435816*fr[15]*nul[64]+0.09943689110435816*fl[15]*nul[64]; 
  incr1[16] = (-0.2870495792324034*fr[8]*nul[71])-0.2870495792324034*fl[8]*nul[71]+0.2983106733130745*fr[2]*nul[71]-0.2983106733130745*fl[2]*nul[71]-0.2870495792324034*fr[16]*nul[67]-0.2870495792324034*fl[16]*nul[67]+0.2983106733130745*fr[6]*nul[67]-0.2983106733130745*fl[6]*nul[67]-0.1657281518405969*fr[8]*nul[65]-0.1657281518405969*fl[8]*nul[65]+0.172229747539442*fr[2]*nul[65]-0.172229747539442*fl[2]*nul[65]-0.1657281518405969*fr[16]*nul[64]-0.1657281518405969*fl[16]*nul[64]+0.172229747539442*fr[6]*nul[64]-0.172229747539442*fl[6]*nul[64]; 
  incr1[17] = 0.1657281518405969*fr[19]*nul[71]+0.1657281518405969*fl[19]*nul[71]-0.172229747539442*fr[10]*nul[71]+0.172229747539442*fl[10]*nul[71]+0.1657281518405969*fr[26]*nul[67]+0.1657281518405969*fl[26]*nul[67]-0.172229747539442*fr[17]*nul[67]+0.172229747539442*fl[17]*nul[67]+0.09568319307746781*fr[19]*nul[65]+0.09568319307746781*fl[19]*nul[65]-0.09943689110435816*fr[10]*nul[65]+0.09943689110435816*fl[10]*nul[65]+0.09568319307746781*fr[26]*nul[64]+0.09568319307746781*fl[26]*nul[64]-0.09943689110435816*fr[17]*nul[64]+0.09943689110435816*fl[17]*nul[64]; 
  incr1[18] = (-0.2870495792324034*fr[11]*nul[71])-0.2870495792324034*fl[11]*nul[71]+0.2983106733130745*fr[4]*nul[71]-0.2983106733130745*fl[4]*nul[71]-0.2870495792324034*fr[18]*nul[67]-0.2870495792324034*fl[18]*nul[67]+0.2983106733130745*fr[9]*nul[67]-0.2983106733130745*fl[9]*nul[67]-0.1657281518405969*fr[11]*nul[65]-0.1657281518405969*fl[11]*nul[65]+0.172229747539442*fr[4]*nul[65]-0.172229747539442*fl[4]*nul[65]-0.1657281518405969*fr[18]*nul[64]-0.1657281518405969*fl[18]*nul[64]+0.172229747539442*fr[9]*nul[64]-0.172229747539442*fl[9]*nul[64]; 
  incr1[19] = (-0.2870495792324034*fr[26]*nul[71])-0.2870495792324034*fl[26]*nul[71]+0.2983106733130745*fr[17]*nul[71]-0.2983106733130745*fl[17]*nul[71]-0.2870495792324034*fr[19]*nul[67]-0.2870495792324034*fl[19]*nul[67]+0.2983106733130745*fr[10]*nul[67]-0.2983106733130745*fl[10]*nul[67]-0.1657281518405969*fr[26]*nul[65]-0.1657281518405969*fl[26]*nul[65]+0.172229747539442*fr[17]*nul[65]-0.172229747539442*fl[17]*nul[65]-0.1657281518405969*fr[19]*nul[64]-0.1657281518405969*fl[19]*nul[64]+0.172229747539442*fr[10]*nul[64]-0.172229747539442*fl[10]*nul[64]; 
  incr1[20] = 0.1657281518405969*fr[22]*nul[71]+0.1657281518405969*fl[22]*nul[71]-0.172229747539442*fr[13]*nul[71]+0.172229747539442*fl[13]*nul[71]+0.1657281518405969*fr[27]*nul[67]+0.1657281518405969*fl[27]*nul[67]-0.172229747539442*fr[20]*nul[67]+0.172229747539442*fl[20]*nul[67]+0.09568319307746781*fr[22]*nul[65]+0.09568319307746781*fl[22]*nul[65]-0.09943689110435816*fr[13]*nul[65]+0.09943689110435816*fl[13]*nul[65]+0.09568319307746781*fr[27]*nul[64]+0.09568319307746781*fl[27]*nul[64]-0.09943689110435816*fr[20]*nul[64]+0.09943689110435816*fl[20]*nul[64]; 
  incr1[21] = (-0.2870495792324034*fr[14]*nul[71])-0.2870495792324034*fl[14]*nul[71]+0.2983106733130745*fr[5]*nul[71]-0.2983106733130745*fl[5]*nul[71]-0.2870495792324034*fr[21]*nul[67]-0.2870495792324034*fl[21]*nul[67]+0.2983106733130745*fr[12]*nul[67]-0.2983106733130745*fl[12]*nul[67]-0.1657281518405969*fr[14]*nul[65]-0.1657281518405969*fl[14]*nul[65]+0.172229747539442*fr[5]*nul[65]-0.172229747539442*fl[5]*nul[65]-0.1657281518405969*fr[21]*nul[64]-0.1657281518405969*fl[21]*nul[64]+0.172229747539442*fr[12]*nul[64]-0.172229747539442*fl[12]*nul[64]; 
  incr1[22] = (-0.2870495792324034*fr[27]*nul[71])-0.2870495792324034*fl[27]*nul[71]+0.2983106733130745*fr[20]*nul[71]-0.2983106733130745*fl[20]*nul[71]-0.2870495792324034*fr[22]*nul[67]-0.2870495792324034*fl[22]*nul[67]+0.2983106733130745*fr[13]*nul[67]-0.2983106733130745*fl[13]*nul[67]-0.1657281518405969*fr[27]*nul[65]-0.1657281518405969*fl[27]*nul[65]+0.172229747539442*fr[20]*nul[65]-0.172229747539442*fl[20]*nul[65]-0.1657281518405969*fr[22]*nul[64]-0.1657281518405969*fl[22]*nul[64]+0.172229747539442*fr[13]*nul[64]-0.172229747539442*fl[13]*nul[64]; 
  incr1[23] = 0.1657281518405969*fr[25]*nul[71]+0.1657281518405969*fl[25]*nul[71]-0.172229747539442*fr[15]*nul[71]+0.172229747539442*fl[15]*nul[71]+0.1657281518405969*fr[29]*nul[67]+0.1657281518405969*fl[29]*nul[67]-0.172229747539442*fr[23]*nul[67]+0.172229747539442*fl[23]*nul[67]+0.09568319307746781*fr[25]*nul[65]+0.09568319307746781*fl[25]*nul[65]-0.09943689110435816*fr[15]*nul[65]+0.09943689110435816*fl[15]*nul[65]+0.09568319307746781*fr[29]*nul[64]+0.09568319307746781*fl[29]*nul[64]-0.09943689110435816*fr[23]*nul[64]+0.09943689110435816*fl[23]*nul[64]; 
  incr1[24] = 0.1657281518405969*fr[31]*nul[71]+0.1657281518405969*fl[31]*nul[71]-0.172229747539442*fr[28]*nul[71]+0.172229747539442*fl[28]*nul[71]+0.1657281518405969*fr[30]*nul[67]+0.1657281518405969*fl[30]*nul[67]-0.172229747539442*fr[24]*nul[67]+0.172229747539442*fl[24]*nul[67]+0.09568319307746781*fr[31]*nul[65]+0.09568319307746781*fl[31]*nul[65]-0.09943689110435816*fr[28]*nul[65]+0.09943689110435816*fl[28]*nul[65]+0.09568319307746781*fr[30]*nul[64]+0.09568319307746781*fl[30]*nul[64]-0.09943689110435816*fr[24]*nul[64]+0.09943689110435816*fl[24]*nul[64]; 
  incr1[25] = (-0.2870495792324034*fr[29]*nul[71])-0.2870495792324034*fl[29]*nul[71]+0.2983106733130745*fr[23]*nul[71]-0.2983106733130745*fl[23]*nul[71]-0.2870495792324034*fr[25]*nul[67]-0.2870495792324034*fl[25]*nul[67]+0.2983106733130745*fr[15]*nul[67]-0.2983106733130745*fl[15]*nul[67]-0.1657281518405969*fr[29]*nul[65]-0.1657281518405969*fl[29]*nul[65]+0.172229747539442*fr[23]*nul[65]-0.172229747539442*fl[23]*nul[65]-0.1657281518405969*fr[25]*nul[64]-0.1657281518405969*fl[25]*nul[64]+0.172229747539442*fr[15]*nul[64]-0.172229747539442*fl[15]*nul[64]; 
  incr1[26] = (-0.2870495792324034*fr[19]*nul[71])-0.2870495792324034*fl[19]*nul[71]+0.2983106733130745*fr[10]*nul[71]-0.2983106733130745*fl[10]*nul[71]-0.2870495792324034*fr[26]*nul[67]-0.2870495792324034*fl[26]*nul[67]+0.2983106733130745*fr[17]*nul[67]-0.2983106733130745*fl[17]*nul[67]-0.1657281518405969*fr[19]*nul[65]-0.1657281518405969*fl[19]*nul[65]+0.172229747539442*fr[10]*nul[65]-0.172229747539442*fl[10]*nul[65]-0.1657281518405969*fr[26]*nul[64]-0.1657281518405969*fl[26]*nul[64]+0.172229747539442*fr[17]*nul[64]-0.172229747539442*fl[17]*nul[64]; 
  incr1[27] = (-0.2870495792324034*fr[22]*nul[71])-0.2870495792324034*fl[22]*nul[71]+0.2983106733130745*fr[13]*nul[71]-0.2983106733130745*fl[13]*nul[71]-0.2870495792324034*fr[27]*nul[67]-0.2870495792324034*fl[27]*nul[67]+0.2983106733130745*fr[20]*nul[67]-0.2983106733130745*fl[20]*nul[67]-0.1657281518405969*fr[22]*nul[65]-0.1657281518405969*fl[22]*nul[65]+0.172229747539442*fr[13]*nul[65]-0.172229747539442*fl[13]*nul[65]-0.1657281518405969*fr[27]*nul[64]-0.1657281518405969*fl[27]*nul[64]+0.172229747539442*fr[20]*nul[64]-0.172229747539442*fl[20]*nul[64]; 
  incr1[28] = 0.1657281518405969*fr[30]*nul[71]+0.1657281518405969*fl[30]*nul[71]-0.172229747539442*fr[24]*nul[71]+0.172229747539442*fl[24]*nul[71]+0.1657281518405969*fr[31]*nul[67]+0.1657281518405969*fl[31]*nul[67]-0.172229747539442*fr[28]*nul[67]+0.172229747539442*fl[28]*nul[67]+0.09568319307746781*fr[30]*nul[65]+0.09568319307746781*fl[30]*nul[65]-0.09943689110435816*fr[24]*nul[65]+0.09943689110435816*fl[24]*nul[65]+0.09568319307746781*fr[31]*nul[64]+0.09568319307746781*fl[31]*nul[64]-0.09943689110435816*fr[28]*nul[64]+0.09943689110435816*fl[28]*nul[64]; 
  incr1[29] = (-0.2870495792324034*fr[25]*nul[71])-0.2870495792324034*fl[25]*nul[71]+0.2983106733130745*fr[15]*nul[71]-0.2983106733130745*fl[15]*nul[71]-0.2870495792324034*fr[29]*nul[67]-0.2870495792324034*fl[29]*nul[67]+0.2983106733130745*fr[23]*nul[67]-0.2983106733130745*fl[23]*nul[67]-0.1657281518405969*fr[25]*nul[65]-0.1657281518405969*fl[25]*nul[65]+0.172229747539442*fr[15]*nul[65]-0.172229747539442*fl[15]*nul[65]-0.1657281518405969*fr[29]*nul[64]-0.1657281518405969*fl[29]*nul[64]+0.172229747539442*fr[23]*nul[64]-0.172229747539442*fl[23]*nul[64]; 
  incr1[30] = (-0.2870495792324034*fr[31]*nul[71])-0.2870495792324034*fl[31]*nul[71]+0.2983106733130745*fr[28]*nul[71]-0.2983106733130745*fl[28]*nul[71]-0.2870495792324034*fr[30]*nul[67]-0.2870495792324034*fl[30]*nul[67]+0.2983106733130745*fr[24]*nul[67]-0.2983106733130745*fl[24]*nul[67]-0.1657281518405969*fr[31]*nul[65]-0.1657281518405969*fl[31]*nul[65]+0.172229747539442*fr[28]*nul[65]-0.172229747539442*fl[28]*nul[65]-0.1657281518405969*fr[30]*nul[64]-0.1657281518405969*fl[30]*nul[64]+0.172229747539442*fr[24]*nul[64]-0.172229747539442*fl[24]*nul[64]; 
  incr1[31] = (-0.2870495792324034*fr[30]*nul[71])-0.2870495792324034*fl[30]*nul[71]+0.2983106733130745*fr[24]*nul[71]-0.2983106733130745*fl[24]*nul[71]-0.2870495792324034*fr[31]*nul[67]-0.2870495792324034*fl[31]*nul[67]+0.2983106733130745*fr[28]*nul[67]-0.2983106733130745*fl[28]*nul[67]-0.1657281518405969*fr[30]*nul[65]-0.1657281518405969*fl[30]*nul[65]+0.172229747539442*fr[24]*nul[65]-0.172229747539442*fl[24]*nul[65]-0.1657281518405969*fr[31]*nul[64]-0.1657281518405969*fl[31]*nul[64]+0.172229747539442*fr[28]*nul[64]-0.172229747539442*fl[28]*nul[64]; 

  incr2[3] = (-0.1530931089239486*fr[7]*nul[71])+0.1530931089239486*fl[7]*nul[71]+0.1325825214724776*fr[1]*nul[71]+0.1325825214724776*fl[1]*nul[71]-0.1530931089239486*fr[3]*nul[67]+0.1530931089239486*fl[3]*nul[67]+0.1325825214724776*fr[0]*nul[67]+0.1325825214724776*fl[0]*nul[67]-0.0883883476483184*fr[7]*nul[65]+0.0883883476483184*fl[7]*nul[65]+0.07654655446197427*fr[1]*nul[65]+0.07654655446197427*fl[1]*nul[65]-0.0883883476483184*fr[3]*nul[64]+0.0883883476483184*fl[3]*nul[64]+0.07654655446197427*fr[0]*nul[64]+0.07654655446197427*fl[0]*nul[64]; 
  incr2[7] = (-0.1530931089239486*fr[3]*nul[71])+0.1530931089239486*fl[3]*nul[71]+0.1325825214724776*fr[0]*nul[71]+0.1325825214724776*fl[0]*nul[71]-0.1530931089239486*fr[7]*nul[67]+0.1530931089239486*fl[7]*nul[67]+0.1325825214724776*fr[1]*nul[67]+0.1325825214724776*fl[1]*nul[67]-0.0883883476483184*fr[3]*nul[65]+0.0883883476483184*fl[3]*nul[65]+0.07654655446197427*fr[0]*nul[65]+0.07654655446197427*fl[0]*nul[65]-0.0883883476483184*fr[7]*nul[64]+0.0883883476483184*fl[7]*nul[64]+0.07654655446197427*fr[1]*nul[64]+0.07654655446197427*fl[1]*nul[64]; 
  incr2[8] = (-0.1530931089239486*fr[16]*nul[71])+0.1530931089239486*fl[16]*nul[71]+0.1325825214724776*fr[6]*nul[71]+0.1325825214724776*fl[6]*nul[71]-0.1530931089239486*fr[8]*nul[67]+0.1530931089239486*fl[8]*nul[67]+0.1325825214724776*fr[2]*nul[67]+0.1325825214724776*fl[2]*nul[67]-0.0883883476483184*fr[16]*nul[65]+0.0883883476483184*fl[16]*nul[65]+0.07654655446197427*fr[6]*nul[65]+0.07654655446197427*fl[6]*nul[65]-0.0883883476483184*fr[8]*nul[64]+0.0883883476483184*fl[8]*nul[64]+0.07654655446197427*fr[2]*nul[64]+0.07654655446197427*fl[2]*nul[64]; 
  incr2[11] = (-0.1530931089239486*fr[18]*nul[71])+0.1530931089239486*fl[18]*nul[71]+0.1325825214724776*fr[9]*nul[71]+0.1325825214724776*fl[9]*nul[71]-0.1530931089239486*fr[11]*nul[67]+0.1530931089239486*fl[11]*nul[67]+0.1325825214724776*fr[4]*nul[67]+0.1325825214724776*fl[4]*nul[67]-0.0883883476483184*fr[18]*nul[65]+0.0883883476483184*fl[18]*nul[65]+0.07654655446197427*fr[9]*nul[65]+0.07654655446197427*fl[9]*nul[65]-0.0883883476483184*fr[11]*nul[64]+0.0883883476483184*fl[11]*nul[64]+0.07654655446197427*fr[4]*nul[64]+0.07654655446197427*fl[4]*nul[64]; 
  incr2[14] = (-0.1530931089239486*fr[21]*nul[71])+0.1530931089239486*fl[21]*nul[71]+0.1325825214724776*fr[12]*nul[71]+0.1325825214724776*fl[12]*nul[71]-0.1530931089239486*fr[14]*nul[67]+0.1530931089239486*fl[14]*nul[67]+0.1325825214724776*fr[5]*nul[67]+0.1325825214724776*fl[5]*nul[67]-0.0883883476483184*fr[21]*nul[65]+0.0883883476483184*fl[21]*nul[65]+0.07654655446197427*fr[12]*nul[65]+0.07654655446197427*fl[12]*nul[65]-0.0883883476483184*fr[14]*nul[64]+0.0883883476483184*fl[14]*nul[64]+0.07654655446197427*fr[5]*nul[64]+0.07654655446197427*fl[5]*nul[64]; 
  incr2[16] = (-0.1530931089239486*fr[8]*nul[71])+0.1530931089239486*fl[8]*nul[71]+0.1325825214724776*fr[2]*nul[71]+0.1325825214724776*fl[2]*nul[71]-0.1530931089239486*fr[16]*nul[67]+0.1530931089239486*fl[16]*nul[67]+0.1325825214724776*fr[6]*nul[67]+0.1325825214724776*fl[6]*nul[67]-0.0883883476483184*fr[8]*nul[65]+0.0883883476483184*fl[8]*nul[65]+0.07654655446197427*fr[2]*nul[65]+0.07654655446197427*fl[2]*nul[65]-0.0883883476483184*fr[16]*nul[64]+0.0883883476483184*fl[16]*nul[64]+0.07654655446197427*fr[6]*nul[64]+0.07654655446197427*fl[6]*nul[64]; 
  incr2[18] = (-0.1530931089239486*fr[11]*nul[71])+0.1530931089239486*fl[11]*nul[71]+0.1325825214724776*fr[4]*nul[71]+0.1325825214724776*fl[4]*nul[71]-0.1530931089239486*fr[18]*nul[67]+0.1530931089239486*fl[18]*nul[67]+0.1325825214724776*fr[9]*nul[67]+0.1325825214724776*fl[9]*nul[67]-0.0883883476483184*fr[11]*nul[65]+0.0883883476483184*fl[11]*nul[65]+0.07654655446197427*fr[4]*nul[65]+0.07654655446197427*fl[4]*nul[65]-0.0883883476483184*fr[18]*nul[64]+0.0883883476483184*fl[18]*nul[64]+0.07654655446197427*fr[9]*nul[64]+0.07654655446197427*fl[9]*nul[64]; 
  incr2[19] = (-0.1530931089239486*fr[26]*nul[71])+0.1530931089239486*fl[26]*nul[71]+0.1325825214724776*fr[17]*nul[71]+0.1325825214724776*fl[17]*nul[71]-0.1530931089239486*fr[19]*nul[67]+0.1530931089239486*fl[19]*nul[67]+0.1325825214724776*fr[10]*nul[67]+0.1325825214724776*fl[10]*nul[67]-0.0883883476483184*fr[26]*nul[65]+0.0883883476483184*fl[26]*nul[65]+0.07654655446197427*fr[17]*nul[65]+0.07654655446197427*fl[17]*nul[65]-0.0883883476483184*fr[19]*nul[64]+0.0883883476483184*fl[19]*nul[64]+0.07654655446197427*fr[10]*nul[64]+0.07654655446197427*fl[10]*nul[64]; 
  incr2[21] = (-0.1530931089239486*fr[14]*nul[71])+0.1530931089239486*fl[14]*nul[71]+0.1325825214724776*fr[5]*nul[71]+0.1325825214724776*fl[5]*nul[71]-0.1530931089239486*fr[21]*nul[67]+0.1530931089239486*fl[21]*nul[67]+0.1325825214724776*fr[12]*nul[67]+0.1325825214724776*fl[12]*nul[67]-0.0883883476483184*fr[14]*nul[65]+0.0883883476483184*fl[14]*nul[65]+0.07654655446197427*fr[5]*nul[65]+0.07654655446197427*fl[5]*nul[65]-0.0883883476483184*fr[21]*nul[64]+0.0883883476483184*fl[21]*nul[64]+0.07654655446197427*fr[12]*nul[64]+0.07654655446197427*fl[12]*nul[64]; 
  incr2[22] = (-0.1530931089239486*fr[27]*nul[71])+0.1530931089239486*fl[27]*nul[71]+0.1325825214724776*fr[20]*nul[71]+0.1325825214724776*fl[20]*nul[71]-0.1530931089239486*fr[22]*nul[67]+0.1530931089239486*fl[22]*nul[67]+0.1325825214724776*fr[13]*nul[67]+0.1325825214724776*fl[13]*nul[67]-0.0883883476483184*fr[27]*nul[65]+0.0883883476483184*fl[27]*nul[65]+0.07654655446197427*fr[20]*nul[65]+0.07654655446197427*fl[20]*nul[65]-0.0883883476483184*fr[22]*nul[64]+0.0883883476483184*fl[22]*nul[64]+0.07654655446197427*fr[13]*nul[64]+0.07654655446197427*fl[13]*nul[64]; 
  incr2[25] = (-0.1530931089239486*fr[29]*nul[71])+0.1530931089239486*fl[29]*nul[71]+0.1325825214724776*fr[23]*nul[71]+0.1325825214724776*fl[23]*nul[71]-0.1530931089239486*fr[25]*nul[67]+0.1530931089239486*fl[25]*nul[67]+0.1325825214724776*fr[15]*nul[67]+0.1325825214724776*fl[15]*nul[67]-0.0883883476483184*fr[29]*nul[65]+0.0883883476483184*fl[29]*nul[65]+0.07654655446197427*fr[23]*nul[65]+0.07654655446197427*fl[23]*nul[65]-0.0883883476483184*fr[25]*nul[64]+0.0883883476483184*fl[25]*nul[64]+0.07654655446197427*fr[15]*nul[64]+0.07654655446197427*fl[15]*nul[64]; 
  incr2[26] = (-0.1530931089239486*fr[19]*nul[71])+0.1530931089239486*fl[19]*nul[71]+0.1325825214724776*fr[10]*nul[71]+0.1325825214724776*fl[10]*nul[71]-0.1530931089239486*fr[26]*nul[67]+0.1530931089239486*fl[26]*nul[67]+0.1325825214724776*fr[17]*nul[67]+0.1325825214724776*fl[17]*nul[67]-0.0883883476483184*fr[19]*nul[65]+0.0883883476483184*fl[19]*nul[65]+0.07654655446197427*fr[10]*nul[65]+0.07654655446197427*fl[10]*nul[65]-0.0883883476483184*fr[26]*nul[64]+0.0883883476483184*fl[26]*nul[64]+0.07654655446197427*fr[17]*nul[64]+0.07654655446197427*fl[17]*nul[64]; 
  incr2[27] = (-0.1530931089239486*fr[22]*nul[71])+0.1530931089239486*fl[22]*nul[71]+0.1325825214724776*fr[13]*nul[71]+0.1325825214724776*fl[13]*nul[71]-0.1530931089239486*fr[27]*nul[67]+0.1530931089239486*fl[27]*nul[67]+0.1325825214724776*fr[20]*nul[67]+0.1325825214724776*fl[20]*nul[67]-0.0883883476483184*fr[22]*nul[65]+0.0883883476483184*fl[22]*nul[65]+0.07654655446197427*fr[13]*nul[65]+0.07654655446197427*fl[13]*nul[65]-0.0883883476483184*fr[27]*nul[64]+0.0883883476483184*fl[27]*nul[64]+0.07654655446197427*fr[20]*nul[64]+0.07654655446197427*fl[20]*nul[64]; 
  incr2[29] = (-0.1530931089239486*fr[25]*nul[71])+0.1530931089239486*fl[25]*nul[71]+0.1325825214724776*fr[15]*nul[71]+0.1325825214724776*fl[15]*nul[71]-0.1530931089239486*fr[29]*nul[67]+0.1530931089239486*fl[29]*nul[67]+0.1325825214724776*fr[23]*nul[67]+0.1325825214724776*fl[23]*nul[67]-0.0883883476483184*fr[25]*nul[65]+0.0883883476483184*fl[25]*nul[65]+0.07654655446197427*fr[15]*nul[65]+0.07654655446197427*fl[15]*nul[65]-0.0883883476483184*fr[29]*nul[64]+0.0883883476483184*fl[29]*nul[64]+0.07654655446197427*fr[23]*nul[64]+0.07654655446197427*fl[23]*nul[64]; 
  incr2[30] = (-0.1530931089239486*fr[31]*nul[71])+0.1530931089239486*fl[31]*nul[71]+0.1325825214724776*fr[28]*nul[71]+0.1325825214724776*fl[28]*nul[71]-0.1530931089239486*fr[30]*nul[67]+0.1530931089239486*fl[30]*nul[67]+0.1325825214724776*fr[24]*nul[67]+0.1325825214724776*fl[24]*nul[67]-0.0883883476483184*fr[31]*nul[65]+0.0883883476483184*fl[31]*nul[65]+0.07654655446197427*fr[28]*nul[65]+0.07654655446197427*fl[28]*nul[65]-0.0883883476483184*fr[30]*nul[64]+0.0883883476483184*fl[30]*nul[64]+0.07654655446197427*fr[24]*nul[64]+0.07654655446197427*fl[24]*nul[64]; 
  incr2[31] = (-0.1530931089239486*fr[30]*nul[71])+0.1530931089239486*fl[30]*nul[71]+0.1325825214724776*fr[24]*nul[71]+0.1325825214724776*fl[24]*nul[71]-0.1530931089239486*fr[31]*nul[67]+0.1530931089239486*fl[31]*nul[67]+0.1325825214724776*fr[28]*nul[67]+0.1325825214724776*fl[28]*nul[67]-0.0883883476483184*fr[30]*nul[65]+0.0883883476483184*fl[30]*nul[65]+0.07654655446197427*fr[24]*nul[65]+0.07654655446197427*fl[24]*nul[65]-0.0883883476483184*fr[31]*nul[64]+0.0883883476483184*fl[31]*nul[64]+0.07654655446197427*fr[28]*nul[64]+0.07654655446197427*fl[28]*nul[64]; 

  outr[0] += incr1[0]*rdxFr; 
  outr[1] += incr1[1]*rdxFr; 
  outr[2] += incr1[2]*rdxFr; 
  outr[3] += incr2[3]*rdxFr+incr1[3]*rdxFr; 
  outr[4] += incr1[4]*rdxFr; 
  outr[5] += incr1[5]*rdxFr; 
  outr[6] += incr1[6]*rdxFr; 
  outr[7] += incr2[7]*rdxFr+incr1[7]*rdxFr; 
  outr[8] += incr2[8]*rdxFr+incr1[8]*rdxFr; 
  outr[9] += incr1[9]*rdxFr; 
  outr[10] += incr1[10]*rdxFr; 
  outr[11] += incr2[11]*rdxFr+incr1[11]*rdxFr; 
  outr[12] += incr1[12]*rdxFr; 
  outr[13] += incr1[13]*rdxFr; 
  outr[14] += incr2[14]*rdxFr+incr1[14]*rdxFr; 
  outr[15] += incr1[15]*rdxFr; 
  outr[16] += incr2[16]*rdxFr+incr1[16]*rdxFr; 
  outr[17] += incr1[17]*rdxFr; 
  outr[18] += incr2[18]*rdxFr+incr1[18]*rdxFr; 
  outr[19] += incr2[19]*rdxFr+incr1[19]*rdxFr; 
  outr[20] += incr1[20]*rdxFr; 
  outr[21] += incr2[21]*rdxFr+incr1[21]*rdxFr; 
  outr[22] += incr2[22]*rdxFr+incr1[22]*rdxFr; 
  outr[23] += incr1[23]*rdxFr; 
  outr[24] += incr1[24]*rdxFr; 
  outr[25] += incr2[25]*rdxFr+incr1[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr+incr1[26]*rdxFr; 
  outr[27] += incr2[27]*rdxFr+incr1[27]*rdxFr; 
  outr[28] += incr1[28]*rdxFr; 
  outr[29] += incr2[29]*rdxFr+incr1[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr+incr1[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr+incr1[31]*rdxFr; 

  outl[0] += -1.0*incr1[0]*rdxFl; 
  outl[1] += -1.0*incr1[1]*rdxFl; 
  outl[2] += -1.0*incr1[2]*rdxFl; 
  outl[3] += incr1[3]*rdxFl-1.0*incr2[3]*rdxFl; 
  outl[4] += -1.0*incr1[4]*rdxFl; 
  outl[5] += -1.0*incr1[5]*rdxFl; 
  outl[6] += -1.0*incr1[6]*rdxFl; 
  outl[7] += incr1[7]*rdxFl-1.0*incr2[7]*rdxFl; 
  outl[8] += incr1[8]*rdxFl-1.0*incr2[8]*rdxFl; 
  outl[9] += -1.0*incr1[9]*rdxFl; 
  outl[10] += -1.0*incr1[10]*rdxFl; 
  outl[11] += incr1[11]*rdxFl-1.0*incr2[11]*rdxFl; 
  outl[12] += -1.0*incr1[12]*rdxFl; 
  outl[13] += -1.0*incr1[13]*rdxFl; 
  outl[14] += incr1[14]*rdxFl-1.0*incr2[14]*rdxFl; 
  outl[15] += -1.0*incr1[15]*rdxFl; 
  outl[16] += incr1[16]*rdxFl-1.0*incr2[16]*rdxFl; 
  outl[17] += -1.0*incr1[17]*rdxFl; 
  outl[18] += incr1[18]*rdxFl-1.0*incr2[18]*rdxFl; 
  outl[19] += incr1[19]*rdxFl-1.0*incr2[19]*rdxFl; 
  outl[20] += -1.0*incr1[20]*rdxFl; 
  outl[21] += incr1[21]*rdxFl-1.0*incr2[21]*rdxFl; 
  outl[22] += incr1[22]*rdxFl-1.0*incr2[22]*rdxFl; 
  outl[23] += -1.0*incr1[23]*rdxFl; 
  outl[24] += -1.0*incr1[24]*rdxFl; 
  outl[25] += incr1[25]*rdxFl-1.0*incr2[25]*rdxFl; 
  outl[26] += incr1[26]*rdxFl-1.0*incr2[26]*rdxFl; 
  outl[27] += incr1[27]*rdxFl-1.0*incr2[27]*rdxFl; 
  outl[28] += -1.0*incr1[28]*rdxFl; 
  outl[29] += incr1[29]*rdxFl-1.0*incr2[29]*rdxFl; 
  outl[30] += incr1[30]*rdxFl-1.0*incr2[30]*rdxFl; 
  outl[31] += incr1[31]*rdxFl-1.0*incr2[31]*rdxFl; 

} 
void ConstDiffusionVarCoeffSurf5xSerP1_X4(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[160]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[3]*dxl[3]); 
  double rdxFr = 4.0/(dxr[3]*dxr[3]); 

  double incr1[32]; 
  double incr2[32]; 

  incr1[0] = 0.09568319307746781*fr[18]*nul[103]+0.09568319307746781*fl[18]*nul[103]-0.09943689110435816*fr[7]*nul[103]+0.09943689110435816*fl[7]*nul[103]+0.09568319307746781*fr[11]*nul[99]+0.09568319307746781*fl[11]*nul[99]-0.09943689110435816*fr[3]*nul[99]+0.09943689110435816*fl[3]*nul[99]+0.09568319307746781*fr[9]*nul[97]+0.09568319307746781*fl[9]*nul[97]-0.09943689110435816*fr[1]*nul[97]+0.09943689110435816*fl[1]*nul[97]+0.09568319307746781*fr[4]*nul[96]+0.09568319307746781*fl[4]*nul[96]-0.09943689110435816*fr[0]*nul[96]+0.09943689110435816*fl[0]*nul[96]; 
  incr1[1] = 0.09568319307746781*fr[11]*nul[103]+0.09568319307746781*fl[11]*nul[103]-0.09943689110435816*fr[3]*nul[103]+0.09943689110435816*fl[3]*nul[103]+0.09568319307746781*fr[18]*nul[99]+0.09568319307746781*fl[18]*nul[99]-0.09943689110435816*fr[7]*nul[99]+0.09943689110435816*fl[7]*nul[99]+0.09568319307746781*fr[4]*nul[97]+0.09568319307746781*fl[4]*nul[97]-0.09943689110435816*fr[0]*nul[97]+0.09943689110435816*fl[0]*nul[97]+0.09568319307746781*fr[9]*nul[96]+0.09568319307746781*fl[9]*nul[96]-0.09943689110435816*fr[1]*nul[96]+0.09943689110435816*fl[1]*nul[96]; 
  incr1[2] = 0.09568319307746781*fr[26]*nul[103]+0.09568319307746781*fl[26]*nul[103]-0.09943689110435816*fr[16]*nul[103]+0.09943689110435816*fl[16]*nul[103]+0.09568319307746781*fr[19]*nul[99]+0.09568319307746781*fl[19]*nul[99]-0.09943689110435816*fr[8]*nul[99]+0.09943689110435816*fl[8]*nul[99]+0.09568319307746781*fr[17]*nul[97]+0.09568319307746781*fl[17]*nul[97]-0.09943689110435816*fr[6]*nul[97]+0.09943689110435816*fl[6]*nul[97]+0.09568319307746781*fr[10]*nul[96]+0.09568319307746781*fl[10]*nul[96]-0.09943689110435816*fr[2]*nul[96]+0.09943689110435816*fl[2]*nul[96]; 
  incr1[3] = 0.09568319307746781*fr[9]*nul[103]+0.09568319307746781*fl[9]*nul[103]-0.09943689110435816*fr[1]*nul[103]+0.09943689110435816*fl[1]*nul[103]+0.09568319307746781*fr[4]*nul[99]+0.09568319307746781*fl[4]*nul[99]-0.09943689110435816*fr[0]*nul[99]+0.09943689110435816*fl[0]*nul[99]+0.09568319307746781*fr[18]*nul[97]+0.09568319307746781*fl[18]*nul[97]-0.09943689110435816*fr[7]*nul[97]+0.09943689110435816*fl[7]*nul[97]+0.09568319307746781*fr[11]*nul[96]+0.09568319307746781*fl[11]*nul[96]-0.09943689110435816*fr[3]*nul[96]+0.09943689110435816*fl[3]*nul[96]; 
  incr1[4] = (-0.1657281518405969*fr[18]*nul[103])-0.1657281518405969*fl[18]*nul[103]+0.172229747539442*fr[7]*nul[103]-0.172229747539442*fl[7]*nul[103]-0.1657281518405969*fr[11]*nul[99]-0.1657281518405969*fl[11]*nul[99]+0.172229747539442*fr[3]*nul[99]-0.172229747539442*fl[3]*nul[99]-0.1657281518405969*fr[9]*nul[97]-0.1657281518405969*fl[9]*nul[97]+0.172229747539442*fr[1]*nul[97]-0.172229747539442*fl[1]*nul[97]-0.1657281518405969*fr[4]*nul[96]-0.1657281518405969*fl[4]*nul[96]+0.172229747539442*fr[0]*nul[96]-0.172229747539442*fl[0]*nul[96]; 
  incr1[5] = 0.09568319307746781*fr[29]*nul[103]+0.09568319307746781*fl[29]*nul[103]-0.09943689110435816*fr[21]*nul[103]+0.09943689110435816*fl[21]*nul[103]+0.09568319307746781*fr[25]*nul[99]+0.09568319307746781*fl[25]*nul[99]-0.09943689110435816*fr[14]*nul[99]+0.09943689110435816*fl[14]*nul[99]+0.09568319307746781*fr[23]*nul[97]+0.09568319307746781*fl[23]*nul[97]-0.09943689110435816*fr[12]*nul[97]+0.09943689110435816*fl[12]*nul[97]+0.09568319307746781*fr[15]*nul[96]+0.09568319307746781*fl[15]*nul[96]-0.09943689110435816*fr[5]*nul[96]+0.09943689110435816*fl[5]*nul[96]; 
  incr1[6] = 0.09568319307746781*fr[19]*nul[103]+0.09568319307746781*fl[19]*nul[103]-0.09943689110435816*fr[8]*nul[103]+0.09943689110435816*fl[8]*nul[103]+0.09568319307746781*fr[26]*nul[99]+0.09568319307746781*fl[26]*nul[99]-0.09943689110435816*fr[16]*nul[99]+0.09943689110435816*fl[16]*nul[99]+0.09568319307746781*fr[10]*nul[97]+0.09568319307746781*fl[10]*nul[97]-0.09943689110435816*fr[2]*nul[97]+0.09943689110435816*fl[2]*nul[97]+0.09568319307746781*fr[17]*nul[96]+0.09568319307746781*fl[17]*nul[96]-0.09943689110435816*fr[6]*nul[96]+0.09943689110435816*fl[6]*nul[96]; 
  incr1[7] = 0.09568319307746781*fr[4]*nul[103]+0.09568319307746781*fl[4]*nul[103]-0.09943689110435816*fr[0]*nul[103]+0.09943689110435816*fl[0]*nul[103]+0.09568319307746781*fr[9]*nul[99]+0.09568319307746781*fl[9]*nul[99]-0.09943689110435816*fr[1]*nul[99]+0.09943689110435816*fl[1]*nul[99]+0.09568319307746781*fr[11]*nul[97]+0.09568319307746781*fl[11]*nul[97]-0.09943689110435816*fr[3]*nul[97]+0.09943689110435816*fl[3]*nul[97]+0.09568319307746781*fr[18]*nul[96]+0.09568319307746781*fl[18]*nul[96]-0.09943689110435816*fr[7]*nul[96]+0.09943689110435816*fl[7]*nul[96]; 
  incr1[8] = 0.09568319307746781*fr[17]*nul[103]+0.09568319307746781*fl[17]*nul[103]-0.09943689110435816*fr[6]*nul[103]+0.09943689110435816*fl[6]*nul[103]+0.09568319307746781*fr[10]*nul[99]+0.09568319307746781*fl[10]*nul[99]-0.09943689110435816*fr[2]*nul[99]+0.09943689110435816*fl[2]*nul[99]+0.09568319307746781*fr[26]*nul[97]+0.09568319307746781*fl[26]*nul[97]-0.09943689110435816*fr[16]*nul[97]+0.09943689110435816*fl[16]*nul[97]+0.09568319307746781*fr[19]*nul[96]+0.09568319307746781*fl[19]*nul[96]-0.09943689110435816*fr[8]*nul[96]+0.09943689110435816*fl[8]*nul[96]; 
  incr1[9] = (-0.1657281518405969*fr[11]*nul[103])-0.1657281518405969*fl[11]*nul[103]+0.172229747539442*fr[3]*nul[103]-0.172229747539442*fl[3]*nul[103]-0.1657281518405969*fr[18]*nul[99]-0.1657281518405969*fl[18]*nul[99]+0.172229747539442*fr[7]*nul[99]-0.172229747539442*fl[7]*nul[99]-0.1657281518405969*fr[4]*nul[97]-0.1657281518405969*fl[4]*nul[97]+0.172229747539442*fr[0]*nul[97]-0.172229747539442*fl[0]*nul[97]-0.1657281518405969*fr[9]*nul[96]-0.1657281518405969*fl[9]*nul[96]+0.172229747539442*fr[1]*nul[96]-0.172229747539442*fl[1]*nul[96]; 
  incr1[10] = (-0.1657281518405969*fr[26]*nul[103])-0.1657281518405969*fl[26]*nul[103]+0.172229747539442*fr[16]*nul[103]-0.172229747539442*fl[16]*nul[103]-0.1657281518405969*fr[19]*nul[99]-0.1657281518405969*fl[19]*nul[99]+0.172229747539442*fr[8]*nul[99]-0.172229747539442*fl[8]*nul[99]-0.1657281518405969*fr[17]*nul[97]-0.1657281518405969*fl[17]*nul[97]+0.172229747539442*fr[6]*nul[97]-0.172229747539442*fl[6]*nul[97]-0.1657281518405969*fr[10]*nul[96]-0.1657281518405969*fl[10]*nul[96]+0.172229747539442*fr[2]*nul[96]-0.172229747539442*fl[2]*nul[96]; 
  incr1[11] = (-0.1657281518405969*fr[9]*nul[103])-0.1657281518405969*fl[9]*nul[103]+0.172229747539442*fr[1]*nul[103]-0.172229747539442*fl[1]*nul[103]-0.1657281518405969*fr[4]*nul[99]-0.1657281518405969*fl[4]*nul[99]+0.172229747539442*fr[0]*nul[99]-0.172229747539442*fl[0]*nul[99]-0.1657281518405969*fr[18]*nul[97]-0.1657281518405969*fl[18]*nul[97]+0.172229747539442*fr[7]*nul[97]-0.172229747539442*fl[7]*nul[97]-0.1657281518405969*fr[11]*nul[96]-0.1657281518405969*fl[11]*nul[96]+0.172229747539442*fr[3]*nul[96]-0.172229747539442*fl[3]*nul[96]; 
  incr1[12] = 0.09568319307746781*fr[25]*nul[103]+0.09568319307746781*fl[25]*nul[103]-0.09943689110435816*fr[14]*nul[103]+0.09943689110435816*fl[14]*nul[103]+0.09568319307746781*fr[29]*nul[99]+0.09568319307746781*fl[29]*nul[99]-0.09943689110435816*fr[21]*nul[99]+0.09943689110435816*fl[21]*nul[99]+0.09568319307746781*fr[15]*nul[97]+0.09568319307746781*fl[15]*nul[97]-0.09943689110435816*fr[5]*nul[97]+0.09943689110435816*fl[5]*nul[97]+0.09568319307746781*fr[23]*nul[96]+0.09568319307746781*fl[23]*nul[96]-0.09943689110435816*fr[12]*nul[96]+0.09943689110435816*fl[12]*nul[96]; 
  incr1[13] = 0.09568319307746781*fr[31]*nul[103]+0.09568319307746781*fl[31]*nul[103]-0.09943689110435816*fr[27]*nul[103]+0.09943689110435816*fl[27]*nul[103]+0.09568319307746781*fr[30]*nul[99]+0.09568319307746781*fl[30]*nul[99]-0.09943689110435816*fr[22]*nul[99]+0.09943689110435816*fl[22]*nul[99]+0.09568319307746781*fr[28]*nul[97]+0.09568319307746781*fl[28]*nul[97]-0.09943689110435816*fr[20]*nul[97]+0.09943689110435816*fl[20]*nul[97]+0.09568319307746781*fr[24]*nul[96]+0.09568319307746781*fl[24]*nul[96]-0.09943689110435816*fr[13]*nul[96]+0.09943689110435816*fl[13]*nul[96]; 
  incr1[14] = 0.09568319307746781*fr[23]*nul[103]+0.09568319307746781*fl[23]*nul[103]-0.09943689110435816*fr[12]*nul[103]+0.09943689110435816*fl[12]*nul[103]+0.09568319307746781*fr[15]*nul[99]+0.09568319307746781*fl[15]*nul[99]-0.09943689110435816*fr[5]*nul[99]+0.09943689110435816*fl[5]*nul[99]+0.09568319307746781*fr[29]*nul[97]+0.09568319307746781*fl[29]*nul[97]-0.09943689110435816*fr[21]*nul[97]+0.09943689110435816*fl[21]*nul[97]+0.09568319307746781*fr[25]*nul[96]+0.09568319307746781*fl[25]*nul[96]-0.09943689110435816*fr[14]*nul[96]+0.09943689110435816*fl[14]*nul[96]; 
  incr1[15] = (-0.1657281518405969*fr[29]*nul[103])-0.1657281518405969*fl[29]*nul[103]+0.172229747539442*fr[21]*nul[103]-0.172229747539442*fl[21]*nul[103]-0.1657281518405969*fr[25]*nul[99]-0.1657281518405969*fl[25]*nul[99]+0.172229747539442*fr[14]*nul[99]-0.172229747539442*fl[14]*nul[99]-0.1657281518405969*fr[23]*nul[97]-0.1657281518405969*fl[23]*nul[97]+0.172229747539442*fr[12]*nul[97]-0.172229747539442*fl[12]*nul[97]-0.1657281518405969*fr[15]*nul[96]-0.1657281518405969*fl[15]*nul[96]+0.172229747539442*fr[5]*nul[96]-0.172229747539442*fl[5]*nul[96]; 
  incr1[16] = 0.09568319307746781*fr[10]*nul[103]+0.09568319307746781*fl[10]*nul[103]-0.09943689110435816*fr[2]*nul[103]+0.09943689110435816*fl[2]*nul[103]+0.09568319307746781*fr[17]*nul[99]+0.09568319307746781*fl[17]*nul[99]-0.09943689110435816*fr[6]*nul[99]+0.09943689110435816*fl[6]*nul[99]+0.09568319307746781*fr[19]*nul[97]+0.09568319307746781*fl[19]*nul[97]-0.09943689110435816*fr[8]*nul[97]+0.09943689110435816*fl[8]*nul[97]+0.09568319307746781*fr[26]*nul[96]+0.09568319307746781*fl[26]*nul[96]-0.09943689110435816*fr[16]*nul[96]+0.09943689110435816*fl[16]*nul[96]; 
  incr1[17] = (-0.1657281518405969*fr[19]*nul[103])-0.1657281518405969*fl[19]*nul[103]+0.172229747539442*fr[8]*nul[103]-0.172229747539442*fl[8]*nul[103]-0.1657281518405969*fr[26]*nul[99]-0.1657281518405969*fl[26]*nul[99]+0.172229747539442*fr[16]*nul[99]-0.172229747539442*fl[16]*nul[99]-0.1657281518405969*fr[10]*nul[97]-0.1657281518405969*fl[10]*nul[97]+0.172229747539442*fr[2]*nul[97]-0.172229747539442*fl[2]*nul[97]-0.1657281518405969*fr[17]*nul[96]-0.1657281518405969*fl[17]*nul[96]+0.172229747539442*fr[6]*nul[96]-0.172229747539442*fl[6]*nul[96]; 
  incr1[18] = (-0.1657281518405969*fr[4]*nul[103])-0.1657281518405969*fl[4]*nul[103]+0.172229747539442*fr[0]*nul[103]-0.172229747539442*fl[0]*nul[103]-0.1657281518405969*fr[9]*nul[99]-0.1657281518405969*fl[9]*nul[99]+0.172229747539442*fr[1]*nul[99]-0.172229747539442*fl[1]*nul[99]-0.1657281518405969*fr[11]*nul[97]-0.1657281518405969*fl[11]*nul[97]+0.172229747539442*fr[3]*nul[97]-0.172229747539442*fl[3]*nul[97]-0.1657281518405969*fr[18]*nul[96]-0.1657281518405969*fl[18]*nul[96]+0.172229747539442*fr[7]*nul[96]-0.172229747539442*fl[7]*nul[96]; 
  incr1[19] = (-0.1657281518405969*fr[17]*nul[103])-0.1657281518405969*fl[17]*nul[103]+0.172229747539442*fr[6]*nul[103]-0.172229747539442*fl[6]*nul[103]-0.1657281518405969*fr[10]*nul[99]-0.1657281518405969*fl[10]*nul[99]+0.172229747539442*fr[2]*nul[99]-0.172229747539442*fl[2]*nul[99]-0.1657281518405969*fr[26]*nul[97]-0.1657281518405969*fl[26]*nul[97]+0.172229747539442*fr[16]*nul[97]-0.172229747539442*fl[16]*nul[97]-0.1657281518405969*fr[19]*nul[96]-0.1657281518405969*fl[19]*nul[96]+0.172229747539442*fr[8]*nul[96]-0.172229747539442*fl[8]*nul[96]; 
  incr1[20] = 0.09568319307746781*fr[30]*nul[103]+0.09568319307746781*fl[30]*nul[103]-0.09943689110435816*fr[22]*nul[103]+0.09943689110435816*fl[22]*nul[103]+0.09568319307746781*fr[31]*nul[99]+0.09568319307746781*fl[31]*nul[99]-0.09943689110435816*fr[27]*nul[99]+0.09943689110435816*fl[27]*nul[99]+0.09568319307746781*fr[24]*nul[97]+0.09568319307746781*fl[24]*nul[97]-0.09943689110435816*fr[13]*nul[97]+0.09943689110435816*fl[13]*nul[97]+0.09568319307746781*fr[28]*nul[96]+0.09568319307746781*fl[28]*nul[96]-0.09943689110435816*fr[20]*nul[96]+0.09943689110435816*fl[20]*nul[96]; 
  incr1[21] = 0.09568319307746781*fr[15]*nul[103]+0.09568319307746781*fl[15]*nul[103]-0.09943689110435816*fr[5]*nul[103]+0.09943689110435816*fl[5]*nul[103]+0.09568319307746781*fr[23]*nul[99]+0.09568319307746781*fl[23]*nul[99]-0.09943689110435816*fr[12]*nul[99]+0.09943689110435816*fl[12]*nul[99]+0.09568319307746781*fr[25]*nul[97]+0.09568319307746781*fl[25]*nul[97]-0.09943689110435816*fr[14]*nul[97]+0.09943689110435816*fl[14]*nul[97]+0.09568319307746781*fr[29]*nul[96]+0.09568319307746781*fl[29]*nul[96]-0.09943689110435816*fr[21]*nul[96]+0.09943689110435816*fl[21]*nul[96]; 
  incr1[22] = 0.09568319307746781*fr[28]*nul[103]+0.09568319307746781*fl[28]*nul[103]-0.09943689110435816*fr[20]*nul[103]+0.09943689110435816*fl[20]*nul[103]+0.09568319307746781*fr[24]*nul[99]+0.09568319307746781*fl[24]*nul[99]-0.09943689110435816*fr[13]*nul[99]+0.09943689110435816*fl[13]*nul[99]+0.09568319307746781*fr[31]*nul[97]+0.09568319307746781*fl[31]*nul[97]-0.09943689110435816*fr[27]*nul[97]+0.09943689110435816*fl[27]*nul[97]+0.09568319307746781*fr[30]*nul[96]+0.09568319307746781*fl[30]*nul[96]-0.09943689110435816*fr[22]*nul[96]+0.09943689110435816*fl[22]*nul[96]; 
  incr1[23] = (-0.1657281518405969*fr[25]*nul[103])-0.1657281518405969*fl[25]*nul[103]+0.172229747539442*fr[14]*nul[103]-0.172229747539442*fl[14]*nul[103]-0.1657281518405969*fr[29]*nul[99]-0.1657281518405969*fl[29]*nul[99]+0.172229747539442*fr[21]*nul[99]-0.172229747539442*fl[21]*nul[99]-0.1657281518405969*fr[15]*nul[97]-0.1657281518405969*fl[15]*nul[97]+0.172229747539442*fr[5]*nul[97]-0.172229747539442*fl[5]*nul[97]-0.1657281518405969*fr[23]*nul[96]-0.1657281518405969*fl[23]*nul[96]+0.172229747539442*fr[12]*nul[96]-0.172229747539442*fl[12]*nul[96]; 
  incr1[24] = (-0.1657281518405969*fr[31]*nul[103])-0.1657281518405969*fl[31]*nul[103]+0.172229747539442*fr[27]*nul[103]-0.172229747539442*fl[27]*nul[103]-0.1657281518405969*fr[30]*nul[99]-0.1657281518405969*fl[30]*nul[99]+0.172229747539442*fr[22]*nul[99]-0.172229747539442*fl[22]*nul[99]-0.1657281518405969*fr[28]*nul[97]-0.1657281518405969*fl[28]*nul[97]+0.172229747539442*fr[20]*nul[97]-0.172229747539442*fl[20]*nul[97]-0.1657281518405969*fr[24]*nul[96]-0.1657281518405969*fl[24]*nul[96]+0.172229747539442*fr[13]*nul[96]-0.172229747539442*fl[13]*nul[96]; 
  incr1[25] = (-0.1657281518405969*fr[23]*nul[103])-0.1657281518405969*fl[23]*nul[103]+0.172229747539442*fr[12]*nul[103]-0.172229747539442*fl[12]*nul[103]-0.1657281518405969*fr[15]*nul[99]-0.1657281518405969*fl[15]*nul[99]+0.172229747539442*fr[5]*nul[99]-0.172229747539442*fl[5]*nul[99]-0.1657281518405969*fr[29]*nul[97]-0.1657281518405969*fl[29]*nul[97]+0.172229747539442*fr[21]*nul[97]-0.172229747539442*fl[21]*nul[97]-0.1657281518405969*fr[25]*nul[96]-0.1657281518405969*fl[25]*nul[96]+0.172229747539442*fr[14]*nul[96]-0.172229747539442*fl[14]*nul[96]; 
  incr1[26] = (-0.1657281518405969*fr[10]*nul[103])-0.1657281518405969*fl[10]*nul[103]+0.172229747539442*fr[2]*nul[103]-0.172229747539442*fl[2]*nul[103]-0.1657281518405969*fr[17]*nul[99]-0.1657281518405969*fl[17]*nul[99]+0.172229747539442*fr[6]*nul[99]-0.172229747539442*fl[6]*nul[99]-0.1657281518405969*fr[19]*nul[97]-0.1657281518405969*fl[19]*nul[97]+0.172229747539442*fr[8]*nul[97]-0.172229747539442*fl[8]*nul[97]-0.1657281518405969*fr[26]*nul[96]-0.1657281518405969*fl[26]*nul[96]+0.172229747539442*fr[16]*nul[96]-0.172229747539442*fl[16]*nul[96]; 
  incr1[27] = 0.09568319307746781*fr[24]*nul[103]+0.09568319307746781*fl[24]*nul[103]-0.09943689110435816*fr[13]*nul[103]+0.09943689110435816*fl[13]*nul[103]+0.09568319307746781*fr[28]*nul[99]+0.09568319307746781*fl[28]*nul[99]-0.09943689110435816*fr[20]*nul[99]+0.09943689110435816*fl[20]*nul[99]+0.09568319307746781*fr[30]*nul[97]+0.09568319307746781*fl[30]*nul[97]-0.09943689110435816*fr[22]*nul[97]+0.09943689110435816*fl[22]*nul[97]+0.09568319307746781*fr[31]*nul[96]+0.09568319307746781*fl[31]*nul[96]-0.09943689110435816*fr[27]*nul[96]+0.09943689110435816*fl[27]*nul[96]; 
  incr1[28] = (-0.1657281518405969*fr[30]*nul[103])-0.1657281518405969*fl[30]*nul[103]+0.172229747539442*fr[22]*nul[103]-0.172229747539442*fl[22]*nul[103]-0.1657281518405969*fr[31]*nul[99]-0.1657281518405969*fl[31]*nul[99]+0.172229747539442*fr[27]*nul[99]-0.172229747539442*fl[27]*nul[99]-0.1657281518405969*fr[24]*nul[97]-0.1657281518405969*fl[24]*nul[97]+0.172229747539442*fr[13]*nul[97]-0.172229747539442*fl[13]*nul[97]-0.1657281518405969*fr[28]*nul[96]-0.1657281518405969*fl[28]*nul[96]+0.172229747539442*fr[20]*nul[96]-0.172229747539442*fl[20]*nul[96]; 
  incr1[29] = (-0.1657281518405969*fr[15]*nul[103])-0.1657281518405969*fl[15]*nul[103]+0.172229747539442*fr[5]*nul[103]-0.172229747539442*fl[5]*nul[103]-0.1657281518405969*fr[23]*nul[99]-0.1657281518405969*fl[23]*nul[99]+0.172229747539442*fr[12]*nul[99]-0.172229747539442*fl[12]*nul[99]-0.1657281518405969*fr[25]*nul[97]-0.1657281518405969*fl[25]*nul[97]+0.172229747539442*fr[14]*nul[97]-0.172229747539442*fl[14]*nul[97]-0.1657281518405969*fr[29]*nul[96]-0.1657281518405969*fl[29]*nul[96]+0.172229747539442*fr[21]*nul[96]-0.172229747539442*fl[21]*nul[96]; 
  incr1[30] = (-0.1657281518405969*fr[28]*nul[103])-0.1657281518405969*fl[28]*nul[103]+0.172229747539442*fr[20]*nul[103]-0.172229747539442*fl[20]*nul[103]-0.1657281518405969*fr[24]*nul[99]-0.1657281518405969*fl[24]*nul[99]+0.172229747539442*fr[13]*nul[99]-0.172229747539442*fl[13]*nul[99]-0.1657281518405969*fr[31]*nul[97]-0.1657281518405969*fl[31]*nul[97]+0.172229747539442*fr[27]*nul[97]-0.172229747539442*fl[27]*nul[97]-0.1657281518405969*fr[30]*nul[96]-0.1657281518405969*fl[30]*nul[96]+0.172229747539442*fr[22]*nul[96]-0.172229747539442*fl[22]*nul[96]; 
  incr1[31] = (-0.1657281518405969*fr[24]*nul[103])-0.1657281518405969*fl[24]*nul[103]+0.172229747539442*fr[13]*nul[103]-0.172229747539442*fl[13]*nul[103]-0.1657281518405969*fr[28]*nul[99]-0.1657281518405969*fl[28]*nul[99]+0.172229747539442*fr[20]*nul[99]-0.172229747539442*fl[20]*nul[99]-0.1657281518405969*fr[30]*nul[97]-0.1657281518405969*fl[30]*nul[97]+0.172229747539442*fr[22]*nul[97]-0.172229747539442*fl[22]*nul[97]-0.1657281518405969*fr[31]*nul[96]-0.1657281518405969*fl[31]*nul[96]+0.172229747539442*fr[27]*nul[96]-0.172229747539442*fl[27]*nul[96]; 

  incr2[4] = (-0.0883883476483184*fr[18]*nul[103])+0.0883883476483184*fl[18]*nul[103]+0.07654655446197427*fr[7]*nul[103]+0.07654655446197427*fl[7]*nul[103]-0.0883883476483184*fr[11]*nul[99]+0.0883883476483184*fl[11]*nul[99]+0.07654655446197427*fr[3]*nul[99]+0.07654655446197427*fl[3]*nul[99]-0.0883883476483184*fr[9]*nul[97]+0.0883883476483184*fl[9]*nul[97]+0.07654655446197427*fr[1]*nul[97]+0.07654655446197427*fl[1]*nul[97]-0.0883883476483184*fr[4]*nul[96]+0.0883883476483184*fl[4]*nul[96]+0.07654655446197427*fr[0]*nul[96]+0.07654655446197427*fl[0]*nul[96]; 
  incr2[9] = (-0.0883883476483184*fr[11]*nul[103])+0.0883883476483184*fl[11]*nul[103]+0.07654655446197427*fr[3]*nul[103]+0.07654655446197427*fl[3]*nul[103]-0.0883883476483184*fr[18]*nul[99]+0.0883883476483184*fl[18]*nul[99]+0.07654655446197427*fr[7]*nul[99]+0.07654655446197427*fl[7]*nul[99]-0.0883883476483184*fr[4]*nul[97]+0.0883883476483184*fl[4]*nul[97]+0.07654655446197427*fr[0]*nul[97]+0.07654655446197427*fl[0]*nul[97]-0.0883883476483184*fr[9]*nul[96]+0.0883883476483184*fl[9]*nul[96]+0.07654655446197427*fr[1]*nul[96]+0.07654655446197427*fl[1]*nul[96]; 
  incr2[10] = (-0.0883883476483184*fr[26]*nul[103])+0.0883883476483184*fl[26]*nul[103]+0.07654655446197427*fr[16]*nul[103]+0.07654655446197427*fl[16]*nul[103]-0.0883883476483184*fr[19]*nul[99]+0.0883883476483184*fl[19]*nul[99]+0.07654655446197427*fr[8]*nul[99]+0.07654655446197427*fl[8]*nul[99]-0.0883883476483184*fr[17]*nul[97]+0.0883883476483184*fl[17]*nul[97]+0.07654655446197427*fr[6]*nul[97]+0.07654655446197427*fl[6]*nul[97]-0.0883883476483184*fr[10]*nul[96]+0.0883883476483184*fl[10]*nul[96]+0.07654655446197427*fr[2]*nul[96]+0.07654655446197427*fl[2]*nul[96]; 
  incr2[11] = (-0.0883883476483184*fr[9]*nul[103])+0.0883883476483184*fl[9]*nul[103]+0.07654655446197427*fr[1]*nul[103]+0.07654655446197427*fl[1]*nul[103]-0.0883883476483184*fr[4]*nul[99]+0.0883883476483184*fl[4]*nul[99]+0.07654655446197427*fr[0]*nul[99]+0.07654655446197427*fl[0]*nul[99]-0.0883883476483184*fr[18]*nul[97]+0.0883883476483184*fl[18]*nul[97]+0.07654655446197427*fr[7]*nul[97]+0.07654655446197427*fl[7]*nul[97]-0.0883883476483184*fr[11]*nul[96]+0.0883883476483184*fl[11]*nul[96]+0.07654655446197427*fr[3]*nul[96]+0.07654655446197427*fl[3]*nul[96]; 
  incr2[15] = (-0.0883883476483184*fr[29]*nul[103])+0.0883883476483184*fl[29]*nul[103]+0.07654655446197427*fr[21]*nul[103]+0.07654655446197427*fl[21]*nul[103]-0.0883883476483184*fr[25]*nul[99]+0.0883883476483184*fl[25]*nul[99]+0.07654655446197427*fr[14]*nul[99]+0.07654655446197427*fl[14]*nul[99]-0.0883883476483184*fr[23]*nul[97]+0.0883883476483184*fl[23]*nul[97]+0.07654655446197427*fr[12]*nul[97]+0.07654655446197427*fl[12]*nul[97]-0.0883883476483184*fr[15]*nul[96]+0.0883883476483184*fl[15]*nul[96]+0.07654655446197427*fr[5]*nul[96]+0.07654655446197427*fl[5]*nul[96]; 
  incr2[17] = (-0.0883883476483184*fr[19]*nul[103])+0.0883883476483184*fl[19]*nul[103]+0.07654655446197427*fr[8]*nul[103]+0.07654655446197427*fl[8]*nul[103]-0.0883883476483184*fr[26]*nul[99]+0.0883883476483184*fl[26]*nul[99]+0.07654655446197427*fr[16]*nul[99]+0.07654655446197427*fl[16]*nul[99]-0.0883883476483184*fr[10]*nul[97]+0.0883883476483184*fl[10]*nul[97]+0.07654655446197427*fr[2]*nul[97]+0.07654655446197427*fl[2]*nul[97]-0.0883883476483184*fr[17]*nul[96]+0.0883883476483184*fl[17]*nul[96]+0.07654655446197427*fr[6]*nul[96]+0.07654655446197427*fl[6]*nul[96]; 
  incr2[18] = (-0.0883883476483184*fr[4]*nul[103])+0.0883883476483184*fl[4]*nul[103]+0.07654655446197427*fr[0]*nul[103]+0.07654655446197427*fl[0]*nul[103]-0.0883883476483184*fr[9]*nul[99]+0.0883883476483184*fl[9]*nul[99]+0.07654655446197427*fr[1]*nul[99]+0.07654655446197427*fl[1]*nul[99]-0.0883883476483184*fr[11]*nul[97]+0.0883883476483184*fl[11]*nul[97]+0.07654655446197427*fr[3]*nul[97]+0.07654655446197427*fl[3]*nul[97]-0.0883883476483184*fr[18]*nul[96]+0.0883883476483184*fl[18]*nul[96]+0.07654655446197427*fr[7]*nul[96]+0.07654655446197427*fl[7]*nul[96]; 
  incr2[19] = (-0.0883883476483184*fr[17]*nul[103])+0.0883883476483184*fl[17]*nul[103]+0.07654655446197427*fr[6]*nul[103]+0.07654655446197427*fl[6]*nul[103]-0.0883883476483184*fr[10]*nul[99]+0.0883883476483184*fl[10]*nul[99]+0.07654655446197427*fr[2]*nul[99]+0.07654655446197427*fl[2]*nul[99]-0.0883883476483184*fr[26]*nul[97]+0.0883883476483184*fl[26]*nul[97]+0.07654655446197427*fr[16]*nul[97]+0.07654655446197427*fl[16]*nul[97]-0.0883883476483184*fr[19]*nul[96]+0.0883883476483184*fl[19]*nul[96]+0.07654655446197427*fr[8]*nul[96]+0.07654655446197427*fl[8]*nul[96]; 
  incr2[23] = (-0.0883883476483184*fr[25]*nul[103])+0.0883883476483184*fl[25]*nul[103]+0.07654655446197427*fr[14]*nul[103]+0.07654655446197427*fl[14]*nul[103]-0.0883883476483184*fr[29]*nul[99]+0.0883883476483184*fl[29]*nul[99]+0.07654655446197427*fr[21]*nul[99]+0.07654655446197427*fl[21]*nul[99]-0.0883883476483184*fr[15]*nul[97]+0.0883883476483184*fl[15]*nul[97]+0.07654655446197427*fr[5]*nul[97]+0.07654655446197427*fl[5]*nul[97]-0.0883883476483184*fr[23]*nul[96]+0.0883883476483184*fl[23]*nul[96]+0.07654655446197427*fr[12]*nul[96]+0.07654655446197427*fl[12]*nul[96]; 
  incr2[24] = (-0.0883883476483184*fr[31]*nul[103])+0.0883883476483184*fl[31]*nul[103]+0.07654655446197427*fr[27]*nul[103]+0.07654655446197427*fl[27]*nul[103]-0.0883883476483184*fr[30]*nul[99]+0.0883883476483184*fl[30]*nul[99]+0.07654655446197427*fr[22]*nul[99]+0.07654655446197427*fl[22]*nul[99]-0.0883883476483184*fr[28]*nul[97]+0.0883883476483184*fl[28]*nul[97]+0.07654655446197427*fr[20]*nul[97]+0.07654655446197427*fl[20]*nul[97]-0.0883883476483184*fr[24]*nul[96]+0.0883883476483184*fl[24]*nul[96]+0.07654655446197427*fr[13]*nul[96]+0.07654655446197427*fl[13]*nul[96]; 
  incr2[25] = (-0.0883883476483184*fr[23]*nul[103])+0.0883883476483184*fl[23]*nul[103]+0.07654655446197427*fr[12]*nul[103]+0.07654655446197427*fl[12]*nul[103]-0.0883883476483184*fr[15]*nul[99]+0.0883883476483184*fl[15]*nul[99]+0.07654655446197427*fr[5]*nul[99]+0.07654655446197427*fl[5]*nul[99]-0.0883883476483184*fr[29]*nul[97]+0.0883883476483184*fl[29]*nul[97]+0.07654655446197427*fr[21]*nul[97]+0.07654655446197427*fl[21]*nul[97]-0.0883883476483184*fr[25]*nul[96]+0.0883883476483184*fl[25]*nul[96]+0.07654655446197427*fr[14]*nul[96]+0.07654655446197427*fl[14]*nul[96]; 
  incr2[26] = (-0.0883883476483184*fr[10]*nul[103])+0.0883883476483184*fl[10]*nul[103]+0.07654655446197427*fr[2]*nul[103]+0.07654655446197427*fl[2]*nul[103]-0.0883883476483184*fr[17]*nul[99]+0.0883883476483184*fl[17]*nul[99]+0.07654655446197427*fr[6]*nul[99]+0.07654655446197427*fl[6]*nul[99]-0.0883883476483184*fr[19]*nul[97]+0.0883883476483184*fl[19]*nul[97]+0.07654655446197427*fr[8]*nul[97]+0.07654655446197427*fl[8]*nul[97]-0.0883883476483184*fr[26]*nul[96]+0.0883883476483184*fl[26]*nul[96]+0.07654655446197427*fr[16]*nul[96]+0.07654655446197427*fl[16]*nul[96]; 
  incr2[28] = (-0.0883883476483184*fr[30]*nul[103])+0.0883883476483184*fl[30]*nul[103]+0.07654655446197427*fr[22]*nul[103]+0.07654655446197427*fl[22]*nul[103]-0.0883883476483184*fr[31]*nul[99]+0.0883883476483184*fl[31]*nul[99]+0.07654655446197427*fr[27]*nul[99]+0.07654655446197427*fl[27]*nul[99]-0.0883883476483184*fr[24]*nul[97]+0.0883883476483184*fl[24]*nul[97]+0.07654655446197427*fr[13]*nul[97]+0.07654655446197427*fl[13]*nul[97]-0.0883883476483184*fr[28]*nul[96]+0.0883883476483184*fl[28]*nul[96]+0.07654655446197427*fr[20]*nul[96]+0.07654655446197427*fl[20]*nul[96]; 
  incr2[29] = (-0.0883883476483184*fr[15]*nul[103])+0.0883883476483184*fl[15]*nul[103]+0.07654655446197427*fr[5]*nul[103]+0.07654655446197427*fl[5]*nul[103]-0.0883883476483184*fr[23]*nul[99]+0.0883883476483184*fl[23]*nul[99]+0.07654655446197427*fr[12]*nul[99]+0.07654655446197427*fl[12]*nul[99]-0.0883883476483184*fr[25]*nul[97]+0.0883883476483184*fl[25]*nul[97]+0.07654655446197427*fr[14]*nul[97]+0.07654655446197427*fl[14]*nul[97]-0.0883883476483184*fr[29]*nul[96]+0.0883883476483184*fl[29]*nul[96]+0.07654655446197427*fr[21]*nul[96]+0.07654655446197427*fl[21]*nul[96]; 
  incr2[30] = (-0.0883883476483184*fr[28]*nul[103])+0.0883883476483184*fl[28]*nul[103]+0.07654655446197427*fr[20]*nul[103]+0.07654655446197427*fl[20]*nul[103]-0.0883883476483184*fr[24]*nul[99]+0.0883883476483184*fl[24]*nul[99]+0.07654655446197427*fr[13]*nul[99]+0.07654655446197427*fl[13]*nul[99]-0.0883883476483184*fr[31]*nul[97]+0.0883883476483184*fl[31]*nul[97]+0.07654655446197427*fr[27]*nul[97]+0.07654655446197427*fl[27]*nul[97]-0.0883883476483184*fr[30]*nul[96]+0.0883883476483184*fl[30]*nul[96]+0.07654655446197427*fr[22]*nul[96]+0.07654655446197427*fl[22]*nul[96]; 
  incr2[31] = (-0.0883883476483184*fr[24]*nul[103])+0.0883883476483184*fl[24]*nul[103]+0.07654655446197427*fr[13]*nul[103]+0.07654655446197427*fl[13]*nul[103]-0.0883883476483184*fr[28]*nul[99]+0.0883883476483184*fl[28]*nul[99]+0.07654655446197427*fr[20]*nul[99]+0.07654655446197427*fl[20]*nul[99]-0.0883883476483184*fr[30]*nul[97]+0.0883883476483184*fl[30]*nul[97]+0.07654655446197427*fr[22]*nul[97]+0.07654655446197427*fl[22]*nul[97]-0.0883883476483184*fr[31]*nul[96]+0.0883883476483184*fl[31]*nul[96]+0.07654655446197427*fr[27]*nul[96]+0.07654655446197427*fl[27]*nul[96]; 

  outr[0] += incr1[0]*rdxFr; 
  outr[1] += incr1[1]*rdxFr; 
  outr[2] += incr1[2]*rdxFr; 
  outr[3] += incr1[3]*rdxFr; 
  outr[4] += incr2[4]*rdxFr+incr1[4]*rdxFr; 
  outr[5] += incr1[5]*rdxFr; 
  outr[6] += incr1[6]*rdxFr; 
  outr[7] += incr1[7]*rdxFr; 
  outr[8] += incr1[8]*rdxFr; 
  outr[9] += incr2[9]*rdxFr+incr1[9]*rdxFr; 
  outr[10] += incr2[10]*rdxFr+incr1[10]*rdxFr; 
  outr[11] += incr2[11]*rdxFr+incr1[11]*rdxFr; 
  outr[12] += incr1[12]*rdxFr; 
  outr[13] += incr1[13]*rdxFr; 
  outr[14] += incr1[14]*rdxFr; 
  outr[15] += incr2[15]*rdxFr+incr1[15]*rdxFr; 
  outr[16] += incr1[16]*rdxFr; 
  outr[17] += incr2[17]*rdxFr+incr1[17]*rdxFr; 
  outr[18] += incr2[18]*rdxFr+incr1[18]*rdxFr; 
  outr[19] += incr2[19]*rdxFr+incr1[19]*rdxFr; 
  outr[20] += incr1[20]*rdxFr; 
  outr[21] += incr1[21]*rdxFr; 
  outr[22] += incr1[22]*rdxFr; 
  outr[23] += incr2[23]*rdxFr+incr1[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr+incr1[24]*rdxFr; 
  outr[25] += incr2[25]*rdxFr+incr1[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr+incr1[26]*rdxFr; 
  outr[27] += incr1[27]*rdxFr; 
  outr[28] += incr2[28]*rdxFr+incr1[28]*rdxFr; 
  outr[29] += incr2[29]*rdxFr+incr1[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr+incr1[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr+incr1[31]*rdxFr; 

  outl[0] += -1.0*incr1[0]*rdxFl; 
  outl[1] += -1.0*incr1[1]*rdxFl; 
  outl[2] += -1.0*incr1[2]*rdxFl; 
  outl[3] += -1.0*incr1[3]*rdxFl; 
  outl[4] += incr1[4]*rdxFl-1.0*incr2[4]*rdxFl; 
  outl[5] += -1.0*incr1[5]*rdxFl; 
  outl[6] += -1.0*incr1[6]*rdxFl; 
  outl[7] += -1.0*incr1[7]*rdxFl; 
  outl[8] += -1.0*incr1[8]*rdxFl; 
  outl[9] += incr1[9]*rdxFl-1.0*incr2[9]*rdxFl; 
  outl[10] += incr1[10]*rdxFl-1.0*incr2[10]*rdxFl; 
  outl[11] += incr1[11]*rdxFl-1.0*incr2[11]*rdxFl; 
  outl[12] += -1.0*incr1[12]*rdxFl; 
  outl[13] += -1.0*incr1[13]*rdxFl; 
  outl[14] += -1.0*incr1[14]*rdxFl; 
  outl[15] += incr1[15]*rdxFl-1.0*incr2[15]*rdxFl; 
  outl[16] += -1.0*incr1[16]*rdxFl; 
  outl[17] += incr1[17]*rdxFl-1.0*incr2[17]*rdxFl; 
  outl[18] += incr1[18]*rdxFl-1.0*incr2[18]*rdxFl; 
  outl[19] += incr1[19]*rdxFl-1.0*incr2[19]*rdxFl; 
  outl[20] += -1.0*incr1[20]*rdxFl; 
  outl[21] += -1.0*incr1[21]*rdxFl; 
  outl[22] += -1.0*incr1[22]*rdxFl; 
  outl[23] += incr1[23]*rdxFl-1.0*incr2[23]*rdxFl; 
  outl[24] += incr1[24]*rdxFl-1.0*incr2[24]*rdxFl; 
  outl[25] += incr1[25]*rdxFl-1.0*incr2[25]*rdxFl; 
  outl[26] += incr1[26]*rdxFl-1.0*incr2[26]*rdxFl; 
  outl[27] += -1.0*incr1[27]*rdxFl; 
  outl[28] += incr1[28]*rdxFl-1.0*incr2[28]*rdxFl; 
  outl[29] += incr1[29]*rdxFl-1.0*incr2[29]*rdxFl; 
  outl[30] += incr1[30]*rdxFl-1.0*incr2[30]*rdxFl; 
  outl[31] += incr1[31]*rdxFl-1.0*incr2[31]*rdxFl; 

} 
void ConstDiffusionVarCoeffSurf5xSerP1_X5(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[5]:      Cell-center coordinates.
  // dx[5]:     Cell spacing.
  // nu[160]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[4]*dxl[4]); 
  double rdxFr = 4.0/(dxr[4]*dxr[4]); 

  double incr1[32]; 
  double incr2[32]; 

  incr1[0] = 0.09568319307746781*fr[21]*nul[135]+0.09568319307746781*fl[21]*nul[135]-0.09943689110435816*fr[7]*nul[135]+0.09943689110435816*fl[7]*nul[135]+0.09568319307746781*fr[14]*nul[131]+0.09568319307746781*fl[14]*nul[131]-0.09943689110435816*fr[3]*nul[131]+0.09943689110435816*fl[3]*nul[131]+0.09568319307746781*fr[12]*nul[129]+0.09568319307746781*fl[12]*nul[129]-0.09943689110435816*fr[1]*nul[129]+0.09943689110435816*fl[1]*nul[129]+0.09568319307746781*fr[5]*nul[128]+0.09568319307746781*fl[5]*nul[128]-0.09943689110435816*fr[0]*nul[128]+0.09943689110435816*fl[0]*nul[128]; 
  incr1[1] = 0.09568319307746781*fr[14]*nul[135]+0.09568319307746781*fl[14]*nul[135]-0.09943689110435816*fr[3]*nul[135]+0.09943689110435816*fl[3]*nul[135]+0.09568319307746781*fr[21]*nul[131]+0.09568319307746781*fl[21]*nul[131]-0.09943689110435816*fr[7]*nul[131]+0.09943689110435816*fl[7]*nul[131]+0.09568319307746781*fr[5]*nul[129]+0.09568319307746781*fl[5]*nul[129]-0.09943689110435816*fr[0]*nul[129]+0.09943689110435816*fl[0]*nul[129]+0.09568319307746781*fr[12]*nul[128]+0.09568319307746781*fl[12]*nul[128]-0.09943689110435816*fr[1]*nul[128]+0.09943689110435816*fl[1]*nul[128]; 
  incr1[2] = 0.09568319307746781*fr[27]*nul[135]+0.09568319307746781*fl[27]*nul[135]-0.09943689110435816*fr[16]*nul[135]+0.09943689110435816*fl[16]*nul[135]+0.09568319307746781*fr[22]*nul[131]+0.09568319307746781*fl[22]*nul[131]-0.09943689110435816*fr[8]*nul[131]+0.09943689110435816*fl[8]*nul[131]+0.09568319307746781*fr[20]*nul[129]+0.09568319307746781*fl[20]*nul[129]-0.09943689110435816*fr[6]*nul[129]+0.09943689110435816*fl[6]*nul[129]+0.09568319307746781*fr[13]*nul[128]+0.09568319307746781*fl[13]*nul[128]-0.09943689110435816*fr[2]*nul[128]+0.09943689110435816*fl[2]*nul[128]; 
  incr1[3] = 0.09568319307746781*fr[12]*nul[135]+0.09568319307746781*fl[12]*nul[135]-0.09943689110435816*fr[1]*nul[135]+0.09943689110435816*fl[1]*nul[135]+0.09568319307746781*fr[5]*nul[131]+0.09568319307746781*fl[5]*nul[131]-0.09943689110435816*fr[0]*nul[131]+0.09943689110435816*fl[0]*nul[131]+0.09568319307746781*fr[21]*nul[129]+0.09568319307746781*fl[21]*nul[129]-0.09943689110435816*fr[7]*nul[129]+0.09943689110435816*fl[7]*nul[129]+0.09568319307746781*fr[14]*nul[128]+0.09568319307746781*fl[14]*nul[128]-0.09943689110435816*fr[3]*nul[128]+0.09943689110435816*fl[3]*nul[128]; 
  incr1[4] = 0.09568319307746781*fr[29]*nul[135]+0.09568319307746781*fl[29]*nul[135]-0.09943689110435816*fr[18]*nul[135]+0.09943689110435816*fl[18]*nul[135]+0.09568319307746781*fr[25]*nul[131]+0.09568319307746781*fl[25]*nul[131]-0.09943689110435816*fr[11]*nul[131]+0.09943689110435816*fl[11]*nul[131]+0.09568319307746781*fr[23]*nul[129]+0.09568319307746781*fl[23]*nul[129]-0.09943689110435816*fr[9]*nul[129]+0.09943689110435816*fl[9]*nul[129]+0.09568319307746781*fr[15]*nul[128]+0.09568319307746781*fl[15]*nul[128]-0.09943689110435816*fr[4]*nul[128]+0.09943689110435816*fl[4]*nul[128]; 
  incr1[5] = (-0.1657281518405969*fr[21]*nul[135])-0.1657281518405969*fl[21]*nul[135]+0.172229747539442*fr[7]*nul[135]-0.172229747539442*fl[7]*nul[135]-0.1657281518405969*fr[14]*nul[131]-0.1657281518405969*fl[14]*nul[131]+0.172229747539442*fr[3]*nul[131]-0.172229747539442*fl[3]*nul[131]-0.1657281518405969*fr[12]*nul[129]-0.1657281518405969*fl[12]*nul[129]+0.172229747539442*fr[1]*nul[129]-0.172229747539442*fl[1]*nul[129]-0.1657281518405969*fr[5]*nul[128]-0.1657281518405969*fl[5]*nul[128]+0.172229747539442*fr[0]*nul[128]-0.172229747539442*fl[0]*nul[128]; 
  incr1[6] = 0.09568319307746781*fr[22]*nul[135]+0.09568319307746781*fl[22]*nul[135]-0.09943689110435816*fr[8]*nul[135]+0.09943689110435816*fl[8]*nul[135]+0.09568319307746781*fr[27]*nul[131]+0.09568319307746781*fl[27]*nul[131]-0.09943689110435816*fr[16]*nul[131]+0.09943689110435816*fl[16]*nul[131]+0.09568319307746781*fr[13]*nul[129]+0.09568319307746781*fl[13]*nul[129]-0.09943689110435816*fr[2]*nul[129]+0.09943689110435816*fl[2]*nul[129]+0.09568319307746781*fr[20]*nul[128]+0.09568319307746781*fl[20]*nul[128]-0.09943689110435816*fr[6]*nul[128]+0.09943689110435816*fl[6]*nul[128]; 
  incr1[7] = 0.09568319307746781*fr[5]*nul[135]+0.09568319307746781*fl[5]*nul[135]-0.09943689110435816*fr[0]*nul[135]+0.09943689110435816*fl[0]*nul[135]+0.09568319307746781*fr[12]*nul[131]+0.09568319307746781*fl[12]*nul[131]-0.09943689110435816*fr[1]*nul[131]+0.09943689110435816*fl[1]*nul[131]+0.09568319307746781*fr[14]*nul[129]+0.09568319307746781*fl[14]*nul[129]-0.09943689110435816*fr[3]*nul[129]+0.09943689110435816*fl[3]*nul[129]+0.09568319307746781*fr[21]*nul[128]+0.09568319307746781*fl[21]*nul[128]-0.09943689110435816*fr[7]*nul[128]+0.09943689110435816*fl[7]*nul[128]; 
  incr1[8] = 0.09568319307746781*fr[20]*nul[135]+0.09568319307746781*fl[20]*nul[135]-0.09943689110435816*fr[6]*nul[135]+0.09943689110435816*fl[6]*nul[135]+0.09568319307746781*fr[13]*nul[131]+0.09568319307746781*fl[13]*nul[131]-0.09943689110435816*fr[2]*nul[131]+0.09943689110435816*fl[2]*nul[131]+0.09568319307746781*fr[27]*nul[129]+0.09568319307746781*fl[27]*nul[129]-0.09943689110435816*fr[16]*nul[129]+0.09943689110435816*fl[16]*nul[129]+0.09568319307746781*fr[22]*nul[128]+0.09568319307746781*fl[22]*nul[128]-0.09943689110435816*fr[8]*nul[128]+0.09943689110435816*fl[8]*nul[128]; 
  incr1[9] = 0.09568319307746781*fr[25]*nul[135]+0.09568319307746781*fl[25]*nul[135]-0.09943689110435816*fr[11]*nul[135]+0.09943689110435816*fl[11]*nul[135]+0.09568319307746781*fr[29]*nul[131]+0.09568319307746781*fl[29]*nul[131]-0.09943689110435816*fr[18]*nul[131]+0.09943689110435816*fl[18]*nul[131]+0.09568319307746781*fr[15]*nul[129]+0.09568319307746781*fl[15]*nul[129]-0.09943689110435816*fr[4]*nul[129]+0.09943689110435816*fl[4]*nul[129]+0.09568319307746781*fr[23]*nul[128]+0.09568319307746781*fl[23]*nul[128]-0.09943689110435816*fr[9]*nul[128]+0.09943689110435816*fl[9]*nul[128]; 
  incr1[10] = 0.09568319307746781*fr[31]*nul[135]+0.09568319307746781*fl[31]*nul[135]-0.09943689110435816*fr[26]*nul[135]+0.09943689110435816*fl[26]*nul[135]+0.09568319307746781*fr[30]*nul[131]+0.09568319307746781*fl[30]*nul[131]-0.09943689110435816*fr[19]*nul[131]+0.09943689110435816*fl[19]*nul[131]+0.09568319307746781*fr[28]*nul[129]+0.09568319307746781*fl[28]*nul[129]-0.09943689110435816*fr[17]*nul[129]+0.09943689110435816*fl[17]*nul[129]+0.09568319307746781*fr[24]*nul[128]+0.09568319307746781*fl[24]*nul[128]-0.09943689110435816*fr[10]*nul[128]+0.09943689110435816*fl[10]*nul[128]; 
  incr1[11] = 0.09568319307746781*fr[23]*nul[135]+0.09568319307746781*fl[23]*nul[135]-0.09943689110435816*fr[9]*nul[135]+0.09943689110435816*fl[9]*nul[135]+0.09568319307746781*fr[15]*nul[131]+0.09568319307746781*fl[15]*nul[131]-0.09943689110435816*fr[4]*nul[131]+0.09943689110435816*fl[4]*nul[131]+0.09568319307746781*fr[29]*nul[129]+0.09568319307746781*fl[29]*nul[129]-0.09943689110435816*fr[18]*nul[129]+0.09943689110435816*fl[18]*nul[129]+0.09568319307746781*fr[25]*nul[128]+0.09568319307746781*fl[25]*nul[128]-0.09943689110435816*fr[11]*nul[128]+0.09943689110435816*fl[11]*nul[128]; 
  incr1[12] = (-0.1657281518405969*fr[14]*nul[135])-0.1657281518405969*fl[14]*nul[135]+0.172229747539442*fr[3]*nul[135]-0.172229747539442*fl[3]*nul[135]-0.1657281518405969*fr[21]*nul[131]-0.1657281518405969*fl[21]*nul[131]+0.172229747539442*fr[7]*nul[131]-0.172229747539442*fl[7]*nul[131]-0.1657281518405969*fr[5]*nul[129]-0.1657281518405969*fl[5]*nul[129]+0.172229747539442*fr[0]*nul[129]-0.172229747539442*fl[0]*nul[129]-0.1657281518405969*fr[12]*nul[128]-0.1657281518405969*fl[12]*nul[128]+0.172229747539442*fr[1]*nul[128]-0.172229747539442*fl[1]*nul[128]; 
  incr1[13] = (-0.1657281518405969*fr[27]*nul[135])-0.1657281518405969*fl[27]*nul[135]+0.172229747539442*fr[16]*nul[135]-0.172229747539442*fl[16]*nul[135]-0.1657281518405969*fr[22]*nul[131]-0.1657281518405969*fl[22]*nul[131]+0.172229747539442*fr[8]*nul[131]-0.172229747539442*fl[8]*nul[131]-0.1657281518405969*fr[20]*nul[129]-0.1657281518405969*fl[20]*nul[129]+0.172229747539442*fr[6]*nul[129]-0.172229747539442*fl[6]*nul[129]-0.1657281518405969*fr[13]*nul[128]-0.1657281518405969*fl[13]*nul[128]+0.172229747539442*fr[2]*nul[128]-0.172229747539442*fl[2]*nul[128]; 
  incr1[14] = (-0.1657281518405969*fr[12]*nul[135])-0.1657281518405969*fl[12]*nul[135]+0.172229747539442*fr[1]*nul[135]-0.172229747539442*fl[1]*nul[135]-0.1657281518405969*fr[5]*nul[131]-0.1657281518405969*fl[5]*nul[131]+0.172229747539442*fr[0]*nul[131]-0.172229747539442*fl[0]*nul[131]-0.1657281518405969*fr[21]*nul[129]-0.1657281518405969*fl[21]*nul[129]+0.172229747539442*fr[7]*nul[129]-0.172229747539442*fl[7]*nul[129]-0.1657281518405969*fr[14]*nul[128]-0.1657281518405969*fl[14]*nul[128]+0.172229747539442*fr[3]*nul[128]-0.172229747539442*fl[3]*nul[128]; 
  incr1[15] = (-0.1657281518405969*fr[29]*nul[135])-0.1657281518405969*fl[29]*nul[135]+0.172229747539442*fr[18]*nul[135]-0.172229747539442*fl[18]*nul[135]-0.1657281518405969*fr[25]*nul[131]-0.1657281518405969*fl[25]*nul[131]+0.172229747539442*fr[11]*nul[131]-0.172229747539442*fl[11]*nul[131]-0.1657281518405969*fr[23]*nul[129]-0.1657281518405969*fl[23]*nul[129]+0.172229747539442*fr[9]*nul[129]-0.172229747539442*fl[9]*nul[129]-0.1657281518405969*fr[15]*nul[128]-0.1657281518405969*fl[15]*nul[128]+0.172229747539442*fr[4]*nul[128]-0.172229747539442*fl[4]*nul[128]; 
  incr1[16] = 0.09568319307746781*fr[13]*nul[135]+0.09568319307746781*fl[13]*nul[135]-0.09943689110435816*fr[2]*nul[135]+0.09943689110435816*fl[2]*nul[135]+0.09568319307746781*fr[20]*nul[131]+0.09568319307746781*fl[20]*nul[131]-0.09943689110435816*fr[6]*nul[131]+0.09943689110435816*fl[6]*nul[131]+0.09568319307746781*fr[22]*nul[129]+0.09568319307746781*fl[22]*nul[129]-0.09943689110435816*fr[8]*nul[129]+0.09943689110435816*fl[8]*nul[129]+0.09568319307746781*fr[27]*nul[128]+0.09568319307746781*fl[27]*nul[128]-0.09943689110435816*fr[16]*nul[128]+0.09943689110435816*fl[16]*nul[128]; 
  incr1[17] = 0.09568319307746781*fr[30]*nul[135]+0.09568319307746781*fl[30]*nul[135]-0.09943689110435816*fr[19]*nul[135]+0.09943689110435816*fl[19]*nul[135]+0.09568319307746781*fr[31]*nul[131]+0.09568319307746781*fl[31]*nul[131]-0.09943689110435816*fr[26]*nul[131]+0.09943689110435816*fl[26]*nul[131]+0.09568319307746781*fr[24]*nul[129]+0.09568319307746781*fl[24]*nul[129]-0.09943689110435816*fr[10]*nul[129]+0.09943689110435816*fl[10]*nul[129]+0.09568319307746781*fr[28]*nul[128]+0.09568319307746781*fl[28]*nul[128]-0.09943689110435816*fr[17]*nul[128]+0.09943689110435816*fl[17]*nul[128]; 
  incr1[18] = 0.09568319307746781*fr[15]*nul[135]+0.09568319307746781*fl[15]*nul[135]-0.09943689110435816*fr[4]*nul[135]+0.09943689110435816*fl[4]*nul[135]+0.09568319307746781*fr[23]*nul[131]+0.09568319307746781*fl[23]*nul[131]-0.09943689110435816*fr[9]*nul[131]+0.09943689110435816*fl[9]*nul[131]+0.09568319307746781*fr[25]*nul[129]+0.09568319307746781*fl[25]*nul[129]-0.09943689110435816*fr[11]*nul[129]+0.09943689110435816*fl[11]*nul[129]+0.09568319307746781*fr[29]*nul[128]+0.09568319307746781*fl[29]*nul[128]-0.09943689110435816*fr[18]*nul[128]+0.09943689110435816*fl[18]*nul[128]; 
  incr1[19] = 0.09568319307746781*fr[28]*nul[135]+0.09568319307746781*fl[28]*nul[135]-0.09943689110435816*fr[17]*nul[135]+0.09943689110435816*fl[17]*nul[135]+0.09568319307746781*fr[24]*nul[131]+0.09568319307746781*fl[24]*nul[131]-0.09943689110435816*fr[10]*nul[131]+0.09943689110435816*fl[10]*nul[131]+0.09568319307746781*fr[31]*nul[129]+0.09568319307746781*fl[31]*nul[129]-0.09943689110435816*fr[26]*nul[129]+0.09943689110435816*fl[26]*nul[129]+0.09568319307746781*fr[30]*nul[128]+0.09568319307746781*fl[30]*nul[128]-0.09943689110435816*fr[19]*nul[128]+0.09943689110435816*fl[19]*nul[128]; 
  incr1[20] = (-0.1657281518405969*fr[22]*nul[135])-0.1657281518405969*fl[22]*nul[135]+0.172229747539442*fr[8]*nul[135]-0.172229747539442*fl[8]*nul[135]-0.1657281518405969*fr[27]*nul[131]-0.1657281518405969*fl[27]*nul[131]+0.172229747539442*fr[16]*nul[131]-0.172229747539442*fl[16]*nul[131]-0.1657281518405969*fr[13]*nul[129]-0.1657281518405969*fl[13]*nul[129]+0.172229747539442*fr[2]*nul[129]-0.172229747539442*fl[2]*nul[129]-0.1657281518405969*fr[20]*nul[128]-0.1657281518405969*fl[20]*nul[128]+0.172229747539442*fr[6]*nul[128]-0.172229747539442*fl[6]*nul[128]; 
  incr1[21] = (-0.1657281518405969*fr[5]*nul[135])-0.1657281518405969*fl[5]*nul[135]+0.172229747539442*fr[0]*nul[135]-0.172229747539442*fl[0]*nul[135]-0.1657281518405969*fr[12]*nul[131]-0.1657281518405969*fl[12]*nul[131]+0.172229747539442*fr[1]*nul[131]-0.172229747539442*fl[1]*nul[131]-0.1657281518405969*fr[14]*nul[129]-0.1657281518405969*fl[14]*nul[129]+0.172229747539442*fr[3]*nul[129]-0.172229747539442*fl[3]*nul[129]-0.1657281518405969*fr[21]*nul[128]-0.1657281518405969*fl[21]*nul[128]+0.172229747539442*fr[7]*nul[128]-0.172229747539442*fl[7]*nul[128]; 
  incr1[22] = (-0.1657281518405969*fr[20]*nul[135])-0.1657281518405969*fl[20]*nul[135]+0.172229747539442*fr[6]*nul[135]-0.172229747539442*fl[6]*nul[135]-0.1657281518405969*fr[13]*nul[131]-0.1657281518405969*fl[13]*nul[131]+0.172229747539442*fr[2]*nul[131]-0.172229747539442*fl[2]*nul[131]-0.1657281518405969*fr[27]*nul[129]-0.1657281518405969*fl[27]*nul[129]+0.172229747539442*fr[16]*nul[129]-0.172229747539442*fl[16]*nul[129]-0.1657281518405969*fr[22]*nul[128]-0.1657281518405969*fl[22]*nul[128]+0.172229747539442*fr[8]*nul[128]-0.172229747539442*fl[8]*nul[128]; 
  incr1[23] = (-0.1657281518405969*fr[25]*nul[135])-0.1657281518405969*fl[25]*nul[135]+0.172229747539442*fr[11]*nul[135]-0.172229747539442*fl[11]*nul[135]-0.1657281518405969*fr[29]*nul[131]-0.1657281518405969*fl[29]*nul[131]+0.172229747539442*fr[18]*nul[131]-0.172229747539442*fl[18]*nul[131]-0.1657281518405969*fr[15]*nul[129]-0.1657281518405969*fl[15]*nul[129]+0.172229747539442*fr[4]*nul[129]-0.172229747539442*fl[4]*nul[129]-0.1657281518405969*fr[23]*nul[128]-0.1657281518405969*fl[23]*nul[128]+0.172229747539442*fr[9]*nul[128]-0.172229747539442*fl[9]*nul[128]; 
  incr1[24] = (-0.1657281518405969*fr[31]*nul[135])-0.1657281518405969*fl[31]*nul[135]+0.172229747539442*fr[26]*nul[135]-0.172229747539442*fl[26]*nul[135]-0.1657281518405969*fr[30]*nul[131]-0.1657281518405969*fl[30]*nul[131]+0.172229747539442*fr[19]*nul[131]-0.172229747539442*fl[19]*nul[131]-0.1657281518405969*fr[28]*nul[129]-0.1657281518405969*fl[28]*nul[129]+0.172229747539442*fr[17]*nul[129]-0.172229747539442*fl[17]*nul[129]-0.1657281518405969*fr[24]*nul[128]-0.1657281518405969*fl[24]*nul[128]+0.172229747539442*fr[10]*nul[128]-0.172229747539442*fl[10]*nul[128]; 
  incr1[25] = (-0.1657281518405969*fr[23]*nul[135])-0.1657281518405969*fl[23]*nul[135]+0.172229747539442*fr[9]*nul[135]-0.172229747539442*fl[9]*nul[135]-0.1657281518405969*fr[15]*nul[131]-0.1657281518405969*fl[15]*nul[131]+0.172229747539442*fr[4]*nul[131]-0.172229747539442*fl[4]*nul[131]-0.1657281518405969*fr[29]*nul[129]-0.1657281518405969*fl[29]*nul[129]+0.172229747539442*fr[18]*nul[129]-0.172229747539442*fl[18]*nul[129]-0.1657281518405969*fr[25]*nul[128]-0.1657281518405969*fl[25]*nul[128]+0.172229747539442*fr[11]*nul[128]-0.172229747539442*fl[11]*nul[128]; 
  incr1[26] = 0.09568319307746781*fr[24]*nul[135]+0.09568319307746781*fl[24]*nul[135]-0.09943689110435816*fr[10]*nul[135]+0.09943689110435816*fl[10]*nul[135]+0.09568319307746781*fr[28]*nul[131]+0.09568319307746781*fl[28]*nul[131]-0.09943689110435816*fr[17]*nul[131]+0.09943689110435816*fl[17]*nul[131]+0.09568319307746781*fr[30]*nul[129]+0.09568319307746781*fl[30]*nul[129]-0.09943689110435816*fr[19]*nul[129]+0.09943689110435816*fl[19]*nul[129]+0.09568319307746781*fr[31]*nul[128]+0.09568319307746781*fl[31]*nul[128]-0.09943689110435816*fr[26]*nul[128]+0.09943689110435816*fl[26]*nul[128]; 
  incr1[27] = (-0.1657281518405969*fr[13]*nul[135])-0.1657281518405969*fl[13]*nul[135]+0.172229747539442*fr[2]*nul[135]-0.172229747539442*fl[2]*nul[135]-0.1657281518405969*fr[20]*nul[131]-0.1657281518405969*fl[20]*nul[131]+0.172229747539442*fr[6]*nul[131]-0.172229747539442*fl[6]*nul[131]-0.1657281518405969*fr[22]*nul[129]-0.1657281518405969*fl[22]*nul[129]+0.172229747539442*fr[8]*nul[129]-0.172229747539442*fl[8]*nul[129]-0.1657281518405969*fr[27]*nul[128]-0.1657281518405969*fl[27]*nul[128]+0.172229747539442*fr[16]*nul[128]-0.172229747539442*fl[16]*nul[128]; 
  incr1[28] = (-0.1657281518405969*fr[30]*nul[135])-0.1657281518405969*fl[30]*nul[135]+0.172229747539442*fr[19]*nul[135]-0.172229747539442*fl[19]*nul[135]-0.1657281518405969*fr[31]*nul[131]-0.1657281518405969*fl[31]*nul[131]+0.172229747539442*fr[26]*nul[131]-0.172229747539442*fl[26]*nul[131]-0.1657281518405969*fr[24]*nul[129]-0.1657281518405969*fl[24]*nul[129]+0.172229747539442*fr[10]*nul[129]-0.172229747539442*fl[10]*nul[129]-0.1657281518405969*fr[28]*nul[128]-0.1657281518405969*fl[28]*nul[128]+0.172229747539442*fr[17]*nul[128]-0.172229747539442*fl[17]*nul[128]; 
  incr1[29] = (-0.1657281518405969*fr[15]*nul[135])-0.1657281518405969*fl[15]*nul[135]+0.172229747539442*fr[4]*nul[135]-0.172229747539442*fl[4]*nul[135]-0.1657281518405969*fr[23]*nul[131]-0.1657281518405969*fl[23]*nul[131]+0.172229747539442*fr[9]*nul[131]-0.172229747539442*fl[9]*nul[131]-0.1657281518405969*fr[25]*nul[129]-0.1657281518405969*fl[25]*nul[129]+0.172229747539442*fr[11]*nul[129]-0.172229747539442*fl[11]*nul[129]-0.1657281518405969*fr[29]*nul[128]-0.1657281518405969*fl[29]*nul[128]+0.172229747539442*fr[18]*nul[128]-0.172229747539442*fl[18]*nul[128]; 
  incr1[30] = (-0.1657281518405969*fr[28]*nul[135])-0.1657281518405969*fl[28]*nul[135]+0.172229747539442*fr[17]*nul[135]-0.172229747539442*fl[17]*nul[135]-0.1657281518405969*fr[24]*nul[131]-0.1657281518405969*fl[24]*nul[131]+0.172229747539442*fr[10]*nul[131]-0.172229747539442*fl[10]*nul[131]-0.1657281518405969*fr[31]*nul[129]-0.1657281518405969*fl[31]*nul[129]+0.172229747539442*fr[26]*nul[129]-0.172229747539442*fl[26]*nul[129]-0.1657281518405969*fr[30]*nul[128]-0.1657281518405969*fl[30]*nul[128]+0.172229747539442*fr[19]*nul[128]-0.172229747539442*fl[19]*nul[128]; 
  incr1[31] = (-0.1657281518405969*fr[24]*nul[135])-0.1657281518405969*fl[24]*nul[135]+0.172229747539442*fr[10]*nul[135]-0.172229747539442*fl[10]*nul[135]-0.1657281518405969*fr[28]*nul[131]-0.1657281518405969*fl[28]*nul[131]+0.172229747539442*fr[17]*nul[131]-0.172229747539442*fl[17]*nul[131]-0.1657281518405969*fr[30]*nul[129]-0.1657281518405969*fl[30]*nul[129]+0.172229747539442*fr[19]*nul[129]-0.172229747539442*fl[19]*nul[129]-0.1657281518405969*fr[31]*nul[128]-0.1657281518405969*fl[31]*nul[128]+0.172229747539442*fr[26]*nul[128]-0.172229747539442*fl[26]*nul[128]; 

  incr2[5] = (-0.0883883476483184*fr[21]*nul[135])+0.0883883476483184*fl[21]*nul[135]+0.07654655446197427*fr[7]*nul[135]+0.07654655446197427*fl[7]*nul[135]-0.0883883476483184*fr[14]*nul[131]+0.0883883476483184*fl[14]*nul[131]+0.07654655446197427*fr[3]*nul[131]+0.07654655446197427*fl[3]*nul[131]-0.0883883476483184*fr[12]*nul[129]+0.0883883476483184*fl[12]*nul[129]+0.07654655446197427*fr[1]*nul[129]+0.07654655446197427*fl[1]*nul[129]-0.0883883476483184*fr[5]*nul[128]+0.0883883476483184*fl[5]*nul[128]+0.07654655446197427*fr[0]*nul[128]+0.07654655446197427*fl[0]*nul[128]; 
  incr2[12] = (-0.0883883476483184*fr[14]*nul[135])+0.0883883476483184*fl[14]*nul[135]+0.07654655446197427*fr[3]*nul[135]+0.07654655446197427*fl[3]*nul[135]-0.0883883476483184*fr[21]*nul[131]+0.0883883476483184*fl[21]*nul[131]+0.07654655446197427*fr[7]*nul[131]+0.07654655446197427*fl[7]*nul[131]-0.0883883476483184*fr[5]*nul[129]+0.0883883476483184*fl[5]*nul[129]+0.07654655446197427*fr[0]*nul[129]+0.07654655446197427*fl[0]*nul[129]-0.0883883476483184*fr[12]*nul[128]+0.0883883476483184*fl[12]*nul[128]+0.07654655446197427*fr[1]*nul[128]+0.07654655446197427*fl[1]*nul[128]; 
  incr2[13] = (-0.0883883476483184*fr[27]*nul[135])+0.0883883476483184*fl[27]*nul[135]+0.07654655446197427*fr[16]*nul[135]+0.07654655446197427*fl[16]*nul[135]-0.0883883476483184*fr[22]*nul[131]+0.0883883476483184*fl[22]*nul[131]+0.07654655446197427*fr[8]*nul[131]+0.07654655446197427*fl[8]*nul[131]-0.0883883476483184*fr[20]*nul[129]+0.0883883476483184*fl[20]*nul[129]+0.07654655446197427*fr[6]*nul[129]+0.07654655446197427*fl[6]*nul[129]-0.0883883476483184*fr[13]*nul[128]+0.0883883476483184*fl[13]*nul[128]+0.07654655446197427*fr[2]*nul[128]+0.07654655446197427*fl[2]*nul[128]; 
  incr2[14] = (-0.0883883476483184*fr[12]*nul[135])+0.0883883476483184*fl[12]*nul[135]+0.07654655446197427*fr[1]*nul[135]+0.07654655446197427*fl[1]*nul[135]-0.0883883476483184*fr[5]*nul[131]+0.0883883476483184*fl[5]*nul[131]+0.07654655446197427*fr[0]*nul[131]+0.07654655446197427*fl[0]*nul[131]-0.0883883476483184*fr[21]*nul[129]+0.0883883476483184*fl[21]*nul[129]+0.07654655446197427*fr[7]*nul[129]+0.07654655446197427*fl[7]*nul[129]-0.0883883476483184*fr[14]*nul[128]+0.0883883476483184*fl[14]*nul[128]+0.07654655446197427*fr[3]*nul[128]+0.07654655446197427*fl[3]*nul[128]; 
  incr2[15] = (-0.0883883476483184*fr[29]*nul[135])+0.0883883476483184*fl[29]*nul[135]+0.07654655446197427*fr[18]*nul[135]+0.07654655446197427*fl[18]*nul[135]-0.0883883476483184*fr[25]*nul[131]+0.0883883476483184*fl[25]*nul[131]+0.07654655446197427*fr[11]*nul[131]+0.07654655446197427*fl[11]*nul[131]-0.0883883476483184*fr[23]*nul[129]+0.0883883476483184*fl[23]*nul[129]+0.07654655446197427*fr[9]*nul[129]+0.07654655446197427*fl[9]*nul[129]-0.0883883476483184*fr[15]*nul[128]+0.0883883476483184*fl[15]*nul[128]+0.07654655446197427*fr[4]*nul[128]+0.07654655446197427*fl[4]*nul[128]; 
  incr2[20] = (-0.0883883476483184*fr[22]*nul[135])+0.0883883476483184*fl[22]*nul[135]+0.07654655446197427*fr[8]*nul[135]+0.07654655446197427*fl[8]*nul[135]-0.0883883476483184*fr[27]*nul[131]+0.0883883476483184*fl[27]*nul[131]+0.07654655446197427*fr[16]*nul[131]+0.07654655446197427*fl[16]*nul[131]-0.0883883476483184*fr[13]*nul[129]+0.0883883476483184*fl[13]*nul[129]+0.07654655446197427*fr[2]*nul[129]+0.07654655446197427*fl[2]*nul[129]-0.0883883476483184*fr[20]*nul[128]+0.0883883476483184*fl[20]*nul[128]+0.07654655446197427*fr[6]*nul[128]+0.07654655446197427*fl[6]*nul[128]; 
  incr2[21] = (-0.0883883476483184*fr[5]*nul[135])+0.0883883476483184*fl[5]*nul[135]+0.07654655446197427*fr[0]*nul[135]+0.07654655446197427*fl[0]*nul[135]-0.0883883476483184*fr[12]*nul[131]+0.0883883476483184*fl[12]*nul[131]+0.07654655446197427*fr[1]*nul[131]+0.07654655446197427*fl[1]*nul[131]-0.0883883476483184*fr[14]*nul[129]+0.0883883476483184*fl[14]*nul[129]+0.07654655446197427*fr[3]*nul[129]+0.07654655446197427*fl[3]*nul[129]-0.0883883476483184*fr[21]*nul[128]+0.0883883476483184*fl[21]*nul[128]+0.07654655446197427*fr[7]*nul[128]+0.07654655446197427*fl[7]*nul[128]; 
  incr2[22] = (-0.0883883476483184*fr[20]*nul[135])+0.0883883476483184*fl[20]*nul[135]+0.07654655446197427*fr[6]*nul[135]+0.07654655446197427*fl[6]*nul[135]-0.0883883476483184*fr[13]*nul[131]+0.0883883476483184*fl[13]*nul[131]+0.07654655446197427*fr[2]*nul[131]+0.07654655446197427*fl[2]*nul[131]-0.0883883476483184*fr[27]*nul[129]+0.0883883476483184*fl[27]*nul[129]+0.07654655446197427*fr[16]*nul[129]+0.07654655446197427*fl[16]*nul[129]-0.0883883476483184*fr[22]*nul[128]+0.0883883476483184*fl[22]*nul[128]+0.07654655446197427*fr[8]*nul[128]+0.07654655446197427*fl[8]*nul[128]; 
  incr2[23] = (-0.0883883476483184*fr[25]*nul[135])+0.0883883476483184*fl[25]*nul[135]+0.07654655446197427*fr[11]*nul[135]+0.07654655446197427*fl[11]*nul[135]-0.0883883476483184*fr[29]*nul[131]+0.0883883476483184*fl[29]*nul[131]+0.07654655446197427*fr[18]*nul[131]+0.07654655446197427*fl[18]*nul[131]-0.0883883476483184*fr[15]*nul[129]+0.0883883476483184*fl[15]*nul[129]+0.07654655446197427*fr[4]*nul[129]+0.07654655446197427*fl[4]*nul[129]-0.0883883476483184*fr[23]*nul[128]+0.0883883476483184*fl[23]*nul[128]+0.07654655446197427*fr[9]*nul[128]+0.07654655446197427*fl[9]*nul[128]; 
  incr2[24] = (-0.0883883476483184*fr[31]*nul[135])+0.0883883476483184*fl[31]*nul[135]+0.07654655446197427*fr[26]*nul[135]+0.07654655446197427*fl[26]*nul[135]-0.0883883476483184*fr[30]*nul[131]+0.0883883476483184*fl[30]*nul[131]+0.07654655446197427*fr[19]*nul[131]+0.07654655446197427*fl[19]*nul[131]-0.0883883476483184*fr[28]*nul[129]+0.0883883476483184*fl[28]*nul[129]+0.07654655446197427*fr[17]*nul[129]+0.07654655446197427*fl[17]*nul[129]-0.0883883476483184*fr[24]*nul[128]+0.0883883476483184*fl[24]*nul[128]+0.07654655446197427*fr[10]*nul[128]+0.07654655446197427*fl[10]*nul[128]; 
  incr2[25] = (-0.0883883476483184*fr[23]*nul[135])+0.0883883476483184*fl[23]*nul[135]+0.07654655446197427*fr[9]*nul[135]+0.07654655446197427*fl[9]*nul[135]-0.0883883476483184*fr[15]*nul[131]+0.0883883476483184*fl[15]*nul[131]+0.07654655446197427*fr[4]*nul[131]+0.07654655446197427*fl[4]*nul[131]-0.0883883476483184*fr[29]*nul[129]+0.0883883476483184*fl[29]*nul[129]+0.07654655446197427*fr[18]*nul[129]+0.07654655446197427*fl[18]*nul[129]-0.0883883476483184*fr[25]*nul[128]+0.0883883476483184*fl[25]*nul[128]+0.07654655446197427*fr[11]*nul[128]+0.07654655446197427*fl[11]*nul[128]; 
  incr2[27] = (-0.0883883476483184*fr[13]*nul[135])+0.0883883476483184*fl[13]*nul[135]+0.07654655446197427*fr[2]*nul[135]+0.07654655446197427*fl[2]*nul[135]-0.0883883476483184*fr[20]*nul[131]+0.0883883476483184*fl[20]*nul[131]+0.07654655446197427*fr[6]*nul[131]+0.07654655446197427*fl[6]*nul[131]-0.0883883476483184*fr[22]*nul[129]+0.0883883476483184*fl[22]*nul[129]+0.07654655446197427*fr[8]*nul[129]+0.07654655446197427*fl[8]*nul[129]-0.0883883476483184*fr[27]*nul[128]+0.0883883476483184*fl[27]*nul[128]+0.07654655446197427*fr[16]*nul[128]+0.07654655446197427*fl[16]*nul[128]; 
  incr2[28] = (-0.0883883476483184*fr[30]*nul[135])+0.0883883476483184*fl[30]*nul[135]+0.07654655446197427*fr[19]*nul[135]+0.07654655446197427*fl[19]*nul[135]-0.0883883476483184*fr[31]*nul[131]+0.0883883476483184*fl[31]*nul[131]+0.07654655446197427*fr[26]*nul[131]+0.07654655446197427*fl[26]*nul[131]-0.0883883476483184*fr[24]*nul[129]+0.0883883476483184*fl[24]*nul[129]+0.07654655446197427*fr[10]*nul[129]+0.07654655446197427*fl[10]*nul[129]-0.0883883476483184*fr[28]*nul[128]+0.0883883476483184*fl[28]*nul[128]+0.07654655446197427*fr[17]*nul[128]+0.07654655446197427*fl[17]*nul[128]; 
  incr2[29] = (-0.0883883476483184*fr[15]*nul[135])+0.0883883476483184*fl[15]*nul[135]+0.07654655446197427*fr[4]*nul[135]+0.07654655446197427*fl[4]*nul[135]-0.0883883476483184*fr[23]*nul[131]+0.0883883476483184*fl[23]*nul[131]+0.07654655446197427*fr[9]*nul[131]+0.07654655446197427*fl[9]*nul[131]-0.0883883476483184*fr[25]*nul[129]+0.0883883476483184*fl[25]*nul[129]+0.07654655446197427*fr[11]*nul[129]+0.07654655446197427*fl[11]*nul[129]-0.0883883476483184*fr[29]*nul[128]+0.0883883476483184*fl[29]*nul[128]+0.07654655446197427*fr[18]*nul[128]+0.07654655446197427*fl[18]*nul[128]; 
  incr2[30] = (-0.0883883476483184*fr[28]*nul[135])+0.0883883476483184*fl[28]*nul[135]+0.07654655446197427*fr[17]*nul[135]+0.07654655446197427*fl[17]*nul[135]-0.0883883476483184*fr[24]*nul[131]+0.0883883476483184*fl[24]*nul[131]+0.07654655446197427*fr[10]*nul[131]+0.07654655446197427*fl[10]*nul[131]-0.0883883476483184*fr[31]*nul[129]+0.0883883476483184*fl[31]*nul[129]+0.07654655446197427*fr[26]*nul[129]+0.07654655446197427*fl[26]*nul[129]-0.0883883476483184*fr[30]*nul[128]+0.0883883476483184*fl[30]*nul[128]+0.07654655446197427*fr[19]*nul[128]+0.07654655446197427*fl[19]*nul[128]; 
  incr2[31] = (-0.0883883476483184*fr[24]*nul[135])+0.0883883476483184*fl[24]*nul[135]+0.07654655446197427*fr[10]*nul[135]+0.07654655446197427*fl[10]*nul[135]-0.0883883476483184*fr[28]*nul[131]+0.0883883476483184*fl[28]*nul[131]+0.07654655446197427*fr[17]*nul[131]+0.07654655446197427*fl[17]*nul[131]-0.0883883476483184*fr[30]*nul[129]+0.0883883476483184*fl[30]*nul[129]+0.07654655446197427*fr[19]*nul[129]+0.07654655446197427*fl[19]*nul[129]-0.0883883476483184*fr[31]*nul[128]+0.0883883476483184*fl[31]*nul[128]+0.07654655446197427*fr[26]*nul[128]+0.07654655446197427*fl[26]*nul[128]; 

  outr[0] += incr1[0]*rdxFr; 
  outr[1] += incr1[1]*rdxFr; 
  outr[2] += incr1[2]*rdxFr; 
  outr[3] += incr1[3]*rdxFr; 
  outr[4] += incr1[4]*rdxFr; 
  outr[5] += incr2[5]*rdxFr+incr1[5]*rdxFr; 
  outr[6] += incr1[6]*rdxFr; 
  outr[7] += incr1[7]*rdxFr; 
  outr[8] += incr1[8]*rdxFr; 
  outr[9] += incr1[9]*rdxFr; 
  outr[10] += incr1[10]*rdxFr; 
  outr[11] += incr1[11]*rdxFr; 
  outr[12] += incr2[12]*rdxFr+incr1[12]*rdxFr; 
  outr[13] += incr2[13]*rdxFr+incr1[13]*rdxFr; 
  outr[14] += incr2[14]*rdxFr+incr1[14]*rdxFr; 
  outr[15] += incr2[15]*rdxFr+incr1[15]*rdxFr; 
  outr[16] += incr1[16]*rdxFr; 
  outr[17] += incr1[17]*rdxFr; 
  outr[18] += incr1[18]*rdxFr; 
  outr[19] += incr1[19]*rdxFr; 
  outr[20] += incr2[20]*rdxFr+incr1[20]*rdxFr; 
  outr[21] += incr2[21]*rdxFr+incr1[21]*rdxFr; 
  outr[22] += incr2[22]*rdxFr+incr1[22]*rdxFr; 
  outr[23] += incr2[23]*rdxFr+incr1[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr+incr1[24]*rdxFr; 
  outr[25] += incr2[25]*rdxFr+incr1[25]*rdxFr; 
  outr[26] += incr1[26]*rdxFr; 
  outr[27] += incr2[27]*rdxFr+incr1[27]*rdxFr; 
  outr[28] += incr2[28]*rdxFr+incr1[28]*rdxFr; 
  outr[29] += incr2[29]*rdxFr+incr1[29]*rdxFr; 
  outr[30] += incr2[30]*rdxFr+incr1[30]*rdxFr; 
  outr[31] += incr2[31]*rdxFr+incr1[31]*rdxFr; 

  outl[0] += -1.0*incr1[0]*rdxFl; 
  outl[1] += -1.0*incr1[1]*rdxFl; 
  outl[2] += -1.0*incr1[2]*rdxFl; 
  outl[3] += -1.0*incr1[3]*rdxFl; 
  outl[4] += -1.0*incr1[4]*rdxFl; 
  outl[5] += incr1[5]*rdxFl-1.0*incr2[5]*rdxFl; 
  outl[6] += -1.0*incr1[6]*rdxFl; 
  outl[7] += -1.0*incr1[7]*rdxFl; 
  outl[8] += -1.0*incr1[8]*rdxFl; 
  outl[9] += -1.0*incr1[9]*rdxFl; 
  outl[10] += -1.0*incr1[10]*rdxFl; 
  outl[11] += -1.0*incr1[11]*rdxFl; 
  outl[12] += incr1[12]*rdxFl-1.0*incr2[12]*rdxFl; 
  outl[13] += incr1[13]*rdxFl-1.0*incr2[13]*rdxFl; 
  outl[14] += incr1[14]*rdxFl-1.0*incr2[14]*rdxFl; 
  outl[15] += incr1[15]*rdxFl-1.0*incr2[15]*rdxFl; 
  outl[16] += -1.0*incr1[16]*rdxFl; 
  outl[17] += -1.0*incr1[17]*rdxFl; 
  outl[18] += -1.0*incr1[18]*rdxFl; 
  outl[19] += -1.0*incr1[19]*rdxFl; 
  outl[20] += incr1[20]*rdxFl-1.0*incr2[20]*rdxFl; 
  outl[21] += incr1[21]*rdxFl-1.0*incr2[21]*rdxFl; 
  outl[22] += incr1[22]*rdxFl-1.0*incr2[22]*rdxFl; 
  outl[23] += incr1[23]*rdxFl-1.0*incr2[23]*rdxFl; 
  outl[24] += incr1[24]*rdxFl-1.0*incr2[24]*rdxFl; 
  outl[25] += incr1[25]*rdxFl-1.0*incr2[25]*rdxFl; 
  outl[26] += -1.0*incr1[26]*rdxFl; 
  outl[27] += incr1[27]*rdxFl-1.0*incr2[27]*rdxFl; 
  outl[28] += incr1[28]*rdxFl-1.0*incr2[28]*rdxFl; 
  outl[29] += incr1[29]*rdxFl-1.0*incr2[29]*rdxFl; 
  outl[30] += incr1[30]*rdxFl-1.0*incr2[30]*rdxFl; 
  outl[31] += incr1[31]*rdxFl-1.0*incr2[31]*rdxFl; 

} 
