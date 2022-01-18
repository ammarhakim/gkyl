#include <ConstDiffusionModDecl.h> 
void ConstDiffusionBoundarySurf3xTensorP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[27]; 
  double incr2[27]; 

  if (edge < 0) {

  incr2[1] = 1.936491673103708*fr[7]-1.5*fr[1]+0.8660254037844386*fr[0]; 
  incr2[4] = 1.936491673103709*fr[11]-1.5*fr[4]+0.8660254037844386*fr[2]; 
  incr2[5] = 1.936491673103709*fr[13]-1.5*fr[5]+0.8660254037844386*fr[3]; 
  incr2[7] = (-7.5*fr[7])+5.809475019311125*fr[1]-3.354101966249685*fr[0]; 
  incr2[10] = 1.936491673103708*fr[17]-1.5*fr[10]+0.8660254037844386*fr[6]; 
  incr2[11] = (-7.5*fr[11])+5.809475019311126*fr[4]-3.354101966249684*fr[2]; 
  incr2[12] = 1.936491673103709*fr[20]-1.5*fr[12]+0.8660254037844387*fr[8]; 
  incr2[13] = (-7.5*fr[13])+5.809475019311126*fr[5]-3.354101966249684*fr[3]; 
  incr2[15] = 1.936491673103709*fr[21]-1.5*fr[15]+0.8660254037844387*fr[9]; 
  incr2[17] = (-7.5*fr[17])+5.809475019311125*fr[10]-3.354101966249685*fr[6]; 
  incr2[18] = 1.936491673103708*fr[23]-1.5*fr[18]+0.8660254037844387*fr[14]; 
  incr2[19] = 1.936491673103708*fr[24]-1.5*fr[19]+0.8660254037844387*fr[16]; 
  incr2[20] = (-7.5*fr[20])+5.809475019311126*fr[12]-3.354101966249685*fr[8]; 
  incr2[21] = (-7.5*fr[21])+5.809475019311126*fr[15]-3.354101966249685*fr[9]; 
  incr2[23] = (-7.5*fr[23])+5.809475019311125*fr[18]-3.354101966249684*fr[14]; 
  incr2[24] = (-7.5*fr[24])+5.809475019311125*fr[19]-3.354101966249684*fr[16]; 
  incr2[25] = 1.936491673103708*fr[26]-1.5*fr[25]+0.8660254037844386*fr[22]; 
  incr2[26] = (-7.5*fr[26])+5.809475019311125*fr[25]-3.354101966249685*fr[22]; 

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
  outr[18] += incr2[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 

  } else {

  incr2[1] = 1.936491673103708*fl[7]+1.5*fl[1]+0.8660254037844386*fl[0]; 
  incr2[4] = 1.936491673103709*fl[11]+1.5*fl[4]+0.8660254037844386*fl[2]; 
  incr2[5] = 1.936491673103709*fl[13]+1.5*fl[5]+0.8660254037844386*fl[3]; 
  incr2[7] = (-7.5*fl[7])-5.809475019311125*fl[1]-3.354101966249685*fl[0]; 
  incr2[10] = 1.936491673103708*fl[17]+1.5*fl[10]+0.8660254037844386*fl[6]; 
  incr2[11] = (-7.5*fl[11])-5.809475019311126*fl[4]-3.354101966249684*fl[2]; 
  incr2[12] = 1.936491673103709*fl[20]+1.5*fl[12]+0.8660254037844387*fl[8]; 
  incr2[13] = (-7.5*fl[13])-5.809475019311126*fl[5]-3.354101966249684*fl[3]; 
  incr2[15] = 1.936491673103709*fl[21]+1.5*fl[15]+0.8660254037844387*fl[9]; 
  incr2[17] = (-7.5*fl[17])-5.809475019311125*fl[10]-3.354101966249685*fl[6]; 
  incr2[18] = 1.936491673103708*fl[23]+1.5*fl[18]+0.8660254037844387*fl[14]; 
  incr2[19] = 1.936491673103708*fl[24]+1.5*fl[19]+0.8660254037844387*fl[16]; 
  incr2[20] = (-7.5*fl[20])-5.809475019311126*fl[12]-3.354101966249685*fl[8]; 
  incr2[21] = (-7.5*fl[21])-5.809475019311126*fl[15]-3.354101966249685*fl[9]; 
  incr2[23] = (-7.5*fl[23])-5.809475019311125*fl[18]-3.354101966249684*fl[14]; 
  incr2[24] = (-7.5*fl[24])-5.809475019311125*fl[19]-3.354101966249684*fl[16]; 
  incr2[25] = 1.936491673103708*fl[26]+1.5*fl[25]+0.8660254037844386*fl[22]; 
  incr2[26] = (-7.5*fl[26])-5.809475019311125*fl[25]-3.354101966249685*fl[22]; 

  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf3xTensorP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[27]; 
  double incr2[27]; 

  if (edge < 0) {

  incr2[2] = 1.936491673103708*fr[8]-1.5*fr[2]+0.8660254037844386*fr[0]; 
  incr2[4] = 1.936491673103709*fr[12]-1.5*fr[4]+0.8660254037844386*fr[1]; 
  incr2[6] = 1.936491673103709*fr[14]-1.5*fr[6]+0.8660254037844386*fr[3]; 
  incr2[8] = (-7.5*fr[8])+5.809475019311125*fr[2]-3.354101966249685*fr[0]; 
  incr2[10] = 1.936491673103708*fr[18]-1.5*fr[10]+0.8660254037844386*fr[5]; 
  incr2[11] = 1.936491673103709*fr[20]-1.5*fr[11]+0.8660254037844387*fr[7]; 
  incr2[12] = (-7.5*fr[12])+5.809475019311126*fr[4]-3.354101966249684*fr[1]; 
  incr2[14] = (-7.5*fr[14])+5.809475019311126*fr[6]-3.354101966249684*fr[3]; 
  incr2[16] = 1.936491673103709*fr[22]-1.5*fr[16]+0.8660254037844387*fr[9]; 
  incr2[17] = 1.936491673103708*fr[23]-1.5*fr[17]+0.8660254037844387*fr[13]; 
  incr2[18] = (-7.5*fr[18])+5.809475019311125*fr[10]-3.354101966249685*fr[5]; 
  incr2[19] = 1.936491673103708*fr[25]-1.5*fr[19]+0.8660254037844387*fr[15]; 
  incr2[20] = (-7.5*fr[20])+5.809475019311126*fr[11]-3.354101966249685*fr[7]; 
  incr2[22] = (-7.5*fr[22])+5.809475019311126*fr[16]-3.354101966249685*fr[9]; 
  incr2[23] = (-7.5*fr[23])+5.809475019311125*fr[17]-3.354101966249684*fr[13]; 
  incr2[24] = 1.936491673103708*fr[26]-1.5*fr[24]+0.8660254037844386*fr[21]; 
  incr2[25] = (-7.5*fr[25])+5.809475019311125*fr[19]-3.354101966249684*fr[15]; 
  incr2[26] = (-7.5*fr[26])+5.809475019311125*fr[24]-3.354101966249685*fr[21]; 

  outr[2] += incr2[2]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 

  } else {

  incr2[2] = 1.936491673103708*fl[8]+1.5*fl[2]+0.8660254037844386*fl[0]; 
  incr2[4] = 1.936491673103709*fl[12]+1.5*fl[4]+0.8660254037844386*fl[1]; 
  incr2[6] = 1.936491673103709*fl[14]+1.5*fl[6]+0.8660254037844386*fl[3]; 
  incr2[8] = (-7.5*fl[8])-5.809475019311125*fl[2]-3.354101966249685*fl[0]; 
  incr2[10] = 1.936491673103708*fl[18]+1.5*fl[10]+0.8660254037844386*fl[5]; 
  incr2[11] = 1.936491673103709*fl[20]+1.5*fl[11]+0.8660254037844387*fl[7]; 
  incr2[12] = (-7.5*fl[12])-5.809475019311126*fl[4]-3.354101966249684*fl[1]; 
  incr2[14] = (-7.5*fl[14])-5.809475019311126*fl[6]-3.354101966249684*fl[3]; 
  incr2[16] = 1.936491673103709*fl[22]+1.5*fl[16]+0.8660254037844387*fl[9]; 
  incr2[17] = 1.936491673103708*fl[23]+1.5*fl[17]+0.8660254037844387*fl[13]; 
  incr2[18] = (-7.5*fl[18])-5.809475019311125*fl[10]-3.354101966249685*fl[5]; 
  incr2[19] = 1.936491673103708*fl[25]+1.5*fl[19]+0.8660254037844387*fl[15]; 
  incr2[20] = (-7.5*fl[20])-5.809475019311126*fl[11]-3.354101966249685*fl[7]; 
  incr2[22] = (-7.5*fl[22])-5.809475019311126*fl[16]-3.354101966249685*fl[9]; 
  incr2[23] = (-7.5*fl[23])-5.809475019311125*fl[17]-3.354101966249684*fl[13]; 
  incr2[24] = 1.936491673103708*fl[26]+1.5*fl[24]+0.8660254037844386*fl[21]; 
  incr2[25] = (-7.5*fl[25])-5.809475019311125*fl[19]-3.354101966249684*fl[15]; 
  incr2[26] = (-7.5*fl[26])-5.809475019311125*fl[24]-3.354101966249685*fl[21]; 

  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul; 

  }

} 
void ConstDiffusionBoundarySurf3xTensorP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[27]; 
  double incr2[27]; 

  if (edge < 0) {

  incr2[3] = 1.936491673103708*fr[9]-1.5*fr[3]+0.8660254037844386*fr[0]; 
  incr2[5] = 1.936491673103709*fr[15]-1.5*fr[5]+0.8660254037844386*fr[1]; 
  incr2[6] = 1.936491673103709*fr[16]-1.5*fr[6]+0.8660254037844386*fr[2]; 
  incr2[9] = (-7.5*fr[9])+5.809475019311125*fr[3]-3.354101966249685*fr[0]; 
  incr2[10] = 1.936491673103708*fr[19]-1.5*fr[10]+0.8660254037844386*fr[4]; 
  incr2[13] = 1.936491673103709*fr[21]-1.5*fr[13]+0.8660254037844387*fr[7]; 
  incr2[14] = 1.936491673103709*fr[22]-1.5*fr[14]+0.8660254037844387*fr[8]; 
  incr2[15] = (-7.5*fr[15])+5.809475019311126*fr[5]-3.354101966249684*fr[1]; 
  incr2[16] = (-7.5*fr[16])+5.809475019311126*fr[6]-3.354101966249684*fr[2]; 
  incr2[17] = 1.936491673103708*fr[24]-1.5*fr[17]+0.8660254037844387*fr[11]; 
  incr2[18] = 1.936491673103708*fr[25]-1.5*fr[18]+0.8660254037844387*fr[12]; 
  incr2[19] = (-7.5*fr[19])+5.809475019311125*fr[10]-3.354101966249685*fr[4]; 
  incr2[21] = (-7.5*fr[21])+5.809475019311126*fr[13]-3.354101966249685*fr[7]; 
  incr2[22] = (-7.5*fr[22])+5.809475019311126*fr[14]-3.354101966249685*fr[8]; 
  incr2[23] = 1.936491673103708*fr[26]-1.5*fr[23]+0.8660254037844386*fr[20]; 
  incr2[24] = (-7.5*fr[24])+5.809475019311125*fr[17]-3.354101966249684*fr[11]; 
  incr2[25] = (-7.5*fr[25])+5.809475019311125*fr[18]-3.354101966249684*fr[12]; 
  incr2[26] = (-7.5*fr[26])+5.809475019311125*fr[23]-3.354101966249685*fr[20]; 

  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 

  } else {

  incr2[3] = 1.936491673103708*fl[9]+1.5*fl[3]+0.8660254037844386*fl[0]; 
  incr2[5] = 1.936491673103709*fl[15]+1.5*fl[5]+0.8660254037844386*fl[1]; 
  incr2[6] = 1.936491673103709*fl[16]+1.5*fl[6]+0.8660254037844386*fl[2]; 
  incr2[9] = (-7.5*fl[9])-5.809475019311125*fl[3]-3.354101966249685*fl[0]; 
  incr2[10] = 1.936491673103708*fl[19]+1.5*fl[10]+0.8660254037844386*fl[4]; 
  incr2[13] = 1.936491673103709*fl[21]+1.5*fl[13]+0.8660254037844387*fl[7]; 
  incr2[14] = 1.936491673103709*fl[22]+1.5*fl[14]+0.8660254037844387*fl[8]; 
  incr2[15] = (-7.5*fl[15])-5.809475019311126*fl[5]-3.354101966249684*fl[1]; 
  incr2[16] = (-7.5*fl[16])-5.809475019311126*fl[6]-3.354101966249684*fl[2]; 
  incr2[17] = 1.936491673103708*fl[24]+1.5*fl[17]+0.8660254037844387*fl[11]; 
  incr2[18] = 1.936491673103708*fl[25]+1.5*fl[18]+0.8660254037844387*fl[12]; 
  incr2[19] = (-7.5*fl[19])-5.809475019311125*fl[10]-3.354101966249685*fl[4]; 
  incr2[21] = (-7.5*fl[21])-5.809475019311126*fl[13]-3.354101966249685*fl[7]; 
  incr2[22] = (-7.5*fl[22])-5.809475019311126*fl[14]-3.354101966249685*fl[8]; 
  incr2[23] = 1.936491673103708*fl[26]+1.5*fl[23]+0.8660254037844386*fl[20]; 
  incr2[24] = (-7.5*fl[24])-5.809475019311125*fl[17]-3.354101966249684*fl[11]; 
  incr2[25] = (-7.5*fl[25])-5.809475019311125*fl[18]-3.354101966249684*fl[12]; 
  incr2[26] = (-7.5*fl[26])-5.809475019311125*fl[23]-3.354101966249685*fl[20]; 

  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xTensorP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[27]; 
  double incr2[27]; 
  double incr3[27]; 
  double incr4[27]; 

  if (edge < 0) {


  incr2[1] = 2.8125*fr[1]-5.083290641897234*fr[7]; 
  incr2[4] = 2.8125*fr[4]-5.083290641897235*fr[11]; 
  incr2[5] = 2.8125*fr[5]-5.083290641897235*fr[13]; 
  incr2[7] = 19.6875*fr[7]-10.89276566120836*fr[1]; 
  incr2[10] = 2.8125*fr[10]-5.083290641897234*fr[17]; 
  incr2[11] = 19.6875*fr[11]-10.89276566120836*fr[4]; 
  incr2[12] = 2.8125*fr[12]-5.083290641897235*fr[20]; 
  incr2[13] = 19.6875*fr[13]-10.89276566120836*fr[5]; 
  incr2[15] = 2.8125*fr[15]-5.083290641897235*fr[21]; 
  incr2[17] = 19.6875*fr[17]-10.89276566120836*fr[10]; 
  incr2[18] = 2.8125*fr[18]-5.083290641897234*fr[23]; 
  incr2[19] = 2.8125*fr[19]-5.083290641897234*fr[24]; 
  incr2[20] = 19.6875*fr[20]-10.89276566120836*fr[12]; 
  incr2[21] = 19.6875*fr[21]-10.89276566120836*fr[15]; 
  incr2[23] = 19.6875*fr[23]-10.89276566120836*fr[18]; 
  incr2[24] = 19.6875*fr[24]-10.89276566120836*fr[19]; 
  incr2[25] = 2.8125*fr[25]-5.083290641897234*fr[26]; 
  incr2[26] = 19.6875*fr[26]-10.89276566120836*fr[25]; 



  outr[1] += -1.0*incr2[1]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[7] += -1.0*incr2[7]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[19] += -1.0*incr2[19]*rdxFnur; 
  outr[20] += -1.0*incr2[20]*rdxFnur; 
  outr[21] += -1.0*incr2[21]*rdxFnur; 
  outr[23] += -1.0*incr2[23]*rdxFnur; 
  outr[24] += -1.0*incr2[24]*rdxFnur; 
  outr[25] += -1.0*incr2[25]*rdxFnur; 
  outr[26] += -1.0*incr2[26]*rdxFnur; 

  } else {


  incr2[1] = 5.083290641897234*fl[7]+2.8125*fl[1]; 
  incr2[4] = 5.083290641897235*fl[11]+2.8125*fl[4]; 
  incr2[5] = 5.083290641897235*fl[13]+2.8125*fl[5]; 
  incr2[7] = 19.6875*fl[7]+10.89276566120836*fl[1]; 
  incr2[10] = 5.083290641897234*fl[17]+2.8125*fl[10]; 
  incr2[11] = 19.6875*fl[11]+10.89276566120836*fl[4]; 
  incr2[12] = 5.083290641897235*fl[20]+2.8125*fl[12]; 
  incr2[13] = 19.6875*fl[13]+10.89276566120836*fl[5]; 
  incr2[15] = 5.083290641897235*fl[21]+2.8125*fl[15]; 
  incr2[17] = 19.6875*fl[17]+10.89276566120836*fl[10]; 
  incr2[18] = 5.083290641897234*fl[23]+2.8125*fl[18]; 
  incr2[19] = 5.083290641897234*fl[24]+2.8125*fl[19]; 
  incr2[20] = 19.6875*fl[20]+10.89276566120836*fl[12]; 
  incr2[21] = 19.6875*fl[21]+10.89276566120836*fl[15]; 
  incr2[23] = 19.6875*fl[23]+10.89276566120836*fl[18]; 
  incr2[24] = 19.6875*fl[24]+10.89276566120836*fl[19]; 
  incr2[25] = 5.083290641897234*fl[26]+2.8125*fl[25]; 
  incr2[26] = 19.6875*fl[26]+10.89276566120836*fl[25]; 



  outl[1] += -1.0*incr2[1]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[7] += -1.0*incr2[7]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xTensorP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[27]; 
  double incr2[27]; 
  double incr3[27]; 
  double incr4[27]; 

  if (edge < 0) {


  incr2[2] = 2.8125*fr[2]-5.083290641897234*fr[8]; 
  incr2[4] = 2.8125*fr[4]-5.083290641897235*fr[12]; 
  incr2[6] = 2.8125*fr[6]-5.083290641897235*fr[14]; 
  incr2[8] = 19.6875*fr[8]-10.89276566120836*fr[2]; 
  incr2[10] = 2.8125*fr[10]-5.083290641897234*fr[18]; 
  incr2[11] = 2.8125*fr[11]-5.083290641897235*fr[20]; 
  incr2[12] = 19.6875*fr[12]-10.89276566120836*fr[4]; 
  incr2[14] = 19.6875*fr[14]-10.89276566120836*fr[6]; 
  incr2[16] = 2.8125*fr[16]-5.083290641897235*fr[22]; 
  incr2[17] = 2.8125*fr[17]-5.083290641897234*fr[23]; 
  incr2[18] = 19.6875*fr[18]-10.89276566120836*fr[10]; 
  incr2[19] = 2.8125*fr[19]-5.083290641897234*fr[25]; 
  incr2[20] = 19.6875*fr[20]-10.89276566120836*fr[11]; 
  incr2[22] = 19.6875*fr[22]-10.89276566120836*fr[16]; 
  incr2[23] = 19.6875*fr[23]-10.89276566120836*fr[17]; 
  incr2[24] = 2.8125*fr[24]-5.083290641897234*fr[26]; 
  incr2[25] = 19.6875*fr[25]-10.89276566120836*fr[19]; 
  incr2[26] = 19.6875*fr[26]-10.89276566120836*fr[24]; 



  outr[2] += -1.0*incr2[2]*rdxFnur; 
  outr[4] += -1.0*incr2[4]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[8] += -1.0*incr2[8]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[11] += -1.0*incr2[11]*rdxFnur; 
  outr[12] += -1.0*incr2[12]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[19] += -1.0*incr2[19]*rdxFnur; 
  outr[20] += -1.0*incr2[20]*rdxFnur; 
  outr[22] += -1.0*incr2[22]*rdxFnur; 
  outr[23] += -1.0*incr2[23]*rdxFnur; 
  outr[24] += -1.0*incr2[24]*rdxFnur; 
  outr[25] += -1.0*incr2[25]*rdxFnur; 
  outr[26] += -1.0*incr2[26]*rdxFnur; 

  } else {


  incr2[2] = 5.083290641897234*fl[8]+2.8125*fl[2]; 
  incr2[4] = 5.083290641897235*fl[12]+2.8125*fl[4]; 
  incr2[6] = 5.083290641897235*fl[14]+2.8125*fl[6]; 
  incr2[8] = 19.6875*fl[8]+10.89276566120836*fl[2]; 
  incr2[10] = 5.083290641897234*fl[18]+2.8125*fl[10]; 
  incr2[11] = 5.083290641897235*fl[20]+2.8125*fl[11]; 
  incr2[12] = 19.6875*fl[12]+10.89276566120836*fl[4]; 
  incr2[14] = 19.6875*fl[14]+10.89276566120836*fl[6]; 
  incr2[16] = 5.083290641897235*fl[22]+2.8125*fl[16]; 
  incr2[17] = 5.083290641897234*fl[23]+2.8125*fl[17]; 
  incr2[18] = 19.6875*fl[18]+10.89276566120836*fl[10]; 
  incr2[19] = 5.083290641897234*fl[25]+2.8125*fl[19]; 
  incr2[20] = 19.6875*fl[20]+10.89276566120836*fl[11]; 
  incr2[22] = 19.6875*fl[22]+10.89276566120836*fl[16]; 
  incr2[23] = 19.6875*fl[23]+10.89276566120836*fl[17]; 
  incr2[24] = 5.083290641897234*fl[26]+2.8125*fl[24]; 
  incr2[25] = 19.6875*fl[25]+10.89276566120836*fl[19]; 
  incr2[26] = 19.6875*fl[26]+10.89276566120836*fl[24]; 



  outl[2] += -1.0*incr2[2]*rdxFnul; 
  outl[4] += -1.0*incr2[4]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[8] += -1.0*incr2[8]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[11] += -1.0*incr2[11]*rdxFnul; 
  outl[12] += -1.0*incr2[12]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[20] += -1.0*incr2[20]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 

  }

} 
void ConstHyperDiffusion4BoundarySurf3xTensorP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[27]; 
  double incr2[27]; 
  double incr3[27]; 
  double incr4[27]; 

  if (edge < 0) {


  incr2[3] = 2.8125*fr[3]-5.083290641897234*fr[9]; 
  incr2[5] = 2.8125*fr[5]-5.083290641897235*fr[15]; 
  incr2[6] = 2.8125*fr[6]-5.083290641897235*fr[16]; 
  incr2[9] = 19.6875*fr[9]-10.89276566120836*fr[3]; 
  incr2[10] = 2.8125*fr[10]-5.083290641897234*fr[19]; 
  incr2[13] = 2.8125*fr[13]-5.083290641897235*fr[21]; 
  incr2[14] = 2.8125*fr[14]-5.083290641897235*fr[22]; 
  incr2[15] = 19.6875*fr[15]-10.89276566120836*fr[5]; 
  incr2[16] = 19.6875*fr[16]-10.89276566120836*fr[6]; 
  incr2[17] = 2.8125*fr[17]-5.083290641897234*fr[24]; 
  incr2[18] = 2.8125*fr[18]-5.083290641897234*fr[25]; 
  incr2[19] = 19.6875*fr[19]-10.89276566120836*fr[10]; 
  incr2[21] = 19.6875*fr[21]-10.89276566120836*fr[13]; 
  incr2[22] = 19.6875*fr[22]-10.89276566120836*fr[14]; 
  incr2[23] = 2.8125*fr[23]-5.083290641897234*fr[26]; 
  incr2[24] = 19.6875*fr[24]-10.89276566120836*fr[17]; 
  incr2[25] = 19.6875*fr[25]-10.89276566120836*fr[18]; 
  incr2[26] = 19.6875*fr[26]-10.89276566120836*fr[23]; 



  outr[3] += -1.0*incr2[3]*rdxFnur; 
  outr[5] += -1.0*incr2[5]*rdxFnur; 
  outr[6] += -1.0*incr2[6]*rdxFnur; 
  outr[9] += -1.0*incr2[9]*rdxFnur; 
  outr[10] += -1.0*incr2[10]*rdxFnur; 
  outr[13] += -1.0*incr2[13]*rdxFnur; 
  outr[14] += -1.0*incr2[14]*rdxFnur; 
  outr[15] += -1.0*incr2[15]*rdxFnur; 
  outr[16] += -1.0*incr2[16]*rdxFnur; 
  outr[17] += -1.0*incr2[17]*rdxFnur; 
  outr[18] += -1.0*incr2[18]*rdxFnur; 
  outr[19] += -1.0*incr2[19]*rdxFnur; 
  outr[21] += -1.0*incr2[21]*rdxFnur; 
  outr[22] += -1.0*incr2[22]*rdxFnur; 
  outr[23] += -1.0*incr2[23]*rdxFnur; 
  outr[24] += -1.0*incr2[24]*rdxFnur; 
  outr[25] += -1.0*incr2[25]*rdxFnur; 
  outr[26] += -1.0*incr2[26]*rdxFnur; 

  } else {


  incr2[3] = 5.083290641897234*fl[9]+2.8125*fl[3]; 
  incr2[5] = 5.083290641897235*fl[15]+2.8125*fl[5]; 
  incr2[6] = 5.083290641897235*fl[16]+2.8125*fl[6]; 
  incr2[9] = 19.6875*fl[9]+10.89276566120836*fl[3]; 
  incr2[10] = 5.083290641897234*fl[19]+2.8125*fl[10]; 
  incr2[13] = 5.083290641897235*fl[21]+2.8125*fl[13]; 
  incr2[14] = 5.083290641897235*fl[22]+2.8125*fl[14]; 
  incr2[15] = 19.6875*fl[15]+10.89276566120836*fl[5]; 
  incr2[16] = 19.6875*fl[16]+10.89276566120836*fl[6]; 
  incr2[17] = 5.083290641897234*fl[24]+2.8125*fl[17]; 
  incr2[18] = 5.083290641897234*fl[25]+2.8125*fl[18]; 
  incr2[19] = 19.6875*fl[19]+10.89276566120836*fl[10]; 
  incr2[21] = 19.6875*fl[21]+10.89276566120836*fl[13]; 
  incr2[22] = 19.6875*fl[22]+10.89276566120836*fl[14]; 
  incr2[23] = 5.083290641897234*fl[26]+2.8125*fl[23]; 
  incr2[24] = 19.6875*fl[24]+10.89276566120836*fl[17]; 
  incr2[25] = 19.6875*fl[25]+10.89276566120836*fl[18]; 
  incr2[26] = 19.6875*fl[26]+10.89276566120836*fl[23]; 



  outl[3] += -1.0*incr2[3]*rdxFnul; 
  outl[5] += -1.0*incr2[5]*rdxFnul; 
  outl[6] += -1.0*incr2[6]*rdxFnul; 
  outl[9] += -1.0*incr2[9]*rdxFnul; 
  outl[10] += -1.0*incr2[10]*rdxFnul; 
  outl[13] += -1.0*incr2[13]*rdxFnul; 
  outl[14] += -1.0*incr2[14]*rdxFnul; 
  outl[15] += -1.0*incr2[15]*rdxFnul; 
  outl[16] += -1.0*incr2[16]*rdxFnul; 
  outl[17] += -1.0*incr2[17]*rdxFnul; 
  outl[18] += -1.0*incr2[18]*rdxFnul; 
  outl[19] += -1.0*incr2[19]*rdxFnul; 
  outl[21] += -1.0*incr2[21]*rdxFnul; 
  outl[22] += -1.0*incr2[22]*rdxFnul; 
  outl[23] += -1.0*incr2[23]*rdxFnul; 
  outl[24] += -1.0*incr2[24]*rdxFnul; 
  outl[25] += -1.0*incr2[25]*rdxFnul; 
  outl[26] += -1.0*incr2[26]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf3xTensorP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[27]; 
  double incr2[27]; 
  double incr3[27]; 
  double incr4[27]; 
  double incr5[27]; 
  double incr6[27]; 

  if (edge < 0) {


  incr2[1] = 19.06233990711463*fr[7]-4.921875*fr[1]; 
  incr2[4] = 19.06233990711463*fr[11]-4.921875*fr[4]; 
  incr2[5] = 19.06233990711463*fr[13]-4.921875*fr[5]; 
  incr2[7] = 19.06233990711463*fr[1]-73.828125*fr[7]; 
  incr2[10] = 19.06233990711463*fr[17]-4.921875*fr[10]; 
  incr2[11] = 19.06233990711463*fr[4]-73.828125*fr[11]; 
  incr2[12] = 19.06233990711463*fr[20]-4.921875*fr[12]; 
  incr2[13] = 19.06233990711463*fr[5]-73.828125*fr[13]; 
  incr2[15] = 19.06233990711463*fr[21]-4.921875*fr[15]; 
  incr2[17] = 19.06233990711463*fr[10]-73.828125*fr[17]; 
  incr2[18] = 19.06233990711463*fr[23]-4.921875*fr[18]; 
  incr2[19] = 19.06233990711463*fr[24]-4.921875*fr[19]; 
  incr2[20] = 19.06233990711463*fr[12]-73.828125*fr[20]; 
  incr2[21] = 19.06233990711463*fr[15]-73.828125*fr[21]; 
  incr2[23] = 19.06233990711463*fr[18]-73.828125*fr[23]; 
  incr2[24] = 19.06233990711463*fr[19]-73.828125*fr[24]; 
  incr2[25] = 19.06233990711463*fr[26]-4.921875*fr[25]; 
  incr2[26] = 19.06233990711463*fr[25]-73.828125*fr[26]; 





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
  outr[18] += incr2[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 

  } else {


  incr2[1] = (-19.06233990711463*fl[7])-4.921875*fl[1]; 
  incr2[4] = (-19.06233990711463*fl[11])-4.921875*fl[4]; 
  incr2[5] = (-19.06233990711463*fl[13])-4.921875*fl[5]; 
  incr2[7] = (-73.828125*fl[7])-19.06233990711463*fl[1]; 
  incr2[10] = (-19.06233990711463*fl[17])-4.921875*fl[10]; 
  incr2[11] = (-73.828125*fl[11])-19.06233990711463*fl[4]; 
  incr2[12] = (-19.06233990711463*fl[20])-4.921875*fl[12]; 
  incr2[13] = (-73.828125*fl[13])-19.06233990711463*fl[5]; 
  incr2[15] = (-19.06233990711463*fl[21])-4.921875*fl[15]; 
  incr2[17] = (-73.828125*fl[17])-19.06233990711463*fl[10]; 
  incr2[18] = (-19.06233990711463*fl[23])-4.921875*fl[18]; 
  incr2[19] = (-19.06233990711463*fl[24])-4.921875*fl[19]; 
  incr2[20] = (-73.828125*fl[20])-19.06233990711463*fl[12]; 
  incr2[21] = (-73.828125*fl[21])-19.06233990711463*fl[15]; 
  incr2[23] = (-73.828125*fl[23])-19.06233990711463*fl[18]; 
  incr2[24] = (-73.828125*fl[24])-19.06233990711463*fl[19]; 
  incr2[25] = (-19.06233990711463*fl[26])-4.921875*fl[25]; 
  incr2[26] = (-73.828125*fl[26])-19.06233990711463*fl[25]; 





  outl[1] += incr2[1]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[7] += incr2[7]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf3xTensorP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[27]; 
  double incr2[27]; 
  double incr3[27]; 
  double incr4[27]; 
  double incr5[27]; 
  double incr6[27]; 

  if (edge < 0) {


  incr2[2] = 19.06233990711463*fr[8]-4.921875*fr[2]; 
  incr2[4] = 19.06233990711463*fr[12]-4.921875*fr[4]; 
  incr2[6] = 19.06233990711463*fr[14]-4.921875*fr[6]; 
  incr2[8] = 19.06233990711463*fr[2]-73.828125*fr[8]; 
  incr2[10] = 19.06233990711463*fr[18]-4.921875*fr[10]; 
  incr2[11] = 19.06233990711463*fr[20]-4.921875*fr[11]; 
  incr2[12] = 19.06233990711463*fr[4]-73.828125*fr[12]; 
  incr2[14] = 19.06233990711463*fr[6]-73.828125*fr[14]; 
  incr2[16] = 19.06233990711463*fr[22]-4.921875*fr[16]; 
  incr2[17] = 19.06233990711463*fr[23]-4.921875*fr[17]; 
  incr2[18] = 19.06233990711463*fr[10]-73.828125*fr[18]; 
  incr2[19] = 19.06233990711463*fr[25]-4.921875*fr[19]; 
  incr2[20] = 19.06233990711463*fr[11]-73.828125*fr[20]; 
  incr2[22] = 19.06233990711463*fr[16]-73.828125*fr[22]; 
  incr2[23] = 19.06233990711463*fr[17]-73.828125*fr[23]; 
  incr2[24] = 19.06233990711463*fr[26]-4.921875*fr[24]; 
  incr2[25] = 19.06233990711463*fr[19]-73.828125*fr[25]; 
  incr2[26] = 19.06233990711463*fr[24]-73.828125*fr[26]; 





  outr[2] += incr2[2]*rdxFnur; 
  outr[4] += incr2[4]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[8] += incr2[8]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[11] += incr2[11]*rdxFnur; 
  outr[12] += incr2[12]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[20] += incr2[20]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 

  } else {


  incr2[2] = (-19.06233990711463*fl[8])-4.921875*fl[2]; 
  incr2[4] = (-19.06233990711463*fl[12])-4.921875*fl[4]; 
  incr2[6] = (-19.06233990711463*fl[14])-4.921875*fl[6]; 
  incr2[8] = (-73.828125*fl[8])-19.06233990711463*fl[2]; 
  incr2[10] = (-19.06233990711463*fl[18])-4.921875*fl[10]; 
  incr2[11] = (-19.06233990711463*fl[20])-4.921875*fl[11]; 
  incr2[12] = (-73.828125*fl[12])-19.06233990711463*fl[4]; 
  incr2[14] = (-73.828125*fl[14])-19.06233990711463*fl[6]; 
  incr2[16] = (-19.06233990711463*fl[22])-4.921875*fl[16]; 
  incr2[17] = (-19.06233990711463*fl[23])-4.921875*fl[17]; 
  incr2[18] = (-73.828125*fl[18])-19.06233990711463*fl[10]; 
  incr2[19] = (-19.06233990711463*fl[25])-4.921875*fl[19]; 
  incr2[20] = (-73.828125*fl[20])-19.06233990711463*fl[11]; 
  incr2[22] = (-73.828125*fl[22])-19.06233990711463*fl[16]; 
  incr2[23] = (-73.828125*fl[23])-19.06233990711463*fl[17]; 
  incr2[24] = (-19.06233990711463*fl[26])-4.921875*fl[24]; 
  incr2[25] = (-73.828125*fl[25])-19.06233990711463*fl[19]; 
  incr2[26] = (-73.828125*fl[26])-19.06233990711463*fl[24]; 





  outl[2] += incr2[2]*rdxFnul; 
  outl[4] += incr2[4]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[8] += incr2[8]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[11] += incr2[11]*rdxFnul; 
  outl[12] += incr2[12]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul; 
  outl[20] += incr2[20]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul; 

  }

} 
void ConstHyperDiffusion6BoundarySurf3xTensorP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
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

  double incr1[27]; 
  double incr2[27]; 
  double incr3[27]; 
  double incr4[27]; 
  double incr5[27]; 
  double incr6[27]; 

  if (edge < 0) {


  incr2[3] = 19.06233990711463*fr[9]-4.921875*fr[3]; 
  incr2[5] = 19.06233990711463*fr[15]-4.921875*fr[5]; 
  incr2[6] = 19.06233990711463*fr[16]-4.921875*fr[6]; 
  incr2[9] = 19.06233990711463*fr[3]-73.828125*fr[9]; 
  incr2[10] = 19.06233990711463*fr[19]-4.921875*fr[10]; 
  incr2[13] = 19.06233990711463*fr[21]-4.921875*fr[13]; 
  incr2[14] = 19.06233990711463*fr[22]-4.921875*fr[14]; 
  incr2[15] = 19.06233990711463*fr[5]-73.828125*fr[15]; 
  incr2[16] = 19.06233990711463*fr[6]-73.828125*fr[16]; 
  incr2[17] = 19.06233990711463*fr[24]-4.921875*fr[17]; 
  incr2[18] = 19.06233990711463*fr[25]-4.921875*fr[18]; 
  incr2[19] = 19.06233990711463*fr[10]-73.828125*fr[19]; 
  incr2[21] = 19.06233990711463*fr[13]-73.828125*fr[21]; 
  incr2[22] = 19.06233990711463*fr[14]-73.828125*fr[22]; 
  incr2[23] = 19.06233990711463*fr[26]-4.921875*fr[23]; 
  incr2[24] = 19.06233990711463*fr[17]-73.828125*fr[24]; 
  incr2[25] = 19.06233990711463*fr[18]-73.828125*fr[25]; 
  incr2[26] = 19.06233990711463*fr[23]-73.828125*fr[26]; 





  outr[3] += incr2[3]*rdxFnur; 
  outr[5] += incr2[5]*rdxFnur; 
  outr[6] += incr2[6]*rdxFnur; 
  outr[9] += incr2[9]*rdxFnur; 
  outr[10] += incr2[10]*rdxFnur; 
  outr[13] += incr2[13]*rdxFnur; 
  outr[14] += incr2[14]*rdxFnur; 
  outr[15] += incr2[15]*rdxFnur; 
  outr[16] += incr2[16]*rdxFnur; 
  outr[17] += incr2[17]*rdxFnur; 
  outr[18] += incr2[18]*rdxFnur; 
  outr[19] += incr2[19]*rdxFnur; 
  outr[21] += incr2[21]*rdxFnur; 
  outr[22] += incr2[22]*rdxFnur; 
  outr[23] += incr2[23]*rdxFnur; 
  outr[24] += incr2[24]*rdxFnur; 
  outr[25] += incr2[25]*rdxFnur; 
  outr[26] += incr2[26]*rdxFnur; 

  } else {


  incr2[3] = (-19.06233990711463*fl[9])-4.921875*fl[3]; 
  incr2[5] = (-19.06233990711463*fl[15])-4.921875*fl[5]; 
  incr2[6] = (-19.06233990711463*fl[16])-4.921875*fl[6]; 
  incr2[9] = (-73.828125*fl[9])-19.06233990711463*fl[3]; 
  incr2[10] = (-19.06233990711463*fl[19])-4.921875*fl[10]; 
  incr2[13] = (-19.06233990711463*fl[21])-4.921875*fl[13]; 
  incr2[14] = (-19.06233990711463*fl[22])-4.921875*fl[14]; 
  incr2[15] = (-73.828125*fl[15])-19.06233990711463*fl[5]; 
  incr2[16] = (-73.828125*fl[16])-19.06233990711463*fl[6]; 
  incr2[17] = (-19.06233990711463*fl[24])-4.921875*fl[17]; 
  incr2[18] = (-19.06233990711463*fl[25])-4.921875*fl[18]; 
  incr2[19] = (-73.828125*fl[19])-19.06233990711463*fl[10]; 
  incr2[21] = (-73.828125*fl[21])-19.06233990711463*fl[13]; 
  incr2[22] = (-73.828125*fl[22])-19.06233990711463*fl[14]; 
  incr2[23] = (-19.06233990711463*fl[26])-4.921875*fl[23]; 
  incr2[24] = (-73.828125*fl[24])-19.06233990711463*fl[17]; 
  incr2[25] = (-73.828125*fl[25])-19.06233990711463*fl[18]; 
  incr2[26] = (-73.828125*fl[26])-19.06233990711463*fl[23]; 





  outl[3] += incr2[3]*rdxFnul; 
  outl[5] += incr2[5]*rdxFnul; 
  outl[6] += incr2[6]*rdxFnul; 
  outl[9] += incr2[9]*rdxFnul; 
  outl[10] += incr2[10]*rdxFnul; 
  outl[13] += incr2[13]*rdxFnul; 
  outl[14] += incr2[14]*rdxFnul; 
  outl[15] += incr2[15]*rdxFnul; 
  outl[16] += incr2[16]*rdxFnul; 
  outl[17] += incr2[17]*rdxFnul; 
  outl[18] += incr2[18]*rdxFnul; 
  outl[19] += incr2[19]*rdxFnul; 
  outl[21] += incr2[21]*rdxFnul; 
  outl[22] += incr2[22]*rdxFnul; 
  outl[23] += incr2[23]*rdxFnul; 
  outl[24] += incr2[24]*rdxFnul; 
  outl[25] += incr2[25]*rdxFnul; 
  outl[26] += incr2[26]*rdxFnul; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf3xTensorP2_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[81]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[0]*dxl[0]); 
  double rdxFr = 4.0/(dxr[0]*dxr[0]); 

  double incr1[27]; 
  double incr2[27]; 

  if (edge < 0) {

  incr2[1] = 1.530931089239485*fr[7]*nul[7]-1.185854122563142*fr[1]*nul[7]+0.6846531968814573*fr[0]*nul[7]+1.185854122563142*nul[1]*fr[7]+0.6846531968814573*nul[0]*fr[7]-0.9185586535436913*fr[1]*nul[1]+0.5303300858899105*fr[0]*nul[1]-0.5303300858899105*nul[0]*fr[1]+0.3061862178478971*fr[0]*nul[0]; 
  incr2[4] = 1.530931089239486*nul[7]*fr[11]+1.185854122563142*nul[1]*fr[11]+0.6846531968814574*nul[0]*fr[11]-1.185854122563142*fr[4]*nul[7]+0.6846531968814573*fr[2]*nul[7]-0.9185586535436913*nul[1]*fr[4]-0.5303300858899105*nul[0]*fr[4]+0.5303300858899105*nul[1]*fr[2]+0.3061862178478971*nul[0]*fr[2]; 
  incr2[5] = 1.530931089239486*nul[7]*fr[13]+1.185854122563142*nul[1]*fr[13]+0.6846531968814574*nul[0]*fr[13]-1.185854122563142*fr[5]*nul[7]+0.6846531968814573*fr[3]*nul[7]-0.9185586535436913*nul[1]*fr[5]-0.5303300858899105*nul[0]*fr[5]+0.5303300858899105*nul[1]*fr[3]+0.3061862178478971*nul[0]*fr[3]; 
  incr2[7] = (-5.929270612815711*fr[7]*nul[7])+4.592793267718456*fr[1]*nul[7]-2.651650429449552*fr[0]*nul[7]-4.592793267718455*nul[1]*fr[7]-2.651650429449552*nul[0]*fr[7]+3.557562367689425*fr[1]*nul[1]-2.053959590644372*fr[0]*nul[1]+2.053959590644372*nul[0]*fr[1]-1.185854122563142*fr[0]*nul[0]; 
  incr2[10] = 1.530931089239485*nul[7]*fr[17]+1.185854122563142*nul[1]*fr[17]+0.6846531968814573*nul[0]*fr[17]-1.185854122563142*nul[7]*fr[10]-0.9185586535436913*nul[1]*fr[10]-0.5303300858899105*nul[0]*fr[10]+0.6846531968814573*fr[6]*nul[7]+0.5303300858899105*nul[1]*fr[6]+0.3061862178478971*nul[0]*fr[6]; 
  incr2[11] = (-5.929270612815711*nul[7]*fr[11])-4.592793267718455*nul[1]*fr[11]-2.651650429449552*nul[0]*fr[11]+4.592793267718458*fr[4]*nul[7]-2.651650429449552*fr[2]*nul[7]+3.557562367689425*nul[1]*fr[4]+2.053959590644372*nul[0]*fr[4]-2.053959590644372*nul[1]*fr[2]-1.185854122563142*nul[0]*fr[2]; 
  incr2[12] = 1.530931089239486*nul[7]*fr[20]+1.185854122563142*nul[1]*fr[20]+0.6846531968814575*nul[0]*fr[20]-1.185854122563142*nul[7]*fr[12]-0.9185586535436913*nul[1]*fr[12]-0.5303300858899105*nul[0]*fr[12]+0.6846531968814574*nul[7]*fr[8]+0.5303300858899104*nul[1]*fr[8]+0.3061862178478971*nul[0]*fr[8]; 
  incr2[13] = (-5.929270612815711*nul[7]*fr[13])-4.592793267718455*nul[1]*fr[13]-2.651650429449552*nul[0]*fr[13]+4.592793267718458*fr[5]*nul[7]-2.651650429449552*fr[3]*nul[7]+3.557562367689425*nul[1]*fr[5]+2.053959590644372*nul[0]*fr[5]-2.053959590644372*nul[1]*fr[3]-1.185854122563142*nul[0]*fr[3]; 
  incr2[15] = 1.530931089239486*nul[7]*fr[21]+1.185854122563142*nul[1]*fr[21]+0.6846531968814575*nul[0]*fr[21]-1.185854122563142*nul[7]*fr[15]-0.9185586535436913*nul[1]*fr[15]-0.5303300858899105*nul[0]*fr[15]+0.6846531968814574*nul[7]*fr[9]+0.5303300858899104*nul[1]*fr[9]+0.3061862178478971*nul[0]*fr[9]; 
  incr2[17] = (-5.929270612815711*nul[7]*fr[17])-4.592793267718455*nul[1]*fr[17]-2.651650429449552*nul[0]*fr[17]+4.592793267718456*nul[7]*fr[10]+3.557562367689425*nul[1]*fr[10]+2.053959590644372*nul[0]*fr[10]-2.651650429449552*fr[6]*nul[7]-2.053959590644372*nul[1]*fr[6]-1.185854122563142*nul[0]*fr[6]; 
  incr2[18] = 1.530931089239485*nul[7]*fr[23]+1.185854122563142*nul[1]*fr[23]+0.6846531968814573*nul[0]*fr[23]-1.185854122563142*nul[7]*fr[18]-0.9185586535436913*nul[1]*fr[18]-0.5303300858899105*nul[0]*fr[18]+0.6846531968814574*nul[7]*fr[14]+0.5303300858899104*nul[1]*fr[14]+0.3061862178478971*nul[0]*fr[14]; 
  incr2[19] = 1.530931089239485*nul[7]*fr[24]+1.185854122563142*nul[1]*fr[24]+0.6846531968814573*nul[0]*fr[24]-1.185854122563142*nul[7]*fr[19]-0.9185586535436913*nul[1]*fr[19]-0.5303300858899105*nul[0]*fr[19]+0.6846531968814574*nul[7]*fr[16]+0.5303300858899104*nul[1]*fr[16]+0.3061862178478971*nul[0]*fr[16]; 
  incr2[20] = (-5.929270612815711*nul[7]*fr[20])-4.592793267718455*nul[1]*fr[20]-2.651650429449552*nul[0]*fr[20]+4.592793267718458*nul[7]*fr[12]+3.557562367689425*nul[1]*fr[12]+2.053959590644372*nul[0]*fr[12]-2.651650429449552*nul[7]*fr[8]-2.053959590644372*nul[1]*fr[8]-1.185854122563142*nul[0]*fr[8]; 
  incr2[21] = (-5.929270612815711*nul[7]*fr[21])-4.592793267718455*nul[1]*fr[21]-2.651650429449552*nul[0]*fr[21]+4.592793267718458*nul[7]*fr[15]+3.557562367689425*nul[1]*fr[15]+2.053959590644372*nul[0]*fr[15]-2.651650429449552*nul[7]*fr[9]-2.053959590644372*nul[1]*fr[9]-1.185854122563142*nul[0]*fr[9]; 
  incr2[23] = (-5.929270612815711*nul[7]*fr[23])-4.592793267718455*nul[1]*fr[23]-2.651650429449552*nul[0]*fr[23]+4.592793267718456*nul[7]*fr[18]+3.557562367689425*nul[1]*fr[18]+2.053959590644372*nul[0]*fr[18]-2.651650429449552*nul[7]*fr[14]-2.053959590644372*nul[1]*fr[14]-1.185854122563142*nul[0]*fr[14]; 
  incr2[24] = (-5.929270612815711*nul[7]*fr[24])-4.592793267718455*nul[1]*fr[24]-2.651650429449552*nul[0]*fr[24]+4.592793267718456*nul[7]*fr[19]+3.557562367689425*nul[1]*fr[19]+2.053959590644372*nul[0]*fr[19]-2.651650429449552*nul[7]*fr[16]-2.053959590644372*nul[1]*fr[16]-1.185854122563142*nul[0]*fr[16]; 
  incr2[25] = 1.530931089239485*nul[7]*fr[26]+1.185854122563142*nul[1]*fr[26]+0.6846531968814573*nul[0]*fr[26]-1.185854122563142*nul[7]*fr[25]-0.9185586535436913*nul[1]*fr[25]-0.5303300858899105*nul[0]*fr[25]+0.6846531968814573*nul[7]*fr[22]+0.5303300858899105*nul[1]*fr[22]+0.3061862178478971*nul[0]*fr[22]; 
  incr2[26] = (-5.929270612815711*nul[7]*fr[26])-4.592793267718455*nul[1]*fr[26]-2.651650429449552*nul[0]*fr[26]+4.592793267718456*nul[7]*fr[25]+3.557562367689425*nul[1]*fr[25]+2.053959590644372*nul[0]*fr[25]-2.651650429449552*nul[7]*fr[22]-2.053959590644372*nul[1]*fr[22]-1.185854122563142*nul[0]*fr[22]; 

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
  outr[18] += incr2[18]*rdxFr; 
  outr[19] += incr2[19]*rdxFr; 
  outr[20] += incr2[20]*rdxFr; 
  outr[21] += incr2[21]*rdxFr; 
  outr[23] += incr2[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr; 
  outr[25] += incr2[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 

  } else {

  incr2[1] = 1.530931089239485*fl[7]*nul[7]+1.185854122563142*fl[1]*nul[7]+0.6846531968814573*fl[0]*nul[7]+1.185854122563142*nul[1]*fl[7]+0.6846531968814573*nul[0]*fl[7]+0.9185586535436913*fl[1]*nul[1]+0.5303300858899105*fl[0]*nul[1]+0.5303300858899105*nul[0]*fl[1]+0.3061862178478971*fl[0]*nul[0]; 
  incr2[4] = 1.530931089239486*nul[7]*fl[11]+1.185854122563142*nul[1]*fl[11]+0.6846531968814574*nul[0]*fl[11]+1.185854122563142*fl[4]*nul[7]+0.6846531968814573*fl[2]*nul[7]+0.9185586535436913*nul[1]*fl[4]+0.5303300858899105*nul[0]*fl[4]+0.5303300858899105*nul[1]*fl[2]+0.3061862178478971*nul[0]*fl[2]; 
  incr2[5] = 1.530931089239486*nul[7]*fl[13]+1.185854122563142*nul[1]*fl[13]+0.6846531968814574*nul[0]*fl[13]+1.185854122563142*fl[5]*nul[7]+0.6846531968814573*fl[3]*nul[7]+0.9185586535436913*nul[1]*fl[5]+0.5303300858899105*nul[0]*fl[5]+0.5303300858899105*nul[1]*fl[3]+0.3061862178478971*nul[0]*fl[3]; 
  incr2[7] = (-5.929270612815711*fl[7]*nul[7])-4.592793267718456*fl[1]*nul[7]-2.651650429449552*fl[0]*nul[7]-4.592793267718455*nul[1]*fl[7]-2.651650429449552*nul[0]*fl[7]-3.557562367689425*fl[1]*nul[1]-2.053959590644372*fl[0]*nul[1]-2.053959590644372*nul[0]*fl[1]-1.185854122563142*fl[0]*nul[0]; 
  incr2[10] = 1.530931089239485*nul[7]*fl[17]+1.185854122563142*nul[1]*fl[17]+0.6846531968814573*nul[0]*fl[17]+1.185854122563142*nul[7]*fl[10]+0.9185586535436913*nul[1]*fl[10]+0.5303300858899105*nul[0]*fl[10]+0.6846531968814573*fl[6]*nul[7]+0.5303300858899105*nul[1]*fl[6]+0.3061862178478971*nul[0]*fl[6]; 
  incr2[11] = (-5.929270612815711*nul[7]*fl[11])-4.592793267718455*nul[1]*fl[11]-2.651650429449552*nul[0]*fl[11]-4.592793267718458*fl[4]*nul[7]-2.651650429449552*fl[2]*nul[7]-3.557562367689425*nul[1]*fl[4]-2.053959590644372*nul[0]*fl[4]-2.053959590644372*nul[1]*fl[2]-1.185854122563142*nul[0]*fl[2]; 
  incr2[12] = 1.530931089239486*nul[7]*fl[20]+1.185854122563142*nul[1]*fl[20]+0.6846531968814575*nul[0]*fl[20]+1.185854122563142*nul[7]*fl[12]+0.9185586535436913*nul[1]*fl[12]+0.5303300858899105*nul[0]*fl[12]+0.6846531968814574*nul[7]*fl[8]+0.5303300858899104*nul[1]*fl[8]+0.3061862178478971*nul[0]*fl[8]; 
  incr2[13] = (-5.929270612815711*nul[7]*fl[13])-4.592793267718455*nul[1]*fl[13]-2.651650429449552*nul[0]*fl[13]-4.592793267718458*fl[5]*nul[7]-2.651650429449552*fl[3]*nul[7]-3.557562367689425*nul[1]*fl[5]-2.053959590644372*nul[0]*fl[5]-2.053959590644372*nul[1]*fl[3]-1.185854122563142*nul[0]*fl[3]; 
  incr2[15] = 1.530931089239486*nul[7]*fl[21]+1.185854122563142*nul[1]*fl[21]+0.6846531968814575*nul[0]*fl[21]+1.185854122563142*nul[7]*fl[15]+0.9185586535436913*nul[1]*fl[15]+0.5303300858899105*nul[0]*fl[15]+0.6846531968814574*nul[7]*fl[9]+0.5303300858899104*nul[1]*fl[9]+0.3061862178478971*nul[0]*fl[9]; 
  incr2[17] = (-5.929270612815711*nul[7]*fl[17])-4.592793267718455*nul[1]*fl[17]-2.651650429449552*nul[0]*fl[17]-4.592793267718456*nul[7]*fl[10]-3.557562367689425*nul[1]*fl[10]-2.053959590644372*nul[0]*fl[10]-2.651650429449552*fl[6]*nul[7]-2.053959590644372*nul[1]*fl[6]-1.185854122563142*nul[0]*fl[6]; 
  incr2[18] = 1.530931089239485*nul[7]*fl[23]+1.185854122563142*nul[1]*fl[23]+0.6846531968814573*nul[0]*fl[23]+1.185854122563142*nul[7]*fl[18]+0.9185586535436913*nul[1]*fl[18]+0.5303300858899105*nul[0]*fl[18]+0.6846531968814574*nul[7]*fl[14]+0.5303300858899104*nul[1]*fl[14]+0.3061862178478971*nul[0]*fl[14]; 
  incr2[19] = 1.530931089239485*nul[7]*fl[24]+1.185854122563142*nul[1]*fl[24]+0.6846531968814573*nul[0]*fl[24]+1.185854122563142*nul[7]*fl[19]+0.9185586535436913*nul[1]*fl[19]+0.5303300858899105*nul[0]*fl[19]+0.6846531968814574*nul[7]*fl[16]+0.5303300858899104*nul[1]*fl[16]+0.3061862178478971*nul[0]*fl[16]; 
  incr2[20] = (-5.929270612815711*nul[7]*fl[20])-4.592793267718455*nul[1]*fl[20]-2.651650429449552*nul[0]*fl[20]-4.592793267718458*nul[7]*fl[12]-3.557562367689425*nul[1]*fl[12]-2.053959590644372*nul[0]*fl[12]-2.651650429449552*nul[7]*fl[8]-2.053959590644372*nul[1]*fl[8]-1.185854122563142*nul[0]*fl[8]; 
  incr2[21] = (-5.929270612815711*nul[7]*fl[21])-4.592793267718455*nul[1]*fl[21]-2.651650429449552*nul[0]*fl[21]-4.592793267718458*nul[7]*fl[15]-3.557562367689425*nul[1]*fl[15]-2.053959590644372*nul[0]*fl[15]-2.651650429449552*nul[7]*fl[9]-2.053959590644372*nul[1]*fl[9]-1.185854122563142*nul[0]*fl[9]; 
  incr2[23] = (-5.929270612815711*nul[7]*fl[23])-4.592793267718455*nul[1]*fl[23]-2.651650429449552*nul[0]*fl[23]-4.592793267718456*nul[7]*fl[18]-3.557562367689425*nul[1]*fl[18]-2.053959590644372*nul[0]*fl[18]-2.651650429449552*nul[7]*fl[14]-2.053959590644372*nul[1]*fl[14]-1.185854122563142*nul[0]*fl[14]; 
  incr2[24] = (-5.929270612815711*nul[7]*fl[24])-4.592793267718455*nul[1]*fl[24]-2.651650429449552*nul[0]*fl[24]-4.592793267718456*nul[7]*fl[19]-3.557562367689425*nul[1]*fl[19]-2.053959590644372*nul[0]*fl[19]-2.651650429449552*nul[7]*fl[16]-2.053959590644372*nul[1]*fl[16]-1.185854122563142*nul[0]*fl[16]; 
  incr2[25] = 1.530931089239485*nul[7]*fl[26]+1.185854122563142*nul[1]*fl[26]+0.6846531968814573*nul[0]*fl[26]+1.185854122563142*nul[7]*fl[25]+0.9185586535436913*nul[1]*fl[25]+0.5303300858899105*nul[0]*fl[25]+0.6846531968814573*nul[7]*fl[22]+0.5303300858899105*nul[1]*fl[22]+0.3061862178478971*nul[0]*fl[22]; 
  incr2[26] = (-5.929270612815711*nul[7]*fl[26])-4.592793267718455*nul[1]*fl[26]-2.651650429449552*nul[0]*fl[26]-4.592793267718456*nul[7]*fl[25]-3.557562367689425*nul[1]*fl[25]-2.053959590644372*nul[0]*fl[25]-2.651650429449552*nul[7]*fl[22]-2.053959590644372*nul[1]*fl[22]-1.185854122563142*nul[0]*fl[22]; 

  outl[1] += -1.0*incr2[1]*rdxFl; 
  outl[4] += -1.0*incr2[4]*rdxFl; 
  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[7] += incr2[7]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[11] += incr2[11]*rdxFl; 
  outl[12] += -1.0*incr2[12]*rdxFl; 
  outl[13] += incr2[13]*rdxFl; 
  outl[15] += -1.0*incr2[15]*rdxFl; 
  outl[17] += incr2[17]*rdxFl; 
  outl[18] += -1.0*incr2[18]*rdxFl; 
  outl[19] += -1.0*incr2[19]*rdxFl; 
  outl[20] += incr2[20]*rdxFl; 
  outl[21] += incr2[21]*rdxFl; 
  outl[23] += incr2[23]*rdxFl; 
  outl[24] += incr2[24]*rdxFl; 
  outl[25] += -1.0*incr2[25]*rdxFl; 
  outl[26] += incr2[26]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf3xTensorP2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[81]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[1]*dxl[1]); 
  double rdxFr = 4.0/(dxr[1]*dxr[1]); 

  double incr1[27]; 
  double incr2[27]; 

  if (edge < 0) {

  incr2[2] = 1.530931089239485*fr[8]*nul[35]-1.185854122563142*fr[2]*nul[35]+0.6846531968814573*fr[0]*nul[35]+1.185854122563142*fr[8]*nul[29]-0.9185586535436913*fr[2]*nul[29]+0.5303300858899105*fr[0]*nul[29]+0.6846531968814573*fr[8]*nul[27]-0.5303300858899105*fr[2]*nul[27]+0.3061862178478971*fr[0]*nul[27]; 
  incr2[4] = 1.530931089239486*fr[12]*nul[35]-1.185854122563142*fr[4]*nul[35]+0.6846531968814573*fr[1]*nul[35]+1.185854122563142*fr[12]*nul[29]-0.9185586535436913*fr[4]*nul[29]+0.5303300858899105*fr[1]*nul[29]+0.6846531968814574*fr[12]*nul[27]-0.5303300858899105*fr[4]*nul[27]+0.3061862178478971*fr[1]*nul[27]; 
  incr2[6] = 1.530931089239486*fr[14]*nul[35]-1.185854122563142*fr[6]*nul[35]+0.6846531968814573*fr[3]*nul[35]+1.185854122563142*fr[14]*nul[29]-0.9185586535436913*fr[6]*nul[29]+0.5303300858899105*fr[3]*nul[29]+0.6846531968814574*fr[14]*nul[27]-0.5303300858899105*fr[6]*nul[27]+0.3061862178478971*fr[3]*nul[27]; 
  incr2[8] = (-5.929270612815711*fr[8]*nul[35])+4.592793267718456*fr[2]*nul[35]-2.651650429449552*fr[0]*nul[35]-4.592793267718455*fr[8]*nul[29]+3.557562367689425*fr[2]*nul[29]-2.053959590644372*fr[0]*nul[29]-2.651650429449552*fr[8]*nul[27]+2.053959590644372*fr[2]*nul[27]-1.185854122563142*fr[0]*nul[27]; 
  incr2[10] = 1.530931089239485*fr[18]*nul[35]-1.185854122563142*fr[10]*nul[35]+0.6846531968814573*fr[5]*nul[35]+1.185854122563142*fr[18]*nul[29]-0.9185586535436913*fr[10]*nul[29]+0.5303300858899105*fr[5]*nul[29]+0.6846531968814573*fr[18]*nul[27]-0.5303300858899105*fr[10]*nul[27]+0.3061862178478971*fr[5]*nul[27]; 
  incr2[11] = 1.530931089239486*fr[20]*nul[35]-1.185854122563142*fr[11]*nul[35]+0.6846531968814574*fr[7]*nul[35]+1.185854122563142*fr[20]*nul[29]-0.9185586535436913*fr[11]*nul[29]+0.5303300858899104*fr[7]*nul[29]+0.6846531968814575*fr[20]*nul[27]-0.5303300858899105*fr[11]*nul[27]+0.3061862178478971*fr[7]*nul[27]; 
  incr2[12] = (-5.929270612815711*fr[12]*nul[35])+4.592793267718458*fr[4]*nul[35]-2.651650429449552*fr[1]*nul[35]-4.592793267718455*fr[12]*nul[29]+3.557562367689425*fr[4]*nul[29]-2.053959590644372*fr[1]*nul[29]-2.651650429449552*fr[12]*nul[27]+2.053959590644372*fr[4]*nul[27]-1.185854122563142*fr[1]*nul[27]; 
  incr2[14] = (-5.929270612815711*fr[14]*nul[35])+4.592793267718458*fr[6]*nul[35]-2.651650429449552*fr[3]*nul[35]-4.592793267718455*fr[14]*nul[29]+3.557562367689425*fr[6]*nul[29]-2.053959590644372*fr[3]*nul[29]-2.651650429449552*fr[14]*nul[27]+2.053959590644372*fr[6]*nul[27]-1.185854122563142*fr[3]*nul[27]; 
  incr2[16] = 1.530931089239486*fr[22]*nul[35]-1.185854122563142*fr[16]*nul[35]+0.6846531968814574*fr[9]*nul[35]+1.185854122563142*fr[22]*nul[29]-0.9185586535436913*fr[16]*nul[29]+0.5303300858899104*fr[9]*nul[29]+0.6846531968814575*fr[22]*nul[27]-0.5303300858899105*fr[16]*nul[27]+0.3061862178478971*fr[9]*nul[27]; 
  incr2[17] = 1.530931089239485*fr[23]*nul[35]-1.185854122563142*fr[17]*nul[35]+0.6846531968814574*fr[13]*nul[35]+1.185854122563142*fr[23]*nul[29]-0.9185586535436913*fr[17]*nul[29]+0.5303300858899104*fr[13]*nul[29]+0.6846531968814573*fr[23]*nul[27]-0.5303300858899105*fr[17]*nul[27]+0.3061862178478971*fr[13]*nul[27]; 
  incr2[18] = (-5.929270612815711*fr[18]*nul[35])+4.592793267718456*fr[10]*nul[35]-2.651650429449552*fr[5]*nul[35]-4.592793267718455*fr[18]*nul[29]+3.557562367689425*fr[10]*nul[29]-2.053959590644372*fr[5]*nul[29]-2.651650429449552*fr[18]*nul[27]+2.053959590644372*fr[10]*nul[27]-1.185854122563142*fr[5]*nul[27]; 
  incr2[19] = 1.530931089239485*fr[25]*nul[35]-1.185854122563142*fr[19]*nul[35]+0.6846531968814574*fr[15]*nul[35]+1.185854122563142*fr[25]*nul[29]-0.9185586535436913*fr[19]*nul[29]+0.5303300858899104*fr[15]*nul[29]+0.6846531968814573*fr[25]*nul[27]-0.5303300858899105*fr[19]*nul[27]+0.3061862178478971*fr[15]*nul[27]; 
  incr2[20] = (-5.929270612815711*fr[20]*nul[35])+4.592793267718458*fr[11]*nul[35]-2.651650429449552*fr[7]*nul[35]-4.592793267718455*fr[20]*nul[29]+3.557562367689425*fr[11]*nul[29]-2.053959590644372*fr[7]*nul[29]-2.651650429449552*fr[20]*nul[27]+2.053959590644372*fr[11]*nul[27]-1.185854122563142*fr[7]*nul[27]; 
  incr2[22] = (-5.929270612815711*fr[22]*nul[35])+4.592793267718458*fr[16]*nul[35]-2.651650429449552*fr[9]*nul[35]-4.592793267718455*fr[22]*nul[29]+3.557562367689425*fr[16]*nul[29]-2.053959590644372*fr[9]*nul[29]-2.651650429449552*fr[22]*nul[27]+2.053959590644372*fr[16]*nul[27]-1.185854122563142*fr[9]*nul[27]; 
  incr2[23] = (-5.929270612815711*fr[23]*nul[35])+4.592793267718456*fr[17]*nul[35]-2.651650429449552*fr[13]*nul[35]-4.592793267718455*fr[23]*nul[29]+3.557562367689425*fr[17]*nul[29]-2.053959590644372*fr[13]*nul[29]-2.651650429449552*fr[23]*nul[27]+2.053959590644372*fr[17]*nul[27]-1.185854122563142*fr[13]*nul[27]; 
  incr2[24] = 1.530931089239485*fr[26]*nul[35]-1.185854122563142*fr[24]*nul[35]+0.6846531968814573*fr[21]*nul[35]+1.185854122563142*fr[26]*nul[29]-0.9185586535436913*fr[24]*nul[29]+0.5303300858899105*fr[21]*nul[29]+0.6846531968814573*fr[26]*nul[27]-0.5303300858899105*fr[24]*nul[27]+0.3061862178478971*fr[21]*nul[27]; 
  incr2[25] = (-5.929270612815711*fr[25]*nul[35])+4.592793267718456*fr[19]*nul[35]-2.651650429449552*fr[15]*nul[35]-4.592793267718455*fr[25]*nul[29]+3.557562367689425*fr[19]*nul[29]-2.053959590644372*fr[15]*nul[29]-2.651650429449552*fr[25]*nul[27]+2.053959590644372*fr[19]*nul[27]-1.185854122563142*fr[15]*nul[27]; 
  incr2[26] = (-5.929270612815711*fr[26]*nul[35])+4.592793267718456*fr[24]*nul[35]-2.651650429449552*fr[21]*nul[35]-4.592793267718455*fr[26]*nul[29]+3.557562367689425*fr[24]*nul[29]-2.053959590644372*fr[21]*nul[29]-2.651650429449552*fr[26]*nul[27]+2.053959590644372*fr[24]*nul[27]-1.185854122563142*fr[21]*nul[27]; 

  outr[2] += incr2[2]*rdxFr; 
  outr[4] += incr2[4]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[8] += incr2[8]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[11] += incr2[11]*rdxFr; 
  outr[12] += incr2[12]*rdxFr; 
  outr[14] += incr2[14]*rdxFr; 
  outr[16] += incr2[16]*rdxFr; 
  outr[17] += incr2[17]*rdxFr; 
  outr[18] += incr2[18]*rdxFr; 
  outr[19] += incr2[19]*rdxFr; 
  outr[20] += incr2[20]*rdxFr; 
  outr[22] += incr2[22]*rdxFr; 
  outr[23] += incr2[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr; 
  outr[25] += incr2[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 

  } else {

  incr2[2] = 1.530931089239485*fl[8]*nul[35]+1.185854122563142*fl[2]*nul[35]+0.6846531968814573*fl[0]*nul[35]+1.185854122563142*fl[8]*nul[29]+0.9185586535436913*fl[2]*nul[29]+0.5303300858899105*fl[0]*nul[29]+0.6846531968814573*fl[8]*nul[27]+0.5303300858899105*fl[2]*nul[27]+0.3061862178478971*fl[0]*nul[27]; 
  incr2[4] = 1.530931089239486*fl[12]*nul[35]+1.185854122563142*fl[4]*nul[35]+0.6846531968814573*fl[1]*nul[35]+1.185854122563142*fl[12]*nul[29]+0.9185586535436913*fl[4]*nul[29]+0.5303300858899105*fl[1]*nul[29]+0.6846531968814574*fl[12]*nul[27]+0.5303300858899105*fl[4]*nul[27]+0.3061862178478971*fl[1]*nul[27]; 
  incr2[6] = 1.530931089239486*fl[14]*nul[35]+1.185854122563142*fl[6]*nul[35]+0.6846531968814573*fl[3]*nul[35]+1.185854122563142*fl[14]*nul[29]+0.9185586535436913*fl[6]*nul[29]+0.5303300858899105*fl[3]*nul[29]+0.6846531968814574*fl[14]*nul[27]+0.5303300858899105*fl[6]*nul[27]+0.3061862178478971*fl[3]*nul[27]; 
  incr2[8] = (-5.929270612815711*fl[8]*nul[35])-4.592793267718456*fl[2]*nul[35]-2.651650429449552*fl[0]*nul[35]-4.592793267718455*fl[8]*nul[29]-3.557562367689425*fl[2]*nul[29]-2.053959590644372*fl[0]*nul[29]-2.651650429449552*fl[8]*nul[27]-2.053959590644372*fl[2]*nul[27]-1.185854122563142*fl[0]*nul[27]; 
  incr2[10] = 1.530931089239485*fl[18]*nul[35]+1.185854122563142*fl[10]*nul[35]+0.6846531968814573*fl[5]*nul[35]+1.185854122563142*fl[18]*nul[29]+0.9185586535436913*fl[10]*nul[29]+0.5303300858899105*fl[5]*nul[29]+0.6846531968814573*fl[18]*nul[27]+0.5303300858899105*fl[10]*nul[27]+0.3061862178478971*fl[5]*nul[27]; 
  incr2[11] = 1.530931089239486*fl[20]*nul[35]+1.185854122563142*fl[11]*nul[35]+0.6846531968814574*fl[7]*nul[35]+1.185854122563142*fl[20]*nul[29]+0.9185586535436913*fl[11]*nul[29]+0.5303300858899104*fl[7]*nul[29]+0.6846531968814575*fl[20]*nul[27]+0.5303300858899105*fl[11]*nul[27]+0.3061862178478971*fl[7]*nul[27]; 
  incr2[12] = (-5.929270612815711*fl[12]*nul[35])-4.592793267718458*fl[4]*nul[35]-2.651650429449552*fl[1]*nul[35]-4.592793267718455*fl[12]*nul[29]-3.557562367689425*fl[4]*nul[29]-2.053959590644372*fl[1]*nul[29]-2.651650429449552*fl[12]*nul[27]-2.053959590644372*fl[4]*nul[27]-1.185854122563142*fl[1]*nul[27]; 
  incr2[14] = (-5.929270612815711*fl[14]*nul[35])-4.592793267718458*fl[6]*nul[35]-2.651650429449552*fl[3]*nul[35]-4.592793267718455*fl[14]*nul[29]-3.557562367689425*fl[6]*nul[29]-2.053959590644372*fl[3]*nul[29]-2.651650429449552*fl[14]*nul[27]-2.053959590644372*fl[6]*nul[27]-1.185854122563142*fl[3]*nul[27]; 
  incr2[16] = 1.530931089239486*fl[22]*nul[35]+1.185854122563142*fl[16]*nul[35]+0.6846531968814574*fl[9]*nul[35]+1.185854122563142*fl[22]*nul[29]+0.9185586535436913*fl[16]*nul[29]+0.5303300858899104*fl[9]*nul[29]+0.6846531968814575*fl[22]*nul[27]+0.5303300858899105*fl[16]*nul[27]+0.3061862178478971*fl[9]*nul[27]; 
  incr2[17] = 1.530931089239485*fl[23]*nul[35]+1.185854122563142*fl[17]*nul[35]+0.6846531968814574*fl[13]*nul[35]+1.185854122563142*fl[23]*nul[29]+0.9185586535436913*fl[17]*nul[29]+0.5303300858899104*fl[13]*nul[29]+0.6846531968814573*fl[23]*nul[27]+0.5303300858899105*fl[17]*nul[27]+0.3061862178478971*fl[13]*nul[27]; 
  incr2[18] = (-5.929270612815711*fl[18]*nul[35])-4.592793267718456*fl[10]*nul[35]-2.651650429449552*fl[5]*nul[35]-4.592793267718455*fl[18]*nul[29]-3.557562367689425*fl[10]*nul[29]-2.053959590644372*fl[5]*nul[29]-2.651650429449552*fl[18]*nul[27]-2.053959590644372*fl[10]*nul[27]-1.185854122563142*fl[5]*nul[27]; 
  incr2[19] = 1.530931089239485*fl[25]*nul[35]+1.185854122563142*fl[19]*nul[35]+0.6846531968814574*fl[15]*nul[35]+1.185854122563142*fl[25]*nul[29]+0.9185586535436913*fl[19]*nul[29]+0.5303300858899104*fl[15]*nul[29]+0.6846531968814573*fl[25]*nul[27]+0.5303300858899105*fl[19]*nul[27]+0.3061862178478971*fl[15]*nul[27]; 
  incr2[20] = (-5.929270612815711*fl[20]*nul[35])-4.592793267718458*fl[11]*nul[35]-2.651650429449552*fl[7]*nul[35]-4.592793267718455*fl[20]*nul[29]-3.557562367689425*fl[11]*nul[29]-2.053959590644372*fl[7]*nul[29]-2.651650429449552*fl[20]*nul[27]-2.053959590644372*fl[11]*nul[27]-1.185854122563142*fl[7]*nul[27]; 
  incr2[22] = (-5.929270612815711*fl[22]*nul[35])-4.592793267718458*fl[16]*nul[35]-2.651650429449552*fl[9]*nul[35]-4.592793267718455*fl[22]*nul[29]-3.557562367689425*fl[16]*nul[29]-2.053959590644372*fl[9]*nul[29]-2.651650429449552*fl[22]*nul[27]-2.053959590644372*fl[16]*nul[27]-1.185854122563142*fl[9]*nul[27]; 
  incr2[23] = (-5.929270612815711*fl[23]*nul[35])-4.592793267718456*fl[17]*nul[35]-2.651650429449552*fl[13]*nul[35]-4.592793267718455*fl[23]*nul[29]-3.557562367689425*fl[17]*nul[29]-2.053959590644372*fl[13]*nul[29]-2.651650429449552*fl[23]*nul[27]-2.053959590644372*fl[17]*nul[27]-1.185854122563142*fl[13]*nul[27]; 
  incr2[24] = 1.530931089239485*fl[26]*nul[35]+1.185854122563142*fl[24]*nul[35]+0.6846531968814573*fl[21]*nul[35]+1.185854122563142*fl[26]*nul[29]+0.9185586535436913*fl[24]*nul[29]+0.5303300858899105*fl[21]*nul[29]+0.6846531968814573*fl[26]*nul[27]+0.5303300858899105*fl[24]*nul[27]+0.3061862178478971*fl[21]*nul[27]; 
  incr2[25] = (-5.929270612815711*fl[25]*nul[35])-4.592793267718456*fl[19]*nul[35]-2.651650429449552*fl[15]*nul[35]-4.592793267718455*fl[25]*nul[29]-3.557562367689425*fl[19]*nul[29]-2.053959590644372*fl[15]*nul[29]-2.651650429449552*fl[25]*nul[27]-2.053959590644372*fl[19]*nul[27]-1.185854122563142*fl[15]*nul[27]; 
  incr2[26] = (-5.929270612815711*fl[26]*nul[35])-4.592793267718456*fl[24]*nul[35]-2.651650429449552*fl[21]*nul[35]-4.592793267718455*fl[26]*nul[29]-3.557562367689425*fl[24]*nul[29]-2.053959590644372*fl[21]*nul[29]-2.651650429449552*fl[26]*nul[27]-2.053959590644372*fl[24]*nul[27]-1.185854122563142*fl[21]*nul[27]; 

  outl[2] += -1.0*incr2[2]*rdxFl; 
  outl[4] += -1.0*incr2[4]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[8] += incr2[8]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[11] += -1.0*incr2[11]*rdxFl; 
  outl[12] += incr2[12]*rdxFl; 
  outl[14] += incr2[14]*rdxFl; 
  outl[16] += -1.0*incr2[16]*rdxFl; 
  outl[17] += -1.0*incr2[17]*rdxFl; 
  outl[18] += incr2[18]*rdxFl; 
  outl[19] += -1.0*incr2[19]*rdxFl; 
  outl[20] += incr2[20]*rdxFl; 
  outl[22] += incr2[22]*rdxFl; 
  outl[23] += incr2[23]*rdxFl; 
  outl[24] += -1.0*incr2[24]*rdxFl; 
  outl[25] += incr2[25]*rdxFl; 
  outl[26] += incr2[26]*rdxFl; 

  }

} 
void ConstDiffusionVarCoeffBoundarySurf3xTensorP2_X3(const double *wl, const double *wr, const double *dxl, const double *dxr, const int *idxl, const int *idxr, const int edge, const double *nul, const double *nur, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:      Cell-center coordinates.
  // dx[3]:     Cell spacing.
  // idx[vdim+cdim]:    current grid index.
  // edge:      =-1 for lower boundary, =1 for upper boundary.
  // nu[81]:     diffusion coefficient.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFl = 4.0/(dxl[2]*dxl[2]); 
  double rdxFr = 4.0/(dxr[2]*dxr[2]); 

  double incr1[27]; 
  double incr2[27]; 

  if (edge < 0) {

  incr2[3] = 1.530931089239485*fr[9]*nul[63]-1.185854122563142*fr[3]*nul[63]+0.6846531968814573*fr[0]*nul[63]+1.185854122563142*fr[9]*nul[57]-0.9185586535436913*fr[3]*nul[57]+0.5303300858899105*fr[0]*nul[57]+0.6846531968814573*fr[9]*nul[54]-0.5303300858899105*fr[3]*nul[54]+0.3061862178478971*fr[0]*nul[54]; 
  incr2[5] = 1.530931089239486*fr[15]*nul[63]-1.185854122563142*fr[5]*nul[63]+0.6846531968814573*fr[1]*nul[63]+1.185854122563142*fr[15]*nul[57]-0.9185586535436913*fr[5]*nul[57]+0.5303300858899105*fr[1]*nul[57]+0.6846531968814574*fr[15]*nul[54]-0.5303300858899105*fr[5]*nul[54]+0.3061862178478971*fr[1]*nul[54]; 
  incr2[6] = 1.530931089239486*fr[16]*nul[63]-1.185854122563142*fr[6]*nul[63]+0.6846531968814573*fr[2]*nul[63]+1.185854122563142*fr[16]*nul[57]-0.9185586535436913*fr[6]*nul[57]+0.5303300858899105*fr[2]*nul[57]+0.6846531968814574*fr[16]*nul[54]-0.5303300858899105*fr[6]*nul[54]+0.3061862178478971*fr[2]*nul[54]; 
  incr2[9] = (-5.929270612815711*fr[9]*nul[63])+4.592793267718456*fr[3]*nul[63]-2.651650429449552*fr[0]*nul[63]-4.592793267718455*fr[9]*nul[57]+3.557562367689425*fr[3]*nul[57]-2.053959590644372*fr[0]*nul[57]-2.651650429449552*fr[9]*nul[54]+2.053959590644372*fr[3]*nul[54]-1.185854122563142*fr[0]*nul[54]; 
  incr2[10] = 1.530931089239485*fr[19]*nul[63]-1.185854122563142*fr[10]*nul[63]+0.6846531968814573*fr[4]*nul[63]+1.185854122563142*fr[19]*nul[57]-0.9185586535436913*fr[10]*nul[57]+0.5303300858899105*fr[4]*nul[57]+0.6846531968814573*fr[19]*nul[54]-0.5303300858899105*fr[10]*nul[54]+0.3061862178478971*fr[4]*nul[54]; 
  incr2[13] = 1.530931089239486*fr[21]*nul[63]-1.185854122563142*fr[13]*nul[63]+0.6846531968814574*fr[7]*nul[63]+1.185854122563142*fr[21]*nul[57]-0.9185586535436913*fr[13]*nul[57]+0.5303300858899104*fr[7]*nul[57]+0.6846531968814575*fr[21]*nul[54]-0.5303300858899105*fr[13]*nul[54]+0.3061862178478971*fr[7]*nul[54]; 
  incr2[14] = 1.530931089239486*fr[22]*nul[63]-1.185854122563142*fr[14]*nul[63]+0.6846531968814574*fr[8]*nul[63]+1.185854122563142*fr[22]*nul[57]-0.9185586535436913*fr[14]*nul[57]+0.5303300858899104*fr[8]*nul[57]+0.6846531968814575*fr[22]*nul[54]-0.5303300858899105*fr[14]*nul[54]+0.3061862178478971*fr[8]*nul[54]; 
  incr2[15] = (-5.929270612815711*fr[15]*nul[63])+4.592793267718458*fr[5]*nul[63]-2.651650429449552*fr[1]*nul[63]-4.592793267718455*fr[15]*nul[57]+3.557562367689425*fr[5]*nul[57]-2.053959590644372*fr[1]*nul[57]-2.651650429449552*fr[15]*nul[54]+2.053959590644372*fr[5]*nul[54]-1.185854122563142*fr[1]*nul[54]; 
  incr2[16] = (-5.929270612815711*fr[16]*nul[63])+4.592793267718458*fr[6]*nul[63]-2.651650429449552*fr[2]*nul[63]-4.592793267718455*fr[16]*nul[57]+3.557562367689425*fr[6]*nul[57]-2.053959590644372*fr[2]*nul[57]-2.651650429449552*fr[16]*nul[54]+2.053959590644372*fr[6]*nul[54]-1.185854122563142*fr[2]*nul[54]; 
  incr2[17] = 1.530931089239485*fr[24]*nul[63]-1.185854122563142*fr[17]*nul[63]+0.6846531968814574*fr[11]*nul[63]+1.185854122563142*fr[24]*nul[57]-0.9185586535436913*fr[17]*nul[57]+0.5303300858899104*fr[11]*nul[57]+0.6846531968814573*fr[24]*nul[54]-0.5303300858899105*fr[17]*nul[54]+0.3061862178478971*fr[11]*nul[54]; 
  incr2[18] = 1.530931089239485*fr[25]*nul[63]-1.185854122563142*fr[18]*nul[63]+0.6846531968814574*fr[12]*nul[63]+1.185854122563142*fr[25]*nul[57]-0.9185586535436913*fr[18]*nul[57]+0.5303300858899104*fr[12]*nul[57]+0.6846531968814573*fr[25]*nul[54]-0.5303300858899105*fr[18]*nul[54]+0.3061862178478971*fr[12]*nul[54]; 
  incr2[19] = (-5.929270612815711*fr[19]*nul[63])+4.592793267718456*fr[10]*nul[63]-2.651650429449552*fr[4]*nul[63]-4.592793267718455*fr[19]*nul[57]+3.557562367689425*fr[10]*nul[57]-2.053959590644372*fr[4]*nul[57]-2.651650429449552*fr[19]*nul[54]+2.053959590644372*fr[10]*nul[54]-1.185854122563142*fr[4]*nul[54]; 
  incr2[21] = (-5.929270612815711*fr[21]*nul[63])+4.592793267718458*fr[13]*nul[63]-2.651650429449552*fr[7]*nul[63]-4.592793267718455*fr[21]*nul[57]+3.557562367689425*fr[13]*nul[57]-2.053959590644372*fr[7]*nul[57]-2.651650429449552*fr[21]*nul[54]+2.053959590644372*fr[13]*nul[54]-1.185854122563142*fr[7]*nul[54]; 
  incr2[22] = (-5.929270612815711*fr[22]*nul[63])+4.592793267718458*fr[14]*nul[63]-2.651650429449552*fr[8]*nul[63]-4.592793267718455*fr[22]*nul[57]+3.557562367689425*fr[14]*nul[57]-2.053959590644372*fr[8]*nul[57]-2.651650429449552*fr[22]*nul[54]+2.053959590644372*fr[14]*nul[54]-1.185854122563142*fr[8]*nul[54]; 
  incr2[23] = 1.530931089239485*fr[26]*nul[63]-1.185854122563142*fr[23]*nul[63]+0.6846531968814573*fr[20]*nul[63]+1.185854122563142*fr[26]*nul[57]-0.9185586535436913*fr[23]*nul[57]+0.5303300858899105*fr[20]*nul[57]+0.6846531968814573*fr[26]*nul[54]-0.5303300858899105*fr[23]*nul[54]+0.3061862178478971*fr[20]*nul[54]; 
  incr2[24] = (-5.929270612815711*fr[24]*nul[63])+4.592793267718456*fr[17]*nul[63]-2.651650429449552*fr[11]*nul[63]-4.592793267718455*fr[24]*nul[57]+3.557562367689425*fr[17]*nul[57]-2.053959590644372*fr[11]*nul[57]-2.651650429449552*fr[24]*nul[54]+2.053959590644372*fr[17]*nul[54]-1.185854122563142*fr[11]*nul[54]; 
  incr2[25] = (-5.929270612815711*fr[25]*nul[63])+4.592793267718456*fr[18]*nul[63]-2.651650429449552*fr[12]*nul[63]-4.592793267718455*fr[25]*nul[57]+3.557562367689425*fr[18]*nul[57]-2.053959590644372*fr[12]*nul[57]-2.651650429449552*fr[25]*nul[54]+2.053959590644372*fr[18]*nul[54]-1.185854122563142*fr[12]*nul[54]; 
  incr2[26] = (-5.929270612815711*fr[26]*nul[63])+4.592793267718456*fr[23]*nul[63]-2.651650429449552*fr[20]*nul[63]-4.592793267718455*fr[26]*nul[57]+3.557562367689425*fr[23]*nul[57]-2.053959590644372*fr[20]*nul[57]-2.651650429449552*fr[26]*nul[54]+2.053959590644372*fr[23]*nul[54]-1.185854122563142*fr[20]*nul[54]; 

  outr[3] += incr2[3]*rdxFr; 
  outr[5] += incr2[5]*rdxFr; 
  outr[6] += incr2[6]*rdxFr; 
  outr[9] += incr2[9]*rdxFr; 
  outr[10] += incr2[10]*rdxFr; 
  outr[13] += incr2[13]*rdxFr; 
  outr[14] += incr2[14]*rdxFr; 
  outr[15] += incr2[15]*rdxFr; 
  outr[16] += incr2[16]*rdxFr; 
  outr[17] += incr2[17]*rdxFr; 
  outr[18] += incr2[18]*rdxFr; 
  outr[19] += incr2[19]*rdxFr; 
  outr[21] += incr2[21]*rdxFr; 
  outr[22] += incr2[22]*rdxFr; 
  outr[23] += incr2[23]*rdxFr; 
  outr[24] += incr2[24]*rdxFr; 
  outr[25] += incr2[25]*rdxFr; 
  outr[26] += incr2[26]*rdxFr; 

  } else {

  incr2[3] = 1.530931089239485*fl[9]*nul[63]+1.185854122563142*fl[3]*nul[63]+0.6846531968814573*fl[0]*nul[63]+1.185854122563142*fl[9]*nul[57]+0.9185586535436913*fl[3]*nul[57]+0.5303300858899105*fl[0]*nul[57]+0.6846531968814573*fl[9]*nul[54]+0.5303300858899105*fl[3]*nul[54]+0.3061862178478971*fl[0]*nul[54]; 
  incr2[5] = 1.530931089239486*fl[15]*nul[63]+1.185854122563142*fl[5]*nul[63]+0.6846531968814573*fl[1]*nul[63]+1.185854122563142*fl[15]*nul[57]+0.9185586535436913*fl[5]*nul[57]+0.5303300858899105*fl[1]*nul[57]+0.6846531968814574*fl[15]*nul[54]+0.5303300858899105*fl[5]*nul[54]+0.3061862178478971*fl[1]*nul[54]; 
  incr2[6] = 1.530931089239486*fl[16]*nul[63]+1.185854122563142*fl[6]*nul[63]+0.6846531968814573*fl[2]*nul[63]+1.185854122563142*fl[16]*nul[57]+0.9185586535436913*fl[6]*nul[57]+0.5303300858899105*fl[2]*nul[57]+0.6846531968814574*fl[16]*nul[54]+0.5303300858899105*fl[6]*nul[54]+0.3061862178478971*fl[2]*nul[54]; 
  incr2[9] = (-5.929270612815711*fl[9]*nul[63])-4.592793267718456*fl[3]*nul[63]-2.651650429449552*fl[0]*nul[63]-4.592793267718455*fl[9]*nul[57]-3.557562367689425*fl[3]*nul[57]-2.053959590644372*fl[0]*nul[57]-2.651650429449552*fl[9]*nul[54]-2.053959590644372*fl[3]*nul[54]-1.185854122563142*fl[0]*nul[54]; 
  incr2[10] = 1.530931089239485*fl[19]*nul[63]+1.185854122563142*fl[10]*nul[63]+0.6846531968814573*fl[4]*nul[63]+1.185854122563142*fl[19]*nul[57]+0.9185586535436913*fl[10]*nul[57]+0.5303300858899105*fl[4]*nul[57]+0.6846531968814573*fl[19]*nul[54]+0.5303300858899105*fl[10]*nul[54]+0.3061862178478971*fl[4]*nul[54]; 
  incr2[13] = 1.530931089239486*fl[21]*nul[63]+1.185854122563142*fl[13]*nul[63]+0.6846531968814574*fl[7]*nul[63]+1.185854122563142*fl[21]*nul[57]+0.9185586535436913*fl[13]*nul[57]+0.5303300858899104*fl[7]*nul[57]+0.6846531968814575*fl[21]*nul[54]+0.5303300858899105*fl[13]*nul[54]+0.3061862178478971*fl[7]*nul[54]; 
  incr2[14] = 1.530931089239486*fl[22]*nul[63]+1.185854122563142*fl[14]*nul[63]+0.6846531968814574*fl[8]*nul[63]+1.185854122563142*fl[22]*nul[57]+0.9185586535436913*fl[14]*nul[57]+0.5303300858899104*fl[8]*nul[57]+0.6846531968814575*fl[22]*nul[54]+0.5303300858899105*fl[14]*nul[54]+0.3061862178478971*fl[8]*nul[54]; 
  incr2[15] = (-5.929270612815711*fl[15]*nul[63])-4.592793267718458*fl[5]*nul[63]-2.651650429449552*fl[1]*nul[63]-4.592793267718455*fl[15]*nul[57]-3.557562367689425*fl[5]*nul[57]-2.053959590644372*fl[1]*nul[57]-2.651650429449552*fl[15]*nul[54]-2.053959590644372*fl[5]*nul[54]-1.185854122563142*fl[1]*nul[54]; 
  incr2[16] = (-5.929270612815711*fl[16]*nul[63])-4.592793267718458*fl[6]*nul[63]-2.651650429449552*fl[2]*nul[63]-4.592793267718455*fl[16]*nul[57]-3.557562367689425*fl[6]*nul[57]-2.053959590644372*fl[2]*nul[57]-2.651650429449552*fl[16]*nul[54]-2.053959590644372*fl[6]*nul[54]-1.185854122563142*fl[2]*nul[54]; 
  incr2[17] = 1.530931089239485*fl[24]*nul[63]+1.185854122563142*fl[17]*nul[63]+0.6846531968814574*fl[11]*nul[63]+1.185854122563142*fl[24]*nul[57]+0.9185586535436913*fl[17]*nul[57]+0.5303300858899104*fl[11]*nul[57]+0.6846531968814573*fl[24]*nul[54]+0.5303300858899105*fl[17]*nul[54]+0.3061862178478971*fl[11]*nul[54]; 
  incr2[18] = 1.530931089239485*fl[25]*nul[63]+1.185854122563142*fl[18]*nul[63]+0.6846531968814574*fl[12]*nul[63]+1.185854122563142*fl[25]*nul[57]+0.9185586535436913*fl[18]*nul[57]+0.5303300858899104*fl[12]*nul[57]+0.6846531968814573*fl[25]*nul[54]+0.5303300858899105*fl[18]*nul[54]+0.3061862178478971*fl[12]*nul[54]; 
  incr2[19] = (-5.929270612815711*fl[19]*nul[63])-4.592793267718456*fl[10]*nul[63]-2.651650429449552*fl[4]*nul[63]-4.592793267718455*fl[19]*nul[57]-3.557562367689425*fl[10]*nul[57]-2.053959590644372*fl[4]*nul[57]-2.651650429449552*fl[19]*nul[54]-2.053959590644372*fl[10]*nul[54]-1.185854122563142*fl[4]*nul[54]; 
  incr2[21] = (-5.929270612815711*fl[21]*nul[63])-4.592793267718458*fl[13]*nul[63]-2.651650429449552*fl[7]*nul[63]-4.592793267718455*fl[21]*nul[57]-3.557562367689425*fl[13]*nul[57]-2.053959590644372*fl[7]*nul[57]-2.651650429449552*fl[21]*nul[54]-2.053959590644372*fl[13]*nul[54]-1.185854122563142*fl[7]*nul[54]; 
  incr2[22] = (-5.929270612815711*fl[22]*nul[63])-4.592793267718458*fl[14]*nul[63]-2.651650429449552*fl[8]*nul[63]-4.592793267718455*fl[22]*nul[57]-3.557562367689425*fl[14]*nul[57]-2.053959590644372*fl[8]*nul[57]-2.651650429449552*fl[22]*nul[54]-2.053959590644372*fl[14]*nul[54]-1.185854122563142*fl[8]*nul[54]; 
  incr2[23] = 1.530931089239485*fl[26]*nul[63]+1.185854122563142*fl[23]*nul[63]+0.6846531968814573*fl[20]*nul[63]+1.185854122563142*fl[26]*nul[57]+0.9185586535436913*fl[23]*nul[57]+0.5303300858899105*fl[20]*nul[57]+0.6846531968814573*fl[26]*nul[54]+0.5303300858899105*fl[23]*nul[54]+0.3061862178478971*fl[20]*nul[54]; 
  incr2[24] = (-5.929270612815711*fl[24]*nul[63])-4.592793267718456*fl[17]*nul[63]-2.651650429449552*fl[11]*nul[63]-4.592793267718455*fl[24]*nul[57]-3.557562367689425*fl[17]*nul[57]-2.053959590644372*fl[11]*nul[57]-2.651650429449552*fl[24]*nul[54]-2.053959590644372*fl[17]*nul[54]-1.185854122563142*fl[11]*nul[54]; 
  incr2[25] = (-5.929270612815711*fl[25]*nul[63])-4.592793267718456*fl[18]*nul[63]-2.651650429449552*fl[12]*nul[63]-4.592793267718455*fl[25]*nul[57]-3.557562367689425*fl[18]*nul[57]-2.053959590644372*fl[12]*nul[57]-2.651650429449552*fl[25]*nul[54]-2.053959590644372*fl[18]*nul[54]-1.185854122563142*fl[12]*nul[54]; 
  incr2[26] = (-5.929270612815711*fl[26]*nul[63])-4.592793267718456*fl[23]*nul[63]-2.651650429449552*fl[20]*nul[63]-4.592793267718455*fl[26]*nul[57]-3.557562367689425*fl[23]*nul[57]-2.053959590644372*fl[20]*nul[57]-2.651650429449552*fl[26]*nul[54]-2.053959590644372*fl[23]*nul[54]-1.185854122563142*fl[20]*nul[54]; 

  outl[3] += -1.0*incr2[3]*rdxFl; 
  outl[5] += -1.0*incr2[5]*rdxFl; 
  outl[6] += -1.0*incr2[6]*rdxFl; 
  outl[9] += incr2[9]*rdxFl; 
  outl[10] += -1.0*incr2[10]*rdxFl; 
  outl[13] += -1.0*incr2[13]*rdxFl; 
  outl[14] += -1.0*incr2[14]*rdxFl; 
  outl[15] += incr2[15]*rdxFl; 
  outl[16] += incr2[16]*rdxFl; 
  outl[17] += -1.0*incr2[17]*rdxFl; 
  outl[18] += -1.0*incr2[18]*rdxFl; 
  outl[19] += incr2[19]*rdxFl; 
  outl[21] += incr2[21]*rdxFl; 
  outl[22] += incr2[22]*rdxFl; 
  outl[23] += -1.0*incr2[23]*rdxFl; 
  outl[24] += incr2[24]*rdxFl; 
  outl[25] += incr2[25]*rdxFl; 
  outl[26] += incr2[26]*rdxFl; 

  }

} 
