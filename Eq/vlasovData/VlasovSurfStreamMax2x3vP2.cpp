#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream2x3vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[0]; 
  double rdxr2 = 2.0/dxvr[0]; 

  double incr[21]; 

  if (wr[2]>0) { 
  incr[0] = wl[2]*(1.118033988749895*fl[16]+0.8660254037844386*fl[1]+0.5*fl[0])+dxvl[2]*(0.25*fl[7]+0.1443375672974065*fl[3]); 
  incr[1] = wl[2]*((-1.936491673103709*fl[16])-1.5*fl[1]-0.8660254037844386*fl[0])+dxvl[2]*((-0.4330127018922193*fl[7])-0.25*fl[3]); 
  incr[2] = 0.1443375672974065*dxvl[2]*fl[8]+wl[2]*(0.8660254037844386*fl[6]+0.5*fl[2]); 
  incr[3] = dxvl[2]*(0.1290994448735806*fl[18]+0.3227486121839515*fl[16]+0.25*fl[1]+0.1443375672974065*fl[0])+wl[2]*(0.8660254037844386*fl[7]+0.5*fl[3]); 
  incr[4] = 0.1443375672974065*dxvl[2]*fl[11]+wl[2]*(0.8660254037844386*fl[9]+0.5*fl[4]); 
  incr[5] = 0.1443375672974065*dxvl[2]*fl[14]+wl[2]*(0.8660254037844386*fl[12]+0.5*fl[5]); 
  incr[6] = wl[2]*((-1.5*fl[6])-0.8660254037844386*fl[2])-0.25*dxvl[2]*fl[8]; 
  incr[7] = dxvl[2]*((-0.223606797749979*fl[18])-0.5590169943749475*fl[16]-0.4330127018922193*fl[1]-0.25*fl[0])+wl[2]*((-1.5*fl[7])-0.8660254037844386*fl[3]); 
  incr[8] = 0.5*wl[2]*fl[8]+dxvl[2]*(0.25*fl[6]+0.1443375672974065*fl[2]); 
  incr[9] = wl[2]*((-1.5*fl[9])-0.8660254037844386*fl[4])-0.25*dxvl[2]*fl[11]; 
  incr[10] = 0.5*wl[2]*fl[10]; 
  incr[11] = 0.5*wl[2]*fl[11]+dxvl[2]*(0.25*fl[9]+0.1443375672974065*fl[4]); 
  incr[12] = wl[2]*((-1.5*fl[12])-0.8660254037844386*fl[5])-0.25*dxvl[2]*fl[14]; 
  incr[13] = 0.5*wl[2]*fl[13]; 
  incr[14] = 0.5*wl[2]*fl[14]+dxvl[2]*(0.25*fl[12]+0.1443375672974065*fl[5]); 
  incr[15] = 0.5*wl[2]*fl[15]; 
  incr[16] = wl[2]*(2.5*fl[16]+1.936491673103709*fl[1]+1.118033988749895*fl[0])+dxvl[2]*(0.5590169943749475*fl[7]+0.3227486121839515*fl[3]); 
  incr[17] = 0.5*wl[2]*fl[17]; 
  incr[18] = 0.5*wl[2]*fl[18]+dxvl[2]*(0.223606797749979*fl[7]+0.1290994448735806*fl[3]); 
  incr[19] = 0.5*wl[2]*fl[19]; 
  incr[20] = 0.5*wl[2]*fl[20]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 
  outr[8] += incr[8]*rdxr2; 
  outr[9] += incr[9]*rdxr2; 
  outr[10] += incr[10]*rdxr2; 
  outr[11] += incr[11]*rdxr2; 
  outr[12] += incr[12]*rdxr2; 
  outr[13] += incr[13]*rdxr2; 
  outr[14] += incr[14]*rdxr2; 
  outr[15] += incr[15]*rdxr2; 
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  outl[8] += -1.0*incr[8]*rdxl2; 
  outl[9] += incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += -1.0*incr[16]*rdxl2; 
  outl[17] += -1.0*incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += -1.0*incr[20]*rdxl2; 
  } else { 
  incr[0] = wr[2]*(1.118033988749895*fr[16]-0.8660254037844386*fr[1]+0.5*fr[0])+dxvr[2]*(0.1443375672974065*fr[3]-0.25*fr[7]); 
  incr[1] = wr[2]*((-1.936491673103709*fr[16])+1.5*fr[1]-0.8660254037844386*fr[0])+dxvr[2]*(0.4330127018922193*fr[7]-0.25*fr[3]); 
  incr[2] = 0.1443375672974065*dxvr[2]*fr[8]+wr[2]*(0.5*fr[2]-0.8660254037844386*fr[6]); 
  incr[3] = dxvr[2]*(0.1290994448735806*fr[18]+0.3227486121839515*fr[16]-0.25*fr[1]+0.1443375672974065*fr[0])+wr[2]*(0.5*fr[3]-0.8660254037844386*fr[7]); 
  incr[4] = 0.1443375672974065*dxvr[2]*fr[11]+wr[2]*(0.5*fr[4]-0.8660254037844386*fr[9]); 
  incr[5] = 0.1443375672974065*dxvr[2]*fr[14]+wr[2]*(0.5*fr[5]-0.8660254037844386*fr[12]); 
  incr[6] = wr[2]*(1.5*fr[6]-0.8660254037844386*fr[2])-0.25*dxvr[2]*fr[8]; 
  incr[7] = dxvr[2]*((-0.223606797749979*fr[18])-0.5590169943749475*fr[16]+0.4330127018922193*fr[1]-0.25*fr[0])+wr[2]*(1.5*fr[7]-0.8660254037844386*fr[3]); 
  incr[8] = 0.5*wr[2]*fr[8]+dxvr[2]*(0.1443375672974065*fr[2]-0.25*fr[6]); 
  incr[9] = wr[2]*(1.5*fr[9]-0.8660254037844386*fr[4])-0.25*dxvr[2]*fr[11]; 
  incr[10] = 0.5*wr[2]*fr[10]; 
  incr[11] = 0.5*wr[2]*fr[11]+dxvr[2]*(0.1443375672974065*fr[4]-0.25*fr[9]); 
  incr[12] = wr[2]*(1.5*fr[12]-0.8660254037844386*fr[5])-0.25*dxvr[2]*fr[14]; 
  incr[13] = 0.5*wr[2]*fr[13]; 
  incr[14] = 0.5*wr[2]*fr[14]+dxvr[2]*(0.1443375672974065*fr[5]-0.25*fr[12]); 
  incr[15] = 0.5*wr[2]*fr[15]; 
  incr[16] = wr[2]*(2.5*fr[16]-1.936491673103709*fr[1]+1.118033988749895*fr[0])+dxvr[2]*(0.3227486121839515*fr[3]-0.5590169943749475*fr[7]); 
  incr[17] = 0.5*wr[2]*fr[17]; 
  incr[18] = 0.5*wr[2]*fr[18]+dxvr[2]*(0.1290994448735806*fr[3]-0.223606797749979*fr[7]); 
  incr[19] = 0.5*wr[2]*fr[19]; 
  incr[20] = 0.5*wr[2]*fr[20]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 
  outr[8] += incr[8]*rdxr2; 
  outr[9] += incr[9]*rdxr2; 
  outr[10] += incr[10]*rdxr2; 
  outr[11] += incr[11]*rdxr2; 
  outr[12] += incr[12]*rdxr2; 
  outr[13] += incr[13]*rdxr2; 
  outr[14] += incr[14]*rdxr2; 
  outr[15] += incr[15]*rdxr2; 
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += incr[1]*rdxl2; 
  outl[2] += -1.0*incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += incr[7]*rdxl2; 
  outl[8] += -1.0*incr[8]*rdxl2; 
  outl[9] += incr[9]*rdxl2; 
  outl[10] += -1.0*incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += incr[12]*rdxl2; 
  outl[13] += -1.0*incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += -1.0*incr[16]*rdxl2; 
  outl[17] += -1.0*incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += -1.0*incr[20]*rdxl2; 
  } 
} 
__host__ __device__ void VlasovSurfStream2x3vMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells.
  double rdxl2 = 2.0/dxvl[1]; 
  double rdxr2 = 2.0/dxvr[1]; 

  double incr[21]; 

  if (wr[3]>0) { 
  incr[0] = wl[3]*(1.118033988749895*fl[17]+0.8660254037844386*fl[2]+0.5*fl[0])+dxvl[3]*(0.25*fl[10]+0.1443375672974065*fl[4]); 
  incr[1] = 0.1443375672974065*dxvl[3]*fl[9]+wl[3]*(0.8660254037844386*fl[6]+0.5*fl[1]); 
  incr[2] = wl[3]*((-1.936491673103709*fl[17])-1.5*fl[2]-0.8660254037844386*fl[0])+dxvl[3]*((-0.4330127018922193*fl[10])-0.25*fl[4]); 
  incr[3] = 0.1443375672974065*dxvl[3]*fl[11]+wl[3]*(0.8660254037844386*fl[8]+0.5*fl[3]); 
  incr[4] = dxvl[3]*(0.1290994448735806*fl[19]+0.3227486121839515*fl[17]+0.25*fl[2]+0.1443375672974065*fl[0])+wl[3]*(0.8660254037844386*fl[10]+0.5*fl[4]); 
  incr[5] = 0.1443375672974065*dxvl[3]*fl[15]+wl[3]*(0.8660254037844386*fl[13]+0.5*fl[5]); 
  incr[6] = wl[3]*((-1.5*fl[6])-0.8660254037844386*fl[1])-0.25*dxvl[3]*fl[9]; 
  incr[7] = 0.5*wl[3]*fl[7]; 
  incr[8] = wl[3]*((-1.5*fl[8])-0.8660254037844386*fl[3])-0.25*dxvl[3]*fl[11]; 
  incr[9] = 0.5*wl[3]*fl[9]+dxvl[3]*(0.25*fl[6]+0.1443375672974065*fl[1]); 
  incr[10] = dxvl[3]*((-0.223606797749979*fl[19])-0.5590169943749475*fl[17]-0.4330127018922193*fl[2]-0.25*fl[0])+wl[3]*((-1.5*fl[10])-0.8660254037844386*fl[4]); 
  incr[11] = 0.5*wl[3]*fl[11]+dxvl[3]*(0.25*fl[8]+0.1443375672974065*fl[3]); 
  incr[12] = 0.5*wl[3]*fl[12]; 
  incr[13] = wl[3]*((-1.5*fl[13])-0.8660254037844386*fl[5])-0.25*dxvl[3]*fl[15]; 
  incr[14] = 0.5*wl[3]*fl[14]; 
  incr[15] = 0.5*wl[3]*fl[15]+dxvl[3]*(0.25*fl[13]+0.1443375672974065*fl[5]); 
  incr[16] = 0.5*wl[3]*fl[16]; 
  incr[17] = wl[3]*(2.5*fl[17]+1.936491673103709*fl[2]+1.118033988749895*fl[0])+dxvl[3]*(0.5590169943749475*fl[10]+0.3227486121839515*fl[4]); 
  incr[18] = 0.5*wl[3]*fl[18]; 
  incr[19] = 0.5*wl[3]*fl[19]+dxvl[3]*(0.223606797749979*fl[10]+0.1290994448735806*fl[4]); 
  incr[20] = 0.5*wl[3]*fl[20]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 
  outr[8] += incr[8]*rdxr2; 
  outr[9] += incr[9]*rdxr2; 
  outr[10] += incr[10]*rdxr2; 
  outr[11] += incr[11]*rdxr2; 
  outr[12] += incr[12]*rdxr2; 
  outr[13] += incr[13]*rdxr2; 
  outr[14] += incr[14]*rdxr2; 
  outr[15] += incr[15]*rdxr2; 
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += -1.0*incr[12]*rdxl2; 
  outl[13] += incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += -1.0*incr[16]*rdxl2; 
  outl[17] += -1.0*incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += -1.0*incr[20]*rdxl2; 
  } else { 
  incr[0] = wr[3]*(1.118033988749895*fr[17]-0.8660254037844386*fr[2]+0.5*fr[0])+dxvr[3]*(0.1443375672974065*fr[4]-0.25*fr[10]); 
  incr[1] = 0.1443375672974065*dxvr[3]*fr[9]+wr[3]*(0.5*fr[1]-0.8660254037844386*fr[6]); 
  incr[2] = wr[3]*((-1.936491673103709*fr[17])+1.5*fr[2]-0.8660254037844386*fr[0])+dxvr[3]*(0.4330127018922193*fr[10]-0.25*fr[4]); 
  incr[3] = 0.1443375672974065*dxvr[3]*fr[11]+wr[3]*(0.5*fr[3]-0.8660254037844386*fr[8]); 
  incr[4] = dxvr[3]*(0.1290994448735806*fr[19]+0.3227486121839515*fr[17]-0.25*fr[2]+0.1443375672974065*fr[0])+wr[3]*(0.5*fr[4]-0.8660254037844386*fr[10]); 
  incr[5] = 0.1443375672974065*dxvr[3]*fr[15]+wr[3]*(0.5*fr[5]-0.8660254037844386*fr[13]); 
  incr[6] = wr[3]*(1.5*fr[6]-0.8660254037844386*fr[1])-0.25*dxvr[3]*fr[9]; 
  incr[7] = 0.5*wr[3]*fr[7]; 
  incr[8] = wr[3]*(1.5*fr[8]-0.8660254037844386*fr[3])-0.25*dxvr[3]*fr[11]; 
  incr[9] = 0.5*wr[3]*fr[9]+dxvr[3]*(0.1443375672974065*fr[1]-0.25*fr[6]); 
  incr[10] = dxvr[3]*((-0.223606797749979*fr[19])-0.5590169943749475*fr[17]+0.4330127018922193*fr[2]-0.25*fr[0])+wr[3]*(1.5*fr[10]-0.8660254037844386*fr[4]); 
  incr[11] = 0.5*wr[3]*fr[11]+dxvr[3]*(0.1443375672974065*fr[3]-0.25*fr[8]); 
  incr[12] = 0.5*wr[3]*fr[12]; 
  incr[13] = wr[3]*(1.5*fr[13]-0.8660254037844386*fr[5])-0.25*dxvr[3]*fr[15]; 
  incr[14] = 0.5*wr[3]*fr[14]; 
  incr[15] = 0.5*wr[3]*fr[15]+dxvr[3]*(0.1443375672974065*fr[5]-0.25*fr[13]); 
  incr[16] = 0.5*wr[3]*fr[16]; 
  incr[17] = wr[3]*(2.5*fr[17]-1.936491673103709*fr[2]+1.118033988749895*fr[0])+dxvr[3]*(0.3227486121839515*fr[4]-0.5590169943749475*fr[10]); 
  incr[18] = 0.5*wr[3]*fr[18]; 
  incr[19] = 0.5*wr[3]*fr[19]+dxvr[3]*(0.1290994448735806*fr[4]-0.223606797749979*fr[10]); 
  incr[20] = 0.5*wr[3]*fr[20]; 

  outr[0] += incr[0]*rdxr2; 
  outr[1] += incr[1]*rdxr2; 
  outr[2] += incr[2]*rdxr2; 
  outr[3] += incr[3]*rdxr2; 
  outr[4] += incr[4]*rdxr2; 
  outr[5] += incr[5]*rdxr2; 
  outr[6] += incr[6]*rdxr2; 
  outr[7] += incr[7]*rdxr2; 
  outr[8] += incr[8]*rdxr2; 
  outr[9] += incr[9]*rdxr2; 
  outr[10] += incr[10]*rdxr2; 
  outr[11] += incr[11]*rdxr2; 
  outr[12] += incr[12]*rdxr2; 
  outr[13] += incr[13]*rdxr2; 
  outr[14] += incr[14]*rdxr2; 
  outr[15] += incr[15]*rdxr2; 
  outr[16] += incr[16]*rdxr2; 
  outr[17] += incr[17]*rdxr2; 
  outr[18] += incr[18]*rdxr2; 
  outr[19] += incr[19]*rdxr2; 
  outr[20] += incr[20]*rdxr2; 

  outl[0] += -1.0*incr[0]*rdxl2; 
  outl[1] += -1.0*incr[1]*rdxl2; 
  outl[2] += incr[2]*rdxl2; 
  outl[3] += -1.0*incr[3]*rdxl2; 
  outl[4] += -1.0*incr[4]*rdxl2; 
  outl[5] += -1.0*incr[5]*rdxl2; 
  outl[6] += incr[6]*rdxl2; 
  outl[7] += -1.0*incr[7]*rdxl2; 
  outl[8] += incr[8]*rdxl2; 
  outl[9] += -1.0*incr[9]*rdxl2; 
  outl[10] += incr[10]*rdxl2; 
  outl[11] += -1.0*incr[11]*rdxl2; 
  outl[12] += -1.0*incr[12]*rdxl2; 
  outl[13] += incr[13]*rdxl2; 
  outl[14] += -1.0*incr[14]*rdxl2; 
  outl[15] += -1.0*incr[15]*rdxl2; 
  outl[16] += -1.0*incr[16]*rdxl2; 
  outl[17] += -1.0*incr[17]*rdxl2; 
  outl[18] += -1.0*incr[18]*rdxl2; 
  outl[19] += -1.0*incr[19]*rdxl2; 
  outl[20] += -1.0*incr[20]*rdxl2; 
  } 
} 
