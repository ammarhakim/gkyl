#include <VlasovModDecl.h> 
void VlasovSurfStream3x3vMax_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[3]/dxv[0]; 
  double wx = w[3]*2/dxv[0]; 

  double incr[7]; 

  if (wx>0) { 
  incr[0] = 0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.2886751345948129*fl[4]*dvx; 
  incr[1] = (-1.5*fl[1]*wx)-0.8660254037844386*fl[0]*wx-0.5*fl[4]*dvx; 
  incr[2] = 0.5*fl[2]*wx; 
  incr[3] = 0.5*fl[3]*wx; 
  incr[4] = 0.5*fl[4]*wx+0.5*fl[1]*dvx+0.2886751345948129*fl[0]*dvx; 
  incr[5] = 0.5*fl[5]*wx; 
  incr[6] = 0.5*fl[6]*wx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  } else { 
  incr[0] = (-0.8660254037844386*fr[1]*wx)+0.5*fr[0]*wx+0.2886751345948129*fr[4]*dvx; 
  incr[1] = 1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx-0.5*fr[4]*dvx; 
  incr[2] = 0.5*fr[2]*wx; 
  incr[3] = 0.5*fr[3]*wx; 
  incr[4] = 0.5*fr[4]*wx-0.5*fr[1]*dvx+0.2886751345948129*fr[0]*dvx; 
  incr[5] = 0.5*fr[5]*wx; 
  incr[6] = 0.5*fr[6]*wx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  } 
} 
void VlasovSurfStream3x3vMax_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[3]/dxv[0]; 
  double wx = w[3]*2/dxv[0]; 

  double incr[28]; 

  if (wx>0) { 
  incr[0] = 1.118033988749895*fl[22]*wx+0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.5*fl[10]*dvx+0.2886751345948129*fl[4]*dvx; 
  incr[1] = (-1.936491673103709*fl[22]*wx)-1.5*fl[1]*wx-0.8660254037844386*fl[0]*wx-0.8660254037844386*fl[10]*dvx-0.5*fl[4]*dvx; 
  incr[2] = 0.8660254037844386*fl[7]*wx+0.5*fl[2]*wx+0.2886751345948129*fl[11]*dvx; 
  incr[3] = 0.8660254037844386*fl[8]*wx+0.5*fl[3]*wx+0.2886751345948129*fl[12]*dvx; 
  incr[4] = 0.8660254037844386*fl[10]*wx+0.5*fl[4]*wx+0.2581988897471612*fl[25]*dvx+0.6454972243679029*fl[22]*dvx+0.5*fl[1]*dvx+0.2886751345948129*fl[0]*dvx; 
  incr[5] = 0.8660254037844386*fl[13]*wx+0.5*fl[5]*wx+0.2886751345948129*fl[16]*dvx; 
  incr[6] = 0.8660254037844386*fl[17]*wx+0.5*fl[6]*wx+0.2886751345948129*fl[20]*dvx; 
  incr[7] = (-1.5*fl[7]*wx)-0.8660254037844386*fl[2]*wx-0.5*fl[11]*dvx; 
  incr[8] = (-1.5*fl[8]*wx)-0.8660254037844386*fl[3]*wx-0.5*fl[12]*dvx; 
  incr[9] = 0.5*fl[9]*wx; 
  incr[10] = (-1.5*fl[10]*wx)-0.8660254037844386*fl[4]*wx-0.4472135954999579*fl[25]*dvx-1.118033988749895*fl[22]*dvx-0.8660254037844386*fl[1]*dvx-0.5*fl[0]*dvx; 
  incr[11] = 0.5*fl[11]*wx+0.5*fl[7]*dvx+0.2886751345948129*fl[2]*dvx; 
  incr[12] = 0.5*fl[12]*wx+0.5*fl[8]*dvx+0.2886751345948129*fl[3]*dvx; 
  incr[13] = (-1.5*fl[13]*wx)-0.8660254037844386*fl[5]*wx-0.5*fl[16]*dvx; 
  incr[14] = 0.5*fl[14]*wx; 
  incr[15] = 0.5*fl[15]*wx; 
  incr[16] = 0.5*fl[16]*wx+0.5*fl[13]*dvx+0.2886751345948129*fl[5]*dvx; 
  incr[17] = (-1.5*fl[17]*wx)-0.8660254037844386*fl[6]*wx-0.5*fl[20]*dvx; 
  incr[18] = 0.5*fl[18]*wx; 
  incr[19] = 0.5*fl[19]*wx; 
  incr[20] = 0.5*fl[20]*wx+0.5*fl[17]*dvx+0.2886751345948129*fl[6]*dvx; 
  incr[21] = 0.5*fl[21]*wx; 
  incr[22] = 2.5*fl[22]*wx+1.936491673103709*fl[1]*wx+1.118033988749895*fl[0]*wx+1.118033988749895*fl[10]*dvx+0.6454972243679029*fl[4]*dvx; 
  incr[23] = 0.5*fl[23]*wx; 
  incr[24] = 0.5*fl[24]*wx; 
  incr[25] = 0.5*fl[25]*wx+0.4472135954999579*fl[10]*dvx+0.2581988897471612*fl[4]*dvx; 
  incr[26] = 0.5*fl[26]*wx; 
  incr[27] = 0.5*fl[27]*wx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += -1.0*incr[19]; 
  outl[20] += -1.0*incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += -1.0*incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += -1.0*incr[26]; 
  outl[27] += -1.0*incr[27]; 
  } else { 
  incr[0] = 1.118033988749895*fr[22]*wx-0.8660254037844386*fr[1]*wx+0.5*fr[0]*wx-0.5*fr[10]*dvx+0.2886751345948129*fr[4]*dvx; 
  incr[1] = (-1.936491673103709*fr[22]*wx)+1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx+0.8660254037844386*fr[10]*dvx-0.5*fr[4]*dvx; 
  incr[2] = (-0.8660254037844386*fr[7]*wx)+0.5*fr[2]*wx+0.2886751345948129*fr[11]*dvx; 
  incr[3] = (-0.8660254037844386*fr[8]*wx)+0.5*fr[3]*wx+0.2886751345948129*fr[12]*dvx; 
  incr[4] = (-0.8660254037844386*fr[10]*wx)+0.5*fr[4]*wx+0.2581988897471612*fr[25]*dvx+0.6454972243679029*fr[22]*dvx-0.5*fr[1]*dvx+0.2886751345948129*fr[0]*dvx; 
  incr[5] = (-0.8660254037844386*fr[13]*wx)+0.5*fr[5]*wx+0.2886751345948129*fr[16]*dvx; 
  incr[6] = (-0.8660254037844386*fr[17]*wx)+0.5*fr[6]*wx+0.2886751345948129*fr[20]*dvx; 
  incr[7] = 1.5*fr[7]*wx-0.8660254037844386*fr[2]*wx-0.5*fr[11]*dvx; 
  incr[8] = 1.5*fr[8]*wx-0.8660254037844386*fr[3]*wx-0.5*fr[12]*dvx; 
  incr[9] = 0.5*fr[9]*wx; 
  incr[10] = 1.5*fr[10]*wx-0.8660254037844386*fr[4]*wx-0.4472135954999579*fr[25]*dvx-1.118033988749895*fr[22]*dvx+0.8660254037844386*fr[1]*dvx-0.5*fr[0]*dvx; 
  incr[11] = 0.5*fr[11]*wx-0.5*fr[7]*dvx+0.2886751345948129*fr[2]*dvx; 
  incr[12] = 0.5*fr[12]*wx-0.5*fr[8]*dvx+0.2886751345948129*fr[3]*dvx; 
  incr[13] = 1.5*fr[13]*wx-0.8660254037844386*fr[5]*wx-0.5*fr[16]*dvx; 
  incr[14] = 0.5*fr[14]*wx; 
  incr[15] = 0.5*fr[15]*wx; 
  incr[16] = 0.5*fr[16]*wx-0.5*fr[13]*dvx+0.2886751345948129*fr[5]*dvx; 
  incr[17] = 1.5*fr[17]*wx-0.8660254037844386*fr[6]*wx-0.5*fr[20]*dvx; 
  incr[18] = 0.5*fr[18]*wx; 
  incr[19] = 0.5*fr[19]*wx; 
  incr[20] = 0.5*fr[20]*wx-0.5*fr[17]*dvx+0.2886751345948129*fr[6]*dvx; 
  incr[21] = 0.5*fr[21]*wx; 
  incr[22] = 2.5*fr[22]*wx-1.936491673103709*fr[1]*wx+1.118033988749895*fr[0]*wx-1.118033988749895*fr[10]*dvx+0.6454972243679029*fr[4]*dvx; 
  incr[23] = 0.5*fr[23]*wx; 
  incr[24] = 0.5*fr[24]*wx; 
  incr[25] = 0.5*fr[25]*wx-0.4472135954999579*fr[10]*dvx+0.2581988897471612*fr[4]*dvx; 
  incr[26] = 0.5*fr[26]*wx; 
  incr[27] = 0.5*fr[27]*wx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += -1.0*incr[19]; 
  outl[20] += -1.0*incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += -1.0*incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += -1.0*incr[26]; 
  outl[27] += -1.0*incr[27]; 
  } 
} 
void VlasovSurfStream3x3vMax_Y_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[4]/dxv[1]; 
  double wx = w[4]*2/dxv[1]; 

  double incr[7]; 

  if (wx>0) { 
  incr[0] = 0.8660254037844386*fl[2]*wx+0.5*fl[0]*wx+0.2886751345948129*fl[5]*dvx; 
  incr[1] = 0.5*fl[1]*wx; 
  incr[2] = (-1.5*fl[2]*wx)-0.8660254037844386*fl[0]*wx-0.5*fl[5]*dvx; 
  incr[3] = 0.5*fl[3]*wx; 
  incr[4] = 0.5*fl[4]*wx; 
  incr[5] = 0.5*fl[5]*wx+0.5*fl[2]*dvx+0.2886751345948129*fl[0]*dvx; 
  incr[6] = 0.5*fl[6]*wx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  } else { 
  incr[0] = (-0.8660254037844386*fr[2]*wx)+0.5*fr[0]*wx+0.2886751345948129*fr[5]*dvx; 
  incr[1] = 0.5*fr[1]*wx; 
  incr[2] = 1.5*fr[2]*wx-0.8660254037844386*fr[0]*wx-0.5*fr[5]*dvx; 
  incr[3] = 0.5*fr[3]*wx; 
  incr[4] = 0.5*fr[4]*wx; 
  incr[5] = 0.5*fr[5]*wx-0.5*fr[2]*dvx+0.2886751345948129*fr[0]*dvx; 
  incr[6] = 0.5*fr[6]*wx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  } 
} 
void VlasovSurfStream3x3vMax_Y_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[4]/dxv[1]; 
  double wx = w[4]*2/dxv[1]; 

  double incr[28]; 

  if (wx>0) { 
  incr[0] = 1.118033988749895*fl[23]*wx+0.8660254037844386*fl[2]*wx+0.5*fl[0]*wx+0.5*fl[14]*dvx+0.2886751345948129*fl[5]*dvx; 
  incr[1] = 0.8660254037844386*fl[7]*wx+0.5*fl[1]*wx+0.2886751345948129*fl[13]*dvx; 
  incr[2] = (-1.936491673103709*fl[23]*wx)-1.5*fl[2]*wx-0.8660254037844386*fl[0]*wx-0.8660254037844386*fl[14]*dvx-0.5*fl[5]*dvx; 
  incr[3] = 0.8660254037844386*fl[9]*wx+0.5*fl[3]*wx+0.2886751345948129*fl[15]*dvx; 
  incr[4] = 0.8660254037844386*fl[11]*wx+0.5*fl[4]*wx+0.2886751345948129*fl[16]*dvx; 
  incr[5] = 0.8660254037844386*fl[14]*wx+0.5*fl[5]*wx+0.2581988897471612*fl[26]*dvx+0.6454972243679029*fl[23]*dvx+0.5*fl[2]*dvx+0.2886751345948129*fl[0]*dvx; 
  incr[6] = 0.8660254037844386*fl[18]*wx+0.5*fl[6]*wx+0.2886751345948129*fl[21]*dvx; 
  incr[7] = (-1.5*fl[7]*wx)-0.8660254037844386*fl[1]*wx-0.5*fl[13]*dvx; 
  incr[8] = 0.5*fl[8]*wx; 
  incr[9] = (-1.5*fl[9]*wx)-0.8660254037844386*fl[3]*wx-0.5*fl[15]*dvx; 
  incr[10] = 0.5*fl[10]*wx; 
  incr[11] = (-1.5*fl[11]*wx)-0.8660254037844386*fl[4]*wx-0.5*fl[16]*dvx; 
  incr[12] = 0.5*fl[12]*wx; 
  incr[13] = 0.5*fl[13]*wx+0.5*fl[7]*dvx+0.2886751345948129*fl[1]*dvx; 
  incr[14] = (-1.5*fl[14]*wx)-0.8660254037844386*fl[5]*wx-0.4472135954999579*fl[26]*dvx-1.118033988749895*fl[23]*dvx-0.8660254037844386*fl[2]*dvx-0.5*fl[0]*dvx; 
  incr[15] = 0.5*fl[15]*wx+0.5*fl[9]*dvx+0.2886751345948129*fl[3]*dvx; 
  incr[16] = 0.5*fl[16]*wx+0.5*fl[11]*dvx+0.2886751345948129*fl[4]*dvx; 
  incr[17] = 0.5*fl[17]*wx; 
  incr[18] = (-1.5*fl[18]*wx)-0.8660254037844386*fl[6]*wx-0.5*fl[21]*dvx; 
  incr[19] = 0.5*fl[19]*wx; 
  incr[20] = 0.5*fl[20]*wx; 
  incr[21] = 0.5*fl[21]*wx+0.5*fl[18]*dvx+0.2886751345948129*fl[6]*dvx; 
  incr[22] = 0.5*fl[22]*wx; 
  incr[23] = 2.5*fl[23]*wx+1.936491673103709*fl[2]*wx+1.118033988749895*fl[0]*wx+1.118033988749895*fl[14]*dvx+0.6454972243679029*fl[5]*dvx; 
  incr[24] = 0.5*fl[24]*wx; 
  incr[25] = 0.5*fl[25]*wx; 
  incr[26] = 0.5*fl[26]*wx+0.4472135954999579*fl[14]*dvx+0.2581988897471612*fl[5]*dvx; 
  incr[27] = 0.5*fl[27]*wx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += -1.0*incr[19]; 
  outl[20] += -1.0*incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += -1.0*incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += -1.0*incr[26]; 
  outl[27] += -1.0*incr[27]; 
  } else { 
  incr[0] = 1.118033988749895*fr[23]*wx-0.8660254037844386*fr[2]*wx+0.5*fr[0]*wx-0.5*fr[14]*dvx+0.2886751345948129*fr[5]*dvx; 
  incr[1] = (-0.8660254037844386*fr[7]*wx)+0.5*fr[1]*wx+0.2886751345948129*fr[13]*dvx; 
  incr[2] = (-1.936491673103709*fr[23]*wx)+1.5*fr[2]*wx-0.8660254037844386*fr[0]*wx+0.8660254037844386*fr[14]*dvx-0.5*fr[5]*dvx; 
  incr[3] = (-0.8660254037844386*fr[9]*wx)+0.5*fr[3]*wx+0.2886751345948129*fr[15]*dvx; 
  incr[4] = (-0.8660254037844386*fr[11]*wx)+0.5*fr[4]*wx+0.2886751345948129*fr[16]*dvx; 
  incr[5] = (-0.8660254037844386*fr[14]*wx)+0.5*fr[5]*wx+0.2581988897471612*fr[26]*dvx+0.6454972243679029*fr[23]*dvx-0.5*fr[2]*dvx+0.2886751345948129*fr[0]*dvx; 
  incr[6] = (-0.8660254037844386*fr[18]*wx)+0.5*fr[6]*wx+0.2886751345948129*fr[21]*dvx; 
  incr[7] = 1.5*fr[7]*wx-0.8660254037844386*fr[1]*wx-0.5*fr[13]*dvx; 
  incr[8] = 0.5*fr[8]*wx; 
  incr[9] = 1.5*fr[9]*wx-0.8660254037844386*fr[3]*wx-0.5*fr[15]*dvx; 
  incr[10] = 0.5*fr[10]*wx; 
  incr[11] = 1.5*fr[11]*wx-0.8660254037844386*fr[4]*wx-0.5*fr[16]*dvx; 
  incr[12] = 0.5*fr[12]*wx; 
  incr[13] = 0.5*fr[13]*wx-0.5*fr[7]*dvx+0.2886751345948129*fr[1]*dvx; 
  incr[14] = 1.5*fr[14]*wx-0.8660254037844386*fr[5]*wx-0.4472135954999579*fr[26]*dvx-1.118033988749895*fr[23]*dvx+0.8660254037844386*fr[2]*dvx-0.5*fr[0]*dvx; 
  incr[15] = 0.5*fr[15]*wx-0.5*fr[9]*dvx+0.2886751345948129*fr[3]*dvx; 
  incr[16] = 0.5*fr[16]*wx-0.5*fr[11]*dvx+0.2886751345948129*fr[4]*dvx; 
  incr[17] = 0.5*fr[17]*wx; 
  incr[18] = 1.5*fr[18]*wx-0.8660254037844386*fr[6]*wx-0.5*fr[21]*dvx; 
  incr[19] = 0.5*fr[19]*wx; 
  incr[20] = 0.5*fr[20]*wx; 
  incr[21] = 0.5*fr[21]*wx-0.5*fr[18]*dvx+0.2886751345948129*fr[6]*dvx; 
  incr[22] = 0.5*fr[22]*wx; 
  incr[23] = 2.5*fr[23]*wx-1.936491673103709*fr[2]*wx+1.118033988749895*fr[0]*wx-1.118033988749895*fr[14]*dvx+0.6454972243679029*fr[5]*dvx; 
  incr[24] = 0.5*fr[24]*wx; 
  incr[25] = 0.5*fr[25]*wx; 
  incr[26] = 0.5*fr[26]*wx-0.4472135954999579*fr[14]*dvx+0.2581988897471612*fr[5]*dvx; 
  incr[27] = 0.5*fr[27]*wx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += -1.0*incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += incr[18]; 
  outl[19] += -1.0*incr[19]; 
  outl[20] += -1.0*incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += -1.0*incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += -1.0*incr[26]; 
  outl[27] += -1.0*incr[27]; 
  } 
} 
void VlasovSurfStream3x3vMax_Z_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[5]/dxv[2]; 
  double wx = w[5]*2/dxv[2]; 

  double incr[7]; 

  if (wx>0) { 
  incr[0] = 0.8660254037844386*fl[3]*wx+0.5*fl[0]*wx+0.2886751345948129*fl[6]*dvx; 
  incr[1] = 0.5*fl[1]*wx; 
  incr[2] = 0.5*fl[2]*wx; 
  incr[3] = (-1.5*fl[3]*wx)-0.8660254037844386*fl[0]*wx-0.5*fl[6]*dvx; 
  incr[4] = 0.5*fl[4]*wx; 
  incr[5] = 0.5*fl[5]*wx; 
  incr[6] = 0.5*fl[6]*wx+0.5*fl[3]*dvx+0.2886751345948129*fl[0]*dvx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  } else { 
  incr[0] = (-0.8660254037844386*fr[3]*wx)+0.5*fr[0]*wx+0.2886751345948129*fr[6]*dvx; 
  incr[1] = 0.5*fr[1]*wx; 
  incr[2] = 0.5*fr[2]*wx; 
  incr[3] = 1.5*fr[3]*wx-0.8660254037844386*fr[0]*wx-0.5*fr[6]*dvx; 
  incr[4] = 0.5*fr[4]*wx; 
  incr[5] = 0.5*fr[5]*wx; 
  incr[6] = 0.5*fr[6]*wx-0.5*fr[3]*dvx+0.2886751345948129*fr[0]*dvx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  } 
} 
void VlasovSurfStream3x3vMax_Z_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[5]/dxv[2]; 
  double wx = w[5]*2/dxv[2]; 

  double incr[28]; 

  if (wx>0) { 
  incr[0] = 1.118033988749895*fl[24]*wx+0.8660254037844386*fl[3]*wx+0.5*fl[0]*wx+0.5*fl[19]*dvx+0.2886751345948129*fl[6]*dvx; 
  incr[1] = 0.8660254037844386*fl[8]*wx+0.5*fl[1]*wx+0.2886751345948129*fl[17]*dvx; 
  incr[2] = 0.8660254037844386*fl[9]*wx+0.5*fl[2]*wx+0.2886751345948129*fl[18]*dvx; 
  incr[3] = (-1.936491673103709*fl[24]*wx)-1.5*fl[3]*wx-0.8660254037844386*fl[0]*wx-0.8660254037844386*fl[19]*dvx-0.5*fl[6]*dvx; 
  incr[4] = 0.8660254037844386*fl[12]*wx+0.5*fl[4]*wx+0.2886751345948129*fl[20]*dvx; 
  incr[5] = 0.8660254037844386*fl[15]*wx+0.5*fl[5]*wx+0.2886751345948129*fl[21]*dvx; 
  incr[6] = 0.8660254037844386*fl[19]*wx+0.5*fl[6]*wx+0.2581988897471612*fl[27]*dvx+0.6454972243679029*fl[24]*dvx+0.5*fl[3]*dvx+0.2886751345948129*fl[0]*dvx; 
  incr[7] = 0.5*fl[7]*wx; 
  incr[8] = (-1.5*fl[8]*wx)-0.8660254037844386*fl[1]*wx-0.5*fl[17]*dvx; 
  incr[9] = (-1.5*fl[9]*wx)-0.8660254037844386*fl[2]*wx-0.5*fl[18]*dvx; 
  incr[10] = 0.5*fl[10]*wx; 
  incr[11] = 0.5*fl[11]*wx; 
  incr[12] = (-1.5*fl[12]*wx)-0.8660254037844386*fl[4]*wx-0.5*fl[20]*dvx; 
  incr[13] = 0.5*fl[13]*wx; 
  incr[14] = 0.5*fl[14]*wx; 
  incr[15] = (-1.5*fl[15]*wx)-0.8660254037844386*fl[5]*wx-0.5*fl[21]*dvx; 
  incr[16] = 0.5*fl[16]*wx; 
  incr[17] = 0.5*fl[17]*wx+0.5*fl[8]*dvx+0.2886751345948129*fl[1]*dvx; 
  incr[18] = 0.5*fl[18]*wx+0.5*fl[9]*dvx+0.2886751345948129*fl[2]*dvx; 
  incr[19] = (-1.5*fl[19]*wx)-0.8660254037844386*fl[6]*wx-0.4472135954999579*fl[27]*dvx-1.118033988749895*fl[24]*dvx-0.8660254037844386*fl[3]*dvx-0.5*fl[0]*dvx; 
  incr[20] = 0.5*fl[20]*wx+0.5*fl[12]*dvx+0.2886751345948129*fl[4]*dvx; 
  incr[21] = 0.5*fl[21]*wx+0.5*fl[15]*dvx+0.2886751345948129*fl[5]*dvx; 
  incr[22] = 0.5*fl[22]*wx; 
  incr[23] = 0.5*fl[23]*wx; 
  incr[24] = 2.5*fl[24]*wx+1.936491673103709*fl[3]*wx+1.118033988749895*fl[0]*wx+1.118033988749895*fl[19]*dvx+0.6454972243679029*fl[6]*dvx; 
  incr[25] = 0.5*fl[25]*wx; 
  incr[26] = 0.5*fl[26]*wx; 
  incr[27] = 0.5*fl[27]*wx+0.4472135954999579*fl[19]*dvx+0.2581988897471612*fl[6]*dvx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += incr[19]; 
  outl[20] += -1.0*incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += -1.0*incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += -1.0*incr[26]; 
  outl[27] += -1.0*incr[27]; 
  } else { 
  incr[0] = 1.118033988749895*fr[24]*wx-0.8660254037844386*fr[3]*wx+0.5*fr[0]*wx-0.5*fr[19]*dvx+0.2886751345948129*fr[6]*dvx; 
  incr[1] = (-0.8660254037844386*fr[8]*wx)+0.5*fr[1]*wx+0.2886751345948129*fr[17]*dvx; 
  incr[2] = (-0.8660254037844386*fr[9]*wx)+0.5*fr[2]*wx+0.2886751345948129*fr[18]*dvx; 
  incr[3] = (-1.936491673103709*fr[24]*wx)+1.5*fr[3]*wx-0.8660254037844386*fr[0]*wx+0.8660254037844386*fr[19]*dvx-0.5*fr[6]*dvx; 
  incr[4] = (-0.8660254037844386*fr[12]*wx)+0.5*fr[4]*wx+0.2886751345948129*fr[20]*dvx; 
  incr[5] = (-0.8660254037844386*fr[15]*wx)+0.5*fr[5]*wx+0.2886751345948129*fr[21]*dvx; 
  incr[6] = (-0.8660254037844386*fr[19]*wx)+0.5*fr[6]*wx+0.2581988897471612*fr[27]*dvx+0.6454972243679029*fr[24]*dvx-0.5*fr[3]*dvx+0.2886751345948129*fr[0]*dvx; 
  incr[7] = 0.5*fr[7]*wx; 
  incr[8] = 1.5*fr[8]*wx-0.8660254037844386*fr[1]*wx-0.5*fr[17]*dvx; 
  incr[9] = 1.5*fr[9]*wx-0.8660254037844386*fr[2]*wx-0.5*fr[18]*dvx; 
  incr[10] = 0.5*fr[10]*wx; 
  incr[11] = 0.5*fr[11]*wx; 
  incr[12] = 1.5*fr[12]*wx-0.8660254037844386*fr[4]*wx-0.5*fr[20]*dvx; 
  incr[13] = 0.5*fr[13]*wx; 
  incr[14] = 0.5*fr[14]*wx; 
  incr[15] = 1.5*fr[15]*wx-0.8660254037844386*fr[5]*wx-0.5*fr[21]*dvx; 
  incr[16] = 0.5*fr[16]*wx; 
  incr[17] = 0.5*fr[17]*wx-0.5*fr[8]*dvx+0.2886751345948129*fr[1]*dvx; 
  incr[18] = 0.5*fr[18]*wx-0.5*fr[9]*dvx+0.2886751345948129*fr[2]*dvx; 
  incr[19] = 1.5*fr[19]*wx-0.8660254037844386*fr[6]*wx-0.4472135954999579*fr[27]*dvx-1.118033988749895*fr[24]*dvx+0.8660254037844386*fr[3]*dvx-0.5*fr[0]*dvx; 
  incr[20] = 0.5*fr[20]*wx-0.5*fr[12]*dvx+0.2886751345948129*fr[4]*dvx; 
  incr[21] = 0.5*fr[21]*wx-0.5*fr[15]*dvx+0.2886751345948129*fr[5]*dvx; 
  incr[22] = 0.5*fr[22]*wx; 
  incr[23] = 0.5*fr[23]*wx; 
  incr[24] = 2.5*fr[24]*wx-1.936491673103709*fr[3]*wx+1.118033988749895*fr[0]*wx-1.118033988749895*fr[19]*dvx+0.6454972243679029*fr[6]*dvx; 
  incr[25] = 0.5*fr[25]*wx; 
  incr[26] = 0.5*fr[26]*wx; 
  incr[27] = 0.5*fr[27]*wx-0.4472135954999579*fr[19]*dvx+0.2581988897471612*fr[6]*dvx; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 
  outr[16] += incr[16]; 
  outr[17] += incr[17]; 
  outr[18] += incr[18]; 
  outr[19] += incr[19]; 
  outr[20] += incr[20]; 
  outr[21] += incr[21]; 
  outr[22] += incr[22]; 
  outr[23] += incr[23]; 
  outr[24] += incr[24]; 
  outr[25] += incr[25]; 
  outr[26] += incr[26]; 
  outr[27] += incr[27]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += -1.0*incr[16]; 
  outl[17] += -1.0*incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += incr[19]; 
  outl[20] += -1.0*incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += -1.0*incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += -1.0*incr[26]; 
  outl[27] += -1.0*incr[27]; 
  } 
} 
