#include <VlasovModDecl.h> 
void VlasovSurfStream1x3vSer_X_P1(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[1]/dxv[0]; 
  double wx = w[1]*2/dxv[0]; 

  double incr[16]; 

  if (wx>0) { 
  incr[0] = 0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.5*fl[5]*dvx+0.2886751345948129*fl[2]*dvx; 
  incr[1] = (-1.5*fl[1]*wx)-0.8660254037844386*fl[0]*wx-0.8660254037844386*fl[5]*dvx-0.5*fl[2]*dvx; 
  incr[2] = 0.8660254037844386*fl[5]*wx+0.5*fl[2]*wx+0.5*fl[1]*dvx+0.2886751345948129*fl[0]*dvx; 
  incr[3] = 0.8660254037844386*fl[6]*wx+0.5*fl[3]*wx+0.5*fl[11]*dvx+0.2886751345948129*fl[7]*dvx; 
  incr[4] = 0.8660254037844386*fl[8]*wx+0.5*fl[4]*wx+0.5*fl[12]*dvx+0.2886751345948129*fl[9]*dvx; 
  incr[5] = (-1.5*fl[5]*wx)-0.8660254037844386*fl[2]*wx-0.8660254037844386*fl[1]*dvx-0.5*fl[0]*dvx; 
  incr[6] = (-1.5*fl[6]*wx)-0.8660254037844386*fl[3]*wx-0.8660254037844386*fl[11]*dvx-0.5*fl[7]*dvx; 
  incr[7] = 0.8660254037844386*fl[11]*wx+0.5*fl[7]*wx+0.5*fl[6]*dvx+0.2886751345948129*fl[3]*dvx; 
  incr[8] = (-1.5*fl[8]*wx)-0.8660254037844386*fl[4]*wx-0.8660254037844386*fl[12]*dvx-0.5*fl[9]*dvx; 
  incr[9] = 0.8660254037844386*fl[12]*wx+0.5*fl[9]*wx+0.5*fl[8]*dvx+0.2886751345948129*fl[4]*dvx; 
  incr[10] = 0.8660254037844386*fl[13]*wx+0.5*fl[10]*wx+0.5*fl[15]*dvx+0.2886751345948129*fl[14]*dvx; 
  incr[11] = (-1.5*fl[11]*wx)-0.8660254037844386*fl[7]*wx-0.8660254037844386*fl[6]*dvx-0.5*fl[3]*dvx; 
  incr[12] = (-1.5*fl[12]*wx)-0.8660254037844386*fl[9]*wx-0.8660254037844386*fl[8]*dvx-0.5*fl[4]*dvx; 
  incr[13] = (-1.5*fl[13]*wx)-0.8660254037844386*fl[10]*wx-0.8660254037844386*fl[15]*dvx-0.5*fl[14]*dvx; 
  incr[14] = 0.8660254037844386*fl[15]*wx+0.5*fl[14]*wx+0.5*fl[13]*dvx+0.2886751345948129*fl[10]*dvx; 
  incr[15] = (-1.5*fl[15]*wx)-0.8660254037844386*fl[14]*wx-0.8660254037844386*fl[13]*dvx-0.5*fl[10]*dvx; 

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

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = (-0.8660254037844386*fr[1]*wx)+0.5*fr[0]*wx-0.5*fr[5]*dvx+0.2886751345948129*fr[2]*dvx; 
  incr[1] = 1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx+0.8660254037844386*fr[5]*dvx-0.5*fr[2]*dvx; 
  incr[2] = (-0.8660254037844386*fr[5]*wx)+0.5*fr[2]*wx-0.5*fr[1]*dvx+0.2886751345948129*fr[0]*dvx; 
  incr[3] = (-0.8660254037844386*fr[6]*wx)+0.5*fr[3]*wx-0.5*fr[11]*dvx+0.2886751345948129*fr[7]*dvx; 
  incr[4] = (-0.8660254037844386*fr[8]*wx)+0.5*fr[4]*wx-0.5*fr[12]*dvx+0.2886751345948129*fr[9]*dvx; 
  incr[5] = 1.5*fr[5]*wx-0.8660254037844386*fr[2]*wx+0.8660254037844386*fr[1]*dvx-0.5*fr[0]*dvx; 
  incr[6] = 1.5*fr[6]*wx-0.8660254037844386*fr[3]*wx+0.8660254037844386*fr[11]*dvx-0.5*fr[7]*dvx; 
  incr[7] = (-0.8660254037844386*fr[11]*wx)+0.5*fr[7]*wx-0.5*fr[6]*dvx+0.2886751345948129*fr[3]*dvx; 
  incr[8] = 1.5*fr[8]*wx-0.8660254037844386*fr[4]*wx+0.8660254037844386*fr[12]*dvx-0.5*fr[9]*dvx; 
  incr[9] = (-0.8660254037844386*fr[12]*wx)+0.5*fr[9]*wx-0.5*fr[8]*dvx+0.2886751345948129*fr[4]*dvx; 
  incr[10] = (-0.8660254037844386*fr[13]*wx)+0.5*fr[10]*wx-0.5*fr[15]*dvx+0.2886751345948129*fr[14]*dvx; 
  incr[11] = 1.5*fr[11]*wx-0.8660254037844386*fr[7]*wx+0.8660254037844386*fr[6]*dvx-0.5*fr[3]*dvx; 
  incr[12] = 1.5*fr[12]*wx-0.8660254037844386*fr[9]*wx+0.8660254037844386*fr[8]*dvx-0.5*fr[4]*dvx; 
  incr[13] = 1.5*fr[13]*wx-0.8660254037844386*fr[10]*wx+0.8660254037844386*fr[15]*dvx-0.5*fr[14]*dvx; 
  incr[14] = (-0.8660254037844386*fr[15]*wx)+0.5*fr[14]*wx-0.5*fr[13]*dvx+0.2886751345948129*fr[10]*dvx; 
  incr[15] = 1.5*fr[15]*wx-0.8660254037844386*fr[14]*wx+0.8660254037844386*fr[13]*dvx-0.5*fr[10]*dvx; 

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

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  } 
} 
void VlasovSurfStream1x3vSer_X_P2(const double *w, const double *dxv, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvx = dxv[1]/dxv[0]; 
  double wx = w[1]*2/dxv[0]; 

  double incr[48]; 

  if (wx>0) { 
  incr[0] = 1.118033988749895*fl[11]*wx+0.8660254037844386*fl[1]*wx+0.5*fl[0]*wx+0.6454972243679029*fl[19]*dvx+0.5*fl[5]*dvx+0.2886751345948129*fl[2]*dvx; 
  incr[1] = (-1.936491673103709*fl[11]*wx)-1.5*fl[1]*wx-0.8660254037844386*fl[0]*wx-1.118033988749895*fl[19]*dvx-0.8660254037844386*fl[5]*dvx-0.5*fl[2]*dvx; 
  incr[2] = 1.118033988749895*fl[19]*wx+0.8660254037844386*fl[5]*wx+0.5*fl[2]*wx+0.4472135954999579*fl[20]*dvx+0.2581988897471612*fl[12]*dvx+0.6454972243679029*fl[11]*dvx+0.5*fl[1]*dvx+0.2886751345948129*fl[0]*dvx; 
  incr[3] = 1.118033988749895*fl[21]*wx+0.8660254037844386*fl[6]*wx+0.5*fl[3]*wx+0.6454972243679029*fl[32]*dvx+0.5*fl[15]*dvx+0.2886751345948129*fl[7]*dvx; 
  incr[4] = 1.118033988749895*fl[25]*wx+0.8660254037844386*fl[8]*wx+0.5*fl[4]*wx+0.6454972243679029*fl[35]*dvx+0.5*fl[16]*dvx+0.2886751345948129*fl[9]*dvx; 
  incr[5] = (-1.936491673103709*fl[19]*wx)-1.5*fl[5]*wx-0.8660254037844386*fl[2]*wx-0.7745966692414833*fl[20]*dvx-0.4472135954999579*fl[12]*dvx-1.118033988749895*fl[11]*dvx-0.8660254037844386*fl[1]*dvx-0.5*fl[0]*dvx; 
  incr[6] = (-1.936491673103709*fl[21]*wx)-1.5*fl[6]*wx-0.8660254037844386*fl[3]*wx-1.118033988749895*fl[32]*dvx-0.8660254037844386*fl[15]*dvx-0.5*fl[7]*dvx; 
  incr[7] = 1.118033988749895*fl[32]*wx+0.8660254037844386*fl[15]*wx+0.5*fl[7]*wx+0.4472135954999579*fl[33]*dvx+0.2581988897471612*fl[22]*dvx+0.6454972243679029*fl[21]*dvx+0.5*fl[6]*dvx+0.2886751345948129*fl[3]*dvx; 
  incr[8] = (-1.936491673103709*fl[25]*wx)-1.5*fl[8]*wx-0.8660254037844386*fl[4]*wx-1.118033988749895*fl[35]*dvx-0.8660254037844386*fl[16]*dvx-0.5*fl[9]*dvx; 
  incr[9] = 1.118033988749895*fl[35]*wx+0.8660254037844386*fl[16]*wx+0.5*fl[9]*wx+0.4472135954999579*fl[36]*dvx+0.2581988897471612*fl[26]*dvx+0.6454972243679029*fl[25]*dvx+0.5*fl[8]*dvx+0.2886751345948129*fl[4]*dvx; 
  incr[10] = 1.118033988749895*fl[37]*wx+0.8660254037844386*fl[17]*wx+0.5*fl[10]*wx+0.6454972243679029*fl[44]*dvx+0.5*fl[31]*dvx+0.2886751345948129*fl[18]*dvx; 
  incr[11] = 2.5*fl[11]*wx+1.936491673103709*fl[1]*wx+1.118033988749895*fl[0]*wx+1.443375672974065*fl[19]*dvx+1.118033988749895*fl[5]*dvx+0.6454972243679029*fl[2]*dvx; 
  incr[12] = 0.8660254037844386*fl[20]*wx+0.5*fl[12]*wx+0.5773502691896258*fl[19]*dvx+0.4472135954999579*fl[5]*dvx+0.2581988897471612*fl[2]*dvx; 
  incr[13] = 0.8660254037844386*fl[23]*wx+0.5*fl[13]*wx+0.5*fl[34]*dvx+0.2886751345948129*fl[24]*dvx; 
  incr[14] = 0.8660254037844386*fl[28]*wx+0.5*fl[14]*wx+0.5*fl[41]*dvx+0.2886751345948129*fl[29]*dvx; 
  incr[15] = (-1.936491673103709*fl[32]*wx)-1.5*fl[15]*wx-0.8660254037844386*fl[7]*wx-0.7745966692414833*fl[33]*dvx-0.4472135954999579*fl[22]*dvx-1.118033988749895*fl[21]*dvx-0.8660254037844386*fl[6]*dvx-0.5*fl[3]*dvx; 
  incr[16] = (-1.936491673103709*fl[35]*wx)-1.5*fl[16]*wx-0.8660254037844386*fl[9]*wx-0.7745966692414833*fl[36]*dvx-0.4472135954999579*fl[26]*dvx-1.118033988749895*fl[25]*dvx-0.8660254037844386*fl[8]*dvx-0.5*fl[4]*dvx; 
  incr[17] = (-1.936491673103709*fl[37]*wx)-1.5*fl[17]*wx-0.8660254037844386*fl[10]*wx-1.118033988749895*fl[44]*dvx-0.8660254037844386*fl[31]*dvx-0.5*fl[18]*dvx; 
  incr[18] = 1.118033988749895*fl[44]*wx+0.8660254037844386*fl[31]*wx+0.5*fl[18]*wx+0.4472135954999579*fl[45]*dvx+0.2581988897471612*fl[38]*dvx+0.6454972243679029*fl[37]*dvx+0.5*fl[17]*dvx+0.2886751345948129*fl[10]*dvx; 
  incr[19] = 2.5*fl[19]*wx+1.936491673103709*fl[5]*wx+1.118033988749895*fl[2]*wx+fl[20]*dvx+0.5773502691896258*fl[12]*dvx+1.443375672974065*fl[11]*dvx+1.118033988749895*fl[1]*dvx+0.6454972243679029*fl[0]*dvx; 
  incr[20] = (-1.5*fl[20]*wx)-0.8660254037844386*fl[12]*wx-1.0*fl[19]*dvx-0.7745966692414833*fl[5]*dvx-0.4472135954999579*fl[2]*dvx; 
  incr[21] = 2.5*fl[21]*wx+1.936491673103709*fl[6]*wx+1.118033988749895*fl[3]*wx+1.443375672974065*fl[32]*dvx+1.118033988749895*fl[15]*dvx+0.6454972243679029*fl[7]*dvx; 
  incr[22] = 0.8660254037844386*fl[33]*wx+0.5*fl[22]*wx+0.5773502691896258*fl[32]*dvx+0.4472135954999579*fl[15]*dvx+0.2581988897471612*fl[7]*dvx; 
  incr[23] = (-1.5*fl[23]*wx)-0.8660254037844386*fl[13]*wx-0.8660254037844386*fl[34]*dvx-0.5*fl[24]*dvx; 
  incr[24] = 0.8660254037844386*fl[34]*wx+0.5*fl[24]*wx+0.5*fl[23]*dvx+0.2886751345948129*fl[13]*dvx; 
  incr[25] = 2.5*fl[25]*wx+1.936491673103709*fl[8]*wx+1.118033988749895*fl[4]*wx+1.443375672974065*fl[35]*dvx+1.118033988749895*fl[16]*dvx+0.6454972243679029*fl[9]*dvx; 
  incr[26] = 0.8660254037844386*fl[36]*wx+0.5*fl[26]*wx+0.5773502691896258*fl[35]*dvx+0.4472135954999579*fl[16]*dvx+0.2581988897471612*fl[9]*dvx; 
  incr[27] = 0.8660254037844386*fl[39]*wx+0.5*fl[27]*wx+0.5*fl[46]*dvx+0.2886751345948129*fl[40]*dvx; 
  incr[28] = (-1.5*fl[28]*wx)-0.8660254037844386*fl[14]*wx-0.8660254037844386*fl[41]*dvx-0.5*fl[29]*dvx; 
  incr[29] = 0.8660254037844386*fl[41]*wx+0.5*fl[29]*wx+0.5*fl[28]*dvx+0.2886751345948129*fl[14]*dvx; 
  incr[30] = 0.8660254037844386*fl[42]*wx+0.5*fl[30]*wx+0.5*fl[47]*dvx+0.2886751345948129*fl[43]*dvx; 
  incr[31] = (-1.936491673103709*fl[44]*wx)-1.5*fl[31]*wx-0.8660254037844386*fl[18]*wx-0.7745966692414833*fl[45]*dvx-0.4472135954999579*fl[38]*dvx-1.118033988749895*fl[37]*dvx-0.8660254037844386*fl[17]*dvx-0.5*fl[10]*dvx; 
  incr[32] = 2.5*fl[32]*wx+1.936491673103709*fl[15]*wx+1.118033988749895*fl[7]*wx+fl[33]*dvx+0.5773502691896258*fl[22]*dvx+1.443375672974065*fl[21]*dvx+1.118033988749895*fl[6]*dvx+0.6454972243679029*fl[3]*dvx; 
  incr[33] = (-1.5*fl[33]*wx)-0.8660254037844386*fl[22]*wx-1.0*fl[32]*dvx-0.7745966692414833*fl[15]*dvx-0.4472135954999579*fl[7]*dvx; 
  incr[34] = (-1.5*fl[34]*wx)-0.8660254037844386*fl[24]*wx-0.8660254037844386*fl[23]*dvx-0.5*fl[13]*dvx; 
  incr[35] = 2.5*fl[35]*wx+1.936491673103709*fl[16]*wx+1.118033988749895*fl[9]*wx+fl[36]*dvx+0.5773502691896258*fl[26]*dvx+1.443375672974065*fl[25]*dvx+1.118033988749895*fl[8]*dvx+0.6454972243679029*fl[4]*dvx; 
  incr[36] = (-1.5*fl[36]*wx)-0.8660254037844386*fl[26]*wx-1.0*fl[35]*dvx-0.7745966692414833*fl[16]*dvx-0.4472135954999579*fl[9]*dvx; 
  incr[37] = 2.5*fl[37]*wx+1.936491673103709*fl[17]*wx+1.118033988749895*fl[10]*wx+1.443375672974065*fl[44]*dvx+1.118033988749895*fl[31]*dvx+0.6454972243679029*fl[18]*dvx; 
  incr[38] = 0.8660254037844386*fl[45]*wx+0.5*fl[38]*wx+0.5773502691896258*fl[44]*dvx+0.4472135954999579*fl[31]*dvx+0.2581988897471612*fl[18]*dvx; 
  incr[39] = (-1.5*fl[39]*wx)-0.8660254037844386*fl[27]*wx-0.8660254037844386*fl[46]*dvx-0.5*fl[40]*dvx; 
  incr[40] = 0.8660254037844386*fl[46]*wx+0.5*fl[40]*wx+0.5*fl[39]*dvx+0.2886751345948129*fl[27]*dvx; 
  incr[41] = (-1.5*fl[41]*wx)-0.8660254037844386*fl[29]*wx-0.8660254037844386*fl[28]*dvx-0.5*fl[14]*dvx; 
  incr[42] = (-1.5*fl[42]*wx)-0.8660254037844386*fl[30]*wx-0.8660254037844386*fl[47]*dvx-0.5*fl[43]*dvx; 
  incr[43] = 0.8660254037844386*fl[47]*wx+0.5*fl[43]*wx+0.5*fl[42]*dvx+0.2886751345948129*fl[30]*dvx; 
  incr[44] = 2.5*fl[44]*wx+1.936491673103709*fl[31]*wx+1.118033988749895*fl[18]*wx+fl[45]*dvx+0.5773502691896258*fl[38]*dvx+1.443375672974065*fl[37]*dvx+1.118033988749895*fl[17]*dvx+0.6454972243679029*fl[10]*dvx; 
  incr[45] = (-1.5*fl[45]*wx)-0.8660254037844386*fl[38]*wx-1.0*fl[44]*dvx-0.7745966692414833*fl[31]*dvx-0.4472135954999579*fl[18]*dvx; 
  incr[46] = (-1.5*fl[46]*wx)-0.8660254037844386*fl[40]*wx-0.8660254037844386*fl[39]*dvx-0.5*fl[27]*dvx; 
  incr[47] = (-1.5*fl[47]*wx)-0.8660254037844386*fl[43]*wx-0.8660254037844386*fl[42]*dvx-0.5*fl[30]*dvx; 

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
  outr[28] += incr[28]; 
  outr[29] += incr[29]; 
  outr[30] += incr[30]; 
  outr[31] += incr[31]; 
  outr[32] += incr[32]; 
  outr[33] += incr[33]; 
  outr[34] += incr[34]; 
  outr[35] += incr[35]; 
  outr[36] += incr[36]; 
  outr[37] += incr[37]; 
  outr[38] += incr[38]; 
  outr[39] += incr[39]; 
  outr[40] += incr[40]; 
  outr[41] += incr[41]; 
  outr[42] += incr[42]; 
  outr[43] += incr[43]; 
  outr[44] += incr[44]; 
  outr[45] += incr[45]; 
  outr[46] += incr[46]; 
  outr[47] += incr[47]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += -1.0*incr[19]; 
  outl[20] += incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += -1.0*incr[26]; 
  outl[27] += -1.0*incr[27]; 
  outl[28] += incr[28]; 
  outl[29] += -1.0*incr[29]; 
  outl[30] += -1.0*incr[30]; 
  outl[31] += incr[31]; 
  outl[32] += -1.0*incr[32]; 
  outl[33] += incr[33]; 
  outl[34] += incr[34]; 
  outl[35] += -1.0*incr[35]; 
  outl[36] += incr[36]; 
  outl[37] += -1.0*incr[37]; 
  outl[38] += -1.0*incr[38]; 
  outl[39] += incr[39]; 
  outl[40] += -1.0*incr[40]; 
  outl[41] += incr[41]; 
  outl[42] += incr[42]; 
  outl[43] += -1.0*incr[43]; 
  outl[44] += -1.0*incr[44]; 
  outl[45] += incr[45]; 
  outl[46] += incr[46]; 
  outl[47] += incr[47]; 
  } else { 
  incr[0] = 1.118033988749895*fr[11]*wx-0.8660254037844386*fr[1]*wx+0.5*fr[0]*wx+0.6454972243679029*fr[19]*dvx-0.5*fr[5]*dvx+0.2886751345948129*fr[2]*dvx; 
  incr[1] = (-1.936491673103709*fr[11]*wx)+1.5*fr[1]*wx-0.8660254037844386*fr[0]*wx-1.118033988749895*fr[19]*dvx+0.8660254037844386*fr[5]*dvx-0.5*fr[2]*dvx; 
  incr[2] = 1.118033988749895*fr[19]*wx-0.8660254037844386*fr[5]*wx+0.5*fr[2]*wx-0.4472135954999579*fr[20]*dvx+0.2581988897471612*fr[12]*dvx+0.6454972243679029*fr[11]*dvx-0.5*fr[1]*dvx+0.2886751345948129*fr[0]*dvx; 
  incr[3] = 1.118033988749895*fr[21]*wx-0.8660254037844386*fr[6]*wx+0.5*fr[3]*wx+0.6454972243679029*fr[32]*dvx-0.5*fr[15]*dvx+0.2886751345948129*fr[7]*dvx; 
  incr[4] = 1.118033988749895*fr[25]*wx-0.8660254037844386*fr[8]*wx+0.5*fr[4]*wx+0.6454972243679029*fr[35]*dvx-0.5*fr[16]*dvx+0.2886751345948129*fr[9]*dvx; 
  incr[5] = (-1.936491673103709*fr[19]*wx)+1.5*fr[5]*wx-0.8660254037844386*fr[2]*wx+0.7745966692414833*fr[20]*dvx-0.4472135954999579*fr[12]*dvx-1.118033988749895*fr[11]*dvx+0.8660254037844386*fr[1]*dvx-0.5*fr[0]*dvx; 
  incr[6] = (-1.936491673103709*fr[21]*wx)+1.5*fr[6]*wx-0.8660254037844386*fr[3]*wx-1.118033988749895*fr[32]*dvx+0.8660254037844386*fr[15]*dvx-0.5*fr[7]*dvx; 
  incr[7] = 1.118033988749895*fr[32]*wx-0.8660254037844386*fr[15]*wx+0.5*fr[7]*wx-0.4472135954999579*fr[33]*dvx+0.2581988897471612*fr[22]*dvx+0.6454972243679029*fr[21]*dvx-0.5*fr[6]*dvx+0.2886751345948129*fr[3]*dvx; 
  incr[8] = (-1.936491673103709*fr[25]*wx)+1.5*fr[8]*wx-0.8660254037844386*fr[4]*wx-1.118033988749895*fr[35]*dvx+0.8660254037844386*fr[16]*dvx-0.5*fr[9]*dvx; 
  incr[9] = 1.118033988749895*fr[35]*wx-0.8660254037844386*fr[16]*wx+0.5*fr[9]*wx-0.4472135954999579*fr[36]*dvx+0.2581988897471612*fr[26]*dvx+0.6454972243679029*fr[25]*dvx-0.5*fr[8]*dvx+0.2886751345948129*fr[4]*dvx; 
  incr[10] = 1.118033988749895*fr[37]*wx-0.8660254037844386*fr[17]*wx+0.5*fr[10]*wx+0.6454972243679029*fr[44]*dvx-0.5*fr[31]*dvx+0.2886751345948129*fr[18]*dvx; 
  incr[11] = 2.5*fr[11]*wx-1.936491673103709*fr[1]*wx+1.118033988749895*fr[0]*wx+1.443375672974065*fr[19]*dvx-1.118033988749895*fr[5]*dvx+0.6454972243679029*fr[2]*dvx; 
  incr[12] = (-0.8660254037844386*fr[20]*wx)+0.5*fr[12]*wx+0.5773502691896258*fr[19]*dvx-0.4472135954999579*fr[5]*dvx+0.2581988897471612*fr[2]*dvx; 
  incr[13] = (-0.8660254037844386*fr[23]*wx)+0.5*fr[13]*wx-0.5*fr[34]*dvx+0.2886751345948129*fr[24]*dvx; 
  incr[14] = (-0.8660254037844386*fr[28]*wx)+0.5*fr[14]*wx-0.5*fr[41]*dvx+0.2886751345948129*fr[29]*dvx; 
  incr[15] = (-1.936491673103709*fr[32]*wx)+1.5*fr[15]*wx-0.8660254037844386*fr[7]*wx+0.7745966692414833*fr[33]*dvx-0.4472135954999579*fr[22]*dvx-1.118033988749895*fr[21]*dvx+0.8660254037844386*fr[6]*dvx-0.5*fr[3]*dvx; 
  incr[16] = (-1.936491673103709*fr[35]*wx)+1.5*fr[16]*wx-0.8660254037844386*fr[9]*wx+0.7745966692414833*fr[36]*dvx-0.4472135954999579*fr[26]*dvx-1.118033988749895*fr[25]*dvx+0.8660254037844386*fr[8]*dvx-0.5*fr[4]*dvx; 
  incr[17] = (-1.936491673103709*fr[37]*wx)+1.5*fr[17]*wx-0.8660254037844386*fr[10]*wx-1.118033988749895*fr[44]*dvx+0.8660254037844386*fr[31]*dvx-0.5*fr[18]*dvx; 
  incr[18] = 1.118033988749895*fr[44]*wx-0.8660254037844386*fr[31]*wx+0.5*fr[18]*wx-0.4472135954999579*fr[45]*dvx+0.2581988897471612*fr[38]*dvx+0.6454972243679029*fr[37]*dvx-0.5*fr[17]*dvx+0.2886751345948129*fr[10]*dvx; 
  incr[19] = 2.5*fr[19]*wx-1.936491673103709*fr[5]*wx+1.118033988749895*fr[2]*wx-1.0*fr[20]*dvx+0.5773502691896258*fr[12]*dvx+1.443375672974065*fr[11]*dvx-1.118033988749895*fr[1]*dvx+0.6454972243679029*fr[0]*dvx; 
  incr[20] = 1.5*fr[20]*wx-0.8660254037844386*fr[12]*wx-1.0*fr[19]*dvx+0.7745966692414833*fr[5]*dvx-0.4472135954999579*fr[2]*dvx; 
  incr[21] = 2.5*fr[21]*wx-1.936491673103709*fr[6]*wx+1.118033988749895*fr[3]*wx+1.443375672974065*fr[32]*dvx-1.118033988749895*fr[15]*dvx+0.6454972243679029*fr[7]*dvx; 
  incr[22] = (-0.8660254037844386*fr[33]*wx)+0.5*fr[22]*wx+0.5773502691896258*fr[32]*dvx-0.4472135954999579*fr[15]*dvx+0.2581988897471612*fr[7]*dvx; 
  incr[23] = 1.5*fr[23]*wx-0.8660254037844386*fr[13]*wx+0.8660254037844386*fr[34]*dvx-0.5*fr[24]*dvx; 
  incr[24] = (-0.8660254037844386*fr[34]*wx)+0.5*fr[24]*wx-0.5*fr[23]*dvx+0.2886751345948129*fr[13]*dvx; 
  incr[25] = 2.5*fr[25]*wx-1.936491673103709*fr[8]*wx+1.118033988749895*fr[4]*wx+1.443375672974065*fr[35]*dvx-1.118033988749895*fr[16]*dvx+0.6454972243679029*fr[9]*dvx; 
  incr[26] = (-0.8660254037844386*fr[36]*wx)+0.5*fr[26]*wx+0.5773502691896258*fr[35]*dvx-0.4472135954999579*fr[16]*dvx+0.2581988897471612*fr[9]*dvx; 
  incr[27] = (-0.8660254037844386*fr[39]*wx)+0.5*fr[27]*wx-0.5*fr[46]*dvx+0.2886751345948129*fr[40]*dvx; 
  incr[28] = 1.5*fr[28]*wx-0.8660254037844386*fr[14]*wx+0.8660254037844386*fr[41]*dvx-0.5*fr[29]*dvx; 
  incr[29] = (-0.8660254037844386*fr[41]*wx)+0.5*fr[29]*wx-0.5*fr[28]*dvx+0.2886751345948129*fr[14]*dvx; 
  incr[30] = (-0.8660254037844386*fr[42]*wx)+0.5*fr[30]*wx-0.5*fr[47]*dvx+0.2886751345948129*fr[43]*dvx; 
  incr[31] = (-1.936491673103709*fr[44]*wx)+1.5*fr[31]*wx-0.8660254037844386*fr[18]*wx+0.7745966692414833*fr[45]*dvx-0.4472135954999579*fr[38]*dvx-1.118033988749895*fr[37]*dvx+0.8660254037844386*fr[17]*dvx-0.5*fr[10]*dvx; 
  incr[32] = 2.5*fr[32]*wx-1.936491673103709*fr[15]*wx+1.118033988749895*fr[7]*wx-1.0*fr[33]*dvx+0.5773502691896258*fr[22]*dvx+1.443375672974065*fr[21]*dvx-1.118033988749895*fr[6]*dvx+0.6454972243679029*fr[3]*dvx; 
  incr[33] = 1.5*fr[33]*wx-0.8660254037844386*fr[22]*wx-1.0*fr[32]*dvx+0.7745966692414833*fr[15]*dvx-0.4472135954999579*fr[7]*dvx; 
  incr[34] = 1.5*fr[34]*wx-0.8660254037844386*fr[24]*wx+0.8660254037844386*fr[23]*dvx-0.5*fr[13]*dvx; 
  incr[35] = 2.5*fr[35]*wx-1.936491673103709*fr[16]*wx+1.118033988749895*fr[9]*wx-1.0*fr[36]*dvx+0.5773502691896258*fr[26]*dvx+1.443375672974065*fr[25]*dvx-1.118033988749895*fr[8]*dvx+0.6454972243679029*fr[4]*dvx; 
  incr[36] = 1.5*fr[36]*wx-0.8660254037844386*fr[26]*wx-1.0*fr[35]*dvx+0.7745966692414833*fr[16]*dvx-0.4472135954999579*fr[9]*dvx; 
  incr[37] = 2.5*fr[37]*wx-1.936491673103709*fr[17]*wx+1.118033988749895*fr[10]*wx+1.443375672974065*fr[44]*dvx-1.118033988749895*fr[31]*dvx+0.6454972243679029*fr[18]*dvx; 
  incr[38] = (-0.8660254037844386*fr[45]*wx)+0.5*fr[38]*wx+0.5773502691896258*fr[44]*dvx-0.4472135954999579*fr[31]*dvx+0.2581988897471612*fr[18]*dvx; 
  incr[39] = 1.5*fr[39]*wx-0.8660254037844386*fr[27]*wx+0.8660254037844386*fr[46]*dvx-0.5*fr[40]*dvx; 
  incr[40] = (-0.8660254037844386*fr[46]*wx)+0.5*fr[40]*wx-0.5*fr[39]*dvx+0.2886751345948129*fr[27]*dvx; 
  incr[41] = 1.5*fr[41]*wx-0.8660254037844386*fr[29]*wx+0.8660254037844386*fr[28]*dvx-0.5*fr[14]*dvx; 
  incr[42] = 1.5*fr[42]*wx-0.8660254037844386*fr[30]*wx+0.8660254037844386*fr[47]*dvx-0.5*fr[43]*dvx; 
  incr[43] = (-0.8660254037844386*fr[47]*wx)+0.5*fr[43]*wx-0.5*fr[42]*dvx+0.2886751345948129*fr[30]*dvx; 
  incr[44] = 2.5*fr[44]*wx-1.936491673103709*fr[31]*wx+1.118033988749895*fr[18]*wx-1.0*fr[45]*dvx+0.5773502691896258*fr[38]*dvx+1.443375672974065*fr[37]*dvx-1.118033988749895*fr[17]*dvx+0.6454972243679029*fr[10]*dvx; 
  incr[45] = 1.5*fr[45]*wx-0.8660254037844386*fr[38]*wx-1.0*fr[44]*dvx+0.7745966692414833*fr[31]*dvx-0.4472135954999579*fr[18]*dvx; 
  incr[46] = 1.5*fr[46]*wx-0.8660254037844386*fr[40]*wx+0.8660254037844386*fr[39]*dvx-0.5*fr[27]*dvx; 
  incr[47] = 1.5*fr[47]*wx-0.8660254037844386*fr[43]*wx+0.8660254037844386*fr[42]*dvx-0.5*fr[30]*dvx; 

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
  outr[28] += incr[28]; 
  outr[29] += incr[29]; 
  outr[30] += incr[30]; 
  outr[31] += incr[31]; 
  outr[32] += incr[32]; 
  outr[33] += incr[33]; 
  outr[34] += incr[34]; 
  outr[35] += incr[35]; 
  outr[36] += incr[36]; 
  outr[37] += incr[37]; 
  outr[38] += incr[38]; 
  outr[39] += incr[39]; 
  outr[40] += incr[40]; 
  outr[41] += incr[41]; 
  outr[42] += incr[42]; 
  outr[43] += incr[43]; 
  outr[44] += incr[44]; 
  outr[45] += incr[45]; 
  outr[46] += incr[46]; 
  outr[47] += incr[47]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += -1.0*incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  outl[16] += incr[16]; 
  outl[17] += incr[17]; 
  outl[18] += -1.0*incr[18]; 
  outl[19] += -1.0*incr[19]; 
  outl[20] += incr[20]; 
  outl[21] += -1.0*incr[21]; 
  outl[22] += -1.0*incr[22]; 
  outl[23] += incr[23]; 
  outl[24] += -1.0*incr[24]; 
  outl[25] += -1.0*incr[25]; 
  outl[26] += -1.0*incr[26]; 
  outl[27] += -1.0*incr[27]; 
  outl[28] += incr[28]; 
  outl[29] += -1.0*incr[29]; 
  outl[30] += -1.0*incr[30]; 
  outl[31] += incr[31]; 
  outl[32] += -1.0*incr[32]; 
  outl[33] += incr[33]; 
  outl[34] += incr[34]; 
  outl[35] += -1.0*incr[35]; 
  outl[36] += incr[36]; 
  outl[37] += -1.0*incr[37]; 
  outl[38] += -1.0*incr[38]; 
  outl[39] += incr[39]; 
  outl[40] += -1.0*incr[40]; 
  outl[41] += incr[41]; 
  outl[42] += incr[42]; 
  outl[43] += -1.0*incr[43]; 
  outl[44] += -1.0*incr[44]; 
  outl[45] += incr[45]; 
  outl[46] += incr[46]; 
  outl[47] += incr[47]; 
  } 
} 
