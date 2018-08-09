#include <VlasovModDecl.h> 
void VlasovSurfStream2x2vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[2]; 
  double wxl = wl[2]; 

  double dvxr = dxvr[2]; 
  double wxr = wr[2]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[5]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[1]+fl[0])*wxl+0.2886751345948129*fl[3]*dvxl; 
  incr[1] = ((-3.0*fl[1])-1.732050807568877*fl[0])*wxl-0.5*fl[3]*dvxl; 
  incr[2] = fl[2]*wxl; 
  incr[3] = fl[3]*wxl+(0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[4] = fl[4]*wxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  } else { 
  incr[0] = (fr[0]-1.732050807568877*fr[1])*wxr+0.2886751345948129*fr[3]*dvxr; 
  incr[1] = (3.0*fr[1]-1.732050807568877*fr[0])*wxr-0.5*fr[3]*dvxr; 
  incr[2] = fr[2]*wxr; 
  incr[3] = fr[3]*wxr+(0.2886751345948129*fr[0]-0.5*fr[1])*dvxr; 
  incr[4] = fr[4]*wxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  } 
} 
void VlasovSurfStream2x2vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[2]; 
  double wxl = wl[2]; 

  double dvxr = dxvr[2]; 
  double wxr = wr[2]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[15]; 

  if (wxr>0) { 
  incr[0] = (2.23606797749979*fl[11]+1.732050807568877*fl[1]+fl[0])*wxl+(0.5*fl[6]+0.2886751345948129*fl[3])*dvxl; 
  incr[1] = ((-3.872983346207417*fl[11])-3.0*fl[1]-1.732050807568877*fl[0])*wxl+((-0.8660254037844386*fl[6])-0.5*fl[3])*dvxl; 
  incr[2] = (1.732050807568877*fl[5]+fl[2])*wxl+0.2886751345948129*fl[7]*dvxl; 
  incr[3] = (1.732050807568877*fl[6]+fl[3])*wxl+(0.2581988897471612*fl[13]+0.6454972243679029*fl[11]+0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[4] = (1.732050807568877*fl[8]+fl[4])*wxl+0.2886751345948129*fl[10]*dvxl; 
  incr[5] = ((-3.0*fl[5])-1.732050807568877*fl[2])*wxl-0.5*fl[7]*dvxl; 
  incr[6] = ((-3.0*fl[6])-1.732050807568877*fl[3])*wxl+((-0.4472135954999579*fl[13])-1.118033988749895*fl[11]-0.8660254037844386*fl[1]-0.5*fl[0])*dvxl; 
  incr[7] = fl[7]*wxl+(0.5*fl[5]+0.2886751345948129*fl[2])*dvxl; 
  incr[8] = ((-3.0*fl[8])-1.732050807568877*fl[4])*wxl-0.5*fl[10]*dvxl; 
  incr[9] = fl[9]*wxl; 
  incr[10] = fl[10]*wxl+(0.5*fl[8]+0.2886751345948129*fl[4])*dvxl; 
  incr[11] = (5.0*fl[11]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*wxl+(1.118033988749895*fl[6]+0.6454972243679029*fl[3])*dvxl; 
  incr[12] = fl[12]*wxl; 
  incr[13] = fl[13]*wxl+(0.4472135954999579*fl[6]+0.2581988897471612*fl[3])*dvxl; 
  incr[14] = fl[14]*wxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  } else { 
  incr[0] = (2.23606797749979*fr[11]-1.732050807568877*fr[1]+fr[0])*wxr+(0.2886751345948129*fr[3]-0.5*fr[6])*dvxr; 
  incr[1] = ((-3.872983346207417*fr[11])+3.0*fr[1]-1.732050807568877*fr[0])*wxr+(0.8660254037844386*fr[6]-0.5*fr[3])*dvxr; 
  incr[2] = (fr[2]-1.732050807568877*fr[5])*wxr+0.2886751345948129*fr[7]*dvxr; 
  incr[3] = (fr[3]-1.732050807568877*fr[6])*wxr+(0.2581988897471612*fr[13]+0.6454972243679029*fr[11]-0.5*fr[1]+0.2886751345948129*fr[0])*dvxr; 
  incr[4] = (fr[4]-1.732050807568877*fr[8])*wxr+0.2886751345948129*fr[10]*dvxr; 
  incr[5] = (3.0*fr[5]-1.732050807568877*fr[2])*wxr-0.5*fr[7]*dvxr; 
  incr[6] = (3.0*fr[6]-1.732050807568877*fr[3])*wxr+((-0.4472135954999579*fr[13])-1.118033988749895*fr[11]+0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[7] = fr[7]*wxr+(0.2886751345948129*fr[2]-0.5*fr[5])*dvxr; 
  incr[8] = (3.0*fr[8]-1.732050807568877*fr[4])*wxr-0.5*fr[10]*dvxr; 
  incr[9] = fr[9]*wxr; 
  incr[10] = fr[10]*wxr+(0.2886751345948129*fr[4]-0.5*fr[8])*dvxr; 
  incr[11] = (5.0*fr[11]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*wxr+(0.6454972243679029*fr[3]-1.118033988749895*fr[6])*dvxr; 
  incr[12] = fr[12]*wxr; 
  incr[13] = fr[13]*wxr+(0.2581988897471612*fr[3]-0.4472135954999579*fr[6])*dvxr; 
  incr[14] = fr[14]*wxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  } 
} 
void VlasovSurfStream2x2vMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[2]; 
  double wxl = wl[2]; 

  double dvxr = dxvr[2]; 
  double wxr = wr[2]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[35]; 

  if (wxr>0) { 
  incr[0] = (2.645751311064591*fl[31]+2.23606797749979*fl[11]+1.732050807568877*fl[1]+fl[0])*wxl+(0.6454972243679028*fl[21]+0.5*fl[6]+0.2886751345948129*fl[3])*dvxl; 
  incr[1] = ((-4.58257569495584*fl[31])-3.872983346207417*fl[11]-3.0*fl[1]-1.732050807568877*fl[0])*wxl+((-1.118033988749895*fl[21])-0.8660254037844386*fl[6]-0.5*fl[3])*dvxl; 
  incr[2] = (2.23606797749979*fl[19]+1.732050807568877*fl[5]+fl[2])*wxl+(0.5*fl[15]+0.2886751345948129*fl[7])*dvxl; 
  incr[3] = (2.23606797749979*fl[21]+1.732050807568877*fl[6]+fl[3])*wxl+(0.7637626158259735*fl[31]+0.447213595499958*fl[23]+0.2581988897471612*fl[13]+0.6454972243679029*fl[11]+0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[4] = (2.23606797749979*fl[25]+1.732050807568877*fl[8]+fl[4])*wxl+(0.5*fl[17]+0.2886751345948129*fl[10])*dvxl; 
  incr[5] = ((-3.872983346207417*fl[19])-3.0*fl[5]-1.732050807568877*fl[2])*wxl+((-0.8660254037844386*fl[15])-0.5*fl[7])*dvxl; 
  incr[6] = ((-3.872983346207417*fl[21])-3.0*fl[6]-1.732050807568877*fl[3])*wxl+((-1.322875655532295*fl[31])-0.7745966692414834*fl[23]-0.4472135954999579*fl[13]-1.118033988749895*fl[11]-0.8660254037844386*fl[1]-0.5*fl[0])*dvxl; 
  incr[7] = (1.732050807568877*fl[15]+fl[7])*wxl+(0.2581988897471611*fl[24]+0.6454972243679028*fl[19]+0.5*fl[5]+0.2886751345948129*fl[2])*dvxl; 
  incr[8] = ((-3.872983346207417*fl[25])-3.0*fl[8]-1.732050807568877*fl[4])*wxl+((-0.8660254037844386*fl[17])-0.5*fl[10])*dvxl; 
  incr[9] = (1.732050807568877*fl[16]+fl[9])*wxl+0.2886751345948129*fl[18]*dvxl; 
  incr[10] = (1.732050807568877*fl[17]+fl[10])*wxl+(0.2581988897471611*fl[27]+0.6454972243679028*fl[25]+0.5*fl[8]+0.2886751345948129*fl[4])*dvxl; 
  incr[11] = (5.916079783099617*fl[31]+5.0*fl[11]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*wxl+(1.443375672974065*fl[21]+1.118033988749895*fl[6]+0.6454972243679029*fl[3])*dvxl; 
  incr[12] = (1.732050807568877*fl[20]+fl[12])*wxl+0.2886751345948129*fl[22]*dvxl; 
  incr[13] = (1.732050807568877*fl[23]+fl[13])*wxl+(0.253546276418555*fl[33]+0.5773502691896257*fl[21]+0.4472135954999579*fl[6]+0.2581988897471612*fl[3])*dvxl; 
  incr[14] = (1.732050807568877*fl[28]+fl[14])*wxl+0.2886751345948129*fl[30]*dvxl; 
  incr[15] = ((-3.0*fl[15])-1.732050807568877*fl[7])*wxl+((-0.447213595499958*fl[24])-1.118033988749895*fl[19]-0.8660254037844386*fl[5]-0.5*fl[2])*dvxl; 
  incr[16] = ((-3.0*fl[16])-1.732050807568877*fl[9])*wxl-0.5*fl[18]*dvxl; 
  incr[17] = ((-3.0*fl[17])-1.732050807568877*fl[10])*wxl+((-0.447213595499958*fl[27])-1.118033988749895*fl[25]-0.8660254037844386*fl[8]-0.5*fl[4])*dvxl; 
  incr[18] = fl[18]*wxl+(0.5*fl[16]+0.2886751345948129*fl[9])*dvxl; 
  incr[19] = (5.0*fl[19]+3.872983346207417*fl[5]+2.23606797749979*fl[2])*wxl+(1.118033988749895*fl[15]+0.6454972243679028*fl[7])*dvxl; 
  incr[20] = ((-3.0*fl[20])-1.732050807568877*fl[12])*wxl-0.5*fl[22]*dvxl; 
  incr[21] = (5.0*fl[21]+3.872983346207417*fl[6]+2.23606797749979*fl[3])*wxl+(1.707825127659933*fl[31]+fl[23]+0.5773502691896257*fl[13]+1.443375672974065*fl[11]+1.118033988749895*fl[1]+0.6454972243679028*fl[0])*dvxl; 
  incr[22] = fl[22]*wxl+(0.5*fl[20]+0.2886751345948129*fl[12])*dvxl; 
  incr[23] = ((-3.0*fl[23])-1.732050807568877*fl[13])*wxl+((-0.4391550328268399*fl[33])-1.0*fl[21]-0.7745966692414834*fl[6]-0.447213595499958*fl[3])*dvxl; 
  incr[24] = fl[24]*wxl+(0.447213595499958*fl[15]+0.2581988897471611*fl[7])*dvxl; 
  incr[25] = (5.0*fl[25]+3.872983346207417*fl[8]+2.23606797749979*fl[4])*wxl+(1.118033988749895*fl[17]+0.6454972243679028*fl[10])*dvxl; 
  incr[26] = fl[26]*wxl; 
  incr[27] = fl[27]*wxl+(0.447213595499958*fl[17]+0.2581988897471611*fl[10])*dvxl; 
  incr[28] = ((-3.0*fl[28])-1.732050807568877*fl[14])*wxl-0.5*fl[30]*dvxl; 
  incr[29] = fl[29]*wxl; 
  incr[30] = fl[30]*wxl+(0.5*fl[28]+0.2886751345948129*fl[14])*dvxl; 
  incr[31] = ((-7.0*fl[31])-5.916079783099617*fl[11]-4.58257569495584*fl[1]-2.645751311064591*fl[0])*wxl+((-1.707825127659933*fl[21])-1.322875655532295*fl[6]-0.7637626158259735*fl[3])*dvxl; 
  incr[32] = fl[32]*wxl; 
  incr[33] = fl[33]*wxl+(0.4391550328268399*fl[23]+0.253546276418555*fl[13])*dvxl; 
  incr[34] = fl[34]*wxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 
  outr[15] += incr[15]*dxr; 
  outr[16] += incr[16]*dxr; 
  outr[17] += incr[17]*dxr; 
  outr[18] += incr[18]*dxr; 
  outr[19] += incr[19]*dxr; 
  outr[20] += incr[20]*dxr; 
  outr[21] += incr[21]*dxr; 
  outr[22] += incr[22]*dxr; 
  outr[23] += incr[23]*dxr; 
  outr[24] += incr[24]*dxr; 
  outr[25] += incr[25]*dxr; 
  outr[26] += incr[26]*dxr; 
  outr[27] += incr[27]*dxr; 
  outr[28] += incr[28]*dxr; 
  outr[29] += incr[29]*dxr; 
  outr[30] += incr[30]*dxr; 
  outr[31] += incr[31]*dxr; 
  outr[32] += incr[32]*dxr; 
  outr[33] += incr[33]*dxr; 
  outr[34] += incr[34]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  outl[16] += incr[16]*dxl; 
  outl[17] += incr[17]*dxl; 
  outl[18] += -1.0*incr[18]*dxl; 
  outl[19] += -1.0*incr[19]*dxl; 
  outl[20] += incr[20]*dxl; 
  outl[21] += -1.0*incr[21]*dxl; 
  outl[22] += -1.0*incr[22]*dxl; 
  outl[23] += incr[23]*dxl; 
  outl[24] += -1.0*incr[24]*dxl; 
  outl[25] += -1.0*incr[25]*dxl; 
  outl[26] += -1.0*incr[26]*dxl; 
  outl[27] += -1.0*incr[27]*dxl; 
  outl[28] += incr[28]*dxl; 
  outl[29] += -1.0*incr[29]*dxl; 
  outl[30] += -1.0*incr[30]*dxl; 
  outl[31] += incr[31]*dxl; 
  outl[32] += -1.0*incr[32]*dxl; 
  outl[33] += -1.0*incr[33]*dxl; 
  outl[34] += -1.0*incr[34]*dxl; 
  } else { 
  incr[0] = ((-2.645751311064591*fr[31])+2.23606797749979*fr[11]-1.732050807568877*fr[1]+fr[0])*wxr+(0.6454972243679028*fr[21]-0.5*fr[6]+0.2886751345948129*fr[3])*dvxr; 
  incr[1] = (4.58257569495584*fr[31]-3.872983346207417*fr[11]+3.0*fr[1]-1.732050807568877*fr[0])*wxr+((-1.118033988749895*fr[21])+0.8660254037844386*fr[6]-0.5*fr[3])*dvxr; 
  incr[2] = (2.23606797749979*fr[19]-1.732050807568877*fr[5]+fr[2])*wxr+(0.2886751345948129*fr[7]-0.5*fr[15])*dvxr; 
  incr[3] = (2.23606797749979*fr[21]-1.732050807568877*fr[6]+fr[3])*wxr+((-0.7637626158259735*fr[31])-0.447213595499958*fr[23]+0.2581988897471612*fr[13]+0.6454972243679029*fr[11]-0.5*fr[1]+0.2886751345948129*fr[0])*dvxr; 
  incr[4] = (2.23606797749979*fr[25]-1.732050807568877*fr[8]+fr[4])*wxr+(0.2886751345948129*fr[10]-0.5*fr[17])*dvxr; 
  incr[5] = ((-3.872983346207417*fr[19])+3.0*fr[5]-1.732050807568877*fr[2])*wxr+(0.8660254037844386*fr[15]-0.5*fr[7])*dvxr; 
  incr[6] = ((-3.872983346207417*fr[21])+3.0*fr[6]-1.732050807568877*fr[3])*wxr+(1.322875655532295*fr[31]+0.7745966692414834*fr[23]-0.4472135954999579*fr[13]-1.118033988749895*fr[11]+0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[7] = (fr[7]-1.732050807568877*fr[15])*wxr+(0.2581988897471611*fr[24]+0.6454972243679028*fr[19]-0.5*fr[5]+0.2886751345948129*fr[2])*dvxr; 
  incr[8] = ((-3.872983346207417*fr[25])+3.0*fr[8]-1.732050807568877*fr[4])*wxr+(0.8660254037844386*fr[17]-0.5*fr[10])*dvxr; 
  incr[9] = (fr[9]-1.732050807568877*fr[16])*wxr+0.2886751345948129*fr[18]*dvxr; 
  incr[10] = (fr[10]-1.732050807568877*fr[17])*wxr+(0.2581988897471611*fr[27]+0.6454972243679028*fr[25]-0.5*fr[8]+0.2886751345948129*fr[4])*dvxr; 
  incr[11] = ((-5.916079783099617*fr[31])+5.0*fr[11]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*wxr+(1.443375672974065*fr[21]-1.118033988749895*fr[6]+0.6454972243679029*fr[3])*dvxr; 
  incr[12] = (fr[12]-1.732050807568877*fr[20])*wxr+0.2886751345948129*fr[22]*dvxr; 
  incr[13] = (fr[13]-1.732050807568877*fr[23])*wxr+(0.253546276418555*fr[33]+0.5773502691896257*fr[21]-0.4472135954999579*fr[6]+0.2581988897471612*fr[3])*dvxr; 
  incr[14] = (fr[14]-1.732050807568877*fr[28])*wxr+0.2886751345948129*fr[30]*dvxr; 
  incr[15] = (3.0*fr[15]-1.732050807568877*fr[7])*wxr+((-0.447213595499958*fr[24])-1.118033988749895*fr[19]+0.8660254037844386*fr[5]-0.5*fr[2])*dvxr; 
  incr[16] = (3.0*fr[16]-1.732050807568877*fr[9])*wxr-0.5*fr[18]*dvxr; 
  incr[17] = (3.0*fr[17]-1.732050807568877*fr[10])*wxr+((-0.447213595499958*fr[27])-1.118033988749895*fr[25]+0.8660254037844386*fr[8]-0.5*fr[4])*dvxr; 
  incr[18] = fr[18]*wxr+(0.2886751345948129*fr[9]-0.5*fr[16])*dvxr; 
  incr[19] = (5.0*fr[19]-3.872983346207417*fr[5]+2.23606797749979*fr[2])*wxr+(0.6454972243679028*fr[7]-1.118033988749895*fr[15])*dvxr; 
  incr[20] = (3.0*fr[20]-1.732050807568877*fr[12])*wxr-0.5*fr[22]*dvxr; 
  incr[21] = (5.0*fr[21]-3.872983346207417*fr[6]+2.23606797749979*fr[3])*wxr+((-1.707825127659933*fr[31])-1.0*fr[23]+0.5773502691896257*fr[13]+1.443375672974065*fr[11]-1.118033988749895*fr[1]+0.6454972243679028*fr[0])*dvxr; 
  incr[22] = fr[22]*wxr+(0.2886751345948129*fr[12]-0.5*fr[20])*dvxr; 
  incr[23] = (3.0*fr[23]-1.732050807568877*fr[13])*wxr+((-0.4391550328268399*fr[33])-1.0*fr[21]+0.7745966692414834*fr[6]-0.447213595499958*fr[3])*dvxr; 
  incr[24] = fr[24]*wxr+(0.2581988897471611*fr[7]-0.447213595499958*fr[15])*dvxr; 
  incr[25] = (5.0*fr[25]-3.872983346207417*fr[8]+2.23606797749979*fr[4])*wxr+(0.6454972243679028*fr[10]-1.118033988749895*fr[17])*dvxr; 
  incr[26] = fr[26]*wxr; 
  incr[27] = fr[27]*wxr+(0.2581988897471611*fr[10]-0.447213595499958*fr[17])*dvxr; 
  incr[28] = (3.0*fr[28]-1.732050807568877*fr[14])*wxr-0.5*fr[30]*dvxr; 
  incr[29] = fr[29]*wxr; 
  incr[30] = fr[30]*wxr+(0.2886751345948129*fr[14]-0.5*fr[28])*dvxr; 
  incr[31] = (7.0*fr[31]-5.916079783099617*fr[11]+4.58257569495584*fr[1]-2.645751311064591*fr[0])*wxr+((-1.707825127659933*fr[21])+1.322875655532295*fr[6]-0.7637626158259735*fr[3])*dvxr; 
  incr[32] = fr[32]*wxr; 
  incr[33] = fr[33]*wxr+(0.253546276418555*fr[13]-0.4391550328268399*fr[23])*dvxr; 
  incr[34] = fr[34]*wxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 
  outr[15] += incr[15]*dxr; 
  outr[16] += incr[16]*dxr; 
  outr[17] += incr[17]*dxr; 
  outr[18] += incr[18]*dxr; 
  outr[19] += incr[19]*dxr; 
  outr[20] += incr[20]*dxr; 
  outr[21] += incr[21]*dxr; 
  outr[22] += incr[22]*dxr; 
  outr[23] += incr[23]*dxr; 
  outr[24] += incr[24]*dxr; 
  outr[25] += incr[25]*dxr; 
  outr[26] += incr[26]*dxr; 
  outr[27] += incr[27]*dxr; 
  outr[28] += incr[28]*dxr; 
  outr[29] += incr[29]*dxr; 
  outr[30] += incr[30]*dxr; 
  outr[31] += incr[31]*dxr; 
  outr[32] += incr[32]*dxr; 
  outr[33] += incr[33]*dxr; 
  outr[34] += incr[34]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  outl[16] += incr[16]*dxl; 
  outl[17] += incr[17]*dxl; 
  outl[18] += -1.0*incr[18]*dxl; 
  outl[19] += -1.0*incr[19]*dxl; 
  outl[20] += incr[20]*dxl; 
  outl[21] += -1.0*incr[21]*dxl; 
  outl[22] += -1.0*incr[22]*dxl; 
  outl[23] += incr[23]*dxl; 
  outl[24] += -1.0*incr[24]*dxl; 
  outl[25] += -1.0*incr[25]*dxl; 
  outl[26] += -1.0*incr[26]*dxl; 
  outl[27] += -1.0*incr[27]*dxl; 
  outl[28] += incr[28]*dxl; 
  outl[29] += -1.0*incr[29]*dxl; 
  outl[30] += -1.0*incr[30]*dxl; 
  outl[31] += incr[31]*dxl; 
  outl[32] += -1.0*incr[32]*dxl; 
  outl[33] += -1.0*incr[33]*dxl; 
  outl[34] += -1.0*incr[34]*dxl; 
  } 
} 
void VlasovSurfStream2x2vMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[3]; 
  double wxl = wl[3]; 

  double dvxr = dxvr[3]; 
  double wxr = wr[3]; 

  double dxl = 1.0/dxvl[1]; 
  double dxr = 1.0/dxvr[1]; 

  double incr[5]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[2]+fl[0])*wxl+0.2886751345948129*fl[4]*dvxl; 
  incr[1] = fl[1]*wxl; 
  incr[2] = ((-3.0*fl[2])-1.732050807568877*fl[0])*wxl-0.5*fl[4]*dvxl; 
  incr[3] = fl[3]*wxl; 
  incr[4] = fl[4]*wxl+(0.5*fl[2]+0.2886751345948129*fl[0])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  } else { 
  incr[0] = (fr[0]-1.732050807568877*fr[2])*wxr+0.2886751345948129*fr[4]*dvxr; 
  incr[1] = fr[1]*wxr; 
  incr[2] = (3.0*fr[2]-1.732050807568877*fr[0])*wxr-0.5*fr[4]*dvxr; 
  incr[3] = fr[3]*wxr; 
  incr[4] = fr[4]*wxr+(0.2886751345948129*fr[0]-0.5*fr[2])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  } 
} 
void VlasovSurfStream2x2vMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[3]; 
  double wxl = wl[3]; 

  double dvxr = dxvr[3]; 
  double wxr = wr[3]; 

  double dxl = 1.0/dxvl[1]; 
  double dxr = 1.0/dxvr[1]; 

  double incr[15]; 

  if (wxr>0) { 
  incr[0] = (2.23606797749979*fl[12]+1.732050807568877*fl[2]+fl[0])*wxl+(0.5*fl[9]+0.2886751345948129*fl[4])*dvxl; 
  incr[1] = (1.732050807568877*fl[5]+fl[1])*wxl+0.2886751345948129*fl[8]*dvxl; 
  incr[2] = ((-3.872983346207417*fl[12])-3.0*fl[2]-1.732050807568877*fl[0])*wxl+((-0.8660254037844386*fl[9])-0.5*fl[4])*dvxl; 
  incr[3] = (1.732050807568877*fl[7]+fl[3])*wxl+0.2886751345948129*fl[10]*dvxl; 
  incr[4] = (1.732050807568877*fl[9]+fl[4])*wxl+(0.2581988897471612*fl[14]+0.6454972243679029*fl[12]+0.5*fl[2]+0.2886751345948129*fl[0])*dvxl; 
  incr[5] = ((-3.0*fl[5])-1.732050807568877*fl[1])*wxl-0.5*fl[8]*dvxl; 
  incr[6] = fl[6]*wxl; 
  incr[7] = ((-3.0*fl[7])-1.732050807568877*fl[3])*wxl-0.5*fl[10]*dvxl; 
  incr[8] = fl[8]*wxl+(0.5*fl[5]+0.2886751345948129*fl[1])*dvxl; 
  incr[9] = ((-3.0*fl[9])-1.732050807568877*fl[4])*wxl+((-0.4472135954999579*fl[14])-1.118033988749895*fl[12]-0.8660254037844386*fl[2]-0.5*fl[0])*dvxl; 
  incr[10] = fl[10]*wxl+(0.5*fl[7]+0.2886751345948129*fl[3])*dvxl; 
  incr[11] = fl[11]*wxl; 
  incr[12] = (5.0*fl[12]+3.872983346207417*fl[2]+2.23606797749979*fl[0])*wxl+(1.118033988749895*fl[9]+0.6454972243679029*fl[4])*dvxl; 
  incr[13] = fl[13]*wxl; 
  incr[14] = fl[14]*wxl+(0.4472135954999579*fl[9]+0.2581988897471612*fl[4])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  outl[9] += incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  } else { 
  incr[0] = (2.23606797749979*fr[12]-1.732050807568877*fr[2]+fr[0])*wxr+(0.2886751345948129*fr[4]-0.5*fr[9])*dvxr; 
  incr[1] = (fr[1]-1.732050807568877*fr[5])*wxr+0.2886751345948129*fr[8]*dvxr; 
  incr[2] = ((-3.872983346207417*fr[12])+3.0*fr[2]-1.732050807568877*fr[0])*wxr+(0.8660254037844386*fr[9]-0.5*fr[4])*dvxr; 
  incr[3] = (fr[3]-1.732050807568877*fr[7])*wxr+0.2886751345948129*fr[10]*dvxr; 
  incr[4] = (fr[4]-1.732050807568877*fr[9])*wxr+(0.2581988897471612*fr[14]+0.6454972243679029*fr[12]-0.5*fr[2]+0.2886751345948129*fr[0])*dvxr; 
  incr[5] = (3.0*fr[5]-1.732050807568877*fr[1])*wxr-0.5*fr[8]*dvxr; 
  incr[6] = fr[6]*wxr; 
  incr[7] = (3.0*fr[7]-1.732050807568877*fr[3])*wxr-0.5*fr[10]*dvxr; 
  incr[8] = fr[8]*wxr+(0.2886751345948129*fr[1]-0.5*fr[5])*dvxr; 
  incr[9] = (3.0*fr[9]-1.732050807568877*fr[4])*wxr+((-0.4472135954999579*fr[14])-1.118033988749895*fr[12]+0.8660254037844386*fr[2]-0.5*fr[0])*dvxr; 
  incr[10] = fr[10]*wxr+(0.2886751345948129*fr[3]-0.5*fr[7])*dvxr; 
  incr[11] = fr[11]*wxr; 
  incr[12] = (5.0*fr[12]-3.872983346207417*fr[2]+2.23606797749979*fr[0])*wxr+(0.6454972243679029*fr[4]-1.118033988749895*fr[9])*dvxr; 
  incr[13] = fr[13]*wxr; 
  incr[14] = fr[14]*wxr+(0.2581988897471612*fr[4]-0.4472135954999579*fr[9])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  outl[9] += incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  } 
} 
void VlasovSurfStream2x2vMax_Y_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[3]; 
  double wxl = wl[3]; 

  double dvxr = dxvr[3]; 
  double wxr = wr[3]; 

  double dxl = 1.0/dxvl[1]; 
  double dxr = 1.0/dxvr[1]; 

  double incr[35]; 

  if (wxr>0) { 
  incr[0] = (2.645751311064591*fl[32]+2.23606797749979*fl[12]+1.732050807568877*fl[2]+fl[0])*wxl+(0.6454972243679028*fl[26]+0.5*fl[9]+0.2886751345948129*fl[4])*dvxl; 
  incr[1] = (2.23606797749979*fl[20]+1.732050807568877*fl[5]+fl[1])*wxl+(0.5*fl[16]+0.2886751345948129*fl[8])*dvxl; 
  incr[2] = ((-4.58257569495584*fl[32])-3.872983346207417*fl[12]-3.0*fl[2]-1.732050807568877*fl[0])*wxl+((-1.118033988749895*fl[26])-0.8660254037844386*fl[9]-0.5*fl[4])*dvxl; 
  incr[3] = (2.23606797749979*fl[22]+1.732050807568877*fl[7]+fl[3])*wxl+(0.5*fl[18]+0.2886751345948129*fl[10])*dvxl; 
  incr[4] = (2.23606797749979*fl[26]+1.732050807568877*fl[9]+fl[4])*wxl+(0.7637626158259735*fl[32]+0.447213595499958*fl[29]+0.2581988897471612*fl[14]+0.6454972243679029*fl[12]+0.5*fl[2]+0.2886751345948129*fl[0])*dvxl; 
  incr[5] = ((-3.872983346207417*fl[20])-3.0*fl[5]-1.732050807568877*fl[1])*wxl+((-0.8660254037844386*fl[16])-0.5*fl[8])*dvxl; 
  incr[6] = (1.732050807568877*fl[15]+fl[6])*wxl+0.2886751345948129*fl[17]*dvxl; 
  incr[7] = ((-3.872983346207417*fl[22])-3.0*fl[7]-1.732050807568877*fl[3])*wxl+((-0.8660254037844386*fl[18])-0.5*fl[10])*dvxl; 
  incr[8] = (1.732050807568877*fl[16]+fl[8])*wxl+(0.2581988897471611*fl[28]+0.6454972243679028*fl[20]+0.5*fl[5]+0.2886751345948129*fl[1])*dvxl; 
  incr[9] = ((-3.872983346207417*fl[26])-3.0*fl[9]-1.732050807568877*fl[4])*wxl+((-1.322875655532295*fl[32])-0.7745966692414834*fl[29]-0.4472135954999579*fl[14]-1.118033988749895*fl[12]-0.8660254037844386*fl[2]-0.5*fl[0])*dvxl; 
  incr[10] = (1.732050807568877*fl[18]+fl[10])*wxl+(0.2581988897471611*fl[30]+0.6454972243679028*fl[22]+0.5*fl[7]+0.2886751345948129*fl[3])*dvxl; 
  incr[11] = (1.732050807568877*fl[19]+fl[11])*wxl+0.2886751345948129*fl[25]*dvxl; 
  incr[12] = (5.916079783099617*fl[32]+5.0*fl[12]+3.872983346207417*fl[2]+2.23606797749979*fl[0])*wxl+(1.443375672974065*fl[26]+1.118033988749895*fl[9]+0.6454972243679029*fl[4])*dvxl; 
  incr[13] = (1.732050807568877*fl[24]+fl[13])*wxl+0.2886751345948129*fl[27]*dvxl; 
  incr[14] = (1.732050807568877*fl[29]+fl[14])*wxl+(0.253546276418555*fl[34]+0.5773502691896257*fl[26]+0.4472135954999579*fl[9]+0.2581988897471612*fl[4])*dvxl; 
  incr[15] = ((-3.0*fl[15])-1.732050807568877*fl[6])*wxl-0.5*fl[17]*dvxl; 
  incr[16] = ((-3.0*fl[16])-1.732050807568877*fl[8])*wxl+((-0.447213595499958*fl[28])-1.118033988749895*fl[20]-0.8660254037844386*fl[5]-0.5*fl[1])*dvxl; 
  incr[17] = fl[17]*wxl+(0.5*fl[15]+0.2886751345948129*fl[6])*dvxl; 
  incr[18] = ((-3.0*fl[18])-1.732050807568877*fl[10])*wxl+((-0.447213595499958*fl[30])-1.118033988749895*fl[22]-0.8660254037844386*fl[7]-0.5*fl[3])*dvxl; 
  incr[19] = ((-3.0*fl[19])-1.732050807568877*fl[11])*wxl-0.5*fl[25]*dvxl; 
  incr[20] = (5.0*fl[20]+3.872983346207417*fl[5]+2.23606797749979*fl[1])*wxl+(1.118033988749895*fl[16]+0.6454972243679028*fl[8])*dvxl; 
  incr[21] = fl[21]*wxl; 
  incr[22] = (5.0*fl[22]+3.872983346207417*fl[7]+2.23606797749979*fl[3])*wxl+(1.118033988749895*fl[18]+0.6454972243679028*fl[10])*dvxl; 
  incr[23] = fl[23]*wxl; 
  incr[24] = ((-3.0*fl[24])-1.732050807568877*fl[13])*wxl-0.5*fl[27]*dvxl; 
  incr[25] = fl[25]*wxl+(0.5*fl[19]+0.2886751345948129*fl[11])*dvxl; 
  incr[26] = (5.0*fl[26]+3.872983346207417*fl[9]+2.23606797749979*fl[4])*wxl+(1.707825127659933*fl[32]+fl[29]+0.5773502691896257*fl[14]+1.443375672974065*fl[12]+1.118033988749895*fl[2]+0.6454972243679028*fl[0])*dvxl; 
  incr[27] = fl[27]*wxl+(0.5*fl[24]+0.2886751345948129*fl[13])*dvxl; 
  incr[28] = fl[28]*wxl+(0.447213595499958*fl[16]+0.2581988897471611*fl[8])*dvxl; 
  incr[29] = ((-3.0*fl[29])-1.732050807568877*fl[14])*wxl+((-0.4391550328268399*fl[34])-1.0*fl[26]-0.7745966692414834*fl[9]-0.447213595499958*fl[4])*dvxl; 
  incr[30] = fl[30]*wxl+(0.447213595499958*fl[18]+0.2581988897471611*fl[10])*dvxl; 
  incr[31] = fl[31]*wxl; 
  incr[32] = ((-7.0*fl[32])-5.916079783099617*fl[12]-4.58257569495584*fl[2]-2.645751311064591*fl[0])*wxl+((-1.707825127659933*fl[26])-1.322875655532295*fl[9]-0.7637626158259735*fl[4])*dvxl; 
  incr[33] = fl[33]*wxl; 
  incr[34] = fl[34]*wxl+(0.4391550328268399*fl[29]+0.253546276418555*fl[14])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 
  outr[15] += incr[15]*dxr; 
  outr[16] += incr[16]*dxr; 
  outr[17] += incr[17]*dxr; 
  outr[18] += incr[18]*dxr; 
  outr[19] += incr[19]*dxr; 
  outr[20] += incr[20]*dxr; 
  outr[21] += incr[21]*dxr; 
  outr[22] += incr[22]*dxr; 
  outr[23] += incr[23]*dxr; 
  outr[24] += incr[24]*dxr; 
  outr[25] += incr[25]*dxr; 
  outr[26] += incr[26]*dxr; 
  outr[27] += incr[27]*dxr; 
  outr[28] += incr[28]*dxr; 
  outr[29] += incr[29]*dxr; 
  outr[30] += incr[30]*dxr; 
  outr[31] += incr[31]*dxr; 
  outr[32] += incr[32]*dxr; 
  outr[33] += incr[33]*dxr; 
  outr[34] += incr[34]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  outl[9] += incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  outl[16] += incr[16]*dxl; 
  outl[17] += -1.0*incr[17]*dxl; 
  outl[18] += incr[18]*dxl; 
  outl[19] += incr[19]*dxl; 
  outl[20] += -1.0*incr[20]*dxl; 
  outl[21] += -1.0*incr[21]*dxl; 
  outl[22] += -1.0*incr[22]*dxl; 
  outl[23] += -1.0*incr[23]*dxl; 
  outl[24] += incr[24]*dxl; 
  outl[25] += -1.0*incr[25]*dxl; 
  outl[26] += -1.0*incr[26]*dxl; 
  outl[27] += -1.0*incr[27]*dxl; 
  outl[28] += -1.0*incr[28]*dxl; 
  outl[29] += incr[29]*dxl; 
  outl[30] += -1.0*incr[30]*dxl; 
  outl[31] += -1.0*incr[31]*dxl; 
  outl[32] += incr[32]*dxl; 
  outl[33] += -1.0*incr[33]*dxl; 
  outl[34] += -1.0*incr[34]*dxl; 
  } else { 
  incr[0] = ((-2.645751311064591*fr[32])+2.23606797749979*fr[12]-1.732050807568877*fr[2]+fr[0])*wxr+(0.6454972243679028*fr[26]-0.5*fr[9]+0.2886751345948129*fr[4])*dvxr; 
  incr[1] = (2.23606797749979*fr[20]-1.732050807568877*fr[5]+fr[1])*wxr+(0.2886751345948129*fr[8]-0.5*fr[16])*dvxr; 
  incr[2] = (4.58257569495584*fr[32]-3.872983346207417*fr[12]+3.0*fr[2]-1.732050807568877*fr[0])*wxr+((-1.118033988749895*fr[26])+0.8660254037844386*fr[9]-0.5*fr[4])*dvxr; 
  incr[3] = (2.23606797749979*fr[22]-1.732050807568877*fr[7]+fr[3])*wxr+(0.2886751345948129*fr[10]-0.5*fr[18])*dvxr; 
  incr[4] = (2.23606797749979*fr[26]-1.732050807568877*fr[9]+fr[4])*wxr+((-0.7637626158259735*fr[32])-0.447213595499958*fr[29]+0.2581988897471612*fr[14]+0.6454972243679029*fr[12]-0.5*fr[2]+0.2886751345948129*fr[0])*dvxr; 
  incr[5] = ((-3.872983346207417*fr[20])+3.0*fr[5]-1.732050807568877*fr[1])*wxr+(0.8660254037844386*fr[16]-0.5*fr[8])*dvxr; 
  incr[6] = (fr[6]-1.732050807568877*fr[15])*wxr+0.2886751345948129*fr[17]*dvxr; 
  incr[7] = ((-3.872983346207417*fr[22])+3.0*fr[7]-1.732050807568877*fr[3])*wxr+(0.8660254037844386*fr[18]-0.5*fr[10])*dvxr; 
  incr[8] = (fr[8]-1.732050807568877*fr[16])*wxr+(0.2581988897471611*fr[28]+0.6454972243679028*fr[20]-0.5*fr[5]+0.2886751345948129*fr[1])*dvxr; 
  incr[9] = ((-3.872983346207417*fr[26])+3.0*fr[9]-1.732050807568877*fr[4])*wxr+(1.322875655532295*fr[32]+0.7745966692414834*fr[29]-0.4472135954999579*fr[14]-1.118033988749895*fr[12]+0.8660254037844386*fr[2]-0.5*fr[0])*dvxr; 
  incr[10] = (fr[10]-1.732050807568877*fr[18])*wxr+(0.2581988897471611*fr[30]+0.6454972243679028*fr[22]-0.5*fr[7]+0.2886751345948129*fr[3])*dvxr; 
  incr[11] = (fr[11]-1.732050807568877*fr[19])*wxr+0.2886751345948129*fr[25]*dvxr; 
  incr[12] = ((-5.916079783099617*fr[32])+5.0*fr[12]-3.872983346207417*fr[2]+2.23606797749979*fr[0])*wxr+(1.443375672974065*fr[26]-1.118033988749895*fr[9]+0.6454972243679029*fr[4])*dvxr; 
  incr[13] = (fr[13]-1.732050807568877*fr[24])*wxr+0.2886751345948129*fr[27]*dvxr; 
  incr[14] = (fr[14]-1.732050807568877*fr[29])*wxr+(0.253546276418555*fr[34]+0.5773502691896257*fr[26]-0.4472135954999579*fr[9]+0.2581988897471612*fr[4])*dvxr; 
  incr[15] = (3.0*fr[15]-1.732050807568877*fr[6])*wxr-0.5*fr[17]*dvxr; 
  incr[16] = (3.0*fr[16]-1.732050807568877*fr[8])*wxr+((-0.447213595499958*fr[28])-1.118033988749895*fr[20]+0.8660254037844386*fr[5]-0.5*fr[1])*dvxr; 
  incr[17] = fr[17]*wxr+(0.2886751345948129*fr[6]-0.5*fr[15])*dvxr; 
  incr[18] = (3.0*fr[18]-1.732050807568877*fr[10])*wxr+((-0.447213595499958*fr[30])-1.118033988749895*fr[22]+0.8660254037844386*fr[7]-0.5*fr[3])*dvxr; 
  incr[19] = (3.0*fr[19]-1.732050807568877*fr[11])*wxr-0.5*fr[25]*dvxr; 
  incr[20] = (5.0*fr[20]-3.872983346207417*fr[5]+2.23606797749979*fr[1])*wxr+(0.6454972243679028*fr[8]-1.118033988749895*fr[16])*dvxr; 
  incr[21] = fr[21]*wxr; 
  incr[22] = (5.0*fr[22]-3.872983346207417*fr[7]+2.23606797749979*fr[3])*wxr+(0.6454972243679028*fr[10]-1.118033988749895*fr[18])*dvxr; 
  incr[23] = fr[23]*wxr; 
  incr[24] = (3.0*fr[24]-1.732050807568877*fr[13])*wxr-0.5*fr[27]*dvxr; 
  incr[25] = fr[25]*wxr+(0.2886751345948129*fr[11]-0.5*fr[19])*dvxr; 
  incr[26] = (5.0*fr[26]-3.872983346207417*fr[9]+2.23606797749979*fr[4])*wxr+((-1.707825127659933*fr[32])-1.0*fr[29]+0.5773502691896257*fr[14]+1.443375672974065*fr[12]-1.118033988749895*fr[2]+0.6454972243679028*fr[0])*dvxr; 
  incr[27] = fr[27]*wxr+(0.2886751345948129*fr[13]-0.5*fr[24])*dvxr; 
  incr[28] = fr[28]*wxr+(0.2581988897471611*fr[8]-0.447213595499958*fr[16])*dvxr; 
  incr[29] = (3.0*fr[29]-1.732050807568877*fr[14])*wxr+((-0.4391550328268399*fr[34])-1.0*fr[26]+0.7745966692414834*fr[9]-0.447213595499958*fr[4])*dvxr; 
  incr[30] = fr[30]*wxr+(0.2581988897471611*fr[10]-0.447213595499958*fr[18])*dvxr; 
  incr[31] = fr[31]*wxr; 
  incr[32] = (7.0*fr[32]-5.916079783099617*fr[12]+4.58257569495584*fr[2]-2.645751311064591*fr[0])*wxr+((-1.707825127659933*fr[26])+1.322875655532295*fr[9]-0.7637626158259735*fr[4])*dvxr; 
  incr[33] = fr[33]*wxr; 
  incr[34] = fr[34]*wxr+(0.253546276418555*fr[14]-0.4391550328268399*fr[29])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 
  outr[7] += incr[7]*dxr; 
  outr[8] += incr[8]*dxr; 
  outr[9] += incr[9]*dxr; 
  outr[10] += incr[10]*dxr; 
  outr[11] += incr[11]*dxr; 
  outr[12] += incr[12]*dxr; 
  outr[13] += incr[13]*dxr; 
  outr[14] += incr[14]*dxr; 
  outr[15] += incr[15]*dxr; 
  outr[16] += incr[16]*dxr; 
  outr[17] += incr[17]*dxr; 
  outr[18] += incr[18]*dxr; 
  outr[19] += incr[19]*dxr; 
  outr[20] += incr[20]*dxr; 
  outr[21] += incr[21]*dxr; 
  outr[22] += incr[22]*dxr; 
  outr[23] += incr[23]*dxr; 
  outr[24] += incr[24]*dxr; 
  outr[25] += incr[25]*dxr; 
  outr[26] += incr[26]*dxr; 
  outr[27] += incr[27]*dxr; 
  outr[28] += incr[28]*dxr; 
  outr[29] += incr[29]*dxr; 
  outr[30] += incr[30]*dxr; 
  outr[31] += incr[31]*dxr; 
  outr[32] += incr[32]*dxr; 
  outr[33] += incr[33]*dxr; 
  outr[34] += incr[34]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  outl[9] += incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  outl[16] += incr[16]*dxl; 
  outl[17] += -1.0*incr[17]*dxl; 
  outl[18] += incr[18]*dxl; 
  outl[19] += incr[19]*dxl; 
  outl[20] += -1.0*incr[20]*dxl; 
  outl[21] += -1.0*incr[21]*dxl; 
  outl[22] += -1.0*incr[22]*dxl; 
  outl[23] += -1.0*incr[23]*dxl; 
  outl[24] += incr[24]*dxl; 
  outl[25] += -1.0*incr[25]*dxl; 
  outl[26] += -1.0*incr[26]*dxl; 
  outl[27] += -1.0*incr[27]*dxl; 
  outl[28] += -1.0*incr[28]*dxl; 
  outl[29] += incr[29]*dxl; 
  outl[30] += -1.0*incr[30]*dxl; 
  outl[31] += -1.0*incr[31]*dxl; 
  outl[32] += incr[32]*dxl; 
  outl[33] += -1.0*incr[33]*dxl; 
  outl[34] += -1.0*incr[34]*dxl; 
  } 
} 
