#include <VlasovModDecl.h> 
void VlasovSurfStream2x3vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[2]; 
  double wxl = wl[2]; 

  double dvxr = dxvr[2]; 
  double wxr = wr[2]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[21]; 

  if (wxr>0) { 
  incr[0] = (2.23606797749979*fl[16]+1.732050807568877*fl[1]+fl[0])*wxl+(0.5*fl[7]+0.2886751345948129*fl[3])*dvxl; 
  incr[1] = ((-3.872983346207417*fl[16])-3.0*fl[1]-1.732050807568877*fl[0])*wxl+((-0.8660254037844386*fl[7])-0.5*fl[3])*dvxl; 
  incr[2] = (1.732050807568877*fl[6]+fl[2])*wxl+0.2886751345948129*fl[8]*dvxl; 
  incr[3] = (1.732050807568877*fl[7]+fl[3])*wxl+(0.2581988897471612*fl[18]+0.6454972243679029*fl[16]+0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[4] = (1.732050807568877*fl[9]+fl[4])*wxl+0.2886751345948129*fl[11]*dvxl; 
  incr[5] = (1.732050807568877*fl[12]+fl[5])*wxl+0.2886751345948129*fl[14]*dvxl; 
  incr[6] = ((-3.0*fl[6])-1.732050807568877*fl[2])*wxl-0.5*fl[8]*dvxl; 
  incr[7] = ((-3.0*fl[7])-1.732050807568877*fl[3])*wxl+((-0.4472135954999579*fl[18])-1.118033988749895*fl[16]-0.8660254037844386*fl[1]-0.5*fl[0])*dvxl; 
  incr[8] = fl[8]*wxl+(0.5*fl[6]+0.2886751345948129*fl[2])*dvxl; 
  incr[9] = ((-3.0*fl[9])-1.732050807568877*fl[4])*wxl-0.5*fl[11]*dvxl; 
  incr[10] = fl[10]*wxl; 
  incr[11] = fl[11]*wxl+(0.5*fl[9]+0.2886751345948129*fl[4])*dvxl; 
  incr[12] = ((-3.0*fl[12])-1.732050807568877*fl[5])*wxl-0.5*fl[14]*dvxl; 
  incr[13] = fl[13]*wxl; 
  incr[14] = fl[14]*wxl+(0.5*fl[12]+0.2886751345948129*fl[5])*dvxl; 
  incr[15] = fl[15]*wxl; 
  incr[16] = (5.0*fl[16]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*wxl+(1.118033988749895*fl[7]+0.6454972243679029*fl[3])*dvxl; 
  incr[17] = fl[17]*wxl; 
  incr[18] = fl[18]*wxl+(0.4472135954999579*fl[7]+0.2581988897471612*fl[3])*dvxl; 
  incr[19] = fl[19]*wxl; 
  incr[20] = fl[20]*wxl; 

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

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  outl[9] += incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += -1.0*incr[15]*dxl; 
  outl[16] += -1.0*incr[16]*dxl; 
  outl[17] += -1.0*incr[17]*dxl; 
  outl[18] += -1.0*incr[18]*dxl; 
  outl[19] += -1.0*incr[19]*dxl; 
  outl[20] += -1.0*incr[20]*dxl; 
  } else { 
  incr[0] = (2.23606797749979*fr[16]-1.732050807568877*fr[1]+fr[0])*wxr+(0.2886751345948129*fr[3]-0.5*fr[7])*dvxr; 
  incr[1] = ((-3.872983346207417*fr[16])+3.0*fr[1]-1.732050807568877*fr[0])*wxr+(0.8660254037844386*fr[7]-0.5*fr[3])*dvxr; 
  incr[2] = (fr[2]-1.732050807568877*fr[6])*wxr+0.2886751345948129*fr[8]*dvxr; 
  incr[3] = (fr[3]-1.732050807568877*fr[7])*wxr+(0.2581988897471612*fr[18]+0.6454972243679029*fr[16]-0.5*fr[1]+0.2886751345948129*fr[0])*dvxr; 
  incr[4] = (fr[4]-1.732050807568877*fr[9])*wxr+0.2886751345948129*fr[11]*dvxr; 
  incr[5] = (fr[5]-1.732050807568877*fr[12])*wxr+0.2886751345948129*fr[14]*dvxr; 
  incr[6] = (3.0*fr[6]-1.732050807568877*fr[2])*wxr-0.5*fr[8]*dvxr; 
  incr[7] = (3.0*fr[7]-1.732050807568877*fr[3])*wxr+((-0.4472135954999579*fr[18])-1.118033988749895*fr[16]+0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[8] = fr[8]*wxr+(0.2886751345948129*fr[2]-0.5*fr[6])*dvxr; 
  incr[9] = (3.0*fr[9]-1.732050807568877*fr[4])*wxr-0.5*fr[11]*dvxr; 
  incr[10] = fr[10]*wxr; 
  incr[11] = fr[11]*wxr+(0.2886751345948129*fr[4]-0.5*fr[9])*dvxr; 
  incr[12] = (3.0*fr[12]-1.732050807568877*fr[5])*wxr-0.5*fr[14]*dvxr; 
  incr[13] = fr[13]*wxr; 
  incr[14] = fr[14]*wxr+(0.2886751345948129*fr[5]-0.5*fr[12])*dvxr; 
  incr[15] = fr[15]*wxr; 
  incr[16] = (5.0*fr[16]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*wxr+(0.6454972243679029*fr[3]-1.118033988749895*fr[7])*dvxr; 
  incr[17] = fr[17]*wxr; 
  incr[18] = fr[18]*wxr+(0.2581988897471612*fr[3]-0.4472135954999579*fr[7])*dvxr; 
  incr[19] = fr[19]*wxr; 
  incr[20] = fr[20]*wxr; 

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

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  outl[9] += incr[9]*dxl; 
  outl[10] += -1.0*incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += -1.0*incr[15]*dxl; 
  outl[16] += -1.0*incr[16]*dxl; 
  outl[17] += -1.0*incr[17]*dxl; 
  outl[18] += -1.0*incr[18]*dxl; 
  outl[19] += -1.0*incr[19]*dxl; 
  outl[20] += -1.0*incr[20]*dxl; 
  } 
} 
void VlasovSurfStream2x3vMax_Y_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[3]; 
  double wxl = wl[3]; 

  double dvxr = dxvr[3]; 
  double wxr = wr[3]; 

  double dxl = 1.0/dxvl[1]; 
  double dxr = 1.0/dxvr[1]; 

  double incr[21]; 

  if (wxr>0) { 
  incr[0] = (2.23606797749979*fl[17]+1.732050807568877*fl[2]+fl[0])*wxl+(0.5*fl[10]+0.2886751345948129*fl[4])*dvxl; 
  incr[1] = (1.732050807568877*fl[6]+fl[1])*wxl+0.2886751345948129*fl[9]*dvxl; 
  incr[2] = ((-3.872983346207417*fl[17])-3.0*fl[2]-1.732050807568877*fl[0])*wxl+((-0.8660254037844386*fl[10])-0.5*fl[4])*dvxl; 
  incr[3] = (1.732050807568877*fl[8]+fl[3])*wxl+0.2886751345948129*fl[11]*dvxl; 
  incr[4] = (1.732050807568877*fl[10]+fl[4])*wxl+(0.2581988897471612*fl[19]+0.6454972243679029*fl[17]+0.5*fl[2]+0.2886751345948129*fl[0])*dvxl; 
  incr[5] = (1.732050807568877*fl[13]+fl[5])*wxl+0.2886751345948129*fl[15]*dvxl; 
  incr[6] = ((-3.0*fl[6])-1.732050807568877*fl[1])*wxl-0.5*fl[9]*dvxl; 
  incr[7] = fl[7]*wxl; 
  incr[8] = ((-3.0*fl[8])-1.732050807568877*fl[3])*wxl-0.5*fl[11]*dvxl; 
  incr[9] = fl[9]*wxl+(0.5*fl[6]+0.2886751345948129*fl[1])*dvxl; 
  incr[10] = ((-3.0*fl[10])-1.732050807568877*fl[4])*wxl+((-0.4472135954999579*fl[19])-1.118033988749895*fl[17]-0.8660254037844386*fl[2]-0.5*fl[0])*dvxl; 
  incr[11] = fl[11]*wxl+(0.5*fl[8]+0.2886751345948129*fl[3])*dvxl; 
  incr[12] = fl[12]*wxl; 
  incr[13] = ((-3.0*fl[13])-1.732050807568877*fl[5])*wxl-0.5*fl[15]*dvxl; 
  incr[14] = fl[14]*wxl; 
  incr[15] = fl[15]*wxl+(0.5*fl[13]+0.2886751345948129*fl[5])*dvxl; 
  incr[16] = fl[16]*wxl; 
  incr[17] = (5.0*fl[17]+3.872983346207417*fl[2]+2.23606797749979*fl[0])*wxl+(1.118033988749895*fl[10]+0.6454972243679029*fl[4])*dvxl; 
  incr[18] = fl[18]*wxl; 
  incr[19] = fl[19]*wxl+(0.4472135954999579*fl[10]+0.2581988897471612*fl[4])*dvxl; 
  incr[20] = fl[20]*wxl; 

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

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += -1.0*incr[15]*dxl; 
  outl[16] += -1.0*incr[16]*dxl; 
  outl[17] += -1.0*incr[17]*dxl; 
  outl[18] += -1.0*incr[18]*dxl; 
  outl[19] += -1.0*incr[19]*dxl; 
  outl[20] += -1.0*incr[20]*dxl; 
  } else { 
  incr[0] = (2.23606797749979*fr[17]-1.732050807568877*fr[2]+fr[0])*wxr+(0.2886751345948129*fr[4]-0.5*fr[10])*dvxr; 
  incr[1] = (fr[1]-1.732050807568877*fr[6])*wxr+0.2886751345948129*fr[9]*dvxr; 
  incr[2] = ((-3.872983346207417*fr[17])+3.0*fr[2]-1.732050807568877*fr[0])*wxr+(0.8660254037844386*fr[10]-0.5*fr[4])*dvxr; 
  incr[3] = (fr[3]-1.732050807568877*fr[8])*wxr+0.2886751345948129*fr[11]*dvxr; 
  incr[4] = (fr[4]-1.732050807568877*fr[10])*wxr+(0.2581988897471612*fr[19]+0.6454972243679029*fr[17]-0.5*fr[2]+0.2886751345948129*fr[0])*dvxr; 
  incr[5] = (fr[5]-1.732050807568877*fr[13])*wxr+0.2886751345948129*fr[15]*dvxr; 
  incr[6] = (3.0*fr[6]-1.732050807568877*fr[1])*wxr-0.5*fr[9]*dvxr; 
  incr[7] = fr[7]*wxr; 
  incr[8] = (3.0*fr[8]-1.732050807568877*fr[3])*wxr-0.5*fr[11]*dvxr; 
  incr[9] = fr[9]*wxr+(0.2886751345948129*fr[1]-0.5*fr[6])*dvxr; 
  incr[10] = (3.0*fr[10]-1.732050807568877*fr[4])*wxr+((-0.4472135954999579*fr[19])-1.118033988749895*fr[17]+0.8660254037844386*fr[2]-0.5*fr[0])*dvxr; 
  incr[11] = fr[11]*wxr+(0.2886751345948129*fr[3]-0.5*fr[8])*dvxr; 
  incr[12] = fr[12]*wxr; 
  incr[13] = (3.0*fr[13]-1.732050807568877*fr[5])*wxr-0.5*fr[15]*dvxr; 
  incr[14] = fr[14]*wxr; 
  incr[15] = fr[15]*wxr+(0.2886751345948129*fr[5]-0.5*fr[13])*dvxr; 
  incr[16] = fr[16]*wxr; 
  incr[17] = (5.0*fr[17]-3.872983346207417*fr[2]+2.23606797749979*fr[0])*wxr+(0.6454972243679029*fr[4]-1.118033988749895*fr[10])*dvxr; 
  incr[18] = fr[18]*wxr; 
  incr[19] = fr[19]*wxr+(0.2581988897471612*fr[4]-0.4472135954999579*fr[10])*dvxr; 
  incr[20] = fr[20]*wxr; 

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

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += -1.0*incr[12]*dxl; 
  outl[13] += incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += -1.0*incr[15]*dxl; 
  outl[16] += -1.0*incr[16]*dxl; 
  outl[17] += -1.0*incr[17]*dxl; 
  outl[18] += -1.0*incr[18]*dxl; 
  outl[19] += -1.0*incr[19]*dxl; 
  outl[20] += -1.0*incr[20]*dxl; 
  } 
} 
