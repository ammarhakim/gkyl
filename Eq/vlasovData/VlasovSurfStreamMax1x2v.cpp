#include <VlasovModDecl.h> 
void VlasovSurfStream1x2vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[4]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[1]+fl[0])*wxl+0.2886751345948129*fl[2]*dvxl; 
  incr[1] = ((-3.0*fl[1])-1.732050807568877*fl[0])*wxl-0.5*fl[2]*dvxl; 
  incr[2] = fl[2]*wxl+(0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[3] = fl[3]*wxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  } else { 
  incr[0] = (fr[0]-1.732050807568877*fr[1])*wxr+0.2886751345948129*fr[2]*dvxr; 
  incr[1] = (3.0*fr[1]-1.732050807568877*fr[0])*wxr-0.5*fr[2]*dvxr; 
  incr[2] = fr[2]*wxr+(0.2886751345948129*fr[0]-0.5*fr[1])*dvxr; 
  incr[3] = fr[3]*wxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  } 
} 
void VlasovSurfStream1x2vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[10]; 

  if (wxr>0) { 
  incr[0] = (2.23606797749979*fl[7]+1.732050807568877*fl[1]+fl[0])*wxl+(0.5*fl[4]+0.2886751345948129*fl[2])*dvxl; 
  incr[1] = ((-3.872983346207417*fl[7])-3.0*fl[1]-1.732050807568877*fl[0])*wxl+((-0.8660254037844386*fl[4])-0.5*fl[2])*dvxl; 
  incr[2] = (1.732050807568877*fl[4]+fl[2])*wxl+(0.2581988897471612*fl[8]+0.6454972243679029*fl[7]+0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[3] = (1.732050807568877*fl[5]+fl[3])*wxl+0.2886751345948129*fl[6]*dvxl; 
  incr[4] = ((-3.0*fl[4])-1.732050807568877*fl[2])*wxl+((-0.4472135954999579*fl[8])-1.118033988749895*fl[7]-0.8660254037844386*fl[1]-0.5*fl[0])*dvxl; 
  incr[5] = ((-3.0*fl[5])-1.732050807568877*fl[3])*wxl-0.5*fl[6]*dvxl; 
  incr[6] = fl[6]*wxl+(0.5*fl[5]+0.2886751345948129*fl[3])*dvxl; 
  incr[7] = (5.0*fl[7]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*wxl+(1.118033988749895*fl[4]+0.6454972243679029*fl[2])*dvxl; 
  incr[8] = fl[8]*wxl+(0.4472135954999579*fl[4]+0.2581988897471612*fl[2])*dvxl; 
  incr[9] = fl[9]*wxl; 

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

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  } else { 
  incr[0] = (2.23606797749979*fr[7]-1.732050807568877*fr[1]+fr[0])*wxr+(0.2886751345948129*fr[2]-0.5*fr[4])*dvxr; 
  incr[1] = ((-3.872983346207417*fr[7])+3.0*fr[1]-1.732050807568877*fr[0])*wxr+(0.8660254037844386*fr[4]-0.5*fr[2])*dvxr; 
  incr[2] = (fr[2]-1.732050807568877*fr[4])*wxr+(0.2581988897471612*fr[8]+0.6454972243679029*fr[7]-0.5*fr[1]+0.2886751345948129*fr[0])*dvxr; 
  incr[3] = (fr[3]-1.732050807568877*fr[5])*wxr+0.2886751345948129*fr[6]*dvxr; 
  incr[4] = (3.0*fr[4]-1.732050807568877*fr[2])*wxr+((-0.4472135954999579*fr[8])-1.118033988749895*fr[7]+0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[5] = (3.0*fr[5]-1.732050807568877*fr[3])*wxr-0.5*fr[6]*dvxr; 
  incr[6] = fr[6]*wxr+(0.2886751345948129*fr[3]-0.5*fr[5])*dvxr; 
  incr[7] = (5.0*fr[7]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*wxr+(0.6454972243679029*fr[2]-1.118033988749895*fr[4])*dvxr; 
  incr[8] = fr[8]*wxr+(0.2581988897471612*fr[2]-0.4472135954999579*fr[4])*dvxr; 
  incr[9] = fr[9]*wxr; 

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

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  } 
} 
void VlasovSurfStream1x2vMax_X_P3(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[20]; 

  if (wxr>0) { 
  incr[0] = (2.645751311064591*fl[17]+2.23606797749979*fl[7]+1.732050807568877*fl[1]+fl[0])*wxl+(0.6454972243679028*fl[11]+0.5*fl[4]+0.2886751345948129*fl[2])*dvxl; 
  incr[1] = ((-4.58257569495584*fl[17])-3.872983346207417*fl[7]-3.0*fl[1]-1.732050807568877*fl[0])*wxl+((-1.118033988749895*fl[11])-0.8660254037844386*fl[4]-0.5*fl[2])*dvxl; 
  incr[2] = (2.23606797749979*fl[11]+1.732050807568877*fl[4]+fl[2])*wxl+(0.7637626158259735*fl[17]+0.447213595499958*fl[12]+0.2581988897471612*fl[8]+0.6454972243679029*fl[7]+0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[3] = (2.23606797749979*fl[13]+1.732050807568877*fl[5]+fl[3])*wxl+(0.5*fl[10]+0.2886751345948129*fl[6])*dvxl; 
  incr[4] = ((-3.872983346207417*fl[11])-3.0*fl[4]-1.732050807568877*fl[2])*wxl+((-1.322875655532295*fl[17])-0.7745966692414834*fl[12]-0.4472135954999579*fl[8]-1.118033988749895*fl[7]-0.8660254037844386*fl[1]-0.5*fl[0])*dvxl; 
  incr[5] = ((-3.872983346207417*fl[13])-3.0*fl[5]-1.732050807568877*fl[3])*wxl+((-0.8660254037844386*fl[10])-0.5*fl[6])*dvxl; 
  incr[6] = (1.732050807568877*fl[10]+fl[6])*wxl+(0.2581988897471611*fl[14]+0.6454972243679028*fl[13]+0.5*fl[5]+0.2886751345948129*fl[3])*dvxl; 
  incr[7] = (5.916079783099617*fl[17]+5.0*fl[7]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*wxl+(1.443375672974065*fl[11]+1.118033988749895*fl[4]+0.6454972243679029*fl[2])*dvxl; 
  incr[8] = (1.732050807568877*fl[12]+fl[8])*wxl+(0.253546276418555*fl[18]+0.5773502691896257*fl[11]+0.4472135954999579*fl[4]+0.2581988897471612*fl[2])*dvxl; 
  incr[9] = (1.732050807568877*fl[15]+fl[9])*wxl+0.2886751345948129*fl[16]*dvxl; 
  incr[10] = ((-3.0*fl[10])-1.732050807568877*fl[6])*wxl+((-0.447213595499958*fl[14])-1.118033988749895*fl[13]-0.8660254037844386*fl[5]-0.5*fl[3])*dvxl; 
  incr[11] = (5.0*fl[11]+3.872983346207417*fl[4]+2.23606797749979*fl[2])*wxl+(1.707825127659933*fl[17]+fl[12]+0.5773502691896257*fl[8]+1.443375672974065*fl[7]+1.118033988749895*fl[1]+0.6454972243679028*fl[0])*dvxl; 
  incr[12] = ((-3.0*fl[12])-1.732050807568877*fl[8])*wxl+((-0.4391550328268399*fl[18])-1.0*fl[11]-0.7745966692414834*fl[4]-0.447213595499958*fl[2])*dvxl; 
  incr[13] = (5.0*fl[13]+3.872983346207417*fl[5]+2.23606797749979*fl[3])*wxl+(1.118033988749895*fl[10]+0.6454972243679028*fl[6])*dvxl; 
  incr[14] = fl[14]*wxl+(0.447213595499958*fl[10]+0.2581988897471611*fl[6])*dvxl; 
  incr[15] = ((-3.0*fl[15])-1.732050807568877*fl[9])*wxl-0.5*fl[16]*dvxl; 
  incr[16] = fl[16]*wxl+(0.5*fl[15]+0.2886751345948129*fl[9])*dvxl; 
  incr[17] = ((-7.0*fl[17])-5.916079783099617*fl[7]-4.58257569495584*fl[1]-2.645751311064591*fl[0])*wxl+((-1.707825127659933*fl[11])-1.322875655532295*fl[4]-0.7637626158259735*fl[2])*dvxl; 
  incr[18] = fl[18]*wxl+(0.4391550328268399*fl[12]+0.253546276418555*fl[8])*dvxl; 
  incr[19] = fl[19]*wxl; 

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

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  outl[16] += -1.0*incr[16]*dxl; 
  outl[17] += incr[17]*dxl; 
  outl[18] += -1.0*incr[18]*dxl; 
  outl[19] += -1.0*incr[19]*dxl; 
  } else { 
  incr[0] = ((-2.645751311064591*fr[17])+2.23606797749979*fr[7]-1.732050807568877*fr[1]+fr[0])*wxr+(0.6454972243679028*fr[11]-0.5*fr[4]+0.2886751345948129*fr[2])*dvxr; 
  incr[1] = (4.58257569495584*fr[17]-3.872983346207417*fr[7]+3.0*fr[1]-1.732050807568877*fr[0])*wxr+((-1.118033988749895*fr[11])+0.8660254037844386*fr[4]-0.5*fr[2])*dvxr; 
  incr[2] = (2.23606797749979*fr[11]-1.732050807568877*fr[4]+fr[2])*wxr+((-0.7637626158259735*fr[17])-0.447213595499958*fr[12]+0.2581988897471612*fr[8]+0.6454972243679029*fr[7]-0.5*fr[1]+0.2886751345948129*fr[0])*dvxr; 
  incr[3] = (2.23606797749979*fr[13]-1.732050807568877*fr[5]+fr[3])*wxr+(0.2886751345948129*fr[6]-0.5*fr[10])*dvxr; 
  incr[4] = ((-3.872983346207417*fr[11])+3.0*fr[4]-1.732050807568877*fr[2])*wxr+(1.322875655532295*fr[17]+0.7745966692414834*fr[12]-0.4472135954999579*fr[8]-1.118033988749895*fr[7]+0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[5] = ((-3.872983346207417*fr[13])+3.0*fr[5]-1.732050807568877*fr[3])*wxr+(0.8660254037844386*fr[10]-0.5*fr[6])*dvxr; 
  incr[6] = (fr[6]-1.732050807568877*fr[10])*wxr+(0.2581988897471611*fr[14]+0.6454972243679028*fr[13]-0.5*fr[5]+0.2886751345948129*fr[3])*dvxr; 
  incr[7] = ((-5.916079783099617*fr[17])+5.0*fr[7]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*wxr+(1.443375672974065*fr[11]-1.118033988749895*fr[4]+0.6454972243679029*fr[2])*dvxr; 
  incr[8] = (fr[8]-1.732050807568877*fr[12])*wxr+(0.253546276418555*fr[18]+0.5773502691896257*fr[11]-0.4472135954999579*fr[4]+0.2581988897471612*fr[2])*dvxr; 
  incr[9] = (fr[9]-1.732050807568877*fr[15])*wxr+0.2886751345948129*fr[16]*dvxr; 
  incr[10] = (3.0*fr[10]-1.732050807568877*fr[6])*wxr+((-0.447213595499958*fr[14])-1.118033988749895*fr[13]+0.8660254037844386*fr[5]-0.5*fr[3])*dvxr; 
  incr[11] = (5.0*fr[11]-3.872983346207417*fr[4]+2.23606797749979*fr[2])*wxr+((-1.707825127659933*fr[17])-1.0*fr[12]+0.5773502691896257*fr[8]+1.443375672974065*fr[7]-1.118033988749895*fr[1]+0.6454972243679028*fr[0])*dvxr; 
  incr[12] = (3.0*fr[12]-1.732050807568877*fr[8])*wxr+((-0.4391550328268399*fr[18])-1.0*fr[11]+0.7745966692414834*fr[4]-0.447213595499958*fr[2])*dvxr; 
  incr[13] = (5.0*fr[13]-3.872983346207417*fr[5]+2.23606797749979*fr[3])*wxr+(0.6454972243679028*fr[6]-1.118033988749895*fr[10])*dvxr; 
  incr[14] = fr[14]*wxr+(0.2581988897471611*fr[6]-0.447213595499958*fr[10])*dvxr; 
  incr[15] = (3.0*fr[15]-1.732050807568877*fr[9])*wxr-0.5*fr[16]*dvxr; 
  incr[16] = fr[16]*wxr+(0.2886751345948129*fr[9]-0.5*fr[15])*dvxr; 
  incr[17] = (7.0*fr[17]-5.916079783099617*fr[7]+4.58257569495584*fr[1]-2.645751311064591*fr[0])*wxr+((-1.707825127659933*fr[11])+1.322875655532295*fr[4]-0.7637626158259735*fr[2])*dvxr; 
  incr[18] = fr[18]*wxr+(0.253546276418555*fr[8]-0.4391550328268399*fr[12])*dvxr; 
  incr[19] = fr[19]*wxr; 

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

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += incr[4]*dxl; 
  outl[5] += incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  outl[7] += -1.0*incr[7]*dxl; 
  outl[8] += -1.0*incr[8]*dxl; 
  outl[9] += -1.0*incr[9]*dxl; 
  outl[10] += incr[10]*dxl; 
  outl[11] += -1.0*incr[11]*dxl; 
  outl[12] += incr[12]*dxl; 
  outl[13] += -1.0*incr[13]*dxl; 
  outl[14] += -1.0*incr[14]*dxl; 
  outl[15] += incr[15]*dxl; 
  outl[16] += -1.0*incr[16]*dxl; 
  outl[17] += incr[17]*dxl; 
  outl[18] += -1.0*incr[18]*dxl; 
  outl[19] += -1.0*incr[19]*dxl; 
  } 
} 
