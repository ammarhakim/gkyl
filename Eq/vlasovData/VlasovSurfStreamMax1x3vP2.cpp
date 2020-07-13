#include <VlasovModDecl.h> 
__host__ __device__ void VlasovSurfStream1x3vMax_X_P2(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[15]; 

  if (wxr>0) { 
  incr[0] = (2.23606797749979*fl[11]+1.732050807568877*fl[1]+fl[0])*wxl+(0.5*fl[5]+0.2886751345948129*fl[2])*dvxl; 
  incr[1] = ((-3.872983346207417*fl[11])-3.0*fl[1]-1.732050807568877*fl[0])*wxl+((-0.8660254037844386*fl[5])-0.5*fl[2])*dvxl; 
  incr[2] = (1.732050807568877*fl[5]+fl[2])*wxl+(0.2581988897471612*fl[12]+0.6454972243679029*fl[11]+0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[3] = (1.732050807568877*fl[6]+fl[3])*wxl+0.2886751345948129*fl[7]*dvxl; 
  incr[4] = (1.732050807568877*fl[8]+fl[4])*wxl+0.2886751345948129*fl[9]*dvxl; 
  incr[5] = ((-3.0*fl[5])-1.732050807568877*fl[2])*wxl+((-0.4472135954999579*fl[12])-1.118033988749895*fl[11]-0.8660254037844386*fl[1]-0.5*fl[0])*dvxl; 
  incr[6] = ((-3.0*fl[6])-1.732050807568877*fl[3])*wxl-0.5*fl[7]*dvxl; 
  incr[7] = fl[7]*wxl+(0.5*fl[6]+0.2886751345948129*fl[3])*dvxl; 
  incr[8] = ((-3.0*fl[8])-1.732050807568877*fl[4])*wxl-0.5*fl[9]*dvxl; 
  incr[9] = fl[9]*wxl+(0.5*fl[8]+0.2886751345948129*fl[4])*dvxl; 
  incr[10] = fl[10]*wxl; 
  incr[11] = (5.0*fl[11]+3.872983346207417*fl[1]+2.23606797749979*fl[0])*wxl+(1.118033988749895*fl[5]+0.6454972243679029*fl[2])*dvxl; 
  incr[12] = fl[12]*wxl+(0.4472135954999579*fl[5]+0.2581988897471612*fl[2])*dvxl; 
  incr[13] = fl[13]*wxl; 
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
  incr[0] = (2.23606797749979*fr[11]-1.732050807568877*fr[1]+fr[0])*wxr+(0.2886751345948129*fr[2]-0.5*fr[5])*dvxr; 
  incr[1] = ((-3.872983346207417*fr[11])+3.0*fr[1]-1.732050807568877*fr[0])*wxr+(0.8660254037844386*fr[5]-0.5*fr[2])*dvxr; 
  incr[2] = (fr[2]-1.732050807568877*fr[5])*wxr+(0.2581988897471612*fr[12]+0.6454972243679029*fr[11]-0.5*fr[1]+0.2886751345948129*fr[0])*dvxr; 
  incr[3] = (fr[3]-1.732050807568877*fr[6])*wxr+0.2886751345948129*fr[7]*dvxr; 
  incr[4] = (fr[4]-1.732050807568877*fr[8])*wxr+0.2886751345948129*fr[9]*dvxr; 
  incr[5] = (3.0*fr[5]-1.732050807568877*fr[2])*wxr+((-0.4472135954999579*fr[12])-1.118033988749895*fr[11]+0.8660254037844386*fr[1]-0.5*fr[0])*dvxr; 
  incr[6] = (3.0*fr[6]-1.732050807568877*fr[3])*wxr-0.5*fr[7]*dvxr; 
  incr[7] = fr[7]*wxr+(0.2886751345948129*fr[3]-0.5*fr[6])*dvxr; 
  incr[8] = (3.0*fr[8]-1.732050807568877*fr[4])*wxr-0.5*fr[9]*dvxr; 
  incr[9] = fr[9]*wxr+(0.2886751345948129*fr[4]-0.5*fr[8])*dvxr; 
  incr[10] = fr[10]*wxr; 
  incr[11] = (5.0*fr[11]-3.872983346207417*fr[1]+2.23606797749979*fr[0])*wxr+(0.6454972243679029*fr[2]-1.118033988749895*fr[5])*dvxr; 
  incr[12] = fr[12]*wxr+(0.2581988897471612*fr[2]-0.4472135954999579*fr[5])*dvxr; 
  incr[13] = fr[13]*wxr; 
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
