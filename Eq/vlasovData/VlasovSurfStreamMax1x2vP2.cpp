#include <VlasovModDecl.h> 
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
