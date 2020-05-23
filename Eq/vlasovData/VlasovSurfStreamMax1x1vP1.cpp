#include <VlasovModDecl.h> 
void VlasovSurfStream1x1vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[1]; 
  double wxl = wl[1]; 

  double dvxr = dxvr[1]; 
  double wxr = wr[1]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[3]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[1]+fl[0])*wxl+0.2886751345948129*fl[2]*dvxl; 
  incr[1] = ((-3.0*fl[1])-1.732050807568877*fl[0])*wxl-0.5*fl[2]*dvxl; 
  incr[2] = fl[2]*wxl+(0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  } else { 
  incr[0] = (fr[0]-1.732050807568877*fr[1])*wxr+0.2886751345948129*fr[2]*dvxr; 
  incr[1] = (3.0*fr[1]-1.732050807568877*fr[0])*wxr-0.5*fr[2]*dvxr; 
  incr[2] = fr[2]*wxr+(0.2886751345948129*fr[0]-0.5*fr[1])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  } 
} 
