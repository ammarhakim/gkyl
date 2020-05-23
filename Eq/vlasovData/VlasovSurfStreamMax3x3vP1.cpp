#include <VlasovModDecl.h> 
void VlasovSurfStream3x3vMax_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[3]; 
  double wxl = wl[3]; 

  double dvxr = dxvr[3]; 
  double wxr = wr[3]; 

  double dxl = 1.0/dxvl[0]; 
  double dxr = 1.0/dxvr[0]; 

  double incr[7]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[1]+fl[0])*wxl+0.2886751345948129*fl[4]*dvxl; 
  incr[1] = ((-3.0*fl[1])-1.732050807568877*fl[0])*wxl-0.5*fl[4]*dvxl; 
  incr[2] = fl[2]*wxl; 
  incr[3] = fl[3]*wxl; 
  incr[4] = fl[4]*wxl+(0.5*fl[1]+0.2886751345948129*fl[0])*dvxl; 
  incr[5] = fl[5]*wxl; 
  incr[6] = fl[6]*wxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  } else { 
  incr[0] = (fr[0]-1.732050807568877*fr[1])*wxr+0.2886751345948129*fr[4]*dvxr; 
  incr[1] = (3.0*fr[1]-1.732050807568877*fr[0])*wxr-0.5*fr[4]*dvxr; 
  incr[2] = fr[2]*wxr; 
  incr[3] = fr[3]*wxr; 
  incr[4] = fr[4]*wxr+(0.2886751345948129*fr[0]-0.5*fr[1])*dvxr; 
  incr[5] = fr[5]*wxr; 
  incr[6] = fr[6]*wxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  } 
} 
void VlasovSurfStream3x3vMax_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[4]; 
  double wxl = wl[4]; 

  double dvxr = dxvr[4]; 
  double wxr = wr[4]; 

  double dxl = 1.0/dxvl[1]; 
  double dxr = 1.0/dxvr[1]; 

  double incr[7]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[2]+fl[0])*wxl+0.2886751345948129*fl[5]*dvxl; 
  incr[1] = fl[1]*wxl; 
  incr[2] = ((-3.0*fl[2])-1.732050807568877*fl[0])*wxl-0.5*fl[5]*dvxl; 
  incr[3] = fl[3]*wxl; 
  incr[4] = fl[4]*wxl; 
  incr[5] = fl[5]*wxl+(0.5*fl[2]+0.2886751345948129*fl[0])*dvxl; 
  incr[6] = fl[6]*wxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  } else { 
  incr[0] = (fr[0]-1.732050807568877*fr[2])*wxr+0.2886751345948129*fr[5]*dvxr; 
  incr[1] = fr[1]*wxr; 
  incr[2] = (3.0*fr[2]-1.732050807568877*fr[0])*wxr-0.5*fr[5]*dvxr; 
  incr[3] = fr[3]*wxr; 
  incr[4] = fr[4]*wxr; 
  incr[5] = fr[5]*wxr+(0.2886751345948129*fr[0]-0.5*fr[2])*dvxr; 
  incr[6] = fr[6]*wxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += incr[2]*dxl; 
  outl[3] += -1.0*incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  } 
} 
void VlasovSurfStream3x3vMax_Z_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double dvxl = dxvl[5]; 
  double wxl = wl[5]; 

  double dvxr = dxvr[5]; 
  double wxr = wr[5]; 

  double dxl = 1.0/dxvl[2]; 
  double dxr = 1.0/dxvr[2]; 

  double incr[7]; 

  if (wxr>0) { 
  incr[0] = (1.732050807568877*fl[3]+fl[0])*wxl+0.2886751345948129*fl[6]*dvxl; 
  incr[1] = fl[1]*wxl; 
  incr[2] = fl[2]*wxl; 
  incr[3] = ((-3.0*fl[3])-1.732050807568877*fl[0])*wxl-0.5*fl[6]*dvxl; 
  incr[4] = fl[4]*wxl; 
  incr[5] = fl[5]*wxl; 
  incr[6] = fl[6]*wxl+(0.5*fl[3]+0.2886751345948129*fl[0])*dvxl; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  } else { 
  incr[0] = (fr[0]-1.732050807568877*fr[3])*wxr+0.2886751345948129*fr[6]*dvxr; 
  incr[1] = fr[1]*wxr; 
  incr[2] = fr[2]*wxr; 
  incr[3] = (3.0*fr[3]-1.732050807568877*fr[0])*wxr-0.5*fr[6]*dvxr; 
  incr[4] = fr[4]*wxr; 
  incr[5] = fr[5]*wxr; 
  incr[6] = fr[6]*wxr+(0.2886751345948129*fr[0]-0.5*fr[3])*dvxr; 

  outr[0] += incr[0]*dxr; 
  outr[1] += incr[1]*dxr; 
  outr[2] += incr[2]*dxr; 
  outr[3] += incr[3]*dxr; 
  outr[4] += incr[4]*dxr; 
  outr[5] += incr[5]*dxr; 
  outr[6] += incr[6]*dxr; 

  outl[0] += -1.0*incr[0]*dxl; 
  outl[1] += -1.0*incr[1]*dxl; 
  outl[2] += -1.0*incr[2]*dxl; 
  outl[3] += incr[3]*dxl; 
  outl[4] += -1.0*incr[4]*dxl; 
  outl[5] += -1.0*incr[5]*dxl; 
  outl[6] += -1.0*incr[6]*dxl; 
  } 
} 
