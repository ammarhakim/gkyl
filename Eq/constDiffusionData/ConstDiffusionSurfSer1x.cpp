#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurf1xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. nu[NDIM]: diffusion coefficient (collisionality). fl/fr: Distribution function in left/right cells 
// outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[2]; 
  double incr2[2]; 
  if ((-0.408248290463863*fr[1])+0.408248290463863*fl[1]+0.3535533905932737*fr[0]+0.3535533905932737*fl[0] >= 0.0) {
  incr1[0] = 1.082531754730548*fr[1]+1.082531754730548*fl[1]-1.125*fr[0]+1.125*fl[0]; 
  incr1[1] = (-1.875*fr[1])-1.875*fl[1]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 

  incr2[1] = (-1.0*fr[1])+fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 

  } else {
  double xBar = (0.3061862178478973*fr[1]+0.3061862178478973*fl[1]+0.5303300858899103*fr[0]-0.5303300858899103*fl[0])/(0.5*(2.449489742783178*fr[1]-2.449489742783178*fl[1])-0.25*(4.898979485566357*fr[1]-4.898979485566357*fl[1]-4.242640687119286*fr[0]-4.242640687119286*fl[0]));
  double xBarSq = xBar*xBar;
  double g1 = (3.0*xBar)/(1.0-1.0*xBarSq)-(1.0*xBar*xBarSq)/(1.0-1.0*xBarSq);
  if (std::abs(g1) > 1.0e-15) {
  incr1[0] = (-(0.25*fr[0]*g1*g1)/std::sinh(g1))-(0.25*fl[0]*g1*g1)/std::sinh(g1); 
  incr1[1] = (0.4330127018922193*fr[0]*g1*g1)/std::sinh(g1)+(0.4330127018922193*fl[0]*g1*g1)/std::sinh(g1); 

  incr2[1] = (0.8660254037844386*fr[0]*g1)/std::sinh(g1)+(0.8660254037844386*fl[0]*g1)/std::sinh(g1); 

  } else {

  incr2[1] = (0.8660254037844386*fr[0]*g1)/std::sinh(g1)+(0.8660254037844386*fl[0]*g1)/std::sinh(g1)+fr[1]-1.0*fl[1]; 

  };

  };

  outr[0] += incr1[0]*rdxSq2nur; 
  outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 

  outl[0] += -1.0*incr1[0]*rdxSq2nul; 
  outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 

} 
