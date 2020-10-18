#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurfPositivity1xSerP1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[1]:      Cell-center coordinates.
  // dx[1]:     Cell spacing.
  // nu[1]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[2]; 
  double incr2[2]; 

  if ( ((-0.408248290463863*fr[1])+0.408248290463863*fl[1]+0.3535533905932737*fr[0]+0.3535533905932737*fl[0]>=0.0) ) {
    incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
    incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 

    incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 

    outr[0] += incr1[0]*rdxFnur; 
    outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 

    outl[0] += -1.0*incr1[0]*rdxFnul; 
    outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  } else {

    double xBar[1];
    xBar[0] = ((-0.53582588123382*fr[1])-0.53582588123382*fl[1]+0.662912607362388*fr[0]-0.662912607362388*fl[0])/(0.5*(0.6123724356957944*fr[1]-0.6123724356957944*fl[1])+3.0*((-0.408248290463863*fr[1])+0.408248290463863*fl[1]+0.3535533905932737*fr[0]+0.3535533905932737*fl[0])); 

    double xBarSq[1];
    xBarSq[0] = xBar[0]*xBar[0]; 

    double g1[1];
    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBarSq[0])-(1.0*xBar[0]*xBarSq[0])/(1.0-1.0*xBarSq[0]); 

    double gBound[1];
    double gBoundP[1];

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (-(0.3061862178478973*g1[0]*fr[1])/std::sinh(g1[0]))+(0.3061862178478973*g1[0]*fl[1])/std::sinh(g1[0])+(0.3535533905932737*fr[0]*g1[0])/std::sinh(g1[0])+(0.3535533905932737*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (-(0.3061862178478973*g1Sq*fr[1])/std::sinh(g1[0]))+(0.3061862178478973*g1Sq*fl[1])/std::sinh(g1[0])+(0.3535533905932737*fr[0]*g1Sq)/std::sinh(g1[0])+(0.3535533905932737*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = (-0.3061862178478973*fr[1])+0.3061862178478973*fl[1]+0.3535533905932737*fr[0]+0.3535533905932737*fl[0]; 
    };

    incr1[0] = -0.7071067811865475*gBoundP[0]; 
    incr1[1] = 1.224744871391589*gBoundP[0]; 

    incr2[1] = 1.224744871391589*gBound[0]; 

    outr[0] += incr1[0]*rdxFnur; 
    outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 

    outl[0] += -1.0*incr1[0]*rdxFnul; 
    outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
  };

} 
