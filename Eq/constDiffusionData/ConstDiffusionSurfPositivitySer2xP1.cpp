#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurfPositivity2xSerP1_diffDirs1_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 

  if ( (0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.1443375672974065*fr[2]-0.1443375672974065*fl[2]-0.2886751345948129*fr[1]+0.2886751345948129*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) && ((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]+0.1443375672974065*fr[2]+0.1443375672974065*fl[2]-0.2886751345948129*fr[1]+0.2886751345948129*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) ) {
    incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
    incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
    incr1[2] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[2]+0.5625*fl[2]; 
    incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 

    incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
    incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 

    outr[0] += incr1[0]*rdxFnur; 
    outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
    outr[2] += incr1[2]*rdxFnur; 
    outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

    outl[0] += -1.0*incr1[0]*rdxFnul; 
    outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
    outl[2] += -1.0*incr1[2]*rdxFnul; 
    outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  } else {

    double xBar[2];
    xBar[0] = (0.21875*fr[3]+0.21875*fl[3]-0.270632938682637*fr[2]+0.270632938682637*fl[2]-0.3788861141556918*fr[1]-0.3788861141556918*fl[1]+0.46875*fr[0]-0.46875*fl[0])/(3.0*(0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.1443375672974065*fr[2]-0.1443375672974065*fl[2]-0.2886751345948129*fr[1]+0.2886751345948129*fl[1]+0.25*fr[0]+0.25*fl[0])+0.5*((-0.25*fr[3])+0.25*fl[3]+0.4330127018922193*fr[1]-0.4330127018922193*fl[1])); 
    xBar[1] = ((-0.21875*fr[3])-0.21875*fl[3]+0.270632938682637*fr[2]-0.270632938682637*fl[2]-0.3788861141556918*fr[1]-0.3788861141556918*fl[1]+0.46875*fr[0]-0.46875*fl[0])/(0.5*(0.25*fr[3]-0.25*fl[3]+0.4330127018922193*fr[1]-0.4330127018922193*fl[1])+3.0*((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]+0.1443375672974065*fr[2]+0.1443375672974065*fl[2]-0.2886751345948129*fr[1]+0.2886751345948129*fl[1]+0.25*fr[0]+0.25*fl[0])); 

    double xBarSq[2];
    xBarSq[0] = xBar[0]*xBar[0]; 
    xBarSq[1] = xBar[1]*xBar[1]; 

    double g1[2];
    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBarSq[0])-(1.0*xBar[0]*xBarSq[0])/(1.0-1.0*xBarSq[0]); 
    g1[1] = (3.0*xBar[1])/(1.0-1.0*xBarSq[1])-(1.0*xBar[1]*xBarSq[1])/(1.0-1.0*xBarSq[1]); 

    double gBound[2];
    double gBoundP[2];

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (0.125*g1[0]*fr[3])/std::sinh(g1[0])-(0.125*g1[0]*fl[3])/std::sinh(g1[0])-(0.1443375672974065*g1[0]*fr[2])/std::sinh(g1[0])-(0.1443375672974065*g1[0]*fl[2])/std::sinh(g1[0])-(0.2165063509461097*g1[0]*fr[1])/std::sinh(g1[0])+(0.2165063509461097*g1[0]*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1[0])/std::sinh(g1[0])+(0.25*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (0.125*g1Sq*fr[3])/std::sinh(g1[0])-(0.125*g1Sq*fl[3])/std::sinh(g1[0])-(0.1443375672974065*g1Sq*fr[2])/std::sinh(g1[0])-(0.1443375672974065*g1Sq*fl[2])/std::sinh(g1[0])-(0.2165063509461097*g1Sq*fr[1])/std::sinh(g1[0])+(0.2165063509461097*g1Sq*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1Sq)/std::sinh(g1[0])+(0.25*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = 0.125*fr[3]-0.125*fl[3]-0.1443375672974065*fr[2]-0.1443375672974065*fl[2]-0.2165063509461097*fr[1]+0.2165063509461097*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (-(0.125*g1[1]*fr[3])/std::sinh(g1[1]))+(0.125*g1[1]*fl[3])/std::sinh(g1[1])+(0.1443375672974065*g1[1]*fr[2])/std::sinh(g1[1])+(0.1443375672974065*g1[1]*fl[2])/std::sinh(g1[1])-(0.2165063509461097*fr[1]*g1[1])/std::sinh(g1[1])+(0.2165063509461097*fl[1]*g1[1])/std::sinh(g1[1])+(0.25*fr[0]*g1[1])/std::sinh(g1[1])+(0.25*fl[0]*g1[1])/std::sinh(g1[1]); 
      gBoundP[1] = (-(0.125*g1Sq*fr[3])/std::sinh(g1[1]))+(0.125*g1Sq*fl[3])/std::sinh(g1[1])+(0.1443375672974065*g1Sq*fr[2])/std::sinh(g1[1])+(0.1443375672974065*g1Sq*fl[2])/std::sinh(g1[1])-(0.2165063509461097*fr[1]*g1Sq)/std::sinh(g1[1])+(0.2165063509461097*fl[1]*g1Sq)/std::sinh(g1[1])+(0.25*fr[0]*g1Sq)/std::sinh(g1[1])+(0.25*fl[0]*g1Sq)/std::sinh(g1[1]); 
    } else {
      gBound[1] = (-0.125*fr[3])+0.125*fl[3]+0.1443375672974065*fr[2]+0.1443375672974065*fl[2]-0.2165063509461097*fr[1]+0.2165063509461097*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    incr1[0] = (-0.5*gBoundP[1])-0.5*gBoundP[0]; 
    incr1[1] = 0.8660254037844386*gBoundP[1]+0.8660254037844386*gBoundP[0]; 
    incr1[2] = 0.8660254037844386*gBoundP[0]-0.8660254037844386*gBoundP[1]; 
    incr1[3] = 1.5*gBoundP[1]-1.5*gBoundP[0]; 

    incr2[1] = 0.8660254037844386*gBound[1]+0.8660254037844386*gBound[0]; 
    incr2[3] = 1.5*gBound[1]-1.5*gBound[0]; 

    outr[0] += incr1[0]*rdxFnur; 
    outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
    outr[2] += incr1[2]*rdxFnur; 
    outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

    outl[0] += -1.0*incr1[0]*rdxFnul; 
    outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
    outl[2] += -1.0*incr1[2]*rdxFnul; 
    outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  };

} 
void ConstDiffusionSurfPositivity2xSerP1_diffDirs12_X1(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[0]*dxl[0]); 
  double rdxFnur = 4.0*nu[0]/(dxr[0]*dxr[0]); 

  double incr1[4]; 
  double incr2[4]; 

  if ( (0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.1443375672974065*fr[2]-0.1443375672974065*fl[2]-0.2886751345948129*fr[1]+0.2886751345948129*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) && ((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]+0.1443375672974065*fr[2]+0.1443375672974065*fl[2]-0.2886751345948129*fr[1]+0.2886751345948129*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) ) {
    incr1[0] = 0.5412658773652741*fr[1]+0.5412658773652741*fl[1]-0.5625*fr[0]+0.5625*fl[0]; 
    incr1[1] = (-0.9375*fr[1])-0.9375*fl[1]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
    incr1[2] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[2]+0.5625*fl[2]; 
    incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[2]-0.9742785792574932*fl[2]; 

    incr2[1] = (-0.5*fr[1])+0.5*fl[1]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
    incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[2]+0.4330127018922193*fl[2]; 

    outr[0] += incr1[0]*rdxFnur; 
    outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
    outr[2] += incr1[2]*rdxFnur; 
    outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

    outl[0] += -1.0*incr1[0]*rdxFnul; 
    outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
    outl[2] += -1.0*incr1[2]*rdxFnul; 
    outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  } else {

    double xBar[2];
    xBar[0] = (0.21875*fr[3]+0.21875*fl[3]-0.270632938682637*fr[2]+0.270632938682637*fl[2]-0.3788861141556918*fr[1]-0.3788861141556918*fl[1]+0.46875*fr[0]-0.46875*fl[0])/(3.0*(0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.1443375672974065*fr[2]-0.1443375672974065*fl[2]-0.2886751345948129*fr[1]+0.2886751345948129*fl[1]+0.25*fr[0]+0.25*fl[0])+0.5*((-0.25*fr[3])+0.25*fl[3]+0.4330127018922193*fr[1]-0.4330127018922193*fl[1])); 
    xBar[1] = ((-0.21875*fr[3])-0.21875*fl[3]+0.270632938682637*fr[2]-0.270632938682637*fl[2]-0.3788861141556918*fr[1]-0.3788861141556918*fl[1]+0.46875*fr[0]-0.46875*fl[0])/(0.5*(0.25*fr[3]-0.25*fl[3]+0.4330127018922193*fr[1]-0.4330127018922193*fl[1])+3.0*((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]+0.1443375672974065*fr[2]+0.1443375672974065*fl[2]-0.2886751345948129*fr[1]+0.2886751345948129*fl[1]+0.25*fr[0]+0.25*fl[0])); 

    double xBarSq[2];
    xBarSq[0] = xBar[0]*xBar[0]; 
    xBarSq[1] = xBar[1]*xBar[1]; 

    double g1[2];
    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBarSq[0])-(1.0*xBar[0]*xBarSq[0])/(1.0-1.0*xBarSq[0]); 
    g1[1] = (3.0*xBar[1])/(1.0-1.0*xBarSq[1])-(1.0*xBar[1]*xBarSq[1])/(1.0-1.0*xBarSq[1]); 

    double gBound[2];
    double gBoundP[2];

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (0.125*g1[0]*fr[3])/std::sinh(g1[0])-(0.125*g1[0]*fl[3])/std::sinh(g1[0])-(0.1443375672974065*g1[0]*fr[2])/std::sinh(g1[0])-(0.1443375672974065*g1[0]*fl[2])/std::sinh(g1[0])-(0.2165063509461097*g1[0]*fr[1])/std::sinh(g1[0])+(0.2165063509461097*g1[0]*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1[0])/std::sinh(g1[0])+(0.25*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (0.125*g1Sq*fr[3])/std::sinh(g1[0])-(0.125*g1Sq*fl[3])/std::sinh(g1[0])-(0.1443375672974065*g1Sq*fr[2])/std::sinh(g1[0])-(0.1443375672974065*g1Sq*fl[2])/std::sinh(g1[0])-(0.2165063509461097*g1Sq*fr[1])/std::sinh(g1[0])+(0.2165063509461097*g1Sq*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1Sq)/std::sinh(g1[0])+(0.25*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = 0.125*fr[3]-0.125*fl[3]-0.1443375672974065*fr[2]-0.1443375672974065*fl[2]-0.2165063509461097*fr[1]+0.2165063509461097*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (-(0.125*g1[1]*fr[3])/std::sinh(g1[1]))+(0.125*g1[1]*fl[3])/std::sinh(g1[1])+(0.1443375672974065*g1[1]*fr[2])/std::sinh(g1[1])+(0.1443375672974065*g1[1]*fl[2])/std::sinh(g1[1])-(0.2165063509461097*fr[1]*g1[1])/std::sinh(g1[1])+(0.2165063509461097*fl[1]*g1[1])/std::sinh(g1[1])+(0.25*fr[0]*g1[1])/std::sinh(g1[1])+(0.25*fl[0]*g1[1])/std::sinh(g1[1]); 
      gBoundP[1] = (-(0.125*g1Sq*fr[3])/std::sinh(g1[1]))+(0.125*g1Sq*fl[3])/std::sinh(g1[1])+(0.1443375672974065*g1Sq*fr[2])/std::sinh(g1[1])+(0.1443375672974065*g1Sq*fl[2])/std::sinh(g1[1])-(0.2165063509461097*fr[1]*g1Sq)/std::sinh(g1[1])+(0.2165063509461097*fl[1]*g1Sq)/std::sinh(g1[1])+(0.25*fr[0]*g1Sq)/std::sinh(g1[1])+(0.25*fl[0]*g1Sq)/std::sinh(g1[1]); 
    } else {
      gBound[1] = (-0.125*fr[3])+0.125*fl[3]+0.1443375672974065*fr[2]+0.1443375672974065*fl[2]-0.2165063509461097*fr[1]+0.2165063509461097*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    incr1[0] = (-0.5*gBoundP[1])-0.5*gBoundP[0]; 
    incr1[1] = 0.8660254037844386*gBoundP[1]+0.8660254037844386*gBoundP[0]; 
    incr1[2] = 0.8660254037844386*gBoundP[0]-0.8660254037844386*gBoundP[1]; 
    incr1[3] = 1.5*gBoundP[1]-1.5*gBoundP[0]; 

    incr2[1] = 0.8660254037844386*gBound[1]+0.8660254037844386*gBound[0]; 
    incr2[3] = 1.5*gBound[1]-1.5*gBound[0]; 

    outr[0] += incr1[0]*rdxFnur; 
    outr[1] += incr2[1]*rdxFnur+incr1[1]*rdxFnur; 
    outr[2] += incr1[2]*rdxFnur; 
    outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

    outl[0] += -1.0*incr1[0]*rdxFnul; 
    outl[1] += incr1[1]*rdxFnul-1.0*incr2[1]*rdxFnul; 
    outl[2] += -1.0*incr1[2]*rdxFnul; 
    outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  };

} 
void ConstDiffusionSurfPositivity2xSerP1_diffDirs12_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[1]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[1]/(dxr[1]*dxr[1]); 

  double incr1[4]; 
  double incr2[4]; 

  if ( (0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]-0.1443375672974065*fr[1]-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) && ((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]+0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) ) {
    incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
    incr1[1] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[1]+0.5625*fl[1]; 
    incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
    incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 

    incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
    incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 

    outr[0] += incr1[0]*rdxFnur; 
    outr[1] += incr1[1]*rdxFnur; 
    outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
    outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

    outl[0] += -1.0*incr1[0]*rdxFnul; 
    outl[1] += -1.0*incr1[1]*rdxFnul; 
    outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
    outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  } else {

    double xBar[2];
    xBar[0] = (0.21875*fr[3]+0.21875*fl[3]-0.3788861141556918*fr[2]-0.3788861141556918*fl[2]-0.270632938682637*fr[1]+0.270632938682637*fl[1]+0.46875*fr[0]-0.46875*fl[0])/(3.0*(0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]-0.1443375672974065*fr[1]-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0])+0.5*((-0.25*fr[3])+0.25*fl[3]+0.4330127018922193*fr[2]-0.4330127018922193*fl[2])); 
    xBar[1] = ((-0.21875*fr[3])-0.21875*fl[3]-0.3788861141556918*fr[2]-0.3788861141556918*fl[2]+0.270632938682637*fr[1]-0.270632938682637*fl[1]+0.46875*fr[0]-0.46875*fl[0])/(0.5*(0.25*fr[3]-0.25*fl[3]+0.4330127018922193*fr[2]-0.4330127018922193*fl[2])+3.0*((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]+0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0])); 

    double xBarSq[2];
    xBarSq[0] = xBar[0]*xBar[0]; 
    xBarSq[1] = xBar[1]*xBar[1]; 

    double g1[2];
    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBarSq[0])-(1.0*xBar[0]*xBarSq[0])/(1.0-1.0*xBarSq[0]); 
    g1[1] = (3.0*xBar[1])/(1.0-1.0*xBarSq[1])-(1.0*xBar[1]*xBarSq[1])/(1.0-1.0*xBarSq[1]); 

    double gBound[2];
    double gBoundP[2];

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (0.125*g1[0]*fr[3])/std::sinh(g1[0])-(0.125*g1[0]*fl[3])/std::sinh(g1[0])-(0.2165063509461097*g1[0]*fr[2])/std::sinh(g1[0])+(0.2165063509461097*g1[0]*fl[2])/std::sinh(g1[0])-(0.1443375672974065*g1[0]*fr[1])/std::sinh(g1[0])-(0.1443375672974065*g1[0]*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1[0])/std::sinh(g1[0])+(0.25*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (0.125*g1Sq*fr[3])/std::sinh(g1[0])-(0.125*g1Sq*fl[3])/std::sinh(g1[0])-(0.2165063509461097*g1Sq*fr[2])/std::sinh(g1[0])+(0.2165063509461097*g1Sq*fl[2])/std::sinh(g1[0])-(0.1443375672974065*g1Sq*fr[1])/std::sinh(g1[0])-(0.1443375672974065*g1Sq*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1Sq)/std::sinh(g1[0])+(0.25*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = 0.125*fr[3]-0.125*fl[3]-0.2165063509461097*fr[2]+0.2165063509461097*fl[2]-0.1443375672974065*fr[1]-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (-(0.125*g1[1]*fr[3])/std::sinh(g1[1]))+(0.125*g1[1]*fl[3])/std::sinh(g1[1])-(0.2165063509461097*g1[1]*fr[2])/std::sinh(g1[1])+(0.2165063509461097*g1[1]*fl[2])/std::sinh(g1[1])+(0.1443375672974065*fr[1]*g1[1])/std::sinh(g1[1])+(0.1443375672974065*fl[1]*g1[1])/std::sinh(g1[1])+(0.25*fr[0]*g1[1])/std::sinh(g1[1])+(0.25*fl[0]*g1[1])/std::sinh(g1[1]); 
      gBoundP[1] = (-(0.125*g1Sq*fr[3])/std::sinh(g1[1]))+(0.125*g1Sq*fl[3])/std::sinh(g1[1])-(0.2165063509461097*g1Sq*fr[2])/std::sinh(g1[1])+(0.2165063509461097*g1Sq*fl[2])/std::sinh(g1[1])+(0.1443375672974065*fr[1]*g1Sq)/std::sinh(g1[1])+(0.1443375672974065*fl[1]*g1Sq)/std::sinh(g1[1])+(0.25*fr[0]*g1Sq)/std::sinh(g1[1])+(0.25*fl[0]*g1Sq)/std::sinh(g1[1]); 
    } else {
      gBound[1] = (-0.125*fr[3])+0.125*fl[3]-0.2165063509461097*fr[2]+0.2165063509461097*fl[2]+0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    incr1[0] = (-0.5*gBoundP[1])-0.5*gBoundP[0]; 
    incr1[1] = 0.8660254037844386*gBoundP[0]-0.8660254037844386*gBoundP[1]; 
    incr1[2] = 0.8660254037844386*gBoundP[1]+0.8660254037844386*gBoundP[0]; 
    incr1[3] = 1.5*gBoundP[1]-1.5*gBoundP[0]; 

    incr2[2] = 0.8660254037844386*gBound[1]+0.8660254037844386*gBound[0]; 
    incr2[3] = 1.5*gBound[1]-1.5*gBound[0]; 

    outr[0] += incr1[0]*rdxFnur; 
    outr[1] += incr1[1]*rdxFnur; 
    outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
    outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

    outl[0] += -1.0*incr1[0]*rdxFnul; 
    outl[1] += -1.0*incr1[1]*rdxFnul; 
    outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
    outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  };

} 
void ConstDiffusionSurfPositivity2xSerP1_diffDirs2_X2(const double *wl, const double *wr, const double *dxl, const double *dxr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dx[2]:     Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxFnul = 4.0*nu[0]/(dxl[1]*dxl[1]); 
  double rdxFnur = 4.0*nu[0]/(dxr[1]*dxr[1]); 

  double incr1[4]; 
  double incr2[4]; 

  if ( (0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]-0.1443375672974065*fr[1]-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) && ((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]+0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) ) {
    incr1[0] = 0.5412658773652741*fr[2]+0.5412658773652741*fl[2]-0.5625*fr[0]+0.5625*fl[0]; 
    incr1[1] = 0.5412658773652741*fr[3]+0.5412658773652741*fl[3]-0.5625*fr[1]+0.5625*fl[1]; 
    incr1[2] = (-0.9375*fr[2])-0.9375*fl[2]+0.9742785792574932*fr[0]-0.9742785792574932*fl[0]; 
    incr1[3] = (-0.9375*fr[3])-0.9375*fl[3]+0.9742785792574932*fr[1]-0.9742785792574932*fl[1]; 

    incr2[2] = (-0.5*fr[2])+0.5*fl[2]+0.4330127018922193*fr[0]+0.4330127018922193*fl[0]; 
    incr2[3] = (-0.5*fr[3])+0.5*fl[3]+0.4330127018922193*fr[1]+0.4330127018922193*fl[1]; 

    outr[0] += incr1[0]*rdxFnur; 
    outr[1] += incr1[1]*rdxFnur; 
    outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
    outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

    outl[0] += -1.0*incr1[0]*rdxFnul; 
    outl[1] += -1.0*incr1[1]*rdxFnul; 
    outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
    outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  } else {

    double xBar[2];
    xBar[0] = (0.21875*fr[3]+0.21875*fl[3]-0.3788861141556918*fr[2]-0.3788861141556918*fl[2]-0.270632938682637*fr[1]+0.270632938682637*fl[1]+0.46875*fr[0]-0.46875*fl[0])/(3.0*(0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]-0.1443375672974065*fr[1]-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0])+0.5*((-0.25*fr[3])+0.25*fl[3]+0.4330127018922193*fr[2]-0.4330127018922193*fl[2])); 
    xBar[1] = ((-0.21875*fr[3])-0.21875*fl[3]-0.3788861141556918*fr[2]-0.3788861141556918*fl[2]+0.270632938682637*fr[1]-0.270632938682637*fl[1]+0.46875*fr[0]-0.46875*fl[0])/(0.5*(0.25*fr[3]-0.25*fl[3]+0.4330127018922193*fr[2]-0.4330127018922193*fl[2])+3.0*((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]+0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0])); 

    double xBarSq[2];
    xBarSq[0] = xBar[0]*xBar[0]; 
    xBarSq[1] = xBar[1]*xBar[1]; 

    double g1[2];
    g1[0] = (3.0*xBar[0])/(1.0-1.0*xBarSq[0])-(1.0*xBar[0]*xBarSq[0])/(1.0-1.0*xBarSq[0]); 
    g1[1] = (3.0*xBar[1])/(1.0-1.0*xBarSq[1])-(1.0*xBar[1]*xBarSq[1])/(1.0-1.0*xBarSq[1]); 

    double gBound[2];
    double gBoundP[2];

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (0.125*g1[0]*fr[3])/std::sinh(g1[0])-(0.125*g1[0]*fl[3])/std::sinh(g1[0])-(0.2165063509461097*g1[0]*fr[2])/std::sinh(g1[0])+(0.2165063509461097*g1[0]*fl[2])/std::sinh(g1[0])-(0.1443375672974065*g1[0]*fr[1])/std::sinh(g1[0])-(0.1443375672974065*g1[0]*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1[0])/std::sinh(g1[0])+(0.25*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (0.125*g1Sq*fr[3])/std::sinh(g1[0])-(0.125*g1Sq*fl[3])/std::sinh(g1[0])-(0.2165063509461097*g1Sq*fr[2])/std::sinh(g1[0])+(0.2165063509461097*g1Sq*fl[2])/std::sinh(g1[0])-(0.1443375672974065*g1Sq*fr[1])/std::sinh(g1[0])-(0.1443375672974065*g1Sq*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1Sq)/std::sinh(g1[0])+(0.25*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = 0.125*fr[3]-0.125*fl[3]-0.2165063509461097*fr[2]+0.2165063509461097*fl[2]-0.1443375672974065*fr[1]-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (-(0.125*g1[1]*fr[3])/std::sinh(g1[1]))+(0.125*g1[1]*fl[3])/std::sinh(g1[1])-(0.2165063509461097*g1[1]*fr[2])/std::sinh(g1[1])+(0.2165063509461097*g1[1]*fl[2])/std::sinh(g1[1])+(0.1443375672974065*fr[1]*g1[1])/std::sinh(g1[1])+(0.1443375672974065*fl[1]*g1[1])/std::sinh(g1[1])+(0.25*fr[0]*g1[1])/std::sinh(g1[1])+(0.25*fl[0]*g1[1])/std::sinh(g1[1]); 
      gBoundP[1] = (-(0.125*g1Sq*fr[3])/std::sinh(g1[1]))+(0.125*g1Sq*fl[3])/std::sinh(g1[1])-(0.2165063509461097*g1Sq*fr[2])/std::sinh(g1[1])+(0.2165063509461097*g1Sq*fl[2])/std::sinh(g1[1])+(0.1443375672974065*fr[1]*g1Sq)/std::sinh(g1[1])+(0.1443375672974065*fl[1]*g1Sq)/std::sinh(g1[1])+(0.25*fr[0]*g1Sq)/std::sinh(g1[1])+(0.25*fl[0]*g1Sq)/std::sinh(g1[1]); 
    } else {
      gBound[1] = (-0.125*fr[3])+0.125*fl[3]-0.2165063509461097*fr[2]+0.2165063509461097*fl[2]+0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    incr1[0] = (-0.5*gBoundP[1])-0.5*gBoundP[0]; 
    incr1[1] = 0.8660254037844386*gBoundP[0]-0.8660254037844386*gBoundP[1]; 
    incr1[2] = 0.8660254037844386*gBoundP[1]+0.8660254037844386*gBoundP[0]; 
    incr1[3] = 1.5*gBoundP[1]-1.5*gBoundP[0]; 

    incr2[2] = 0.8660254037844386*gBound[1]+0.8660254037844386*gBound[0]; 
    incr2[3] = 1.5*gBound[1]-1.5*gBound[0]; 

    outr[0] += incr1[0]*rdxFnur; 
    outr[1] += incr1[1]*rdxFnur; 
    outr[2] += incr2[2]*rdxFnur+incr1[2]*rdxFnur; 
    outr[3] += incr2[3]*rdxFnur+incr1[3]*rdxFnur; 

    outl[0] += -1.0*incr1[0]*rdxFnul; 
    outl[1] += -1.0*incr1[1]*rdxFnul; 
    outl[2] += incr1[2]*rdxFnul-1.0*incr2[2]*rdxFnul; 
    outl[3] += incr1[3]*rdxFnul-1.0*incr2[3]*rdxFnul; 
  };

} 
