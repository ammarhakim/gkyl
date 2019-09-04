#include <ConstDiffusionModDecl.h> 
void ConstDiffusionSurfPositivity2xSer_X_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[0]/(dxvl[0]*dxvl[0]); 
  double rdxSq2nur = 2.0*nu[0]/(dxvr[0]*dxvr[0]); 

  double incr1[4]; 
  double incr2[4]; 

  if ( (0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.1443375672974065*fr[2]-0.1443375672974065*fl[2]-0.2886751345948129*fr[1]+0.2886751345948129*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) && ((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]+0.1443375672974065*fr[2]+0.1443375672974065*fl[2]-0.2886751345948129*fr[1]+0.2886751345948129*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) ) {
    incr1[0] = 1.082531754730548*fr[1]+1.082531754730548*fl[1]-1.125*fr[0]+1.125*fl[0]; 
    incr1[1] = (-1.875*fr[1])-1.875*fl[1]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 
    incr1[2] = 1.082531754730548*fr[3]+1.082531754730548*fl[3]-1.125*fr[2]+1.125*fl[2]; 
    incr1[3] = (-1.875*fr[3])-1.875*fl[3]+1.948557158514986*fr[2]-1.948557158514986*fl[2]; 

    incr2[1] = (-1.0*fr[1])+fl[1]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
    incr2[3] = (-1.0*fr[3])+fl[3]+0.8660254037844386*fr[2]+0.8660254037844386*fl[2]; 

    outr[0] += incr1[0]*rdxSq2nur; 
    outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 
    outr[2] += incr1[2]*rdxSq2nur; 
    outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 

    outl[0] += -1.0*incr1[0]*rdxSq2nul; 
    outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 
    outl[2] += -1.0*incr1[2]*rdxSq2nul; 
    outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  } else {

    double xBar[2];
    xBar[0] = ((-0.125*fr[3])-0.125*fl[3]-0.2165063509461096*fr[2]+0.2165063509461096*fl[2]+0.2165063509461096*fr[1]+0.2165063509461096*fl[1]+0.375*fr[0]-0.375*fl[0])/(0.5*(1.732050807568877*(fr[1]-0.5773502691896258*fr[3])-1.732050807568877*(fl[1]-0.5773502691896258*fl[3]))-0.25*(3.464101615137754*(fr[1]-0.5773502691896258*fr[3])-3.464101615137754*(fl[1]-0.5773502691896258*fl[3])-3.0*(fr[0]-0.5773502691896258*fr[2])-3.0*(fl[0]-0.5773502691896258*fl[2]))); 
    xBar[1] = (0.125*fr[3]+0.125*fl[3]+0.2165063509461096*fr[2]-0.2165063509461096*fl[2]+0.2165063509461096*fr[1]+0.2165063509461096*fl[1]+0.375*fr[0]-0.375*fl[0])/(0.5*(1.732050807568877*(0.5773502691896258*fr[3]+fr[1])-1.732050807568877*(0.5773502691896258*fl[3]+fl[1]))-0.25*(3.464101615137754*(0.5773502691896258*fr[3]+fr[1])-3.464101615137754*(0.5773502691896258*fl[3]+fl[1])-3.0*(0.5773502691896258*fr[2]+fr[0])-3.0*(0.5773502691896258*fl[2]+fl[0]))); 

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
      gBound[0] = (-(0.1443375672974065*g1[0]*fr[2])/std::sinh(g1[0]))-(0.1443375672974065*g1[0]*fl[2])/std::sinh(g1[0])+(0.25*fr[0]*g1[0])/std::sinh(g1[0])+(0.25*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (-(0.1443375672974065*g1Sq*fr[2])/std::sinh(g1[0]))-(0.1443375672974065*g1Sq*fl[2])/std::sinh(g1[0])+(0.25*fr[0]*g1Sq)/std::sinh(g1[0])+(0.25*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = (-0.1443375672974065*fr[2])-0.1443375672974065*fl[2]+0.25*fr[0]+0.25*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (0.1443375672974065*g1[1]*fr[2])/std::sinh(g1[1])+(0.1443375672974065*g1[1]*fl[2])/std::sinh(g1[1])+(0.25*fr[0]*g1[1])/std::sinh(g1[1])+(0.25*fl[0]*g1[1])/std::sinh(g1[1]); 
      gBoundP[1] = (0.1443375672974065*g1Sq*fr[2])/std::sinh(g1[1])+(0.1443375672974065*g1Sq*fl[2])/std::sinh(g1[1])+(0.25*fr[0]*g1Sq)/std::sinh(g1[1])+(0.25*fl[0]*g1Sq)/std::sinh(g1[1]); 
    } else {
      gBound[1] = 0.1443375672974065*fr[2]+0.1443375672974065*fl[2]+0.25*fr[0]+0.25*fl[0]; 
    };

    incr1[0] = (-0.5*gBoundP[1])-0.5*gBoundP[0]; 
    incr1[1] = 0.8660254037844386*gBoundP[1]+0.8660254037844386*gBoundP[0]; 
    incr1[2] = 0.8660254037844386*gBoundP[1]-0.8660254037844386*gBoundP[0]; 
    incr1[3] = 1.5*gBoundP[0]-1.5*gBoundP[1]; 

    incr2[1] = 1.732050807568877*gBound[1]+1.732050807568877*gBound[0]; 
    incr2[3] = 3.0*gBound[0]-3.0*gBound[1]; 

    outr[0] += incr1[0]*rdxSq2nur; 
    outr[1] += incr2[1]*rdxSq2nur+incr1[1]*rdxSq2nur; 
    outr[2] += incr1[2]*rdxSq2nur; 
    outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 

    outl[0] += -1.0*incr1[0]*rdxSq2nul; 
    outl[1] += incr1[1]*rdxSq2nul-1.0*incr2[1]*rdxSq2nul; 
    outl[2] += -1.0*incr1[2]*rdxSq2nul; 
    outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  };

} 
void ConstDiffusionSurfPositivity2xSer_Y_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *nu, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:      Cell-center coordinates.
  // dxv[2]:    Cell spacing.
  // nu[2]:     diffusion coefficient (collisionality).
  // fl/fr:     Distribution function in left/right cells.
  // outl/outr: Incremented distribution function in left/right cells 
  double rdxSq2nul = 2.0*nu[1]/(dxvl[1]*dxvl[1]); 
  double rdxSq2nur = 2.0*nu[1]/(dxvr[1]*dxvr[1]); 

  double incr1[4]; 
  double incr2[4]; 

  if ( (0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]-0.1443375672974065*fr[1]-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) && ((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]+0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) ) {
    incr1[0] = 1.082531754730548*fr[2]+1.082531754730548*fl[2]-1.125*fr[0]+1.125*fl[0]; 
    incr1[1] = 1.082531754730548*fr[3]+1.082531754730548*fl[3]-1.125*fr[1]+1.125*fl[1]; 
    incr1[2] = (-1.875*fr[2])-1.875*fl[2]+1.948557158514986*fr[0]-1.948557158514986*fl[0]; 
    incr1[3] = (-1.875*fr[3])-1.875*fl[3]+1.948557158514986*fr[1]-1.948557158514986*fl[1]; 

    incr2[2] = (-1.0*fr[2])+fl[2]+0.8660254037844386*fr[0]+0.8660254037844386*fl[0]; 
    incr2[3] = (-1.0*fr[3])+fl[3]+0.8660254037844386*fr[1]+0.8660254037844386*fl[1]; 

    outr[0] += incr1[0]*rdxSq2nur; 
    outr[1] += incr1[1]*rdxSq2nur; 
    outr[2] += incr2[2]*rdxSq2nur+incr1[2]*rdxSq2nur; 
    outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 

    outl[0] += -1.0*incr1[0]*rdxSq2nul; 
    outl[1] += -1.0*incr1[1]*rdxSq2nul; 
    outl[2] += incr1[2]*rdxSq2nul-1.0*incr2[2]*rdxSq2nul; 
    outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  } else {

    double xBar[2];
    xBar[0] = ((-0.125*fr[3])-0.125*fl[3]+0.2165063509461096*fr[2]+0.2165063509461096*fl[2]-0.2165063509461096*fr[1]+0.2165063509461096*fl[1]+0.375*fr[0]-0.375*fl[0])/(0.5*(1.732050807568877*(fr[2]-0.5773502691896258*fr[3])-1.732050807568877*(fl[2]-0.5773502691896258*fl[3]))-0.25*(3.464101615137754*(fr[2]-0.5773502691896258*fr[3])-3.464101615137754*(fl[2]-0.5773502691896258*fl[3])-3.0*(fr[0]-0.5773502691896258*fr[1])-3.0*(fl[0]-0.5773502691896258*fl[1]))); 
    xBar[1] = (0.125*fr[3]+0.125*fl[3]+0.2165063509461096*fr[2]+0.2165063509461096*fl[2]+0.2165063509461096*fr[1]-0.2165063509461096*fl[1]+0.375*fr[0]-0.375*fl[0])/(0.5*(1.732050807568877*(0.5773502691896258*fr[3]+fr[2])-1.732050807568877*(0.5773502691896258*fl[3]+fl[2]))-0.25*(3.464101615137754*(0.5773502691896258*fr[3]+fr[2])-3.464101615137754*(0.5773502691896258*fl[3]+fl[2])-3.0*(0.5773502691896258*fr[1]+fr[0])-3.0*(0.5773502691896258*fl[1]+fl[0]))); 

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
      gBound[0] = (-(0.1443375672974065*g1[0]*fr[1])/std::sinh(g1[0]))-(0.1443375672974065*g1[0]*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1[0])/std::sinh(g1[0])+(0.25*fl[0]*g1[0])/std::sinh(g1[0]); 
      gBoundP[0] = (-(0.1443375672974065*g1Sq*fr[1])/std::sinh(g1[0]))-(0.1443375672974065*g1Sq*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1Sq)/std::sinh(g1[0])+(0.25*fl[0]*g1Sq)/std::sinh(g1[0]); 
    } else {
      gBound[0] = (-0.1443375672974065*fr[1])-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (0.1443375672974065*fr[1]*g1[1])/std::sinh(g1[1])+(0.1443375672974065*fl[1]*g1[1])/std::sinh(g1[1])+(0.25*fr[0]*g1[1])/std::sinh(g1[1])+(0.25*fl[0]*g1[1])/std::sinh(g1[1]); 
      gBoundP[1] = (0.1443375672974065*fr[1]*g1Sq)/std::sinh(g1[1])+(0.1443375672974065*fl[1]*g1Sq)/std::sinh(g1[1])+(0.25*fr[0]*g1Sq)/std::sinh(g1[1])+(0.25*fl[0]*g1Sq)/std::sinh(g1[1]); 
    } else {
      gBound[1] = 0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    incr1[0] = (-0.5*gBoundP[1])-0.5*gBoundP[0]; 
    incr1[1] = 0.8660254037844386*gBoundP[1]-0.8660254037844386*gBoundP[0]; 
    incr1[2] = 0.8660254037844386*gBoundP[1]+0.8660254037844386*gBoundP[0]; 
    incr1[3] = 1.5*gBoundP[0]-1.5*gBoundP[1]; 

    incr2[2] = 1.732050807568877*gBound[1]+1.732050807568877*gBound[0]; 
    incr2[3] = 3.0*gBound[0]-3.0*gBound[1]; 

    outr[0] += incr1[0]*rdxSq2nur; 
    outr[1] += incr1[1]*rdxSq2nur; 
    outr[2] += incr2[2]*rdxSq2nur+incr1[2]*rdxSq2nur; 
    outr[3] += incr2[3]*rdxSq2nur+incr1[3]*rdxSq2nur; 

    outl[0] += -1.0*incr1[0]*rdxSq2nul; 
    outl[1] += -1.0*incr1[1]*rdxSq2nul; 
    outl[2] += incr1[2]*rdxSq2nul-1.0*incr2[2]*rdxSq2nul; 
    outl[3] += incr1[3]*rdxSq2nul-1.0*incr2[3]*rdxSq2nul; 
  };

} 
