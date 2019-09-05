#include <DistFuncMomentCalcModDecl.h> 
#include <cmath>
void GkM0Star1x1vSer_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[2]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 1.0*intFac*(wr[1]-wl[1]); 
 
  out[0] += ((-0.408248290463863*fr[2])+0.408248290463863*fl[2]+0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*dS; 
  out[1] += ((-0.408248290463863*fr[3])+0.408248290463863*fl[3]+0.3535533905932737*fr[1]+0.3535533905932737*fl[1])*dS; 
 
} 
 
void GkM0StarPositivity1x1vSer_VX(const double intFac, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *fl, const double *fr, double *out) 
{ 
  // intFac:  =2pi/m for gyrokinetics (not used in Vlasov). 
  // w[NDIM]: Cell-center coordinates. 
  // dxv[2]:  cell length in each direciton. 
  // fl/fr:   Distribution function in left/right cells 
  // out:     Increment to M_0^star from this cell surface. 
 
  const double dS = 1.0*intFac*(wr[1]-wl[1]); 
 
  if ( (0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]-0.1443375672974065*fr[1]-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) && ((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]+0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) ) {
    out[0] += ((-0.408248290463863*fr[2])+0.408248290463863*fl[2]+0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*dS; 
    out[1] += ((-0.408248290463863*fr[3])+0.408248290463863*fl[3]+0.3535533905932737*fr[1]+0.3535533905932737*fl[1])*dS; 
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

    if (std::abs(g1[0]) > 1.0e-15) {
      double g1Sq = g1[0]*g1[0];
      gBound[0] = (-(0.1443375672974065*g1[0]*fr[1])/std::sinh(g1[0]))-(0.1443375672974065*g1[0]*fl[1])/std::sinh(g1[0])+(0.25*fr[0]*g1[0])/std::sinh(g1[0])+(0.25*fl[0]*g1[0])/std::sinh(g1[0]); 
    } else {
      gBound[0] = (-0.1443375672974065*fr[1])-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    if (std::abs(g1[1]) > 1.0e-15) {
      double g1Sq = g1[1]*g1[1];
      gBound[1] = (0.1443375672974065*fr[1]*g1[1])/std::sinh(g1[1])+(0.1443375672974065*fl[1]*g1[1])/std::sinh(g1[1])+(0.25*fr[0]*g1[1])/std::sinh(g1[1])+(0.25*fl[0]*g1[1])/std::sinh(g1[1]); 
    } else {
      gBound[1] = 0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]; 
    };

    out[0] += (0.7071067811865475*gBound[1]+0.7071067811865475*gBound[0])*dS; 
    out[1] += (1.224744871391589*gBound[0]-1.224744871391589*gBound[1])*dS; 
  };
 
} 
 
void GkM1iM2Star1x1vSer(const double *w, const double *dxv, const double intFac, const double m_, const double *Bmag, const double *f, double *outM1i, double *outM2) 
{ 
  // w[2]:    Cell-center coordinates. 
  // dxv[2]:  Cell length in each direciton. 
  // intFac:  =2pi/m for gyrokinetics. 
  // m_:      mass. 
  // Bmag[2]: Magnetic field magnitude. 
  // f:       Distribution function. 
  // outM1i:  Contribution to M_1^star from this cell. 
  // outM2:   Contribution to M_2^star from this cell. 
 
  const double volFact = 0.5*dxv[1]; 
  double wvSq[1]; 
  wvSq[0]  = w[1]*w[1]; 
  double dvSq[1]; 
  dvSq[0] = dxv[1]*dxv[1]; 
 
  outM1i[0] += 1.414213562373095*f[0]*w[1]*volFact; 
  outM1i[1] += 1.414213562373095*f[1]*w[1]*volFact; 
 
  outM2[0] += (0.408248290463863*dxv[1]*w[1]*f[2]+1.414213562373095*f[0]*wvSq[0])*volFact; 
  outM2[1] += (0.408248290463863*dxv[1]*w[1]*f[3]+1.414213562373095*f[1]*wvSq[0])*volFact; 
 
} 
void GkBoundaryIntegral1x1vSer_F_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[4]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*intFac; 
 
  if (atLower) {
 
    out[0] += (1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x1vSer_F_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*intFac; 
 
  if (atLower) {
 
    out[0] += ((-1.58113883008419*fIn[5])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[1] += ((-1.58113883008419*fIn[7])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS; 
 
  } else {
 
    out[0] += (1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x1vSer_F_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[12]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*intFac; 
 
  if (atLower) {
 
    out[0] += (1.870828693386971*fIn[9]-1.58113883008419*fIn[5]+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.870828693386971*fIn[11]-1.58113883008419*fIn[7]+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS; 
    out[3] += (1.224744871391589*fIn[10]-0.7071067811865475*fIn[8])*dS; 
 
  } else {
 
    out[0] += (1.870828693386971*fIn[9]+1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS; 
    out[1] += (1.870828693386971*fIn[11]+1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS; 
    out[3] += (1.224744871391589*fIn[10]+0.7071067811865475*fIn[8])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x1vSer_vF_VX_P1(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[4]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*intFac; 
 
  if (atLower) {
 
    out[0] += (1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS*vBoundary+(0.6123724356957944*dxv[1]*fIn[2]-0.3535533905932737*fIn[0]*dxv[1])*dS; 
    out[1] += (1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS*vBoundary+(0.6123724356957944*dxv[1]*fIn[3]-0.3535533905932737*dxv[1]*fIn[1])*dS; 
 
  } else {
 
    out[0] += (1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS*vBoundary+((-0.6123724356957944*dxv[1]*fIn[2])-0.3535533905932737*fIn[0]*dxv[1])*dS; 
    out[1] += (1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS*vBoundary+((-0.6123724356957944*dxv[1]*fIn[3])-0.3535533905932737*dxv[1]*fIn[1])*dS; 
 
  }
 
} 
 
void GkBoundaryIntegral1x1vSer_vF_VX_P2(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[8]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*intFac; 
 
  if (atLower) {
 
    out[0] += ((-1.58113883008419*fIn[5])+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += ((-1.58113883008419*fIn[7])+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS*vBoundary; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += (1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS*vBoundary; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS*vBoundary; 
 
  }
 
} 
 
void GkBoundaryIntegral1x1vSer_vF_VX_P3(const bool atLower, const double intFac, const double vBoundary, const double *dxv, const double *fIn, double *out) 
{ 
  // atLower:   =true(false) if in cell at lower(upper) velocity boundary. 
  // intFac:    =2pi/m or 4pi/m for GkLBO (not used for Vlasov). 
  // vBoundary: velocity at the boundary of the velocity grid. 
  // dxv[2]:    cell length in each direciton. 
  // fIn[12]:    distribution function at velocity boundaries. 
  // out:       int dS of f|^(vmax)_(vmin) or vf^(vmax)_(vmin). 
 
  const double dS = 1.0*intFac; 
 
  if (atLower) {
 
    out[0] += (1.870828693386971*fIn[9]-1.58113883008419*fIn[5]+1.224744871391589*fIn[2]-0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += (1.870828693386971*fIn[11]-1.58113883008419*fIn[7]+1.224744871391589*fIn[3]-0.7071067811865475*fIn[1])*dS*vBoundary; 
    out[2] += (1.224744871391589*fIn[6]-0.7071067811865475*fIn[4])*dS*vBoundary; 
    out[3] += (1.224744871391589*fIn[10]-0.7071067811865475*fIn[8])*dS*vBoundary; 
 
  } else {
 
    out[0] += (1.870828693386971*fIn[9]+1.58113883008419*fIn[5]+1.224744871391589*fIn[2]+0.7071067811865475*fIn[0])*dS*vBoundary; 
    out[1] += (1.870828693386971*fIn[11]+1.58113883008419*fIn[7]+1.224744871391589*fIn[3]+0.7071067811865475*fIn[1])*dS*vBoundary; 
    out[2] += (1.224744871391589*fIn[6]+0.7071067811865475*fIn[4])*dS*vBoundary; 
    out[3] += (1.224744871391589*fIn[10]+0.7071067811865475*fIn[8])*dS*vBoundary; 
 
  }
 
} 
 
