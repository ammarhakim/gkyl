#include <GkLBOModDecl.h> 
double GkLBOconstNuSurfPositivity1x1vSer_Vpar_P1(const double m_, const double *cflRateByDirL, const double *cflRateByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // m_:              species mass. 
  // cflRateByDir[2]: CFL rate in each direction. 
  // w[2]:            Cell-center coordinates. 
  // dxv[2]:          Cell spacing. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:       maximum midpoint value of v-u. 
  // nuUSum[2]:       sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:    sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:           Distribution function in left/right cells 
  // outl/outr:       Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[2]; 
  alphaDrSurf[0] = 0.7071067811865475*((2.0*wl[1]+dxvl[1])*nuSum-1.414213562373095*nuUSum[0]); 
  alphaDrSurf[1] = -1.0*nuUSum[1]; 

  double rCtrlL[2], rCtrlR[2];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rCtrlL[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rCtrlR[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fCtrlL[2], fCtrlR[2];  // fCtrl = anti-limited f evaluated at each control node on vx surface 
  // control node [x] = [-1/3] 
  fCtrlL[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rCtrlL[0],1.0); 
  fCtrlR[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rCtrlR[0],-1.0); 
  // control node [x] = [1/3] 
  fCtrlL[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rCtrlR[1],-1.0); 
  double fL_AL[2], fR_AL[2];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.7071067811865475*(fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 1.224744871391589*(fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.7071067811865475*(fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 1.224744871391589*(fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[2], fhatAL[2]; 
  alphaQuad = 0.7071067811865475*alphaDrSurf[1]-0.7071067811865475*alphaDrSurf[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = -0.5*(1.414213562373095*fL_AL[1]-1.414213562373095*fL_AL[0]); 
  } else {
  fhatALQuad[0] = -0.5*(1.414213562373095*fR_AL[1]-1.414213562373095*fR_AL[0]); 
  } 
  alphaQuad = -0.7071067811865475*(alphaDrSurf[1]+alphaDrSurf[0]); 
  if(alphaQuad > 0) {
  fhatALQuad[1] = 0.7071067811865476*(fL_AL[1]+fL_AL[0]); 
  } else {
  fhatALQuad[1] = 0.7071067811865476*(fR_AL[1]+fR_AL[0]); 
  } 
  fhatAL[0] = 0.7071067811865475*(fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.7071067811865475*(fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double Gdiff[4]; 
  double Ghat[4]; 
  double incr2[4]; 

  if ( (0.1666666666666667*fr[3]-0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]-0.1443375672974065*fr[1]-0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) && ((-0.1666666666666667*fr[3])+0.1666666666666667*fl[3]-0.2886751345948129*fr[2]+0.2886751345948129*fl[2]+0.1443375672974065*fr[1]+0.1443375672974065*fl[1]+0.25*fr[0]+0.25*fl[0]>=0.0) ) {
    incr2[2] = (nuVtSqSum[1]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])))*rdvSq4R; 
    incr2[3] = (nuVtSqSum[0]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])))*rdvSq4R; 

    Gdiff[0] = nuVtSqSum[1]*((-1.530931089239486*(fr[3]+fl[3]))+1.590990257669731*fr[1]-1.590990257669731*fl[1])+nuVtSqSum[0]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0]); 
    Gdiff[1] = nuVtSqSum[0]*((-1.530931089239486*(fr[3]+fl[3]))+1.590990257669731*fr[1]-1.590990257669731*fl[1])+nuVtSqSum[1]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0]); 

    Ghat[0] = Gdiff[0]*rdv+alphaDrSurf[1]*fhatAL[1]+alphaDrSurf[0]*fhatAL[0]; 
    Ghat[1] = Gdiff[1]*rdv+alphaDrSurf[0]*fhatAL[1]+fhatAL[0]*alphaDrSurf[1]; 
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

    incr2[2] = ((1.060660171779821*gBound[1]-1.060660171779821*gBound[0])*nuVtSqSum[1]+0.6123724356957944*nuVtSqSum[0]*(gBound[1]+gBound[0]))*rdvSq4R; 
    incr2[3] = (0.6123724356957944*(gBound[1]+gBound[0])*nuVtSqSum[1]+nuVtSqSum[0]*(1.060660171779821*gBound[1]-1.060660171779821*gBound[0]))*rdvSq4R; 

    Gdiff[0] = (1.224744871391589*gBoundP[1]-1.224744871391589*gBoundP[0])*nuVtSqSum[1]+0.7071067811865475*nuVtSqSum[0]*(gBoundP[1]+gBoundP[0]); 
    Gdiff[1] = 0.7071067811865475*(gBoundP[1]+gBoundP[0])*nuVtSqSum[1]+nuVtSqSum[0]*(1.224744871391589*gBoundP[1]-1.224744871391589*gBoundP[0]); 

    Ghat[0] = Gdiff[0]*rdv+alphaDrSurf[1]*fhatAL[1]+alphaDrSurf[0]*fhatAL[0]; 
    Ghat[1] = Gdiff[1]*rdv+alphaDrSurf[0]*fhatAL[1]+fhatAL[0]*alphaDrSurf[1]; 
  };

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]*rdv2R; 
  incr1[1] = -0.5*Ghat[1]*rdv2R; 
  incr1[2] = 0.8660254037844386*Ghat[0]*rdv2R; 
  incr1[3] = 0.8660254037844386*Ghat[1]*rdv2R; 

  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPosP[2], outrPosP[2]; 
  double outlPosM[2], outrPosM[2]; 
  fluxFracL = cflRateByDirL[0] == 0. ? 1.0 : cflRateByDirL[1]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] == 0. ? 1.0 : cflRateByDirR[1]/cflRateByDirR[0]; 
  outlPosP[0] = 0.1666666666666667*(incr2[3]-1.0*incr1[3]-1.732050807568877*incr2[2]+1.732050807568877*(incr1[2]+incr1[1])-3.0*incr1[0]); 
  outlPosP[1] = -0.1666666666666667*(incr2[3]-1.0*incr1[3]+1.732050807568877*incr2[2]-1.732050807568877*incr1[2]+1.732050807568877*incr1[1]+3.0*incr1[0]); 
  outlPosM[0] = -0.1666666666666667*(incr2[3]-1.0*incr1[3]-1.732050807568877*incr2[2]+1.732050807568877*incr1[2]-1.732050807568877*incr1[1]+3.0*incr1[0]); 
  outlPosM[1] = 0.1666666666666667*(incr2[3]-1.0*incr1[3]+1.732050807568877*incr2[2]-1.732050807568877*(incr1[2]+incr1[1])-3.0*incr1[0]); 
  outrPosP[0] = -0.1666666666666667*(incr2[3]+incr1[3]-1.732050807568877*(incr2[2]+incr1[2])+1.732050807568877*incr1[1]-3.0*incr1[0]); 
  outrPosP[1] = 0.1666666666666667*(incr2[3]+incr1[3]+1.732050807568877*(incr2[2]+incr1[2]+incr1[1])+3.0*incr1[0]); 
  outrPosM[0] = 0.1666666666666667*(incr2[3]+incr1[3]-1.732050807568877*(incr2[2]+incr1[2]+incr1[1])+3.0*incr1[0]); 
  outrPosM[1] = -0.1666666666666667*(incr2[3]+incr1[3]+1.732050807568877*(incr2[2]+incr1[2])-1.732050807568877*incr1[1]-3.0*incr1[0]); 
  if(outlPosP[0] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.1666666666666667*(fl[3]-1.732050807568877*fl[2]+1.732050807568877*fl[1]-3.0*fl[0]))/dtApprox/outlPosP[0]); 
  else limFac = 1.0; 
  if(outrPosM[0] < 0.) limFac = std::min(limFac, -fluxFracR*(0.1666666666666667*(fr[3]-1.732050807568877*(fr[2]+fr[1])+3.0*fr[0]))/dtApprox/outrPosM[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPosP[0] *= limFac; 
  outlPosM[0] *= limFac; 
  outrPosP[0] *= limFac; 
  outrPosM[0] *= limFac; 
  if(outlPosP[1] < 0.) limFac = std::min(1.0, -fluxFracL*(0.1666666666666667*(fl[3]+1.732050807568877*(fl[2]+fl[1])+3.0*fl[0]))/dtApprox/outlPosP[1]); 
  else limFac = 1.0; 
  if(outrPosM[1] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.1666666666666667*(fr[3]+1.732050807568877*fr[2]-1.732050807568877*fr[1]-3.0*fr[0]))/dtApprox/outrPosM[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPosP[1] *= limFac; 
  outlPosM[1] *= limFac; 
  outrPosP[1] *= limFac; 
  outrPosM[1] *= limFac; 
  outr[0] += 0.5*(outrPosP[1]+outrPosM[1]+outrPosP[0]+outrPosM[0]); 
  outr[1] += 0.5*(1.732050807568877*(outrPosP[1]+outrPosM[1])-1.732050807568877*(outrPosP[0]+outrPosM[0])); 
  outr[2] += 0.5*(1.732050807568877*outrPosP[1]-1.732050807568877*outrPosM[1]+1.732050807568877*outrPosP[0]-1.732050807568877*outrPosM[0]); 
  outr[3] += 0.5*(3.0*outrPosP[1]-3.0*(outrPosM[1]+outrPosP[0])+3.0*outrPosM[0]); 

  outl[0] += 0.5*(outlPosP[1]+outlPosM[1]+outlPosP[0]+outlPosM[0]); 
  outl[1] += 0.5*(1.732050807568877*(outlPosP[1]+outlPosM[1])-1.732050807568877*(outlPosP[0]+outlPosM[0])); 
  outl[2] += 0.5*(1.732050807568877*outlPosP[1]-1.732050807568877*outlPosM[1]+1.732050807568877*outlPosP[0]-1.732050807568877*outlPosM[0]); 
  outl[3] += 0.5*(3.0*outlPosP[1]-3.0*(outlPosM[1]+outlPosP[0])+3.0*outlPosM[0]); 

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
