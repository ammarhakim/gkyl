#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf1x1vSer_Vpar_P1(const double m_, const double *cflRateByDirL, const double *cflRateByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
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

  double fUpwindQuad[2];
  fUpwindQuad[0] = 0.5*(0.8660254037844386*fr[3]-0.8660254037844386*(fl[3]+fr[2])+0.8660254037844386*fl[2]-0.5*(fr[1]+fl[1])+0.5*(fr[0]+fl[0]))-0.5*((-0.8660254037844386*(fr[3]+fl[3]))+0.8660254037844386*(fr[2]+fl[2])+0.5*fr[1]-0.5*(fl[1]+fr[0])+0.5*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[1] = 0.5*((-0.8660254037844386*(fr[3]+fr[2]))+0.8660254037844386*(fl[3]+fl[2])+0.5*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.8660254037844386*(fr[3]+fl[3]+fr[2]+fl[2])-0.5*(fr[1]+fr[0])+0.5*(fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 

  double fUpwind[2];
  fUpwind[0] = 0.7071067811865475*(fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.7071067811865475*(fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  double Gdiff[4]; 
  double Ghat[4]; 
  double incr2[4]; 

  incr2[2] = (nuVtSqSum[1]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])))*rdvSq4R; 
  incr2[3] = (nuVtSqSum[0]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])))*rdvSq4R; 

  Gdiff[0] = nuVtSqSum[1]*((-1.530931089239486*(fr[3]+fl[3]))+1.590990257669731*fr[1]-1.590990257669731*fl[1])+nuVtSqSum[0]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0]); 
  Gdiff[1] = nuVtSqSum[0]*((-1.530931089239486*(fr[3]+fl[3]))+1.590990257669731*fr[1]-1.590990257669731*fl[1])+nuVtSqSum[1]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0]); 

  Ghat[0] = Gdiff[0]*rdv+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]*rdv2R; 
  incr1[1] = -0.5*Ghat[1]*rdv2R; 
  incr1[2] = 0.8660254037844386*Ghat[0]*rdv2R; 
  incr1[3] = 0.8660254037844386*Ghat[1]*rdv2R; 

  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[2], outrPos[2]; 
  fluxFracL = cflRateByDirL[0] == 0. ? 0.5 : cflRateByDirL[1]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] == 0. ? 0.5 : cflRateByDirR[1]/cflRateByDirR[0]; 
  outlPos[0] = 0.1666666666666667*(incr2[3]-1.0*incr1[3]-1.732050807568877*incr2[2]+1.732050807568877*(incr1[2]+incr1[1])-3.0*incr1[0]); 
  outlPos[1] = -0.1666666666666667*(incr2[3]-1.0*incr1[3]+1.732050807568877*incr2[2]-1.732050807568877*incr1[2]+1.732050807568877*incr1[1]+3.0*incr1[0]); 
  outrPos[0] = 0.1666666666666667*(incr2[3]+incr1[3]-1.732050807568877*(incr2[2]+incr1[2]+incr1[1])+3.0*incr1[0]); 
  outrPos[1] = -0.1666666666666667*(incr2[3]+incr1[3]+1.732050807568877*(incr2[2]+incr1[2])-1.732050807568877*incr1[1]-3.0*incr1[0]); 
  if(outlPos[0] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.1666666666666667*(fl[3]-1.732050807568877*fl[2]+1.732050807568877*fl[1]-3.0*fl[0]))/dtApprox/outlPos[0]); 
  else limFac = 1.0; 
  if(outrPos[0] < 0.) limFac = std::min(limFac, -fluxFracR*(0.1666666666666667*(fr[3]-1.732050807568877*(fr[2]+fr[1])+3.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[0] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[1] < 0.) limFac = std::min(1.0, -fluxFracL*(0.1666666666666667*(fl[3]+1.732050807568877*(fl[2]+fl[1])+3.0*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[1] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.1666666666666667*(fr[3]+1.732050807568877*fr[2]-1.732050807568877*fr[1]-3.0*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outrPos[1] *= limFac; 
  outr[0] += 0.5*outrPos[1]+0.5*outrPos[0]; 
  outr[1] += 0.8660254037844386*outrPos[1]-0.8660254037844386*outrPos[0]; 
  outr[2] += (-0.8660254037844386*outrPos[1])-0.8660254037844386*outrPos[0]; 
  outr[3] += 1.5*outrPos[0]-1.5*outrPos[1]; 

  outl[0] += 0.5*outlPos[1]+0.5*outlPos[0]; 
  outl[1] += 0.8660254037844386*outlPos[1]-0.8660254037844386*outlPos[0]; 
  outl[2] += 0.8660254037844386*outlPos[1]+0.8660254037844386*outlPos[0]; 
  outl[3] += 1.5*outlPos[1]-1.5*outlPos[0]; 

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
