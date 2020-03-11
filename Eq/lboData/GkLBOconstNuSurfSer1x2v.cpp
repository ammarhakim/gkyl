#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf1x2vSer_Vpar_P1(const double m_, const double *cflRateByDirL, const double *cflRateByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // m_:              species mass. 
  // cflRateByDir[3]: CFL rate in each direction. 
  // w[3]:            Cell-center coordinates. 
  // dxv[3]:          Cell spacing. 
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

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (2.0*wl[1]+dxvl[1])*nuSum-1.414213562373095*nuUSum[0]; 
  alphaDrSurf[1] = -1.414213562373095*nuUSum[1]; 

  double fUpwindQuad[4];
  fUpwindQuad[0] = 0.5*((-0.6123724356957944*fr[7])+0.6123724356957944*(fl[7]+fr[6])-0.6123724356957944*fl[6]+0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[1] = 0.5*(0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6])-0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[2] = 0.5*(0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*fl[6]-0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[3] = 0.5*((-0.6123724356957944*(fr[7]+fr[6]))+0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 

  double fUpwind[4];
  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  double Gdiff[8]; 
  double Ghat[8]; 
  double incr2[8]; 

  incr2[2] = (nuVtSqSum[1]*((-0.3535533905932737*fr[4])+0.3535533905932737*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])))*rdvSq4R; 
  incr2[4] = (nuVtSqSum[0]*((-0.3535533905932737*fr[4])+0.3535533905932737*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])))*rdvSq4R; 
  incr2[6] = (nuVtSqSum[1]*((-0.3535533905932737*fr[7])+0.3535533905932737*fl[7]+0.3061862178478971*(fr[5]+fl[5]))+nuVtSqSum[0]*((-0.3535533905932737*fr[6])+0.3535533905932737*fl[6]+0.3061862178478971*(fr[3]+fl[3])))*rdvSq4R; 
  incr2[7] = (nuVtSqSum[0]*((-0.3535533905932737*fr[7])+0.3535533905932737*fl[7]+0.3061862178478971*(fr[5]+fl[5]))+nuVtSqSum[1]*((-0.3535533905932737*fr[6])+0.3535533905932737*fl[6]+0.3061862178478971*(fr[3]+fl[3])))*rdvSq4R; 

  Gdiff[0] = nuVtSqSum[1]*((-1.530931089239486*(fr[4]+fl[4]))+1.590990257669731*fr[1]-1.590990257669731*fl[1])+nuVtSqSum[0]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0]); 
  Gdiff[1] = nuVtSqSum[0]*((-1.530931089239486*(fr[4]+fl[4]))+1.590990257669731*fr[1]-1.590990257669731*fl[1])+nuVtSqSum[1]*((-1.530931089239486*(fr[2]+fl[2]))+1.590990257669731*fr[0]-1.590990257669731*fl[0]); 
  Gdiff[3] = nuVtSqSum[1]*((-1.530931089239486*(fr[7]+fl[7]))+1.590990257669731*fr[5]-1.590990257669731*fl[5])+nuVtSqSum[0]*((-1.530931089239486*(fr[6]+fl[6]))+1.590990257669731*fr[3]-1.590990257669731*fl[3]); 
  Gdiff[5] = nuVtSqSum[0]*((-1.530931089239486*(fr[7]+fl[7]))+1.590990257669731*fr[5]-1.590990257669731*fl[5])+nuVtSqSum[1]*((-1.530931089239486*(fr[6]+fl[6]))+1.590990257669731*fr[3]-1.590990257669731*fl[3]); 

  Ghat[0] = Gdiff[0]*rdv+0.7071067811865475*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = Gdiff[1]*rdv+0.7071067811865475*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[3] = Gdiff[3]*rdv+0.7071067811865475*(alphaDrSurf[1]*fUpwind[3]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[5] = Gdiff[5]*rdv+0.7071067811865475*(alphaDrSurf[0]*fUpwind[3]+alphaDrSurf[1]*fUpwind[2]); 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]*rdv2R; 
  incr1[1] = -0.5*Ghat[1]*rdv2R; 
  incr1[2] = 0.8660254037844386*Ghat[0]*rdv2R; 
  incr1[3] = -0.5*Ghat[3]*rdv2R; 
  incr1[4] = 0.8660254037844386*Ghat[1]*rdv2R; 
  incr1[5] = -0.5*Ghat[5]*rdv2R; 
  incr1[6] = 0.8660254037844386*Ghat[3]*rdv2R; 
  incr1[7] = 0.8660254037844386*Ghat[5]*rdv2R; 

  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[4], outrPos[4]; 
  fluxFracL = cflRateByDirL[0] == 0. ? 0.3333333333333333 : cflRateByDirL[1]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] == 0. ? 0.3333333333333333 : cflRateByDirR[1]/cflRateByDirR[0]; 
  outlPos[0] = -0.02777777777777778*(2.449489742783178*incr2[7]-2.449489742783178*incr1[7]-4.242640687119286*incr2[6]+4.242640687119286*(incr1[6]+incr1[5])-4.242640687119286*incr2[4]+4.242640687119286*incr1[4]-7.348469228349534*incr1[3]+7.348469228349534*incr2[2]-7.348469228349534*(incr1[2]+incr1[1])+12.72792206135786*incr1[0]); 
  outlPos[1] = 0.02777777777777778*(2.449489742783178*incr2[7]-2.449489742783178*incr1[7]+4.242640687119286*incr2[6]-4.242640687119286*incr1[6]+4.242640687119286*incr1[5]-4.242640687119286*incr2[4]+4.242640687119286*incr1[4]+7.348469228349534*incr1[3]-7.348469228349534*incr2[2]+7.348469228349534*incr1[2]-7.348469228349534*incr1[1]-12.72792206135786*incr1[0]); 
  outlPos[2] = 0.02777777777777778*(2.449489742783178*incr2[7]-2.449489742783178*incr1[7]-4.242640687119286*incr2[6]+4.242640687119286*(incr1[6]+incr1[5]+incr2[4])-4.242640687119286*incr1[4]-7.348469228349534*(incr1[3]+incr2[2])+7.348469228349534*(incr1[2]+incr1[1])-12.72792206135786*incr1[0]); 
  outlPos[3] = -0.02777777777777778*(2.449489742783178*incr2[7]-2.449489742783178*incr1[7]+4.242640687119286*incr2[6]-4.242640687119286*incr1[6]+4.242640687119286*(incr1[5]+incr2[4])-4.242640687119286*incr1[4]+7.348469228349534*(incr1[3]+incr2[2])-7.348469228349534*incr1[2]+7.348469228349534*incr1[1]+12.72792206135786*incr1[0]); 
  outrPos[0] = -0.02777777777777778*(2.449489742783178*(incr2[7]+incr1[7])-4.242640687119286*(incr2[6]+incr1[6]+incr1[5]+incr2[4]+incr1[4])+7.348469228349534*(incr1[3]+incr2[2]+incr1[2]+incr1[1])-12.72792206135786*incr1[0]); 
  outrPos[1] = 0.02777777777777778*(2.449489742783178*(incr2[7]+incr1[7])+4.242640687119286*(incr2[6]+incr1[6])-4.242640687119286*(incr1[5]+incr2[4]+incr1[4])-7.348469228349534*(incr1[3]+incr2[2]+incr1[2])+7.348469228349534*incr1[1]+12.72792206135786*incr1[0]); 
  outrPos[2] = 0.02777777777777778*(2.449489742783178*(incr2[7]+incr1[7])-4.242640687119286*(incr2[6]+incr1[6]+incr1[5])+4.242640687119286*(incr2[4]+incr1[4])+7.348469228349534*incr1[3]-7.348469228349534*(incr2[2]+incr1[2]+incr1[1])+12.72792206135786*incr1[0]); 
  outrPos[3] = -0.02777777777777778*(2.449489742783178*(incr2[7]+incr1[7])+4.242640687119286*(incr2[6]+incr1[6])-4.242640687119286*incr1[5]+4.242640687119286*(incr2[4]+incr1[4])-7.348469228349534*incr1[3]+7.348469228349534*(incr2[2]+incr1[2])-7.348469228349534*incr1[1]-12.72792206135786*incr1[0]); 
  if(outlPos[0] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*fl[5]-4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]+12.72792206135786*fl[0]))/dtApprox/outlPos[0]); 
  else limFac = 1.0; 
  if(outrPos[0] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[0] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[1] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5])-4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])-12.72792206135786*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[1] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outrPos[1] *= limFac; 
  if(outlPos[2] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 1.0; 
  if(outrPos[2] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5])+4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])+12.72792206135786*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outrPos[2] *= limFac; 
  if(outlPos[3] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[3] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*fr[5]+4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]-12.72792206135786*fr[0]))/dtApprox/outrPos[3]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outrPos[3] *= limFac; 
  outr[0] += 0.3535533905932737*outrPos[3]+0.3535533905932737*outrPos[2]+0.3535533905932737*outrPos[1]+0.3535533905932737*outrPos[0]; 
  outr[1] += 0.6123724356957944*outrPos[3]-0.6123724356957944*outrPos[2]+0.6123724356957944*outrPos[1]-0.6123724356957944*outrPos[0]; 
  outr[2] += (-0.6123724356957944*outrPos[3])-0.6123724356957944*outrPos[2]-0.6123724356957944*outrPos[1]-0.6123724356957944*outrPos[0]; 
  outr[3] += 0.6123724356957944*outrPos[3]+0.6123724356957944*outrPos[2]-0.6123724356957944*outrPos[1]-0.6123724356957944*outrPos[0]; 
  outr[4] += (-1.060660171779821*outrPos[3])+1.060660171779821*outrPos[2]-1.060660171779821*outrPos[1]+1.060660171779821*outrPos[0]; 
  outr[5] += 1.060660171779821*outrPos[3]-1.060660171779821*outrPos[2]-1.060660171779821*outrPos[1]+1.060660171779821*outrPos[0]; 
  outr[6] += (-1.060660171779821*outrPos[3])-1.060660171779821*outrPos[2]+1.060660171779821*outrPos[1]+1.060660171779821*outrPos[0]; 
  outr[7] += (-1.837117307087383*outrPos[3])+1.837117307087383*outrPos[2]+1.837117307087383*outrPos[1]-1.837117307087383*outrPos[0]; 

  outl[0] += 0.3535533905932737*outlPos[3]+0.3535533905932737*outlPos[2]+0.3535533905932737*outlPos[1]+0.3535533905932737*outlPos[0]; 
  outl[1] += 0.6123724356957944*outlPos[3]-0.6123724356957944*outlPos[2]+0.6123724356957944*outlPos[1]-0.6123724356957944*outlPos[0]; 
  outl[2] += 0.6123724356957944*outlPos[3]+0.6123724356957944*outlPos[2]+0.6123724356957944*outlPos[1]+0.6123724356957944*outlPos[0]; 
  outl[3] += 0.6123724356957944*outlPos[3]+0.6123724356957944*outlPos[2]-0.6123724356957944*outlPos[1]-0.6123724356957944*outlPos[0]; 
  outl[4] += 1.060660171779821*outlPos[3]-1.060660171779821*outlPos[2]+1.060660171779821*outlPos[1]-1.060660171779821*outlPos[0]; 
  outl[5] += 1.060660171779821*outlPos[3]-1.060660171779821*outlPos[2]-1.060660171779821*outlPos[1]+1.060660171779821*outlPos[0]; 
  outl[6] += 1.060660171779821*outlPos[3]+1.060660171779821*outlPos[2]-1.060660171779821*outlPos[1]-1.060660171779821*outlPos[0]; 
  outl[7] += 1.837117307087383*outlPos[3]-1.837117307087383*outlPos[2]-1.837117307087383*outlPos[1]+1.837117307087383*outlPos[0]; 

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
double GkLBOconstNuSurf1x2vSer_Mu_P1(const double m_, const double *cflRateByDirL, const double *cflRateByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // m_:              species mass. 
  // cflRateByDir[3]: CFL rate in each direction. 
  // w[3]:            Cell-center coordinates. 
  // dxv[3]:          Cell spacing. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:       maximum midpoint value of v-u. 
  // nuUSum[2]:       sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:    sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:           Distribution function in left/right cells 
  // outl/outr:       Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (4.0*wl[2]+2.0*dxvl[2])*nuSum; 

  double fUpwindQuad[4];
  fUpwindQuad[0] = 0.5*((-0.6123724356957944*fr[7])+0.6123724356957944*(fl[7]+fr[6]+fr[5])-0.6123724356957944*(fl[6]+fl[5])+0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2]+fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6]+fr[5]+fl[5])-0.3535533905932737*fr[4]+0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*(fr[2]+fr[1])-0.3535533905932737*(fl[2]+fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]); 
  fUpwindQuad[1] = 0.5*(0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6]+fr[5])+0.6123724356957944*fl[5]-0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2])+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.6123724356957944*(fr[5]+fl[5])+0.3535533905932737*fr[4]-0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*fr[2]-0.3535533905932737*(fl[2]+fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaDrSurf[0]); 
  fUpwindQuad[2] = 0.5*(0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*(fl[6]+fr[5])-0.6123724356957944*fl[5]-0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])-0.6123724356957944*(fr[5]+fl[5])+0.3535533905932737*fr[4]-0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*fr[2]+0.3535533905932737*(fl[2]+fr[1])-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]); 
  fUpwindQuad[3] = 0.5*((-0.6123724356957944*(fr[7]+fr[6]+fr[5]))+0.6123724356957944*(fl[7]+fl[6]+fl[5])+0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.3535533905932737*fr[4]+0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*(fr[2]+fr[1]+fr[0])+0.3535533905932737*(fl[2]+fl[1]+fl[0]))*sgn(alphaDrSurf[0]); 

  double fUpwind[4];
  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  double diffFac[2]; 
  diffFac[0] = (BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 
  diffFac[1] = (BmagInv[0]*nuVtSqSum[1]+nuVtSqSum[0]*BmagInv[1])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 

  double Gdiff[8]; 
  double Ghat[8]; 
  double incr2[8]; 

  incr2[3] = (diffFac[1]*(0.3535533905932737*fl[5]-0.3535533905932737*fr[5])+diffFac[0]*(0.3535533905932737*fl[3]-0.3535533905932737*fr[3])+0.3061862178478971*(diffFac[1]*(fr[1]+fl[1])+diffFac[0]*(fr[0]+fl[0])))*rdvSq4R; 
  incr2[5] = (diffFac[0]*(0.3535533905932737*fl[5]-0.3535533905932737*fr[5])+diffFac[1]*(0.3535533905932737*fl[3]-0.3535533905932737*fr[3])+0.3061862178478971*(diffFac[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*diffFac[1]))*rdvSq4R; 
  incr2[6] = (diffFac[1]*(0.3535533905932737*fl[7]-0.3535533905932737*fr[7])+diffFac[0]*(0.3535533905932737*fl[6]-0.3535533905932737*fr[6])+0.3061862178478971*(diffFac[1]*(fr[4]+fl[4])+diffFac[0]*(fr[2]+fl[2])))*rdvSq4R; 
  incr2[7] = (diffFac[0]*(0.3535533905932737*fl[7]-0.3535533905932737*fr[7])+diffFac[1]*(0.3535533905932737*fl[6]-0.3535533905932737*fr[6])+0.3061862178478971*(diffFac[0]*(fr[4]+fl[4])+diffFac[1]*(fr[2]+fl[2])))*rdvSq4R; 

  Gdiff[0] = (-1.530931089239486*(diffFac[1]*(fr[5]+fl[5])+diffFac[0]*(fr[3]+fl[3])))+diffFac[1]*(1.590990257669731*fr[1]-1.590990257669731*fl[1])+diffFac[0]*(1.590990257669731*fr[0]-1.590990257669731*fl[0]); 
  Gdiff[1] = (-1.530931089239486*(diffFac[0]*(fr[5]+fl[5])+diffFac[1]*(fr[3]+fl[3])))+diffFac[0]*(1.590990257669731*fr[1]-1.590990257669731*fl[1])+(1.590990257669731*fr[0]-1.590990257669731*fl[0])*diffFac[1]; 
  Gdiff[2] = (-1.530931089239486*(diffFac[1]*(fr[7]+fl[7])+diffFac[0]*(fr[6]+fl[6])))+diffFac[1]*(1.590990257669731*fr[4]-1.590990257669731*fl[4])+diffFac[0]*(1.590990257669731*fr[2]-1.590990257669731*fl[2]); 
  Gdiff[4] = (-1.530931089239486*(diffFac[0]*(fr[7]+fl[7])+diffFac[1]*(fr[6]+fl[6])))+diffFac[0]*(1.590990257669731*fr[4]-1.590990257669731*fl[4])+diffFac[1]*(1.590990257669731*fr[2]-1.590990257669731*fl[2]); 

  Ghat[0] = Gdiff[0]*rdv+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]; 
  Ghat[2] = Gdiff[2]*rdv+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[4] = Gdiff[4]*rdv+0.7071067811865475*alphaDrSurf[0]*fUpwind[3]; 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]*rdv2R; 
  incr1[1] = -0.5*Ghat[1]*rdv2R; 
  incr1[2] = -0.5*Ghat[2]*rdv2R; 
  incr1[3] = 0.8660254037844386*Ghat[0]*rdv2R; 
  incr1[4] = -0.5*Ghat[4]*rdv2R; 
  incr1[5] = 0.8660254037844386*Ghat[1]*rdv2R; 
  incr1[6] = 0.8660254037844386*Ghat[2]*rdv2R; 
  incr1[7] = 0.8660254037844386*Ghat[4]*rdv2R; 

  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[4], outrPos[4]; 
  fluxFracL = cflRateByDirL[0] == 0. ? 0.3333333333333333 : cflRateByDirL[2]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] == 0. ? 0.3333333333333333 : cflRateByDirR[2]/cflRateByDirR[0]; 
  outlPos[0] = -0.02777777777777778*(2.449489742783178*incr2[7]-2.449489742783178*incr1[7]-4.242640687119286*incr2[6]+4.242640687119286*incr1[6]-4.242640687119286*incr2[5]+4.242640687119286*(incr1[5]+incr1[4])+7.348469228349534*incr2[3]-7.348469228349534*(incr1[3]+incr1[2]+incr1[1])+12.72792206135786*incr1[0]); 
  outlPos[1] = 0.02777777777777778*(2.449489742783178*incr2[7]-2.449489742783178*incr1[7]+4.242640687119286*incr2[6]-4.242640687119286*(incr1[6]+incr2[5])+4.242640687119286*(incr1[5]+incr1[4])-7.348469228349534*incr2[3]+7.348469228349534*(incr1[3]+incr1[2])-7.348469228349534*incr1[1]-12.72792206135786*incr1[0]); 
  outlPos[2] = 0.02777777777777778*(2.449489742783178*incr2[7]-2.449489742783178*incr1[7]-4.242640687119286*incr2[6]+4.242640687119286*(incr1[6]+incr2[5])-4.242640687119286*incr1[5]+4.242640687119286*incr1[4]-7.348469228349534*incr2[3]+7.348469228349534*incr1[3]-7.348469228349534*incr1[2]+7.348469228349534*incr1[1]-12.72792206135786*incr1[0]); 
  outlPos[3] = -0.02777777777777778*(2.449489742783178*incr2[7]-2.449489742783178*incr1[7]+4.242640687119286*incr2[6]-4.242640687119286*incr1[6]+4.242640687119286*incr2[5]-4.242640687119286*incr1[5]+4.242640687119286*incr1[4]+7.348469228349534*incr2[3]-7.348469228349534*incr1[3]+7.348469228349534*(incr1[2]+incr1[1])+12.72792206135786*incr1[0]); 
  outrPos[0] = -0.02777777777777778*(2.449489742783178*(incr2[7]+incr1[7])-4.242640687119286*(incr2[6]+incr1[6]+incr2[5]+incr1[5]+incr1[4])+7.348469228349534*(incr2[3]+incr1[3]+incr1[2]+incr1[1])-12.72792206135786*incr1[0]); 
  outrPos[1] = 0.02777777777777778*(2.449489742783178*(incr2[7]+incr1[7])+4.242640687119286*(incr2[6]+incr1[6])-4.242640687119286*(incr2[5]+incr1[5]+incr1[4])-7.348469228349534*(incr2[3]+incr1[3]+incr1[2])+7.348469228349534*incr1[1]+12.72792206135786*incr1[0]); 
  outrPos[2] = 0.02777777777777778*(2.449489742783178*(incr2[7]+incr1[7])-4.242640687119286*(incr2[6]+incr1[6])+4.242640687119286*(incr2[5]+incr1[5])-4.242640687119286*incr1[4]-7.348469228349534*(incr2[3]+incr1[3])+7.348469228349534*incr1[2]-7.348469228349534*incr1[1]+12.72792206135786*incr1[0]); 
  outrPos[3] = -0.02777777777777778*(2.449489742783178*(incr2[7]+incr1[7])+4.242640687119286*(incr2[6]+incr1[6]+incr2[5]+incr1[5])-4.242640687119286*incr1[4]+7.348469228349534*(incr2[3]+incr1[3])-7.348469228349534*(incr1[2]+incr1[1])-12.72792206135786*incr1[0]); 
  if(outlPos[0] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*(fl[6]+fl[5])+4.242640687119286*fl[4]+7.348469228349534*fl[3]-7.348469228349534*(fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPos[0]); 
  else limFac = 1.0; 
  if(outrPos[0] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[5]+fr[4])+7.348469228349534*(fr[3]+fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[0] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[1] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*fl[5]+4.242640687119286*fl[4]-7.348469228349534*fl[3]+7.348469228349534*fl[2]-7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[1] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*(fr[5]+fr[4])-7.348469228349534*(fr[3]+fr[2])+7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outrPos[1] *= limFac; 
  if(outlPos[2] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*(fl[5]+fl[4])-7.348469228349534*(fl[3]+fl[2])+7.348469228349534*fl[1]-12.72792206135786*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 1.0; 
  if(outrPos[2] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*fr[5]-4.242640687119286*fr[4]-7.348469228349534*fr[3]+7.348469228349534*fr[2]-7.348469228349534*fr[1]+12.72792206135786*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outrPos[2] *= limFac; 
  if(outlPos[3] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[5]+fl[4])+7.348469228349534*(fl[3]+fl[2]+fl[1])+12.72792206135786*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[3] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(2.449489742783178*fr[7]+4.242640687119286*(fr[6]+fr[5])-4.242640687119286*fr[4]+7.348469228349534*fr[3]-7.348469228349534*(fr[2]+fr[1])-12.72792206135786*fr[0]))/dtApprox/outrPos[3]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outrPos[3] *= limFac; 
  outr[0] += 0.3535533905932737*outrPos[3]+0.3535533905932737*outrPos[2]+0.3535533905932737*outrPos[1]+0.3535533905932737*outrPos[0]; 
  outr[1] += 0.6123724356957944*outrPos[3]-0.6123724356957944*outrPos[2]+0.6123724356957944*outrPos[1]-0.6123724356957944*outrPos[0]; 
  outr[2] += 0.6123724356957944*outrPos[3]+0.6123724356957944*outrPos[2]-0.6123724356957944*outrPos[1]-0.6123724356957944*outrPos[0]; 
  outr[3] += (-0.6123724356957944*outrPos[3])-0.6123724356957944*outrPos[2]-0.6123724356957944*outrPos[1]-0.6123724356957944*outrPos[0]; 
  outr[4] += 1.060660171779821*outrPos[3]-1.060660171779821*outrPos[2]-1.060660171779821*outrPos[1]+1.060660171779821*outrPos[0]; 
  outr[5] += (-1.060660171779821*outrPos[3])+1.060660171779821*outrPos[2]-1.060660171779821*outrPos[1]+1.060660171779821*outrPos[0]; 
  outr[6] += (-1.060660171779821*outrPos[3])-1.060660171779821*outrPos[2]+1.060660171779821*outrPos[1]+1.060660171779821*outrPos[0]; 
  outr[7] += (-1.837117307087383*outrPos[3])+1.837117307087383*outrPos[2]+1.837117307087383*outrPos[1]-1.837117307087383*outrPos[0]; 

  outl[0] += 0.3535533905932737*outlPos[3]+0.3535533905932737*outlPos[2]+0.3535533905932737*outlPos[1]+0.3535533905932737*outlPos[0]; 
  outl[1] += 0.6123724356957944*outlPos[3]-0.6123724356957944*outlPos[2]+0.6123724356957944*outlPos[1]-0.6123724356957944*outlPos[0]; 
  outl[2] += 0.6123724356957944*outlPos[3]+0.6123724356957944*outlPos[2]-0.6123724356957944*outlPos[1]-0.6123724356957944*outlPos[0]; 
  outl[3] += 0.6123724356957944*outlPos[3]+0.6123724356957944*outlPos[2]+0.6123724356957944*outlPos[1]+0.6123724356957944*outlPos[0]; 
  outl[4] += 1.060660171779821*outlPos[3]-1.060660171779821*outlPos[2]-1.060660171779821*outlPos[1]+1.060660171779821*outlPos[0]; 
  outl[5] += 1.060660171779821*outlPos[3]-1.060660171779821*outlPos[2]+1.060660171779821*outlPos[1]-1.060660171779821*outlPos[0]; 
  outl[6] += 1.060660171779821*outlPos[3]+1.060660171779821*outlPos[2]-1.060660171779821*outlPos[1]-1.060660171779821*outlPos[0]; 
  outl[7] += 1.837117307087383*outlPos[3]-1.837117307087383*outlPos[2]-1.837117307087383*outlPos[1]+1.837117307087383*outlPos[0]; 

  return std::abs(2.0*wl[2]); 
} 
