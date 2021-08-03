#include <VmLBOModDecl.h> 
double VmLBOconstNuUpwindSurf1x2vTensor_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[4]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = 4.0/(dxvl[1]*dxvl[1]); 
  double rdvSq4R = 4.0/(dxvr[1]*dxvr[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (2.0*wl[1]+dxvl[1])*nuSum-1.414213562373095*sumNuUx[0]; 
  alphaDrSurf[1] = -1.414213562373095*sumNuUx[1]; 

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


  incr2[2] = -0.1767766952966368*(2.0*nuVtSqSum[1]*fr[4]-2.0*nuVtSqSum[1]*fl[4]+2.0*nuVtSqSum[0]*fr[2]-2.0*nuVtSqSum[0]*fl[2]+((-1.732050807568877*fr[1])-1.732050807568877*fl[1])*nuVtSqSum[1]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[0]); 
  incr2[4] = -0.1767766952966368*(2.0*nuVtSqSum[0]*fr[4]-2.0*nuVtSqSum[0]*fl[4]+2.0*nuVtSqSum[1]*fr[2]-2.0*nuVtSqSum[1]*fl[2]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0]*fr[1]-1.732050807568877*nuVtSqSum[0]*fl[1]); 
  incr2[6] = -0.1767766952966368*(2.0*nuVtSqSum[1]*fr[7]-2.0*nuVtSqSum[1]*fl[7]+2.0*nuVtSqSum[0]*fr[6]-2.0*nuVtSqSum[0]*fl[6]-1.732050807568877*nuVtSqSum[1]*fr[5]-1.732050807568877*nuVtSqSum[1]*fl[5]-1.732050807568877*nuVtSqSum[0]*fr[3]-1.732050807568877*nuVtSqSum[0]*fl[3]); 
  incr2[7] = -0.1767766952966368*(2.0*nuVtSqSum[0]*fr[7]-2.0*nuVtSqSum[0]*fl[7]+2.0*nuVtSqSum[1]*fr[6]-2.0*nuVtSqSum[1]*fl[6]-1.732050807568877*nuVtSqSum[0]*fr[5]-1.732050807568877*nuVtSqSum[0]*fl[5]-1.732050807568877*nuVtSqSum[1]*fr[3]-1.732050807568877*nuVtSqSum[1]*fl[3]); 


  Gdiff[0] = -0.0883883476483184*(1.732050807568877*(5.0*nuVtSqSum[1]*fr[4]+5.0*nuVtSqSum[1]*fl[4]+5.0*nuVtSqSum[0]*fr[2]+5.0*nuVtSqSum[0]*fl[2])+(9.0*fl[1]-9.0*fr[1])*nuVtSqSum[1]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[0]); 
  Gdiff[1] = -0.0883883476483184*(1.732050807568877*(5.0*nuVtSqSum[0]*fr[4]+5.0*nuVtSqSum[0]*fl[4]+5.0*nuVtSqSum[1]*fr[2]+5.0*nuVtSqSum[1]*fl[2])+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[1]-9.0*nuVtSqSum[0]*fr[1]+9.0*nuVtSqSum[0]*fl[1]); 
  Gdiff[3] = -0.0883883476483184*(1.732050807568877*(5.0*nuVtSqSum[1]*fr[7]+5.0*nuVtSqSum[1]*fl[7]+5.0*nuVtSqSum[0]*fr[6]+5.0*nuVtSqSum[0]*fl[6])-9.0*nuVtSqSum[1]*fr[5]+9.0*nuVtSqSum[1]*fl[5]-9.0*nuVtSqSum[0]*fr[3]+9.0*nuVtSqSum[0]*fl[3]); 
  Gdiff[5] = -0.0883883476483184*(1.732050807568877*(5.0*nuVtSqSum[0]*fr[7]+5.0*nuVtSqSum[0]*fl[7]+5.0*nuVtSqSum[1]*fr[6]+5.0*nuVtSqSum[1]*fl[6])-9.0*nuVtSqSum[0]*fr[5]+9.0*nuVtSqSum[0]*fl[5]-9.0*nuVtSqSum[1]*fr[3]+9.0*nuVtSqSum[1]*fl[3]); 

  Ghat[0] = Gdiff[0]*rdv2L+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[3] = Gdiff[3]*rdv2L+0.7071067811865475*alphaDrSurf[1]*fUpwind[3]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[5] = Gdiff[5]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fUpwind[3]+0.7071067811865475*alphaDrSurf[1]*fUpwind[2]; 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[1]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[3]; 
  incr1[7] = 0.8660254037844386*Ghat[5]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 

  return std::abs(wl[1]-(0.7071067811865475*sumNuUx[0])/nuSum); 
} 
double VmLBOconstNuUpwindSurf1x2vTensor_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[4]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUy = &nuUSum[2]; 

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (2.0*wl[2]+dxvl[2])*nuSum-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf[1] = -1.414213562373095*sumNuUy[1]; 

  double fUpwindQuad[4];
  fUpwindQuad[0] = 0.5*((-0.6123724356957944*fr[7])+0.6123724356957944*(fl[7]+fr[6]+fr[5])-0.6123724356957944*(fl[6]+fl[5])+0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2]+fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6]+fr[5]+fl[5])-0.3535533905932737*fr[4]+0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*(fr[2]+fr[1])-0.3535533905932737*(fl[2]+fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[1] = 0.5*(0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6]+fr[5])+0.6123724356957944*fl[5]-0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2])+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.6123724356957944*(fr[5]+fl[5])+0.3535533905932737*fr[4]-0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*fr[2]-0.3535533905932737*(fl[2]+fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[2] = 0.5*(0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*(fl[6]+fr[5])-0.6123724356957944*fl[5]-0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])-0.6123724356957944*(fr[5]+fl[5])+0.3535533905932737*fr[4]-0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*fr[2]+0.3535533905932737*(fl[2]+fr[1])-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[3] = 0.5*((-0.6123724356957944*(fr[7]+fr[6]+fr[5]))+0.6123724356957944*(fl[7]+fl[6]+fl[5])+0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.3535533905932737*fr[4]+0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*(fr[2]+fr[1]+fr[0])+0.3535533905932737*(fl[2]+fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 

  double fUpwind[4];
  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  double Gdiff[8]; 
  double Ghat[8]; 
  double incr2[8]; 


  incr2[3] = -0.1767766952966368*(2.0*nuVtSqSum[1]*fr[5]-2.0*nuVtSqSum[1]*fl[5]+2.0*nuVtSqSum[0]*fr[3]-2.0*nuVtSqSum[0]*fl[3]+((-1.732050807568877*fr[1])-1.732050807568877*fl[1])*nuVtSqSum[1]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[0]); 
  incr2[5] = -0.1767766952966368*(2.0*nuVtSqSum[0]*fr[5]-2.0*nuVtSqSum[0]*fl[5]+2.0*nuVtSqSum[1]*fr[3]-2.0*nuVtSqSum[1]*fl[3]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0]*fr[1]-1.732050807568877*nuVtSqSum[0]*fl[1]); 
  incr2[6] = -0.1767766952966368*(2.0*nuVtSqSum[1]*fr[7]-2.0*nuVtSqSum[1]*fl[7]+2.0*nuVtSqSum[0]*fr[6]-2.0*nuVtSqSum[0]*fl[6]-1.732050807568877*nuVtSqSum[1]*fr[4]-1.732050807568877*nuVtSqSum[1]*fl[4]-1.732050807568877*nuVtSqSum[0]*fr[2]-1.732050807568877*nuVtSqSum[0]*fl[2]); 
  incr2[7] = -0.1767766952966368*(2.0*nuVtSqSum[0]*fr[7]-2.0*nuVtSqSum[0]*fl[7]+2.0*nuVtSqSum[1]*fr[6]-2.0*nuVtSqSum[1]*fl[6]-1.732050807568877*nuVtSqSum[0]*fr[4]-1.732050807568877*nuVtSqSum[0]*fl[4]-1.732050807568877*nuVtSqSum[1]*fr[2]-1.732050807568877*nuVtSqSum[1]*fl[2]); 


  Gdiff[0] = -0.0883883476483184*(1.732050807568877*(5.0*nuVtSqSum[1]*fr[5]+5.0*nuVtSqSum[1]*fl[5]+5.0*nuVtSqSum[0]*fr[3]+5.0*nuVtSqSum[0]*fl[3])+(9.0*fl[1]-9.0*fr[1])*nuVtSqSum[1]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[0]); 
  Gdiff[1] = -0.0883883476483184*(1.732050807568877*(5.0*nuVtSqSum[0]*fr[5]+5.0*nuVtSqSum[0]*fl[5]+5.0*nuVtSqSum[1]*fr[3]+5.0*nuVtSqSum[1]*fl[3])+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[1]-9.0*nuVtSqSum[0]*fr[1]+9.0*nuVtSqSum[0]*fl[1]); 
  Gdiff[2] = -0.0883883476483184*(1.732050807568877*(5.0*nuVtSqSum[1]*fr[7]+5.0*nuVtSqSum[1]*fl[7]+5.0*nuVtSqSum[0]*fr[6]+5.0*nuVtSqSum[0]*fl[6])-9.0*nuVtSqSum[1]*fr[4]+9.0*nuVtSqSum[1]*fl[4]-9.0*nuVtSqSum[0]*fr[2]+9.0*nuVtSqSum[0]*fl[2]); 
  Gdiff[4] = -0.0883883476483184*(1.732050807568877*(5.0*nuVtSqSum[0]*fr[7]+5.0*nuVtSqSum[0]*fl[7]+5.0*nuVtSqSum[1]*fr[6]+5.0*nuVtSqSum[1]*fl[6])-9.0*nuVtSqSum[0]*fr[4]+9.0*nuVtSqSum[0]*fl[4]-9.0*nuVtSqSum[1]*fr[2]+9.0*nuVtSqSum[1]*fl[2]); 

  Ghat[0] = Gdiff[0]*rdv2L+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = Gdiff[2]*rdv2L+0.7071067811865475*alphaDrSurf[1]*fUpwind[3]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[4] = Gdiff[4]*rdv2L+0.7071067811865475*alphaDrSurf[0]*fUpwind[3]+0.7071067811865475*alphaDrSurf[1]*fUpwind[2]; 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = 0.8660254037844386*Ghat[2]; 
  incr1[7] = 0.8660254037844386*Ghat[4]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 

  return std::abs(wl[2]-(0.7071067811865475*sumNuUy[0])/nuSum); 
} 
