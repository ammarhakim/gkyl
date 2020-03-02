#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf1x2vSer_Vpar_P1(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[1]; 
  double rdv2L = 2.0/dxvl[1]; 
  double rdv2R = 2.0/dxvr[1]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (2.0*wl[1]+dxvl[1])*nuSum-1.414213562373095*nuUSum[0]; 
  alphaDrSurf[1] = -1.414213562373095*nuUSum[1]; 

  double fUpwindQuad[4];
  fUpwindQuad[0] = 0.5*(0.7905694150420947*(fr[8]+fl[8])-0.6123724356957944*fr[7]+0.6123724356957944*(fl[7]+fr[6])-0.6123724356957944*fl[6]+0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*((-0.7905694150420947*fr[8])+0.7905694150420947*fl[8]+0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[1] = 0.5*(0.7905694150420947*(fr[8]+fl[8])+0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6])-0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]-0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.7905694150420947*fr[8])+0.7905694150420947*fl[8]-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])+0.3535533905932737*fr[3]-0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[2] = 0.5*(0.7905694150420947*(fr[8]+fl[8])+0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*fl[6]-0.3535533905932737*(fr[5]+fl[5])+0.6123724356957944*fr[4]-0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*((-0.7905694150420947*fr[8])+0.7905694150420947*fl[8]-0.6123724356957944*(fr[7]+fl[7])+0.6123724356957944*(fr[6]+fl[6])+0.3535533905932737*fr[5]-0.3535533905932737*fl[5]-0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[3] = 0.5*(0.7905694150420947*(fr[8]+fl[8])-0.6123724356957944*(fr[7]+fr[6])+0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*(fr[5]+fl[5])-0.6123724356957944*fr[4]+0.6123724356957944*fl[4]+0.3535533905932737*(fr[3]+fl[3])-0.6123724356957944*fr[2]+0.6123724356957944*fl[2]+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.7905694150420947*fr[8])+0.7905694150420947*fl[8]+0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6])-0.3535533905932737*fr[5]+0.3535533905932737*fl[5]+0.6123724356957944*(fr[4]+fl[4])-0.3535533905932737*fr[3]+0.3535533905932737*fl[3]+0.6123724356957944*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 

  double fUpwind[4];
  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  double Gdiff[9]; 
  double Ghat[9]; 
  double incr2[9]; 

  incr2[2] = nuVtSqSum[0]*(0.2995357736356374*(fr[8]+fl[8])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[1]*((-0.430893194785552*fr[4])+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[4] = nuVtSqSum[1]*(0.2995357736356374*(fr[8]+fl[8])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[0]*((-0.430893194785552*fr[4])+0.430893194785552*fl[4]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[6] = nuVtSqSum[1]*((-0.430893194785552*fr[7])+0.430893194785552*fl[7]+0.3061862178478971*(fr[5]+fl[5]))+nuVtSqSum[0]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 
  incr2[7] = nuVtSqSum[0]*((-0.430893194785552*fr[7])+0.430893194785552*fl[7]+0.3061862178478971*(fr[5]+fl[5]))+nuVtSqSum[1]*((-0.430893194785552*fr[6])+0.430893194785552*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 
  incr2[8] = nuVtSqSum[0]*((-1.160097062884178*(fr[8]+fl[8]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0]))+nuVtSqSum[1]*(1.668842167398551*fr[4]-1.668842167398551*fl[4]-1.185854122563142*(fr[1]+fl[1])); 

  Gdiff[0] = nuVtSqSum[0]*(1.897366596101028*fr[8]-1.897366596101028*fl[8]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+nuVtSqSum[1]*((-3.368048396326869*(fr[4]+fl[4]))+2.651650429449552*fr[1]-2.651650429449552*fl[1]); 
  Gdiff[1] = nuVtSqSum[1]*(1.897366596101028*fr[8]-1.897366596101028*fl[8]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+nuVtSqSum[0]*((-3.368048396326869*(fr[4]+fl[4]))+2.651650429449552*fr[1]-2.651650429449552*fl[1]); 
  Gdiff[3] = nuVtSqSum[1]*((-3.368048396326869*(fr[7]+fl[7]))+2.651650429449552*fr[5]-2.651650429449552*fl[5])+nuVtSqSum[0]*((-3.368048396326869*(fr[6]+fl[6]))+2.651650429449552*fr[3]-2.651650429449552*fl[3]); 
  Gdiff[5] = nuVtSqSum[0]*((-3.368048396326869*(fr[7]+fl[7]))+2.651650429449552*fr[5]-2.651650429449552*fl[5])+nuVtSqSum[1]*((-3.368048396326869*(fr[6]+fl[6]))+2.651650429449552*fr[3]-2.651650429449552*fl[3]); 

  Ghat[0] = Gdiff[0]*rdv+0.7071067811865475*(alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = Gdiff[1]*rdv+0.7071067811865475*(alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[3] = Gdiff[3]*rdv+0.7071067811865475*(alphaDrSurf[1]*fUpwind[3]+alphaDrSurf[0]*fUpwind[2]); 
  Ghat[5] = Gdiff[5]*rdv+0.7071067811865475*(alphaDrSurf[0]*fUpwind[3]+alphaDrSurf[1]*fUpwind[2]); 

  double incr1[9]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[1]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[3]; 
  incr1[7] = 0.8660254037844386*Ghat[5]; 
  incr1[8] = -1.118033988749895*Ghat[0]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr2[8]*rdvSq4R+incr1[8]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 
  outl[8] += incr2[8]*rdvSq4L-1.0*incr1[8]*rdv2L; 

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
double GkLBOconstNuSurf1x2vSer_Mu_P1(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[3]:         Cell-center coordinates. 
  // dxv[3]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[2]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[5]; 
  alphaDrSurf[0] = (4.0*wl[2]+2.0*dxvl[2])*nuSum; 

  double fUpwindQuad[5];
  fUpwindQuad[0] = 0.5*((-0.6123724356957944*fr[7])+0.6123724356957944*(fl[7]+fr[6]+fr[5])-0.6123724356957944*(fl[6]+fl[5])+0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2]+fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7])-0.6123724356957944*(fr[6]+fl[6]+fr[5]+fl[5])-0.3535533905932737*fr[4]+0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*(fr[2]+fr[1])-0.3535533905932737*(fl[2]+fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]); 
  fUpwindQuad[1] = 0.5*(0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fl[7]+fl[6]+fr[5])+0.6123724356957944*fl[5]-0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]-0.3535533905932737*(fr[2]+fl[2])+0.3535533905932737*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]))+0.6123724356957944*(fr[5]+fl[5])+0.3535533905932737*fr[4]-0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])+0.3535533905932737*fr[2]-0.3535533905932737*(fl[2]+fr[1]+fr[0])+0.3535533905932737*(fl[1]+fl[0]))*sgn(alphaDrSurf[0]); 
  fUpwindQuad[2] = 0.5*(0.6123724356957944*fr[7]-0.6123724356957944*(fl[7]+fr[6])+0.6123724356957944*(fl[6]+fr[5])-0.6123724356957944*fl[5]-0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[2]+fl[2])-0.3535533905932737*(fr[1]+fl[1])+0.3535533905932737*(fr[0]+fl[0]))-0.5*((-0.6123724356957944*(fr[7]+fl[7]))+0.6123724356957944*(fr[6]+fl[6])-0.6123724356957944*(fr[5]+fl[5])+0.3535533905932737*fr[4]-0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*fr[2]+0.3535533905932737*(fl[2]+fr[1])-0.3535533905932737*(fl[1]+fr[0])+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]); 
  fUpwindQuad[3] = 0.5*((-0.6123724356957944*(fr[7]+fr[6]+fr[5]))+0.6123724356957944*(fl[7]+fl[6]+fl[5])+0.3535533905932737*(fr[4]+fl[4])-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.6123724356957944*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.3535533905932737*fr[4]+0.3535533905932737*fl[4]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*(fr[2]+fr[1]+fr[0])+0.3535533905932737*(fl[2]+fl[1]+fl[0]))*sgn(alphaDrSurf[0]); 
  fUpwindQuad[4] = 0.5*((-0.3952847075210473*(fr[8]+fl[8]))-0.6123724356957944*fr[3]+0.6123724356957944*fl[3]+0.3535533905932737*(fr[0]+fl[0]))-0.5*(0.3952847075210473*fr[8]-0.3952847075210473*fl[8]+0.6123724356957944*(fr[3]+fl[3])-0.3535533905932737*fr[0]+0.3535533905932737*fl[0])*sgn(alphaDrSurf[0]); 

  double fUpwind[5];
  fUpwind[0] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.5*(fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.5*(fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.5*(fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[4] = -0.4472135954999579*(4.0*fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 

  double diffFac[2]; 
  diffFac[0] = (BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 
  diffFac[1] = (BmagInv[0]*nuVtSqSum[1]+nuVtSqSum[0]*BmagInv[1])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 

  double Gdiff[9]; 
  double Ghat[9]; 
  double incr2[9]; 

  incr2[3] = diffFac[1]*(0.3535533905932737*fl[5]-0.3535533905932737*fr[5])+diffFac[0]*(0.3535533905932737*fl[3]-0.3535533905932737*fr[3])+0.3061862178478971*(diffFac[1]*(fr[1]+fl[1])+diffFac[0]*(fr[0]+fl[0])); 
  incr2[5] = diffFac[0]*(0.3535533905932737*fl[5]-0.3535533905932737*fr[5])+diffFac[1]*(0.3535533905932737*fl[3]-0.3535533905932737*fr[3])+0.3061862178478971*(diffFac[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*diffFac[1]); 
  incr2[6] = diffFac[1]*(0.3535533905932737*fl[7]-0.3535533905932737*fr[7])+diffFac[0]*(0.3535533905932737*fl[6]-0.3535533905932737*fr[6])+0.3061862178478971*(diffFac[1]*(fr[4]+fl[4])+diffFac[0]*(fr[2]+fl[2])); 
  incr2[7] = diffFac[0]*(0.3535533905932737*fl[7]-0.3535533905932737*fr[7])+diffFac[1]*(0.3535533905932737*fl[6]-0.3535533905932737*fr[6])+0.3061862178478971*(diffFac[0]*(fr[4]+fl[4])+diffFac[1]*(fr[2]+fl[2])); 

  Gdiff[0] = (-1.530931089239486*(diffFac[1]*(fr[5]+fl[5])+diffFac[0]*(fr[3]+fl[3])))+diffFac[1]*(1.590990257669731*fr[1]-1.590990257669731*fl[1])+diffFac[0]*(1.590990257669731*fr[0]-1.590990257669731*fl[0]); 
  Gdiff[1] = (-1.530931089239486*(diffFac[0]*(fr[5]+fl[5])+diffFac[1]*(fr[3]+fl[3])))+diffFac[0]*(1.590990257669731*fr[1]-1.590990257669731*fl[1])+(1.590990257669731*fr[0]-1.590990257669731*fl[0])*diffFac[1]; 
  Gdiff[2] = (-1.530931089239486*(diffFac[1]*(fr[7]+fl[7])+diffFac[0]*(fr[6]+fl[6])))+diffFac[1]*(1.590990257669731*fr[4]-1.590990257669731*fl[4])+diffFac[0]*(1.590990257669731*fr[2]-1.590990257669731*fl[2]); 
  Gdiff[4] = (-1.530931089239486*(diffFac[0]*(fr[7]+fl[7])+diffFac[1]*(fr[6]+fl[6])))+diffFac[0]*(1.590990257669731*fr[4]-1.590990257669731*fl[4])+diffFac[1]*(1.590990257669731*fr[2]-1.590990257669731*fl[2]); 
  Gdiff[8] = diffFac[0]*(1.590990257669731*fr[8]-1.590990257669731*fl[8]); 

  Ghat[0] = Gdiff[0]*rdv+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]; 
  Ghat[2] = Gdiff[2]*rdv+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[4] = Gdiff[4]*rdv+0.7071067811865475*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[8] = Gdiff[8]*rdv+0.7071067811865475*alphaDrSurf[0]*fUpwind[4]; 

  double incr1[9]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = 0.8660254037844386*Ghat[2]; 
  incr1[7] = 0.8660254037844386*Ghat[4]; 
  incr1[8] = -0.5*Ghat[8]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr2[5]*rdvSq4R+incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += incr1[5]*rdv2L-1.0*incr2[5]*rdvSq4L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 

  return std::abs(2.0*wl[2]); 
} 
