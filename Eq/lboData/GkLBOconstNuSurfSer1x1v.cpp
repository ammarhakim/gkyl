#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf1x1vSer_Vpar_P1(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
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

  double alphaDrSurf[2]; 
  alphaDrSurf[0] = 0.7071067811865475*((2.0*wl[1]+dxvl[1])*nuSum-1.414213562373095*nuUSum[0]); 
  alphaDrSurf[1] = -1.0*nuUSum[1]; 

  double fUpwindQuad[2];
  fUpwindQuad[0] = 0.5*(1.118033988749895*(fr[4]+fl[4])+0.8660254037844386*fr[3]-0.8660254037844386*(fl[3]+fr[2])+0.8660254037844386*fl[2]-0.5*(fr[1]+fl[1])+0.5*(fr[0]+fl[0]))-0.5*((-1.118033988749895*fr[4])+1.118033988749895*fl[4]-0.8660254037844386*(fr[3]+fl[3])+0.8660254037844386*(fr[2]+fl[2])+0.5*fr[1]-0.5*(fl[1]+fr[0])+0.5*fl[0])*sgn(alphaDrSurf[0]-alphaDrSurf[1]); 
  fUpwindQuad[1] = 0.5*(1.118033988749895*(fr[4]+fl[4])-0.8660254037844386*(fr[3]+fr[2])+0.8660254037844386*(fl[3]+fl[2])+0.5*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-1.118033988749895*fr[4])+1.118033988749895*fl[4]+0.8660254037844386*(fr[3]+fl[3]+fr[2]+fl[2])-0.5*(fr[1]+fr[0])+0.5*(fl[1]+fl[0]))*sgn(alphaDrSurf[1]+alphaDrSurf[0]); 

  double fUpwind[2];
  fUpwind[0] = 0.7071067811865475*(fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.7071067811865475*(fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  double Gdiff[5]; 
  double Ghat[5]; 
  double incr2[5]; 

  incr2[2] = nuVtSqSum[0]*(0.2995357736356374*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[1]*((-0.430893194785552*fr[3])+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[3] = nuVtSqSum[1]*(0.2995357736356374*(fr[4]+fl[4])-0.430893194785552*fr[2]+0.430893194785552*fl[2]+0.3061862178478971*(fr[0]+fl[0]))+nuVtSqSum[0]*((-0.430893194785552*fr[3])+0.430893194785552*fl[3]+0.3061862178478971*(fr[1]+fl[1])); 
  incr2[4] = nuVtSqSum[0]*((-1.160097062884178*(fr[4]+fl[4]))+1.668842167398551*fr[2]-1.668842167398551*fl[2]-1.185854122563142*(fr[0]+fl[0]))+nuVtSqSum[1]*(1.668842167398551*fr[3]-1.668842167398551*fl[3]-1.185854122563142*(fr[1]+fl[1])); 

  Gdiff[0] = nuVtSqSum[0]*(1.897366596101028*fr[4]-1.897366596101028*fl[4]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+nuVtSqSum[1]*((-3.368048396326869*(fr[3]+fl[3]))+2.651650429449552*fr[1]-2.651650429449552*fl[1]); 
  Gdiff[1] = nuVtSqSum[1]*(1.897366596101028*fr[4]-1.897366596101028*fl[4]-3.368048396326869*(fr[2]+fl[2])+2.651650429449552*fr[0]-2.651650429449552*fl[0])+nuVtSqSum[0]*((-3.368048396326869*(fr[3]+fl[3]))+2.651650429449552*fr[1]-2.651650429449552*fl[1]); 

  Ghat[0] = Gdiff[0]*rdv+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]; 

  double incr1[5]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 
  incr1[4] = -1.118033988749895*Ghat[0]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += incr2[4]*rdvSq4L-1.0*incr1[4]*rdv2L; 

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
