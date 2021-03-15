#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf2x2vSer_Vpar_P1(const double m_, const double cflL, const double cflR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[8]; 
  alphaDrSurf[0] = 0.7071067811865475*((4.0*wl[2]+2.0*dxvl[2])*nuSum-2.0*nuUSum[0]); 
  alphaDrSurf[1] = -1.414213562373095*nuUSum[1]; 
  alphaDrSurf[2] = -1.414213562373095*nuUSum[2]; 
  alphaDrSurf[4] = -1.414213562373095*nuUSum[3]; 

  double fUpwindQuad[8];
  fUpwindQuad[0] = 0.5*(0.4330127018922193*fr[15]-0.4330127018922193*(fl[15]+fr[14]+fr[13])+0.4330127018922193*(fl[14]+fl[13])-0.25*(fr[12]+fl[12])-0.4330127018922193*fr[11]+0.4330127018922193*(fl[11]+fr[10])-0.4330127018922193*fl[10]+0.25*(fr[9]+fl[9]+fr[8]+fl[8])+0.4330127018922193*(fr[7]+fr[6])-0.4330127018922193*(fl[7]+fl[6])+0.25*(fr[5]+fl[5])-0.25*(fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]-0.25*(fr[2]+fl[2]+fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]))+0.4330127018922193*(fr[14]+fl[14]+fr[13]+fl[13])+0.25*fr[12]-0.25*fl[12]+0.4330127018922193*(fr[11]+fl[11])-0.4330127018922193*(fr[10]+fl[10])-0.25*(fr[9]+fr[8])+0.25*(fl[9]+fl[8])-0.4330127018922193*(fr[7]+fl[7]+fr[6]+fl[6])-0.25*fr[5]+0.25*(fl[5]+fr[4])-0.25*fl[4]+0.4330127018922193*(fr[3]+fl[3])+0.25*(fr[2]+fr[1])-0.25*(fl[2]+fl[1]+fr[0])+0.25*fl[0])*sgn(alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[1] = 0.5*((-0.4330127018922193*(fr[15]+fr[14]))+0.4330127018922193*(fl[15]+fl[14]+fr[13])-0.4330127018922193*fl[13]+0.25*(fr[12]+fl[12])+0.4330127018922193*(fr[11]+fr[10])-0.4330127018922193*(fl[11]+fl[10])+0.25*(fr[9]+fl[9])-0.25*(fr[8]+fl[8])+0.4330127018922193*fr[7]-0.4330127018922193*(fl[7]+fr[6])+0.4330127018922193*fl[6]-0.25*(fr[5]+fl[5]+fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]-0.25*(fr[2]+fl[2])+0.25*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14])-0.4330127018922193*(fr[13]+fl[13])-0.25*fr[12]+0.25*fl[12]-0.4330127018922193*(fr[11]+fl[11]+fr[10]+fl[10])-0.25*fr[9]+0.25*(fl[9]+fr[8])-0.25*fl[8]-0.4330127018922193*(fr[7]+fl[7])+0.4330127018922193*(fr[6]+fl[6])+0.25*(fr[5]+fr[4])-0.25*(fl[5]+fl[4])+0.4330127018922193*(fr[3]+fl[3])+0.25*fr[2]-0.25*(fl[2]+fr[1]+fr[0])+0.25*(fl[1]+fl[0]))*sgn((-alphaDrSurf[4])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[2] = 0.5*((-0.4330127018922193*fr[15])+0.4330127018922193*(fl[15]+fr[14])-0.4330127018922193*(fl[14]+fr[13])+0.4330127018922193*fl[13]+0.25*(fr[12]+fl[12])+0.4330127018922193*(fr[11]+fr[10])-0.4330127018922193*(fl[11]+fl[10])-0.25*(fr[9]+fl[9])+0.25*(fr[8]+fl[8])-0.4330127018922193*fr[7]+0.4330127018922193*(fl[7]+fr[6])-0.4330127018922193*fl[6]-0.25*(fr[5]+fl[5]+fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]+0.25*(fr[2]+fl[2])-0.25*(fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15])-0.4330127018922193*(fr[14]+fl[14])+0.4330127018922193*(fr[13]+fl[13])-0.25*fr[12]+0.25*fl[12]-0.4330127018922193*(fr[11]+fl[11]+fr[10]+fl[10])+0.25*fr[9]-0.25*(fl[9]+fr[8])+0.25*fl[8]+0.4330127018922193*(fr[7]+fl[7])-0.4330127018922193*(fr[6]+fl[6])+0.25*(fr[5]+fr[4])-0.25*(fl[5]+fl[4])+0.4330127018922193*(fr[3]+fl[3])-0.25*fr[2]+0.25*(fl[2]+fr[1])-0.25*(fl[1]+fr[0])+0.25*fl[0])*sgn((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[3] = 0.5*(0.4330127018922193*(fr[15]+fr[14]+fr[13])-0.4330127018922193*(fl[15]+fl[14]+fl[13])-0.25*(fr[12]+fl[12])-0.4330127018922193*fr[11]+0.4330127018922193*(fl[11]+fr[10])-0.4330127018922193*fl[10]-0.25*(fr[9]+fl[9]+fr[8]+fl[8])-0.4330127018922193*(fr[7]+fr[6])+0.4330127018922193*(fl[7]+fl[6])+0.25*(fr[5]+fl[5])-0.25*(fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]+0.25*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]))+0.25*fr[12]-0.25*fl[12]+0.4330127018922193*(fr[11]+fl[11])-0.4330127018922193*(fr[10]+fl[10])+0.25*(fr[9]+fr[8])-0.25*(fl[9]+fl[8])+0.4330127018922193*(fr[7]+fl[7]+fr[6]+fl[6])-0.25*fr[5]+0.25*(fl[5]+fr[4])-0.25*fl[4]+0.4330127018922193*(fr[3]+fl[3])-0.25*(fr[2]+fr[1]+fr[0])+0.25*(fl[2]+fl[1]+fl[0]))*sgn(alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[4] = 0.5*((-0.4330127018922193*fr[15])+0.4330127018922193*(fl[15]+fr[14]+fr[13])-0.4330127018922193*(fl[14]+fl[13])+0.25*(fr[12]+fl[12])-0.4330127018922193*(fr[11]+fr[10])+0.4330127018922193*(fl[11]+fl[10])-0.25*(fr[9]+fl[9]+fr[8]+fl[8])+0.4330127018922193*(fr[7]+fr[6])-0.4330127018922193*(fl[7]+fl[6])+0.25*(fr[5]+fl[5]+fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]-0.25*(fr[2]+fl[2]+fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15])-0.4330127018922193*(fr[14]+fl[14]+fr[13]+fl[13])-0.25*fr[12]+0.25*fl[12]+0.4330127018922193*(fr[11]+fl[11]+fr[10]+fl[10])+0.25*(fr[9]+fr[8])-0.25*(fl[9]+fl[8])-0.4330127018922193*(fr[7]+fl[7]+fr[6]+fl[6])-0.25*(fr[5]+fr[4])+0.25*(fl[5]+fl[4])+0.4330127018922193*(fr[3]+fl[3])+0.25*(fr[2]+fr[1])-0.25*(fl[2]+fl[1]+fr[0])+0.25*fl[0])*sgn(alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[5] = 0.5*(0.4330127018922193*(fr[15]+fr[14])-0.4330127018922193*(fl[15]+fl[14]+fr[13])+0.4330127018922193*fl[13]-0.25*(fr[12]+fl[12])+0.4330127018922193*fr[11]-0.4330127018922193*(fl[11]+fr[10])+0.4330127018922193*fl[10]-0.25*(fr[9]+fl[9])+0.25*(fr[8]+fl[8])+0.4330127018922193*fr[7]-0.4330127018922193*(fl[7]+fr[6])+0.4330127018922193*fl[6]-0.25*(fr[5]+fl[5])+0.25*(fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]-0.25*(fr[2]+fl[2])+0.25*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]))+0.4330127018922193*(fr[13]+fl[13])+0.25*fr[12]-0.25*fl[12]-0.4330127018922193*(fr[11]+fl[11])+0.4330127018922193*(fr[10]+fl[10])+0.25*fr[9]-0.25*(fl[9]+fr[8])+0.25*fl[8]-0.4330127018922193*(fr[7]+fl[7])+0.4330127018922193*(fr[6]+fl[6])+0.25*fr[5]-0.25*(fl[5]+fr[4])+0.25*fl[4]+0.4330127018922193*(fr[3]+fl[3])+0.25*fr[2]-0.25*(fl[2]+fr[1]+fr[0])+0.25*(fl[1]+fl[0]))*sgn((-alphaDrSurf[4])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[6] = 0.5*(0.4330127018922193*fr[15]-0.4330127018922193*(fl[15]+fr[14])+0.4330127018922193*(fl[14]+fr[13])-0.4330127018922193*fl[13]-0.25*(fr[12]+fl[12])+0.4330127018922193*fr[11]-0.4330127018922193*(fl[11]+fr[10])+0.4330127018922193*fl[10]+0.25*(fr[9]+fl[9])-0.25*(fr[8]+fl[8])-0.4330127018922193*fr[7]+0.4330127018922193*(fl[7]+fr[6])-0.4330127018922193*fl[6]-0.25*(fr[5]+fl[5])+0.25*(fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]+0.25*(fr[2]+fl[2])-0.25*(fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]))+0.4330127018922193*(fr[14]+fl[14])-0.4330127018922193*(fr[13]+fl[13])+0.25*fr[12]-0.25*fl[12]-0.4330127018922193*(fr[11]+fl[11])+0.4330127018922193*(fr[10]+fl[10])-0.25*fr[9]+0.25*(fl[9]+fr[8])-0.25*fl[8]+0.4330127018922193*(fr[7]+fl[7])-0.4330127018922193*(fr[6]+fl[6])+0.25*fr[5]-0.25*(fl[5]+fr[4])+0.25*fl[4]+0.4330127018922193*(fr[3]+fl[3])-0.25*fr[2]+0.25*(fl[2]+fr[1])-0.25*(fl[1]+fr[0])+0.25*fl[0])*sgn((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[7] = 0.5*((-0.4330127018922193*(fr[15]+fr[14]+fr[13]))+0.4330127018922193*(fl[15]+fl[14]+fl[13])+0.25*(fr[12]+fl[12])-0.4330127018922193*(fr[11]+fr[10])+0.4330127018922193*(fl[11]+fl[10])+0.25*(fr[9]+fl[9]+fr[8]+fl[8])-0.4330127018922193*(fr[7]+fr[6])+0.4330127018922193*(fl[7]+fl[6])+0.25*(fr[5]+fl[5]+fr[4]+fl[4])-0.4330127018922193*fr[3]+0.4330127018922193*fl[3]+0.25*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13])-0.25*fr[12]+0.25*fl[12]+0.4330127018922193*(fr[11]+fl[11]+fr[10]+fl[10])-0.25*(fr[9]+fr[8])+0.25*(fl[9]+fl[8])+0.4330127018922193*(fr[7]+fl[7]+fr[6]+fl[6])-0.25*(fr[5]+fr[4])+0.25*(fl[5]+fl[4])+0.4330127018922193*(fr[3]+fl[3])-0.25*(fr[2]+fr[1]+fr[0])+0.25*(fl[2]+fl[1]+fl[0]))*sgn(alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]); 

  double fUpwind[8];
  fUpwind[0] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[5] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[6] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  double Gdiff[16]; 
  double Ghat[16]; 
  double incr2[16]; 

  incr2[3] = nuVtSqSum[3]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[2]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[6] = nuVtSqSum[2]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[3]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[7] = nuVtSqSum[1]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[0]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[10] = nuVtSqSum[3]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[2]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[1]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[0]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[11] = nuVtSqSum[0]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[1]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[2]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0]))*nuVtSqSum[3]; 
  incr2[13] = nuVtSqSum[2]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[3]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[0]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[1]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[14] = nuVtSqSum[1]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[0]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[3]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[2]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[15] = nuVtSqSum[0]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[1]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[2]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[3]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 

  Gdiff[0] = nuVtSqSum[3]*((-0.5412658773652741*(fr[11]+fl[11]))+0.5625*fr[5]-0.5625*fl[5])+nuVtSqSum[2]*((-0.5412658773652741*(fr[7]+fl[7]))+0.5625*fr[2]-0.5625*fl[2])+nuVtSqSum[1]*((-0.5412658773652741*(fr[6]+fl[6]))+0.5625*fr[1]-0.5625*fl[1])+nuVtSqSum[0]*((-0.5412658773652741*(fr[3]+fl[3]))+0.5625*fr[0]-0.5625*fl[0]); 
  Gdiff[1] = nuVtSqSum[2]*((-0.5412658773652741*(fr[11]+fl[11]))+0.5625*fr[5]-0.5625*fl[5])+nuVtSqSum[3]*((-0.5412658773652741*(fr[7]+fl[7]))+0.5625*fr[2]-0.5625*fl[2])+nuVtSqSum[0]*((-0.5412658773652741*(fr[6]+fl[6]))+0.5625*fr[1]-0.5625*fl[1])+nuVtSqSum[1]*((-0.5412658773652741*(fr[3]+fl[3]))+0.5625*fr[0]-0.5625*fl[0]); 
  Gdiff[2] = nuVtSqSum[1]*((-0.5412658773652741*(fr[11]+fl[11]))+0.5625*fr[5]-0.5625*fl[5])+nuVtSqSum[0]*((-0.5412658773652741*(fr[7]+fl[7]))+0.5625*fr[2]-0.5625*fl[2])+nuVtSqSum[3]*((-0.5412658773652741*(fr[6]+fl[6]))+0.5625*fr[1]-0.5625*fl[1])+nuVtSqSum[2]*((-0.5412658773652741*(fr[3]+fl[3]))+0.5625*fr[0]-0.5625*fl[0]); 
  Gdiff[4] = nuVtSqSum[3]*((-0.5412658773652741*(fr[15]+fl[15]))+0.5625*fr[12]-0.5625*fl[12])+nuVtSqSum[2]*((-0.5412658773652741*(fr[14]+fl[14]))+0.5625*fr[9]-0.5625*fl[9])+nuVtSqSum[1]*((-0.5412658773652741*(fr[13]+fl[13]))+0.5625*fr[8]-0.5625*fl[8])+nuVtSqSum[0]*((-0.5412658773652741*(fr[10]+fl[10]))+0.5625*fr[4]-0.5625*fl[4]); 
  Gdiff[5] = nuVtSqSum[0]*((-0.5412658773652741*(fr[11]+fl[11]))+0.5625*fr[5]-0.5625*fl[5])+nuVtSqSum[1]*((-0.5412658773652741*(fr[7]+fl[7]))+0.5625*fr[2]-0.5625*fl[2])+nuVtSqSum[2]*((-0.5412658773652741*(fr[6]+fl[6]))+0.5625*fr[1]-0.5625*fl[1])+((-0.5412658773652741*(fr[3]+fl[3]))+0.5625*fr[0]-0.5625*fl[0])*nuVtSqSum[3]; 
  Gdiff[8] = nuVtSqSum[2]*((-0.5412658773652741*(fr[15]+fl[15]))+0.5625*fr[12]-0.5625*fl[12])+nuVtSqSum[3]*((-0.5412658773652741*(fr[14]+fl[14]))+0.5625*fr[9]-0.5625*fl[9])+nuVtSqSum[0]*((-0.5412658773652741*(fr[13]+fl[13]))+0.5625*fr[8]-0.5625*fl[8])+nuVtSqSum[1]*((-0.5412658773652741*(fr[10]+fl[10]))+0.5625*fr[4]-0.5625*fl[4]); 
  Gdiff[9] = nuVtSqSum[1]*((-0.5412658773652741*(fr[15]+fl[15]))+0.5625*fr[12]-0.5625*fl[12])+nuVtSqSum[0]*((-0.5412658773652741*(fr[14]+fl[14]))+0.5625*fr[9]-0.5625*fl[9])+nuVtSqSum[3]*((-0.5412658773652741*(fr[13]+fl[13]))+0.5625*fr[8]-0.5625*fl[8])+nuVtSqSum[2]*((-0.5412658773652741*(fr[10]+fl[10]))+0.5625*fr[4]-0.5625*fl[4]); 
  Gdiff[12] = nuVtSqSum[0]*((-0.5412658773652741*(fr[15]+fl[15]))+0.5625*fr[12]-0.5625*fl[12])+nuVtSqSum[1]*((-0.5412658773652741*(fr[14]+fl[14]))+0.5625*fr[9]-0.5625*fl[9])+nuVtSqSum[2]*((-0.5412658773652741*(fr[13]+fl[13]))+0.5625*fr[8]-0.5625*fl[8])+nuVtSqSum[3]*((-0.5412658773652741*(fr[10]+fl[10]))+0.5625*fr[4]-0.5625*fl[4]); 

  Ghat[0] = Gdiff[0]*rdv2L+0.5*(alphaDrSurf[4]*fUpwind[4]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = Gdiff[1]*rdv2L+0.5*(alphaDrSurf[2]*fUpwind[4]+fUpwind[2]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = Gdiff[2]*rdv2L+0.5*(alphaDrSurf[1]*fUpwind[4]+fUpwind[1]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[4] = Gdiff[4]*rdv2L+0.5*(alphaDrSurf[4]*fUpwind[7]+alphaDrSurf[2]*fUpwind[6]+alphaDrSurf[1]*fUpwind[5]+alphaDrSurf[0]*fUpwind[3]); 
  Ghat[5] = Gdiff[5]*rdv2L+0.5*(alphaDrSurf[0]*fUpwind[4]+fUpwind[0]*alphaDrSurf[4]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[8] = Gdiff[8]*rdv2L+0.5*(alphaDrSurf[2]*fUpwind[7]+alphaDrSurf[4]*fUpwind[6]+alphaDrSurf[0]*fUpwind[5]+alphaDrSurf[1]*fUpwind[3]); 
  Ghat[9] = Gdiff[9]*rdv2L+0.5*(alphaDrSurf[1]*fUpwind[7]+alphaDrSurf[0]*fUpwind[6]+alphaDrSurf[4]*fUpwind[5]+alphaDrSurf[2]*fUpwind[3]); 
  Ghat[12] = Gdiff[12]*rdv2L+0.5*(alphaDrSurf[0]*fUpwind[7]+alphaDrSurf[1]*fUpwind[6]+alphaDrSurf[2]*fUpwind[5]+fUpwind[3]*alphaDrSurf[4]); 

  double incr1[16]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[1]; 
  incr1[7] = 0.8660254037844386*Ghat[2]; 
  incr1[8] = -0.5*Ghat[8]; 
  incr1[9] = -0.5*Ghat[9]; 
  incr1[10] = 0.8660254037844386*Ghat[4]; 
  incr1[11] = 0.8660254037844386*Ghat[5]; 
  incr1[12] = -0.5*Ghat[12]; 
  incr1[13] = 0.8660254037844386*Ghat[8]; 
  incr1[14] = 0.8660254037844386*Ghat[9]; 
  incr1[15] = 0.8660254037844386*Ghat[12]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 
  outr[4] += incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr2[6]*rdvSq4R+incr1[6]*rdv2R; 
  outr[7] += incr2[7]*rdvSq4R+incr1[7]*rdv2R; 
  outr[8] += incr1[8]*rdv2R; 
  outr[9] += incr1[9]*rdv2R; 
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr2[11]*rdvSq4R+incr1[11]*rdv2R; 
  outr[12] += incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 
  outl[4] += -1.0*incr1[4]*rdv2L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += incr1[6]*rdv2L-1.0*incr2[6]*rdvSq4L; 
  outl[7] += incr1[7]*rdv2L-1.0*incr2[7]*rdvSq4L; 
  outl[8] += -1.0*incr1[8]*rdv2L; 
  outl[9] += -1.0*incr1[9]*rdv2L; 
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += incr1[11]*rdv2L-1.0*incr2[11]*rdvSq4L; 
  outl[12] += -1.0*incr1[12]*rdv2L; 
  outl[13] += incr1[13]*rdv2L-1.0*incr2[13]*rdvSq4L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr1[15]*rdv2L-1.0*incr2[15]*rdvSq4L; 

  return std::abs(wl[2]-(0.5*nuUSum[0])/nuSum); 
} 
double GkLBOconstNuSurf2x2vSer_Mu_P1(const double m_, const double cflL, const double cflR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[8]; 
  alphaDrSurf[0] = (5.656854249492382*wl[3]+2.828427124746191*dxvl[3])*nuSum; 

  double fUpwindQuad[8];
  fUpwindQuad[0] = 0.5*(0.4330127018922193*fr[15]-0.4330127018922193*(fl[15]+fr[14]+fr[13]+fr[12])+0.4330127018922193*(fl[14]+fl[13]+fl[12])-0.25*(fr[11]+fl[11])+0.4330127018922193*(fr[10]+fr[9]+fr[8])-0.4330127018922193*(fl[10]+fl[9]+fl[8])+0.25*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3]+fr[2]+fl[2]+fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]))+0.4330127018922193*(fr[14]+fl[14]+fr[13]+fl[13]+fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9]+fr[8]+fl[8])-0.25*(fr[7]+fr[6]+fr[5])+0.25*(fl[7]+fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])+0.25*(fr[3]+fr[2]+fr[1])-0.25*(fl[3]+fl[2]+fl[1]+fr[0])+0.25*fl[0])*sgn(alphaDrSurf[0]); 
  fUpwindQuad[1] = 0.5*((-0.4330127018922193*(fr[15]+fr[14]))+0.4330127018922193*(fl[15]+fl[14]+fr[13]+fr[12])-0.4330127018922193*(fl[13]+fl[12])+0.25*(fr[11]+fl[11])+0.4330127018922193*(fr[10]+fr[9])-0.4330127018922193*(fl[10]+fl[9]+fr[8])+0.4330127018922193*fl[8]+0.25*(fr[7]+fl[7])-0.25*(fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3]+fr[2]+fl[2])+0.25*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14])-0.4330127018922193*(fr[13]+fl[13]+fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9])+0.4330127018922193*(fr[8]+fl[8])-0.25*fr[7]+0.25*(fl[7]+fr[6]+fr[5])-0.25*(fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])+0.25*(fr[3]+fr[2])-0.25*(fl[3]+fl[2]+fr[1]+fr[0])+0.25*(fl[1]+fl[0]))*sgn(alphaDrSurf[0]); 
  fUpwindQuad[2] = 0.5*((-0.4330127018922193*fr[15])+0.4330127018922193*(fl[15]+fr[14])-0.4330127018922193*(fl[14]+fr[13])+0.4330127018922193*(fl[13]+fr[12])-0.4330127018922193*fl[12]+0.25*(fr[11]+fl[11])+0.4330127018922193*fr[10]-0.4330127018922193*(fl[10]+fr[9])+0.4330127018922193*(fl[9]+fr[8])-0.4330127018922193*fl[8]-0.25*(fr[7]+fl[7])+0.25*(fr[6]+fl[6])-0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3])+0.25*(fr[2]+fl[2])-0.25*(fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15])-0.4330127018922193*(fr[14]+fl[14])+0.4330127018922193*(fr[13]+fl[13])-0.4330127018922193*(fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10])+0.4330127018922193*(fr[9]+fl[9])-0.4330127018922193*(fr[8]+fl[8])+0.25*fr[7]-0.25*(fl[7]+fr[6])+0.25*(fl[6]+fr[5])-0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])+0.25*fr[3]-0.25*(fl[3]+fr[2])+0.25*(fl[2]+fr[1])-0.25*(fl[1]+fr[0])+0.25*fl[0])*sgn(alphaDrSurf[0]); 
  fUpwindQuad[3] = 0.5*(0.4330127018922193*(fr[15]+fr[14]+fr[13])-0.4330127018922193*(fl[15]+fl[14]+fl[13]+fr[12])+0.4330127018922193*fl[12]-0.25*(fr[11]+fl[11])+0.4330127018922193*fr[10]-0.4330127018922193*(fl[10]+fr[9]+fr[8])+0.4330127018922193*(fl[9]+fl[8])-0.25*(fr[7]+fl[7]+fr[6]+fl[6])+0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3])+0.25*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]))+0.4330127018922193*(fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10])+0.4330127018922193*(fr[9]+fl[9]+fr[8]+fl[8])+0.25*(fr[7]+fr[6])-0.25*(fl[7]+fl[6]+fr[5])+0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])+0.25*fr[3]-0.25*(fl[3]+fr[2]+fr[1]+fr[0])+0.25*(fl[2]+fl[1]+fl[0]))*sgn(alphaDrSurf[0]); 
  fUpwindQuad[4] = 0.5*((-0.4330127018922193*fr[15])+0.4330127018922193*(fl[15]+fr[14]+fr[13])-0.4330127018922193*(fl[14]+fl[13]+fr[12])+0.4330127018922193*fl[12]+0.25*(fr[11]+fl[11])-0.4330127018922193*fr[10]+0.4330127018922193*(fl[10]+fr[9]+fr[8])-0.4330127018922193*(fl[9]+fl[8])-0.25*(fr[7]+fl[7]+fr[6]+fl[6])+0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3])-0.25*(fr[2]+fl[2]+fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15])-0.4330127018922193*(fr[14]+fl[14]+fr[13]+fl[13])+0.4330127018922193*(fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10])-0.4330127018922193*(fr[9]+fl[9]+fr[8]+fl[8])+0.25*(fr[7]+fr[6])-0.25*(fl[7]+fl[6]+fr[5])+0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])-0.25*fr[3]+0.25*(fl[3]+fr[2]+fr[1])-0.25*(fl[2]+fl[1]+fr[0])+0.25*fl[0])*sgn(alphaDrSurf[0]); 
  fUpwindQuad[5] = 0.5*(0.4330127018922193*(fr[15]+fr[14])-0.4330127018922193*(fl[15]+fl[14]+fr[13])+0.4330127018922193*(fl[13]+fr[12])-0.4330127018922193*fl[12]-0.25*(fr[11]+fl[11])-0.4330127018922193*fr[10]+0.4330127018922193*(fl[10]+fr[9])-0.4330127018922193*(fl[9]+fr[8])+0.4330127018922193*fl[8]-0.25*(fr[7]+fl[7])+0.25*(fr[6]+fl[6])-0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3])-0.25*(fr[2]+fl[2])+0.25*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]))+0.4330127018922193*(fr[13]+fl[13])-0.4330127018922193*(fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10])-0.4330127018922193*(fr[9]+fl[9])+0.4330127018922193*(fr[8]+fl[8])+0.25*fr[7]-0.25*(fl[7]+fr[6])+0.25*(fl[6]+fr[5])-0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])-0.25*fr[3]+0.25*(fl[3]+fr[2])-0.25*(fl[2]+fr[1]+fr[0])+0.25*(fl[1]+fl[0]))*sgn(alphaDrSurf[0]); 
  fUpwindQuad[6] = 0.5*(0.4330127018922193*fr[15]-0.4330127018922193*(fl[15]+fr[14])+0.4330127018922193*(fl[14]+fr[13]+fr[12])-0.4330127018922193*(fl[13]+fl[12])-0.25*(fr[11]+fl[11])-0.4330127018922193*(fr[10]+fr[9])+0.4330127018922193*(fl[10]+fl[9]+fr[8])-0.4330127018922193*fl[8]+0.25*(fr[7]+fl[7])-0.25*(fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3]+fr[2]+fl[2])-0.25*(fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]))+0.4330127018922193*(fr[14]+fl[14])-0.4330127018922193*(fr[13]+fl[13]+fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9])-0.4330127018922193*(fr[8]+fl[8])-0.25*fr[7]+0.25*(fl[7]+fr[6]+fr[5])-0.25*(fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])-0.25*(fr[3]+fr[2])+0.25*(fl[3]+fl[2]+fr[1])-0.25*(fl[1]+fr[0])+0.25*fl[0])*sgn(alphaDrSurf[0]); 
  fUpwindQuad[7] = 0.5*((-0.4330127018922193*(fr[15]+fr[14]+fr[13]+fr[12]))+0.4330127018922193*(fl[15]+fl[14]+fl[13]+fl[12])+0.25*(fr[11]+fl[11])-0.4330127018922193*(fr[10]+fr[9]+fr[8])+0.4330127018922193*(fl[10]+fl[9]+fl[8])+0.25*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3]+fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]+fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9]+fr[8]+fl[8])-0.25*(fr[7]+fr[6]+fr[5])+0.25*(fl[7]+fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])-0.25*(fr[3]+fr[2]+fr[1]+fr[0])+0.25*(fl[3]+fl[2]+fl[1]+fl[0]))*sgn(alphaDrSurf[0]); 

  double fUpwind[8];
  fUpwind[0] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[5] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[6] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  double diffFac[4]; 
  diffFac[0] = (BmagInv[3]*nuVtSqSum[3]+BmagInv[2]*nuVtSqSum[2]+BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*(wl[3]+0.5*dxvl[3])*m_; 
  diffFac[1] = (BmagInv[2]*nuVtSqSum[3]+nuVtSqSum[2]*BmagInv[3]+BmagInv[0]*nuVtSqSum[1]+nuVtSqSum[0]*BmagInv[1])*(wl[3]+0.5*dxvl[3])*m_; 
  diffFac[2] = (BmagInv[1]*nuVtSqSum[3]+nuVtSqSum[1]*BmagInv[3]+BmagInv[0]*nuVtSqSum[2]+nuVtSqSum[0]*BmagInv[2])*(wl[3]+0.5*dxvl[3])*m_; 
  diffFac[3] = (BmagInv[0]*nuVtSqSum[3]+nuVtSqSum[0]*BmagInv[3]+BmagInv[1]*nuVtSqSum[2]+nuVtSqSum[1]*BmagInv[2])*(wl[3]+0.5*dxvl[3])*m_; 

  double Gdiff[16]; 
  double Ghat[16]; 
  double incr2[16]; 

  incr2[4] = diffFac[3]*(0.25*fl[12]-0.25*fr[12])+diffFac[2]*(0.25*fl[9]-0.25*fr[9])+diffFac[1]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[3]*(fr[5]+fl[5])+diffFac[0]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*(diffFac[2]*(fr[2]+fl[2])+diffFac[1]*(fr[1]+fl[1])+diffFac[0]*(fr[0]+fl[0])); 
  incr2[8] = diffFac[2]*(0.25*fl[12]-0.25*fr[12])+diffFac[3]*(0.25*fl[9]-0.25*fr[9])+diffFac[0]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[2]*(fr[5]+fl[5])+diffFac[1]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*((fr[2]+fl[2])*diffFac[3]+diffFac[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*diffFac[1]); 
  incr2[9] = diffFac[1]*(0.25*fl[12]-0.25*fr[12])+diffFac[0]*(0.25*fl[9]-0.25*fr[9])+diffFac[3]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[1]*(fr[5]+fl[5])+diffFac[2]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*((fr[1]+fl[1])*diffFac[3]+diffFac[0]*(fr[2]+fl[2])+(fr[0]+fl[0])*diffFac[2]); 
  incr2[10] = diffFac[3]*(0.25*fl[15]-0.25*fr[15])+diffFac[2]*(0.25*fl[14]-0.25*fr[14])+diffFac[1]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[3]*(fr[11]+fl[11])+diffFac[0]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[2]*(fr[7]+fl[7])+diffFac[1]*(fr[6]+fl[6])+diffFac[0]*(fr[3]+fl[3])); 
  incr2[12] = diffFac[0]*(0.25*fl[12]-0.25*fr[12])+diffFac[1]*(0.25*fl[9]-0.25*fr[9])+diffFac[2]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[0]*(fr[5]+fl[5])+diffFac[3]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*((fr[0]+fl[0])*diffFac[3]+diffFac[1]*(fr[2]+fl[2])+(fr[1]+fl[1])*diffFac[2]); 
  incr2[13] = diffFac[2]*(0.25*fl[15]-0.25*fr[15])+diffFac[3]*(0.25*fl[14]-0.25*fr[14])+diffFac[0]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[2]*(fr[11]+fl[11])+diffFac[1]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[3]*(fr[7]+fl[7])+diffFac[0]*(fr[6]+fl[6])+diffFac[1]*(fr[3]+fl[3])); 
  incr2[14] = diffFac[1]*(0.25*fl[15]-0.25*fr[15])+diffFac[0]*(0.25*fl[14]-0.25*fr[14])+diffFac[3]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[1]*(fr[11]+fl[11])+diffFac[2]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[0]*(fr[7]+fl[7])+diffFac[3]*(fr[6]+fl[6])+diffFac[2]*(fr[3]+fl[3])); 
  incr2[15] = diffFac[0]*(0.25*fl[15]-0.25*fr[15])+diffFac[1]*(0.25*fl[14]-0.25*fr[14])+diffFac[2]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[0]*(fr[11]+fl[11])+diffFac[3]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[1]*(fr[7]+fl[7])+diffFac[2]*(fr[6]+fl[6])+diffFac[3]*(fr[3]+fl[3])); 

  Gdiff[0] = (-0.5412658773652741*(diffFac[3]*(fr[12]+fl[12])+diffFac[2]*(fr[9]+fl[9])+diffFac[1]*(fr[8]+fl[8])))+diffFac[3]*(0.5625*fr[5]-0.5625*fl[5])-0.5412658773652741*diffFac[0]*(fr[4]+fl[4])+diffFac[2]*(0.5625*fr[2]-0.5625*fl[2])+diffFac[1]*(0.5625*fr[1]-0.5625*fl[1])+diffFac[0]*(0.5625*fr[0]-0.5625*fl[0]); 
  Gdiff[1] = (-0.5412658773652741*(diffFac[2]*(fr[12]+fl[12])+diffFac[3]*(fr[9]+fl[9])+diffFac[0]*(fr[8]+fl[8])))+diffFac[2]*(0.5625*fr[5]-0.5625*fl[5])-0.5412658773652741*diffFac[1]*(fr[4]+fl[4])+(0.5625*fr[2]-0.5625*fl[2])*diffFac[3]+diffFac[0]*(0.5625*fr[1]-0.5625*fl[1])+(0.5625*fr[0]-0.5625*fl[0])*diffFac[1]; 
  Gdiff[2] = (-0.5412658773652741*(diffFac[1]*(fr[12]+fl[12])+diffFac[0]*(fr[9]+fl[9])+diffFac[3]*(fr[8]+fl[8])))+diffFac[1]*(0.5625*fr[5]-0.5625*fl[5])-0.5412658773652741*diffFac[2]*(fr[4]+fl[4])+(0.5625*fr[1]-0.5625*fl[1])*diffFac[3]+diffFac[0]*(0.5625*fr[2]-0.5625*fl[2])+(0.5625*fr[0]-0.5625*fl[0])*diffFac[2]; 
  Gdiff[3] = (-0.5412658773652741*(diffFac[3]*(fr[15]+fl[15])+diffFac[2]*(fr[14]+fl[14])+diffFac[1]*(fr[13]+fl[13])))+diffFac[3]*(0.5625*fr[11]-0.5625*fl[11])-0.5412658773652741*diffFac[0]*(fr[10]+fl[10])+diffFac[2]*(0.5625*fr[7]-0.5625*fl[7])+diffFac[1]*(0.5625*fr[6]-0.5625*fl[6])+diffFac[0]*(0.5625*fr[3]-0.5625*fl[3]); 
  Gdiff[5] = (-0.5412658773652741*(diffFac[0]*(fr[12]+fl[12])+diffFac[1]*(fr[9]+fl[9])+diffFac[2]*(fr[8]+fl[8])))+diffFac[0]*(0.5625*fr[5]-0.5625*fl[5])+diffFac[3]*((-0.5412658773652741*(fr[4]+fl[4]))+0.5625*fr[0]-0.5625*fl[0])+diffFac[1]*(0.5625*fr[2]-0.5625*fl[2])+(0.5625*fr[1]-0.5625*fl[1])*diffFac[2]; 
  Gdiff[6] = (-0.5412658773652741*(diffFac[2]*(fr[15]+fl[15])+diffFac[3]*(fr[14]+fl[14])+diffFac[0]*(fr[13]+fl[13])))+diffFac[2]*(0.5625*fr[11]-0.5625*fl[11])-0.5412658773652741*diffFac[1]*(fr[10]+fl[10])+diffFac[3]*(0.5625*fr[7]-0.5625*fl[7])+diffFac[0]*(0.5625*fr[6]-0.5625*fl[6])+diffFac[1]*(0.5625*fr[3]-0.5625*fl[3]); 
  Gdiff[7] = (-0.5412658773652741*(diffFac[1]*(fr[15]+fl[15])+diffFac[0]*(fr[14]+fl[14])+diffFac[3]*(fr[13]+fl[13])))+diffFac[1]*(0.5625*fr[11]-0.5625*fl[11])-0.5412658773652741*diffFac[2]*(fr[10]+fl[10])+diffFac[0]*(0.5625*fr[7]-0.5625*fl[7])+diffFac[3]*(0.5625*fr[6]-0.5625*fl[6])+diffFac[2]*(0.5625*fr[3]-0.5625*fl[3]); 
  Gdiff[11] = (-0.5412658773652741*(diffFac[0]*(fr[15]+fl[15])+diffFac[1]*(fr[14]+fl[14])+diffFac[2]*(fr[13]+fl[13])))+diffFac[0]*(0.5625*fr[11]-0.5625*fl[11])-0.5412658773652741*diffFac[3]*(fr[10]+fl[10])+diffFac[1]*(0.5625*fr[7]-0.5625*fl[7])+diffFac[2]*(0.5625*fr[6]-0.5625*fl[6])+diffFac[3]*(0.5625*fr[3]-0.5625*fl[3]); 

  Ghat[0] = Gdiff[0]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[1]; 
  Ghat[2] = Gdiff[2]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[3] = Gdiff[3]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[5] = Gdiff[5]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[4]; 
  Ghat[6] = Gdiff[6]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[5]; 
  Ghat[7] = Gdiff[7]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[6]; 
  Ghat[11] = Gdiff[11]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[7]; 

  double incr1[16]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[0]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = -0.5*Ghat[6]; 
  incr1[7] = -0.5*Ghat[7]; 
  incr1[8] = 0.8660254037844386*Ghat[1]; 
  incr1[9] = 0.8660254037844386*Ghat[2]; 
  incr1[10] = 0.8660254037844386*Ghat[3]; 
  incr1[11] = -0.5*Ghat[11]; 
  incr1[12] = 0.8660254037844386*Ghat[5]; 
  incr1[13] = 0.8660254037844386*Ghat[6]; 
  incr1[14] = 0.8660254037844386*Ghat[7]; 
  incr1[15] = 0.8660254037844386*Ghat[11]; 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr1[2]*rdv2R; 
  outr[3] += incr1[3]*rdv2R; 
  outr[4] += incr2[4]*rdvSq4R+incr1[4]*rdv2R; 
  outr[5] += incr1[5]*rdv2R; 
  outr[6] += incr1[6]*rdv2R; 
  outr[7] += incr1[7]*rdv2R; 
  outr[8] += incr2[8]*rdvSq4R+incr1[8]*rdv2R; 
  outr[9] += incr2[9]*rdvSq4R+incr1[9]*rdv2R; 
  outr[10] += incr2[10]*rdvSq4R+incr1[10]*rdv2R; 
  outr[11] += incr1[11]*rdv2R; 
  outr[12] += incr2[12]*rdvSq4R+incr1[12]*rdv2R; 
  outr[13] += incr2[13]*rdvSq4R+incr1[13]*rdv2R; 
  outr[14] += incr2[14]*rdvSq4R+incr1[14]*rdv2R; 
  outr[15] += incr2[15]*rdvSq4R+incr1[15]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += -1.0*incr1[2]*rdv2L; 
  outl[3] += -1.0*incr1[3]*rdv2L; 
  outl[4] += incr1[4]*rdv2L-1.0*incr2[4]*rdvSq4L; 
  outl[5] += -1.0*incr1[5]*rdv2L; 
  outl[6] += -1.0*incr1[6]*rdv2L; 
  outl[7] += -1.0*incr1[7]*rdv2L; 
  outl[8] += incr1[8]*rdv2L-1.0*incr2[8]*rdvSq4L; 
  outl[9] += incr1[9]*rdv2L-1.0*incr2[9]*rdvSq4L; 
  outl[10] += incr1[10]*rdv2L-1.0*incr2[10]*rdvSq4L; 
  outl[11] += -1.0*incr1[11]*rdv2L; 
  outl[12] += incr1[12]*rdv2L-1.0*incr2[12]*rdvSq4L; 
  outl[13] += incr1[13]*rdv2L-1.0*incr2[13]*rdvSq4L; 
  outl[14] += incr1[14]*rdv2L-1.0*incr2[14]*rdvSq4L; 
  outl[15] += incr1[15]*rdv2L-1.0*incr2[15]*rdvSq4L; 

  return std::abs(2.0*wl[3]); 
} 
