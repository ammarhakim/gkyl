#include <GkLBOModDecl.h> 
double GkLBOconstNuSurf2x2vSer_Vpar_P1(const double m_, const double *cflRateByDirL, const double *cflRateByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // m_:              species mass. 
  // cflRateByDir[4]: CFL rate in each direction. 
  // w[4]:            Cell-center coordinates. 
  // dxv[4]:          Cell spacing. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:       maximum midpoint value of v-u. 
  // nuUSum[4]:       sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:    sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:           Distribution function in left/right cells 
  // outl/outr:       Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[2]; 
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

  incr2[3] = (nuVtSqSum[3]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[2]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])))*rdvSq4R; 
  incr2[6] = (nuVtSqSum[2]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[3]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])))*rdvSq4R; 
  incr2[7] = (nuVtSqSum[1]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[0]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])))*rdvSq4R; 
  incr2[10] = (nuVtSqSum[3]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[2]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[1]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[0]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])))*rdvSq4R; 
  incr2[11] = (nuVtSqSum[0]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[1]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[2]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0]))*nuVtSqSum[3])*rdvSq4R; 
  incr2[13] = (nuVtSqSum[2]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[3]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[0]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[1]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])))*rdvSq4R; 
  incr2[14] = (nuVtSqSum[1]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[0]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[3]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[2]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])))*rdvSq4R; 
  incr2[15] = (nuVtSqSum[0]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[1]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[2]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[3]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])))*rdvSq4R; 

  Gdiff[0] = nuVtSqSum[3]*((-1.082531754730548*(fr[11]+fl[11]))+1.125*fr[5]-1.125*fl[5])+nuVtSqSum[2]*((-1.082531754730548*(fr[7]+fl[7]))+1.125*fr[2]-1.125*fl[2])+nuVtSqSum[1]*((-1.082531754730548*(fr[6]+fl[6]))+1.125*fr[1]-1.125*fl[1])+nuVtSqSum[0]*((-1.082531754730548*(fr[3]+fl[3]))+1.125*fr[0]-1.125*fl[0]); 
  Gdiff[1] = nuVtSqSum[2]*((-1.082531754730548*(fr[11]+fl[11]))+1.125*fr[5]-1.125*fl[5])+nuVtSqSum[3]*((-1.082531754730548*(fr[7]+fl[7]))+1.125*fr[2]-1.125*fl[2])+nuVtSqSum[0]*((-1.082531754730548*(fr[6]+fl[6]))+1.125*fr[1]-1.125*fl[1])+nuVtSqSum[1]*((-1.082531754730548*(fr[3]+fl[3]))+1.125*fr[0]-1.125*fl[0]); 
  Gdiff[2] = nuVtSqSum[1]*((-1.082531754730548*(fr[11]+fl[11]))+1.125*fr[5]-1.125*fl[5])+nuVtSqSum[0]*((-1.082531754730548*(fr[7]+fl[7]))+1.125*fr[2]-1.125*fl[2])+nuVtSqSum[3]*((-1.082531754730548*(fr[6]+fl[6]))+1.125*fr[1]-1.125*fl[1])+nuVtSqSum[2]*((-1.082531754730548*(fr[3]+fl[3]))+1.125*fr[0]-1.125*fl[0]); 
  Gdiff[4] = nuVtSqSum[3]*((-1.082531754730548*(fr[15]+fl[15]))+1.125*fr[12]-1.125*fl[12])+nuVtSqSum[2]*((-1.082531754730548*(fr[14]+fl[14]))+1.125*fr[9]-1.125*fl[9])+nuVtSqSum[1]*((-1.082531754730548*(fr[13]+fl[13]))+1.125*fr[8]-1.125*fl[8])+nuVtSqSum[0]*((-1.082531754730548*(fr[10]+fl[10]))+1.125*fr[4]-1.125*fl[4]); 
  Gdiff[5] = nuVtSqSum[0]*((-1.082531754730548*(fr[11]+fl[11]))+1.125*fr[5]-1.125*fl[5])+nuVtSqSum[1]*((-1.082531754730548*(fr[7]+fl[7]))+1.125*fr[2]-1.125*fl[2])+nuVtSqSum[2]*((-1.082531754730548*(fr[6]+fl[6]))+1.125*fr[1]-1.125*fl[1])+((-1.082531754730548*(fr[3]+fl[3]))+1.125*fr[0]-1.125*fl[0])*nuVtSqSum[3]; 
  Gdiff[8] = nuVtSqSum[2]*((-1.082531754730548*(fr[15]+fl[15]))+1.125*fr[12]-1.125*fl[12])+nuVtSqSum[3]*((-1.082531754730548*(fr[14]+fl[14]))+1.125*fr[9]-1.125*fl[9])+nuVtSqSum[0]*((-1.082531754730548*(fr[13]+fl[13]))+1.125*fr[8]-1.125*fl[8])+nuVtSqSum[1]*((-1.082531754730548*(fr[10]+fl[10]))+1.125*fr[4]-1.125*fl[4]); 
  Gdiff[9] = nuVtSqSum[1]*((-1.082531754730548*(fr[15]+fl[15]))+1.125*fr[12]-1.125*fl[12])+nuVtSqSum[0]*((-1.082531754730548*(fr[14]+fl[14]))+1.125*fr[9]-1.125*fl[9])+nuVtSqSum[3]*((-1.082531754730548*(fr[13]+fl[13]))+1.125*fr[8]-1.125*fl[8])+nuVtSqSum[2]*((-1.082531754730548*(fr[10]+fl[10]))+1.125*fr[4]-1.125*fl[4]); 
  Gdiff[12] = nuVtSqSum[0]*((-1.082531754730548*(fr[15]+fl[15]))+1.125*fr[12]-1.125*fl[12])+nuVtSqSum[1]*((-1.082531754730548*(fr[14]+fl[14]))+1.125*fr[9]-1.125*fl[9])+nuVtSqSum[2]*((-1.082531754730548*(fr[13]+fl[13]))+1.125*fr[8]-1.125*fl[8])+nuVtSqSum[3]*((-1.082531754730548*(fr[10]+fl[10]))+1.125*fr[4]-1.125*fl[4]); 

  Ghat[0] = Gdiff[0]*rdv+0.5*(alphaDrSurf[4]*fUpwind[4]+alphaDrSurf[2]*fUpwind[2]+alphaDrSurf[1]*fUpwind[1]+alphaDrSurf[0]*fUpwind[0]); 
  Ghat[1] = Gdiff[1]*rdv+0.5*(alphaDrSurf[2]*fUpwind[4]+fUpwind[2]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[1]+fUpwind[0]*alphaDrSurf[1]); 
  Ghat[2] = Gdiff[2]*rdv+0.5*(alphaDrSurf[1]*fUpwind[4]+fUpwind[1]*alphaDrSurf[4]+alphaDrSurf[0]*fUpwind[2]+fUpwind[0]*alphaDrSurf[2]); 
  Ghat[4] = Gdiff[4]*rdv+0.5*(alphaDrSurf[4]*fUpwind[7]+alphaDrSurf[2]*fUpwind[6]+alphaDrSurf[1]*fUpwind[5]+alphaDrSurf[0]*fUpwind[3]); 
  Ghat[5] = Gdiff[5]*rdv+0.5*(alphaDrSurf[0]*fUpwind[4]+fUpwind[0]*alphaDrSurf[4]+alphaDrSurf[1]*fUpwind[2]+fUpwind[1]*alphaDrSurf[2]); 
  Ghat[8] = Gdiff[8]*rdv+0.5*(alphaDrSurf[2]*fUpwind[7]+alphaDrSurf[4]*fUpwind[6]+alphaDrSurf[0]*fUpwind[5]+alphaDrSurf[1]*fUpwind[3]); 
  Ghat[9] = Gdiff[9]*rdv+0.5*(alphaDrSurf[1]*fUpwind[7]+alphaDrSurf[0]*fUpwind[6]+alphaDrSurf[4]*fUpwind[5]+alphaDrSurf[2]*fUpwind[3]); 
  Ghat[12] = Gdiff[12]*rdv+0.5*(alphaDrSurf[0]*fUpwind[7]+alphaDrSurf[1]*fUpwind[6]+alphaDrSurf[2]*fUpwind[5]+fUpwind[3]*alphaDrSurf[4]); 

  double incr1[16]; 
  incr1[0] = -0.5*Ghat[0]*rdv2R; 
  incr1[1] = -0.5*Ghat[1]*rdv2R; 
  incr1[2] = -0.5*Ghat[2]*rdv2R; 
  incr1[3] = 0.8660254037844386*Ghat[0]*rdv2R; 
  incr1[4] = -0.5*Ghat[4]*rdv2R; 
  incr1[5] = -0.5*Ghat[5]*rdv2R; 
  incr1[6] = 0.8660254037844386*Ghat[1]*rdv2R; 
  incr1[7] = 0.8660254037844386*Ghat[2]*rdv2R; 
  incr1[8] = -0.5*Ghat[8]*rdv2R; 
  incr1[9] = -0.5*Ghat[9]*rdv2R; 
  incr1[10] = 0.8660254037844386*Ghat[4]*rdv2R; 
  incr1[11] = 0.8660254037844386*Ghat[5]*rdv2R; 
  incr1[12] = -0.5*Ghat[12]*rdv2R; 
  incr1[13] = 0.8660254037844386*Ghat[8]*rdv2R; 
  incr1[14] = 0.8660254037844386*Ghat[9]*rdv2R; 
  incr1[15] = 0.8660254037844386*Ghat[12]*rdv2R; 

  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[8], outrPos[8]; 
  fluxFracL = cflRateByDirL[0] == 0. ? 0.25 : cflRateByDirL[1]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] == 0. ? 0.25 : cflRateByDirR[1]/cflRateByDirR[0]; 
  outlPos[0] = 0.02777777777777778*(incr2[15]-1.0*incr1[15]-1.732050807568877*incr2[14]+1.732050807568877*incr1[14]-1.732050807568877*incr2[13]+1.732050807568877*(incr1[13]+incr1[12])-1.732050807568877*incr2[11]+1.732050807568877*incr1[11]+3.0*incr2[10]-3.0*(incr1[10]+incr1[9]+incr1[8])+3.0*incr2[7]-3.0*incr1[7]+3.0*incr2[6]-3.0*(incr1[6]+incr1[5])+5.196152422706631*incr1[4]-5.196152422706631*incr2[3]+5.196152422706631*(incr1[3]+incr1[2]+incr1[1])-9.0*incr1[0]); 
  outlPos[1] = -0.02777777777777778*(incr2[15]-1.0*incr1[15]+1.732050807568877*incr2[14]-1.732050807568877*(incr1[14]+incr2[13])+1.732050807568877*(incr1[13]+incr1[12])-1.732050807568877*incr2[11]+1.732050807568877*incr1[11]-3.0*incr2[10]+3.0*(incr1[10]+incr1[9])-3.0*(incr1[8]+incr2[7])+3.0*(incr1[7]+incr2[6])-3.0*(incr1[6]+incr1[5])-5.196152422706631*incr1[4]+5.196152422706631*incr2[3]-5.196152422706631*(incr1[3]+incr1[2])+5.196152422706631*incr1[1]+9.0*incr1[0]); 
  outlPos[2] = -0.02777777777777778*(incr2[15]-1.0*incr1[15]-1.732050807568877*incr2[14]+1.732050807568877*(incr1[14]+incr2[13])-1.732050807568877*incr1[13]+1.732050807568877*incr1[12]-1.732050807568877*incr2[11]+1.732050807568877*incr1[11]-3.0*incr2[10]+3.0*incr1[10]-3.0*incr1[9]+3.0*(incr1[8]+incr2[7])-3.0*(incr1[7]+incr2[6])+3.0*incr1[6]-3.0*incr1[5]-5.196152422706631*incr1[4]+5.196152422706631*incr2[3]-5.196152422706631*incr1[3]+5.196152422706631*incr1[2]-5.196152422706631*incr1[1]+9.0*incr1[0]); 
  outlPos[3] = 0.02777777777777778*(incr2[15]-1.0*incr1[15]+1.732050807568877*incr2[14]-1.732050807568877*incr1[14]+1.732050807568877*incr2[13]-1.732050807568877*incr1[13]+1.732050807568877*incr1[12]-1.732050807568877*incr2[11]+1.732050807568877*incr1[11]+3.0*incr2[10]-3.0*incr1[10]+3.0*(incr1[9]+incr1[8])-3.0*incr2[7]+3.0*incr1[7]-3.0*incr2[6]+3.0*incr1[6]-3.0*incr1[5]+5.196152422706631*incr1[4]-5.196152422706631*incr2[3]+5.196152422706631*incr1[3]-5.196152422706631*(incr1[2]+incr1[1])-9.0*incr1[0]); 
  outlPos[4] = -0.02777777777777778*(incr2[15]-1.0*incr1[15]-1.732050807568877*incr2[14]+1.732050807568877*incr1[14]-1.732050807568877*incr2[13]+1.732050807568877*(incr1[13]+incr1[12]+incr2[11])-1.732050807568877*incr1[11]+3.0*incr2[10]-3.0*(incr1[10]+incr1[9]+incr1[8]+incr2[7])+3.0*incr1[7]-3.0*incr2[6]+3.0*(incr1[6]+incr1[5])+5.196152422706631*(incr1[4]+incr2[3])-5.196152422706631*(incr1[3]+incr1[2]+incr1[1])+9.0*incr1[0]); 
  outlPos[5] = 0.02777777777777778*(incr2[15]-1.0*incr1[15]+1.732050807568877*incr2[14]-1.732050807568877*(incr1[14]+incr2[13])+1.732050807568877*(incr1[13]+incr1[12]+incr2[11])-1.732050807568877*incr1[11]-3.0*incr2[10]+3.0*(incr1[10]+incr1[9])-3.0*incr1[8]+3.0*incr2[7]-3.0*(incr1[7]+incr2[6])+3.0*(incr1[6]+incr1[5])-5.196152422706631*(incr1[4]+incr2[3])+5.196152422706631*(incr1[3]+incr1[2])-5.196152422706631*incr1[1]-9.0*incr1[0]); 
  outlPos[6] = 0.02777777777777778*(incr2[15]-1.0*incr1[15]-1.732050807568877*incr2[14]+1.732050807568877*(incr1[14]+incr2[13])-1.732050807568877*incr1[13]+1.732050807568877*(incr1[12]+incr2[11])-1.732050807568877*incr1[11]-3.0*incr2[10]+3.0*incr1[10]-3.0*incr1[9]+3.0*incr1[8]-3.0*incr2[7]+3.0*(incr1[7]+incr2[6])-3.0*incr1[6]+3.0*incr1[5]-5.196152422706631*(incr1[4]+incr2[3])+5.196152422706631*incr1[3]-5.196152422706631*incr1[2]+5.196152422706631*incr1[1]-9.0*incr1[0]); 
  outlPos[7] = -0.02777777777777778*(incr2[15]-1.0*incr1[15]+1.732050807568877*incr2[14]-1.732050807568877*incr1[14]+1.732050807568877*incr2[13]-1.732050807568877*incr1[13]+1.732050807568877*(incr1[12]+incr2[11])-1.732050807568877*incr1[11]+3.0*incr2[10]-3.0*incr1[10]+3.0*(incr1[9]+incr1[8]+incr2[7])-3.0*incr1[7]+3.0*incr2[6]-3.0*incr1[6]+3.0*incr1[5]+5.196152422706631*(incr1[4]+incr2[3])-5.196152422706631*incr1[3]+5.196152422706631*(incr1[2]+incr1[1])+9.0*incr1[0]); 
  outrPos[0] = 0.02777777777777778*(incr2[15]+incr1[15]-1.732050807568877*(incr2[14]+incr1[14]+incr2[13]+incr1[13]+incr1[12]+incr2[11]+incr1[11])+3.0*(incr2[10]+incr1[10]+incr1[9]+incr1[8]+incr2[7]+incr1[7]+incr2[6]+incr1[6]+incr1[5])-5.196152422706631*(incr1[4]+incr2[3]+incr1[3]+incr1[2]+incr1[1])+9.0*incr1[0]); 
  outrPos[1] = -0.02777777777777778*(incr2[15]+incr1[15]+1.732050807568877*(incr2[14]+incr1[14])-1.732050807568877*(incr2[13]+incr1[13]+incr1[12]+incr2[11]+incr1[11])-3.0*(incr2[10]+incr1[10]+incr1[9])+3.0*incr1[8]-3.0*(incr2[7]+incr1[7])+3.0*(incr2[6]+incr1[6]+incr1[5])+5.196152422706631*(incr1[4]+incr2[3]+incr1[3]+incr1[2])-5.196152422706631*incr1[1]-9.0*incr1[0]); 
  outrPos[2] = -0.02777777777777778*(incr2[15]+incr1[15]-1.732050807568877*(incr2[14]+incr1[14])+1.732050807568877*(incr2[13]+incr1[13])-1.732050807568877*(incr1[12]+incr2[11]+incr1[11])-3.0*(incr2[10]+incr1[10])+3.0*incr1[9]-3.0*incr1[8]+3.0*(incr2[7]+incr1[7])-3.0*(incr2[6]+incr1[6])+3.0*incr1[5]+5.196152422706631*(incr1[4]+incr2[3]+incr1[3])-5.196152422706631*incr1[2]+5.196152422706631*incr1[1]-9.0*incr1[0]); 
  outrPos[3] = 0.02777777777777778*(incr2[15]+incr1[15]+1.732050807568877*(incr2[14]+incr1[14]+incr2[13]+incr1[13])-1.732050807568877*(incr1[12]+incr2[11]+incr1[11])+3.0*(incr2[10]+incr1[10])-3.0*(incr1[9]+incr1[8]+incr2[7]+incr1[7]+incr2[6]+incr1[6])+3.0*incr1[5]-5.196152422706631*(incr1[4]+incr2[3]+incr1[3])+5.196152422706631*(incr1[2]+incr1[1])+9.0*incr1[0]); 
  outrPos[4] = -0.02777777777777778*(incr2[15]+incr1[15]-1.732050807568877*(incr2[14]+incr1[14]+incr2[13]+incr1[13]+incr1[12])+1.732050807568877*(incr2[11]+incr1[11])+3.0*(incr2[10]+incr1[10]+incr1[9]+incr1[8])-3.0*(incr2[7]+incr1[7]+incr2[6]+incr1[6]+incr1[5])-5.196152422706631*incr1[4]+5.196152422706631*(incr2[3]+incr1[3]+incr1[2]+incr1[1])-9.0*incr1[0]); 
  outrPos[5] = 0.02777777777777778*(incr2[15]+incr1[15]+1.732050807568877*(incr2[14]+incr1[14])-1.732050807568877*(incr2[13]+incr1[13]+incr1[12])+1.732050807568877*(incr2[11]+incr1[11])-3.0*(incr2[10]+incr1[10]+incr1[9])+3.0*(incr1[8]+incr2[7]+incr1[7])-3.0*(incr2[6]+incr1[6]+incr1[5])+5.196152422706631*incr1[4]-5.196152422706631*(incr2[3]+incr1[3]+incr1[2])+5.196152422706631*incr1[1]+9.0*incr1[0]); 
  outrPos[6] = 0.02777777777777778*(incr2[15]+incr1[15]-1.732050807568877*(incr2[14]+incr1[14])+1.732050807568877*(incr2[13]+incr1[13])-1.732050807568877*incr1[12]+1.732050807568877*(incr2[11]+incr1[11])-3.0*(incr2[10]+incr1[10])+3.0*incr1[9]-3.0*(incr1[8]+incr2[7]+incr1[7])+3.0*(incr2[6]+incr1[6])-3.0*incr1[5]+5.196152422706631*incr1[4]-5.196152422706631*(incr2[3]+incr1[3])+5.196152422706631*incr1[2]-5.196152422706631*incr1[1]+9.0*incr1[0]); 
  outrPos[7] = -0.02777777777777778*(incr2[15]+incr1[15]+1.732050807568877*(incr2[14]+incr1[14]+incr2[13]+incr1[13])-1.732050807568877*incr1[12]+1.732050807568877*(incr2[11]+incr1[11])+3.0*(incr2[10]+incr1[10])-3.0*(incr1[9]+incr1[8])+3.0*(incr2[7]+incr1[7]+incr2[6]+incr1[6])-3.0*incr1[5]-5.196152422706631*incr1[4]+5.196152422706631*(incr2[3]+incr1[3])-5.196152422706631*(incr1[2]+incr1[1])-9.0*incr1[0]); 
  if(outlPos[0] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]-1.732050807568877*(fl[14]+fl[13])+1.732050807568877*fl[12]-1.732050807568877*fl[11]+3.0*fl[10]-3.0*(fl[9]+fl[8])+3.0*(fl[7]+fl[6])-3.0*fl[5]+5.196152422706631*fl[4]-5.196152422706631*fl[3]+5.196152422706631*(fl[2]+fl[1])-9.0*fl[0]))/dtApprox/outlPos[0]); 
  else limFac = 1.0; 
  if(outrPos[0] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]-1.732050807568877*(fr[14]+fr[13]+fr[12]+fr[11])+3.0*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-5.196152422706631*(fr[4]+fr[3]+fr[2]+fr[1])+9.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[0] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[1] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*fl[14]-1.732050807568877*fl[13]+1.732050807568877*fl[12]-1.732050807568877*fl[11]-3.0*fl[10]+3.0*fl[9]-3.0*(fl[8]+fl[7])+3.0*fl[6]-3.0*fl[5]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*fl[2]+5.196152422706631*fl[1]+9.0*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[1] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*fr[14]-1.732050807568877*(fr[13]+fr[12]+fr[11])-3.0*(fr[10]+fr[9])+3.0*fr[8]-3.0*fr[7]+3.0*(fr[6]+fr[5])+5.196152422706631*(fr[4]+fr[3]+fr[2])-5.196152422706631*fr[1]-9.0*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outrPos[1] *= limFac; 
  if(outlPos[2] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]-1.732050807568877*fl[14]+1.732050807568877*(fl[13]+fl[12])-1.732050807568877*fl[11]-3.0*(fl[10]+fl[9])+3.0*(fl[8]+fl[7])-3.0*(fl[6]+fl[5])-5.196152422706631*fl[4]+5.196152422706631*(fl[3]+fl[2])-5.196152422706631*fl[1]+9.0*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 1.0; 
  if(outrPos[2] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]-1.732050807568877*fr[14]+1.732050807568877*fr[13]-1.732050807568877*(fr[12]+fr[11])-3.0*fr[10]+3.0*fr[9]-3.0*fr[8]+3.0*fr[7]-3.0*fr[6]+3.0*fr[5]+5.196152422706631*(fr[4]+fr[3])-5.196152422706631*fr[2]+5.196152422706631*fr[1]-9.0*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outrPos[2] *= limFac; 
  if(outlPos[3] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12])-1.732050807568877*fl[11]+3.0*(fl[10]+fl[9]+fl[8])-3.0*(fl[7]+fl[6]+fl[5])+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[2]+fl[1])-9.0*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[3] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]+fr[13])-1.732050807568877*(fr[12]+fr[11])+3.0*fr[10]-3.0*(fr[9]+fr[8]+fr[7]+fr[6])+3.0*fr[5]-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*(fr[2]+fr[1])+9.0*fr[0]))/dtApprox/outrPos[3]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outrPos[3] *= limFac; 
  if(outlPos[4] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]-1.732050807568877*(fl[14]+fl[13])+1.732050807568877*(fl[12]+fl[11])+3.0*fl[10]-3.0*(fl[9]+fl[8]+fl[7]+fl[6])+3.0*fl[5]+5.196152422706631*(fl[4]+fl[3])-5.196152422706631*(fl[2]+fl[1])+9.0*fl[0]))/dtApprox/outlPos[4]); 
  else limFac = 1.0; 
  if(outrPos[4] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]-1.732050807568877*(fr[14]+fr[13]+fr[12])+1.732050807568877*fr[11]+3.0*(fr[10]+fr[9]+fr[8])-3.0*(fr[7]+fr[6]+fr[5])-5.196152422706631*fr[4]+5.196152422706631*(fr[3]+fr[2]+fr[1])-9.0*fr[0]))/dtApprox/outrPos[4]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[4] *= limFac; 
  outrPos[4] *= limFac; 
  if(outlPos[5] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*fl[14]-1.732050807568877*fl[13]+1.732050807568877*(fl[12]+fl[11])-3.0*fl[10]+3.0*fl[9]-3.0*fl[8]+3.0*fl[7]-3.0*fl[6]+3.0*fl[5]-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[2]-5.196152422706631*fl[1]-9.0*fl[0]))/dtApprox/outlPos[5]); 
  else limFac = 1.0; 
  if(outrPos[5] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*fr[14]-1.732050807568877*(fr[13]+fr[12])+1.732050807568877*fr[11]-3.0*(fr[10]+fr[9])+3.0*(fr[8]+fr[7])-3.0*(fr[6]+fr[5])+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[2])+5.196152422706631*fr[1]+9.0*fr[0]))/dtApprox/outrPos[5]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[5] *= limFac; 
  outrPos[5] *= limFac; 
  if(outlPos[6] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]-1.732050807568877*fl[14]+1.732050807568877*(fl[13]+fl[12]+fl[11])-3.0*(fl[10]+fl[9])+3.0*fl[8]-3.0*fl[7]+3.0*(fl[6]+fl[5])-5.196152422706631*(fl[4]+fl[3]+fl[2])+5.196152422706631*fl[1]-9.0*fl[0]))/dtApprox/outlPos[6]); 
  else limFac = 1.0; 
  if(outrPos[6] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]-1.732050807568877*fr[14]+1.732050807568877*fr[13]-1.732050807568877*fr[12]+1.732050807568877*fr[11]-3.0*fr[10]+3.0*fr[9]-3.0*(fr[8]+fr[7])+3.0*fr[6]-3.0*fr[5]+5.196152422706631*fr[4]-5.196152422706631*fr[3]+5.196152422706631*fr[2]-5.196152422706631*fr[1]+9.0*fr[0]))/dtApprox/outrPos[6]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[6] *= limFac; 
  outrPos[6] *= limFac; 
  if(outlPos[7] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]+fl[11])+3.0*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+5.196152422706631*(fl[4]+fl[3]+fl[2]+fl[1])+9.0*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[7] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]+fr[13])-1.732050807568877*fr[12]+1.732050807568877*fr[11]+3.0*fr[10]-3.0*(fr[9]+fr[8])+3.0*(fr[7]+fr[6])-3.0*fr[5]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*(fr[2]+fr[1])-9.0*fr[0]))/dtApprox/outrPos[7]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outrPos[7] *= limFac; 
  outr[0] += 0.25*outrPos[7]+0.25*outrPos[6]+0.25*outrPos[5]+0.25*outrPos[4]+0.25*outrPos[3]+0.25*outrPos[2]+0.25*outrPos[1]+0.25*outrPos[0]; 
  outr[1] += 0.4330127018922193*outrPos[7]-0.4330127018922193*outrPos[6]+0.4330127018922193*outrPos[5]-0.4330127018922193*outrPos[4]+0.4330127018922193*outrPos[3]-0.4330127018922193*outrPos[2]+0.4330127018922193*outrPos[1]-0.4330127018922193*outrPos[0]; 
  outr[2] += 0.4330127018922193*outrPos[7]+0.4330127018922193*outrPos[6]-0.4330127018922193*outrPos[5]-0.4330127018922193*outrPos[4]+0.4330127018922193*outrPos[3]+0.4330127018922193*outrPos[2]-0.4330127018922193*outrPos[1]-0.4330127018922193*outrPos[0]; 
  outr[3] += (-0.4330127018922193*outrPos[7])-0.4330127018922193*outrPos[6]-0.4330127018922193*outrPos[5]-0.4330127018922193*outrPos[4]-0.4330127018922193*outrPos[3]-0.4330127018922193*outrPos[2]-0.4330127018922193*outrPos[1]-0.4330127018922193*outrPos[0]; 
  outr[4] += 0.4330127018922193*outrPos[7]+0.4330127018922193*outrPos[6]+0.4330127018922193*outrPos[5]+0.4330127018922193*outrPos[4]-0.4330127018922193*outrPos[3]-0.4330127018922193*outrPos[2]-0.4330127018922193*outrPos[1]-0.4330127018922193*outrPos[0]; 
  outr[5] += 0.75*outrPos[7]-0.75*outrPos[6]-0.75*outrPos[5]+0.75*outrPos[4]+0.75*outrPos[3]-0.75*outrPos[2]-0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[6] += (-0.75*outrPos[7])+0.75*outrPos[6]-0.75*outrPos[5]+0.75*outrPos[4]-0.75*outrPos[3]+0.75*outrPos[2]-0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[7] += (-0.75*outrPos[7])-0.75*outrPos[6]+0.75*outrPos[5]+0.75*outrPos[4]-0.75*outrPos[3]-0.75*outrPos[2]+0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[8] += 0.75*outrPos[7]-0.75*outrPos[6]+0.75*outrPos[5]-0.75*outrPos[4]-0.75*outrPos[3]+0.75*outrPos[2]-0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[9] += 0.75*outrPos[7]+0.75*outrPos[6]-0.75*outrPos[5]-0.75*outrPos[4]-0.75*outrPos[3]-0.75*outrPos[2]+0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[10] += (-0.75*outrPos[7])-0.75*outrPos[6]-0.75*outrPos[5]-0.75*outrPos[4]+0.75*outrPos[3]+0.75*outrPos[2]+0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[11] += (-1.299038105676658*outrPos[7])+1.299038105676658*outrPos[6]+1.299038105676658*outrPos[5]-1.299038105676658*outrPos[4]-1.299038105676658*outrPos[3]+1.299038105676658*outrPos[2]+1.299038105676658*outrPos[1]-1.299038105676658*outrPos[0]; 
  outr[12] += 1.299038105676658*outrPos[7]-1.299038105676658*outrPos[6]-1.299038105676658*outrPos[5]+1.299038105676658*outrPos[4]-1.299038105676658*outrPos[3]+1.299038105676658*outrPos[2]+1.299038105676658*outrPos[1]-1.299038105676658*outrPos[0]; 
  outr[13] += (-1.299038105676658*outrPos[7])+1.299038105676658*outrPos[6]-1.299038105676658*outrPos[5]+1.299038105676658*outrPos[4]+1.299038105676658*outrPos[3]-1.299038105676658*outrPos[2]+1.299038105676658*outrPos[1]-1.299038105676658*outrPos[0]; 
  outr[14] += (-1.299038105676658*outrPos[7])-1.299038105676658*outrPos[6]+1.299038105676658*outrPos[5]+1.299038105676658*outrPos[4]+1.299038105676658*outrPos[3]+1.299038105676658*outrPos[2]-1.299038105676658*outrPos[1]-1.299038105676658*outrPos[0]; 
  outr[15] += (-2.25*outrPos[7])+2.25*outrPos[6]+2.25*outrPos[5]-2.25*outrPos[4]+2.25*outrPos[3]-2.25*outrPos[2]-2.25*outrPos[1]+2.25*outrPos[0]; 

  outl[0] += 0.25*outlPos[7]+0.25*outlPos[6]+0.25*outlPos[5]+0.25*outlPos[4]+0.25*outlPos[3]+0.25*outlPos[2]+0.25*outlPos[1]+0.25*outlPos[0]; 
  outl[1] += 0.4330127018922193*outlPos[7]-0.4330127018922193*outlPos[6]+0.4330127018922193*outlPos[5]-0.4330127018922193*outlPos[4]+0.4330127018922193*outlPos[3]-0.4330127018922193*outlPos[2]+0.4330127018922193*outlPos[1]-0.4330127018922193*outlPos[0]; 
  outl[2] += 0.4330127018922193*outlPos[7]+0.4330127018922193*outlPos[6]-0.4330127018922193*outlPos[5]-0.4330127018922193*outlPos[4]+0.4330127018922193*outlPos[3]+0.4330127018922193*outlPos[2]-0.4330127018922193*outlPos[1]-0.4330127018922193*outlPos[0]; 
  outl[3] += 0.4330127018922193*outlPos[7]+0.4330127018922193*outlPos[6]+0.4330127018922193*outlPos[5]+0.4330127018922193*outlPos[4]+0.4330127018922193*outlPos[3]+0.4330127018922193*outlPos[2]+0.4330127018922193*outlPos[1]+0.4330127018922193*outlPos[0]; 
  outl[4] += 0.4330127018922193*outlPos[7]+0.4330127018922193*outlPos[6]+0.4330127018922193*outlPos[5]+0.4330127018922193*outlPos[4]-0.4330127018922193*outlPos[3]-0.4330127018922193*outlPos[2]-0.4330127018922193*outlPos[1]-0.4330127018922193*outlPos[0]; 
  outl[5] += 0.75*outlPos[7]-0.75*outlPos[6]-0.75*outlPos[5]+0.75*outlPos[4]+0.75*outlPos[3]-0.75*outlPos[2]-0.75*outlPos[1]+0.75*outlPos[0]; 
  outl[6] += 0.75*outlPos[7]-0.75*outlPos[6]+0.75*outlPos[5]-0.75*outlPos[4]+0.75*outlPos[3]-0.75*outlPos[2]+0.75*outlPos[1]-0.75*outlPos[0]; 
  outl[7] += 0.75*outlPos[7]+0.75*outlPos[6]-0.75*outlPos[5]-0.75*outlPos[4]+0.75*outlPos[3]+0.75*outlPos[2]-0.75*outlPos[1]-0.75*outlPos[0]; 
  outl[8] += 0.75*outlPos[7]-0.75*outlPos[6]+0.75*outlPos[5]-0.75*outlPos[4]-0.75*outlPos[3]+0.75*outlPos[2]-0.75*outlPos[1]+0.75*outlPos[0]; 
  outl[9] += 0.75*outlPos[7]+0.75*outlPos[6]-0.75*outlPos[5]-0.75*outlPos[4]-0.75*outlPos[3]-0.75*outlPos[2]+0.75*outlPos[1]+0.75*outlPos[0]; 
  outl[10] += 0.75*outlPos[7]+0.75*outlPos[6]+0.75*outlPos[5]+0.75*outlPos[4]-0.75*outlPos[3]-0.75*outlPos[2]-0.75*outlPos[1]-0.75*outlPos[0]; 
  outl[11] += 1.299038105676658*outlPos[7]-1.299038105676658*outlPos[6]-1.299038105676658*outlPos[5]+1.299038105676658*outlPos[4]+1.299038105676658*outlPos[3]-1.299038105676658*outlPos[2]-1.299038105676658*outlPos[1]+1.299038105676658*outlPos[0]; 
  outl[12] += 1.299038105676658*outlPos[7]-1.299038105676658*outlPos[6]-1.299038105676658*outlPos[5]+1.299038105676658*outlPos[4]-1.299038105676658*outlPos[3]+1.299038105676658*outlPos[2]+1.299038105676658*outlPos[1]-1.299038105676658*outlPos[0]; 
  outl[13] += 1.299038105676658*outlPos[7]-1.299038105676658*outlPos[6]+1.299038105676658*outlPos[5]-1.299038105676658*outlPos[4]-1.299038105676658*outlPos[3]+1.299038105676658*outlPos[2]-1.299038105676658*outlPos[1]+1.299038105676658*outlPos[0]; 
  outl[14] += 1.299038105676658*outlPos[7]+1.299038105676658*outlPos[6]-1.299038105676658*outlPos[5]-1.299038105676658*outlPos[4]-1.299038105676658*outlPos[3]-1.299038105676658*outlPos[2]+1.299038105676658*outlPos[1]+1.299038105676658*outlPos[0]; 
  outl[15] += 2.25*outlPos[7]-2.25*outlPos[6]-2.25*outlPos[5]+2.25*outlPos[4]-2.25*outlPos[3]+2.25*outlPos[2]+2.25*outlPos[1]-2.25*outlPos[0]; 

  return std::abs(wl[2]-(0.5*nuUSum[0])/nuSum); 
} 
double GkLBOconstNuSurf2x2vSer_Mu_P1(const double m_, const double *cflRateByDirL, const double *cflRateByDirR, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double dtApprox, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // m_:              species mass. 
  // cflRateByDir[4]: CFL rate in each direction. 
  // w[4]:            Cell-center coordinates. 
  // dxv[4]:          Cell spacing. 
  // nuSum:           collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:       maximum midpoint value of v-u. 
  // nuUSum[4]:       sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:    sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:           Distribution function in left/right cells 
  // outl/outr:       Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[3]; 
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

  incr2[4] = (diffFac[3]*(0.25*fl[12]-0.25*fr[12])+diffFac[2]*(0.25*fl[9]-0.25*fr[9])+diffFac[1]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[3]*(fr[5]+fl[5])+diffFac[0]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*(diffFac[2]*(fr[2]+fl[2])+diffFac[1]*(fr[1]+fl[1])+diffFac[0]*(fr[0]+fl[0])))*rdvSq4R; 
  incr2[8] = (diffFac[2]*(0.25*fl[12]-0.25*fr[12])+diffFac[3]*(0.25*fl[9]-0.25*fr[9])+diffFac[0]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[2]*(fr[5]+fl[5])+diffFac[1]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*((fr[2]+fl[2])*diffFac[3]+diffFac[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*diffFac[1]))*rdvSq4R; 
  incr2[9] = (diffFac[1]*(0.25*fl[12]-0.25*fr[12])+diffFac[0]*(0.25*fl[9]-0.25*fr[9])+diffFac[3]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[1]*(fr[5]+fl[5])+diffFac[2]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*((fr[1]+fl[1])*diffFac[3]+diffFac[0]*(fr[2]+fl[2])+(fr[0]+fl[0])*diffFac[2]))*rdvSq4R; 
  incr2[10] = (diffFac[3]*(0.25*fl[15]-0.25*fr[15])+diffFac[2]*(0.25*fl[14]-0.25*fr[14])+diffFac[1]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[3]*(fr[11]+fl[11])+diffFac[0]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[2]*(fr[7]+fl[7])+diffFac[1]*(fr[6]+fl[6])+diffFac[0]*(fr[3]+fl[3])))*rdvSq4R; 
  incr2[12] = (diffFac[0]*(0.25*fl[12]-0.25*fr[12])+diffFac[1]*(0.25*fl[9]-0.25*fr[9])+diffFac[2]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[0]*(fr[5]+fl[5])+diffFac[3]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*((fr[0]+fl[0])*diffFac[3]+diffFac[1]*(fr[2]+fl[2])+(fr[1]+fl[1])*diffFac[2]))*rdvSq4R; 
  incr2[13] = (diffFac[2]*(0.25*fl[15]-0.25*fr[15])+diffFac[3]*(0.25*fl[14]-0.25*fr[14])+diffFac[0]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[2]*(fr[11]+fl[11])+diffFac[1]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[3]*(fr[7]+fl[7])+diffFac[0]*(fr[6]+fl[6])+diffFac[1]*(fr[3]+fl[3])))*rdvSq4R; 
  incr2[14] = (diffFac[1]*(0.25*fl[15]-0.25*fr[15])+diffFac[0]*(0.25*fl[14]-0.25*fr[14])+diffFac[3]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[1]*(fr[11]+fl[11])+diffFac[2]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[0]*(fr[7]+fl[7])+diffFac[3]*(fr[6]+fl[6])+diffFac[2]*(fr[3]+fl[3])))*rdvSq4R; 
  incr2[15] = (diffFac[0]*(0.25*fl[15]-0.25*fr[15])+diffFac[1]*(0.25*fl[14]-0.25*fr[14])+diffFac[2]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[0]*(fr[11]+fl[11])+diffFac[3]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[1]*(fr[7]+fl[7])+diffFac[2]*(fr[6]+fl[6])+diffFac[3]*(fr[3]+fl[3])))*rdvSq4R; 

  Gdiff[0] = (-1.082531754730548*(diffFac[3]*(fr[12]+fl[12])+diffFac[2]*(fr[9]+fl[9])+diffFac[1]*(fr[8]+fl[8])))+diffFac[3]*(1.125*fr[5]-1.125*fl[5])-1.082531754730548*diffFac[0]*(fr[4]+fl[4])+diffFac[2]*(1.125*fr[2]-1.125*fl[2])+diffFac[1]*(1.125*fr[1]-1.125*fl[1])+diffFac[0]*(1.125*fr[0]-1.125*fl[0]); 
  Gdiff[1] = (-1.082531754730548*(diffFac[2]*(fr[12]+fl[12])+diffFac[3]*(fr[9]+fl[9])+diffFac[0]*(fr[8]+fl[8])))+diffFac[2]*(1.125*fr[5]-1.125*fl[5])-1.082531754730548*diffFac[1]*(fr[4]+fl[4])+(1.125*fr[2]-1.125*fl[2])*diffFac[3]+diffFac[0]*(1.125*fr[1]-1.125*fl[1])+(1.125*fr[0]-1.125*fl[0])*diffFac[1]; 
  Gdiff[2] = (-1.082531754730548*(diffFac[1]*(fr[12]+fl[12])+diffFac[0]*(fr[9]+fl[9])+diffFac[3]*(fr[8]+fl[8])))+diffFac[1]*(1.125*fr[5]-1.125*fl[5])-1.082531754730548*diffFac[2]*(fr[4]+fl[4])+(1.125*fr[1]-1.125*fl[1])*diffFac[3]+diffFac[0]*(1.125*fr[2]-1.125*fl[2])+(1.125*fr[0]-1.125*fl[0])*diffFac[2]; 
  Gdiff[3] = (-1.082531754730548*(diffFac[3]*(fr[15]+fl[15])+diffFac[2]*(fr[14]+fl[14])+diffFac[1]*(fr[13]+fl[13])))+diffFac[3]*(1.125*fr[11]-1.125*fl[11])-1.082531754730548*diffFac[0]*(fr[10]+fl[10])+diffFac[2]*(1.125*fr[7]-1.125*fl[7])+diffFac[1]*(1.125*fr[6]-1.125*fl[6])+diffFac[0]*(1.125*fr[3]-1.125*fl[3]); 
  Gdiff[5] = (-1.082531754730548*(diffFac[0]*(fr[12]+fl[12])+diffFac[1]*(fr[9]+fl[9])+diffFac[2]*(fr[8]+fl[8])))+diffFac[0]*(1.125*fr[5]-1.125*fl[5])+diffFac[3]*((-1.082531754730548*(fr[4]+fl[4]))+1.125*fr[0]-1.125*fl[0])+diffFac[1]*(1.125*fr[2]-1.125*fl[2])+(1.125*fr[1]-1.125*fl[1])*diffFac[2]; 
  Gdiff[6] = (-1.082531754730548*(diffFac[2]*(fr[15]+fl[15])+diffFac[3]*(fr[14]+fl[14])+diffFac[0]*(fr[13]+fl[13])))+diffFac[2]*(1.125*fr[11]-1.125*fl[11])-1.082531754730548*diffFac[1]*(fr[10]+fl[10])+diffFac[3]*(1.125*fr[7]-1.125*fl[7])+diffFac[0]*(1.125*fr[6]-1.125*fl[6])+diffFac[1]*(1.125*fr[3]-1.125*fl[3]); 
  Gdiff[7] = (-1.082531754730548*(diffFac[1]*(fr[15]+fl[15])+diffFac[0]*(fr[14]+fl[14])+diffFac[3]*(fr[13]+fl[13])))+diffFac[1]*(1.125*fr[11]-1.125*fl[11])-1.082531754730548*diffFac[2]*(fr[10]+fl[10])+diffFac[0]*(1.125*fr[7]-1.125*fl[7])+diffFac[3]*(1.125*fr[6]-1.125*fl[6])+diffFac[2]*(1.125*fr[3]-1.125*fl[3]); 
  Gdiff[11] = (-1.082531754730548*(diffFac[0]*(fr[15]+fl[15])+diffFac[1]*(fr[14]+fl[14])+diffFac[2]*(fr[13]+fl[13])))+diffFac[0]*(1.125*fr[11]-1.125*fl[11])-1.082531754730548*diffFac[3]*(fr[10]+fl[10])+diffFac[1]*(1.125*fr[7]-1.125*fl[7])+diffFac[2]*(1.125*fr[6]-1.125*fl[6])+diffFac[3]*(1.125*fr[3]-1.125*fl[3]); 

  Ghat[0] = Gdiff[0]*rdv+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv+0.5*alphaDrSurf[0]*fUpwind[1]; 
  Ghat[2] = Gdiff[2]*rdv+0.5*alphaDrSurf[0]*fUpwind[2]; 
  Ghat[3] = Gdiff[3]*rdv+0.5*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[5] = Gdiff[5]*rdv+0.5*alphaDrSurf[0]*fUpwind[4]; 
  Ghat[6] = Gdiff[6]*rdv+0.5*alphaDrSurf[0]*fUpwind[5]; 
  Ghat[7] = Gdiff[7]*rdv+0.5*alphaDrSurf[0]*fUpwind[6]; 
  Ghat[11] = Gdiff[11]*rdv+0.5*alphaDrSurf[0]*fUpwind[7]; 

  double incr1[16]; 
  incr1[0] = -0.5*Ghat[0]*rdv2R; 
  incr1[1] = -0.5*Ghat[1]*rdv2R; 
  incr1[2] = -0.5*Ghat[2]*rdv2R; 
  incr1[3] = -0.5*Ghat[3]*rdv2R; 
  incr1[4] = 0.8660254037844386*Ghat[0]*rdv2R; 
  incr1[5] = -0.5*Ghat[5]*rdv2R; 
  incr1[6] = -0.5*Ghat[6]*rdv2R; 
  incr1[7] = -0.5*Ghat[7]*rdv2R; 
  incr1[8] = 0.8660254037844386*Ghat[1]*rdv2R; 
  incr1[9] = 0.8660254037844386*Ghat[2]*rdv2R; 
  incr1[10] = 0.8660254037844386*Ghat[3]*rdv2R; 
  incr1[11] = -0.5*Ghat[11]*rdv2R; 
  incr1[12] = 0.8660254037844386*Ghat[5]*rdv2R; 
  incr1[13] = 0.8660254037844386*Ghat[6]*rdv2R; 
  incr1[14] = 0.8660254037844386*Ghat[7]*rdv2R; 
  incr1[15] = 0.8660254037844386*Ghat[11]*rdv2R; 

  double fluxFracL, fluxFracR, limFac=1.0;
  double outlPos[8], outrPos[8]; 
  fluxFracL = cflRateByDirL[0] == 0. ? 0.25 : cflRateByDirL[2]/cflRateByDirL[0]; 
  fluxFracR = cflRateByDirR[0] == 0. ? 0.25 : cflRateByDirR[2]/cflRateByDirR[0]; 
  outlPos[0] = 0.02777777777777778*(incr2[15]-1.0*incr1[15]-1.732050807568877*incr2[14]+1.732050807568877*incr1[14]-1.732050807568877*incr2[13]+1.732050807568877*incr1[13]-1.732050807568877*incr2[12]+1.732050807568877*(incr1[12]+incr1[11])+3.0*incr2[10]-3.0*incr1[10]+3.0*incr2[9]-3.0*incr1[9]+3.0*incr2[8]-3.0*(incr1[8]+incr1[7]+incr1[6]+incr1[5])-5.196152422706631*incr2[4]+5.196152422706631*(incr1[4]+incr1[3]+incr1[2]+incr1[1])-9.0*incr1[0]); 
  outlPos[1] = -0.02777777777777778*(incr2[15]-1.0*incr1[15]+1.732050807568877*incr2[14]-1.732050807568877*(incr1[14]+incr2[13])+1.732050807568877*incr1[13]-1.732050807568877*incr2[12]+1.732050807568877*(incr1[12]+incr1[11])-3.0*incr2[10]+3.0*incr1[10]-3.0*incr2[9]+3.0*(incr1[9]+incr2[8])-3.0*incr1[8]+3.0*incr1[7]-3.0*(incr1[6]+incr1[5])+5.196152422706631*incr2[4]-5.196152422706631*(incr1[4]+incr1[3]+incr1[2])+5.196152422706631*incr1[1]+9.0*incr1[0]); 
  outlPos[2] = -0.02777777777777778*(incr2[15]-1.0*incr1[15]-1.732050807568877*incr2[14]+1.732050807568877*(incr1[14]+incr2[13])-1.732050807568877*(incr1[13]+incr2[12])+1.732050807568877*(incr1[12]+incr1[11])-3.0*incr2[10]+3.0*(incr1[10]+incr2[9])-3.0*(incr1[9]+incr2[8])+3.0*incr1[8]-3.0*incr1[7]+3.0*incr1[6]-3.0*incr1[5]+5.196152422706631*incr2[4]-5.196152422706631*(incr1[4]+incr1[3])+5.196152422706631*incr1[2]-5.196152422706631*incr1[1]+9.0*incr1[0]); 
  outlPos[3] = 0.02777777777777778*(incr2[15]-1.0*incr1[15]+1.732050807568877*incr2[14]-1.732050807568877*incr1[14]+1.732050807568877*incr2[13]-1.732050807568877*(incr1[13]+incr2[12])+1.732050807568877*(incr1[12]+incr1[11])+3.0*incr2[10]-3.0*(incr1[10]+incr2[9])+3.0*incr1[9]-3.0*incr2[8]+3.0*(incr1[8]+incr1[7]+incr1[6])-3.0*incr1[5]-5.196152422706631*incr2[4]+5.196152422706631*(incr1[4]+incr1[3])-5.196152422706631*(incr1[2]+incr1[1])-9.0*incr1[0]); 
  outlPos[4] = -0.02777777777777778*(incr2[15]-1.0*incr1[15]-1.732050807568877*incr2[14]+1.732050807568877*incr1[14]-1.732050807568877*incr2[13]+1.732050807568877*(incr1[13]+incr2[12])-1.732050807568877*incr1[12]+1.732050807568877*incr1[11]+3.0*incr2[10]-3.0*(incr1[10]+incr2[9])+3.0*incr1[9]-3.0*incr2[8]+3.0*incr1[8]-3.0*(incr1[7]+incr1[6])+3.0*incr1[5]+5.196152422706631*incr2[4]-5.196152422706631*incr1[4]+5.196152422706631*incr1[3]-5.196152422706631*(incr1[2]+incr1[1])+9.0*incr1[0]); 
  outlPos[5] = 0.02777777777777778*(incr2[15]-1.0*incr1[15]+1.732050807568877*incr2[14]-1.732050807568877*(incr1[14]+incr2[13])+1.732050807568877*(incr1[13]+incr2[12])-1.732050807568877*incr1[12]+1.732050807568877*incr1[11]-3.0*incr2[10]+3.0*(incr1[10]+incr2[9])-3.0*(incr1[9]+incr2[8])+3.0*(incr1[8]+incr1[7])-3.0*incr1[6]+3.0*incr1[5]-5.196152422706631*incr2[4]+5.196152422706631*incr1[4]-5.196152422706631*incr1[3]+5.196152422706631*incr1[2]-5.196152422706631*incr1[1]-9.0*incr1[0]); 
  outlPos[6] = 0.02777777777777778*(incr2[15]-1.0*incr1[15]-1.732050807568877*incr2[14]+1.732050807568877*(incr1[14]+incr2[13])-1.732050807568877*incr1[13]+1.732050807568877*incr2[12]-1.732050807568877*incr1[12]+1.732050807568877*incr1[11]-3.0*incr2[10]+3.0*incr1[10]-3.0*incr2[9]+3.0*(incr1[9]+incr2[8])-3.0*(incr1[8]+incr1[7])+3.0*(incr1[6]+incr1[5])-5.196152422706631*incr2[4]+5.196152422706631*incr1[4]-5.196152422706631*(incr1[3]+incr1[2])+5.196152422706631*incr1[1]-9.0*incr1[0]); 
  outlPos[7] = -0.02777777777777778*(incr2[15]-1.0*incr1[15]+1.732050807568877*incr2[14]-1.732050807568877*incr1[14]+1.732050807568877*incr2[13]-1.732050807568877*incr1[13]+1.732050807568877*incr2[12]-1.732050807568877*incr1[12]+1.732050807568877*incr1[11]+3.0*incr2[10]-3.0*incr1[10]+3.0*incr2[9]-3.0*incr1[9]+3.0*incr2[8]-3.0*incr1[8]+3.0*(incr1[7]+incr1[6]+incr1[5])+5.196152422706631*incr2[4]-5.196152422706631*incr1[4]+5.196152422706631*(incr1[3]+incr1[2]+incr1[1])+9.0*incr1[0]); 
  outrPos[0] = 0.02777777777777778*(incr2[15]+incr1[15]-1.732050807568877*(incr2[14]+incr1[14]+incr2[13]+incr1[13]+incr2[12]+incr1[12]+incr1[11])+3.0*(incr2[10]+incr1[10]+incr2[9]+incr1[9]+incr2[8]+incr1[8]+incr1[7]+incr1[6]+incr1[5])-5.196152422706631*(incr2[4]+incr1[4]+incr1[3]+incr1[2]+incr1[1])+9.0*incr1[0]); 
  outrPos[1] = -0.02777777777777778*(incr2[15]+incr1[15]+1.732050807568877*(incr2[14]+incr1[14])-1.732050807568877*(incr2[13]+incr1[13]+incr2[12]+incr1[12]+incr1[11])-3.0*(incr2[10]+incr1[10]+incr2[9]+incr1[9])+3.0*(incr2[8]+incr1[8])-3.0*incr1[7]+3.0*(incr1[6]+incr1[5])+5.196152422706631*(incr2[4]+incr1[4]+incr1[3]+incr1[2])-5.196152422706631*incr1[1]-9.0*incr1[0]); 
  outrPos[2] = -0.02777777777777778*(incr2[15]+incr1[15]-1.732050807568877*(incr2[14]+incr1[14])+1.732050807568877*(incr2[13]+incr1[13])-1.732050807568877*(incr2[12]+incr1[12]+incr1[11])-3.0*(incr2[10]+incr1[10])+3.0*(incr2[9]+incr1[9])-3.0*(incr2[8]+incr1[8])+3.0*incr1[7]-3.0*incr1[6]+3.0*incr1[5]+5.196152422706631*(incr2[4]+incr1[4]+incr1[3])-5.196152422706631*incr1[2]+5.196152422706631*incr1[1]-9.0*incr1[0]); 
  outrPos[3] = 0.02777777777777778*(incr2[15]+incr1[15]+1.732050807568877*(incr2[14]+incr1[14]+incr2[13]+incr1[13])-1.732050807568877*(incr2[12]+incr1[12]+incr1[11])+3.0*(incr2[10]+incr1[10])-3.0*(incr2[9]+incr1[9]+incr2[8]+incr1[8]+incr1[7]+incr1[6])+3.0*incr1[5]-5.196152422706631*(incr2[4]+incr1[4]+incr1[3])+5.196152422706631*(incr1[2]+incr1[1])+9.0*incr1[0]); 
  outrPos[4] = -0.02777777777777778*(incr2[15]+incr1[15]-1.732050807568877*(incr2[14]+incr1[14]+incr2[13]+incr1[13])+1.732050807568877*(incr2[12]+incr1[12])-1.732050807568877*incr1[11]+3.0*(incr2[10]+incr1[10])-3.0*(incr2[9]+incr1[9]+incr2[8]+incr1[8])+3.0*(incr1[7]+incr1[6])-3.0*incr1[5]+5.196152422706631*(incr2[4]+incr1[4])-5.196152422706631*incr1[3]+5.196152422706631*(incr1[2]+incr1[1])-9.0*incr1[0]); 
  outrPos[5] = 0.02777777777777778*(incr2[15]+incr1[15]+1.732050807568877*(incr2[14]+incr1[14])-1.732050807568877*(incr2[13]+incr1[13])+1.732050807568877*(incr2[12]+incr1[12])-1.732050807568877*incr1[11]-3.0*(incr2[10]+incr1[10])+3.0*(incr2[9]+incr1[9])-3.0*(incr2[8]+incr1[8]+incr1[7])+3.0*incr1[6]-3.0*incr1[5]-5.196152422706631*(incr2[4]+incr1[4])+5.196152422706631*incr1[3]-5.196152422706631*incr1[2]+5.196152422706631*incr1[1]+9.0*incr1[0]); 
  outrPos[6] = 0.02777777777777778*(incr2[15]+incr1[15]-1.732050807568877*(incr2[14]+incr1[14])+1.732050807568877*(incr2[13]+incr1[13]+incr2[12]+incr1[12])-1.732050807568877*incr1[11]-3.0*(incr2[10]+incr1[10]+incr2[9]+incr1[9])+3.0*(incr2[8]+incr1[8]+incr1[7])-3.0*(incr1[6]+incr1[5])-5.196152422706631*(incr2[4]+incr1[4])+5.196152422706631*(incr1[3]+incr1[2])-5.196152422706631*incr1[1]+9.0*incr1[0]); 
  outrPos[7] = -0.02777777777777778*(incr2[15]+incr1[15]+1.732050807568877*(incr2[14]+incr1[14]+incr2[13]+incr1[13]+incr2[12]+incr1[12])-1.732050807568877*incr1[11]+3.0*(incr2[10]+incr1[10]+incr2[9]+incr1[9]+incr2[8]+incr1[8])-3.0*(incr1[7]+incr1[6]+incr1[5])+5.196152422706631*(incr2[4]+incr1[4])-5.196152422706631*(incr1[3]+incr1[2]+incr1[1])-9.0*incr1[0]); 
  if(outlPos[0] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]-1.732050807568877*(fl[14]+fl[13]+fl[12])+1.732050807568877*fl[11]+3.0*(fl[10]+fl[9]+fl[8])-3.0*(fl[7]+fl[6]+fl[5])-5.196152422706631*fl[4]+5.196152422706631*(fl[3]+fl[2]+fl[1])-9.0*fl[0]))/dtApprox/outlPos[0]); 
  else limFac = 1.0; 
  if(outrPos[0] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]-1.732050807568877*(fr[14]+fr[13]+fr[12]+fr[11])+3.0*(fr[10]+fr[9]+fr[8]+fr[7]+fr[6]+fr[5])-5.196152422706631*(fr[4]+fr[3]+fr[2]+fr[1])+9.0*fr[0]))/dtApprox/outrPos[0]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[0] *= limFac; 
  outrPos[0] *= limFac; 
  if(outlPos[1] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*fl[14]-1.732050807568877*(fl[13]+fl[12])+1.732050807568877*fl[11]-3.0*(fl[10]+fl[9])+3.0*(fl[8]+fl[7])-3.0*(fl[6]+fl[5])+5.196152422706631*fl[4]-5.196152422706631*(fl[3]+fl[2])+5.196152422706631*fl[1]+9.0*fl[0]))/dtApprox/outlPos[1]); 
  else limFac = 1.0; 
  if(outrPos[1] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*fr[14]-1.732050807568877*(fr[13]+fr[12]+fr[11])-3.0*(fr[10]+fr[9])+3.0*fr[8]-3.0*fr[7]+3.0*(fr[6]+fr[5])+5.196152422706631*(fr[4]+fr[3]+fr[2])-5.196152422706631*fr[1]-9.0*fr[0]))/dtApprox/outrPos[1]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[1] *= limFac; 
  outrPos[1] *= limFac; 
  if(outlPos[2] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]-1.732050807568877*fl[14]+1.732050807568877*fl[13]-1.732050807568877*fl[12]+1.732050807568877*fl[11]-3.0*fl[10]+3.0*fl[9]-3.0*(fl[8]+fl[7])+3.0*fl[6]-3.0*fl[5]+5.196152422706631*fl[4]-5.196152422706631*fl[3]+5.196152422706631*fl[2]-5.196152422706631*fl[1]+9.0*fl[0]))/dtApprox/outlPos[2]); 
  else limFac = 1.0; 
  if(outrPos[2] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]-1.732050807568877*fr[14]+1.732050807568877*fr[13]-1.732050807568877*(fr[12]+fr[11])-3.0*fr[10]+3.0*fr[9]-3.0*fr[8]+3.0*fr[7]-3.0*fr[6]+3.0*fr[5]+5.196152422706631*(fr[4]+fr[3])-5.196152422706631*fr[2]+5.196152422706631*fr[1]-9.0*fr[0]))/dtApprox/outrPos[2]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[2] *= limFac; 
  outrPos[2] *= limFac; 
  if(outlPos[3] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13])-1.732050807568877*fl[12]+1.732050807568877*fl[11]+3.0*fl[10]-3.0*(fl[9]+fl[8])+3.0*(fl[7]+fl[6])-3.0*fl[5]-5.196152422706631*fl[4]+5.196152422706631*fl[3]-5.196152422706631*(fl[2]+fl[1])-9.0*fl[0]))/dtApprox/outlPos[3]); 
  else limFac = 1.0; 
  if(outrPos[3] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]+fr[13])-1.732050807568877*(fr[12]+fr[11])+3.0*fr[10]-3.0*(fr[9]+fr[8]+fr[7]+fr[6])+3.0*fr[5]-5.196152422706631*(fr[4]+fr[3])+5.196152422706631*(fr[2]+fr[1])+9.0*fr[0]))/dtApprox/outrPos[3]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[3] *= limFac; 
  outrPos[3] *= limFac; 
  if(outlPos[4] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]-1.732050807568877*(fl[14]+fl[13])+1.732050807568877*(fl[12]+fl[11])+3.0*fl[10]-3.0*(fl[9]+fl[8]+fl[7]+fl[6])+3.0*fl[5]+5.196152422706631*(fl[4]+fl[3])-5.196152422706631*(fl[2]+fl[1])+9.0*fl[0]))/dtApprox/outlPos[4]); 
  else limFac = 1.0; 
  if(outrPos[4] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]-1.732050807568877*(fr[14]+fr[13])+1.732050807568877*fr[12]-1.732050807568877*fr[11]+3.0*fr[10]-3.0*(fr[9]+fr[8])+3.0*(fr[7]+fr[6])-3.0*fr[5]+5.196152422706631*fr[4]-5.196152422706631*fr[3]+5.196152422706631*(fr[2]+fr[1])-9.0*fr[0]))/dtApprox/outrPos[4]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[4] *= limFac; 
  outrPos[4] *= limFac; 
  if(outlPos[5] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]+1.732050807568877*fl[14]-1.732050807568877*fl[13]+1.732050807568877*(fl[12]+fl[11])-3.0*fl[10]+3.0*fl[9]-3.0*fl[8]+3.0*fl[7]-3.0*fl[6]+3.0*fl[5]-5.196152422706631*(fl[4]+fl[3])+5.196152422706631*fl[2]-5.196152422706631*fl[1]-9.0*fl[0]))/dtApprox/outlPos[5]); 
  else limFac = 1.0; 
  if(outrPos[5] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]+1.732050807568877*fr[14]-1.732050807568877*fr[13]+1.732050807568877*fr[12]-1.732050807568877*fr[11]-3.0*fr[10]+3.0*fr[9]-3.0*(fr[8]+fr[7])+3.0*fr[6]-3.0*fr[5]-5.196152422706631*fr[4]+5.196152422706631*fr[3]-5.196152422706631*fr[2]+5.196152422706631*fr[1]+9.0*fr[0]))/dtApprox/outrPos[5]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[5] *= limFac; 
  outrPos[5] *= limFac; 
  if(outlPos[6] < 0.) limFac = std::min(1.0, -fluxFracL*(-0.02777777777777778*(fl[15]-1.732050807568877*fl[14]+1.732050807568877*(fl[13]+fl[12]+fl[11])-3.0*(fl[10]+fl[9])+3.0*fl[8]-3.0*fl[7]+3.0*(fl[6]+fl[5])-5.196152422706631*(fl[4]+fl[3]+fl[2])+5.196152422706631*fl[1]-9.0*fl[0]))/dtApprox/outlPos[6]); 
  else limFac = 1.0; 
  if(outrPos[6] < 0.) limFac = std::min(limFac, -fluxFracR*(0.02777777777777778*(fr[15]-1.732050807568877*fr[14]+1.732050807568877*(fr[13]+fr[12])-1.732050807568877*fr[11]-3.0*(fr[10]+fr[9])+3.0*(fr[8]+fr[7])-3.0*(fr[6]+fr[5])-5.196152422706631*fr[4]+5.196152422706631*(fr[3]+fr[2])-5.196152422706631*fr[1]+9.0*fr[0]))/dtApprox/outrPos[6]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[6] *= limFac; 
  outrPos[6] *= limFac; 
  if(outlPos[7] < 0.) limFac = std::min(1.0, -fluxFracL*(0.02777777777777778*(fl[15]+1.732050807568877*(fl[14]+fl[13]+fl[12]+fl[11])+3.0*(fl[10]+fl[9]+fl[8]+fl[7]+fl[6]+fl[5])+5.196152422706631*(fl[4]+fl[3]+fl[2]+fl[1])+9.0*fl[0]))/dtApprox/outlPos[7]); 
  else limFac = 1.0; 
  if(outrPos[7] < 0.) limFac = std::min(limFac, -fluxFracR*(-0.02777777777777778*(fr[15]+1.732050807568877*(fr[14]+fr[13]+fr[12])-1.732050807568877*fr[11]+3.0*(fr[10]+fr[9]+fr[8])-3.0*(fr[7]+fr[6]+fr[5])+5.196152422706631*fr[4]-5.196152422706631*(fr[3]+fr[2]+fr[1])-9.0*fr[0]))/dtApprox/outrPos[7]); 
  if(limFac < 0.) limFac = 0.; 
  outlPos[7] *= limFac; 
  outrPos[7] *= limFac; 
  outr[0] += 0.25*outrPos[7]+0.25*outrPos[6]+0.25*outrPos[5]+0.25*outrPos[4]+0.25*outrPos[3]+0.25*outrPos[2]+0.25*outrPos[1]+0.25*outrPos[0]; 
  outr[1] += 0.4330127018922193*outrPos[7]-0.4330127018922193*outrPos[6]+0.4330127018922193*outrPos[5]-0.4330127018922193*outrPos[4]+0.4330127018922193*outrPos[3]-0.4330127018922193*outrPos[2]+0.4330127018922193*outrPos[1]-0.4330127018922193*outrPos[0]; 
  outr[2] += 0.4330127018922193*outrPos[7]+0.4330127018922193*outrPos[6]-0.4330127018922193*outrPos[5]-0.4330127018922193*outrPos[4]+0.4330127018922193*outrPos[3]+0.4330127018922193*outrPos[2]-0.4330127018922193*outrPos[1]-0.4330127018922193*outrPos[0]; 
  outr[3] += 0.4330127018922193*outrPos[7]+0.4330127018922193*outrPos[6]+0.4330127018922193*outrPos[5]+0.4330127018922193*outrPos[4]-0.4330127018922193*outrPos[3]-0.4330127018922193*outrPos[2]-0.4330127018922193*outrPos[1]-0.4330127018922193*outrPos[0]; 
  outr[4] += (-0.4330127018922193*outrPos[7])-0.4330127018922193*outrPos[6]-0.4330127018922193*outrPos[5]-0.4330127018922193*outrPos[4]-0.4330127018922193*outrPos[3]-0.4330127018922193*outrPos[2]-0.4330127018922193*outrPos[1]-0.4330127018922193*outrPos[0]; 
  outr[5] += 0.75*outrPos[7]-0.75*outrPos[6]-0.75*outrPos[5]+0.75*outrPos[4]+0.75*outrPos[3]-0.75*outrPos[2]-0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[6] += 0.75*outrPos[7]-0.75*outrPos[6]+0.75*outrPos[5]-0.75*outrPos[4]-0.75*outrPos[3]+0.75*outrPos[2]-0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[7] += 0.75*outrPos[7]+0.75*outrPos[6]-0.75*outrPos[5]-0.75*outrPos[4]-0.75*outrPos[3]-0.75*outrPos[2]+0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[8] += (-0.75*outrPos[7])+0.75*outrPos[6]-0.75*outrPos[5]+0.75*outrPos[4]-0.75*outrPos[3]+0.75*outrPos[2]-0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[9] += (-0.75*outrPos[7])-0.75*outrPos[6]+0.75*outrPos[5]+0.75*outrPos[4]-0.75*outrPos[3]-0.75*outrPos[2]+0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[10] += (-0.75*outrPos[7])-0.75*outrPos[6]-0.75*outrPos[5]-0.75*outrPos[4]+0.75*outrPos[3]+0.75*outrPos[2]+0.75*outrPos[1]+0.75*outrPos[0]; 
  outr[11] += 1.299038105676658*outrPos[7]-1.299038105676658*outrPos[6]-1.299038105676658*outrPos[5]+1.299038105676658*outrPos[4]-1.299038105676658*outrPos[3]+1.299038105676658*outrPos[2]+1.299038105676658*outrPos[1]-1.299038105676658*outrPos[0]; 
  outr[12] += (-1.299038105676658*outrPos[7])+1.299038105676658*outrPos[6]+1.299038105676658*outrPos[5]-1.299038105676658*outrPos[4]-1.299038105676658*outrPos[3]+1.299038105676658*outrPos[2]+1.299038105676658*outrPos[1]-1.299038105676658*outrPos[0]; 
  outr[13] += (-1.299038105676658*outrPos[7])+1.299038105676658*outrPos[6]-1.299038105676658*outrPos[5]+1.299038105676658*outrPos[4]+1.299038105676658*outrPos[3]-1.299038105676658*outrPos[2]+1.299038105676658*outrPos[1]-1.299038105676658*outrPos[0]; 
  outr[14] += (-1.299038105676658*outrPos[7])-1.299038105676658*outrPos[6]+1.299038105676658*outrPos[5]+1.299038105676658*outrPos[4]+1.299038105676658*outrPos[3]+1.299038105676658*outrPos[2]-1.299038105676658*outrPos[1]-1.299038105676658*outrPos[0]; 
  outr[15] += (-2.25*outrPos[7])+2.25*outrPos[6]+2.25*outrPos[5]-2.25*outrPos[4]+2.25*outrPos[3]-2.25*outrPos[2]-2.25*outrPos[1]+2.25*outrPos[0]; 

  outl[0] += 0.25*outlPos[7]+0.25*outlPos[6]+0.25*outlPos[5]+0.25*outlPos[4]+0.25*outlPos[3]+0.25*outlPos[2]+0.25*outlPos[1]+0.25*outlPos[0]; 
  outl[1] += 0.4330127018922193*outlPos[7]-0.4330127018922193*outlPos[6]+0.4330127018922193*outlPos[5]-0.4330127018922193*outlPos[4]+0.4330127018922193*outlPos[3]-0.4330127018922193*outlPos[2]+0.4330127018922193*outlPos[1]-0.4330127018922193*outlPos[0]; 
  outl[2] += 0.4330127018922193*outlPos[7]+0.4330127018922193*outlPos[6]-0.4330127018922193*outlPos[5]-0.4330127018922193*outlPos[4]+0.4330127018922193*outlPos[3]+0.4330127018922193*outlPos[2]-0.4330127018922193*outlPos[1]-0.4330127018922193*outlPos[0]; 
  outl[3] += 0.4330127018922193*outlPos[7]+0.4330127018922193*outlPos[6]+0.4330127018922193*outlPos[5]+0.4330127018922193*outlPos[4]-0.4330127018922193*outlPos[3]-0.4330127018922193*outlPos[2]-0.4330127018922193*outlPos[1]-0.4330127018922193*outlPos[0]; 
  outl[4] += 0.4330127018922193*outlPos[7]+0.4330127018922193*outlPos[6]+0.4330127018922193*outlPos[5]+0.4330127018922193*outlPos[4]+0.4330127018922193*outlPos[3]+0.4330127018922193*outlPos[2]+0.4330127018922193*outlPos[1]+0.4330127018922193*outlPos[0]; 
  outl[5] += 0.75*outlPos[7]-0.75*outlPos[6]-0.75*outlPos[5]+0.75*outlPos[4]+0.75*outlPos[3]-0.75*outlPos[2]-0.75*outlPos[1]+0.75*outlPos[0]; 
  outl[6] += 0.75*outlPos[7]-0.75*outlPos[6]+0.75*outlPos[5]-0.75*outlPos[4]-0.75*outlPos[3]+0.75*outlPos[2]-0.75*outlPos[1]+0.75*outlPos[0]; 
  outl[7] += 0.75*outlPos[7]+0.75*outlPos[6]-0.75*outlPos[5]-0.75*outlPos[4]-0.75*outlPos[3]-0.75*outlPos[2]+0.75*outlPos[1]+0.75*outlPos[0]; 
  outl[8] += 0.75*outlPos[7]-0.75*outlPos[6]+0.75*outlPos[5]-0.75*outlPos[4]+0.75*outlPos[3]-0.75*outlPos[2]+0.75*outlPos[1]-0.75*outlPos[0]; 
  outl[9] += 0.75*outlPos[7]+0.75*outlPos[6]-0.75*outlPos[5]-0.75*outlPos[4]+0.75*outlPos[3]+0.75*outlPos[2]-0.75*outlPos[1]-0.75*outlPos[0]; 
  outl[10] += 0.75*outlPos[7]+0.75*outlPos[6]+0.75*outlPos[5]+0.75*outlPos[4]-0.75*outlPos[3]-0.75*outlPos[2]-0.75*outlPos[1]-0.75*outlPos[0]; 
  outl[11] += 1.299038105676658*outlPos[7]-1.299038105676658*outlPos[6]-1.299038105676658*outlPos[5]+1.299038105676658*outlPos[4]-1.299038105676658*outlPos[3]+1.299038105676658*outlPos[2]+1.299038105676658*outlPos[1]-1.299038105676658*outlPos[0]; 
  outl[12] += 1.299038105676658*outlPos[7]-1.299038105676658*outlPos[6]-1.299038105676658*outlPos[5]+1.299038105676658*outlPos[4]+1.299038105676658*outlPos[3]-1.299038105676658*outlPos[2]-1.299038105676658*outlPos[1]+1.299038105676658*outlPos[0]; 
  outl[13] += 1.299038105676658*outlPos[7]-1.299038105676658*outlPos[6]+1.299038105676658*outlPos[5]-1.299038105676658*outlPos[4]-1.299038105676658*outlPos[3]+1.299038105676658*outlPos[2]-1.299038105676658*outlPos[1]+1.299038105676658*outlPos[0]; 
  outl[14] += 1.299038105676658*outlPos[7]+1.299038105676658*outlPos[6]-1.299038105676658*outlPos[5]-1.299038105676658*outlPos[4]-1.299038105676658*outlPos[3]-1.299038105676658*outlPos[2]+1.299038105676658*outlPos[1]+1.299038105676658*outlPos[0]; 
  outl[15] += 2.25*outlPos[7]-2.25*outlPos[6]-2.25*outlPos[5]+2.25*outlPos[4]-2.25*outlPos[3]+2.25*outlPos[2]+2.25*outlPos[1]-2.25*outlPos[0]; 

  return std::abs(2.0*wl[3]); 
} 
