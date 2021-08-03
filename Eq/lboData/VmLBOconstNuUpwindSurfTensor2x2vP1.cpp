#include <VmLBOModDecl.h> 
double VmLBOconstNuUpwindSurf2x2vTensor_VX_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[8]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[2]; 
  double rdv2R = 2.0/dxvr[2]; 
  double rdvSq4L = 4.0/(dxvl[2]*dxvl[2]); 
  double rdvSq4R = 4.0/(dxvr[2]*dxvr[2]); 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf[8]; 
  alphaDrSurf[0] = 0.7071067811865475*((4.0*wl[2]+2.0*dxvl[2])*nuSum-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = -1.414213562373095*sumNuUx[1]; 
  alphaDrSurf[2] = -1.414213562373095*sumNuUx[2]; 
  alphaDrSurf[4] = -1.414213562373095*sumNuUx[3]; 

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


  incr2[3] = -0.125*(2.0*nuVtSqSum[3]*fr[11]-2.0*nuVtSqSum[3]*fl[11]+2.0*nuVtSqSum[2]*fr[7]-2.0*nuVtSqSum[2]*fl[7]+2.0*nuVtSqSum[1]*fr[6]-2.0*nuVtSqSum[1]*fl[6]-1.732050807568877*nuVtSqSum[3]*fr[5]-1.732050807568877*nuVtSqSum[3]*fl[5]+2.0*nuVtSqSum[0]*fr[3]-2.0*nuVtSqSum[0]*fl[3]+((-1.732050807568877*fr[2])-1.732050807568877*fl[2])*nuVtSqSum[2]+((-1.732050807568877*fr[1])-1.732050807568877*fl[1])*nuVtSqSum[1]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[0]); 
  incr2[6] = -0.125*(2.0*nuVtSqSum[2]*fr[11]-2.0*nuVtSqSum[2]*fl[11]+2.0*nuVtSqSum[3]*fr[7]-2.0*nuVtSqSum[3]*fl[7]+2.0*nuVtSqSum[0]*fr[6]-2.0*nuVtSqSum[0]*fl[6]-1.732050807568877*nuVtSqSum[2]*fr[5]-1.732050807568877*nuVtSqSum[2]*fl[5]+((-1.732050807568877*fr[2])-1.732050807568877*fl[2])*nuVtSqSum[3]+2.0*nuVtSqSum[1]*fr[3]-2.0*nuVtSqSum[1]*fl[3]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0]*fr[1]-1.732050807568877*nuVtSqSum[0]*fl[1]); 
  incr2[7] = -0.125*(2.0*nuVtSqSum[1]*fr[11]-2.0*nuVtSqSum[1]*fl[11]+2.0*nuVtSqSum[0]*fr[7]-2.0*nuVtSqSum[0]*fl[7]+2.0*nuVtSqSum[3]*fr[6]-2.0*nuVtSqSum[3]*fl[6]-1.732050807568877*nuVtSqSum[1]*fr[5]-1.732050807568877*nuVtSqSum[1]*fl[5]+((-1.732050807568877*fr[1])-1.732050807568877*fl[1])*nuVtSqSum[3]+2.0*nuVtSqSum[2]*fr[3]-2.0*nuVtSqSum[2]*fl[3]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[2]-1.732050807568877*nuVtSqSum[0]*fr[2]-1.732050807568877*nuVtSqSum[0]*fl[2]); 
  incr2[10] = -0.125*(2.0*nuVtSqSum[3]*fr[15]-2.0*nuVtSqSum[3]*fl[15]+2.0*nuVtSqSum[2]*fr[14]-2.0*nuVtSqSum[2]*fl[14]+2.0*nuVtSqSum[1]*fr[13]-2.0*nuVtSqSum[1]*fl[13]-1.732050807568877*nuVtSqSum[3]*fr[12]-1.732050807568877*nuVtSqSum[3]*fl[12]+2.0*nuVtSqSum[0]*fr[10]-2.0*nuVtSqSum[0]*fl[10]-1.732050807568877*nuVtSqSum[2]*fr[9]-1.732050807568877*nuVtSqSum[2]*fl[9]-1.732050807568877*nuVtSqSum[1]*fr[8]-1.732050807568877*nuVtSqSum[1]*fl[8]-1.732050807568877*nuVtSqSum[0]*fr[4]-1.732050807568877*nuVtSqSum[0]*fl[4]); 
  incr2[11] = -0.125*(2.0*nuVtSqSum[0]*fr[11]-2.0*nuVtSqSum[0]*fl[11]+2.0*nuVtSqSum[1]*fr[7]-2.0*nuVtSqSum[1]*fl[7]+2.0*nuVtSqSum[2]*fr[6]-2.0*nuVtSqSum[2]*fl[6]-1.732050807568877*nuVtSqSum[0]*fr[5]-1.732050807568877*nuVtSqSum[0]*fl[5]+(2.0*fr[3]-2.0*fl[3]-1.732050807568877*fr[0]-1.732050807568877*fl[0])*nuVtSqSum[3]+((-1.732050807568877*fr[1])-1.732050807568877*fl[1])*nuVtSqSum[2]-1.732050807568877*nuVtSqSum[1]*fr[2]-1.732050807568877*nuVtSqSum[1]*fl[2]); 
  incr2[13] = -0.125*(2.0*nuVtSqSum[2]*fr[15]-2.0*nuVtSqSum[2]*fl[15]+2.0*nuVtSqSum[3]*fr[14]-2.0*nuVtSqSum[3]*fl[14]+2.0*nuVtSqSum[0]*fr[13]-2.0*nuVtSqSum[0]*fl[13]-1.732050807568877*nuVtSqSum[2]*fr[12]-1.732050807568877*nuVtSqSum[2]*fl[12]+2.0*nuVtSqSum[1]*fr[10]-2.0*nuVtSqSum[1]*fl[10]-1.732050807568877*nuVtSqSum[3]*fr[9]-1.732050807568877*nuVtSqSum[3]*fl[9]-1.732050807568877*nuVtSqSum[0]*fr[8]-1.732050807568877*nuVtSqSum[0]*fl[8]-1.732050807568877*nuVtSqSum[1]*fr[4]-1.732050807568877*nuVtSqSum[1]*fl[4]); 
  incr2[14] = -0.125*(2.0*nuVtSqSum[1]*fr[15]-2.0*nuVtSqSum[1]*fl[15]+2.0*nuVtSqSum[0]*fr[14]-2.0*nuVtSqSum[0]*fl[14]+2.0*nuVtSqSum[3]*fr[13]-2.0*nuVtSqSum[3]*fl[13]-1.732050807568877*nuVtSqSum[1]*fr[12]-1.732050807568877*nuVtSqSum[1]*fl[12]+2.0*nuVtSqSum[2]*fr[10]-2.0*nuVtSqSum[2]*fl[10]-1.732050807568877*nuVtSqSum[0]*fr[9]-1.732050807568877*nuVtSqSum[0]*fl[9]-1.732050807568877*nuVtSqSum[3]*fr[8]-1.732050807568877*nuVtSqSum[3]*fl[8]-1.732050807568877*nuVtSqSum[2]*fr[4]-1.732050807568877*nuVtSqSum[2]*fl[4]); 
  incr2[15] = -0.125*(2.0*nuVtSqSum[0]*fr[15]-2.0*nuVtSqSum[0]*fl[15]+2.0*nuVtSqSum[1]*fr[14]-2.0*nuVtSqSum[1]*fl[14]+2.0*nuVtSqSum[2]*fr[13]-2.0*nuVtSqSum[2]*fl[13]-1.732050807568877*nuVtSqSum[0]*fr[12]-1.732050807568877*nuVtSqSum[0]*fl[12]+2.0*nuVtSqSum[3]*fr[10]-2.0*nuVtSqSum[3]*fl[10]-1.732050807568877*nuVtSqSum[1]*fr[9]-1.732050807568877*nuVtSqSum[1]*fl[9]-1.732050807568877*nuVtSqSum[2]*fr[8]-1.732050807568877*nuVtSqSum[2]*fl[8]-1.732050807568877*nuVtSqSum[3]*fr[4]-1.732050807568877*nuVtSqSum[3]*fl[4]); 


  Gdiff[0] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[3]*fr[11]+5.0*nuVtSqSum[3]*fl[11]+5.0*nuVtSqSum[2]*fr[7]+5.0*nuVtSqSum[2]*fl[7]+5.0*nuVtSqSum[1]*fr[6]+5.0*nuVtSqSum[1]*fl[6]+5.0*nuVtSqSum[0]*fr[3]+5.0*nuVtSqSum[0]*fl[3])-9.0*nuVtSqSum[3]*fr[5]+9.0*nuVtSqSum[3]*fl[5]+(9.0*fl[2]-9.0*fr[2])*nuVtSqSum[2]+(9.0*fl[1]-9.0*fr[1])*nuVtSqSum[1]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[0]); 
  Gdiff[1] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[2]*fr[11]+5.0*nuVtSqSum[2]*fl[11]+5.0*nuVtSqSum[3]*fr[7]+5.0*nuVtSqSum[3]*fl[7]+5.0*nuVtSqSum[0]*fr[6]+5.0*nuVtSqSum[0]*fl[6]+5.0*nuVtSqSum[1]*fr[3]+5.0*nuVtSqSum[1]*fl[3])-9.0*nuVtSqSum[2]*fr[5]+9.0*nuVtSqSum[2]*fl[5]+(9.0*fl[2]-9.0*fr[2])*nuVtSqSum[3]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[1]-9.0*nuVtSqSum[0]*fr[1]+9.0*nuVtSqSum[0]*fl[1]); 
  Gdiff[2] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[1]*fr[11]+5.0*nuVtSqSum[1]*fl[11]+5.0*nuVtSqSum[0]*fr[7]+5.0*nuVtSqSum[0]*fl[7]+5.0*nuVtSqSum[3]*fr[6]+5.0*nuVtSqSum[3]*fl[6]+5.0*nuVtSqSum[2]*fr[3]+5.0*nuVtSqSum[2]*fl[3])-9.0*nuVtSqSum[1]*fr[5]+9.0*nuVtSqSum[1]*fl[5]+(9.0*fl[1]-9.0*fr[1])*nuVtSqSum[3]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[2]-9.0*nuVtSqSum[0]*fr[2]+9.0*nuVtSqSum[0]*fl[2]); 
  Gdiff[4] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[3]*fr[15]+5.0*nuVtSqSum[3]*fl[15]+5.0*nuVtSqSum[2]*fr[14]+5.0*nuVtSqSum[2]*fl[14]+5.0*nuVtSqSum[1]*fr[13]+5.0*nuVtSqSum[1]*fl[13]+5.0*nuVtSqSum[0]*fr[10]+5.0*nuVtSqSum[0]*fl[10])-9.0*nuVtSqSum[3]*fr[12]+9.0*nuVtSqSum[3]*fl[12]-9.0*nuVtSqSum[2]*fr[9]+9.0*nuVtSqSum[2]*fl[9]-9.0*nuVtSqSum[1]*fr[8]+9.0*nuVtSqSum[1]*fl[8]-9.0*nuVtSqSum[0]*fr[4]+9.0*nuVtSqSum[0]*fl[4]); 
  Gdiff[5] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[0]*fr[11]+5.0*nuVtSqSum[0]*fl[11]+5.0*nuVtSqSum[1]*fr[7]+5.0*nuVtSqSum[1]*fl[7]+5.0*nuVtSqSum[2]*fr[6]+5.0*nuVtSqSum[2]*fl[6]+(5.0*fr[3]+5.0*fl[3])*nuVtSqSum[3])-9.0*nuVtSqSum[0]*fr[5]+9.0*nuVtSqSum[0]*fl[5]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[3]+(9.0*fl[1]-9.0*fr[1])*nuVtSqSum[2]-9.0*nuVtSqSum[1]*fr[2]+9.0*nuVtSqSum[1]*fl[2]); 
  Gdiff[8] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[2]*fr[15]+5.0*nuVtSqSum[2]*fl[15]+5.0*nuVtSqSum[3]*fr[14]+5.0*nuVtSqSum[3]*fl[14]+5.0*nuVtSqSum[0]*fr[13]+5.0*nuVtSqSum[0]*fl[13]+5.0*nuVtSqSum[1]*fr[10]+5.0*nuVtSqSum[1]*fl[10])-9.0*nuVtSqSum[2]*fr[12]+9.0*nuVtSqSum[2]*fl[12]-9.0*nuVtSqSum[3]*fr[9]+9.0*nuVtSqSum[3]*fl[9]-9.0*nuVtSqSum[0]*fr[8]+9.0*nuVtSqSum[0]*fl[8]-9.0*nuVtSqSum[1]*fr[4]+9.0*nuVtSqSum[1]*fl[4]); 
  Gdiff[9] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[1]*fr[15]+5.0*nuVtSqSum[1]*fl[15]+5.0*nuVtSqSum[0]*fr[14]+5.0*nuVtSqSum[0]*fl[14]+5.0*nuVtSqSum[3]*fr[13]+5.0*nuVtSqSum[3]*fl[13]+5.0*nuVtSqSum[2]*fr[10]+5.0*nuVtSqSum[2]*fl[10])-9.0*nuVtSqSum[1]*fr[12]+9.0*nuVtSqSum[1]*fl[12]-9.0*nuVtSqSum[0]*fr[9]+9.0*nuVtSqSum[0]*fl[9]-9.0*nuVtSqSum[3]*fr[8]+9.0*nuVtSqSum[3]*fl[8]-9.0*nuVtSqSum[2]*fr[4]+9.0*nuVtSqSum[2]*fl[4]); 
  Gdiff[12] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[0]*fr[15]+5.0*nuVtSqSum[0]*fl[15]+5.0*nuVtSqSum[1]*fr[14]+5.0*nuVtSqSum[1]*fl[14]+5.0*nuVtSqSum[2]*fr[13]+5.0*nuVtSqSum[2]*fl[13]+5.0*nuVtSqSum[3]*fr[10]+5.0*nuVtSqSum[3]*fl[10])-9.0*nuVtSqSum[0]*fr[12]+9.0*nuVtSqSum[0]*fl[12]-9.0*nuVtSqSum[1]*fr[9]+9.0*nuVtSqSum[1]*fl[9]-9.0*nuVtSqSum[2]*fr[8]+9.0*nuVtSqSum[2]*fl[8]-9.0*nuVtSqSum[3]*fr[4]+9.0*nuVtSqSum[3]*fl[4]); 

  Ghat[0] = Gdiff[0]*rdv2L+0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[2]*fUpwind[2]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+0.5*alphaDrSurf[2]*fUpwind[4]+0.5*fUpwind[2]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = Gdiff[2]*rdv2L+0.5*alphaDrSurf[1]*fUpwind[4]+0.5*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[2]+0.5*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[4] = Gdiff[4]*rdv2L+0.5*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[2]*fUpwind[6]+0.5*alphaDrSurf[1]*fUpwind[5]+0.5*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[5] = Gdiff[5]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.5*alphaDrSurf[1]*fUpwind[2]+0.5*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[8] = Gdiff[8]*rdv2L+0.5*alphaDrSurf[2]*fUpwind[7]+0.5*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[5]+0.5*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[9] = Gdiff[9]*rdv2L+0.5*alphaDrSurf[1]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[6]+0.5*alphaDrSurf[4]*fUpwind[5]+0.5*alphaDrSurf[2]*fUpwind[3]; 
  Ghat[12] = Gdiff[12]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[7]+0.5*alphaDrSurf[1]*fUpwind[6]+0.5*alphaDrSurf[2]*fUpwind[5]+0.5*fUpwind[3]*alphaDrSurf[4]; 

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

  return std::abs(wl[2]-(0.5*sumNuUx[0])/nuSum); 
} 
double VmLBOconstNuUpwindSurf2x2vTensor_VY_P1(const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:     maximum midpoint value of v-u. 
  // nuUSum[8]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:         Distribution function in left/right cells 
  // outl/outr:     Incremented distribution function in left/right cells 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = 4.0/(dxvl[3]*dxvl[3]); 
  double rdvSq4R = 4.0/(dxvr[3]*dxvr[3]); 

  const double *sumNuUy = &nuUSum[4]; 

  double alphaDrSurf[8]; 
  alphaDrSurf[0] = 0.7071067811865475*((4.0*wl[3]+2.0*dxvl[3])*nuSum-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = -1.414213562373095*sumNuUy[1]; 
  alphaDrSurf[2] = -1.414213562373095*sumNuUy[2]; 
  alphaDrSurf[4] = -1.414213562373095*sumNuUy[3]; 

  double fUpwindQuad[8];
  fUpwindQuad[0] = 0.5*(0.4330127018922193*fr[15]-0.4330127018922193*(fl[15]+fr[14]+fr[13]+fr[12])+0.4330127018922193*(fl[14]+fl[13]+fl[12])-0.25*(fr[11]+fl[11])+0.4330127018922193*(fr[10]+fr[9]+fr[8])-0.4330127018922193*(fl[10]+fl[9]+fl[8])+0.25*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3]+fr[2]+fl[2]+fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]))+0.4330127018922193*(fr[14]+fl[14]+fr[13]+fl[13]+fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9]+fr[8]+fl[8])-0.25*(fr[7]+fr[6]+fr[5])+0.25*(fl[7]+fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])+0.25*(fr[3]+fr[2]+fr[1])-0.25*(fl[3]+fl[2]+fl[1]+fr[0])+0.25*fl[0])*sgn(alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[1] = 0.5*((-0.4330127018922193*(fr[15]+fr[14]))+0.4330127018922193*(fl[15]+fl[14]+fr[13]+fr[12])-0.4330127018922193*(fl[13]+fl[12])+0.25*(fr[11]+fl[11])+0.4330127018922193*(fr[10]+fr[9])-0.4330127018922193*(fl[10]+fl[9]+fr[8])+0.4330127018922193*fl[8]+0.25*(fr[7]+fl[7])-0.25*(fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3]+fr[2]+fl[2])+0.25*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14])-0.4330127018922193*(fr[13]+fl[13]+fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9])+0.4330127018922193*(fr[8]+fl[8])-0.25*fr[7]+0.25*(fl[7]+fr[6]+fr[5])-0.25*(fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])+0.25*(fr[3]+fr[2])-0.25*(fl[3]+fl[2]+fr[1]+fr[0])+0.25*(fl[1]+fl[0]))*sgn((-alphaDrSurf[4])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[2] = 0.5*((-0.4330127018922193*fr[15])+0.4330127018922193*(fl[15]+fr[14])-0.4330127018922193*(fl[14]+fr[13])+0.4330127018922193*(fl[13]+fr[12])-0.4330127018922193*fl[12]+0.25*(fr[11]+fl[11])+0.4330127018922193*fr[10]-0.4330127018922193*(fl[10]+fr[9])+0.4330127018922193*(fl[9]+fr[8])-0.4330127018922193*fl[8]-0.25*(fr[7]+fl[7])+0.25*(fr[6]+fl[6])-0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3])+0.25*(fr[2]+fl[2])-0.25*(fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15])-0.4330127018922193*(fr[14]+fl[14])+0.4330127018922193*(fr[13]+fl[13])-0.4330127018922193*(fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10])+0.4330127018922193*(fr[9]+fl[9])-0.4330127018922193*(fr[8]+fl[8])+0.25*fr[7]-0.25*(fl[7]+fr[6])+0.25*(fl[6]+fr[5])-0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])+0.25*fr[3]-0.25*(fl[3]+fr[2])+0.25*(fl[2]+fr[1])-0.25*(fl[1]+fr[0])+0.25*fl[0])*sgn((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[3] = 0.5*(0.4330127018922193*(fr[15]+fr[14]+fr[13])-0.4330127018922193*(fl[15]+fl[14]+fl[13]+fr[12])+0.4330127018922193*fl[12]-0.25*(fr[11]+fl[11])+0.4330127018922193*fr[10]-0.4330127018922193*(fl[10]+fr[9]+fr[8])+0.4330127018922193*(fl[9]+fl[8])-0.25*(fr[7]+fl[7]+fr[6]+fl[6])+0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]-0.25*(fr[3]+fl[3])+0.25*(fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]))+0.4330127018922193*(fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]-0.4330127018922193*(fr[10]+fl[10])+0.4330127018922193*(fr[9]+fl[9]+fr[8]+fl[8])+0.25*(fr[7]+fr[6])-0.25*(fl[7]+fl[6]+fr[5])+0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])+0.25*fr[3]-0.25*(fl[3]+fr[2]+fr[1]+fr[0])+0.25*(fl[2]+fl[1]+fl[0]))*sgn(alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[4] = 0.5*((-0.4330127018922193*fr[15])+0.4330127018922193*(fl[15]+fr[14]+fr[13])-0.4330127018922193*(fl[14]+fl[13]+fr[12])+0.4330127018922193*fl[12]+0.25*(fr[11]+fl[11])-0.4330127018922193*fr[10]+0.4330127018922193*(fl[10]+fr[9]+fr[8])-0.4330127018922193*(fl[9]+fl[8])-0.25*(fr[7]+fl[7]+fr[6]+fl[6])+0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3])-0.25*(fr[2]+fl[2]+fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15])-0.4330127018922193*(fr[14]+fl[14]+fr[13]+fl[13])+0.4330127018922193*(fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10])-0.4330127018922193*(fr[9]+fl[9]+fr[8]+fl[8])+0.25*(fr[7]+fr[6])-0.25*(fl[7]+fl[6]+fr[5])+0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])-0.25*fr[3]+0.25*(fl[3]+fr[2]+fr[1])-0.25*(fl[2]+fl[1]+fr[0])+0.25*fl[0])*sgn(alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[5] = 0.5*(0.4330127018922193*(fr[15]+fr[14])-0.4330127018922193*(fl[15]+fl[14]+fr[13])+0.4330127018922193*(fl[13]+fr[12])-0.4330127018922193*fl[12]-0.25*(fr[11]+fl[11])-0.4330127018922193*fr[10]+0.4330127018922193*(fl[10]+fr[9])-0.4330127018922193*(fl[9]+fr[8])+0.4330127018922193*fl[8]-0.25*(fr[7]+fl[7])+0.25*(fr[6]+fl[6])-0.25*(fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3])-0.25*(fr[2]+fl[2])+0.25*(fr[1]+fl[1]+fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]))+0.4330127018922193*(fr[13]+fl[13])-0.4330127018922193*(fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10])-0.4330127018922193*(fr[9]+fl[9])+0.4330127018922193*(fr[8]+fl[8])+0.25*fr[7]-0.25*(fl[7]+fr[6])+0.25*(fl[6]+fr[5])-0.25*fl[5]+0.4330127018922193*(fr[4]+fl[4])-0.25*fr[3]+0.25*(fl[3]+fr[2])-0.25*(fl[2]+fr[1]+fr[0])+0.25*(fl[1]+fl[0]))*sgn((-alphaDrSurf[4])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[6] = 0.5*(0.4330127018922193*fr[15]-0.4330127018922193*(fl[15]+fr[14])+0.4330127018922193*(fl[14]+fr[13]+fr[12])-0.4330127018922193*(fl[13]+fl[12])-0.25*(fr[11]+fl[11])-0.4330127018922193*(fr[10]+fr[9])+0.4330127018922193*(fl[10]+fl[9]+fr[8])-0.4330127018922193*fl[8]+0.25*(fr[7]+fl[7])-0.25*(fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3]+fr[2]+fl[2])-0.25*(fr[1]+fl[1])+0.25*(fr[0]+fl[0]))-0.5*((-0.4330127018922193*(fr[15]+fl[15]))+0.4330127018922193*(fr[14]+fl[14])-0.4330127018922193*(fr[13]+fl[13]+fr[12]+fl[12])+0.25*fr[11]-0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9])-0.4330127018922193*(fr[8]+fl[8])-0.25*fr[7]+0.25*(fl[7]+fr[6]+fr[5])-0.25*(fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])-0.25*(fr[3]+fr[2])+0.25*(fl[3]+fl[2]+fr[1])-0.25*(fl[1]+fr[0])+0.25*fl[0])*sgn((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0]); 
  fUpwindQuad[7] = 0.5*((-0.4330127018922193*(fr[15]+fr[14]+fr[13]+fr[12]))+0.4330127018922193*(fl[15]+fl[14]+fl[13]+fl[12])+0.25*(fr[11]+fl[11])-0.4330127018922193*(fr[10]+fr[9]+fr[8])+0.4330127018922193*(fl[10]+fl[9]+fl[8])+0.25*(fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5])-0.4330127018922193*fr[4]+0.4330127018922193*fl[4]+0.25*(fr[3]+fl[3]+fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]))-0.5*(0.4330127018922193*(fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]+fr[12]+fl[12])-0.25*fr[11]+0.25*fl[11]+0.4330127018922193*(fr[10]+fl[10]+fr[9]+fl[9]+fr[8]+fl[8])-0.25*(fr[7]+fr[6]+fr[5])+0.25*(fl[7]+fl[6]+fl[5])+0.4330127018922193*(fr[4]+fl[4])-0.25*(fr[3]+fr[2]+fr[1]+fr[0])+0.25*(fl[3]+fl[2]+fl[1]+fl[0]))*sgn(alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]); 

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


  incr2[4] = -0.125*(2.0*nuVtSqSum[3]*fr[12]-2.0*nuVtSqSum[3]*fl[12]+2.0*nuVtSqSum[2]*fr[9]-2.0*nuVtSqSum[2]*fl[9]+2.0*nuVtSqSum[1]*fr[8]-2.0*nuVtSqSum[1]*fl[8]-1.732050807568877*nuVtSqSum[3]*fr[5]-1.732050807568877*nuVtSqSum[3]*fl[5]+2.0*nuVtSqSum[0]*fr[4]-2.0*nuVtSqSum[0]*fl[4]+((-1.732050807568877*fr[2])-1.732050807568877*fl[2])*nuVtSqSum[2]+((-1.732050807568877*fr[1])-1.732050807568877*fl[1])*nuVtSqSum[1]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[0]); 
  incr2[8] = -0.125*(2.0*nuVtSqSum[2]*fr[12]-2.0*nuVtSqSum[2]*fl[12]+2.0*nuVtSqSum[3]*fr[9]-2.0*nuVtSqSum[3]*fl[9]+2.0*nuVtSqSum[0]*fr[8]-2.0*nuVtSqSum[0]*fl[8]-1.732050807568877*nuVtSqSum[2]*fr[5]-1.732050807568877*nuVtSqSum[2]*fl[5]+2.0*nuVtSqSum[1]*fr[4]-2.0*nuVtSqSum[1]*fl[4]+((-1.732050807568877*fr[2])-1.732050807568877*fl[2])*nuVtSqSum[3]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[1]-1.732050807568877*nuVtSqSum[0]*fr[1]-1.732050807568877*nuVtSqSum[0]*fl[1]); 
  incr2[9] = -0.125*(2.0*nuVtSqSum[1]*fr[12]-2.0*nuVtSqSum[1]*fl[12]+2.0*nuVtSqSum[0]*fr[9]-2.0*nuVtSqSum[0]*fl[9]+2.0*nuVtSqSum[3]*fr[8]-2.0*nuVtSqSum[3]*fl[8]-1.732050807568877*nuVtSqSum[1]*fr[5]-1.732050807568877*nuVtSqSum[1]*fl[5]+2.0*nuVtSqSum[2]*fr[4]-2.0*nuVtSqSum[2]*fl[4]+((-1.732050807568877*fr[1])-1.732050807568877*fl[1])*nuVtSqSum[3]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[2]-1.732050807568877*nuVtSqSum[0]*fr[2]-1.732050807568877*nuVtSqSum[0]*fl[2]); 
  incr2[10] = -0.125*(2.0*nuVtSqSum[3]*fr[15]-2.0*nuVtSqSum[3]*fl[15]+2.0*nuVtSqSum[2]*fr[14]-2.0*nuVtSqSum[2]*fl[14]+2.0*nuVtSqSum[1]*fr[13]-2.0*nuVtSqSum[1]*fl[13]-1.732050807568877*nuVtSqSum[3]*fr[11]-1.732050807568877*nuVtSqSum[3]*fl[11]+2.0*nuVtSqSum[0]*fr[10]-2.0*nuVtSqSum[0]*fl[10]-1.732050807568877*nuVtSqSum[2]*fr[7]-1.732050807568877*nuVtSqSum[2]*fl[7]-1.732050807568877*nuVtSqSum[1]*fr[6]-1.732050807568877*nuVtSqSum[1]*fl[6]-1.732050807568877*nuVtSqSum[0]*fr[3]-1.732050807568877*nuVtSqSum[0]*fl[3]); 
  incr2[12] = -0.125*(2.0*nuVtSqSum[0]*fr[12]-2.0*nuVtSqSum[0]*fl[12]+2.0*nuVtSqSum[1]*fr[9]-2.0*nuVtSqSum[1]*fl[9]+2.0*nuVtSqSum[2]*fr[8]-2.0*nuVtSqSum[2]*fl[8]-1.732050807568877*nuVtSqSum[0]*fr[5]-1.732050807568877*nuVtSqSum[0]*fl[5]+2.0*nuVtSqSum[3]*fr[4]-2.0*nuVtSqSum[3]*fl[4]+((-1.732050807568877*fr[0])-1.732050807568877*fl[0])*nuVtSqSum[3]+((-1.732050807568877*fr[1])-1.732050807568877*fl[1])*nuVtSqSum[2]-1.732050807568877*nuVtSqSum[1]*fr[2]-1.732050807568877*nuVtSqSum[1]*fl[2]); 
  incr2[13] = -0.125*(2.0*nuVtSqSum[2]*fr[15]-2.0*nuVtSqSum[2]*fl[15]+2.0*nuVtSqSum[3]*fr[14]-2.0*nuVtSqSum[3]*fl[14]+2.0*nuVtSqSum[0]*fr[13]-2.0*nuVtSqSum[0]*fl[13]-1.732050807568877*nuVtSqSum[2]*fr[11]-1.732050807568877*nuVtSqSum[2]*fl[11]+2.0*nuVtSqSum[1]*fr[10]-2.0*nuVtSqSum[1]*fl[10]-1.732050807568877*nuVtSqSum[3]*fr[7]-1.732050807568877*nuVtSqSum[3]*fl[7]-1.732050807568877*nuVtSqSum[0]*fr[6]-1.732050807568877*nuVtSqSum[0]*fl[6]-1.732050807568877*nuVtSqSum[1]*fr[3]-1.732050807568877*nuVtSqSum[1]*fl[3]); 
  incr2[14] = -0.125*(2.0*nuVtSqSum[1]*fr[15]-2.0*nuVtSqSum[1]*fl[15]+2.0*nuVtSqSum[0]*fr[14]-2.0*nuVtSqSum[0]*fl[14]+2.0*nuVtSqSum[3]*fr[13]-2.0*nuVtSqSum[3]*fl[13]-1.732050807568877*nuVtSqSum[1]*fr[11]-1.732050807568877*nuVtSqSum[1]*fl[11]+2.0*nuVtSqSum[2]*fr[10]-2.0*nuVtSqSum[2]*fl[10]-1.732050807568877*nuVtSqSum[0]*fr[7]-1.732050807568877*nuVtSqSum[0]*fl[7]-1.732050807568877*nuVtSqSum[3]*fr[6]-1.732050807568877*nuVtSqSum[3]*fl[6]-1.732050807568877*nuVtSqSum[2]*fr[3]-1.732050807568877*nuVtSqSum[2]*fl[3]); 
  incr2[15] = -0.125*(2.0*nuVtSqSum[0]*fr[15]-2.0*nuVtSqSum[0]*fl[15]+2.0*nuVtSqSum[1]*fr[14]-2.0*nuVtSqSum[1]*fl[14]+2.0*nuVtSqSum[2]*fr[13]-2.0*nuVtSqSum[2]*fl[13]-1.732050807568877*nuVtSqSum[0]*fr[11]-1.732050807568877*nuVtSqSum[0]*fl[11]+2.0*nuVtSqSum[3]*fr[10]-2.0*nuVtSqSum[3]*fl[10]-1.732050807568877*nuVtSqSum[1]*fr[7]-1.732050807568877*nuVtSqSum[1]*fl[7]-1.732050807568877*nuVtSqSum[2]*fr[6]-1.732050807568877*nuVtSqSum[2]*fl[6]+((-1.732050807568877*fr[3])-1.732050807568877*fl[3])*nuVtSqSum[3]); 


  Gdiff[0] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[3]*fr[12]+5.0*nuVtSqSum[3]*fl[12]+5.0*nuVtSqSum[2]*fr[9]+5.0*nuVtSqSum[2]*fl[9]+5.0*nuVtSqSum[1]*fr[8]+5.0*nuVtSqSum[1]*fl[8]+5.0*nuVtSqSum[0]*fr[4]+5.0*nuVtSqSum[0]*fl[4])-9.0*nuVtSqSum[3]*fr[5]+9.0*nuVtSqSum[3]*fl[5]+(9.0*fl[2]-9.0*fr[2])*nuVtSqSum[2]+(9.0*fl[1]-9.0*fr[1])*nuVtSqSum[1]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[0]); 
  Gdiff[1] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[2]*fr[12]+5.0*nuVtSqSum[2]*fl[12]+5.0*nuVtSqSum[3]*fr[9]+5.0*nuVtSqSum[3]*fl[9]+5.0*nuVtSqSum[0]*fr[8]+5.0*nuVtSqSum[0]*fl[8]+5.0*nuVtSqSum[1]*fr[4]+5.0*nuVtSqSum[1]*fl[4])-9.0*nuVtSqSum[2]*fr[5]+9.0*nuVtSqSum[2]*fl[5]+(9.0*fl[2]-9.0*fr[2])*nuVtSqSum[3]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[1]-9.0*nuVtSqSum[0]*fr[1]+9.0*nuVtSqSum[0]*fl[1]); 
  Gdiff[2] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[1]*fr[12]+5.0*nuVtSqSum[1]*fl[12]+5.0*nuVtSqSum[0]*fr[9]+5.0*nuVtSqSum[0]*fl[9]+5.0*nuVtSqSum[3]*fr[8]+5.0*nuVtSqSum[3]*fl[8]+5.0*nuVtSqSum[2]*fr[4]+5.0*nuVtSqSum[2]*fl[4])-9.0*nuVtSqSum[1]*fr[5]+9.0*nuVtSqSum[1]*fl[5]+(9.0*fl[1]-9.0*fr[1])*nuVtSqSum[3]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[2]-9.0*nuVtSqSum[0]*fr[2]+9.0*nuVtSqSum[0]*fl[2]); 
  Gdiff[3] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[3]*fr[15]+5.0*nuVtSqSum[3]*fl[15]+5.0*nuVtSqSum[2]*fr[14]+5.0*nuVtSqSum[2]*fl[14]+5.0*nuVtSqSum[1]*fr[13]+5.0*nuVtSqSum[1]*fl[13]+5.0*nuVtSqSum[0]*fr[10]+5.0*nuVtSqSum[0]*fl[10])-9.0*nuVtSqSum[3]*fr[11]+9.0*nuVtSqSum[3]*fl[11]-9.0*nuVtSqSum[2]*fr[7]+9.0*nuVtSqSum[2]*fl[7]-9.0*nuVtSqSum[1]*fr[6]+9.0*nuVtSqSum[1]*fl[6]-9.0*nuVtSqSum[0]*fr[3]+9.0*nuVtSqSum[0]*fl[3]); 
  Gdiff[5] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[0]*fr[12]+5.0*nuVtSqSum[0]*fl[12]+5.0*nuVtSqSum[1]*fr[9]+5.0*nuVtSqSum[1]*fl[9]+5.0*nuVtSqSum[2]*fr[8]+5.0*nuVtSqSum[2]*fl[8]+5.0*nuVtSqSum[3]*fr[4]+5.0*nuVtSqSum[3]*fl[4])-9.0*nuVtSqSum[0]*fr[5]+9.0*nuVtSqSum[0]*fl[5]+(9.0*fl[0]-9.0*fr[0])*nuVtSqSum[3]+(9.0*fl[1]-9.0*fr[1])*nuVtSqSum[2]-9.0*nuVtSqSum[1]*fr[2]+9.0*nuVtSqSum[1]*fl[2]); 
  Gdiff[6] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[2]*fr[15]+5.0*nuVtSqSum[2]*fl[15]+5.0*nuVtSqSum[3]*fr[14]+5.0*nuVtSqSum[3]*fl[14]+5.0*nuVtSqSum[0]*fr[13]+5.0*nuVtSqSum[0]*fl[13]+5.0*nuVtSqSum[1]*fr[10]+5.0*nuVtSqSum[1]*fl[10])-9.0*nuVtSqSum[2]*fr[11]+9.0*nuVtSqSum[2]*fl[11]-9.0*nuVtSqSum[3]*fr[7]+9.0*nuVtSqSum[3]*fl[7]-9.0*nuVtSqSum[0]*fr[6]+9.0*nuVtSqSum[0]*fl[6]-9.0*nuVtSqSum[1]*fr[3]+9.0*nuVtSqSum[1]*fl[3]); 
  Gdiff[7] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[1]*fr[15]+5.0*nuVtSqSum[1]*fl[15]+5.0*nuVtSqSum[0]*fr[14]+5.0*nuVtSqSum[0]*fl[14]+5.0*nuVtSqSum[3]*fr[13]+5.0*nuVtSqSum[3]*fl[13]+5.0*nuVtSqSum[2]*fr[10]+5.0*nuVtSqSum[2]*fl[10])-9.0*nuVtSqSum[1]*fr[11]+9.0*nuVtSqSum[1]*fl[11]-9.0*nuVtSqSum[0]*fr[7]+9.0*nuVtSqSum[0]*fl[7]-9.0*nuVtSqSum[3]*fr[6]+9.0*nuVtSqSum[3]*fl[6]-9.0*nuVtSqSum[2]*fr[3]+9.0*nuVtSqSum[2]*fl[3]); 
  Gdiff[11] = -0.0625*(1.732050807568877*(5.0*nuVtSqSum[0]*fr[15]+5.0*nuVtSqSum[0]*fl[15]+5.0*nuVtSqSum[1]*fr[14]+5.0*nuVtSqSum[1]*fl[14]+5.0*nuVtSqSum[2]*fr[13]+5.0*nuVtSqSum[2]*fl[13]+5.0*nuVtSqSum[3]*fr[10]+5.0*nuVtSqSum[3]*fl[10])-9.0*nuVtSqSum[0]*fr[11]+9.0*nuVtSqSum[0]*fl[11]-9.0*nuVtSqSum[1]*fr[7]+9.0*nuVtSqSum[1]*fl[7]-9.0*nuVtSqSum[2]*fr[6]+9.0*nuVtSqSum[2]*fl[6]+(9.0*fl[3]-9.0*fr[3])*nuVtSqSum[3]); 

  Ghat[0] = Gdiff[0]*rdv2L+0.5*alphaDrSurf[4]*fUpwind[4]+0.5*alphaDrSurf[2]*fUpwind[2]+0.5*alphaDrSurf[1]*fUpwind[1]+0.5*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2L+0.5*alphaDrSurf[2]*fUpwind[4]+0.5*fUpwind[2]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[1]+0.5*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = Gdiff[2]*rdv2L+0.5*alphaDrSurf[1]*fUpwind[4]+0.5*fUpwind[1]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fUpwind[2]+0.5*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = Gdiff[3]*rdv2L+0.5*alphaDrSurf[4]*fUpwind[7]+0.5*alphaDrSurf[2]*fUpwind[6]+0.5*alphaDrSurf[1]*fUpwind[5]+0.5*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[5] = Gdiff[5]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[4]+0.5*fUpwind[0]*alphaDrSurf[4]+0.5*alphaDrSurf[1]*fUpwind[2]+0.5*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[6] = Gdiff[6]*rdv2L+0.5*alphaDrSurf[2]*fUpwind[7]+0.5*alphaDrSurf[4]*fUpwind[6]+0.5*alphaDrSurf[0]*fUpwind[5]+0.5*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[7] = Gdiff[7]*rdv2L+0.5*alphaDrSurf[1]*fUpwind[7]+0.5*alphaDrSurf[0]*fUpwind[6]+0.5*alphaDrSurf[4]*fUpwind[5]+0.5*alphaDrSurf[2]*fUpwind[3]; 
  Ghat[11] = Gdiff[11]*rdv2L+0.5*alphaDrSurf[0]*fUpwind[7]+0.5*alphaDrSurf[1]*fUpwind[6]+0.5*alphaDrSurf[2]*fUpwind[5]+0.5*fUpwind[3]*alphaDrSurf[4]; 

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

  return std::abs(wl[3]-(0.5*sumNuUy[0])/nuSum); 
} 
