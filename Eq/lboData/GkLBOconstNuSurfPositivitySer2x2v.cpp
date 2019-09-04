#include <GkLBOModDecl.h> 
double GkLBOconstNuSurfPositivity2x2vSer_Vpar_P1(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
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

  double f0Quad[8]; 
  double f1Quad[8]; 
  double limQuad[8]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = (-0.3535533905932737*alphaDrSurf[4])+0.3535533905932737*(alphaDrSurf[2]+alphaDrSurf[1])-0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[0] = 0.1767766952966368*((fr[12]-1.0*(fl[12]+fr[9])+fl[9]-1.0*fr[8]+fl[8]-1.0*fr[5]+fl[5]+fr[4]-1.0*fl[4]+fr[2]-1.0*fl[2]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12])+fr[9]+fl[9]+fr[8]+fl[8]+fr[5]+fl[5]-1.0*(fr[4]+fl[4]+fr[2]+fl[2]+fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[0] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14]+fr[13]+fl[13]+fr[11]+fl[11])+fr[10]+fl[10]+fr[7]+fl[7]+fr[6]+fl[6]-1.0*(fr[3]+fl[3]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*fl[14]+fr[13]-1.0*fl[13]+fr[11]-1.0*(fl[11]+fr[10])+fl[10]-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]+fr[3]-1.0*fl[3]); 
  limQuad[0] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = 0.3535533905932737*(alphaDrSurf[4]+alphaDrSurf[2])-0.3535533905932737*(alphaDrSurf[1]+alphaDrSurf[0]); 
  f0Quad[1] = -0.1767766952966368*((fr[12]-1.0*fl[12]+fr[9]-1.0*(fl[9]+fr[8])+fl[8]-1.0*fr[5]+fl[5]-1.0*fr[4]+fl[4]-1.0*fr[2]+fl[2]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12]+fr[9]+fl[9])+fr[8]+fl[8]+fr[5]+fl[5]+fr[4]+fl[4]+fr[2]+fl[2]-1.0*(fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[1] = -0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]-1.0*(fr[13]+fl[13]+fr[11]+fl[11]+fr[10]+fl[10]+fr[7]+fl[7])+fr[6]+fl[6]+fr[3]+fl[3])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]+fr[13]-1.0*fl[13]+fr[11]-1.0*fl[11]+fr[10]-1.0*fl[10]+fr[7]-1.0*(fl[7]+fr[6])+fl[6]-1.0*fr[3]+fl[3]); 
  limQuad[1] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = 0.3535533905932737*alphaDrSurf[4]-0.3535533905932737*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[1]-0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[2] = -0.1767766952966368*((fr[12]-1.0*(fl[12]+fr[9])+fl[9]+fr[8]-1.0*(fl[8]+fr[5])+fl[5]-1.0*fr[4]+fl[4]+fr[2]-1.0*(fl[2]+fr[1])+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12])+fr[9]+fl[9]-1.0*(fr[8]+fl[8])+fr[5]+fl[5]+fr[4]+fl[4]-1.0*(fr[2]+fl[2])+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[2] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14])+fr[13]+fl[13]-1.0*(fr[11]+fl[11]+fr[10]+fl[10])+fr[7]+fl[7]-1.0*(fr[6]+fl[6])+fr[3]+fl[3])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*(fl[14]+fr[13])+fl[13]+fr[11]-1.0*fl[11]+fr[10]-1.0*(fl[10]+fr[7])+fl[7]+fr[6]-1.0*(fl[6]+fr[3])+fl[3]); 
  limQuad[2] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = -0.3535533905932737*(alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]); 
  f0Quad[3] = 0.1767766952966368*((fr[12]-1.0*fl[12]+fr[9]-1.0*fl[9]+fr[8]-1.0*(fl[8]+fr[5])+fl[5]+fr[4]-1.0*(fl[4]+fr[2])+fl[2]-1.0*fr[1]+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12]+fr[9]+fl[9]+fr[8]+fl[8])+fr[5]+fl[5]-1.0*(fr[4]+fl[4])+fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[3] = 0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]-1.0*(fr[11]+fl[11])+fr[10]+fl[10]-1.0*(fr[7]+fl[7]+fr[6]+fl[6]+fr[3]+fl[3]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]-1.0*fr[13]+fl[13]+fr[11]-1.0*(fl[11]+fr[10])+fl[10]+fr[7]-1.0*fl[7]+fr[6]-1.0*fl[6]+fr[3]-1.0*fl[3]); 
  limQuad[3] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = (-0.3535533905932737*alphaDrSurf[4])+0.3535533905932737*(alphaDrSurf[2]+alphaDrSurf[1])-0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[4] = -0.1767766952966368*((fr[12]-1.0*(fl[12]+fr[9])+fl[9]-1.0*fr[8]+fl[8]+fr[5]-1.0*fl[5]+fr[4]-1.0*(fl[4]+fr[2])+fl[2]-1.0*fr[1]+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12])+fr[9]+fl[9]+fr[8]+fl[8]-1.0*(fr[5]+fl[5]+fr[4]+fl[4])+fr[2]+fl[2]+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[4] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14]+fr[13]+fl[13])+fr[11]+fl[11]+fr[10]+fl[10]-1.0*(fr[7]+fl[7]+fr[6]+fl[6])+fr[3]+fl[3])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*fl[14]+fr[13]-1.0*(fl[13]+fr[11])+fl[11]-1.0*fr[10]+fl[10]+fr[7]-1.0*fl[7]+fr[6]-1.0*(fl[6]+fr[3])+fl[3]); 
  limQuad[4] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = 0.3535533905932737*(alphaDrSurf[4]+alphaDrSurf[2])-0.3535533905932737*(alphaDrSurf[1]+alphaDrSurf[0]); 
  f0Quad[5] = 0.1767766952966368*((fr[12]-1.0*fl[12]+fr[9]-1.0*(fl[9]+fr[8])+fl[8]+fr[5]-1.0*(fl[5]+fr[4])+fl[4]+fr[2]-1.0*(fl[2]+fr[1])+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12]+fr[9]+fl[9])+fr[8]+fl[8]-1.0*(fr[5]+fl[5])+fr[4]+fl[4]-1.0*(fr[2]+fl[2])+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[5] = 0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]-1.0*(fr[13]+fl[13])+fr[11]+fl[11]-1.0*(fr[10]+fl[10])+fr[7]+fl[7]-1.0*(fr[6]+fl[6]+fr[3]+fl[3]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]+fr[13]-1.0*(fl[13]+fr[11])+fl[11]+fr[10]-1.0*(fl[10]+fr[7])+fl[7]+fr[6]-1.0*fl[6]+fr[3]-1.0*fl[3]); 
  limQuad[5] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = 0.3535533905932737*alphaDrSurf[4]-0.3535533905932737*alphaDrSurf[2]+0.3535533905932737*alphaDrSurf[1]-0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[6] = 0.1767766952966368*((fr[12]-1.0*(fl[12]+fr[9])+fl[9]+fr[8]-1.0*fl[8]+fr[5]-1.0*(fl[5]+fr[4])+fl[4]-1.0*fr[2]+fl[2]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12])+fr[9]+fl[9]-1.0*(fr[8]+fl[8]+fr[5]+fl[5])+fr[4]+fl[4]+fr[2]+fl[2]-1.0*(fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[6] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14])+fr[13]+fl[13]+fr[11]+fl[11]-1.0*(fr[10]+fl[10]+fr[7]+fl[7])+fr[6]+fl[6]-1.0*(fr[3]+fl[3]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*(fl[14]+fr[13])+fl[13]-1.0*fr[11]+fl[11]+fr[10]-1.0*fl[10]+fr[7]-1.0*(fl[7]+fr[6])+fl[6]+fr[3]-1.0*fl[3]); 
  limQuad[6] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = -0.3535533905932737*(alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0]); 
  f0Quad[7] = -0.1767766952966368*((fr[12]-1.0*fl[12]+fr[9]-1.0*fl[9]+fr[8]-1.0*fl[8]+fr[5]-1.0*fl[5]+fr[4]-1.0*fl[4]+fr[2]-1.0*fl[2]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[12]+fl[12]+fr[9]+fl[9]+fr[8]+fl[8]+fr[5]+fl[5]+fr[4]+fl[4]+fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[7] = -0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]+fr[11]+fl[11]+fr[10]+fl[10]+fr[7]+fl[7]+fr[6]+fl[6]+fr[3]+fl[3])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]-1.0*fr[13]+fl[13]-1.0*fr[11]+fl[11]-1.0*fr[10]+fl[10]-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]-1.0*fr[3]+fl[3]); 
  limQuad[7] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  double fhat[16]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[3] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[4] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[5] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[6] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[7] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[8] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[9] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[10] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[11] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[12] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[13] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[14] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[15] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[8];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrl[0] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[14]+fhat[13]+fhat[11]+3.0*fhat[3])+9.0*(fhat[10]+fhat[7]+fhat[6])))/(36.0*EPSILON-1.732050807568877*(fhat[12]+3.0*(fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[9]+fhat[8]+fhat[5]+3.0*fhat[0])); 
  rCtrl[1] = (3.0*fhat[15]+5.196152422706631*(fhat[14]-1.0*(fhat[13]+fhat[11])+3.0*fhat[3])+9.0*(fhat[6]-1.0*(fhat[10]+fhat[7])))/(36.0*EPSILON+1.732050807568877*(fhat[12]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[2])))+3.0*(fhat[9]-1.0*(fhat[8]+fhat[5]-3.0*fhat[0]))); 
  rCtrl[2] = (3.0*fhat[15]-5.196152422706631*(fhat[14]-1.0*fhat[13]+fhat[11]-3.0*fhat[3])+9.0*((-1.0*fhat[10])+fhat[7]-1.0*fhat[6]))/(36.0*EPSILON+1.732050807568877*(fhat[12]+3.0*((-1.0*fhat[4])+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*fhat[9])+fhat[8]-1.0*fhat[5]+3.0*fhat[0])); 
  rCtrl[3] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[14]+fhat[13]-1.0*(fhat[11]+3.0*fhat[3]))+9.0*(fhat[10]-1.0*(fhat[7]+fhat[6]))))/(36.0*EPSILON-1.732050807568877*(fhat[12]+3.0*(fhat[4]-1.0*(fhat[2]+fhat[1])))+3.0*((-1.0*(fhat[9]+fhat[8]))+fhat[5]+3.0*fhat[0])); 
  rCtrl[4] = (3.0*fhat[15]-5.196152422706631*(fhat[14]+fhat[13]-1.0*(fhat[11]+3.0*fhat[3]))+9.0*(fhat[10]-1.0*(fhat[7]+fhat[6])))/(36.0*EPSILON+1.732050807568877*(fhat[12]+3.0*(fhat[4]-1.0*(fhat[2]+fhat[1])))+3.0*((-1.0*(fhat[9]+fhat[8]))+fhat[5]+3.0*fhat[0])); 
  rCtrl[5] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[14]-1.0*fhat[13]+fhat[11]-3.0*fhat[3])+9.0*((-1.0*fhat[10])+fhat[7]-1.0*fhat[6])))/(36.0*EPSILON-1.732050807568877*(fhat[12]+3.0*((-1.0*fhat[4])+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*fhat[9])+fhat[8]-1.0*fhat[5]+3.0*fhat[0])); 
  rCtrl[6] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[14]-1.0*(fhat[13]+fhat[11])+3.0*fhat[3])+9.0*(fhat[6]-1.0*(fhat[10]+fhat[7]))))/(36.0*EPSILON-1.732050807568877*(fhat[12]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[2])))+3.0*(fhat[9]-1.0*(fhat[8]+fhat[5]-3.0*fhat[0]))); 
  rCtrl[7] = (3.0*fhat[15]+5.196152422706631*(fhat[14]+fhat[13]+fhat[11]+3.0*fhat[3])+9.0*(fhat[10]+fhat[7]+fhat[6]))/(36.0*EPSILON+1.732050807568877*(fhat[12]+3.0*(fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[9]+fhat[8]+fhat[5]+3.0*fhat[0])); 
  double fhatCtrl[8];  // fhatCtrl = anti-limited fhat evaluated at each control node on vx surface 
  fhatCtrl[0] = -0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[4]+fhat[2]+fhat[1]))-3.0*(fhat[9]+fhat[8]+fhat[5]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[2])))+3.0*(fhat[9]-1.0*(fhat[8]+fhat[5])+3.0*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*((-1.0*fhat[4])+fhat[2]-1.0*fhat[1]))-3.0*(fhat[9]-1.0*fhat[8]+fhat[5]-3.0*fhat[0]))*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = -0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[4]-1.0*(fhat[2]+fhat[1])))+3.0*(fhat[9]+fhat[8]-1.0*(fhat[5]+3.0*fhat[0])))*limTheta(rCtrl[3],-1.0); 
  fhatCtrl[4] = 0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[4]-1.0*(fhat[2]+fhat[1])))-3.0*(fhat[9]+fhat[8]-1.0*(fhat[5]+3.0*fhat[0])))*limTheta(rCtrl[4],-1.0); 
  fhatCtrl[5] = -0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*((-1.0*fhat[4])+fhat[2]-1.0*fhat[1]))+3.0*(fhat[9]-1.0*fhat[8]+fhat[5]-3.0*fhat[0]))*limTheta(rCtrl[5],-1.0); 
  fhatCtrl[6] = -0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[1]-1.0*(fhat[4]+fhat[2])))-3.0*(fhat[9]-1.0*(fhat[8]+fhat[5])+3.0*fhat[0]))*limTheta(rCtrl[6],-1.0); 
  fhatCtrl[7] = 0.02777777777777778*(1.732050807568877*(fhat[12]+3.0*(fhat[4]+fhat[2]+fhat[1]))+3.0*(fhat[9]+fhat[8]+fhat[5]+3.0*fhat[0]))*limTheta(rCtrl[7],-1.0); 
  double fhatAL[8];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.3535533905932737*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.6123724356957944*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 1.060660171779821*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[5] = 1.060660171779821*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[6] = 1.060660171779821*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 1.837117307087383*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[8]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.25*(1.414213562373095*(1.0*(fhatAL[6]+fhatAL[4])+fhatAL[5]+fhatAL[0])-1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*fhatAL[1]-0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[2]))+1.414213562373095*((-0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[6]))-1.0*fhatAL[5]+fhatAL[0])), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[7]+fhatAL[2])-0.5773502691896258*(fhatAL[3]+fhatAL[1]))+1.414213562373095*((-1.0*(fhatAL[6]+fhatAL[4]))+fhatAL[5]+fhatAL[0])), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[2]-1.0*fhatAL[7])-0.5773502691896258*fhatAL[3]+0.5773502691896258*fhatAL[1])+1.414213562373095*(0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[6])-1.0*fhatAL[5]+fhatAL[0])), limQuad[3])); 
  fhatALQuad[4] = std::max(0.0, std::min(0.25*(2.449489742783178*((-0.5773502691896258*(fhatAL[2]-1.0*fhatAL[7]))+0.5773502691896258*fhatAL[3]-0.5773502691896258*fhatAL[1])+1.414213562373095*((-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[4]))-1.0*fhatAL[5]+fhatAL[0])), limQuad[4])); 
  fhatALQuad[5] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[3]+fhatAL[1])-0.5773502691896258*(fhatAL[7]+fhatAL[2]))+1.414213562373095*((-1.0*(fhatAL[6]+fhatAL[4]))+fhatAL[5]+fhatAL[0])), limQuad[5])); 
  fhatALQuad[6] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[2])-0.5773502691896258*fhatAL[1])+1.414213562373095*(0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[4])-1.0*fhatAL[5]+fhatAL[0])), limQuad[6])); 
  fhatALQuad[7] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])+1.414213562373095*(1.0*(fhatAL[6]+fhatAL[4])+fhatAL[5]+fhatAL[0])), limQuad[7])); 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double Ghat[16]; 
  for(unsigned short int i=0; i<16; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.082531754730548*nuVtSqSum[3]*fr[11])-1.082531754730548*nuVtSqSum[3]*fl[11]-1.082531754730548*nuVtSqSum[2]*fr[7]-1.082531754730548*nuVtSqSum[2]*fl[7]-1.082531754730548*nuVtSqSum[1]*fr[6]-1.082531754730548*nuVtSqSum[1]*fl[6]+1.125*nuVtSqSum[3]*fr[5]-1.125*nuVtSqSum[3]*fl[5]-1.082531754730548*nuVtSqSum[0]*fr[3]-1.082531754730548*nuVtSqSum[0]*fl[3]+1.125*fr[2]*nuVtSqSum[2]-1.125*fl[2]*nuVtSqSum[2]+1.125*fr[1]*nuVtSqSum[1]-1.125*fl[1]*nuVtSqSum[1]+1.125*fr[0]*nuVtSqSum[0]-1.125*fl[0]*nuVtSqSum[0])*rdv+0.5*alphaDrSurf[4]*fhatAL[4]+0.5*alphaDrSurf[2]*fhatAL[2]+0.5*alphaDrSurf[1]*fhatAL[1]+0.5*alphaDrSurf[0]*fhatAL[0]; 
  Ghat[1] = ((-1.082531754730548*nuVtSqSum[2]*fr[11])-1.082531754730548*nuVtSqSum[2]*fl[11]-1.082531754730548*nuVtSqSum[3]*fr[7]-1.082531754730548*nuVtSqSum[3]*fl[7]-1.082531754730548*nuVtSqSum[0]*fr[6]-1.082531754730548*nuVtSqSum[0]*fl[6]+1.125*nuVtSqSum[2]*fr[5]-1.125*nuVtSqSum[2]*fl[5]+1.125*fr[2]*nuVtSqSum[3]-1.125*fl[2]*nuVtSqSum[3]-1.082531754730548*nuVtSqSum[1]*fr[3]-1.082531754730548*nuVtSqSum[1]*fl[3]+1.125*fr[0]*nuVtSqSum[1]-1.125*fl[0]*nuVtSqSum[1]+1.125*nuVtSqSum[0]*fr[1]-1.125*nuVtSqSum[0]*fl[1])*rdv+0.5*alphaDrSurf[2]*fhatAL[4]+0.5*fhatAL[2]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fhatAL[1]+0.5*fhatAL[0]*alphaDrSurf[1]; 
  Ghat[2] = ((-1.082531754730548*nuVtSqSum[1]*fr[11])-1.082531754730548*nuVtSqSum[1]*fl[11]-1.082531754730548*nuVtSqSum[0]*fr[7]-1.082531754730548*nuVtSqSum[0]*fl[7]-1.082531754730548*nuVtSqSum[3]*fr[6]-1.082531754730548*nuVtSqSum[3]*fl[6]+1.125*nuVtSqSum[1]*fr[5]-1.125*nuVtSqSum[1]*fl[5]+1.125*fr[1]*nuVtSqSum[3]-1.125*fl[1]*nuVtSqSum[3]-1.082531754730548*nuVtSqSum[2]*fr[3]-1.082531754730548*nuVtSqSum[2]*fl[3]+1.125*fr[0]*nuVtSqSum[2]-1.125*fl[0]*nuVtSqSum[2]+1.125*nuVtSqSum[0]*fr[2]-1.125*nuVtSqSum[0]*fl[2])*rdv+0.5*alphaDrSurf[1]*fhatAL[4]+0.5*fhatAL[1]*alphaDrSurf[4]+0.5*alphaDrSurf[0]*fhatAL[2]+0.5*fhatAL[0]*alphaDrSurf[2]; 
  Ghat[4] = ((-1.082531754730548*nuVtSqSum[3]*fr[15])-1.082531754730548*nuVtSqSum[3]*fl[15]-1.082531754730548*nuVtSqSum[2]*fr[14]-1.082531754730548*nuVtSqSum[2]*fl[14]-1.082531754730548*nuVtSqSum[1]*fr[13]-1.082531754730548*nuVtSqSum[1]*fl[13]+1.125*nuVtSqSum[3]*fr[12]-1.125*nuVtSqSum[3]*fl[12]-1.082531754730548*nuVtSqSum[0]*fr[10]-1.082531754730548*nuVtSqSum[0]*fl[10]+1.125*nuVtSqSum[2]*fr[9]-1.125*nuVtSqSum[2]*fl[9]+1.125*nuVtSqSum[1]*fr[8]-1.125*nuVtSqSum[1]*fl[8]+1.125*nuVtSqSum[0]*fr[4]-1.125*nuVtSqSum[0]*fl[4])*rdv+0.5*alphaDrSurf[4]*fhatAL[7]+0.5*alphaDrSurf[2]*fhatAL[6]+0.5*alphaDrSurf[1]*fhatAL[5]+0.5*alphaDrSurf[0]*fhatAL[3]; 
  Ghat[5] = ((-1.082531754730548*nuVtSqSum[0]*fr[11])-1.082531754730548*nuVtSqSum[0]*fl[11]-1.082531754730548*nuVtSqSum[1]*fr[7]-1.082531754730548*nuVtSqSum[1]*fl[7]-1.082531754730548*nuVtSqSum[2]*fr[6]-1.082531754730548*nuVtSqSum[2]*fl[6]+1.125*nuVtSqSum[0]*fr[5]-1.125*nuVtSqSum[0]*fl[5]-1.082531754730548*fr[3]*nuVtSqSum[3]-1.082531754730548*fl[3]*nuVtSqSum[3]+1.125*fr[0]*nuVtSqSum[3]-1.125*fl[0]*nuVtSqSum[3]+1.125*fr[1]*nuVtSqSum[2]-1.125*fl[1]*nuVtSqSum[2]+1.125*nuVtSqSum[1]*fr[2]-1.125*nuVtSqSum[1]*fl[2])*rdv+0.5*alphaDrSurf[0]*fhatAL[4]+0.5*fhatAL[0]*alphaDrSurf[4]+0.5*alphaDrSurf[1]*fhatAL[2]+0.5*fhatAL[1]*alphaDrSurf[2]; 
  Ghat[8] = ((-1.082531754730548*nuVtSqSum[2]*fr[15])-1.082531754730548*nuVtSqSum[2]*fl[15]-1.082531754730548*nuVtSqSum[3]*fr[14]-1.082531754730548*nuVtSqSum[3]*fl[14]-1.082531754730548*nuVtSqSum[0]*fr[13]-1.082531754730548*nuVtSqSum[0]*fl[13]+1.125*nuVtSqSum[2]*fr[12]-1.125*nuVtSqSum[2]*fl[12]-1.082531754730548*nuVtSqSum[1]*fr[10]-1.082531754730548*nuVtSqSum[1]*fl[10]+1.125*nuVtSqSum[3]*fr[9]-1.125*nuVtSqSum[3]*fl[9]+1.125*nuVtSqSum[0]*fr[8]-1.125*nuVtSqSum[0]*fl[8]+1.125*nuVtSqSum[1]*fr[4]-1.125*nuVtSqSum[1]*fl[4])*rdv+0.5*alphaDrSurf[2]*fhatAL[7]+0.5*alphaDrSurf[4]*fhatAL[6]+0.5*alphaDrSurf[0]*fhatAL[5]+0.5*alphaDrSurf[1]*fhatAL[3]; 
  Ghat[9] = ((-1.082531754730548*nuVtSqSum[1]*fr[15])-1.082531754730548*nuVtSqSum[1]*fl[15]-1.082531754730548*nuVtSqSum[0]*fr[14]-1.082531754730548*nuVtSqSum[0]*fl[14]-1.082531754730548*nuVtSqSum[3]*fr[13]-1.082531754730548*nuVtSqSum[3]*fl[13]+1.125*nuVtSqSum[1]*fr[12]-1.125*nuVtSqSum[1]*fl[12]-1.082531754730548*nuVtSqSum[2]*fr[10]-1.082531754730548*nuVtSqSum[2]*fl[10]+1.125*nuVtSqSum[0]*fr[9]-1.125*nuVtSqSum[0]*fl[9]+1.125*nuVtSqSum[3]*fr[8]-1.125*nuVtSqSum[3]*fl[8]+1.125*nuVtSqSum[2]*fr[4]-1.125*nuVtSqSum[2]*fl[4])*rdv+0.5*alphaDrSurf[1]*fhatAL[7]+0.5*alphaDrSurf[0]*fhatAL[6]+0.5*alphaDrSurf[4]*fhatAL[5]+0.5*alphaDrSurf[2]*fhatAL[3]; 
  Ghat[12] = ((-1.082531754730548*nuVtSqSum[0]*fr[15])-1.082531754730548*nuVtSqSum[0]*fl[15]-1.082531754730548*nuVtSqSum[1]*fr[14]-1.082531754730548*nuVtSqSum[1]*fl[14]-1.082531754730548*nuVtSqSum[2]*fr[13]-1.082531754730548*nuVtSqSum[2]*fl[13]+1.125*nuVtSqSum[0]*fr[12]-1.125*nuVtSqSum[0]*fl[12]-1.082531754730548*nuVtSqSum[3]*fr[10]-1.082531754730548*nuVtSqSum[3]*fl[10]+1.125*nuVtSqSum[1]*fr[9]-1.125*nuVtSqSum[1]*fl[9]+1.125*nuVtSqSum[2]*fr[8]-1.125*nuVtSqSum[2]*fl[8]+1.125*nuVtSqSum[3]*fr[4]-1.125*nuVtSqSum[3]*fl[4])*rdv+0.5*alphaDrSurf[0]*fhatAL[7]+0.5*alphaDrSurf[1]*fhatAL[6]+0.5*alphaDrSurf[2]*fhatAL[5]+0.5*fhatAL[3]*alphaDrSurf[4]; 

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

  double incr2[16]; 
  incr2[3] = nuVtSqSum[3]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[2]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[1]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[6] = nuVtSqSum[2]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[3]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[0]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[7] = nuVtSqSum[1]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[0]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[3]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+nuVtSqSum[2]*((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0])); 
  incr2[10] = nuVtSqSum[3]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[2]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[1]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[0]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[11] = nuVtSqSum[0]*((-0.25*fr[11])+0.25*fl[11]+0.2165063509461096*(fr[5]+fl[5]))+nuVtSqSum[1]*((-0.25*fr[7])+0.25*fl[7]+0.2165063509461096*(fr[2]+fl[2]))+nuVtSqSum[2]*((-0.25*fr[6])+0.25*fl[6]+0.2165063509461096*(fr[1]+fl[1]))+((-0.25*fr[3])+0.25*fl[3]+0.2165063509461096*(fr[0]+fl[0]))*nuVtSqSum[3]; 
  incr2[13] = nuVtSqSum[2]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[3]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[0]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[1]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[14] = nuVtSqSum[1]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[0]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[3]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[2]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 
  incr2[15] = nuVtSqSum[0]*((-0.25*fr[15])+0.25*fl[15]+0.2165063509461096*(fr[12]+fl[12]))+nuVtSqSum[1]*((-0.25*fr[14])+0.25*fl[14]+0.2165063509461096*(fr[9]+fl[9]))+nuVtSqSum[2]*((-0.25*fr[13])+0.25*fl[13]+0.2165063509461096*(fr[8]+fl[8]))+nuVtSqSum[3]*((-0.25*fr[10])+0.25*fl[10]+0.2165063509461096*(fr[4]+fl[4])); 

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
double GkLBOconstNuSurfPositivity2x2vSer_Mu_P1(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:        collisionalities added (self and cross species collisionalities). 
  // vMuMidMax:    maximum midpoint value of v-u. 
  // nuUSum[4]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fr:        Distribution function in left/right cells 
  // outl/outr:    Incremented distribution function in left/right cells 
  double rdv = 1.0/dxvl[3]; 
  double rdv2L = 2.0/dxvl[3]; 
  double rdv2R = 2.0/dxvr[3]; 
  double rdvSq4L = rdv2L*rdv2L; 
  double rdvSq4R = rdv2R*rdv2R; 

  double alphaDrSurf[8]; 
  alphaDrSurf[0] = 1.414213562373095*(4.0*wl[3]+2.0*dxvl[3])*nuSum; 

  double f0Quad[8]; 
  double f1Quad[8]; 
  double limQuad[8]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = -0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[0] = 0.1767766952966368*((fr[11]-1.0*(fl[11]+fr[7])+fl[7]-1.0*fr[6]+fl[6]-1.0*fr[5]+fl[5]+fr[3]-1.0*fl[3]+fr[2]-1.0*fl[2]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[11]+fl[11])+fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5]-1.0*(fr[3]+fl[3]+fr[2]+fl[2]+fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[0] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14]+fr[13]+fl[13]+fr[12]+fl[12])+fr[10]+fl[10]+fr[9]+fl[9]+fr[8]+fl[8]-1.0*(fr[4]+fl[4]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*fl[14]+fr[13]-1.0*fl[13]+fr[12]-1.0*(fl[12]+fr[10])+fl[10]-1.0*fr[9]+fl[9]-1.0*fr[8]+fl[8]+fr[4]-1.0*fl[4]); 
  limQuad[0] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = -0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[1] = -0.1767766952966368*((fr[11]-1.0*fl[11]+fr[7]-1.0*(fl[7]+fr[6])+fl[6]-1.0*fr[5]+fl[5]-1.0*fr[3]+fl[3]-1.0*fr[2]+fl[2]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[11]+fl[11]+fr[7]+fl[7])+fr[6]+fl[6]+fr[5]+fl[5]+fr[3]+fl[3]+fr[2]+fl[2]-1.0*(fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[1] = -0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]-1.0*(fr[13]+fl[13]+fr[12]+fl[12]+fr[10]+fl[10]+fr[9]+fl[9])+fr[8]+fl[8]+fr[4]+fl[4])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]+fr[13]-1.0*fl[13]+fr[12]-1.0*fl[12]+fr[10]-1.0*fl[10]+fr[9]-1.0*(fl[9]+fr[8])+fl[8]-1.0*fr[4]+fl[4]); 
  limQuad[1] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = -0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[2] = -0.1767766952966368*((fr[11]-1.0*(fl[11]+fr[7])+fl[7]+fr[6]-1.0*(fl[6]+fr[5])+fl[5]-1.0*fr[3]+fl[3]+fr[2]-1.0*(fl[2]+fr[1])+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[11]+fl[11])+fr[7]+fl[7]-1.0*(fr[6]+fl[6])+fr[5]+fl[5]+fr[3]+fl[3]-1.0*(fr[2]+fl[2])+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[2] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14])+fr[13]+fl[13]-1.0*(fr[12]+fl[12]+fr[10]+fl[10])+fr[9]+fl[9]-1.0*(fr[8]+fl[8])+fr[4]+fl[4])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*(fl[14]+fr[13])+fl[13]+fr[12]-1.0*fl[12]+fr[10]-1.0*(fl[10]+fr[9])+fl[9]+fr[8]-1.0*(fl[8]+fr[4])+fl[4]); 
  limQuad[2] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = -0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[3] = 0.1767766952966368*((fr[11]-1.0*fl[11]+fr[7]-1.0*fl[7]+fr[6]-1.0*(fl[6]+fr[5])+fl[5]+fr[3]-1.0*(fl[3]+fr[2])+fl[2]-1.0*fr[1]+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[11]+fl[11]+fr[7]+fl[7]+fr[6]+fl[6])+fr[5]+fl[5]-1.0*(fr[3]+fl[3])+fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[3] = 0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]-1.0*(fr[12]+fl[12])+fr[10]+fl[10]-1.0*(fr[9]+fl[9]+fr[8]+fl[8]+fr[4]+fl[4]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]-1.0*fr[13]+fl[13]+fr[12]-1.0*(fl[12]+fr[10])+fl[10]+fr[9]-1.0*fl[9]+fr[8]-1.0*fl[8]+fr[4]-1.0*fl[4]); 
  limQuad[3] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = -0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[4] = -0.1767766952966368*((fr[11]-1.0*(fl[11]+fr[7])+fl[7]-1.0*fr[6]+fl[6]+fr[5]-1.0*fl[5]+fr[3]-1.0*(fl[3]+fr[2])+fl[2]-1.0*fr[1]+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[11]+fl[11])+fr[7]+fl[7]+fr[6]+fl[6]-1.0*(fr[5]+fl[5]+fr[3]+fl[3])+fr[2]+fl[2]+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[4] = -0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14]+fr[13]+fl[13])+fr[12]+fl[12]+fr[10]+fl[10]-1.0*(fr[9]+fl[9]+fr[8]+fl[8])+fr[4]+fl[4])*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*fl[14]+fr[13]-1.0*(fl[13]+fr[12])+fl[12]-1.0*fr[10]+fl[10]+fr[9]-1.0*fl[9]+fr[8]-1.0*(fl[8]+fr[4])+fl[4]); 
  limQuad[4] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = -0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[5] = 0.1767766952966368*((fr[11]-1.0*fl[11]+fr[7]-1.0*(fl[7]+fr[6])+fl[6]+fr[5]-1.0*(fl[5]+fr[3])+fl[3]+fr[2]-1.0*(fl[2]+fr[1])+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[11]+fl[11]+fr[7]+fl[7])+fr[6]+fl[6]-1.0*(fr[5]+fl[5])+fr[3]+fl[3]-1.0*(fr[2]+fl[2])+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[5] = 0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]-1.0*(fr[13]+fl[13])+fr[12]+fl[12]-1.0*(fr[10]+fl[10])+fr[9]+fl[9]-1.0*(fr[8]+fl[8]+fr[4]+fl[4]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]+fr[13]-1.0*(fl[13]+fr[12])+fl[12]+fr[10]-1.0*(fl[10]+fr[9])+fl[9]+fr[8]-1.0*fl[8]+fr[4]-1.0*fl[4]); 
  limQuad[5] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = -0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[6] = 0.1767766952966368*((fr[11]-1.0*(fl[11]+fr[7])+fl[7]+fr[6]-1.0*fl[6]+fr[5]-1.0*(fl[5]+fr[3])+fl[3]-1.0*fr[2]+fl[2]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[11]+fl[11])+fr[7]+fl[7]-1.0*(fr[6]+fl[6]+fr[5]+fl[5])+fr[3]+fl[3]+fr[2]+fl[2]-1.0*(fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[6] = 0.1767766952966368*((fr[15]+fl[15]-1.0*(fr[14]+fl[14])+fr[13]+fl[13]+fr[12]+fl[12]-1.0*(fr[10]+fl[10]+fr[9]+fl[9])+fr[8]+fl[8]-1.0*(fr[4]+fl[4]))*sgn(alphaQuad)-1.0*fr[15]+fl[15]+fr[14]-1.0*(fl[14]+fr[13])+fl[13]-1.0*fr[12]+fl[12]+fr[10]-1.0*fl[10]+fr[9]-1.0*(fl[9]+fr[8])+fl[8]+fr[4]-1.0*fl[4]); 
  limQuad[6] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  alphaQuad = -0.3535533905932737*alphaDrSurf[0]; 
  f0Quad[7] = -0.1767766952966368*((fr[11]-1.0*fl[11]+fr[7]-1.0*fl[7]+fr[6]-1.0*fl[6]+fr[5]-1.0*fl[5]+fr[3]-1.0*fl[3]+fr[2]-1.0*fl[2]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[11]+fl[11]+fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5]+fr[3]+fl[3]+fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[7] = -0.1767766952966368*((fr[15]+fl[15]+fr[14]+fl[14]+fr[13]+fl[13]+fr[12]+fl[12]+fr[10]+fl[10]+fr[9]+fl[9]+fr[8]+fl[8]+fr[4]+fl[4])*sgn(alphaQuad)-1.0*fr[15]+fl[15]-1.0*fr[14]+fl[14]-1.0*fr[13]+fl[13]-1.0*fr[12]+fl[12]-1.0*fr[10]+fl[10]-1.0*fr[9]+fl[9]-1.0*fr[8]+fl[8]-1.0*fr[4]+fl[4]); 
  limQuad[7] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.25; 
  double fhat[16]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*f0Quad[4]+f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4])+f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[3] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]+f0Quad[5]+f0Quad[4]-1.0*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[5] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]+f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[6] = 0.3535533905932737*(f0Quad[7]-1.0*f0Quad[6]+f0Quad[5]-1.0*(f0Quad[4]+f0Quad[3])+f0Quad[2]-1.0*f0Quad[1]+f0Quad[0]); 
  fhat[7] = 0.3535533905932737*(f0Quad[7]+f0Quad[6]-1.0*(f0Quad[5]+f0Quad[4]+f0Quad[3]+f0Quad[2])+f0Quad[1]+f0Quad[0]); 
  fhat[8] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*f1Quad[4]+f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[9] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4])+f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[10] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]+f1Quad[5]+f1Quad[4]-1.0*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0])); 
  fhat[11] = 0.3535533905932737*(f0Quad[7]-1.0*(f0Quad[6]+f0Quad[5])+f0Quad[4]-1.0*f0Quad[3]+f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[12] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]+f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  fhat[13] = 0.3535533905932737*(f1Quad[7]-1.0*f1Quad[6]+f1Quad[5]-1.0*(f1Quad[4]+f1Quad[3])+f1Quad[2]-1.0*f1Quad[1]+f1Quad[0]); 
  fhat[14] = 0.3535533905932737*(f1Quad[7]+f1Quad[6]-1.0*(f1Quad[5]+f1Quad[4]+f1Quad[3]+f1Quad[2])+f1Quad[1]+f1Quad[0]); 
  fhat[15] = 0.3535533905932737*(f1Quad[7]-1.0*(f1Quad[6]+f1Quad[5])+f1Quad[4]-1.0*f1Quad[3]+f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[8];  // rCtrl=f1/f0 at each control node in dimensions other than vy 
  rCtrl[0] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[14]+fhat[13]+fhat[12]+3.0*fhat[4])+9.0*(fhat[10]+fhat[9]+fhat[8])))/(36.0*EPSILON-1.732050807568877*(fhat[11]+3.0*(fhat[3]+fhat[2]+fhat[1]))+3.0*(fhat[7]+fhat[6]+fhat[5]+3.0*fhat[0])); 
  rCtrl[1] = (3.0*fhat[15]+5.196152422706631*(fhat[14]-1.0*(fhat[13]+fhat[12])+3.0*fhat[4])+9.0*(fhat[8]-1.0*(fhat[10]+fhat[9])))/(36.0*EPSILON+1.732050807568877*(fhat[11]+3.0*(fhat[1]-1.0*(fhat[3]+fhat[2])))+3.0*(fhat[7]-1.0*(fhat[6]+fhat[5]-3.0*fhat[0]))); 
  rCtrl[2] = (3.0*fhat[15]-5.196152422706631*(fhat[14]-1.0*fhat[13]+fhat[12]-3.0*fhat[4])+9.0*((-1.0*fhat[10])+fhat[9]-1.0*fhat[8]))/(36.0*EPSILON+1.732050807568877*(fhat[11]+3.0*((-1.0*fhat[3])+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*fhat[7])+fhat[6]-1.0*fhat[5]+3.0*fhat[0])); 
  rCtrl[3] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[14]+fhat[13]-1.0*(fhat[12]+3.0*fhat[4]))+9.0*(fhat[10]-1.0*(fhat[9]+fhat[8]))))/(36.0*EPSILON-1.732050807568877*(fhat[11]+3.0*(fhat[3]-1.0*(fhat[2]+fhat[1])))+3.0*((-1.0*(fhat[7]+fhat[6]))+fhat[5]+3.0*fhat[0])); 
  rCtrl[4] = (3.0*fhat[15]-5.196152422706631*(fhat[14]+fhat[13]-1.0*(fhat[12]+3.0*fhat[4]))+9.0*(fhat[10]-1.0*(fhat[9]+fhat[8])))/(36.0*EPSILON+1.732050807568877*(fhat[11]+3.0*(fhat[3]-1.0*(fhat[2]+fhat[1])))+3.0*((-1.0*(fhat[7]+fhat[6]))+fhat[5]+3.0*fhat[0])); 
  rCtrl[5] = -(1.0*(3.0*fhat[15]+5.196152422706631*(fhat[14]-1.0*fhat[13]+fhat[12]-3.0*fhat[4])+9.0*((-1.0*fhat[10])+fhat[9]-1.0*fhat[8])))/(36.0*EPSILON-1.732050807568877*(fhat[11]+3.0*((-1.0*fhat[3])+fhat[2]-1.0*fhat[1]))+3.0*((-1.0*fhat[7])+fhat[6]-1.0*fhat[5]+3.0*fhat[0])); 
  rCtrl[6] = -(1.0*(3.0*fhat[15]-5.196152422706631*(fhat[14]-1.0*(fhat[13]+fhat[12])+3.0*fhat[4])+9.0*(fhat[8]-1.0*(fhat[10]+fhat[9]))))/(36.0*EPSILON-1.732050807568877*(fhat[11]+3.0*(fhat[1]-1.0*(fhat[3]+fhat[2])))+3.0*(fhat[7]-1.0*(fhat[6]+fhat[5]-3.0*fhat[0]))); 
  rCtrl[7] = (3.0*fhat[15]+5.196152422706631*(fhat[14]+fhat[13]+fhat[12]+3.0*fhat[4])+9.0*(fhat[10]+fhat[9]+fhat[8]))/(36.0*EPSILON+1.732050807568877*(fhat[11]+3.0*(fhat[3]+fhat[2]+fhat[1]))+3.0*(fhat[7]+fhat[6]+fhat[5]+3.0*fhat[0])); 
  double fhatCtrl[8];  // fhatCtrl = anti-limited fhat evaluated at each control node on vy surface 
  fhatCtrl[0] = -0.02777777777777778*(1.732050807568877*(fhat[11]+3.0*(fhat[3]+fhat[2]+fhat[1]))-3.0*(fhat[7]+fhat[6]+fhat[5]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.02777777777777778*(1.732050807568877*(fhat[11]+3.0*(fhat[1]-1.0*(fhat[3]+fhat[2])))+3.0*(fhat[7]-1.0*(fhat[6]+fhat[5])+3.0*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.02777777777777778*(1.732050807568877*(fhat[11]+3.0*((-1.0*fhat[3])+fhat[2]-1.0*fhat[1]))-3.0*(fhat[7]-1.0*fhat[6]+fhat[5]-3.0*fhat[0]))*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = -0.02777777777777778*(1.732050807568877*(fhat[11]+3.0*(fhat[3]-1.0*(fhat[2]+fhat[1])))+3.0*(fhat[7]+fhat[6]-1.0*(fhat[5]+3.0*fhat[0])))*limTheta(rCtrl[3],-1.0); 
  fhatCtrl[4] = 0.02777777777777778*(1.732050807568877*(fhat[11]+3.0*(fhat[3]-1.0*(fhat[2]+fhat[1])))-3.0*(fhat[7]+fhat[6]-1.0*(fhat[5]+3.0*fhat[0])))*limTheta(rCtrl[4],-1.0); 
  fhatCtrl[5] = -0.02777777777777778*(1.732050807568877*(fhat[11]+3.0*((-1.0*fhat[3])+fhat[2]-1.0*fhat[1]))+3.0*(fhat[7]-1.0*fhat[6]+fhat[5]-3.0*fhat[0]))*limTheta(rCtrl[5],-1.0); 
  fhatCtrl[6] = -0.02777777777777778*(1.732050807568877*(fhat[11]+3.0*(fhat[1]-1.0*(fhat[3]+fhat[2])))-3.0*(fhat[7]-1.0*(fhat[6]+fhat[5])+3.0*fhat[0]))*limTheta(rCtrl[6],-1.0); 
  fhatCtrl[7] = 0.02777777777777778*(1.732050807568877*(fhat[11]+3.0*(fhat[3]+fhat[2]+fhat[1]))+3.0*(fhat[7]+fhat[6]+fhat[5]+3.0*fhat[0]))*limTheta(rCtrl[7],-1.0); 
  double fhatAL[8];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.3535533905932737*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.6123724356957944*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*fhatCtrl[4]+fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4])+fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 0.6123724356957944*(fhatCtrl[7]+fhatCtrl[6]+fhatCtrl[5]+fhatCtrl[4]-1.0*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[4] = 1.060660171779821*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]+fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  fhatAL[5] = 1.060660171779821*(fhatCtrl[7]-1.0*fhatCtrl[6]+fhatCtrl[5]-1.0*(fhatCtrl[4]+fhatCtrl[3])+fhatCtrl[2]-1.0*fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[6] = 1.060660171779821*(fhatCtrl[7]+fhatCtrl[6]-1.0*(fhatCtrl[5]+fhatCtrl[4]+fhatCtrl[3]+fhatCtrl[2])+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[7] = 1.837117307087383*(fhatCtrl[7]-1.0*(fhatCtrl[6]+fhatCtrl[5])+fhatCtrl[4]-1.0*fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[8]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.25*(1.414213562373095*(1.0*(fhatAL[6]+fhatAL[4])+fhatAL[5]+fhatAL[0])-1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*fhatAL[1]-0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[2]))+1.414213562373095*((-0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[6]))-1.0*fhatAL[5]+fhatAL[0])), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[7]+fhatAL[2])-0.5773502691896258*(fhatAL[3]+fhatAL[1]))+1.414213562373095*((-1.0*(fhatAL[6]+fhatAL[4]))+fhatAL[5]+fhatAL[0])), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[2]-1.0*fhatAL[7])-0.5773502691896258*fhatAL[3]+0.5773502691896258*fhatAL[1])+1.414213562373095*(0.5773502691896258*(1.732050807568877*fhatAL[4]-1.732050807568877*fhatAL[6])-1.0*fhatAL[5]+fhatAL[0])), limQuad[3])); 
  fhatALQuad[4] = std::max(0.0, std::min(0.25*(2.449489742783178*((-0.5773502691896258*(fhatAL[2]-1.0*fhatAL[7]))+0.5773502691896258*fhatAL[3]-0.5773502691896258*fhatAL[1])+1.414213562373095*((-0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[4]))-1.0*fhatAL[5]+fhatAL[0])), limQuad[4])); 
  fhatALQuad[5] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*(fhatAL[3]+fhatAL[1])-0.5773502691896258*(fhatAL[7]+fhatAL[2]))+1.414213562373095*((-1.0*(fhatAL[6]+fhatAL[4]))+fhatAL[5]+fhatAL[0])), limQuad[5])); 
  fhatALQuad[6] = std::max(0.0, std::min(0.25*(2.449489742783178*(0.5773502691896258*((-1.0*fhatAL[7])+fhatAL[3]+fhatAL[2])-0.5773502691896258*fhatAL[1])+1.414213562373095*(0.5773502691896258*(1.732050807568877*fhatAL[6]-1.732050807568877*fhatAL[4])-1.0*fhatAL[5]+fhatAL[0])), limQuad[6])); 
  fhatALQuad[7] = std::max(0.0, std::min(0.25*(1.414213562373095*(fhatAL[7]+fhatAL[3]+fhatAL[2]+fhatAL[1])+1.414213562373095*(1.0*(fhatAL[6]+fhatAL[4])+fhatAL[5]+fhatAL[0])), limQuad[7])); 
  fhatAL[0] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*fhatALQuad[4]+fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4])+fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]+fhatALQuad[5]+fhatALQuad[4]-1.0*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[4] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]+fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 
  fhatAL[5] = 0.3535533905932737*(fhatALQuad[7]-1.0*fhatALQuad[6]+fhatALQuad[5]-1.0*(fhatALQuad[4]+fhatALQuad[3])+fhatALQuad[2]-1.0*fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[6] = 0.3535533905932737*(fhatALQuad[7]+fhatALQuad[6]-1.0*(fhatALQuad[5]+fhatALQuad[4]+fhatALQuad[3]+fhatALQuad[2])+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[7] = 0.3535533905932737*(fhatALQuad[7]-1.0*(fhatALQuad[6]+fhatALQuad[5])+fhatALQuad[4]-1.0*fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double diffFac[4]; 
  diffFac[0] = (BmagInv[3]*nuVtSqSum[3]+BmagInv[2]*nuVtSqSum[2]+BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*(wl[3]+0.5*dxvl[3])*m_; 
  diffFac[1] = (BmagInv[2]*nuVtSqSum[3]+nuVtSqSum[2]*BmagInv[3]+BmagInv[0]*nuVtSqSum[1]+nuVtSqSum[0]*BmagInv[1])*(wl[3]+0.5*dxvl[3])*m_; 
  diffFac[2] = (BmagInv[1]*nuVtSqSum[3]+nuVtSqSum[1]*BmagInv[3]+BmagInv[0]*nuVtSqSum[2]+nuVtSqSum[0]*BmagInv[2])*(wl[3]+0.5*dxvl[3])*m_; 
  diffFac[3] = (BmagInv[0]*nuVtSqSum[3]+nuVtSqSum[0]*BmagInv[3]+BmagInv[1]*nuVtSqSum[2]+nuVtSqSum[1]*BmagInv[2])*(wl[3]+0.5*dxvl[3])*m_; 

  double Ghat[16]; 
  for(unsigned short int i=0; i<16; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.082531754730548*diffFac[3]*fr[12])-1.082531754730548*diffFac[3]*fl[12]-1.082531754730548*diffFac[2]*fr[9]-1.082531754730548*diffFac[2]*fl[9]-1.082531754730548*diffFac[1]*fr[8]-1.082531754730548*diffFac[1]*fl[8]+1.125*diffFac[3]*fr[5]-1.125*diffFac[3]*fl[5]-1.082531754730548*diffFac[0]*fr[4]-1.082531754730548*diffFac[0]*fl[4]+1.125*diffFac[2]*fr[2]-1.125*diffFac[2]*fl[2]+1.125*diffFac[1]*fr[1]-1.125*diffFac[1]*fl[1]+1.125*diffFac[0]*fr[0]-1.125*diffFac[0]*fl[0])*rdv+0.5*alphaDrSurf[0]*fhatAL[0]; 
  Ghat[1] = ((-1.082531754730548*diffFac[2]*fr[12])-1.082531754730548*diffFac[2]*fl[12]-1.082531754730548*diffFac[3]*fr[9]-1.082531754730548*diffFac[3]*fl[9]-1.082531754730548*diffFac[0]*fr[8]-1.082531754730548*diffFac[0]*fl[8]+1.125*diffFac[2]*fr[5]-1.125*diffFac[2]*fl[5]-1.082531754730548*diffFac[1]*fr[4]-1.082531754730548*diffFac[1]*fl[4]+1.125*fr[2]*diffFac[3]-1.125*fl[2]*diffFac[3]+1.125*diffFac[0]*fr[1]-1.125*diffFac[0]*fl[1]+1.125*fr[0]*diffFac[1]-1.125*fl[0]*diffFac[1])*rdv+0.5*alphaDrSurf[0]*fhatAL[1]; 
  Ghat[2] = ((-1.082531754730548*diffFac[1]*fr[12])-1.082531754730548*diffFac[1]*fl[12]-1.082531754730548*diffFac[0]*fr[9]-1.082531754730548*diffFac[0]*fl[9]-1.082531754730548*diffFac[3]*fr[8]-1.082531754730548*diffFac[3]*fl[8]+1.125*diffFac[1]*fr[5]-1.125*diffFac[1]*fl[5]-1.082531754730548*diffFac[2]*fr[4]-1.082531754730548*diffFac[2]*fl[4]+1.125*fr[1]*diffFac[3]-1.125*fl[1]*diffFac[3]+1.125*diffFac[0]*fr[2]-1.125*diffFac[0]*fl[2]+1.125*fr[0]*diffFac[2]-1.125*fl[0]*diffFac[2])*rdv+0.5*alphaDrSurf[0]*fhatAL[2]; 
  Ghat[3] = ((-1.082531754730548*diffFac[3]*fr[15])-1.082531754730548*diffFac[3]*fl[15]-1.082531754730548*diffFac[2]*fr[14]-1.082531754730548*diffFac[2]*fl[14]-1.082531754730548*diffFac[1]*fr[13]-1.082531754730548*diffFac[1]*fl[13]+1.125*diffFac[3]*fr[11]-1.125*diffFac[3]*fl[11]-1.082531754730548*diffFac[0]*fr[10]-1.082531754730548*diffFac[0]*fl[10]+1.125*diffFac[2]*fr[7]-1.125*diffFac[2]*fl[7]+1.125*diffFac[1]*fr[6]-1.125*diffFac[1]*fl[6]+1.125*diffFac[0]*fr[3]-1.125*diffFac[0]*fl[3])*rdv+0.5*alphaDrSurf[0]*fhatAL[3]; 
  Ghat[5] = ((-1.082531754730548*diffFac[0]*fr[12])-1.082531754730548*diffFac[0]*fl[12]-1.082531754730548*diffFac[1]*fr[9]-1.082531754730548*diffFac[1]*fl[9]-1.082531754730548*diffFac[2]*fr[8]-1.082531754730548*diffFac[2]*fl[8]+1.125*diffFac[0]*fr[5]-1.125*diffFac[0]*fl[5]-1.082531754730548*diffFac[3]*fr[4]-1.082531754730548*diffFac[3]*fl[4]+1.125*fr[0]*diffFac[3]-1.125*fl[0]*diffFac[3]+1.125*diffFac[1]*fr[2]-1.125*diffFac[1]*fl[2]+1.125*fr[1]*diffFac[2]-1.125*fl[1]*diffFac[2])*rdv+0.5*alphaDrSurf[0]*fhatAL[4]; 
  Ghat[6] = ((-1.082531754730548*diffFac[2]*fr[15])-1.082531754730548*diffFac[2]*fl[15]-1.082531754730548*diffFac[3]*fr[14]-1.082531754730548*diffFac[3]*fl[14]-1.082531754730548*diffFac[0]*fr[13]-1.082531754730548*diffFac[0]*fl[13]+1.125*diffFac[2]*fr[11]-1.125*diffFac[2]*fl[11]-1.082531754730548*diffFac[1]*fr[10]-1.082531754730548*diffFac[1]*fl[10]+1.125*diffFac[3]*fr[7]-1.125*diffFac[3]*fl[7]+1.125*diffFac[0]*fr[6]-1.125*diffFac[0]*fl[6]+1.125*diffFac[1]*fr[3]-1.125*diffFac[1]*fl[3])*rdv+0.5*alphaDrSurf[0]*fhatAL[5]; 
  Ghat[7] = ((-1.082531754730548*diffFac[1]*fr[15])-1.082531754730548*diffFac[1]*fl[15]-1.082531754730548*diffFac[0]*fr[14]-1.082531754730548*diffFac[0]*fl[14]-1.082531754730548*diffFac[3]*fr[13]-1.082531754730548*diffFac[3]*fl[13]+1.125*diffFac[1]*fr[11]-1.125*diffFac[1]*fl[11]-1.082531754730548*diffFac[2]*fr[10]-1.082531754730548*diffFac[2]*fl[10]+1.125*diffFac[0]*fr[7]-1.125*diffFac[0]*fl[7]+1.125*diffFac[3]*fr[6]-1.125*diffFac[3]*fl[6]+1.125*diffFac[2]*fr[3]-1.125*diffFac[2]*fl[3])*rdv+0.5*alphaDrSurf[0]*fhatAL[6]; 
  Ghat[11] = ((-1.082531754730548*diffFac[0]*fr[15])-1.082531754730548*diffFac[0]*fl[15]-1.082531754730548*diffFac[1]*fr[14]-1.082531754730548*diffFac[1]*fl[14]-1.082531754730548*diffFac[2]*fr[13]-1.082531754730548*diffFac[2]*fl[13]+1.125*diffFac[0]*fr[11]-1.125*diffFac[0]*fl[11]-1.082531754730548*diffFac[3]*fr[10]-1.082531754730548*diffFac[3]*fl[10]+1.125*diffFac[1]*fr[7]-1.125*diffFac[1]*fl[7]+1.125*diffFac[2]*fr[6]-1.125*diffFac[2]*fl[6]+1.125*diffFac[3]*fr[3]-1.125*diffFac[3]*fl[3])*rdv+0.5*alphaDrSurf[0]*fhatAL[7]; 

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

  double incr2[16]; 
  incr2[4] = diffFac[3]*(0.25*fl[12]-0.25*fr[12])+diffFac[2]*(0.25*fl[9]-0.25*fr[9])+diffFac[1]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[3]*(fr[5]+fl[5])+diffFac[0]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*(diffFac[2]*(fr[2]+fl[2])+diffFac[1]*(fr[1]+fl[1])+diffFac[0]*(fr[0]+fl[0])); 
  incr2[8] = diffFac[2]*(0.25*fl[12]-0.25*fr[12])+diffFac[3]*(0.25*fl[9]-0.25*fr[9])+diffFac[0]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[2]*(fr[5]+fl[5])+diffFac[1]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*((fr[2]+fl[2])*diffFac[3]+diffFac[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*diffFac[1]); 
  incr2[9] = diffFac[1]*(0.25*fl[12]-0.25*fr[12])+diffFac[0]*(0.25*fl[9]-0.25*fr[9])+diffFac[3]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[1]*(fr[5]+fl[5])+diffFac[2]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*((fr[1]+fl[1])*diffFac[3]+diffFac[0]*(fr[2]+fl[2])+(fr[0]+fl[0])*diffFac[2]); 
  incr2[10] = diffFac[3]*(0.25*fl[15]-0.25*fr[15])+diffFac[2]*(0.25*fl[14]-0.25*fr[14])+diffFac[1]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[3]*(fr[11]+fl[11])+diffFac[0]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[2]*(fr[7]+fl[7])+diffFac[1]*(fr[6]+fl[6])+diffFac[0]*(fr[3]+fl[3])); 
  incr2[12] = diffFac[0]*(0.25*fl[12]-0.25*fr[12])+diffFac[1]*(0.25*fl[9]-0.25*fr[9])+diffFac[2]*(0.25*fl[8]-0.25*fr[8])+0.2165063509461096*diffFac[0]*(fr[5]+fl[5])+diffFac[3]*(0.25*fl[4]-0.25*fr[4])+0.2165063509461096*((fr[0]+fl[0])*diffFac[3]+diffFac[1]*(fr[2]+fl[2])+(fr[1]+fl[1])*diffFac[2]); 
  incr2[13] = diffFac[2]*(0.25*fl[15]-0.25*fr[15])+diffFac[3]*(0.25*fl[14]-0.25*fr[14])+diffFac[0]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[2]*(fr[11]+fl[11])+diffFac[1]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[3]*(fr[7]+fl[7])+diffFac[0]*(fr[6]+fl[6])+diffFac[1]*(fr[3]+fl[3])); 
  incr2[14] = diffFac[1]*(0.25*fl[15]-0.25*fr[15])+diffFac[0]*(0.25*fl[14]-0.25*fr[14])+diffFac[3]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[1]*(fr[11]+fl[11])+diffFac[2]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[0]*(fr[7]+fl[7])+diffFac[3]*(fr[6]+fl[6])+diffFac[2]*(fr[3]+fl[3])); 
  incr2[15] = diffFac[0]*(0.25*fl[15]-0.25*fr[15])+diffFac[1]*(0.25*fl[14]-0.25*fr[14])+diffFac[2]*(0.25*fl[13]-0.25*fr[13])+0.2165063509461096*diffFac[0]*(fr[11]+fl[11])+diffFac[3]*(0.25*fl[10]-0.25*fr[10])+0.2165063509461096*(diffFac[1]*(fr[7]+fl[7])+diffFac[2]*(fr[6]+fl[6])+diffFac[3]*(fr[3]+fl[3])); 

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
