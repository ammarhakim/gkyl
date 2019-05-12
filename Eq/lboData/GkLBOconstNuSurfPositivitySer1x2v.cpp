#include <GkLBOModDecl.h> 
double GkLBOconstNuSurfPositivity1x2vSer_Vpar_P1(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
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

  double f0Quad[4]; 
  double f1Quad[4]; 
  double limQuad[4]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.5*alphaDrSurf[1]-0.5*alphaDrSurf[0]; 
  f0Quad[0] = -0.25*((fr[5]-1.0*(fl[5]+fr[3])+fl[3]-1.0*fr[1]+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5])+fr[3]+fl[3]+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[0] = -0.25*((fr[7]+fl[7]-1.0*(fr[6]+fl[6]+fr[4]+fl[4])+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[6]-1.0*fl[6]+fr[4]-1.0*(fl[4]+fr[2])+fl[2]); 
  limQuad[0] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.3535533905932737; 
  alphaQuad = -0.5*(alphaDrSurf[1]+alphaDrSurf[0]); 
  f0Quad[1] = 0.25*((fr[5]-1.0*fl[5]+fr[3]-1.0*(fl[3]+fr[1])+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5]+fr[3]+fl[3])+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[1] = 0.25*((fr[7]+fl[7]+fr[6]+fl[6]-1.0*(fr[4]+fl[4]+fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]+fr[4]-1.0*fl[4]+fr[2]-1.0*fl[2]); 
  limQuad[1] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.3535533905932737; 
  alphaQuad = 0.5*alphaDrSurf[1]-0.5*alphaDrSurf[0]; 
  f0Quad[2] = 0.25*((fr[5]-1.0*(fl[5]+fr[3])+fl[3]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5])+fr[3]+fl[3]-1.0*(fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[2] = 0.25*((fr[7]+fl[7]-1.0*(fr[6]+fl[6])+fr[4]+fl[4]-1.0*(fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[6]-1.0*(fl[6]+fr[4])+fl[4]+fr[2]-1.0*fl[2]); 
  limQuad[2] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.3535533905932737; 
  alphaQuad = -0.5*(alphaDrSurf[1]+alphaDrSurf[0]); 
  f0Quad[3] = -0.25*((fr[5]-1.0*fl[5]+fr[3]-1.0*fl[3]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5]+fr[3]+fl[3]+fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[3] = -0.25*((fr[7]+fl[7]+fr[6]+fl[6]+fr[4]+fl[4]+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]-1.0*fr[4]+fl[4]-1.0*fr[2]+fl[2]); 
  limQuad[3] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.3535533905932737; 
  double fhat[8]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.5*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.5*(f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.5*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[3] = 0.5*(f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.5*(f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[5] = 0.5*(f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[6] = 0.5*(f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[7] = 0.5*(f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[4];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrl[0] = (1.414213562373095*(3.0*fhat[7]-5.196152422706631*(fhat[6]+fhat[4])+9.0*fhat[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fhat[5]-3.0*(fhat[3]+fhat[1])+5.196152422706631*fhat[0])); 
  rCtrl[1] = -(1.414213562373095*(3.0*fhat[7]+5.196152422706631*fhat[6]-1.0*(5.196152422706631*fhat[4]+9.0*fhat[2])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fhat[5])+3.0*(fhat[1]-1.0*fhat[3])+5.196152422706631*fhat[0])); 
  rCtrl[2] = -(1.414213562373095*(3.0*fhat[7]+5.196152422706631*(fhat[4]-1.0*fhat[6])-9.0*fhat[2]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fhat[5])+3.0*(fhat[3]-1.0*fhat[1])+5.196152422706631*fhat[0])); 
  rCtrl[3] = (1.414213562373095*(3.0*fhat[7]+5.196152422706631*(fhat[6]+fhat[4])+9.0*fhat[2]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fhat[5]+3.0*(fhat[3]+fhat[1])+5.196152422706631*fhat[0])); 
  double fhatCtrl[4];  // fhatCtrl = anti-limited fhat evaluated at each control node on vx surface 
  fhatCtrl[0] = 0.06804138174397717*(1.732050807568877*fhat[5]-3.0*(fhat[3]+fhat[1])+5.196152422706631*fhat[0])*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = -0.06804138174397717*(1.732050807568877*fhat[5]+3.0*fhat[3]-1.0*(3.0*fhat[1]+5.196152422706631*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = -0.06804138174397717*(1.732050807568877*fhat[5]+3.0*(fhat[1]-1.0*fhat[3])-5.196152422706631*fhat[0])*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = 0.06804138174397717*(1.732050807568877*fhat[5]+3.0*(fhat[3]+fhat[1])+5.196152422706631*fhat[0])*limTheta(rCtrl[3],-1.0); 
  double fhatAL[4];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.5*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.8660254037844386*(fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.8660254037844386*(fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 1.5*(fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[4]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.5*((-0.5773502691896258*(1.732050807568877*fhatAL[1]-1.732050807568877*fhatAL[3]))-1.0*fhatAL[2]+fhatAL[0]), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.5*(0.5773502691896258*(1.732050807568877*fhatAL[1]-1.732050807568877*fhatAL[3])-1.0*fhatAL[2]+fhatAL[0]), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.5*((-1.0*(fhatAL[3]+fhatAL[1]))+fhatAL[2]+fhatAL[0]), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.5*(1.0*(fhatAL[3]+fhatAL[1])+fhatAL[2]+fhatAL[0]), limQuad[3])); 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double Ghat[8]; 
  for(unsigned short int i=0; i<8; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.530931089239486*nuVtSqSum[1]*fr[4])-1.530931089239486*nuVtSqSum[1]*fl[4]-1.530931089239486*nuVtSqSum[0]*fr[2]-1.530931089239486*nuVtSqSum[0]*fl[2]+1.590990257669731*fr[1]*nuVtSqSum[1]-1.590990257669731*fl[1]*nuVtSqSum[1]+1.590990257669731*fr[0]*nuVtSqSum[0]-1.590990257669731*fl[0]*nuVtSqSum[0])*rdv+0.7071067811865475*alphaDrSurf[1]*fhatAL[1]+0.7071067811865475*alphaDrSurf[0]*fhatAL[0]; 
  Ghat[1] = ((-1.530931089239486*nuVtSqSum[0]*fr[4])-1.530931089239486*nuVtSqSum[0]*fl[4]-1.530931089239486*nuVtSqSum[1]*fr[2]-1.530931089239486*nuVtSqSum[1]*fl[2]+1.590990257669731*fr[0]*nuVtSqSum[1]-1.590990257669731*fl[0]*nuVtSqSum[1]+1.590990257669731*nuVtSqSum[0]*fr[1]-1.590990257669731*nuVtSqSum[0]*fl[1])*rdv+0.7071067811865475*alphaDrSurf[0]*fhatAL[1]+0.7071067811865475*fhatAL[0]*alphaDrSurf[1]; 
  Ghat[3] = ((-1.530931089239486*nuVtSqSum[1]*fr[7])-1.530931089239486*nuVtSqSum[1]*fl[7]-1.530931089239486*nuVtSqSum[0]*fr[6]-1.530931089239486*nuVtSqSum[0]*fl[6]+1.590990257669731*nuVtSqSum[1]*fr[5]-1.590990257669731*nuVtSqSum[1]*fl[5]+1.590990257669731*nuVtSqSum[0]*fr[3]-1.590990257669731*nuVtSqSum[0]*fl[3])*rdv+0.7071067811865475*alphaDrSurf[1]*fhatAL[3]+0.7071067811865475*alphaDrSurf[0]*fhatAL[2]; 
  Ghat[5] = ((-1.530931089239486*nuVtSqSum[0]*fr[7])-1.530931089239486*nuVtSqSum[0]*fl[7]-1.530931089239486*nuVtSqSum[1]*fr[6]-1.530931089239486*nuVtSqSum[1]*fl[6]+1.590990257669731*nuVtSqSum[0]*fr[5]-1.590990257669731*nuVtSqSum[0]*fl[5]+1.590990257669731*nuVtSqSum[1]*fr[3]-1.590990257669731*nuVtSqSum[1]*fl[3])*rdv+0.7071067811865475*alphaDrSurf[0]*fhatAL[3]+0.7071067811865475*alphaDrSurf[1]*fhatAL[2]; 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = -0.5*Ghat[3]; 
  incr1[4] = 0.8660254037844386*Ghat[1]; 
  incr1[5] = -0.5*Ghat[5]; 
  incr1[6] = 0.8660254037844386*Ghat[3]; 
  incr1[7] = 0.8660254037844386*Ghat[5]; 

  double incr2[8]; 
  incr2[2] = nuVtSqSum[1]*((-0.3535533905932737*fr[4])+0.3535533905932737*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[4] = nuVtSqSum[0]*((-0.3535533905932737*fr[4])+0.3535533905932737*fl[4]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[6] = nuVtSqSum[1]*((-0.3535533905932737*fr[7])+0.3535533905932737*fl[7]+0.3061862178478971*(fr[5]+fl[5]))+nuVtSqSum[0]*((-0.3535533905932737*fr[6])+0.3535533905932737*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 
  incr2[7] = nuVtSqSum[0]*((-0.3535533905932737*fr[7])+0.3535533905932737*fl[7]+0.3061862178478971*(fr[5]+fl[5]))+nuVtSqSum[1]*((-0.3535533905932737*fr[6])+0.3535533905932737*fl[6]+0.3061862178478971*(fr[3]+fl[3])); 

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

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
double GkLBOconstNuSurfPositivity1x2vSer_Mu_P1(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
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

  double alphaDrSurf[4]; 
  alphaDrSurf[0] = (4.0*wl[2]+2.0*dxvl[2])*nuSum; 

  double f0Quad[4]; 
  double f1Quad[4]; 
  double limQuad[4]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  f0Quad[0] = -0.25*((fr[4]-1.0*(fl[4]+fr[2])+fl[2]-1.0*fr[1]+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[4]+fl[4])+fr[2]+fl[2]+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[0] = -0.25*((fr[7]+fl[7]-1.0*(fr[6]+fl[6]+fr[5]+fl[5])+fr[3]+fl[3])*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[6]-1.0*fl[6]+fr[5]-1.0*(fl[5]+fr[3])+fl[3]); 
  limQuad[0] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.3535533905932737; 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  f0Quad[1] = 0.25*((fr[4]-1.0*fl[4]+fr[2]-1.0*(fl[2]+fr[1])+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[4]+fl[4]+fr[2]+fl[2])+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[1] = 0.25*((fr[7]+fl[7]+fr[6]+fl[6]-1.0*(fr[5]+fl[5]+fr[3]+fl[3]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]+fr[5]-1.0*fl[5]+fr[3]-1.0*fl[3]); 
  limQuad[1] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.3535533905932737; 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  f0Quad[2] = 0.25*((fr[4]-1.0*(fl[4]+fr[2])+fl[2]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[4]+fl[4])+fr[2]+fl[2]-1.0*(fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[2] = 0.25*((fr[7]+fl[7]-1.0*(fr[6]+fl[6])+fr[5]+fl[5]-1.0*(fr[3]+fl[3]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[6]-1.0*(fl[6]+fr[5])+fl[5]+fr[3]-1.0*fl[3]); 
  limQuad[2] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.3535533905932737; 
  alphaQuad = -0.5*alphaDrSurf[0]; 
  f0Quad[3] = -0.25*((fr[4]-1.0*fl[4]+fr[2]-1.0*fl[2]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[4]+fl[4]+fr[2]+fl[2]+fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[3] = -0.25*((fr[7]+fl[7]+fr[6]+fl[6]+fr[5]+fl[5]+fr[3]+fl[3])*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]-1.0*fr[5]+fl[5]-1.0*fr[3]+fl[3]); 
  limQuad[3] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.3535533905932737; 
  double fhat[8]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.5*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.5*(f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.5*(f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[3] = 0.5*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[4] = 0.5*(f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[5] = 0.5*(f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[6] = 0.5*(f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[7] = 0.5*(f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[4];  // rCtrl=f1/f0 at each control node in dimensions other than vy 
  rCtrl[0] = (1.414213562373095*(3.0*fhat[7]-5.196152422706631*(fhat[6]+fhat[5])+9.0*fhat[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fhat[4]-3.0*(fhat[2]+fhat[1])+5.196152422706631*fhat[0])); 
  rCtrl[1] = -(1.414213562373095*(3.0*fhat[7]+5.196152422706631*fhat[6]-1.0*(5.196152422706631*fhat[5]+9.0*fhat[3])))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fhat[4])+3.0*(fhat[1]-1.0*fhat[2])+5.196152422706631*fhat[0])); 
  rCtrl[2] = -(1.414213562373095*(3.0*fhat[7]+5.196152422706631*(fhat[5]-1.0*fhat[6])-9.0*fhat[3]))/(20.78460969082652*EPSILON+1.414213562373095*((-1.732050807568877*fhat[4])+3.0*(fhat[2]-1.0*fhat[1])+5.196152422706631*fhat[0])); 
  rCtrl[3] = (1.414213562373095*(3.0*fhat[7]+5.196152422706631*(fhat[6]+fhat[5])+9.0*fhat[3]))/(20.78460969082652*EPSILON+1.414213562373095*(1.732050807568877*fhat[4]+3.0*(fhat[2]+fhat[1])+5.196152422706631*fhat[0])); 
  double fhatCtrl[4];  // fhatCtrl = anti-limited fhat evaluated at each control node on vy surface 
  fhatCtrl[0] = 0.06804138174397717*(1.732050807568877*fhat[4]-3.0*(fhat[2]+fhat[1])+5.196152422706631*fhat[0])*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = -0.06804138174397717*(1.732050807568877*fhat[4]+3.0*fhat[2]-1.0*(3.0*fhat[1]+5.196152422706631*fhat[0]))*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = -0.06804138174397717*(1.732050807568877*fhat[4]+3.0*(fhat[1]-1.0*fhat[2])-5.196152422706631*fhat[0])*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = 0.06804138174397717*(1.732050807568877*fhat[4]+3.0*(fhat[2]+fhat[1])+5.196152422706631*fhat[0])*limTheta(rCtrl[3],-1.0); 
  double fhatAL[4];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.5*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.8660254037844386*(fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.8660254037844386*(fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 1.5*(fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[4]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.5*((-0.5773502691896258*(1.732050807568877*fhatAL[1]-1.732050807568877*fhatAL[3]))-1.0*fhatAL[2]+fhatAL[0]), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.5*(0.5773502691896258*(1.732050807568877*fhatAL[1]-1.732050807568877*fhatAL[3])-1.0*fhatAL[2]+fhatAL[0]), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.5*((-1.0*(fhatAL[3]+fhatAL[1]))+fhatAL[2]+fhatAL[0]), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.5*(1.0*(fhatAL[3]+fhatAL[1])+fhatAL[2]+fhatAL[0]), limQuad[3])); 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  double diffFac[2]; 
  diffFac[0] = (BmagInv[1]*nuVtSqSum[1]+BmagInv[0]*nuVtSqSum[0])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 
  diffFac[1] = (BmagInv[0]*nuVtSqSum[1]+nuVtSqSum[0]*BmagInv[1])*(1.414213562373095*wl[2]+0.7071067811865475*dxvl[2])*m_; 

  double Ghat[8]; 
  for(unsigned short int i=0; i<8; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.530931089239486*diffFac[1]*fr[5])-1.530931089239486*diffFac[1]*fl[5]-1.530931089239486*diffFac[0]*fr[3]-1.530931089239486*diffFac[0]*fl[3]+1.590990257669731*diffFac[1]*fr[1]-1.590990257669731*diffFac[1]*fl[1]+1.590990257669731*diffFac[0]*fr[0]-1.590990257669731*diffFac[0]*fl[0])*rdv+0.7071067811865475*alphaDrSurf[0]*fhatAL[0]; 
  Ghat[1] = ((-1.530931089239486*diffFac[0]*fr[5])-1.530931089239486*diffFac[0]*fl[5]-1.530931089239486*diffFac[1]*fr[3]-1.530931089239486*diffFac[1]*fl[3]+1.590990257669731*diffFac[0]*fr[1]-1.590990257669731*diffFac[0]*fl[1]+1.590990257669731*fr[0]*diffFac[1]-1.590990257669731*fl[0]*diffFac[1])*rdv+0.7071067811865475*alphaDrSurf[0]*fhatAL[1]; 
  Ghat[2] = ((-1.530931089239486*diffFac[1]*fr[7])-1.530931089239486*diffFac[1]*fl[7]-1.530931089239486*diffFac[0]*fr[6]-1.530931089239486*diffFac[0]*fl[6]+1.590990257669731*diffFac[1]*fr[4]-1.590990257669731*diffFac[1]*fl[4]+1.590990257669731*diffFac[0]*fr[2]-1.590990257669731*diffFac[0]*fl[2])*rdv+0.7071067811865475*alphaDrSurf[0]*fhatAL[2]; 
  Ghat[4] = ((-1.530931089239486*diffFac[0]*fr[7])-1.530931089239486*diffFac[0]*fl[7]-1.530931089239486*diffFac[1]*fr[6]-1.530931089239486*diffFac[1]*fl[6]+1.590990257669731*diffFac[0]*fr[4]-1.590990257669731*diffFac[0]*fl[4]+1.590990257669731*diffFac[1]*fr[2]-1.590990257669731*diffFac[1]*fl[2])*rdv+0.7071067811865475*alphaDrSurf[0]*fhatAL[3]; 

  double incr1[8]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = -0.5*Ghat[2]; 
  incr1[3] = 0.8660254037844386*Ghat[0]; 
  incr1[4] = -0.5*Ghat[4]; 
  incr1[5] = 0.8660254037844386*Ghat[1]; 
  incr1[6] = 0.8660254037844386*Ghat[2]; 
  incr1[7] = 0.8660254037844386*Ghat[4]; 

  double incr2[8]; 
  incr2[3] = diffFac[1]*(0.3535533905932737*fl[5]-0.3535533905932737*fr[5])+diffFac[0]*(0.3535533905932737*fl[3]-0.3535533905932737*fr[3])+0.3061862178478971*(diffFac[1]*(fr[1]+fl[1])+diffFac[0]*(fr[0]+fl[0])); 
  incr2[5] = diffFac[0]*(0.3535533905932737*fl[5]-0.3535533905932737*fr[5])+diffFac[1]*(0.3535533905932737*fl[3]-0.3535533905932737*fr[3])+0.3061862178478971*(diffFac[0]*(fr[1]+fl[1])+(fr[0]+fl[0])*diffFac[1]); 
  incr2[6] = diffFac[1]*(0.3535533905932737*fl[7]-0.3535533905932737*fr[7])+diffFac[0]*(0.3535533905932737*fl[6]-0.3535533905932737*fr[6])+0.3061862178478971*(diffFac[1]*(fr[4]+fl[4])+diffFac[0]*(fr[2]+fl[2])); 
  incr2[7] = diffFac[0]*(0.3535533905932737*fl[7]-0.3535533905932737*fr[7])+diffFac[1]*(0.3535533905932737*fl[6]-0.3535533905932737*fr[6])+0.3061862178478971*(diffFac[0]*(fr[4]+fl[4])+diffFac[1]*(fr[2]+fl[2])); 

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

  return std::abs(2.0*wl[2]); 
} 
