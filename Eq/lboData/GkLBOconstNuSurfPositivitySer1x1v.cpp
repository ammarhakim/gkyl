#include <GkLBOModDecl.h> 
double GkLBOconstNuSurfPositivity1x1vSer_Vpar_P1(const double m_, const double cfll, const double cflr, const double *wl, const double *wr, const double *dxvl, const double *dxvr, const double *BmagInv, const double nuSum, const double vMuMidMax, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fr, double *outl, double *outr) 
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

  double f0Quad[2]; 
  double f1Quad[2]; 
  double limQuad[2]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.7071067811865475*alphaDrSurf[1]-0.7071067811865475*alphaDrSurf[0]; 
  f0Quad[0] = 0.25*((1.414213562373095*fr[1]-1.414213562373095*(fl[1]+fr[0])+1.414213562373095*fl[0])*sgn(alphaQuad)-1.414213562373095*(fr[1]+fl[1])+1.414213562373095*(fr[0]+fl[0])); 
  f1Quad[0] = 0.25*((1.414213562373095*(fr[3]+fl[3])-1.414213562373095*(fr[2]+fl[2]))*sgn(alphaQuad)-1.414213562373095*fr[3]+1.414213562373095*(fl[3]+fr[2])-1.414213562373095*fl[2]); 
  limQuad[0] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.5; 
  alphaQuad = -0.7071067811865475*(alphaDrSurf[1]+alphaDrSurf[0]); 
  f0Quad[1] = -0.25*((1.414213562373095*fr[1]-1.414213562373095*fl[1]+1.414213562373095*fr[0]-1.414213562373095*fl[0])*sgn(alphaQuad)-1.414213562373095*(fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[1] = -0.25*(1.414213562373095*(fr[3]+fl[3]+fr[2]+fl[2])*sgn(alphaQuad)-1.414213562373095*fr[3]+1.414213562373095*fl[3]-1.414213562373095*fr[2]+1.414213562373095*fl[2]); 
  limQuad[1] = .5*(fl[0]/cfll+fr[0]/cflr + sgn(alphaQuad)*(fl[0]/cfll-fr[0]/cflr))*0.5; 
  double fhat[4]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.7071067811865475*(f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.7071067811865475*(f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.7071067811865475*(f1Quad[1]+f1Quad[0]); 
  fhat[3] = 0.7071067811865475*(f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[2];  // rCtrl=f1/f0 at each control node in dimensions other than vx 
  rCtrl[0] = -(1.0*(1.732050807568877*fhat[3]-3.0*fhat[2]))/(3.464101615137754*EPSILON-1.0*fhat[1]+1.732050807568877*fhat[0]); 
  rCtrl[1] = (1.732050807568877*fhat[3]+3.0*fhat[2])/(3.464101615137754*EPSILON+fhat[1]+1.732050807568877*fhat[0]); 
  double fhatCtrl[2];  // fhatCtrl = anti-limited fhat evaluated at each control node on vx surface 
  fhatCtrl[0] = -0.2886751345948129*(fhat[1]-1.732050807568877*fhat[0])*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.2886751345948129*(fhat[1]+1.732050807568877*fhat[0])*limTheta(rCtrl[1],-1.0); 
  double fhatAL[2];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.7071067811865475*(fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 1.224744871391589*(fhatCtrl[1]-1.0*fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[2]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.5*(1.414213562373095*fhatAL[0]-1.414213562373095*fhatAL[1]), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.7071067811865476*(fhatAL[1]+fhatAL[0]), limQuad[1])); 
  fhatAL[0] = 0.7071067811865475*(fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.7071067811865475*(fhatALQuad[1]-1.0*fhatALQuad[0]); 

  // begin surface update 
 
  double Ghat[4]; 
  for(unsigned short int i=0; i<4; ++i){ 
    Ghat[i]=0.0; 
  }; 

  Ghat[0] = ((-1.530931089239486*nuVtSqSum[1]*fr[3])-1.530931089239486*nuVtSqSum[1]*fl[3]-1.530931089239486*nuVtSqSum[0]*fr[2]-1.530931089239486*nuVtSqSum[0]*fl[2]+1.590990257669731*fr[1]*nuVtSqSum[1]-1.590990257669731*fl[1]*nuVtSqSum[1]+1.590990257669731*fr[0]*nuVtSqSum[0]-1.590990257669731*fl[0]*nuVtSqSum[0])*rdv+alphaDrSurf[1]*fhatAL[1]+alphaDrSurf[0]*fhatAL[0]; 
  Ghat[1] = ((-1.530931089239486*nuVtSqSum[0]*fr[3])-1.530931089239486*nuVtSqSum[0]*fl[3]-1.530931089239486*nuVtSqSum[1]*fr[2]-1.530931089239486*nuVtSqSum[1]*fl[2]+1.590990257669731*fr[0]*nuVtSqSum[1]-1.590990257669731*fl[0]*nuVtSqSum[1]+1.590990257669731*nuVtSqSum[0]*fr[1]-1.590990257669731*nuVtSqSum[0]*fl[1])*rdv+alphaDrSurf[0]*fhatAL[1]+fhatAL[0]*alphaDrSurf[1]; 

  double incr1[4]; 
  incr1[0] = -0.5*Ghat[0]; 
  incr1[1] = -0.5*Ghat[1]; 
  incr1[2] = 0.8660254037844386*Ghat[0]; 
  incr1[3] = 0.8660254037844386*Ghat[1]; 

  double incr2[4]; 
  incr2[2] = nuVtSqSum[1]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[0]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 
  incr2[3] = nuVtSqSum[0]*((-0.3535533905932737*fr[3])+0.3535533905932737*fl[3]+0.3061862178478971*(fr[1]+fl[1]))+nuVtSqSum[1]*((-0.3535533905932737*fr[2])+0.3535533905932737*fl[2]+0.3061862178478971*(fr[0]+fl[0])); 

  outr[0] += incr1[0]*rdv2R; 
  outr[1] += incr1[1]*rdv2R; 
  outr[2] += incr2[2]*rdvSq4R+incr1[2]*rdv2R; 
  outr[3] += incr2[3]*rdvSq4R+incr1[3]*rdv2R; 

  outl[0] += -1.0*incr1[0]*rdv2L; 
  outl[1] += -1.0*incr1[1]*rdv2L; 
  outl[2] += incr1[2]*rdv2L-1.0*incr2[2]*rdvSq4L; 
  outl[3] += incr1[3]*rdv2L-1.0*incr2[3]*rdvSq4L; 

  return std::abs(wl[1]-(0.7071067811865475*nuUSum[0])/nuSum); 
} 
