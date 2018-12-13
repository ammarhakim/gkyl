#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity2x0vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.125*BmagInv[0]*(1.732050807568877*Phi[2]-3.0*Phi[3])*dfac_y; 

  double alpha[2]; 
  alpha[0] = -0.3535533905932737*BmagInv[0]*(1.732050807568877*Phi[2]-3.0*Phi[3])*dfac_y; 
  double f0Quad[2]; 
  double f1Quad[2]; 
  double limQuad[2]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  f0Quad[0] = 0.25*((1.414213562373095*fr[2]-1.414213562373095*(fl[2]+fr[0])+1.414213562373095*fl[0])*sgn(alphaQuad)-1.414213562373095*(fr[2]+fl[2])+1.414213562373095*(fr[0]+fl[0])); 
  f1Quad[0] = 0.25*((1.414213562373095*(fr[3]+fl[3])-1.414213562373095*(fr[1]+fl[1]))*sgn(alphaQuad)-1.414213562373095*fr[3]+1.414213562373095*(fl[3]+fr[1])-1.414213562373095*fl[1]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.5; 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  f0Quad[1] = -0.25*((1.414213562373095*fr[2]-1.414213562373095*fl[2]+1.414213562373095*fr[0]-1.414213562373095*fl[0])*sgn(alphaQuad)-1.414213562373095*(fr[2]+fl[2]+fr[0]+fl[0])); 
  f1Quad[1] = -0.25*(1.414213562373095*(fr[3]+fl[3]+fr[1]+fl[1])*sgn(alphaQuad)-1.414213562373095*fr[3]+1.414213562373095*fl[3]-1.414213562373095*fr[1]+1.414213562373095*fl[1]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.5; 
  double fhat[4]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.7071067811865475*(f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.7071067811865475*(f1Quad[1]+f1Quad[0]); 
  fhat[2] = 0.7071067811865475*(f0Quad[1]-1.0*f0Quad[0]); 
  fhat[3] = 0.7071067811865475*(f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[2];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrl[0] = -(1.0*(1.732050807568877*fhat[3]-3.0*fhat[1]))/(1.732050807568877*(2.0*EPSILON+fhat[0])-1.0*fhat[2]); 
  rCtrl[1] = (1.732050807568877*fhat[3]+3.0*fhat[1])/(1.732050807568877*(2.0*EPSILON+fhat[0])+fhat[2]); 
  double fhatCtrl[2];  // fhatCtrl = anti-limited fhat evaluated at each control node on x surface 
  fhatCtrl[0] = 0.2886751345948129*(1.732050807568877*fhat[0]-1.0*fhat[2])*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.2886751345948129*(fhat[2]+1.732050807568877*fhat[0])*limTheta(rCtrl[1],-1.0); 
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
 
  incr[0] = 0.5*alpha[0]*fhatAL[0]*dfac_x; 
  incr[1] = -0.8660254037844386*alpha[0]*fhatAL[0]*dfac_x; 
  incr[2] = 0.5*alpha[0]*fhatAL[1]*dfac_x; 
  incr[3] = -0.8660254037844386*alpha[0]*fhatAL[1]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity2x0vSer_Y_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.125*BmagInv[0]*(1.732050807568877*Phi[1]-3.0*Phi[3])*dfac_x; 

  double alpha[2]; 
  alpha[0] = 0.3535533905932737*BmagInv[0]*(1.732050807568877*Phi[1]-3.0*Phi[3])*dfac_x; 
  double f0Quad[2]; 
  double f1Quad[2]; 
  double limQuad[2]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  f0Quad[0] = 0.25*((1.414213562373095*fr[1]-1.414213562373095*(fl[1]+fr[0])+1.414213562373095*fl[0])*sgn(alphaQuad)-1.414213562373095*(fr[1]+fl[1])+1.414213562373095*(fr[0]+fl[0])); 
  f1Quad[0] = 0.25*((1.414213562373095*(fr[3]+fl[3])-1.414213562373095*(fr[2]+fl[2]))*sgn(alphaQuad)-1.414213562373095*fr[3]+1.414213562373095*(fl[3]+fr[2])-1.414213562373095*fl[2]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.5; 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  f0Quad[1] = -0.25*((1.414213562373095*fr[1]-1.414213562373095*fl[1]+1.414213562373095*fr[0]-1.414213562373095*fl[0])*sgn(alphaQuad)-1.414213562373095*(fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[1] = -0.25*(1.414213562373095*(fr[3]+fl[3]+fr[2]+fl[2])*sgn(alphaQuad)-1.414213562373095*fr[3]+1.414213562373095*fl[3]-1.414213562373095*fr[2]+1.414213562373095*fl[2]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.5; 
  double fhat[4]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.7071067811865475*(f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.7071067811865475*(f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.7071067811865475*(f1Quad[1]+f1Quad[0]); 
  fhat[3] = 0.7071067811865475*(f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[2];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrl[0] = -(1.0*(1.732050807568877*fhat[3]-3.0*fhat[2]))/(1.732050807568877*(2.0*EPSILON+fhat[0])-1.0*fhat[1]); 
  rCtrl[1] = (1.732050807568877*fhat[3]+3.0*fhat[2])/(1.732050807568877*(2.0*EPSILON+fhat[0])+fhat[1]); 
  double fhatCtrl[2];  // fhatCtrl = anti-limited fhat evaluated at each control node on y surface 
  fhatCtrl[0] = 0.2886751345948129*(1.732050807568877*fhat[0]-1.0*fhat[1])*limTheta(rCtrl[0],-1.0); 
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
 
  incr[0] = 0.5*alpha[0]*fhatAL[0]*dfac_y; 
  incr[1] = 0.5*alpha[0]*fhatAL[1]*dfac_y; 
  incr[2] = -0.8660254037844386*alpha[0]*fhatAL[0]*dfac_y; 
  incr[3] = -0.8660254037844386*alpha[0]*fhatAL[1]*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  return std::abs(alpha0); 
} 
