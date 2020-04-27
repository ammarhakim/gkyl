#include <IncompEulerModDecl.h> 
double IncompEulerSurfPositivity2xSer_X_P1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // H/f: Input Hamiltonian/distribution function.
  // out: Incremented output.

  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // Surface-averaged phase velocity in this direction.
  double alpha0 = -0.25*(3.0*Phi[3]-1.732050807568877*Phi[2])*dfac_y; 

  double alpha[2]; 
  alpha[0] = -0.5*(4.242640687119286*Phi[3]-2.449489742783178*Phi[2])*dfac_y; 
  double f0Quad[2]; 
  double f1Quad[2]; 
  double limQuad[2]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
     f0Quad[0] = -0.7071067811865475*(fl[2]-1.0*fl[0]); 
     f1Quad[0] = 0.7071067811865475*(fl[3]-1.0*fl[1]); 
     limQuad[0] = fl[0]/cflL*0.5; 
  } else {
     f0Quad[0] = -0.7071067811865475*(fr[2]-1.0*fr[0]); 
     f1Quad[0] = -0.7071067811865475*(fr[3]-1.0*fr[1]); 
     limQuad[0] = fr[0]/cflR*0.5; 
  } 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
     f0Quad[1] = 0.7071067811865475*(fl[2]+fl[0]); 
     f1Quad[1] = -0.7071067811865475*(fl[3]+fl[1]); 
     limQuad[1] = fl[0]/cflL*0.5; 
  } else {
     f0Quad[1] = 0.7071067811865475*(fr[2]+fr[0]); 
     f1Quad[1] = 0.7071067811865475*(fr[3]+fr[1]); 
     limQuad[1] = fr[0]/cflR*0.5; 
  } 
  double fhat[4]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.7071067811865475*(f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.7071067811865475*(f1Quad[1]+f1Quad[0]); 
  fhat[2] = 0.7071067811865475*(f0Quad[1]-1.0*f0Quad[0]); 
  fhat[3] = 0.7071067811865475*(f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[2];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrl[0] = -(1.0*(1.732050807568877*fhat[3]-3.0*fhat[1]))/(3.464101615137754*EPSILON-1.0*fhat[2]+1.732050807568877*fhat[0]); 
  rCtrl[1] = (1.732050807568877*fhat[3]+3.0*fhat[1])/(3.464101615137754*EPSILON+fhat[2]+1.732050807568877*fhat[0]); 
  double fhatCtrl[2];  // fhatCtrl = anti-limited fhat evaluated at each control node on x surface 
  fhatCtrl[0] = -0.2886751345948129*(fhat[2]-1.732050807568877*fhat[0])*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.2886751345948129*(fhat[2]+1.732050807568877*fhat[0])*limTheta(rCtrl[1],-1.0); 
  double fhatAL[2];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.7071067811865475*(fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 1.224744871391589*(fhatCtrl[1]-1.0*fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[2]; 
  fhatALQuad[0] = std::max(0., std::min(0.5*(1.414213562373095*fhatAL[0]-1.414213562373095*fhatAL[1]), limQuad[0])); 
  fhatALQuad[1] = std::max(0., std::min(0.7071067811865476*(fhatAL[1]+fhatAL[0]), limQuad[1])); 
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
double IncompEulerSurfPositivity2xSer_Y_P1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // H/f: Input Hamiltonian/distribution function.
  // out: Incremented output.

  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double incr[4]; 
  // Surface-averaged phase velocity in this direction.
  double alpha0 = 0.25*(3.0*Phi[3]-1.732050807568877*Phi[1])*dfac_x; 

  double alpha[2]; 
  alpha[0] = 0.5*(4.242640687119286*Phi[3]-2.449489742783178*Phi[1])*dfac_x; 
  double f0Quad[2]; 
  double f1Quad[2]; 
  double limQuad[2]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
     f0Quad[0] = -0.7071067811865475*(fl[1]-1.0*fl[0]); 
     f1Quad[0] = 0.7071067811865475*(fl[3]-1.0*fl[2]); 
     limQuad[0] = fl[0]/cflL*0.5; 
  } else {
     f0Quad[0] = -0.7071067811865475*(fr[1]-1.0*fr[0]); 
     f1Quad[0] = -0.7071067811865475*(fr[3]-1.0*fr[2]); 
     limQuad[0] = fr[0]/cflR*0.5; 
  } 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
     f0Quad[1] = 0.7071067811865475*(fl[1]+fl[0]); 
     f1Quad[1] = -0.7071067811865475*(fl[3]+fl[2]); 
     limQuad[1] = fl[0]/cflL*0.5; 
  } else {
     f0Quad[1] = 0.7071067811865475*(fr[1]+fr[0]); 
     f1Quad[1] = 0.7071067811865475*(fr[3]+fr[2]); 
     limQuad[1] = fr[0]/cflR*0.5; 
  } 
  double fhat[4]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.7071067811865475*(f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.7071067811865475*(f0Quad[1]-1.0*f0Quad[0]); 
  fhat[2] = 0.7071067811865475*(f1Quad[1]+f1Quad[0]); 
  fhat[3] = 0.7071067811865475*(f1Quad[1]-1.0*f1Quad[0]); 
  double rCtrl[2];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrl[0] = -(1.0*(1.732050807568877*fhat[3]-3.0*fhat[2]))/(3.464101615137754*EPSILON-1.0*fhat[1]+1.732050807568877*fhat[0]); 
  rCtrl[1] = (1.732050807568877*fhat[3]+3.0*fhat[2])/(3.464101615137754*EPSILON+fhat[1]+1.732050807568877*fhat[0]); 
  double fhatCtrl[2];  // fhatCtrl = anti-limited fhat evaluated at each control node on y surface 
  fhatCtrl[0] = -0.2886751345948129*(fhat[1]-1.732050807568877*fhat[0])*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = 0.2886751345948129*(fhat[1]+1.732050807568877*fhat[0])*limTheta(rCtrl[1],-1.0); 
  double fhatAL[2];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.7071067811865475*(fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 1.224744871391589*(fhatCtrl[1]-1.0*fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[2]; 
  fhatALQuad[0] = std::max(0., std::min(0.5*(1.414213562373095*fhatAL[0]-1.414213562373095*fhatAL[1]), limQuad[0])); 
  fhatALQuad[1] = std::max(0., std::min(0.7071067811865476*(fhatAL[1]+fhatAL[0]), limQuad[1])); 
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
