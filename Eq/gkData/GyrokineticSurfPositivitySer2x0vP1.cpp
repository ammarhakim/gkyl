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
  double rCtrlL[2], rCtrlR[2];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(1.732050807568877*(2.0*EPSILON+fl[0])-1.0*fl[2]); 
  rCtrlL[1] = (1.732050807568877*fl[3]+3.0*fl[1])/(1.732050807568877*(2.0*EPSILON+fl[0])+fl[2]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(1.732050807568877*(2.0*EPSILON+fr[0])-1.0*fr[2]); 
  rCtrlR[1] = (1.732050807568877*fr[3]+3.0*fr[1])/(1.732050807568877*(2.0*EPSILON+fr[0])+fr[2]); 
  double fCtrlL[2], fCtrlR[2];  // fCtrl = anti-limited f evaluated at each control node on x surface 
  fCtrlL[0] = 0.2886751345948129*(1.732050807568877*fl[0]-1.0*fl[2])*limTheta(rCtrlL[0],1.0); 
  fCtrlL[1] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[0] = 0.2886751345948129*(1.732050807568877*fr[0]-1.0*fr[2])*limTheta(rCtrlR[0],-1.0); 
  fCtrlR[1] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rCtrlR[1],-1.0); 
  double fL_AL[2], fR_AL[2];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.7071067811865475*(fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 1.224744871391589*(fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.7071067811865475*(fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 1.224744871391589*(fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[2], fhatAL[2]; 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = std::max(0., std::min(0.5*(1.414213562373095*fL_AL[0]-1.414213562373095*fL_AL[1]), fl[0]/cflL*0.5)); 
  } else {
  fhatALQuad[0] = std::max(0., std::min(0.5*(1.414213562373095*fR_AL[0]-1.414213562373095*fR_AL[1]), fr[0]/cflR*0.5)); 
  } 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., std::min(0.7071067811865476*(fL_AL[1]+fL_AL[0]), fl[0]/cflL*0.5)); 
  } else {
  fhatALQuad[1] = std::max(0., std::min(0.7071067811865476*(fR_AL[1]+fR_AL[0]), fr[0]/cflR*0.5)); 
  } 
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
  double rCtrlL[2], rCtrlR[2];  // rCtrl=f1/f0 at each control node in dimensions other than y 
  rCtrlL[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(1.732050807568877*(2.0*EPSILON+fl[0])-1.0*fl[1]); 
  rCtrlL[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(1.732050807568877*(2.0*EPSILON+fl[0])+fl[1]); 
  rCtrlR[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(1.732050807568877*(2.0*EPSILON+fr[0])-1.0*fr[1]); 
  rCtrlR[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(1.732050807568877*(2.0*EPSILON+fr[0])+fr[1]); 
  double fCtrlL[2], fCtrlR[2];  // fCtrl = anti-limited f evaluated at each control node on y surface 
  fCtrlL[0] = 0.2886751345948129*(1.732050807568877*fl[0]-1.0*fl[1])*limTheta(rCtrlL[0],1.0); 
  fCtrlL[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rCtrlL[1],1.0); 
  fCtrlR[0] = 0.2886751345948129*(1.732050807568877*fr[0]-1.0*fr[1])*limTheta(rCtrlR[0],-1.0); 
  fCtrlR[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rCtrlR[1],-1.0); 
  double fL_AL[2], fR_AL[2];  // f_AL = mode coefficients of anti-limited f on surface 
  fL_AL[0] = 0.7071067811865475*(fCtrlL[1]+fCtrlL[0]); 
  fL_AL[1] = 1.224744871391589*(fCtrlL[1]-1.0*fCtrlL[0]); 
  fR_AL[0] = 0.7071067811865475*(fCtrlR[1]+fCtrlR[0]); 
  fR_AL[1] = 1.224744871391589*(fCtrlR[1]-1.0*fCtrlR[0]); 
  double alphaQuad; 
  // determine upwinding and enforce limiters at each surface quadrature node 
  double fhatALQuad[2], fhatAL[2]; 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[0] = std::max(0., std::min(0.5*(1.414213562373095*fL_AL[0]-1.414213562373095*fL_AL[1]), fl[0]/cflL*0.5)); 
  } else {
  fhatALQuad[0] = std::max(0., std::min(0.5*(1.414213562373095*fR_AL[0]-1.414213562373095*fR_AL[1]), fr[0]/cflR*0.5)); 
  } 
  alphaQuad = 0.7071067811865475*alpha[0]; 
  if(alphaQuad > 0) {
  fhatALQuad[1] = std::max(0., std::min(0.7071067811865476*(fL_AL[1]+fL_AL[0]), fl[0]/cflL*0.5)); 
  } else {
  fhatALQuad[1] = std::max(0., std::min(0.7071067811865476*(fR_AL[1]+fR_AL[0]), fr[0]/cflR*0.5)); 
  } 
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
