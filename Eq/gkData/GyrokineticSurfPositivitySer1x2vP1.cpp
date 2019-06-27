#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[4]; 
  alpha[0] = 1.414213562373095*Gradpar[0]*wv; 
  double f0Quad[4]; 
  double f1Quad[4]; 
  double limQuad[4]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.5*alpha[0]; 
  f0Quad[0] = -0.25*((fr[6]-1.0*(fl[6]+fr[3])+fl[3]-1.0*fr[2]+fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[6]+fl[6])+fr[3]+fl[3]+fr[2]+fl[2]-1.0*(fr[0]+fl[0])); 
  f1Quad[0] = -0.25*((fr[7]+fl[7]-1.0*(fr[5]+fl[5]+fr[4]+fl[4])+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[5]-1.0*fl[5]+fr[4]-1.0*(fl[4]+fr[1])+fl[1]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*alpha[0]; 
  f0Quad[1] = 0.25*((fr[6]-1.0*fl[6]+fr[3]-1.0*(fl[3]+fr[2])+fl[2]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[6]+fl[6]+fr[3]+fl[3])+fr[2]+fl[2]+fr[0]+fl[0]); 
  f1Quad[1] = 0.25*((fr[7]+fl[7]+fr[5]+fl[5]-1.0*(fr[4]+fl[4]+fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[5]+fl[5]+fr[4]-1.0*fl[4]+fr[1]-1.0*fl[1]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*alpha[0]; 
  f0Quad[2] = 0.25*((fr[6]-1.0*(fl[6]+fr[3])+fl[3]+fr[2]-1.0*(fl[2]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[6]+fl[6])+fr[3]+fl[3]-1.0*(fr[2]+fl[2])+fr[0]+fl[0]); 
  f1Quad[2] = 0.25*((fr[7]+fl[7]-1.0*(fr[5]+fl[5])+fr[4]+fl[4]-1.0*(fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[5]-1.0*(fl[5]+fr[4])+fl[4]+fr[1]-1.0*fl[1]); 
  limQuad[2] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*alpha[0]; 
  f0Quad[3] = -0.25*((fr[6]-1.0*fl[6]+fr[3]-1.0*fl[3]+fr[2]-1.0*fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[6]+fl[6]+fr[3]+fl[3]+fr[2]+fl[2]+fr[0]+fl[0])); 
  f1Quad[3] = -0.25*((fr[7]+fl[7]+fr[5]+fl[5]+fr[4]+fl[4]+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[5]+fl[5]-1.0*fr[4]+fl[4]-1.0*fr[1]+fl[1]); 
  limQuad[3] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  double fhat[8]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.5*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.5*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[2] = 0.5*(f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[3] = 0.5*(f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.5*(f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[5] = 0.5*(f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[6] = 0.5*(f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[7] = 0.5*(f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[4];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[7]+3.0*fhat[1])-3.0*(fhat[5]+fhat[4])))/(12.0*EPSILON+1.414213562373095*(fhat[6]-1.732050807568877*(fhat[3]+fhat[2])+3.0*fhat[0])); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[7]-3.0*fhat[1])+3.0*(fhat[5]-1.0*fhat[4])))/(12.0*EPSILON-1.414213562373095*(fhat[6]+1.732050807568877*(fhat[3]-1.0*fhat[2])-3.0*fhat[0])); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[7]-3.0*fhat[1])+3.0*(fhat[4]-1.0*fhat[5])))/(12.0*EPSILON+1.414213562373095*((-1.0*fhat[6])+1.732050807568877*(fhat[3]-1.0*fhat[2])+3.0*fhat[0])); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[7]+3.0*fhat[1])+3.0*(fhat[5]+fhat[4])))/(12.0*EPSILON+1.414213562373095*(fhat[6]+1.732050807568877*(fhat[3]+fhat[2])+3.0*fhat[0])); 
  double fhatCtrl[4];  // fhatCtrl = anti-limited fhat evaluated at each control node on x surface 
  fhatCtrl[0] = -0.1178511301977579*(1.732050807568877*(fhat[3]+fhat[2])-1.0*(fhat[6]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = -0.1178511301977579*(fhat[6]+1.732050807568877*(fhat[3]-1.0*fhat[2])-3.0*fhat[0])*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.1178511301977579*((-1.0*fhat[6])+1.732050807568877*(fhat[3]-1.0*fhat[2])+3.0*fhat[0])*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = 0.1178511301977579*(fhat[6]+1.732050807568877*(fhat[3]+fhat[2])+3.0*fhat[0])*limTheta(rCtrl[3],-1.0); 
  double fhatAL[4];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.5*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.8660254037844386*(fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.8660254037844386*(fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 1.5*(fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[4]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.5*(fhatAL[3]-1.0*(fhatAL[2]+fhatAL[1])+fhatAL[0]), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.5*((-1.0*fhatAL[3])+1.732050807568877*(0.5773502691896258*fhatAL[1]-0.5773502691896258*fhatAL[2])+fhatAL[0]), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.5*((-1.0*fhatAL[3])+1.732050807568877*(0.5773502691896258*fhatAL[2]-0.5773502691896258*fhatAL[1])+fhatAL[0]), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.5*(fhatAL[3]+1.0*(fhatAL[2]+fhatAL[1])+fhatAL[0]), limQuad[3])); 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.3535533905932737*alpha[0]*fhatAL[0]*dfac_x; 
  incr[1] = -0.6123724356957944*alpha[0]*fhatAL[0]*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*fhatAL[1]*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*fhatAL[2]*dfac_x; 
  incr[4] = -0.6123724356957944*alpha[0]*fhatAL[1]*dfac_x; 
  incr[5] = -0.6123724356957944*alpha[0]*fhatAL[2]*dfac_x; 
  incr[6] = 0.3535533905932737*alpha[0]*fhatAL[3]*dfac_x; 
  incr[7] = -0.6123724356957944*alpha[0]*fhatAL[3]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[4]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  double f0Quad[4]; 
  double f1Quad[4]; 
  double limQuad[4]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.5*alpha[0]; 
  f0Quad[0] = -0.25*((fr[5]-1.0*(fl[5]+fr[3])+fl[3]-1.0*fr[1]+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5])+fr[3]+fl[3]+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[0] = -0.25*((fr[7]+fl[7]-1.0*(fr[6]+fl[6]+fr[4]+fl[4])+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[6]-1.0*fl[6]+fr[4]-1.0*(fl[4]+fr[2])+fl[2]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*alpha[0]; 
  f0Quad[1] = 0.25*((fr[5]-1.0*fl[5]+fr[3]-1.0*(fl[3]+fr[1])+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5]+fr[3]+fl[3])+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[1] = 0.25*((fr[7]+fl[7]+fr[6]+fl[6]-1.0*(fr[4]+fl[4]+fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]+fr[4]-1.0*fl[4]+fr[2]-1.0*fl[2]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*alpha[0]; 
  f0Quad[2] = 0.25*((fr[5]-1.0*(fl[5]+fr[3])+fl[3]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5])+fr[3]+fl[3]-1.0*(fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[2] = 0.25*((fr[7]+fl[7]-1.0*(fr[6]+fl[6])+fr[4]+fl[4]-1.0*(fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[6]-1.0*(fl[6]+fr[4])+fl[4]+fr[2]-1.0*fl[2]); 
  limQuad[2] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*alpha[0]; 
  f0Quad[3] = -0.25*((fr[5]-1.0*fl[5]+fr[3]-1.0*fl[3]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5]+fr[3]+fl[3]+fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[3] = -0.25*((fr[7]+fl[7]+fr[6]+fl[6]+fr[4]+fl[4]+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]-1.0*fr[4]+fl[4]-1.0*fr[2]+fl[2]); 
  limQuad[3] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
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
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[7]+3.0*fhat[2])-3.0*(fhat[6]+fhat[4])))/(12.0*EPSILON+1.414213562373095*(fhat[5]-1.732050807568877*(fhat[3]+fhat[1])+3.0*fhat[0])); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[7]-3.0*fhat[2])+3.0*(fhat[6]-1.0*fhat[4])))/(12.0*EPSILON-1.414213562373095*(fhat[5]+1.732050807568877*(fhat[3]-1.0*fhat[1])-3.0*fhat[0])); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[7]-3.0*fhat[2])+3.0*(fhat[4]-1.0*fhat[6])))/(12.0*EPSILON+1.414213562373095*((-1.0*fhat[5])+1.732050807568877*(fhat[3]-1.0*fhat[1])+3.0*fhat[0])); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[7]+3.0*fhat[2])+3.0*(fhat[6]+fhat[4])))/(12.0*EPSILON+1.414213562373095*(fhat[5]+1.732050807568877*(fhat[3]+fhat[1])+3.0*fhat[0])); 
  double fhatCtrl[4];  // fhatCtrl = anti-limited fhat evaluated at each control node on vx surface 
  fhatCtrl[0] = -0.1178511301977579*(1.732050807568877*(fhat[3]+fhat[1])-1.0*(fhat[5]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = -0.1178511301977579*(fhat[5]+1.732050807568877*(fhat[3]-1.0*fhat[1])-3.0*fhat[0])*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.1178511301977579*((-1.0*fhat[5])+1.732050807568877*(fhat[3]-1.0*fhat[1])+3.0*fhat[0])*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = 0.1178511301977579*(fhat[5]+1.732050807568877*(fhat[3]+fhat[1])+3.0*fhat[0])*limTheta(rCtrl[3],-1.0); 
  double fhatAL[4];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.5*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.8660254037844386*(fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.8660254037844386*(fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 1.5*(fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[4]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.5*(fhatAL[3]-1.0*(fhatAL[2]+fhatAL[1])+fhatAL[0]), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.5*((-1.0*fhatAL[3])+1.732050807568877*(0.5773502691896258*fhatAL[1]-0.5773502691896258*fhatAL[2])+fhatAL[0]), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.5*((-1.0*fhatAL[3])+1.732050807568877*(0.5773502691896258*fhatAL[2]-0.5773502691896258*fhatAL[1])+fhatAL[0]), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.5*(fhatAL[3]+1.0*(fhatAL[2]+fhatAL[1])+fhatAL[0]), limQuad[3])); 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.3535533905932737*alpha[0]*fhatAL[0]*dfac_v; 
  incr[1] = 0.3535533905932737*alpha[0]*fhatAL[1]*dfac_v; 
  incr[2] = -0.6123724356957944*alpha[0]*fhatAL[0]*dfac_v; 
  incr[3] = 0.3535533905932737*alpha[0]*fhatAL[2]*dfac_v; 
  incr[4] = -0.6123724356957944*alpha[0]*fhatAL[1]*dfac_v; 
  incr[5] = 0.3535533905932737*alpha[0]*fhatAL[3]*dfac_v; 
  incr[6] = -0.6123724356957944*alpha[0]*fhatAL[2]*dfac_v; 
  incr[7] = -0.6123724356957944*alpha[0]*fhatAL[3]*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.1767766952966368*wv*(1.732050807568877*(Bmag[1]*(3.0*BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x*m_*wv+2.0*Gradpar[1]*q_)-1.0*(3.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x*m_*wv+2.0*Gradpar[0]*q_)))/q_; 

  double alpha[4]; 
  alpha[0] = -(0.5*(2.449489742783178*(Bmag[1]*(3.0*BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x*m_*wv2+2.0*Gradpar[1]*q_*wv)+1.414213562373095*((-3.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x*m_*wv2)-2.0*Gradpar[0]*q_*wv)))/q_; 
  alpha[1] = (0.5*Bmag[1]*(1.414213562373095*((-3.0*BmagInv[1]*geoY[1])-1.0*BmagInv[0]*geoY[0])+2.449489742783178*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))*dfac_x*m_*wv)/(dfac_v*q_); 
  double f0Quad[4]; 
  double f1Quad[4]; 
  double limQuad[4]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.5*alpha[0]-0.5*alpha[1]; 
  f0Quad[0] = -0.25*((fr[6]-1.0*(fl[6]+fr[3])+fl[3]-1.0*fr[2]+fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[6]+fl[6])+fr[3]+fl[3]+fr[2]+fl[2]-1.0*(fr[0]+fl[0])); 
  f1Quad[0] = -0.25*((fr[7]+fl[7]-1.0*(fr[5]+fl[5]+fr[4]+fl[4])+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[5]-1.0*fl[5]+fr[4]-1.0*(fl[4]+fr[1])+fl[1]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*(alpha[1]+alpha[0]); 
  f0Quad[1] = 0.25*((fr[6]-1.0*fl[6]+fr[3]-1.0*(fl[3]+fr[2])+fl[2]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[6]+fl[6]+fr[3]+fl[3])+fr[2]+fl[2]+fr[0]+fl[0]); 
  f1Quad[1] = 0.25*((fr[7]+fl[7]+fr[5]+fl[5]-1.0*(fr[4]+fl[4]+fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[5]+fl[5]+fr[4]-1.0*fl[4]+fr[1]-1.0*fl[1]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*alpha[0]-0.5*alpha[1]; 
  f0Quad[2] = 0.25*((fr[6]-1.0*(fl[6]+fr[3])+fl[3]+fr[2]-1.0*(fl[2]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[6]+fl[6])+fr[3]+fl[3]-1.0*(fr[2]+fl[2])+fr[0]+fl[0]); 
  f1Quad[2] = 0.25*((fr[7]+fl[7]-1.0*(fr[5]+fl[5])+fr[4]+fl[4]-1.0*(fr[1]+fl[1]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[5]-1.0*(fl[5]+fr[4])+fl[4]+fr[1]-1.0*fl[1]); 
  limQuad[2] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*(alpha[1]+alpha[0]); 
  f0Quad[3] = -0.25*((fr[6]-1.0*fl[6]+fr[3]-1.0*fl[3]+fr[2]-1.0*fl[2]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[6]+fl[6]+fr[3]+fl[3]+fr[2]+fl[2]+fr[0]+fl[0])); 
  f1Quad[3] = -0.25*((fr[7]+fl[7]+fr[5]+fl[5]+fr[4]+fl[4]+fr[1]+fl[1])*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[5]+fl[5]-1.0*fr[4]+fl[4]-1.0*fr[1]+fl[1]); 
  limQuad[3] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  double fhat[8]; // (volume) mode coefficients of fhat 
  fhat[0] = 0.5*(f0Quad[3]+f0Quad[2]+f0Quad[1]+f0Quad[0]); 
  fhat[1] = 0.5*(f1Quad[3]+f1Quad[2]+f1Quad[1]+f1Quad[0]); 
  fhat[2] = 0.5*(f0Quad[3]-1.0*f0Quad[2]+f0Quad[1]-1.0*f0Quad[0]); 
  fhat[3] = 0.5*(f0Quad[3]+f0Quad[2]-1.0*(f0Quad[1]+f0Quad[0])); 
  fhat[4] = 0.5*(f1Quad[3]-1.0*f1Quad[2]+f1Quad[1]-1.0*f1Quad[0]); 
  fhat[5] = 0.5*(f1Quad[3]+f1Quad[2]-1.0*(f1Quad[1]+f1Quad[0])); 
  fhat[6] = 0.5*(f0Quad[3]-1.0*(f0Quad[2]+f0Quad[1])+f0Quad[0]); 
  fhat[7] = 0.5*(f1Quad[3]-1.0*(f1Quad[2]+f1Quad[1])+f1Quad[0]); 
  double rCtrl[4];  // rCtrl=f1/f0 at each control node in dimensions other than x 
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[7]+3.0*fhat[1])-3.0*(fhat[5]+fhat[4])))/(12.0*EPSILON+1.414213562373095*(fhat[6]-1.732050807568877*(fhat[3]+fhat[2])+3.0*fhat[0])); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[7]-3.0*fhat[1])+3.0*(fhat[5]-1.0*fhat[4])))/(12.0*EPSILON-1.414213562373095*(fhat[6]+1.732050807568877*(fhat[3]-1.0*fhat[2])-3.0*fhat[0])); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[7]-3.0*fhat[1])+3.0*(fhat[4]-1.0*fhat[5])))/(12.0*EPSILON+1.414213562373095*((-1.0*fhat[6])+1.732050807568877*(fhat[3]-1.0*fhat[2])+3.0*fhat[0])); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[7]+3.0*fhat[1])+3.0*(fhat[5]+fhat[4])))/(12.0*EPSILON+1.414213562373095*(fhat[6]+1.732050807568877*(fhat[3]+fhat[2])+3.0*fhat[0])); 
  double fhatCtrl[4];  // fhatCtrl = anti-limited fhat evaluated at each control node on x surface 
  fhatCtrl[0] = -0.1178511301977579*(1.732050807568877*(fhat[3]+fhat[2])-1.0*(fhat[6]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = -0.1178511301977579*(fhat[6]+1.732050807568877*(fhat[3]-1.0*fhat[2])-3.0*fhat[0])*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.1178511301977579*((-1.0*fhat[6])+1.732050807568877*(fhat[3]-1.0*fhat[2])+3.0*fhat[0])*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = 0.1178511301977579*(fhat[6]+1.732050807568877*(fhat[3]+fhat[2])+3.0*fhat[0])*limTheta(rCtrl[3],-1.0); 
  double fhatAL[4];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.5*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.8660254037844386*(fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.8660254037844386*(fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 1.5*(fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[4]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.5*(fhatAL[3]-1.0*(fhatAL[2]+fhatAL[1])+fhatAL[0]), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.5*((-1.0*fhatAL[3])+1.732050807568877*(0.5773502691896258*fhatAL[1]-0.5773502691896258*fhatAL[2])+fhatAL[0]), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.5*((-1.0*fhatAL[3])+1.732050807568877*(0.5773502691896258*fhatAL[2]-0.5773502691896258*fhatAL[1])+fhatAL[0]), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.5*(fhatAL[3]+1.0*(fhatAL[2]+fhatAL[1])+fhatAL[0]), limQuad[3])); 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.3535533905932737*(alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_x; 
  incr[1] = -0.6123724356957944*(alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_x; 
  incr[2] = 0.3535533905932737*(alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_x; 
  incr[3] = 0.3535533905932737*(alpha[1]*fhatAL[3]+alpha[0]*fhatAL[2])*dfac_x; 
  incr[4] = -0.6123724356957944*(alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_x; 
  incr[5] = -0.6123724356957944*(alpha[1]*fhatAL[3]+alpha[0]*fhatAL[2])*dfac_x; 
  incr[6] = 0.3535533905932737*(alpha[0]*fhatAL[3]+alpha[1]*fhatAL[2])*dfac_x; 
  incr[7] = -0.6123724356957944*(alpha[0]*fhatAL[3]+alpha[1]*fhatAL[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*dfac_x*(Bmag[1]*wm+Phi[1]*q_)*(3.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_v*dfac_x*m_*wv-1.0*(3.464101615137754*Gradpar[0]*dfac_v*q_+3.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x*m_)))/(dfac_v*m_*q_); 

  double alpha[4]; 
  alpha[0] = -(0.5*(dfac_v*(3.464101615137754*Gradpar[0]*dfac_x*(Bmag[1]*q_*wm+Phi[1]*q2)-3.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_*(Bmag[1]*wm+Phi[1]*q_)*wv)+3.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_*(Bmag[1]*wm+Phi[1]*q_)))/(dfac_v*m_*q_); 
  alpha[1] = -(0.5*(dfac_v*(3.464101615137754*Gradpar[1]*dfac_x*(Bmag[1]*q_*wm+Phi[1]*q2)-3.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_*(Bmag[1]*wm+Phi[1]*q_)*wv)+3.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_*(Bmag[1]*wm+Phi[1]*q_)))/(dfac_v*m_*q_); 
  alpha[2] = (0.5*Bmag[1]*(1.732050807568877*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_*(dfac_v*wv-1.0)-2.0*Gradpar[0]*dfac_v*dfac_x*q_))/(dfac_m*dfac_v*m_*q_); 
  alpha[3] = (0.5*Bmag[1]*(1.732050807568877*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_*(dfac_v*wv-1.0)-2.0*Gradpar[1]*dfac_v*dfac_x*q_))/(dfac_m*dfac_v*m_*q_); 
  double f0Quad[4]; 
  double f1Quad[4]; 
  double limQuad[4]; 
  double alphaQuad; 
  // determine upwinding at each surface quadrature node 
  alphaQuad = 0.5*alpha[3]-0.5*(alpha[2]+alpha[1])+0.5*alpha[0]; 
  f0Quad[0] = -0.25*((fr[5]-1.0*(fl[5]+fr[3])+fl[3]-1.0*fr[1]+fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5])+fr[3]+fl[3]+fr[1]+fl[1]-1.0*(fr[0]+fl[0])); 
  f1Quad[0] = -0.25*((fr[7]+fl[7]-1.0*(fr[6]+fl[6]+fr[4]+fl[4])+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[6]-1.0*fl[6]+fr[4]-1.0*(fl[4]+fr[2])+fl[2]); 
  limQuad[0] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*(alpha[1]+alpha[0])-0.5*(alpha[3]+alpha[2]); 
  f0Quad[1] = 0.25*((fr[5]-1.0*fl[5]+fr[3]-1.0*(fl[3]+fr[1])+fl[1]-1.0*fr[0]+fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5]+fr[3]+fl[3])+fr[1]+fl[1]+fr[0]+fl[0]); 
  f1Quad[1] = 0.25*((fr[7]+fl[7]+fr[6]+fl[6]-1.0*(fr[4]+fl[4]+fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]+fr[4]-1.0*fl[4]+fr[2]-1.0*fl[2]); 
  limQuad[1] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = (-0.5*alpha[3])+0.5*alpha[2]-0.5*alpha[1]+0.5*alpha[0]; 
  f0Quad[2] = 0.25*((fr[5]-1.0*(fl[5]+fr[3])+fl[3]+fr[1]-1.0*(fl[1]+fr[0])+fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5])+fr[3]+fl[3]-1.0*(fr[1]+fl[1])+fr[0]+fl[0]); 
  f1Quad[2] = 0.25*((fr[7]+fl[7]-1.0*(fr[6]+fl[6])+fr[4]+fl[4]-1.0*(fr[2]+fl[2]))*sgn(alphaQuad)-1.0*fr[7]+fl[7]+fr[6]-1.0*(fl[6]+fr[4])+fl[4]+fr[2]-1.0*fl[2]); 
  limQuad[2] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
  alphaQuad = 0.5*(alpha[3]+alpha[2]+alpha[1]+alpha[0]); 
  f0Quad[3] = -0.25*((fr[5]-1.0*fl[5]+fr[3]-1.0*fl[3]+fr[1]-1.0*fl[1]+fr[0]-1.0*fl[0])*sgn(alphaQuad)-1.0*(fr[5]+fl[5]+fr[3]+fl[3]+fr[1]+fl[1]+fr[0]+fl[0])); 
  f1Quad[3] = -0.25*((fr[7]+fl[7]+fr[6]+fl[6]+fr[4]+fl[4]+fr[2]+fl[2])*sgn(alphaQuad)-1.0*fr[7]+fl[7]-1.0*fr[6]+fl[6]-1.0*fr[4]+fl[4]-1.0*fr[2]+fl[2]); 
  limQuad[3] = .5*(fl[0]/cflL+fr[0]/cflR + sgn(alphaQuad)*(fl[0]/cflL-fr[0]/cflR))*0.3535533905932737; 
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
  rCtrl[0] = (1.414213562373095*(1.732050807568877*(fhat[7]+3.0*fhat[2])-3.0*(fhat[6]+fhat[4])))/(12.0*EPSILON+1.414213562373095*(fhat[5]-1.732050807568877*(fhat[3]+fhat[1])+3.0*fhat[0])); 
  rCtrl[1] = -(1.414213562373095*(1.732050807568877*(fhat[7]-3.0*fhat[2])+3.0*(fhat[6]-1.0*fhat[4])))/(12.0*EPSILON-1.414213562373095*(fhat[5]+1.732050807568877*(fhat[3]-1.0*fhat[1])-3.0*fhat[0])); 
  rCtrl[2] = -(1.414213562373095*(1.732050807568877*(fhat[7]-3.0*fhat[2])+3.0*(fhat[4]-1.0*fhat[6])))/(12.0*EPSILON+1.414213562373095*((-1.0*fhat[5])+1.732050807568877*(fhat[3]-1.0*fhat[1])+3.0*fhat[0])); 
  rCtrl[3] = (1.414213562373095*(1.732050807568877*(fhat[7]+3.0*fhat[2])+3.0*(fhat[6]+fhat[4])))/(12.0*EPSILON+1.414213562373095*(fhat[5]+1.732050807568877*(fhat[3]+fhat[1])+3.0*fhat[0])); 
  double fhatCtrl[4];  // fhatCtrl = anti-limited fhat evaluated at each control node on vx surface 
  fhatCtrl[0] = -0.1178511301977579*(1.732050807568877*(fhat[3]+fhat[1])-1.0*(fhat[5]+3.0*fhat[0]))*limTheta(rCtrl[0],-1.0); 
  fhatCtrl[1] = -0.1178511301977579*(fhat[5]+1.732050807568877*(fhat[3]-1.0*fhat[1])-3.0*fhat[0])*limTheta(rCtrl[1],-1.0); 
  fhatCtrl[2] = 0.1178511301977579*((-1.0*fhat[5])+1.732050807568877*(fhat[3]-1.0*fhat[1])+3.0*fhat[0])*limTheta(rCtrl[2],-1.0); 
  fhatCtrl[3] = 0.1178511301977579*(fhat[5]+1.732050807568877*(fhat[3]+fhat[1])+3.0*fhat[0])*limTheta(rCtrl[3],-1.0); 
  double fhatAL[4];  // fhatAL = mode coefficients of anti-limited f on surface 
  fhatAL[0] = 0.5*(fhatCtrl[3]+fhatCtrl[2]+fhatCtrl[1]+fhatCtrl[0]); 
  fhatAL[1] = 0.8660254037844386*(fhatCtrl[3]-1.0*fhatCtrl[2]+fhatCtrl[1]-1.0*fhatCtrl[0]); 
  fhatAL[2] = 0.8660254037844386*(fhatCtrl[3]+fhatCtrl[2]-1.0*(fhatCtrl[1]+fhatCtrl[0])); 
  fhatAL[3] = 1.5*(fhatCtrl[3]-1.0*(fhatCtrl[2]+fhatCtrl[1])+fhatCtrl[0]); 
  // enforce limiters at surface quadrature nodes 
  double fhatALQuad[4]; 
  fhatALQuad[0] = std::max(0.0, std::min(0.5*(fhatAL[3]-1.0*(fhatAL[2]+fhatAL[1])+fhatAL[0]), limQuad[0])); 
  fhatALQuad[1] = std::max(0.0, std::min(0.5*((-1.0*fhatAL[3])+1.732050807568877*(0.5773502691896258*fhatAL[1]-0.5773502691896258*fhatAL[2])+fhatAL[0]), limQuad[1])); 
  fhatALQuad[2] = std::max(0.0, std::min(0.5*((-1.0*fhatAL[3])+1.732050807568877*(0.5773502691896258*fhatAL[2]-0.5773502691896258*fhatAL[1])+fhatAL[0]), limQuad[2])); 
  fhatALQuad[3] = std::max(0.0, std::min(0.5*(fhatAL[3]+1.0*(fhatAL[2]+fhatAL[1])+fhatAL[0]), limQuad[3])); 
  fhatAL[0] = 0.5*(fhatALQuad[3]+fhatALQuad[2]+fhatALQuad[1]+fhatALQuad[0]); 
  fhatAL[1] = 0.5*(fhatALQuad[3]-1.0*fhatALQuad[2]+fhatALQuad[1]-1.0*fhatALQuad[0]); 
  fhatAL[2] = 0.5*(fhatALQuad[3]+fhatALQuad[2]-1.0*(fhatALQuad[1]+fhatALQuad[0])); 
  fhatAL[3] = 0.5*(fhatALQuad[3]-1.0*(fhatALQuad[2]+fhatALQuad[1])+fhatALQuad[0]); 

  // begin surface update 
 
  incr[0] = 0.3535533905932737*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[2] = -0.6123724356957944*(alpha[3]*fhatAL[3]+alpha[2]*fhatAL[2]+alpha[1]*fhatAL[1]+alpha[0]*fhatAL[0])*dfac_v; 
  incr[3] = 0.3535533905932737*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[4] = -0.6123724356957944*(alpha[2]*fhatAL[3]+fhatAL[2]*alpha[3]+alpha[0]*fhatAL[1]+fhatAL[0]*alpha[1])*dfac_v; 
  incr[5] = 0.3535533905932737*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 
  incr[6] = -0.6123724356957944*(alpha[1]*fhatAL[3]+fhatAL[1]*alpha[3]+alpha[0]*fhatAL[2]+fhatAL[0]*alpha[2])*dfac_v; 
  incr[7] = -0.6123724356957944*(alpha[0]*fhatAL[3]+fhatAL[0]*alpha[3]+alpha[1]*fhatAL[2]+fhatAL[1]*alpha[2])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  return std::abs(alpha0); 
} 
