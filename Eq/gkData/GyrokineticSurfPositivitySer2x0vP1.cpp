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
  double alpha0 = 0.125*BmagInv[0]*(3.0*Phi[3]-1.732050807568877*Phi[2])*dfac_y; 

  double alpha[2]; 
  alpha[0] = 0.3535533905932737*BmagInv[0]*(3.0*Phi[3]-1.732050807568877*Phi[2])*dfac_y; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  double rVal[2];  // rVal=f1/f0 at each control node in dimensions other than x 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each control node on x surface 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflL); 
  fqVal[1] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflL); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[1] = -0.8660254037844386*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[2] = 0.5*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[3] = -0.8660254037844386*alpha[0]*fhatALVal[1]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[2];  // rVal=f1/f0 at each control node in dimensions other than x 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each control node on x surface 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflR); 
  fqVal[1] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflR); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[1] = -0.8660254037844386*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[2] = 0.5*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[3] = -0.8660254037844386*alpha[0]*fhatALVal[1]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
  if(0.7071067811865475*alpha[0] > 0) {
  fupwindQuad[0] = 0.5*(fl[1]+fl[0])-0.5*(fl[3]+fl[2]); 
  } else {
  fupwindQuad[0] = 0.5*fr[3]-0.5*(fr[2]+fr[1])+0.5*fr[0]; 
  }
  if(0.7071067811865475*alpha[0] > 0) {
  fupwindQuad[1] = 0.5*fl[3]-0.5*(fl[2]+fl[1])+0.5*fl[0]; 
  } else {
  fupwindQuad[1] = 0.5*(fr[1]+fr[0])-0.5*(fr[3]+fr[2]); 
  }
  if(0.7071067811865475*alpha[0] > 0) {
  fupwindQuad[2] = 0.5*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[2] = (-0.5*fr[3])+0.5*fr[2]-0.5*fr[1]+0.5*fr[0]; 
  }
  if(0.7071067811865475*alpha[0] > 0) {
  fupwindQuad[3] = (-0.5*fl[3])+0.5*fl[2]-0.5*fl[1]+0.5*fl[0]; 
  } else {
  fupwindQuad[3] = 0.5*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*(fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.5*(fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  double rVal[2];  // rVal=f1/f0 at each control node in dimensions other than x 
  rVal[0] = -(1.0*(1.732050807568877*fupwind[3]-3.0*fupwind[1]))/(3.464101615137754*EPSILON-1.0*fupwind[2]+1.732050807568877*fupwind[0]); 
  rVal[1] = (1.732050807568877*fupwind[3]+3.0*fupwind[1])/(3.464101615137754*EPSILON+fupwind[2]+1.732050807568877*fupwind[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each control node on x surface 
  fqVal[0] = -0.2886751345948129*(fupwind[2]-1.732050807568877*fupwind[0])*limTheta(rVal[0],-1.0,cflR); 
  fqVal[1] = 0.2886751345948129*(fupwind[2]+1.732050807568877*fupwind[0])*limTheta(rVal[1],-1.0,cflR); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[1] = -0.8660254037844386*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[2] = 0.5*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[3] = -0.8660254037844386*alpha[0]*fhatALVal[1]*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
#endif
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
  double alpha0 = -0.125*BmagInv[0]*(3.0*Phi[3]-1.732050807568877*Phi[1])*dfac_x; 

  double alpha[2]; 
  alpha[0] = -0.3535533905932737*BmagInv[0]*(3.0*Phi[3]-1.732050807568877*Phi[1])*dfac_x; 
#if upwindType == SURFAVG 
  if (alpha0>0) { 
  double rVal[2];  // rVal=f1/f0 at each control node in dimensions other than y 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each control node on y surface 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflL); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflL); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*alpha[0]*fhatALVal[0]*dfac_y; 
  incr[1] = 0.5*alpha[0]*fhatALVal[1]*dfac_y; 
  incr[2] = -0.8660254037844386*alpha[0]*fhatALVal[0]*dfac_y; 
  incr[3] = -0.8660254037844386*alpha[0]*fhatALVal[1]*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[2];  // rVal=f1/f0 at each control node in dimensions other than y 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each control node on y surface 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflR); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflR); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*alpha[0]*fhatALVal[0]*dfac_y; 
  incr[1] = 0.5*alpha[0]*fhatALVal[1]*dfac_y; 
  incr[2] = -0.8660254037844386*alpha[0]*fhatALVal[0]*dfac_y; 
  incr[3] = -0.8660254037844386*alpha[0]*fhatALVal[1]*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
#elif upwindType == QUAD 
double fupwind[4];
double fupwindQuad[4];
  if(0.7071067811865475*alpha[0] > 0) {
  fupwindQuad[0] = (-0.5*fl[3])+0.5*fl[2]-0.5*fl[1]+0.5*fl[0]; 
  } else {
  fupwindQuad[0] = 0.5*fr[3]-0.5*(fr[2]+fr[1])+0.5*fr[0]; 
  }
  if(0.7071067811865475*alpha[0] > 0) {
  fupwindQuad[1] = 0.5*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else {
  fupwindQuad[1] = 0.5*(fr[1]+fr[0])-0.5*(fr[3]+fr[2]); 
  }
  if(0.7071067811865475*alpha[0] > 0) {
  fupwindQuad[2] = 0.5*fl[3]-0.5*(fl[2]+fl[1])+0.5*fl[0]; 
  } else {
  fupwindQuad[2] = (-0.5*fr[3])+0.5*fr[2]-0.5*fr[1]+0.5*fr[0]; 
  }
  if(0.7071067811865475*alpha[0] > 0) {
  fupwindQuad[3] = 0.5*(fl[1]+fl[0])-0.5*(fl[3]+fl[2]); 
  } else {
  fupwindQuad[3] = 0.5*(fr[3]+fr[2]+fr[1]+fr[0]); 
  }
  fupwind[0] = 0.5*(fupwindQuad[3]+fupwindQuad[2]+fupwindQuad[1]+fupwindQuad[0]); 
  fupwind[1] = 0.5*(fupwindQuad[3]-1.0*fupwindQuad[2]+fupwindQuad[1]-1.0*fupwindQuad[0]); 
  fupwind[2] = 0.5*(fupwindQuad[3]+fupwindQuad[2]-1.0*(fupwindQuad[1]+fupwindQuad[0])); 
  fupwind[3] = 0.5*(fupwindQuad[3]-1.0*(fupwindQuad[2]+fupwindQuad[1])+fupwindQuad[0]); 
  double rVal[2];  // rVal=f1/f0 at each control node in dimensions other than y 
  rVal[0] = -(1.0*(1.732050807568877*fupwind[3]-3.0*fupwind[2]))/(3.464101615137754*EPSILON-1.0*fupwind[1]+1.732050807568877*fupwind[0]); 
  rVal[1] = (1.732050807568877*fupwind[3]+3.0*fupwind[2])/(3.464101615137754*EPSILON+fupwind[1]+1.732050807568877*fupwind[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each control node on y surface 
  fqVal[0] = -0.2886751345948129*(fupwind[1]-1.732050807568877*fupwind[0])*limTheta(rVal[0],-1.0,cflR); 
  fqVal[1] = 0.2886751345948129*(fupwind[1]+1.732050807568877*fupwind[0])*limTheta(rVal[1],-1.0,cflR); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*alpha[0]*fhatALVal[0]*dfac_y; 
  incr[1] = 0.5*alpha[0]*fhatALVal[1]*dfac_y; 
  incr[2] = -0.8660254037844386*alpha[0]*fhatALVal[0]*dfac_y; 
  incr[3] = -0.8660254037844386*alpha[0]*fhatALVal[1]*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
#endif
  return std::abs(alpha0); 
} 
