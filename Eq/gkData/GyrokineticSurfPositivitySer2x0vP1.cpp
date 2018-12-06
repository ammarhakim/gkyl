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
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than x 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each control node on x surface 
  // determine upwinding at each surface control node 
  if(0.7071067811865475*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[3]-5.196152422706631*fl[1]))/(6.0*EPSILON-1.732050807568877*fl[2]+3.0*fl[0]); 
  fqVal[0] = -0.1666666666666667*(1.732050807568877*fl[2]-3.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[3]-5.196152422706631*fr[1]))/(6.0*EPSILON-1.732050807568877*fr[2]+3.0*fr[0]); 
  fqVal[0] = -0.1666666666666667*(1.732050807568877*fr[2]-3.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.7071067811865475*alpha[0] > 0) {
  rVal = (3.0*fl[3]+5.196152422706631*fl[1])/(6.0*EPSILON+1.732050807568877*fl[2]+3.0*fl[0]); 
  fqVal[1] = 0.1666666666666667*(1.732050807568877*fl[2]+3.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[3]+5.196152422706631*fr[1])/(6.0*EPSILON+1.732050807568877*fr[2]+3.0*fr[0]); 
  fqVal[1] = 0.1666666666666667*(1.732050807568877*fr[2]+3.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
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
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than y 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each control node on y surface 
  // determine upwinding at each surface control node 
  if(0.7071067811865475*alpha[0] > 0) {
  rVal = -(1.0*(3.0*fl[3]-5.196152422706631*fl[2]))/(6.0*EPSILON-1.732050807568877*fl[1]+3.0*fl[0]); 
  fqVal[0] = -0.1666666666666667*(1.732050807568877*fl[1]-3.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = -(1.0*(3.0*fr[3]-5.196152422706631*fr[2]))/(6.0*EPSILON-1.732050807568877*fr[1]+3.0*fr[0]); 
  fqVal[0] = -0.1666666666666667*(1.732050807568877*fr[1]-3.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
  if(0.7071067811865475*alpha[0] > 0) {
  rVal = (3.0*fl[3]+5.196152422706631*fl[2])/(6.0*EPSILON+1.732050807568877*fl[1]+3.0*fl[0]); 
  fqVal[1] = 0.1666666666666667*(1.732050807568877*fl[1]+3.0*fl[0])*limTheta(rVal,1.0,cflL); 
  } else {
  rVal = (3.0*fr[3]+5.196152422706631*fr[2])/(6.0*EPSILON+1.732050807568877*fr[1]+3.0*fr[0]); 
  fqVal[1] = 0.1666666666666667*(1.732050807568877*fr[1]+3.0*fr[0])*limTheta(rVal,-1.0,cflR); 
  }
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
  return std::abs(alpha0); 
} 
