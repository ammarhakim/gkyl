#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[4]; 
  alpha[0] = 1.414213562373095*Gradpar[0]*wv; 
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than x 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each control node on x surface 
  // determine upwinding at each surface control node 
  if(0.5*alpha[0] > 0) {
  rVal = (2.449489742783178*fl[7]-4.242640687119286*(fl[5]+fl[4])+7.348469228349534*fl[1])/(12.0*EPSILON+1.414213562373095*fl[6]-2.449489742783178*(fl[3]+fl[2])+4.242640687119286*fl[0]); 
  fqVal[0] = fmin(0.08333333333333333*(1.414213562373095*fl[6]-2.449489742783178*(fl[3]+fl[2])+4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = (2.449489742783178*fr[7]-4.242640687119286*(fr[5]+fr[4])+7.348469228349534*fr[1])/(12.0*EPSILON+1.414213562373095*fr[6]-2.449489742783178*(fr[3]+fr[2])+4.242640687119286*fr[0]); 
  fqVal[0] = fmin(0.08333333333333333*(1.414213562373095*fr[6]-2.449489742783178*(fr[3]+fr[2])+4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if(0.5*alpha[0] > 0) {
  rVal = -(1.0*(2.449489742783178*fl[7]+4.242640687119286*fl[5]-4.242640687119286*fl[4]-7.348469228349534*fl[1]))/(12.0*EPSILON-1.414213562373095*fl[6]-2.449489742783178*fl[3]+2.449489742783178*fl[2]+4.242640687119286*fl[0]); 
  fqVal[1] = fmin(-0.08333333333333333*(1.414213562373095*fl[6]+2.449489742783178*fl[3]-2.449489742783178*fl[2]-4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = -(1.0*(2.449489742783178*fr[7]+4.242640687119286*fr[5]-4.242640687119286*fr[4]-7.348469228349534*fr[1]))/(12.0*EPSILON-1.414213562373095*fr[6]-2.449489742783178*fr[3]+2.449489742783178*fr[2]+4.242640687119286*fr[0]); 
  fqVal[1] = fmin(-0.08333333333333333*(1.414213562373095*fr[6]+2.449489742783178*fr[3]-2.449489742783178*fr[2]-4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if(0.5*alpha[0] > 0) {
  rVal = -(1.0*(2.449489742783178*fl[7]-4.242640687119286*fl[5]+4.242640687119286*fl[4]-7.348469228349534*fl[1]))/(12.0*EPSILON-1.414213562373095*fl[6]+2.449489742783178*fl[3]-2.449489742783178*fl[2]+4.242640687119286*fl[0]); 
  fqVal[2] = fmin(-0.08333333333333333*(1.414213562373095*fl[6]-2.449489742783178*fl[3]+2.449489742783178*fl[2]-4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = -(1.0*(2.449489742783178*fr[7]-4.242640687119286*fr[5]+4.242640687119286*fr[4]-7.348469228349534*fr[1]))/(12.0*EPSILON-1.414213562373095*fr[6]+2.449489742783178*fr[3]-2.449489742783178*fr[2]+4.242640687119286*fr[0]); 
  fqVal[2] = fmin(-0.08333333333333333*(1.414213562373095*fr[6]-2.449489742783178*fr[3]+2.449489742783178*fr[2]-4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if(0.5*alpha[0] > 0) {
  rVal = (2.449489742783178*fl[7]+4.242640687119286*(fl[5]+fl[4])+7.348469228349534*fl[1])/(12.0*EPSILON+1.414213562373095*fl[6]+2.449489742783178*(fl[3]+fl[2])+4.242640687119286*fl[0]); 
  fqVal[3] = fmin(0.08333333333333333*(1.414213562373095*fl[6]+2.449489742783178*(fl[3]+fl[2])+4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = (2.449489742783178*fr[7]+4.242640687119286*(fr[5]+fr[4])+7.348469228349534*fr[1])/(12.0*EPSILON+1.414213562373095*fr[6]+2.449489742783178*(fr[3]+fr[2])+4.242640687119286*fr[0]); 
  fqVal[3] = fmin(0.08333333333333333*(1.414213562373095*fr[6]+2.449489742783178*(fr[3]+fr[2])+4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[1] = -0.6123724356957944*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*fhatALVal[2]*dfac_x; 
  incr[4] = -0.6123724356957944*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[5] = -0.6123724356957944*alpha[0]*fhatALVal[2]*dfac_x; 
  incr[6] = 0.3535533905932737*alpha[0]*fhatALVal[3]*dfac_x; 
  incr[7] = -0.6123724356957944*alpha[0]*fhatALVal[3]*dfac_x; 

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
double GyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[4]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than vx 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each control node on vx surface 
  // determine upwinding at each surface control node 
  if(0.5*alpha[0] > 0) {
  rVal = (2.449489742783178*fl[7]-4.242640687119286*(fl[6]+fl[4])+7.348469228349534*fl[2])/(12.0*EPSILON+1.414213562373095*fl[5]-2.449489742783178*(fl[3]+fl[1])+4.242640687119286*fl[0]); 
  fqVal[0] = fmin(0.08333333333333333*(1.414213562373095*fl[5]-2.449489742783178*(fl[3]+fl[1])+4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = (2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[4])+7.348469228349534*fr[2])/(12.0*EPSILON+1.414213562373095*fr[5]-2.449489742783178*(fr[3]+fr[1])+4.242640687119286*fr[0]); 
  fqVal[0] = fmin(0.08333333333333333*(1.414213562373095*fr[5]-2.449489742783178*(fr[3]+fr[1])+4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if(0.5*alpha[0] > 0) {
  rVal = -(1.0*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*fl[4]-7.348469228349534*fl[2]))/(12.0*EPSILON-1.414213562373095*fl[5]-2.449489742783178*fl[3]+2.449489742783178*fl[1]+4.242640687119286*fl[0]); 
  fqVal[1] = fmin(-0.08333333333333333*(1.414213562373095*fl[5]+2.449489742783178*fl[3]-2.449489742783178*fl[1]-4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = -(1.0*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*fr[4]-7.348469228349534*fr[2]))/(12.0*EPSILON-1.414213562373095*fr[5]-2.449489742783178*fr[3]+2.449489742783178*fr[1]+4.242640687119286*fr[0]); 
  fqVal[1] = fmin(-0.08333333333333333*(1.414213562373095*fr[5]+2.449489742783178*fr[3]-2.449489742783178*fr[1]-4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if(0.5*alpha[0] > 0) {
  rVal = -(1.0*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*fl[4]-7.348469228349534*fl[2]))/(12.0*EPSILON-1.414213562373095*fl[5]+2.449489742783178*fl[3]-2.449489742783178*fl[1]+4.242640687119286*fl[0]); 
  fqVal[2] = fmin(-0.08333333333333333*(1.414213562373095*fl[5]-2.449489742783178*fl[3]+2.449489742783178*fl[1]-4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = -(1.0*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*fr[4]-7.348469228349534*fr[2]))/(12.0*EPSILON-1.414213562373095*fr[5]+2.449489742783178*fr[3]-2.449489742783178*fr[1]+4.242640687119286*fr[0]); 
  fqVal[2] = fmin(-0.08333333333333333*(1.414213562373095*fr[5]-2.449489742783178*fr[3]+2.449489742783178*fr[1]-4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if(0.5*alpha[0] > 0) {
  rVal = (2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[4])+7.348469228349534*fl[2])/(12.0*EPSILON+1.414213562373095*fl[5]+2.449489742783178*(fl[3]+fl[1])+4.242640687119286*fl[0]); 
  fqVal[3] = fmin(0.08333333333333333*(1.414213562373095*fl[5]+2.449489742783178*(fl[3]+fl[1])+4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = (2.449489742783178*fr[7]+4.242640687119286*(fr[6]+fr[4])+7.348469228349534*fr[2])/(12.0*EPSILON+1.414213562373095*fr[5]+2.449489742783178*(fr[3]+fr[1])+4.242640687119286*fr[0]); 
  fqVal[3] = fmin(0.08333333333333333*(1.414213562373095*fr[5]+2.449489742783178*(fr[3]+fr[1])+4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*alpha[0]*fhatALVal[0]*dfac_v; 
  incr[1] = 0.3535533905932737*alpha[0]*fhatALVal[1]*dfac_v; 
  incr[2] = -0.6123724356957944*alpha[0]*fhatALVal[0]*dfac_v; 
  incr[3] = 0.3535533905932737*alpha[0]*fhatALVal[2]*dfac_v; 
  incr[4] = -0.6123724356957944*alpha[0]*fhatALVal[1]*dfac_v; 
  incr[5] = 0.3535533905932737*alpha[0]*fhatALVal[3]*dfac_v; 
  incr[6] = -0.6123724356957944*alpha[0]*fhatALVal[2]*dfac_v; 
  incr[7] = -0.6123724356957944*alpha[0]*fhatALVal[3]*dfac_v; 

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
double GyrokineticSurfPositivity1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.3535533905932737*(1.732050807568877*Gradpar[1]-1.0*Gradpar[0])*wv; 

  double alpha[4]; 
  alpha[0] = (1.414213562373095*Gradpar[0]-2.449489742783178*Gradpar[1])*wv; 
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than x 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each control node on x surface 
  // determine upwinding at each surface control node 
  if(0.5*alpha[0] > 0) {
  rVal = (2.449489742783178*fl[7]-4.242640687119286*(fl[5]+fl[4])+7.348469228349534*fl[1])/(12.0*EPSILON+1.414213562373095*fl[6]-2.449489742783178*(fl[3]+fl[2])+4.242640687119286*fl[0]); 
  fqVal[0] = fmin(0.08333333333333333*(1.414213562373095*fl[6]-2.449489742783178*(fl[3]+fl[2])+4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = (2.449489742783178*fr[7]-4.242640687119286*(fr[5]+fr[4])+7.348469228349534*fr[1])/(12.0*EPSILON+1.414213562373095*fr[6]-2.449489742783178*(fr[3]+fr[2])+4.242640687119286*fr[0]); 
  fqVal[0] = fmin(0.08333333333333333*(1.414213562373095*fr[6]-2.449489742783178*(fr[3]+fr[2])+4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if(0.5*alpha[0] > 0) {
  rVal = -(1.0*(2.449489742783178*fl[7]+4.242640687119286*fl[5]-4.242640687119286*fl[4]-7.348469228349534*fl[1]))/(12.0*EPSILON-1.414213562373095*fl[6]-2.449489742783178*fl[3]+2.449489742783178*fl[2]+4.242640687119286*fl[0]); 
  fqVal[1] = fmin(-0.08333333333333333*(1.414213562373095*fl[6]+2.449489742783178*fl[3]-2.449489742783178*fl[2]-4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = -(1.0*(2.449489742783178*fr[7]+4.242640687119286*fr[5]-4.242640687119286*fr[4]-7.348469228349534*fr[1]))/(12.0*EPSILON-1.414213562373095*fr[6]-2.449489742783178*fr[3]+2.449489742783178*fr[2]+4.242640687119286*fr[0]); 
  fqVal[1] = fmin(-0.08333333333333333*(1.414213562373095*fr[6]+2.449489742783178*fr[3]-2.449489742783178*fr[2]-4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if(0.5*alpha[0] > 0) {
  rVal = -(1.0*(2.449489742783178*fl[7]-4.242640687119286*fl[5]+4.242640687119286*fl[4]-7.348469228349534*fl[1]))/(12.0*EPSILON-1.414213562373095*fl[6]+2.449489742783178*fl[3]-2.449489742783178*fl[2]+4.242640687119286*fl[0]); 
  fqVal[2] = fmin(-0.08333333333333333*(1.414213562373095*fl[6]-2.449489742783178*fl[3]+2.449489742783178*fl[2]-4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = -(1.0*(2.449489742783178*fr[7]-4.242640687119286*fr[5]+4.242640687119286*fr[4]-7.348469228349534*fr[1]))/(12.0*EPSILON-1.414213562373095*fr[6]+2.449489742783178*fr[3]-2.449489742783178*fr[2]+4.242640687119286*fr[0]); 
  fqVal[2] = fmin(-0.08333333333333333*(1.414213562373095*fr[6]-2.449489742783178*fr[3]+2.449489742783178*fr[2]-4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if(0.5*alpha[0] > 0) {
  rVal = (2.449489742783178*fl[7]+4.242640687119286*(fl[5]+fl[4])+7.348469228349534*fl[1])/(12.0*EPSILON+1.414213562373095*fl[6]+2.449489742783178*(fl[3]+fl[2])+4.242640687119286*fl[0]); 
  fqVal[3] = fmin(0.08333333333333333*(1.414213562373095*fl[6]+2.449489742783178*(fl[3]+fl[2])+4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = (2.449489742783178*fr[7]+4.242640687119286*(fr[5]+fr[4])+7.348469228349534*fr[1])/(12.0*EPSILON+1.414213562373095*fr[6]+2.449489742783178*(fr[3]+fr[2])+4.242640687119286*fr[0]); 
  fqVal[3] = fmin(0.08333333333333333*(1.414213562373095*fr[6]+2.449489742783178*(fr[3]+fr[2])+4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[1] = -0.6123724356957944*alpha[0]*fhatALVal[0]*dfac_x; 
  incr[2] = 0.3535533905932737*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[3] = 0.3535533905932737*alpha[0]*fhatALVal[2]*dfac_x; 
  incr[4] = -0.6123724356957944*alpha[0]*fhatALVal[1]*dfac_x; 
  incr[5] = -0.6123724356957944*alpha[0]*fhatALVal[2]*dfac_x; 
  incr[6] = 0.3535533905932737*alpha[0]*fhatALVal[3]*dfac_x; 
  incr[7] = -0.6123724356957944*alpha[0]*fhatALVal[3]*dfac_x; 

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
double GyrokineticSurfPositivity1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cflL, const double cflR, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*Gradpar[0]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 

  double alpha[4]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 
  alpha[1] = -(1.732050807568877*Gradpar[1]*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 
  alpha[2] = -(1.0*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[3] = -(1.0*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  double rVal;  // rVal=f1/f0 at each control node in dimensions other than vx 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each control node on vx surface 
  // determine upwinding at each surface control node 
  if(0.1666666666666667*alpha[3]-0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) {
  rVal = (2.449489742783178*fl[7]-4.242640687119286*(fl[6]+fl[4])+7.348469228349534*fl[2])/(12.0*EPSILON+1.414213562373095*fl[5]-2.449489742783178*(fl[3]+fl[1])+4.242640687119286*fl[0]); 
  fqVal[0] = fmin(0.08333333333333333*(1.414213562373095*fl[5]-2.449489742783178*(fl[3]+fl[1])+4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = (2.449489742783178*fr[7]-4.242640687119286*(fr[6]+fr[4])+7.348469228349534*fr[2])/(12.0*EPSILON+1.414213562373095*fr[5]-2.449489742783178*(fr[3]+fr[1])+4.242640687119286*fr[0]); 
  fqVal[0] = fmin(0.08333333333333333*(1.414213562373095*fr[5]-2.449489742783178*(fr[3]+fr[1])+4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if((-0.1666666666666667*alpha[3])-0.2886751345948129*alpha[2]+0.2886751345948129*alpha[1]+0.5*alpha[0] > 0) {
  rVal = -(1.0*(2.449489742783178*fl[7]+4.242640687119286*fl[6]-4.242640687119286*fl[4]-7.348469228349534*fl[2]))/(12.0*EPSILON-1.414213562373095*fl[5]-2.449489742783178*fl[3]+2.449489742783178*fl[1]+4.242640687119286*fl[0]); 
  fqVal[1] = fmin(-0.08333333333333333*(1.414213562373095*fl[5]+2.449489742783178*fl[3]-2.449489742783178*fl[1]-4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = -(1.0*(2.449489742783178*fr[7]+4.242640687119286*fr[6]-4.242640687119286*fr[4]-7.348469228349534*fr[2]))/(12.0*EPSILON-1.414213562373095*fr[5]-2.449489742783178*fr[3]+2.449489742783178*fr[1]+4.242640687119286*fr[0]); 
  fqVal[1] = fmin(-0.08333333333333333*(1.414213562373095*fr[5]+2.449489742783178*fr[3]-2.449489742783178*fr[1]-4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if((-0.1666666666666667*alpha[3])+0.2886751345948129*alpha[2]-0.2886751345948129*alpha[1]+0.5*alpha[0] > 0) {
  rVal = -(1.0*(2.449489742783178*fl[7]-4.242640687119286*fl[6]+4.242640687119286*fl[4]-7.348469228349534*fl[2]))/(12.0*EPSILON-1.414213562373095*fl[5]+2.449489742783178*fl[3]-2.449489742783178*fl[1]+4.242640687119286*fl[0]); 
  fqVal[2] = fmin(-0.08333333333333333*(1.414213562373095*fl[5]-2.449489742783178*fl[3]+2.449489742783178*fl[1]-4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = -(1.0*(2.449489742783178*fr[7]-4.242640687119286*fr[6]+4.242640687119286*fr[4]-7.348469228349534*fr[2]))/(12.0*EPSILON-1.414213562373095*fr[5]+2.449489742783178*fr[3]-2.449489742783178*fr[1]+4.242640687119286*fr[0]); 
  fqVal[2] = fmin(-0.08333333333333333*(1.414213562373095*fr[5]-2.449489742783178*fr[3]+2.449489742783178*fr[1]-4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  if(0.1666666666666667*alpha[3]+0.2886751345948129*(alpha[2]+alpha[1])+0.5*alpha[0] > 0) {
  rVal = (2.449489742783178*fl[7]+4.242640687119286*(fl[6]+fl[4])+7.348469228349534*fl[2])/(12.0*EPSILON+1.414213562373095*fl[5]+2.449489742783178*(fl[3]+fl[1])+4.242640687119286*fl[0]); 
  fqVal[3] = fmin(0.08333333333333333*(1.414213562373095*fl[5]+2.449489742783178*(fl[3]+fl[1])+4.242640687119286*fl[0])*limTheta(rVal,1.0,cflL), fl[0]*cflL); 
  } else {
  rVal = (2.449489742783178*fr[7]+4.242640687119286*(fr[6]+fr[4])+7.348469228349534*fr[2])/(12.0*EPSILON+1.414213562373095*fr[5]+2.449489742783178*(fr[3]+fr[1])+4.242640687119286*fr[0]); 
  fqVal[3] = fmin(0.08333333333333333*(1.414213562373095*fr[5]+2.449489742783178*(fr[3]+fr[1])+4.242640687119286*fr[0])*limTheta(rVal,-1.0,cflR), fr[0]*cflR); 
  }
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.3535533905932737*(alpha[3]*fhatALVal[3]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.3535533905932737*(alpha[2]*fhatALVal[3]+fhatALVal[2]*alpha[3]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.6123724356957944*(alpha[3]*fhatALVal[3]+alpha[2]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.3535533905932737*(alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_v; 
  incr[4] = -0.6123724356957944*(alpha[2]*fhatALVal[3]+fhatALVal[2]*alpha[3]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[5] = 0.3535533905932737*(alpha[0]*fhatALVal[3]+fhatALVal[0]*alpha[3]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_v; 
  incr[6] = -0.6123724356957944*(alpha[1]*fhatALVal[3]+fhatALVal[1]*alpha[3]+alpha[0]*fhatALVal[2]+fhatALVal[0]*alpha[2])*dfac_v; 
  incr[7] = -0.6123724356957944*(alpha[0]*fhatALVal[3]+fhatALVal[0]*alpha[3]+alpha[1]*fhatALVal[2]+fhatALVal[1]*alpha[2])*dfac_v; 

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
