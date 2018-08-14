#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.7071067811865475*Gradpar[0]*wv; 

  double alpha[4]; 
  alpha[0] = 1.414213562373095*Gradpar[0]*wv; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[2] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[1]+fhatALVal[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fhatALVal[1]+1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[3]+fhatALVal[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(3.0*fhatALVal[3]+1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[2] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[1]-1.0*fhatALVal[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fhatALVal[1]-1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[3]-1.0*fhatALVal[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(3.0*fhatALVal[3]-1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.8660254037844386*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[4]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[2]+fhatALVal[0])*dfac_v; 
  incr[1] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[3]+fhatALVal[1])*dfac_v; 
  incr[2] = -0.25*alpha[0]*(3.0*fhatALVal[2]+1.732050807568877*fhatALVal[0])*dfac_v; 
  incr[3] = -0.25*alpha[0]*(3.0*fhatALVal[3]+1.732050807568877*fhatALVal[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[2]-1.0*fhatALVal[0])*dfac_v; 
  incr[1] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[3]-1.0*fhatALVal[1])*dfac_v; 
  incr[2] = 0.25*alpha[0]*(3.0*fhatALVal[2]-1.732050807568877*fhatALVal[0])*dfac_v; 
  incr[3] = 0.25*alpha[0]*(3.0*fhatALVal[3]-1.732050807568877*fhatALVal[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*(2.449489742783178*geoY[0]*Apar[1]*dfac_x+2.0*Gradpar[0])*wv; 

  double alpha[4]; 
  alpha[0] = 1.732050807568877*geoY[0]*Apar[1]*dfac_x*wv+1.414213562373095*Gradpar[0]*wv; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[2] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[1]+fhatALVal[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fhatALVal[1]+1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fhatALVal[3]+fhatALVal[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(3.0*fhatALVal[3]+1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[2] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[1]-1.0*fhatALVal[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fhatALVal[1]-1.732050807568877*fhatALVal[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fhatALVal[3]-1.0*fhatALVal[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(3.0*fhatALVal[3]-1.732050807568877*fhatALVal[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(4.242640687119286*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x+3.464101615137754*Gradpar[0]*Phi[1]*dfac_x+2.828427124746191*dApardt[0])*q_)/m_; 

  double alpha[4]; 
  alpha[0] = (-(2.121320343559643*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x*q_)/m_)-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*(1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*alpha[1]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]+1.732050807568877*alpha[1]*fhatALVal[1]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = -0.25*(3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]+1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*(1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]-1.0*alpha[1]*fhatALVal[1]-1.0*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]-1.0*alpha[0]*fhatALVal[1]-1.0*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(3.0*alpha[1]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]-1.732050807568877*alpha[1]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.25*(3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[1]-1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(((7.348469228349534*Bmag[1]*BmagInv[1]-4.242640687119286*BmagInv[0]*Bmag[1])*geoY[1]-4.242640687119286*geoY[0]*Bmag[1]*BmagInv[1]+2.449489742783178*BmagInv[0]*geoY[0]*Bmag[1])*dfac_x*m_*wv2+(4.898979485566357*Gradpar[1]-2.828427124746191*Gradpar[0])*q_*wv))/q_; 

  double alpha[4]; 
  alpha[0] = (-(3.674234614174767*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv2)/q_)+(2.121320343559643*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv2)/q_+(2.121320343559643*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv2)/q_-(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv2)/q_-2.449489742783178*Gradpar[1]*wv+1.414213562373095*Gradpar[0]*wv; 
  alpha[2] = (-(2.121320343559642*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_))+(1.224744871391589*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(1.224744871391589*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv)/(dfac_v*q_)-(0.7071067811865475*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[2] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*(1.732050807568877*alpha[2]*fhatALVal[3]+alpha[2]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = -0.25*(3.0*alpha[2]*fhatALVal[3]+1.732050807568877*alpha[2]*fhatALVal[2]+3.0*alpha[0]*fhatALVal[1]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+alpha[0]*fhatALVal[2]+1.732050807568877*fhatALVal[1]*alpha[2]+fhatALVal[0]*alpha[2])*dfac_x; 
  incr[3] = -0.25*(3.0*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]+3.0*fhatALVal[1]*alpha[2]+1.732050807568877*fhatALVal[0]*alpha[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[2] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*(1.732050807568877*alpha[2]*fhatALVal[3]-1.0*alpha[2]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]-1.0*alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = 0.25*(3.0*alpha[2]*fhatALVal[3]-1.732050807568877*alpha[2]*fhatALVal[2]+3.0*alpha[0]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = -0.25*(1.732050807568877*alpha[0]*fhatALVal[3]-1.0*alpha[0]*fhatALVal[2]+1.732050807568877*fhatALVal[1]*alpha[2]-1.0*fhatALVal[0]*alpha[2])*dfac_x; 
  incr[3] = 0.25*(3.0*alpha[0]*fhatALVal[3]-1.732050807568877*alpha[0]*fhatALVal[2]+3.0*fhatALVal[1]*alpha[2]-1.732050807568877*fhatALVal[0]*alpha[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.25*((3.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]+3.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_v*dfac_x*dfac_x*m_*wv-3.464101615137754*Gradpar[0]*Phi[1]*dfac_v*dfac_x*q_+((-3.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1])-3.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_x*dfac_x*m_))/(dfac_v*m_); 

  double alpha[4]; 
  alpha[0] = 1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*wv-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[1] = 1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*wv-(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x*q_)/m_-(1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*(1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*alpha[1]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]+1.732050807568877*alpha[1]*fhatALVal[1]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = -0.25*(3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]+1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*(1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]-1.0*alpha[1]*fhatALVal[1]-1.0*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]-1.0*alpha[0]*fhatALVal[1]-1.0*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(3.0*alpha[1]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]-1.732050807568877*alpha[1]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.25*(3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[1]-1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(((7.348469228349534*Bmag[1]*BmagInv[1]-4.242640687119286*BmagInv[0]*Bmag[1])*geoY[1]-4.242640687119286*geoY[0]*Bmag[1]*BmagInv[1]+2.449489742783178*BmagInv[0]*geoY[0]*Bmag[1])*dfac_x*m_*wv2+((((5.196152422706631*Apar[0]-9.0*Apar[1])*Bmag[1]*BmagInv[1]+(5.196152422706631*BmagInv[0]*Apar[1]-3.0*Apar[0]*BmagInv[0])*Bmag[1]+6.0*Apar[1])*geoY[1]+(5.196152422706631*geoY[0]*Apar[1]-3.0*Apar[0]*geoY[0])*Bmag[1]*BmagInv[1]+(1.732050807568877*Apar[0]*BmagInv[0]*geoY[0]-3.0*BmagInv[0]*geoY[0]*Apar[1])*Bmag[1]-3.464101615137754*geoY[0]*Apar[1])*dfac_x+4.898979485566357*Gradpar[1]-2.828427124746191*Gradpar[0])*q_*wv))/q_; 

  double alpha[4]; 
  alpha[0] = (-(3.674234614174767*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv2)/q_)+(2.121320343559643*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv2)/q_+(2.121320343559643*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv2)/q_-(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv2)/q_+4.5*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*wv-2.598076211353316*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*wv-2.598076211353316*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*dfac_x*wv+1.5*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*wv-3.0*Apar[1]*geoY[1]*dfac_x*wv-2.598076211353316*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*dfac_x*wv+1.5*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*wv+1.5*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*dfac_x*wv-0.8660254037844386*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*wv+1.732050807568877*geoY[0]*Apar[1]*dfac_x*wv-2.449489742783178*Gradpar[1]*wv+1.414213562373095*Gradpar[0]*wv; 
  alpha[2] = (-(2.121320343559642*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_))+(1.224744871391589*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(1.224744871391589*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv)/(dfac_v*q_)-(0.7071067811865475*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[2] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*(1.732050807568877*alpha[2]*fhatALVal[3]+alpha[2]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = -0.25*(3.0*alpha[2]*fhatALVal[3]+1.732050807568877*alpha[2]*fhatALVal[2]+3.0*alpha[0]*fhatALVal[1]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = 0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+alpha[0]*fhatALVal[2]+1.732050807568877*fhatALVal[1]*alpha[2]+fhatALVal[0]*alpha[2])*dfac_x; 
  incr[3] = -0.25*(3.0*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]+3.0*fhatALVal[1]*alpha[2]+1.732050807568877*fhatALVal[0]*alpha[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[2] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[0]); 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*(1.732050807568877*alpha[2]*fhatALVal[3]-1.0*alpha[2]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]-1.0*alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = 0.25*(3.0*alpha[2]*fhatALVal[3]-1.732050807568877*alpha[2]*fhatALVal[2]+3.0*alpha[0]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = -0.25*(1.732050807568877*alpha[0]*fhatALVal[3]-1.0*alpha[0]*fhatALVal[2]+1.732050807568877*fhatALVal[1]*alpha[2]-1.0*fhatALVal[0]*alpha[2])*dfac_x; 
  incr[3] = 0.25*(3.0*alpha[0]*fhatALVal[3]-1.732050807568877*alpha[0]*fhatALVal[2]+3.0*fhatALVal[1]*alpha[2]-1.732050807568877*fhatALVal[0]*alpha[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurfPositivity1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*((6.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]+6.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_v*dfac_x*dfac_x*m_*wv+(((4.242640687119286*Apar[0]*Bmag[1]*BmagInv[1]+4.242640687119286*BmagInv[0]*Apar[1]*Bmag[1])*Phi[1]*geoY[1]+(4.242640687119286*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]+4.242640687119286*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]-8.485281374238571*geoY[0]*Apar[1])*Phi[1])*dfac_v*dfac_x*dfac_x-6.928203230275509*Gradpar[0]*Phi[1]*dfac_v*dfac_x-5.656854249492382*dApardt[0]*dfac_v)*q_+((-6.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1])-6.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_x*dfac_x*m_))/(dfac_v*m_); 

  double alpha[4]; 
  alpha[0] = 1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*wv+(1.060660171779821*Apar[0]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*BmagInv[0]*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(2.121320343559643*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(1.414213562373096*dApardt[0]*q_)/m_-(1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[1] = 1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*wv+(1.909188309203679*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_-(2.121320343559643*Apar[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x*q_)/m_-(1.414213562373096*dApardt[1]*q_)/m_-(1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  if (alpha0>0) { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[3] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[2],1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[3],1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = 0.25*(1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]+alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]+alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*alpha[1]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]+1.732050807568877*alpha[1]*fhatALVal[1]+1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = -0.25*(3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]+1.732050807568877*alpha[0]*fhatALVal[1]+1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[4];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  rVal[2] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[3] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[4];  // fqVal = anti-limited f evaluated at each node 
  double cflVal = std::abs((alpha0*dt)/dxv[1]); 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cflVal); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cflVal); 
  fqVal[2] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[2],-1.0,cflVal); 
  fqVal[3] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[3],-1.0,cflVal); 
  double fhatALVal[4];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.5*(fqVal[3]+fqVal[2]+fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 0.8660254037844386*(fqVal[3]-1.0*fqVal[2]+fqVal[1]-1.0*fqVal[0]); 
  fhatALVal[2] = 0.8660254037844386*(fqVal[3]+fqVal[2]-1.0*(fqVal[1]+fqVal[0])); 
  fhatALVal[3] = 1.5*(fqVal[3]-1.0*(fqVal[2]+fqVal[1])+fqVal[0]); 
  incr[0] = -0.25*(1.732050807568877*alpha[1]*fhatALVal[3]+1.732050807568877*alpha[0]*fhatALVal[2]-1.0*alpha[1]*fhatALVal[1]-1.0*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*alpha[0]*fhatALVal[3]+1.732050807568877*alpha[1]*fhatALVal[2]-1.0*alpha[0]*fhatALVal[1]-1.0*fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(3.0*alpha[1]*fhatALVal[3]+3.0*alpha[0]*fhatALVal[2]-1.732050807568877*alpha[1]*fhatALVal[1]-1.732050807568877*alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = 0.25*(3.0*alpha[0]*fhatALVal[3]+3.0*alpha[1]*fhatALVal[2]-1.732050807568877*alpha[0]*fhatALVal[1]-1.732050807568877*fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } 
  return std::abs(alpha0); 
} 
