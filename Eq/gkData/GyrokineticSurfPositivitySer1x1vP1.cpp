#include <GyrokineticModDecl.h> 
double GyrokineticSurfPositivity1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*Gradpar[0]*wv; 

  double alpha[2]; 
  alpha[0] = Gradpar[0]*wv; 
  if (alpha0>0) { 
  double rVal[2];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cfl); 
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
  double rVal[2];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cfl); 
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
  return std::abs(alpha0); 
} 
double GyrokineticSurfPositivity1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.4330127018922193*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[2]; 
  alpha[0] = -(1.224744871391589*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  if (alpha0>0) { 
  double rVal[2];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cfl); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*alpha[0]*fhatALVal[0]*dfac_v; 
  incr[1] = 0.5*alpha[0]*fhatALVal[1]*dfac_v; 
  incr[2] = -0.8660254037844386*alpha[0]*fhatALVal[0]*dfac_v; 
  incr[3] = -0.8660254037844386*alpha[0]*fhatALVal[1]*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[2];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cfl); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*alpha[0]*fhatALVal[0]*dfac_v; 
  incr[1] = 0.5*alpha[0]*fhatALVal[1]*dfac_v; 
  incr[2] = -0.8660254037844386*alpha[0]*fhatALVal[0]*dfac_v; 
  incr[3] = -0.8660254037844386*alpha[0]*fhatALVal[1]*dfac_v; 

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
double GyrokineticSurfPositivity1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.1767766952966368*wv*(Bmag[1]*((5.196152422706631*BmagInv[1]-3.0*BmagInv[0])*geoY[1]+geoY[0]*(1.732050807568877*BmagInv[0]-3.0*BmagInv[1]))*dfac_x*m_*wv+2.0*(1.732050807568877*Gradpar[1]-1.0*Gradpar[0])*q_))/q_; 

  double alpha[2]; 
  alpha[0] = -(0.5*(Bmag[1]*((5.196152422706631*BmagInv[1]-3.0*BmagInv[0])*geoY[1]+geoY[0]*(1.732050807568877*BmagInv[0]-3.0*BmagInv[1]))*dfac_x*m_*wv2+(3.464101615137754*Gradpar[1]-2.0*Gradpar[0])*q_*wv))/q_; 
  alpha[1] = -(0.5*Bmag[1]*((3.0*BmagInv[1]-1.732050807568877*BmagInv[0])*geoY[1]+geoY[0]*(BmagInv[0]-1.732050807568877*BmagInv[1]))*dfac_x*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  double rVal[2];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[1]))/(3.464101615137754*EPSILON-1.0*fl[2]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[1])/(3.464101615137754*EPSILON+fl[2]+1.732050807568877*fl[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fl[2]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = 0.2886751345948129*(fl[2]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cfl); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = -0.8660254037844386*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = 0.5*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 
  incr[3] = -0.8660254037844386*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[2];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[1]))/(3.464101615137754*EPSILON-1.0*fr[2]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[1])/(3.464101615137754*EPSILON+fr[2]+1.732050807568877*fr[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fr[2]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = 0.2886751345948129*(fr[2]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cfl); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[1] = -0.8660254037844386*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_x; 
  incr[2] = 0.5*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 
  incr[3] = -0.8660254037844386*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_x; 

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
double GyrokineticSurfPositivity1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*Phi[1]*dfac_x*(3.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_v*dfac_x*m_*wv-1.0*(3.464101615137754*Gradpar[0]*dfac_v*q_+3.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x*m_)))/(dfac_v*m_); 

  double alpha[2]; 
  alpha[0] = (0.25*Phi[1]*(dfac_v*(4.242640687119286*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_*wv-4.898979485566357*Gradpar[0]*dfac_x*q_)-4.242640687119286*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_))/(dfac_v*m_); 
  alpha[1] = (0.25*Phi[1]*(dfac_v*(4.242640687119286*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_*wv-4.898979485566357*Gradpar[1]*dfac_x*q_)-4.242640687119286*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_))/(dfac_v*m_); 
  if (alpha0>0) { 
  double rVal[2];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fl[3]-3.0*fl[2]))/(3.464101615137754*EPSILON-1.0*fl[1]+1.732050807568877*fl[0]); 
  rVal[1] = (1.732050807568877*fl[3]+3.0*fl[2])/(3.464101615137754*EPSILON+fl[1]+1.732050807568877*fl[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fl[1]-1.732050807568877*fl[0])*limTheta(rVal[0],1.0,cfl); 
  fqVal[1] = 0.2886751345948129*(fl[1]+1.732050807568877*fl[0])*limTheta(rVal[1],1.0,cfl); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.5*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.8660254037844386*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = -0.8660254037844386*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  double rVal[2];  // rVal=f1/f0 at each node 
  rVal[0] = -(1.0*(1.732050807568877*fr[3]-3.0*fr[2]))/(3.464101615137754*EPSILON-1.0*fr[1]+1.732050807568877*fr[0]); 
  rVal[1] = (1.732050807568877*fr[3]+3.0*fr[2])/(3.464101615137754*EPSILON+fr[1]+1.732050807568877*fr[0]); 
  double fqVal[2];  // fqVal = anti-limited f evaluated at each node 
  fqVal[0] = -0.2886751345948129*(fr[1]-1.732050807568877*fr[0])*limTheta(rVal[0],-1.0,cfl); 
  fqVal[1] = 0.2886751345948129*(fr[1]+1.732050807568877*fr[0])*limTheta(rVal[1],-1.0,cfl); 
  double fhatALVal[2];  // fhatALVal = mode coefficients of anti-limited f 
  fhatALVal[0] = 0.7071067811865475*(fqVal[1]+fqVal[0]); 
  fhatALVal[1] = 1.224744871391589*(fqVal[1]-1.0*fqVal[0]); 
  incr[0] = 0.5*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[1] = 0.5*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 
  incr[2] = -0.8660254037844386*(alpha[1]*fhatALVal[1]+alpha[0]*fhatALVal[0])*dfac_v; 
  incr[3] = -0.8660254037844386*(alpha[0]*fhatALVal[1]+fhatALVal[0]*alpha[1])*dfac_v; 

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
