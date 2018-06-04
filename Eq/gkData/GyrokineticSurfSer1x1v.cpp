#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = wv;

  if (alpha0>0) { 
  incr[0] = 0.5*(1.732050807568877*fl[1]+fl[0])*dfac_x*wv; 
  incr[1] = -0.5*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x*wv; 
  incr[2] = 0.5*(1.732050807568877*fl[3]+fl[2])*dfac_x*wv; 
  incr[3] = -0.5*(3.0*fl[3]+1.732050807568877*fl[2])*dfac_x*wv; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.5*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x*wv; 
  incr[1] = 0.5*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x*wv; 
  incr[2] = -0.5*(1.732050807568877*fr[3]-1.0*fr[2])*dfac_x*wv; 
  incr[3] = 0.5*(3.0*fr[3]-1.732050807568877*fr[2])*dfac_x*wv; 

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
double GyrokineticSurf1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.224744871391589*Phi[1]*dfac_x*q_)/m_;

  if (alpha0>0) { 
  incr[0] = -(0.3535533905932737*Phi[1]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v*dfac_x*q_)/m_; 
  incr[1] = -(0.3535533905932737*Phi[1]*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_v*dfac_x*q_)/m_; 
  incr[2] = (1.060660171779821*Phi[1]*(1.732050807568877*fl[2]+fl[0])*dfac_v*dfac_x*q_)/m_; 
  incr[3] = (1.060660171779821*Phi[1]*(1.732050807568877*fl[3]+fl[1])*dfac_v*dfac_x*q_)/m_; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = (0.3535533905932737*Phi[1]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v*dfac_x*q_)/m_; 
  incr[1] = (0.3535533905932737*Phi[1]*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_v*dfac_x*q_)/m_; 
  incr[2] = -(1.060660171779821*Phi[1]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v*dfac_x*q_)/m_; 
  incr[3] = -(1.060660171779821*Phi[1]*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_v*dfac_x*q_)/m_; 

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
double EmGyrokineticSurf1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = wv;

  if (alpha0>0) { 
  incr[0] = 0.5*(1.732050807568877*fl[1]+fl[0])*dfac_x*wv; 
  incr[1] = -0.5*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x*wv; 
  incr[2] = 0.5*(1.732050807568877*fl[3]+fl[2])*dfac_x*wv; 
  incr[3] = -0.5*(3.0*fl[3]+1.732050807568877*fl[2])*dfac_x*wv; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.5*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x*wv; 
  incr[1] = 0.5*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x*wv; 
  incr[2] = -0.5*(1.732050807568877*fr[3]-1.0*fr[2])*dfac_x*wv; 
  incr[3] = 0.5*(3.0*fr[3]-1.732050807568877*fr[2])*dfac_x*wv; 

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
double EmGyrokineticSurf1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.224744871391589*Phi[1]*dfac_x*q_)/m_;

  if (alpha0>0) { 
  incr[0] = -(0.3535533905932737*Phi[1]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v*dfac_x*q_)/m_; 
  incr[1] = -(0.3535533905932737*Phi[1]*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_v*dfac_x*q_)/m_; 
  incr[2] = (1.060660171779821*Phi[1]*(1.732050807568877*fl[2]+fl[0])*dfac_v*dfac_x*q_)/m_; 
  incr[3] = (1.060660171779821*Phi[1]*(1.732050807568877*fl[3]+fl[1])*dfac_v*dfac_x*q_)/m_; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = (0.3535533905932737*Phi[1]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v*dfac_x*q_)/m_; 
  incr[1] = (0.3535533905932737*Phi[1]*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_v*dfac_x*q_)/m_; 
  incr[2] = -(1.060660171779821*Phi[1]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v*dfac_x*q_)/m_; 
  incr[3] = -(1.060660171779821*Phi[1]*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_v*dfac_x*q_)/m_; 

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
double GyrokineticSurf1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = wv;

  if (alpha0>0) { 
  incr[0] = 0.5*(1.732050807568877*fl[1]+fl[0])*dfac_x*wv; 
  incr[1] = -0.5*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x*wv; 
  incr[2] = 0.5*(1.732050807568877*fl[3]+fl[2])*dfac_x*wv; 
  incr[3] = -0.5*(3.0*fl[3]+1.732050807568877*fl[2])*dfac_x*wv; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.5*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x*wv; 
  incr[1] = 0.5*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x*wv; 
  incr[2] = -0.5*(1.732050807568877*fr[3]-1.0*fr[2])*dfac_x*wv; 
  incr[3] = 0.5*(3.0*fr[3]-1.732050807568877*fr[2])*dfac_x*wv; 

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
double GyrokineticSurf1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.224744871391589*Phi[1]*dfac_x*q_)/m_;

  if (alpha0>0) { 
  incr[0] = -(0.3535533905932737*Phi[1]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v*dfac_x*q_)/m_; 
  incr[1] = -(0.3535533905932737*Phi[1]*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_v*dfac_x*q_)/m_; 
  incr[2] = (1.060660171779821*Phi[1]*(1.732050807568877*fl[2]+fl[0])*dfac_v*dfac_x*q_)/m_; 
  incr[3] = (1.060660171779821*Phi[1]*(1.732050807568877*fl[3]+fl[1])*dfac_v*dfac_x*q_)/m_; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = (0.3535533905932737*Phi[1]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v*dfac_x*q_)/m_; 
  incr[1] = (0.3535533905932737*Phi[1]*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_v*dfac_x*q_)/m_; 
  incr[2] = -(1.060660171779821*Phi[1]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v*dfac_x*q_)/m_; 
  incr[3] = -(1.060660171779821*Phi[1]*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_v*dfac_x*q_)/m_; 

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
double EmGyrokineticSurf1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = wv;

  if (alpha0>0) { 
  incr[0] = 0.5*(1.732050807568877*fl[1]+fl[0])*dfac_x*wv; 
  incr[1] = -0.5*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x*wv; 
  incr[2] = 0.5*(1.732050807568877*fl[3]+fl[2])*dfac_x*wv; 
  incr[3] = -0.5*(3.0*fl[3]+1.732050807568877*fl[2])*dfac_x*wv; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.5*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x*wv; 
  incr[1] = 0.5*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x*wv; 
  incr[2] = -0.5*(1.732050807568877*fr[3]-1.0*fr[2])*dfac_x*wv; 
  incr[3] = 0.5*(3.0*fr[3]-1.732050807568877*fr[2])*dfac_x*wv; 

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
double EmGyrokineticSurf1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(1.224744871391589*Phi[1]*dfac_x*q_)/m_;

  if (alpha0>0) { 
  incr[0] = -(0.3535533905932737*Phi[1]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v*dfac_x*q_)/m_; 
  incr[1] = -(0.3535533905932737*Phi[1]*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_v*dfac_x*q_)/m_; 
  incr[2] = (1.060660171779821*Phi[1]*(1.732050807568877*fl[2]+fl[0])*dfac_v*dfac_x*q_)/m_; 
  incr[3] = (1.060660171779821*Phi[1]*(1.732050807568877*fl[3]+fl[1])*dfac_v*dfac_x*q_)/m_; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = (0.3535533905932737*Phi[1]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v*dfac_x*q_)/m_; 
  incr[1] = (0.3535533905932737*Phi[1]*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_v*dfac_x*q_)/m_; 
  incr[2] = -(1.060660171779821*Phi[1]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v*dfac_x*q_)/m_; 
  incr[3] = -(1.060660171779821*Phi[1]*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_v*dfac_x*q_)/m_; 

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
double dAdtSurf1x1vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
  double dfac_v = 2.0/dxv[1]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.7071067811865475*dApardt[0]*q_)/m_;

  if (alpha0>0) { 
  incr[0] = -(0.3535533905932737*(1.732050807568877*dApardt[1]*fl[3]+1.732050807568877*dApardt[0]*fl[2]+dApardt[1]*fl[1]+dApardt[0]*fl[0])*dfac_v*q_)/m_; 
  incr[1] = -(0.3535533905932737*(1.732050807568877*dApardt[0]*fl[3]+1.732050807568877*dApardt[1]*fl[2]+dApardt[0]*fl[1]+fl[0]*dApardt[1])*dfac_v*q_)/m_; 
  incr[2] = (0.3535533905932737*(3.0*dApardt[1]*fl[3]+3.0*dApardt[0]*fl[2]+1.732050807568877*dApardt[1]*fl[1]+1.732050807568877*dApardt[0]*fl[0])*dfac_v*q_)/m_; 
  incr[3] = (0.3535533905932737*(3.0*dApardt[0]*fl[3]+3.0*dApardt[1]*fl[2]+1.732050807568877*dApardt[0]*fl[1]+1.732050807568877*fl[0]*dApardt[1])*dfac_v*q_)/m_; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = (0.3535533905932737*(1.732050807568877*dApardt[1]*fr[3]+1.732050807568877*dApardt[0]*fr[2]-1.0*dApardt[1]*fr[1]-1.0*dApardt[0]*fr[0])*dfac_v*q_)/m_; 
  incr[1] = (0.3535533905932737*(1.732050807568877*dApardt[0]*fr[3]+1.732050807568877*dApardt[1]*fr[2]-1.0*dApardt[0]*fr[1]-1.0*fr[0]*dApardt[1])*dfac_v*q_)/m_; 
  incr[2] = -(0.3535533905932737*(3.0*dApardt[1]*fr[3]+3.0*dApardt[0]*fr[2]-1.732050807568877*dApardt[1]*fr[1]-1.732050807568877*dApardt[0]*fr[0])*dfac_v*q_)/m_; 
  incr[3] = -(0.3535533905932737*(3.0*dApardt[0]*fr[3]+3.0*dApardt[1]*fr[2]-1.732050807568877*dApardt[0]*fr[1]-1.732050807568877*fr[0]*dApardt[1])*dfac_v*q_)/m_; 

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
