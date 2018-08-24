#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = 0.5*wv; 

  double alpha[8]; 
  alpha[0] = 2.828427124746191*wv; 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_x; 

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
  } else { 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.6123724356957944*Phi[1]*dfac_x*q_)/m_; 

  double alpha[8]; 
  alpha[0] = -(3.464101615137754*Phi[1]*dfac_x*q_)/m_; 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[4]+fl[1])*dfac_v; 
  incr[2] = -0.1767766952966368*alpha[0]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[6]+fl[3])*dfac_v; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[1])*dfac_v; 
  incr[5] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[7]+fl[5])*dfac_v; 
  incr[6] = -0.1767766952966368*alpha[0]*(3.0*fl[6]+1.732050807568877*fl[3])*dfac_v; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[5])*dfac_v; 

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
  } else { 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  incr[1] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[1])*dfac_v; 
  incr[2] = 0.1767766952966368*alpha[0]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[6]-1.0*fr[3])*dfac_v; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[1])*dfac_v; 
  incr[5] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[5])*dfac_v; 
  incr[6] = 0.1767766952966368*alpha[0]*(3.0*fr[6]-1.732050807568877*fr[3])*dfac_v; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[5])*dfac_v; 

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
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = 0.5*wv; 

  double alpha[8]; 
  alpha[0] = 2.828427124746191*wv; 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[4]+fl[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[5]+fl[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(3.0*fl[4]+1.732050807568877*fl[2])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[3])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[7]+fl[6])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[6])*dfac_x; 

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
  } else { 
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[4]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(3.0*fr[4]-1.732050807568877*fr[2])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[3])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[6])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[6])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.25*(2.449489742783178*Bmag[1]*dfac_x*wm+2.449489742783178*Phi[1]*dfac_x*q_))/m_; 

  double alpha[8]; 
  alpha[0] = (-(3.464101615137754*Bmag[1]*dfac_x*wm)/m_)-(3.464101615137754*Phi[1]*dfac_x*q_)/m_; 
  alpha[3] = -(2.0*Bmag[1]*dfac_x)/(dfac_m*m_); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*alpha[3]*fl[6]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*alpha[3]*fl[7]+alpha[3]*fl[5]+1.732050807568877*alpha[0]*fl[4]+alpha[0]*fl[1])*dfac_v; 
  incr[2] = -0.1767766952966368*(3.0*alpha[3]*fl[6]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[3] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[6]+alpha[0]*fl[3]+1.732050807568877*fl[2]*alpha[3]+fl[0]*alpha[3])*dfac_v; 
  incr[4] = -0.1767766952966368*(3.0*alpha[3]*fl[7]+1.732050807568877*alpha[3]*fl[5]+3.0*alpha[0]*fl[4]+1.732050807568877*alpha[0]*fl[1])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*alpha[0]*fl[7]+alpha[0]*fl[5]+1.732050807568877*alpha[3]*fl[4]+fl[1]*alpha[3])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*alpha[0]*fl[6]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*alpha[0]*fl[7]+1.732050807568877*alpha[0]*fl[5]+3.0*alpha[3]*fl[4]+1.732050807568877*fl[1]*alpha[3])*dfac_v; 

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
  } else { 
  incr[0] = -0.1767766952966368*(1.732050807568877*alpha[3]*fr[6]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*alpha[3]*fr[7]-1.0*alpha[3]*fr[5]+1.732050807568877*alpha[0]*fr[4]-1.0*alpha[0]*fr[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(3.0*alpha[3]*fr[6]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[3] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[6]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_v; 
  incr[4] = 0.1767766952966368*(3.0*alpha[3]*fr[7]-1.732050807568877*alpha[3]*fr[5]+3.0*alpha[0]*fr[4]-1.732050807568877*alpha[0]*fr[1])*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*alpha[0]*fr[7]-1.0*alpha[0]*fr[5]+1.732050807568877*alpha[3]*fr[4]-1.0*fr[1]*alpha[3])*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*alpha[0]*fr[6]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*alpha[0]*fr[7]-1.732050807568877*alpha[0]*fr[5]+3.0*alpha[3]*fr[4]-1.732050807568877*fr[1]*alpha[3])*dfac_v; 

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
  } 
  return std::abs(alpha0); 
} 
