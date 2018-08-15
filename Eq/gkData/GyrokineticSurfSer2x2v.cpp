#include <GyrokineticModDecl.h> 
double GyrokineticSurf2x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.5*(2.0*BdriftX[0]*m_*wv2+(3.0*BmagInv[0]*Phi[3]-1.732050807568877*BmagInv[0]*Phi[2])*dfac_y*q_))/q_; 

  double alpha[16]; 
  alpha[0] = (2.0*BdriftX[0]*m_*wv2)/q_+3.0*BmagInv[0]*Phi[3]*dfac_y-1.732050807568877*BmagInv[0]*Phi[2]*dfac_y; 
  alpha[3] = (1.154700538379252*BdriftX[0]*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[3]*fl[6]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[1]+alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.125*(3.0*alpha[3]*fl[6]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.125*(1.732050807568877*alpha[3]*fl[11]+alpha[3]*fl[7]+1.732050807568877*alpha[0]*fl[5]+alpha[0]*fl[2])*dfac_x; 
  incr[3] = 0.125*(1.732050807568877*alpha[0]*fl[6]+alpha[0]*fl[3]+1.732050807568877*fl[1]*alpha[3]+fl[0]*alpha[3])*dfac_x; 
  incr[4] = 0.125*(1.732050807568877*alpha[3]*fl[13]+alpha[3]*fl[10]+1.732050807568877*alpha[0]*fl[8]+alpha[0]*fl[4])*dfac_x; 
  incr[5] = -0.125*(3.0*alpha[3]*fl[11]+1.732050807568877*alpha[3]*fl[7]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[2])*dfac_x; 
  incr[6] = -0.125*(3.0*alpha[0]*fl[6]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[1]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_x; 
  incr[7] = 0.125*(1.732050807568877*alpha[0]*fl[11]+alpha[0]*fl[7]+1.732050807568877*alpha[3]*fl[5]+fl[2]*alpha[3])*dfac_x; 
  incr[8] = -0.125*(3.0*alpha[3]*fl[13]+1.732050807568877*alpha[3]*fl[10]+3.0*alpha[0]*fl[8]+1.732050807568877*alpha[0]*fl[4])*dfac_x; 
  incr[9] = 0.125*(1.732050807568877*alpha[3]*fl[15]+alpha[3]*fl[14]+1.732050807568877*alpha[0]*fl[12]+alpha[0]*fl[9])*dfac_x; 
  incr[10] = 0.125*(1.732050807568877*alpha[0]*fl[13]+alpha[0]*fl[10]+1.732050807568877*alpha[3]*fl[8]+alpha[3]*fl[4])*dfac_x; 
  incr[11] = -0.125*(3.0*alpha[0]*fl[11]+1.732050807568877*alpha[0]*fl[7]+3.0*alpha[3]*fl[5]+1.732050807568877*fl[2]*alpha[3])*dfac_x; 
  incr[12] = -0.125*(3.0*alpha[3]*fl[15]+1.732050807568877*alpha[3]*fl[14]+3.0*alpha[0]*fl[12]+1.732050807568877*alpha[0]*fl[9])*dfac_x; 
  incr[13] = -0.125*(3.0*alpha[0]*fl[13]+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[3]*fl[8]+1.732050807568877*alpha[3]*fl[4])*dfac_x; 
  incr[14] = 0.125*(1.732050807568877*alpha[0]*fl[15]+alpha[0]*fl[14]+1.732050807568877*alpha[3]*fl[12]+alpha[3]*fl[9])*dfac_x; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+1.732050807568877*alpha[0]*fl[14]+3.0*alpha[3]*fl[12]+1.732050807568877*alpha[3]*fl[9])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[3]*fr[6]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[1]-1.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.125*(3.0*alpha[3]*fr[6]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.125*(1.732050807568877*alpha[3]*fr[11]-1.0*alpha[3]*fr[7]+1.732050807568877*alpha[0]*fr[5]-1.0*alpha[0]*fr[2])*dfac_x; 
  incr[3] = -0.125*(1.732050807568877*alpha[0]*fr[6]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[1]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_x; 
  incr[4] = -0.125*(1.732050807568877*alpha[3]*fr[13]-1.0*alpha[3]*fr[10]+1.732050807568877*alpha[0]*fr[8]-1.0*alpha[0]*fr[4])*dfac_x; 
  incr[5] = 0.125*(3.0*alpha[3]*fr[11]-1.732050807568877*alpha[3]*fr[7]+3.0*alpha[0]*fr[5]-1.732050807568877*alpha[0]*fr[2])*dfac_x; 
  incr[6] = 0.125*(3.0*alpha[0]*fr[6]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[1]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_x; 
  incr[7] = -0.125*(1.732050807568877*alpha[0]*fr[11]-1.0*alpha[0]*fr[7]+1.732050807568877*alpha[3]*fr[5]-1.0*fr[2]*alpha[3])*dfac_x; 
  incr[8] = 0.125*(3.0*alpha[3]*fr[13]-1.732050807568877*alpha[3]*fr[10]+3.0*alpha[0]*fr[8]-1.732050807568877*alpha[0]*fr[4])*dfac_x; 
  incr[9] = -0.125*(1.732050807568877*alpha[3]*fr[15]-1.0*alpha[3]*fr[14]+1.732050807568877*alpha[0]*fr[12]-1.0*alpha[0]*fr[9])*dfac_x; 
  incr[10] = -0.125*(1.732050807568877*alpha[0]*fr[13]-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[3]*fr[8]-1.0*alpha[3]*fr[4])*dfac_x; 
  incr[11] = 0.125*(3.0*alpha[0]*fr[11]-1.732050807568877*alpha[0]*fr[7]+3.0*alpha[3]*fr[5]-1.732050807568877*fr[2]*alpha[3])*dfac_x; 
  incr[12] = 0.125*(3.0*alpha[3]*fr[15]-1.732050807568877*alpha[3]*fr[14]+3.0*alpha[0]*fr[12]-1.732050807568877*alpha[0]*fr[9])*dfac_x; 
  incr[13] = 0.125*(3.0*alpha[0]*fr[13]-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[3]*fr[8]-1.732050807568877*alpha[3]*fr[4])*dfac_x; 
  incr[14] = -0.125*(1.732050807568877*alpha[0]*fr[15]-1.0*alpha[0]*fr[14]+1.732050807568877*alpha[3]*fr[12]-1.0*alpha[3]*fr[9])*dfac_x; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]-1.732050807568877*alpha[0]*fr[14]+3.0*alpha[3]*fr[12]-1.732050807568877*alpha[3]*fr[9])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf2x2vSer_Y_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.5*(2.0*BdriftY[0]*m_*wv2+(1.732050807568877*BmagInv[0]*Phi[1]-3.0*BmagInv[0]*Phi[3])*dfac_x*q_))/q_; 

  double alpha[16]; 
  alpha[0] = (2.0*BdriftY[0]*m_*wv2)/q_-3.0*BmagInv[0]*Phi[3]*dfac_x+1.732050807568877*BmagInv[0]*Phi[1]*dfac_x; 
  alpha[3] = (1.154700538379252*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[3]*fl[7]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[0]*fl[0])*dfac_y; 
  incr[1] = 0.125*(1.732050807568877*alpha[3]*fl[11]+alpha[3]*fl[6]+1.732050807568877*alpha[0]*fl[5]+alpha[0]*fl[1])*dfac_y; 
  incr[2] = -0.125*(3.0*alpha[3]*fl[7]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[0]*fl[0])*dfac_y; 
  incr[3] = 0.125*(1.732050807568877*alpha[0]*fl[7]+alpha[0]*fl[3]+1.732050807568877*fl[2]*alpha[3]+fl[0]*alpha[3])*dfac_y; 
  incr[4] = 0.125*(1.732050807568877*alpha[3]*fl[14]+alpha[3]*fl[10]+1.732050807568877*alpha[0]*fl[9]+alpha[0]*fl[4])*dfac_y; 
  incr[5] = -0.125*(3.0*alpha[3]*fl[11]+1.732050807568877*alpha[3]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[1])*dfac_y; 
  incr[6] = 0.125*(1.732050807568877*alpha[0]*fl[11]+alpha[0]*fl[6]+1.732050807568877*alpha[3]*fl[5]+fl[1]*alpha[3])*dfac_y; 
  incr[7] = -0.125*(3.0*alpha[0]*fl[7]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_y; 
  incr[8] = 0.125*(1.732050807568877*alpha[3]*fl[15]+alpha[3]*fl[13]+1.732050807568877*alpha[0]*fl[12]+alpha[0]*fl[8])*dfac_y; 
  incr[9] = -0.125*(3.0*alpha[3]*fl[14]+1.732050807568877*alpha[3]*fl[10]+3.0*alpha[0]*fl[9]+1.732050807568877*alpha[0]*fl[4])*dfac_y; 
  incr[10] = 0.125*(1.732050807568877*alpha[0]*fl[14]+alpha[0]*fl[10]+1.732050807568877*alpha[3]*fl[9]+alpha[3]*fl[4])*dfac_y; 
  incr[11] = -0.125*(3.0*alpha[0]*fl[11]+1.732050807568877*alpha[0]*fl[6]+3.0*alpha[3]*fl[5]+1.732050807568877*fl[1]*alpha[3])*dfac_y; 
  incr[12] = -0.125*(3.0*alpha[3]*fl[15]+1.732050807568877*alpha[3]*fl[13]+3.0*alpha[0]*fl[12]+1.732050807568877*alpha[0]*fl[8])*dfac_y; 
  incr[13] = 0.125*(1.732050807568877*alpha[0]*fl[15]+alpha[0]*fl[13]+1.732050807568877*alpha[3]*fl[12]+alpha[3]*fl[8])*dfac_y; 
  incr[14] = -0.125*(3.0*alpha[0]*fl[14]+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[3]*fl[9]+1.732050807568877*alpha[3]*fl[4])*dfac_y; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+1.732050807568877*alpha[0]*fl[13]+3.0*alpha[3]*fl[12]+1.732050807568877*alpha[3]*fl[8])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[3]*fr[7]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[0]*fr[0])*dfac_y; 
  incr[1] = -0.125*(1.732050807568877*alpha[3]*fr[11]-1.0*alpha[3]*fr[6]+1.732050807568877*alpha[0]*fr[5]-1.0*alpha[0]*fr[1])*dfac_y; 
  incr[2] = 0.125*(3.0*alpha[3]*fr[7]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[0]*fr[0])*dfac_y; 
  incr[3] = -0.125*(1.732050807568877*alpha[0]*fr[7]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_y; 
  incr[4] = -0.125*(1.732050807568877*alpha[3]*fr[14]-1.0*alpha[3]*fr[10]+1.732050807568877*alpha[0]*fr[9]-1.0*alpha[0]*fr[4])*dfac_y; 
  incr[5] = 0.125*(3.0*alpha[3]*fr[11]-1.732050807568877*alpha[3]*fr[6]+3.0*alpha[0]*fr[5]-1.732050807568877*alpha[0]*fr[1])*dfac_y; 
  incr[6] = -0.125*(1.732050807568877*alpha[0]*fr[11]-1.0*alpha[0]*fr[6]+1.732050807568877*alpha[3]*fr[5]-1.0*fr[1]*alpha[3])*dfac_y; 
  incr[7] = 0.125*(3.0*alpha[0]*fr[7]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_y; 
  incr[8] = -0.125*(1.732050807568877*alpha[3]*fr[15]-1.0*alpha[3]*fr[13]+1.732050807568877*alpha[0]*fr[12]-1.0*alpha[0]*fr[8])*dfac_y; 
  incr[9] = 0.125*(3.0*alpha[3]*fr[14]-1.732050807568877*alpha[3]*fr[10]+3.0*alpha[0]*fr[9]-1.732050807568877*alpha[0]*fr[4])*dfac_y; 
  incr[10] = -0.125*(1.732050807568877*alpha[0]*fr[14]-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[3]*fr[9]-1.0*alpha[3]*fr[4])*dfac_y; 
  incr[11] = 0.125*(3.0*alpha[0]*fr[11]-1.732050807568877*alpha[0]*fr[6]+3.0*alpha[3]*fr[5]-1.732050807568877*fr[1]*alpha[3])*dfac_y; 
  incr[12] = 0.125*(3.0*alpha[3]*fr[15]-1.732050807568877*alpha[3]*fr[13]+3.0*alpha[0]*fr[12]-1.732050807568877*alpha[0]*fr[8])*dfac_y; 
  incr[13] = -0.125*(1.732050807568877*alpha[0]*fr[15]-1.0*alpha[0]*fr[13]+1.732050807568877*alpha[3]*fr[12]-1.0*alpha[3]*fr[8])*dfac_y; 
  incr[14] = 0.125*(3.0*alpha[0]*fr[14]-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[3]*fr[9]-1.732050807568877*alpha[3]*fr[4])*dfac_y; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]-1.732050807568877*alpha[0]*fr[13]+3.0*alpha[3]*fr[12]-1.732050807568877*alpha[3]*fr[8])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf2x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.5*((1.732050807568877*BdriftY[0]*Phi[2]*dfac_v*dfac_y+1.732050807568877*BdriftX[0]*Phi[1]*dfac_v*dfac_x)*wv-1.732050807568877*BdriftY[0]*Phi[2]*dfac_y-1.732050807568877*BdriftX[0]*Phi[1]*dfac_x))/dfac_v; 

  double alpha[16]; 
  alpha[0] = (-1.732050807568877*BdriftY[0]*Phi[2]*dfac_y*wv)-1.732050807568877*BdriftX[0]*Phi[1]*dfac_x*wv+(1.732050807568877*BdriftY[0]*Phi[2]*dfac_y)/dfac_v+(1.732050807568877*BdriftX[0]*Phi[1]*dfac_x)/dfac_v; 
  alpha[1] = (1.732050807568877*BdriftY[0]*Phi[3]*dfac_y)/dfac_v-1.732050807568877*BdriftY[0]*Phi[3]*dfac_y*wv; 
  alpha[2] = (1.732050807568877*BdriftX[0]*Phi[3]*dfac_x)/dfac_v-1.732050807568877*BdriftX[0]*Phi[3]*dfac_x*wv; 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[2]*fl[7]+1.732050807568877*alpha[1]*fl[6]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.125*(1.732050807568877*alpha[2]*fl[11]+1.732050807568877*alpha[0]*fl[6]+alpha[2]*fl[5]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.125*(1.732050807568877*alpha[1]*fl[11]+1.732050807568877*alpha[0]*fl[7]+alpha[1]*fl[5]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.125*(3.0*alpha[2]*fl[7]+3.0*alpha[1]*fl[6]+3.0*alpha[0]*fl[3]+1.732050807568877*alpha[2]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[4] = 0.125*(1.732050807568877*alpha[2]*fl[14]+1.732050807568877*alpha[1]*fl[13]+1.732050807568877*alpha[0]*fl[10]+alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4])*dfac_v; 
  incr[5] = 0.125*(1.732050807568877*alpha[0]*fl[11]+1.732050807568877*alpha[1]*fl[7]+1.732050807568877*alpha[2]*fl[6]+alpha[0]*fl[5]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.125*(3.0*alpha[2]*fl[11]+3.0*alpha[0]*fl[6]+1.732050807568877*alpha[2]*fl[5]+3.0*alpha[1]*fl[3]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 
  incr[7] = -0.125*(3.0*alpha[1]*fl[11]+3.0*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[5]+3.0*alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+1.732050807568877*fl[0]*alpha[2])*dfac_v; 
  incr[8] = 0.125*(1.732050807568877*alpha[2]*fl[15]+1.732050807568877*alpha[0]*fl[13]+alpha[2]*fl[12]+1.732050807568877*alpha[1]*fl[10]+alpha[0]*fl[8]+alpha[1]*fl[4])*dfac_v; 
  incr[9] = 0.125*(1.732050807568877*alpha[1]*fl[15]+1.732050807568877*alpha[0]*fl[14]+alpha[1]*fl[12]+1.732050807568877*alpha[2]*fl[10]+alpha[0]*fl[9]+alpha[2]*fl[4])*dfac_v; 
  incr[10] = -0.125*(3.0*alpha[2]*fl[14]+3.0*alpha[1]*fl[13]+3.0*alpha[0]*fl[10]+1.732050807568877*alpha[2]*fl[9]+1.732050807568877*alpha[1]*fl[8]+1.732050807568877*alpha[0]*fl[4])*dfac_v; 
  incr[11] = -0.125*(3.0*alpha[0]*fl[11]+3.0*alpha[1]*fl[7]+3.0*alpha[2]*fl[6]+1.732050807568877*alpha[0]*fl[5]+1.732050807568877*alpha[1]*fl[2]+1.732050807568877*fl[1]*alpha[2])*dfac_v; 
  incr[12] = 0.125*(1.732050807568877*alpha[0]*fl[15]+1.732050807568877*alpha[1]*fl[14]+1.732050807568877*alpha[2]*fl[13]+alpha[0]*fl[12]+alpha[1]*fl[9]+alpha[2]*fl[8])*dfac_v; 
  incr[13] = -0.125*(3.0*alpha[2]*fl[15]+3.0*alpha[0]*fl[13]+1.732050807568877*alpha[2]*fl[12]+3.0*alpha[1]*fl[10]+1.732050807568877*alpha[0]*fl[8]+1.732050807568877*alpha[1]*fl[4])*dfac_v; 
  incr[14] = -0.125*(3.0*alpha[1]*fl[15]+3.0*alpha[0]*fl[14]+1.732050807568877*alpha[1]*fl[12]+3.0*alpha[2]*fl[10]+1.732050807568877*alpha[0]*fl[9]+1.732050807568877*alpha[2]*fl[4])*dfac_v; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+3.0*alpha[1]*fl[14]+3.0*alpha[2]*fl[13]+1.732050807568877*alpha[0]*fl[12]+1.732050807568877*alpha[1]*fl[9]+1.732050807568877*alpha[2]*fl[8])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[2]*fr[7]+1.732050807568877*alpha[1]*fr[6]+1.732050807568877*alpha[0]*fr[3]-1.0*alpha[2]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.125*(1.732050807568877*alpha[2]*fr[11]+1.732050807568877*alpha[0]*fr[6]-1.0*alpha[2]*fr[5]+1.732050807568877*alpha[1]*fr[3]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = -0.125*(1.732050807568877*alpha[1]*fr[11]+1.732050807568877*alpha[0]*fr[7]-1.0*alpha[1]*fr[5]+1.732050807568877*alpha[2]*fr[3]-1.0*alpha[0]*fr[2]-1.0*fr[0]*alpha[2])*dfac_v; 
  incr[3] = 0.125*(3.0*alpha[2]*fr[7]+3.0*alpha[1]*fr[6]+3.0*alpha[0]*fr[3]-1.732050807568877*alpha[2]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[4] = -0.125*(1.732050807568877*alpha[2]*fr[14]+1.732050807568877*alpha[1]*fr[13]+1.732050807568877*alpha[0]*fr[10]-1.0*alpha[2]*fr[9]-1.0*alpha[1]*fr[8]-1.0*alpha[0]*fr[4])*dfac_v; 
  incr[5] = -0.125*(1.732050807568877*alpha[0]*fr[11]+1.732050807568877*alpha[1]*fr[7]+1.732050807568877*alpha[2]*fr[6]-1.0*alpha[0]*fr[5]-1.0*alpha[1]*fr[2]-1.0*fr[1]*alpha[2])*dfac_v; 
  incr[6] = 0.125*(3.0*alpha[2]*fr[11]+3.0*alpha[0]*fr[6]-1.732050807568877*alpha[2]*fr[5]+3.0*alpha[1]*fr[3]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 
  incr[7] = 0.125*(3.0*alpha[1]*fr[11]+3.0*alpha[0]*fr[7]-1.732050807568877*alpha[1]*fr[5]+3.0*alpha[2]*fr[3]-1.732050807568877*alpha[0]*fr[2]-1.732050807568877*fr[0]*alpha[2])*dfac_v; 
  incr[8] = -0.125*(1.732050807568877*alpha[2]*fr[15]+1.732050807568877*alpha[0]*fr[13]-1.0*alpha[2]*fr[12]+1.732050807568877*alpha[1]*fr[10]-1.0*alpha[0]*fr[8]-1.0*alpha[1]*fr[4])*dfac_v; 
  incr[9] = -0.125*(1.732050807568877*alpha[1]*fr[15]+1.732050807568877*alpha[0]*fr[14]-1.0*alpha[1]*fr[12]+1.732050807568877*alpha[2]*fr[10]-1.0*alpha[0]*fr[9]-1.0*alpha[2]*fr[4])*dfac_v; 
  incr[10] = 0.125*(3.0*alpha[2]*fr[14]+3.0*alpha[1]*fr[13]+3.0*alpha[0]*fr[10]-1.732050807568877*alpha[2]*fr[9]-1.732050807568877*alpha[1]*fr[8]-1.732050807568877*alpha[0]*fr[4])*dfac_v; 
  incr[11] = 0.125*(3.0*alpha[0]*fr[11]+3.0*alpha[1]*fr[7]+3.0*alpha[2]*fr[6]-1.732050807568877*alpha[0]*fr[5]-1.732050807568877*alpha[1]*fr[2]-1.732050807568877*fr[1]*alpha[2])*dfac_v; 
  incr[12] = -0.125*(1.732050807568877*alpha[0]*fr[15]+1.732050807568877*alpha[1]*fr[14]+1.732050807568877*alpha[2]*fr[13]-1.0*alpha[0]*fr[12]-1.0*alpha[1]*fr[9]-1.0*alpha[2]*fr[8])*dfac_v; 
  incr[13] = 0.125*(3.0*alpha[2]*fr[15]+3.0*alpha[0]*fr[13]-1.732050807568877*alpha[2]*fr[12]+3.0*alpha[1]*fr[10]-1.732050807568877*alpha[0]*fr[8]-1.732050807568877*alpha[1]*fr[4])*dfac_v; 
  incr[14] = 0.125*(3.0*alpha[1]*fr[15]+3.0*alpha[0]*fr[14]-1.732050807568877*alpha[1]*fr[12]+3.0*alpha[2]*fr[10]-1.732050807568877*alpha[0]*fr[9]-1.732050807568877*alpha[2]*fr[4])*dfac_v; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]+3.0*alpha[1]*fr[14]+3.0*alpha[2]*fr[13]-1.732050807568877*alpha[0]*fr[12]-1.732050807568877*alpha[1]*fr[9]-1.732050807568877*alpha[2]*fr[8])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.5*(2.0*BdriftX[0]*m_*wv2+((3.464101615137754*Apar[2]-6.0*Apar[3])*dfac_y-1.732050807568877*BdriftX[0]*Apar[1]+Apar[0]*BdriftX[0])*q_*wv+(3.0*BmagInv[0]*Phi[3]-1.732050807568877*BmagInv[0]*Phi[2])*dfac_y*q_))/q_; 

  double alpha[16]; 
  alpha[0] = (2.0*BdriftX[0]*m_*wv2)/q_-6.0*Apar[3]*dfac_y*wv+3.464101615137754*Apar[2]*dfac_y*wv-1.732050807568877*BdriftX[0]*Apar[1]*wv+Apar[0]*BdriftX[0]*wv+3.0*BmagInv[0]*Phi[3]*dfac_y-1.732050807568877*BmagInv[0]*Phi[2]*dfac_y; 
  alpha[2] = 1.0*BdriftX[0]*Apar[2]*wv-1.732050807568878*BdriftX[0]*Apar[3]*wv; 
  alpha[3] = (1.154700538379252*BdriftX[0]*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[3]*fl[6]+1.732050807568877*alpha[2]*fl[5]+alpha[3]*fl[3]+alpha[2]*fl[2]+1.732050807568877*alpha[0]*fl[1]+alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.125*(3.0*alpha[3]*fl[6]+3.0*alpha[2]*fl[5]+1.732050807568877*alpha[3]*fl[3]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.125*(1.732050807568877*alpha[3]*fl[11]+alpha[3]*fl[7]+1.732050807568877*alpha[0]*fl[5]+alpha[0]*fl[2]+1.732050807568877*fl[1]*alpha[2]+fl[0]*alpha[2])*dfac_x; 
  incr[3] = 0.125*(1.732050807568877*alpha[2]*fl[11]+alpha[2]*fl[7]+1.732050807568877*alpha[0]*fl[6]+alpha[0]*fl[3]+1.732050807568877*fl[1]*alpha[3]+fl[0]*alpha[3])*dfac_x; 
  incr[4] = 0.125*(1.732050807568877*alpha[3]*fl[13]+1.732050807568877*alpha[2]*fl[12]+alpha[3]*fl[10]+alpha[2]*fl[9]+1.732050807568877*alpha[0]*fl[8]+alpha[0]*fl[4])*dfac_x; 
  incr[5] = -0.125*(3.0*alpha[3]*fl[11]+1.732050807568877*alpha[3]*fl[7]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[2]+3.0*fl[1]*alpha[2]+1.732050807568877*fl[0]*alpha[2])*dfac_x; 
  incr[6] = -0.125*(3.0*alpha[2]*fl[11]+1.732050807568877*alpha[2]*fl[7]+3.0*alpha[0]*fl[6]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[1]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_x; 
  incr[7] = 0.125*(1.732050807568877*alpha[0]*fl[11]+alpha[0]*fl[7]+1.732050807568877*alpha[2]*fl[6]+1.732050807568877*alpha[3]*fl[5]+alpha[2]*fl[3]+fl[2]*alpha[3])*dfac_x; 
  incr[8] = -0.125*(3.0*alpha[3]*fl[13]+3.0*alpha[2]*fl[12]+1.732050807568877*alpha[3]*fl[10]+1.732050807568877*alpha[2]*fl[9]+3.0*alpha[0]*fl[8]+1.732050807568877*alpha[0]*fl[4])*dfac_x; 
  incr[9] = 0.125*(1.732050807568877*alpha[3]*fl[15]+alpha[3]*fl[14]+1.732050807568877*alpha[0]*fl[12]+alpha[0]*fl[9]+1.732050807568877*alpha[2]*fl[8]+alpha[2]*fl[4])*dfac_x; 
  incr[10] = 0.125*(1.732050807568877*alpha[2]*fl[15]+alpha[2]*fl[14]+1.732050807568877*alpha[0]*fl[13]+alpha[0]*fl[10]+1.732050807568877*alpha[3]*fl[8]+alpha[3]*fl[4])*dfac_x; 
  incr[11] = -0.125*(3.0*alpha[0]*fl[11]+1.732050807568877*alpha[0]*fl[7]+3.0*alpha[2]*fl[6]+3.0*alpha[3]*fl[5]+1.732050807568877*alpha[2]*fl[3]+1.732050807568877*fl[2]*alpha[3])*dfac_x; 
  incr[12] = -0.125*(3.0*alpha[3]*fl[15]+1.732050807568877*alpha[3]*fl[14]+3.0*alpha[0]*fl[12]+1.732050807568877*alpha[0]*fl[9]+3.0*alpha[2]*fl[8]+1.732050807568877*alpha[2]*fl[4])*dfac_x; 
  incr[13] = -0.125*(3.0*alpha[2]*fl[15]+1.732050807568877*alpha[2]*fl[14]+3.0*alpha[0]*fl[13]+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[3]*fl[8]+1.732050807568877*alpha[3]*fl[4])*dfac_x; 
  incr[14] = 0.125*(1.732050807568877*alpha[0]*fl[15]+alpha[0]*fl[14]+1.732050807568877*alpha[2]*fl[13]+1.732050807568877*alpha[3]*fl[12]+alpha[2]*fl[10]+alpha[3]*fl[9])*dfac_x; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+1.732050807568877*alpha[0]*fl[14]+3.0*alpha[2]*fl[13]+3.0*alpha[3]*fl[12]+1.732050807568877*alpha[2]*fl[10]+1.732050807568877*alpha[3]*fl[9])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[3]*fr[6]+1.732050807568877*alpha[2]*fr[5]-1.0*alpha[3]*fr[3]-1.0*alpha[2]*fr[2]+1.732050807568877*alpha[0]*fr[1]-1.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.125*(3.0*alpha[3]*fr[6]+3.0*alpha[2]*fr[5]-1.732050807568877*alpha[3]*fr[3]-1.732050807568877*alpha[2]*fr[2]+3.0*alpha[0]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.125*(1.732050807568877*alpha[3]*fr[11]-1.0*alpha[3]*fr[7]+1.732050807568877*alpha[0]*fr[5]-1.0*alpha[0]*fr[2]+1.732050807568877*fr[1]*alpha[2]-1.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = -0.125*(1.732050807568877*alpha[2]*fr[11]-1.0*alpha[2]*fr[7]+1.732050807568877*alpha[0]*fr[6]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[1]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_x; 
  incr[4] = -0.125*(1.732050807568877*alpha[3]*fr[13]+1.732050807568877*alpha[2]*fr[12]-1.0*alpha[3]*fr[10]-1.0*alpha[2]*fr[9]+1.732050807568877*alpha[0]*fr[8]-1.0*alpha[0]*fr[4])*dfac_x; 
  incr[5] = 0.125*(3.0*alpha[3]*fr[11]-1.732050807568877*alpha[3]*fr[7]+3.0*alpha[0]*fr[5]-1.732050807568877*alpha[0]*fr[2]+3.0*fr[1]*alpha[2]-1.732050807568877*fr[0]*alpha[2])*dfac_x; 
  incr[6] = 0.125*(3.0*alpha[2]*fr[11]-1.732050807568877*alpha[2]*fr[7]+3.0*alpha[0]*fr[6]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[1]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_x; 
  incr[7] = -0.125*(1.732050807568877*alpha[0]*fr[11]-1.0*alpha[0]*fr[7]+1.732050807568877*alpha[2]*fr[6]+1.732050807568877*alpha[3]*fr[5]-1.0*alpha[2]*fr[3]-1.0*fr[2]*alpha[3])*dfac_x; 
  incr[8] = 0.125*(3.0*alpha[3]*fr[13]+3.0*alpha[2]*fr[12]-1.732050807568877*alpha[3]*fr[10]-1.732050807568877*alpha[2]*fr[9]+3.0*alpha[0]*fr[8]-1.732050807568877*alpha[0]*fr[4])*dfac_x; 
  incr[9] = -0.125*(1.732050807568877*alpha[3]*fr[15]-1.0*alpha[3]*fr[14]+1.732050807568877*alpha[0]*fr[12]-1.0*alpha[0]*fr[9]+1.732050807568877*alpha[2]*fr[8]-1.0*alpha[2]*fr[4])*dfac_x; 
  incr[10] = -0.125*(1.732050807568877*alpha[2]*fr[15]-1.0*alpha[2]*fr[14]+1.732050807568877*alpha[0]*fr[13]-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[3]*fr[8]-1.0*alpha[3]*fr[4])*dfac_x; 
  incr[11] = 0.125*(3.0*alpha[0]*fr[11]-1.732050807568877*alpha[0]*fr[7]+3.0*alpha[2]*fr[6]+3.0*alpha[3]*fr[5]-1.732050807568877*alpha[2]*fr[3]-1.732050807568877*fr[2]*alpha[3])*dfac_x; 
  incr[12] = 0.125*(3.0*alpha[3]*fr[15]-1.732050807568877*alpha[3]*fr[14]+3.0*alpha[0]*fr[12]-1.732050807568877*alpha[0]*fr[9]+3.0*alpha[2]*fr[8]-1.732050807568877*alpha[2]*fr[4])*dfac_x; 
  incr[13] = 0.125*(3.0*alpha[2]*fr[15]-1.732050807568877*alpha[2]*fr[14]+3.0*alpha[0]*fr[13]-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[3]*fr[8]-1.732050807568877*alpha[3]*fr[4])*dfac_x; 
  incr[14] = -0.125*(1.732050807568877*alpha[0]*fr[15]-1.0*alpha[0]*fr[14]+1.732050807568877*alpha[2]*fr[13]+1.732050807568877*alpha[3]*fr[12]-1.0*alpha[2]*fr[10]-1.0*alpha[3]*fr[9])*dfac_x; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]-1.732050807568877*alpha[0]*fr[14]+3.0*alpha[2]*fr[13]+3.0*alpha[3]*fr[12]-1.732050807568877*alpha[2]*fr[10]-1.732050807568877*alpha[3]*fr[9])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_Y_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.5*(2.0*BdriftY[0]*m_*wv2+((6.0*Apar[3]-3.464101615137754*Apar[1])*dfac_x-1.732050807568877*BdriftY[0]*Apar[2]+Apar[0]*BdriftY[0])*q_*wv+(1.732050807568877*BmagInv[0]*Phi[1]-3.0*BmagInv[0]*Phi[3])*dfac_x*q_))/q_; 

  double alpha[16]; 
  alpha[0] = (2.0*BdriftY[0]*m_*wv2)/q_+6.0*Apar[3]*dfac_x*wv-3.464101615137754*Apar[1]*dfac_x*wv-1.732050807568877*BdriftY[0]*Apar[2]*wv+Apar[0]*BdriftY[0]*wv-3.0*BmagInv[0]*Phi[3]*dfac_x+1.732050807568877*BmagInv[0]*Phi[1]*dfac_x; 
  alpha[1] = 1.0*BdriftY[0]*Apar[1]*wv-1.732050807568878*BdriftY[0]*Apar[3]*wv; 
  alpha[3] = (1.154700538379252*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[3]*fl[7]+1.732050807568877*alpha[1]*fl[5]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_y; 
  incr[1] = 0.125*(1.732050807568877*alpha[3]*fl[11]+alpha[3]*fl[6]+1.732050807568877*alpha[0]*fl[5]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_y; 
  incr[2] = -0.125*(3.0*alpha[3]*fl[7]+3.0*alpha[1]*fl[5]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_y; 
  incr[3] = 0.125*(1.732050807568877*alpha[1]*fl[11]+1.732050807568877*alpha[0]*fl[7]+alpha[1]*fl[6]+alpha[0]*fl[3]+1.732050807568877*fl[2]*alpha[3]+fl[0]*alpha[3])*dfac_y; 
  incr[4] = 0.125*(1.732050807568877*alpha[3]*fl[14]+1.732050807568877*alpha[1]*fl[12]+alpha[3]*fl[10]+1.732050807568877*alpha[0]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4])*dfac_y; 
  incr[5] = -0.125*(3.0*alpha[3]*fl[11]+1.732050807568877*alpha[3]*fl[6]+3.0*alpha[0]*fl[5]+3.0*alpha[1]*fl[2]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_y; 
  incr[6] = 0.125*(1.732050807568877*alpha[0]*fl[11]+1.732050807568877*alpha[1]*fl[7]+alpha[0]*fl[6]+1.732050807568877*alpha[3]*fl[5]+alpha[1]*fl[3]+fl[1]*alpha[3])*dfac_y; 
  incr[7] = -0.125*(3.0*alpha[1]*fl[11]+3.0*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[6]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_y; 
  incr[8] = 0.125*(1.732050807568877*alpha[3]*fl[15]+alpha[3]*fl[13]+1.732050807568877*alpha[0]*fl[12]+1.732050807568877*alpha[1]*fl[9]+alpha[0]*fl[8]+alpha[1]*fl[4])*dfac_y; 
  incr[9] = -0.125*(3.0*alpha[3]*fl[14]+3.0*alpha[1]*fl[12]+1.732050807568877*alpha[3]*fl[10]+3.0*alpha[0]*fl[9]+1.732050807568877*alpha[1]*fl[8]+1.732050807568877*alpha[0]*fl[4])*dfac_y; 
  incr[10] = 0.125*(1.732050807568877*alpha[1]*fl[15]+1.732050807568877*alpha[0]*fl[14]+alpha[1]*fl[13]+alpha[0]*fl[10]+1.732050807568877*alpha[3]*fl[9]+alpha[3]*fl[4])*dfac_y; 
  incr[11] = -0.125*(3.0*alpha[0]*fl[11]+3.0*alpha[1]*fl[7]+1.732050807568877*alpha[0]*fl[6]+3.0*alpha[3]*fl[5]+1.732050807568877*alpha[1]*fl[3]+1.732050807568877*fl[1]*alpha[3])*dfac_y; 
  incr[12] = -0.125*(3.0*alpha[3]*fl[15]+1.732050807568877*alpha[3]*fl[13]+3.0*alpha[0]*fl[12]+3.0*alpha[1]*fl[9]+1.732050807568877*alpha[0]*fl[8]+1.732050807568877*alpha[1]*fl[4])*dfac_y; 
  incr[13] = 0.125*(1.732050807568877*alpha[0]*fl[15]+1.732050807568877*alpha[1]*fl[14]+alpha[0]*fl[13]+1.732050807568877*alpha[3]*fl[12]+alpha[1]*fl[10]+alpha[3]*fl[8])*dfac_y; 
  incr[14] = -0.125*(3.0*alpha[1]*fl[15]+3.0*alpha[0]*fl[14]+1.732050807568877*alpha[1]*fl[13]+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[3]*fl[9]+1.732050807568877*alpha[3]*fl[4])*dfac_y; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+3.0*alpha[1]*fl[14]+1.732050807568877*alpha[0]*fl[13]+3.0*alpha[3]*fl[12]+1.732050807568877*alpha[1]*fl[10]+1.732050807568877*alpha[3]*fl[8])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[3]*fr[7]+1.732050807568877*alpha[1]*fr[5]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_y; 
  incr[1] = -0.125*(1.732050807568877*alpha[3]*fr[11]-1.0*alpha[3]*fr[6]+1.732050807568877*alpha[0]*fr[5]+1.732050807568877*alpha[1]*fr[2]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_y; 
  incr[2] = 0.125*(3.0*alpha[3]*fr[7]+3.0*alpha[1]*fr[5]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_y; 
  incr[3] = -0.125*(1.732050807568877*alpha[1]*fr[11]+1.732050807568877*alpha[0]*fr[7]-1.0*alpha[1]*fr[6]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_y; 
  incr[4] = -0.125*(1.732050807568877*alpha[3]*fr[14]+1.732050807568877*alpha[1]*fr[12]-1.0*alpha[3]*fr[10]+1.732050807568877*alpha[0]*fr[9]-1.0*alpha[1]*fr[8]-1.0*alpha[0]*fr[4])*dfac_y; 
  incr[5] = 0.125*(3.0*alpha[3]*fr[11]-1.732050807568877*alpha[3]*fr[6]+3.0*alpha[0]*fr[5]+3.0*alpha[1]*fr[2]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_y; 
  incr[6] = -0.125*(1.732050807568877*alpha[0]*fr[11]+1.732050807568877*alpha[1]*fr[7]-1.0*alpha[0]*fr[6]+1.732050807568877*alpha[3]*fr[5]-1.0*alpha[1]*fr[3]-1.0*fr[1]*alpha[3])*dfac_y; 
  incr[7] = 0.125*(3.0*alpha[1]*fr[11]+3.0*alpha[0]*fr[7]-1.732050807568877*alpha[1]*fr[6]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_y; 
  incr[8] = -0.125*(1.732050807568877*alpha[3]*fr[15]-1.0*alpha[3]*fr[13]+1.732050807568877*alpha[0]*fr[12]+1.732050807568877*alpha[1]*fr[9]-1.0*alpha[0]*fr[8]-1.0*alpha[1]*fr[4])*dfac_y; 
  incr[9] = 0.125*(3.0*alpha[3]*fr[14]+3.0*alpha[1]*fr[12]-1.732050807568877*alpha[3]*fr[10]+3.0*alpha[0]*fr[9]-1.732050807568877*alpha[1]*fr[8]-1.732050807568877*alpha[0]*fr[4])*dfac_y; 
  incr[10] = -0.125*(1.732050807568877*alpha[1]*fr[15]+1.732050807568877*alpha[0]*fr[14]-1.0*alpha[1]*fr[13]-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[3]*fr[9]-1.0*alpha[3]*fr[4])*dfac_y; 
  incr[11] = 0.125*(3.0*alpha[0]*fr[11]+3.0*alpha[1]*fr[7]-1.732050807568877*alpha[0]*fr[6]+3.0*alpha[3]*fr[5]-1.732050807568877*alpha[1]*fr[3]-1.732050807568877*fr[1]*alpha[3])*dfac_y; 
  incr[12] = 0.125*(3.0*alpha[3]*fr[15]-1.732050807568877*alpha[3]*fr[13]+3.0*alpha[0]*fr[12]+3.0*alpha[1]*fr[9]-1.732050807568877*alpha[0]*fr[8]-1.732050807568877*alpha[1]*fr[4])*dfac_y; 
  incr[13] = -0.125*(1.732050807568877*alpha[0]*fr[15]+1.732050807568877*alpha[1]*fr[14]-1.0*alpha[0]*fr[13]+1.732050807568877*alpha[3]*fr[12]-1.0*alpha[1]*fr[10]-1.0*alpha[3]*fr[8])*dfac_y; 
  incr[14] = 0.125*(3.0*alpha[1]*fr[15]+3.0*alpha[0]*fr[14]-1.732050807568877*alpha[1]*fr[13]-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[3]*fr[9]-1.732050807568877*alpha[3]*fr[4])*dfac_y; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]+3.0*alpha[1]*fr[14]-1.732050807568877*alpha[0]*fr[13]+3.0*alpha[3]*fr[12]-1.732050807568877*alpha[1]*fr[10]-1.732050807568877*alpha[3]*fr[8])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*((3.464101615137754*BdriftY[0]*Phi[2]*dfac_v*dfac_y+3.464101615137754*BdriftX[0]*Phi[1]*dfac_v*dfac_x)*m_*wv+(((6.0*Phi[1]*Apar[2]-6.0*Apar[1]*Phi[2])*dfac_v*dfac_x+(1.732050807568877*BdriftY[0]*Apar[1]*Phi[3]+1.732050807568877*Apar[0]*BdriftY[0]*Phi[2])*dfac_v)*dfac_y+(1.732050807568877*BdriftX[0]*Apar[2]*Phi[3]+1.732050807568877*Apar[0]*BdriftX[0]*Phi[1])*dfac_v*dfac_x+4.0*dApardt[0]*dfac_v)*q_+((-3.464101615137754*BdriftY[0]*Phi[2]*dfac_y)-3.464101615137754*BdriftX[0]*Phi[1]*dfac_x)*m_))/(dfac_v*m_); 

  double alpha[16]; 
  alpha[0] = (-1.732050807568877*BdriftY[0]*Phi[2]*dfac_y*wv)-1.732050807568877*BdriftX[0]*Phi[1]*dfac_x*wv+(3.0*Apar[1]*Phi[2]*dfac_x*dfac_y*q_)/m_-(3.0*Phi[1]*Apar[2]*dfac_x*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[1]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[0]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[2]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*Apar[0]*BdriftX[0]*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[0]*q_)/m_+(1.732050807568877*BdriftY[0]*Phi[2]*dfac_y)/dfac_v+(1.732050807568877*BdriftX[0]*Phi[1]*dfac_x)/dfac_v; 
  alpha[1] = (-1.732050807568877*BdriftY[0]*Phi[3]*dfac_y*wv)+(3.0*Apar[1]*Phi[3]*dfac_x*dfac_y*q_)/m_-(3.0*Phi[1]*Apar[3]*dfac_x*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[0]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[1]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[3]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[1]*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[1]*q_)/m_+(1.732050807568877*BdriftY[0]*Phi[3]*dfac_y)/dfac_v; 
  alpha[2] = (-1.732050807568877*BdriftX[0]*Phi[3]*dfac_x*wv)-(3.0*Apar[2]*Phi[3]*dfac_x*dfac_y*q_)/m_+(3.0*Phi[2]*Apar[3]*dfac_x*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[3]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[2]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftX[0]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Phi[1]*Apar[2]*dfac_x*q_)/m_-(2.0*dApardt[2]*q_)/m_+(1.732050807568877*BdriftX[0]*Phi[3]*dfac_x)/dfac_v; 
  alpha[5] = (-(0.8660254037844386*BdriftY[0]*Apar[2]*Phi[3]*dfac_y*q_)/m_)-(0.8660254037844386*BdriftY[0]*Phi[2]*Apar[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[1]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Phi[1]*Apar[3]*dfac_x*q_)/m_-(2.0*dApardt[3]*q_)/m_; 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[5]*fl[11]+1.732050807568877*alpha[2]*fl[7]+1.732050807568877*alpha[1]*fl[6]+alpha[5]*fl[5]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.125*(1.732050807568877*alpha[2]*fl[11]+1.732050807568877*alpha[5]*fl[7]+1.732050807568877*alpha[0]*fl[6]+alpha[2]*fl[5]+fl[2]*alpha[5]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.125*(1.732050807568877*alpha[1]*fl[11]+1.732050807568877*alpha[0]*fl[7]+1.732050807568877*alpha[5]*fl[6]+alpha[1]*fl[5]+fl[1]*alpha[5]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.125*(3.0*alpha[5]*fl[11]+3.0*alpha[2]*fl[7]+3.0*alpha[1]*fl[6]+1.732050807568877*alpha[5]*fl[5]+3.0*alpha[0]*fl[3]+1.732050807568877*alpha[2]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[4] = 0.125*(1.732050807568877*alpha[5]*fl[15]+1.732050807568877*alpha[2]*fl[14]+1.732050807568877*alpha[1]*fl[13]+alpha[5]*fl[12]+1.732050807568877*alpha[0]*fl[10]+alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4])*dfac_v; 
  incr[5] = 0.125*(1.732050807568877*alpha[0]*fl[11]+1.732050807568877*alpha[1]*fl[7]+1.732050807568877*alpha[2]*fl[6]+alpha[0]*fl[5]+1.732050807568877*fl[3]*alpha[5]+fl[0]*alpha[5]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.125*(3.0*alpha[2]*fl[11]+3.0*alpha[5]*fl[7]+3.0*alpha[0]*fl[6]+1.732050807568877*alpha[2]*fl[5]+1.732050807568877*fl[2]*alpha[5]+3.0*alpha[1]*fl[3]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 
  incr[7] = -0.125*(3.0*alpha[1]*fl[11]+3.0*alpha[0]*fl[7]+3.0*alpha[5]*fl[6]+1.732050807568877*alpha[1]*fl[5]+1.732050807568877*fl[1]*alpha[5]+3.0*alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+1.732050807568877*fl[0]*alpha[2])*dfac_v; 
  incr[8] = 0.125*(1.732050807568877*alpha[2]*fl[15]+1.732050807568877*alpha[5]*fl[14]+1.732050807568877*alpha[0]*fl[13]+alpha[2]*fl[12]+1.732050807568877*alpha[1]*fl[10]+alpha[5]*fl[9]+alpha[0]*fl[8]+alpha[1]*fl[4])*dfac_v; 
  incr[9] = 0.125*(1.732050807568877*alpha[1]*fl[15]+1.732050807568877*alpha[0]*fl[14]+1.732050807568877*alpha[5]*fl[13]+alpha[1]*fl[12]+1.732050807568877*alpha[2]*fl[10]+alpha[0]*fl[9]+alpha[5]*fl[8]+alpha[2]*fl[4])*dfac_v; 
  incr[10] = -0.125*(3.0*alpha[5]*fl[15]+3.0*alpha[2]*fl[14]+3.0*alpha[1]*fl[13]+1.732050807568877*alpha[5]*fl[12]+3.0*alpha[0]*fl[10]+1.732050807568877*alpha[2]*fl[9]+1.732050807568877*alpha[1]*fl[8]+1.732050807568877*alpha[0]*fl[4])*dfac_v; 
  incr[11] = -0.125*(3.0*alpha[0]*fl[11]+3.0*alpha[1]*fl[7]+3.0*alpha[2]*fl[6]+1.732050807568877*alpha[0]*fl[5]+3.0*fl[3]*alpha[5]+1.732050807568877*fl[0]*alpha[5]+1.732050807568877*alpha[1]*fl[2]+1.732050807568877*fl[1]*alpha[2])*dfac_v; 
  incr[12] = 0.125*(1.732050807568877*alpha[0]*fl[15]+1.732050807568877*alpha[1]*fl[14]+1.732050807568877*alpha[2]*fl[13]+alpha[0]*fl[12]+1.732050807568877*alpha[5]*fl[10]+alpha[1]*fl[9]+alpha[2]*fl[8]+fl[4]*alpha[5])*dfac_v; 
  incr[13] = -0.125*(3.0*alpha[2]*fl[15]+3.0*alpha[5]*fl[14]+3.0*alpha[0]*fl[13]+1.732050807568877*alpha[2]*fl[12]+3.0*alpha[1]*fl[10]+1.732050807568877*alpha[5]*fl[9]+1.732050807568877*alpha[0]*fl[8]+1.732050807568877*alpha[1]*fl[4])*dfac_v; 
  incr[14] = -0.125*(3.0*alpha[1]*fl[15]+3.0*alpha[0]*fl[14]+3.0*alpha[5]*fl[13]+1.732050807568877*alpha[1]*fl[12]+3.0*alpha[2]*fl[10]+1.732050807568877*alpha[0]*fl[9]+1.732050807568877*alpha[5]*fl[8]+1.732050807568877*alpha[2]*fl[4])*dfac_v; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+3.0*alpha[1]*fl[14]+3.0*alpha[2]*fl[13]+1.732050807568877*alpha[0]*fl[12]+3.0*alpha[5]*fl[10]+1.732050807568877*alpha[1]*fl[9]+1.732050807568877*alpha[2]*fl[8]+1.732050807568877*fl[4]*alpha[5])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[5]*fr[11]+1.732050807568877*alpha[2]*fr[7]+1.732050807568877*alpha[1]*fr[6]-1.0*alpha[5]*fr[5]+1.732050807568877*alpha[0]*fr[3]-1.0*alpha[2]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.125*(1.732050807568877*alpha[2]*fr[11]+1.732050807568877*alpha[5]*fr[7]+1.732050807568877*alpha[0]*fr[6]-1.0*alpha[2]*fr[5]-1.0*fr[2]*alpha[5]+1.732050807568877*alpha[1]*fr[3]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = -0.125*(1.732050807568877*alpha[1]*fr[11]+1.732050807568877*alpha[0]*fr[7]+1.732050807568877*alpha[5]*fr[6]-1.0*alpha[1]*fr[5]-1.0*fr[1]*alpha[5]+1.732050807568877*alpha[2]*fr[3]-1.0*alpha[0]*fr[2]-1.0*fr[0]*alpha[2])*dfac_v; 
  incr[3] = 0.125*(3.0*alpha[5]*fr[11]+3.0*alpha[2]*fr[7]+3.0*alpha[1]*fr[6]-1.732050807568877*alpha[5]*fr[5]+3.0*alpha[0]*fr[3]-1.732050807568877*alpha[2]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[4] = -0.125*(1.732050807568877*alpha[5]*fr[15]+1.732050807568877*alpha[2]*fr[14]+1.732050807568877*alpha[1]*fr[13]-1.0*alpha[5]*fr[12]+1.732050807568877*alpha[0]*fr[10]-1.0*alpha[2]*fr[9]-1.0*alpha[1]*fr[8]-1.0*alpha[0]*fr[4])*dfac_v; 
  incr[5] = -0.125*(1.732050807568877*alpha[0]*fr[11]+1.732050807568877*alpha[1]*fr[7]+1.732050807568877*alpha[2]*fr[6]-1.0*alpha[0]*fr[5]+1.732050807568877*fr[3]*alpha[5]-1.0*fr[0]*alpha[5]-1.0*alpha[1]*fr[2]-1.0*fr[1]*alpha[2])*dfac_v; 
  incr[6] = 0.125*(3.0*alpha[2]*fr[11]+3.0*alpha[5]*fr[7]+3.0*alpha[0]*fr[6]-1.732050807568877*alpha[2]*fr[5]-1.732050807568877*fr[2]*alpha[5]+3.0*alpha[1]*fr[3]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 
  incr[7] = 0.125*(3.0*alpha[1]*fr[11]+3.0*alpha[0]*fr[7]+3.0*alpha[5]*fr[6]-1.732050807568877*alpha[1]*fr[5]-1.732050807568877*fr[1]*alpha[5]+3.0*alpha[2]*fr[3]-1.732050807568877*alpha[0]*fr[2]-1.732050807568877*fr[0]*alpha[2])*dfac_v; 
  incr[8] = -0.125*(1.732050807568877*alpha[2]*fr[15]+1.732050807568877*alpha[5]*fr[14]+1.732050807568877*alpha[0]*fr[13]-1.0*alpha[2]*fr[12]+1.732050807568877*alpha[1]*fr[10]-1.0*alpha[5]*fr[9]-1.0*alpha[0]*fr[8]-1.0*alpha[1]*fr[4])*dfac_v; 
  incr[9] = -0.125*(1.732050807568877*alpha[1]*fr[15]+1.732050807568877*alpha[0]*fr[14]+1.732050807568877*alpha[5]*fr[13]-1.0*alpha[1]*fr[12]+1.732050807568877*alpha[2]*fr[10]-1.0*alpha[0]*fr[9]-1.0*alpha[5]*fr[8]-1.0*alpha[2]*fr[4])*dfac_v; 
  incr[10] = 0.125*(3.0*alpha[5]*fr[15]+3.0*alpha[2]*fr[14]+3.0*alpha[1]*fr[13]-1.732050807568877*alpha[5]*fr[12]+3.0*alpha[0]*fr[10]-1.732050807568877*alpha[2]*fr[9]-1.732050807568877*alpha[1]*fr[8]-1.732050807568877*alpha[0]*fr[4])*dfac_v; 
  incr[11] = 0.125*(3.0*alpha[0]*fr[11]+3.0*alpha[1]*fr[7]+3.0*alpha[2]*fr[6]-1.732050807568877*alpha[0]*fr[5]+3.0*fr[3]*alpha[5]-1.732050807568877*fr[0]*alpha[5]-1.732050807568877*alpha[1]*fr[2]-1.732050807568877*fr[1]*alpha[2])*dfac_v; 
  incr[12] = -0.125*(1.732050807568877*alpha[0]*fr[15]+1.732050807568877*alpha[1]*fr[14]+1.732050807568877*alpha[2]*fr[13]-1.0*alpha[0]*fr[12]+1.732050807568877*alpha[5]*fr[10]-1.0*alpha[1]*fr[9]-1.0*alpha[2]*fr[8]-1.0*fr[4]*alpha[5])*dfac_v; 
  incr[13] = 0.125*(3.0*alpha[2]*fr[15]+3.0*alpha[5]*fr[14]+3.0*alpha[0]*fr[13]-1.732050807568877*alpha[2]*fr[12]+3.0*alpha[1]*fr[10]-1.732050807568877*alpha[5]*fr[9]-1.732050807568877*alpha[0]*fr[8]-1.732050807568877*alpha[1]*fr[4])*dfac_v; 
  incr[14] = 0.125*(3.0*alpha[1]*fr[15]+3.0*alpha[0]*fr[14]+3.0*alpha[5]*fr[13]-1.732050807568877*alpha[1]*fr[12]+3.0*alpha[2]*fr[10]-1.732050807568877*alpha[0]*fr[9]-1.732050807568877*alpha[5]*fr[8]-1.732050807568877*alpha[2]*fr[4])*dfac_v; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]+3.0*alpha[1]*fr[14]+3.0*alpha[2]*fr[13]-1.732050807568877*alpha[0]*fr[12]+3.0*alpha[5]*fr[10]-1.732050807568877*alpha[1]*fr[9]-1.732050807568877*alpha[2]*fr[8]-1.732050807568877*fr[4]*alpha[5])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf2x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.5*((3.464101615137754*BdriftX[1]-2.0*BdriftX[0])*m_*wv2+((5.196152422706631*BmagInv[1]-3.0*BmagInv[0])*Phi[3]+(1.732050807568877*BmagInv[0]-3.0*BmagInv[1])*Phi[2])*dfac_y*q_))/q_; 

  double alpha[16]; 
  alpha[0] = (-(3.464101615137754*BdriftX[1]*m_*wv2)/q_)+(2.0*BdriftX[0]*m_*wv2)/q_-5.196152422706631*BmagInv[1]*Phi[3]*dfac_y+3.0*BmagInv[0]*Phi[3]*dfac_y+3.0*BmagInv[1]*Phi[2]*dfac_y-1.732050807568877*BmagInv[0]*Phi[2]*dfac_y; 
  alpha[3] = (1.154700538379251*BdriftX[0]*m_*wv)/(dfac_v*q_)-(2.0*BdriftX[1]*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[3]*fl[6]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[1]+alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.125*(3.0*alpha[3]*fl[6]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.125*(1.732050807568877*alpha[3]*fl[11]+alpha[3]*fl[7]+1.732050807568877*alpha[0]*fl[5]+alpha[0]*fl[2])*dfac_x; 
  incr[3] = 0.125*(1.732050807568877*alpha[0]*fl[6]+alpha[0]*fl[3]+1.732050807568877*fl[1]*alpha[3]+fl[0]*alpha[3])*dfac_x; 
  incr[4] = 0.125*(1.732050807568877*alpha[3]*fl[13]+alpha[3]*fl[10]+1.732050807568877*alpha[0]*fl[8]+alpha[0]*fl[4])*dfac_x; 
  incr[5] = -0.125*(3.0*alpha[3]*fl[11]+1.732050807568877*alpha[3]*fl[7]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[2])*dfac_x; 
  incr[6] = -0.125*(3.0*alpha[0]*fl[6]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[1]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_x; 
  incr[7] = 0.125*(1.732050807568877*alpha[0]*fl[11]+alpha[0]*fl[7]+1.732050807568877*alpha[3]*fl[5]+fl[2]*alpha[3])*dfac_x; 
  incr[8] = -0.125*(3.0*alpha[3]*fl[13]+1.732050807568877*alpha[3]*fl[10]+3.0*alpha[0]*fl[8]+1.732050807568877*alpha[0]*fl[4])*dfac_x; 
  incr[9] = 0.125*(1.732050807568877*alpha[3]*fl[15]+alpha[3]*fl[14]+1.732050807568877*alpha[0]*fl[12]+alpha[0]*fl[9])*dfac_x; 
  incr[10] = 0.125*(1.732050807568877*alpha[0]*fl[13]+alpha[0]*fl[10]+1.732050807568877*alpha[3]*fl[8]+alpha[3]*fl[4])*dfac_x; 
  incr[11] = -0.125*(3.0*alpha[0]*fl[11]+1.732050807568877*alpha[0]*fl[7]+3.0*alpha[3]*fl[5]+1.732050807568877*fl[2]*alpha[3])*dfac_x; 
  incr[12] = -0.125*(3.0*alpha[3]*fl[15]+1.732050807568877*alpha[3]*fl[14]+3.0*alpha[0]*fl[12]+1.732050807568877*alpha[0]*fl[9])*dfac_x; 
  incr[13] = -0.125*(3.0*alpha[0]*fl[13]+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[3]*fl[8]+1.732050807568877*alpha[3]*fl[4])*dfac_x; 
  incr[14] = 0.125*(1.732050807568877*alpha[0]*fl[15]+alpha[0]*fl[14]+1.732050807568877*alpha[3]*fl[12]+alpha[3]*fl[9])*dfac_x; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+1.732050807568877*alpha[0]*fl[14]+3.0*alpha[3]*fl[12]+1.732050807568877*alpha[3]*fl[9])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[3]*fr[6]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[1]-1.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.125*(3.0*alpha[3]*fr[6]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.125*(1.732050807568877*alpha[3]*fr[11]-1.0*alpha[3]*fr[7]+1.732050807568877*alpha[0]*fr[5]-1.0*alpha[0]*fr[2])*dfac_x; 
  incr[3] = -0.125*(1.732050807568877*alpha[0]*fr[6]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[1]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_x; 
  incr[4] = -0.125*(1.732050807568877*alpha[3]*fr[13]-1.0*alpha[3]*fr[10]+1.732050807568877*alpha[0]*fr[8]-1.0*alpha[0]*fr[4])*dfac_x; 
  incr[5] = 0.125*(3.0*alpha[3]*fr[11]-1.732050807568877*alpha[3]*fr[7]+3.0*alpha[0]*fr[5]-1.732050807568877*alpha[0]*fr[2])*dfac_x; 
  incr[6] = 0.125*(3.0*alpha[0]*fr[6]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[1]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_x; 
  incr[7] = -0.125*(1.732050807568877*alpha[0]*fr[11]-1.0*alpha[0]*fr[7]+1.732050807568877*alpha[3]*fr[5]-1.0*fr[2]*alpha[3])*dfac_x; 
  incr[8] = 0.125*(3.0*alpha[3]*fr[13]-1.732050807568877*alpha[3]*fr[10]+3.0*alpha[0]*fr[8]-1.732050807568877*alpha[0]*fr[4])*dfac_x; 
  incr[9] = -0.125*(1.732050807568877*alpha[3]*fr[15]-1.0*alpha[3]*fr[14]+1.732050807568877*alpha[0]*fr[12]-1.0*alpha[0]*fr[9])*dfac_x; 
  incr[10] = -0.125*(1.732050807568877*alpha[0]*fr[13]-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[3]*fr[8]-1.0*alpha[3]*fr[4])*dfac_x; 
  incr[11] = 0.125*(3.0*alpha[0]*fr[11]-1.732050807568877*alpha[0]*fr[7]+3.0*alpha[3]*fr[5]-1.732050807568877*fr[2]*alpha[3])*dfac_x; 
  incr[12] = 0.125*(3.0*alpha[3]*fr[15]-1.732050807568877*alpha[3]*fr[14]+3.0*alpha[0]*fr[12]-1.732050807568877*alpha[0]*fr[9])*dfac_x; 
  incr[13] = 0.125*(3.0*alpha[0]*fr[13]-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[3]*fr[8]-1.732050807568877*alpha[3]*fr[4])*dfac_x; 
  incr[14] = -0.125*(1.732050807568877*alpha[0]*fr[15]-1.0*alpha[0]*fr[14]+1.732050807568877*alpha[3]*fr[12]-1.0*alpha[3]*fr[9])*dfac_x; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]-1.732050807568877*alpha[0]*fr[14]+3.0*alpha[3]*fr[12]-1.732050807568877*alpha[3]*fr[9])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf2x2vSer_Y_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.5*(2.0*BdriftY[0]*m_*wv2+1.732050807568877*BmagInv[0]*Bmag[1]*dfac_x*wm+(1.732050807568877*BmagInv[0]*Phi[1]-3.0*BmagInv[0]*Phi[3])*dfac_x*q_))/q_; 

  double alpha[16]; 
  alpha[0] = (2.0*BdriftY[0]*m_*wv2)/q_+(1.732050807568877*BmagInv[0]*Bmag[1]*dfac_x*wm)/q_-3.0*BmagInv[0]*Phi[3]*dfac_x+1.732050807568877*BmagInv[0]*Phi[1]*dfac_x; 
  alpha[1] = (2.0*BdriftY[1]*m_*wv2)/q_+(1.732050807568877*Bmag[1]*BmagInv[1]*dfac_x*wm)/q_-3.0*BmagInv[1]*Phi[3]*dfac_x+1.732050807568877*BmagInv[1]*Phi[1]*dfac_x; 
  alpha[3] = (1.154700538379252*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alpha[4] = (BmagInv[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alpha[6] = (1.154700538379252*BdriftY[1]*m_*wv)/(dfac_v*q_); 
  alpha[8] = (Bmag[1]*BmagInv[1]*dfac_x)/(dfac_m*q_); 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[8]*fl[12]+1.732050807568877*alpha[6]*fl[11]+1.732050807568877*alpha[4]*fl[9]+alpha[8]*fl[8]+1.732050807568877*alpha[3]*fl[7]+alpha[6]*fl[6]+1.732050807568877*alpha[1]*fl[5]+alpha[4]*fl[4]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_y; 
  incr[1] = 0.125*(1.732050807568877*alpha[4]*fl[12]+1.732050807568877*alpha[3]*fl[11]+1.732050807568877*alpha[8]*fl[9]+alpha[4]*fl[8]+fl[4]*alpha[8]+1.732050807568877*alpha[6]*fl[7]+alpha[3]*fl[6]+fl[3]*alpha[6]+1.732050807568877*alpha[0]*fl[5]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_y; 
  incr[2] = -0.125*(3.0*alpha[8]*fl[12]+3.0*alpha[6]*fl[11]+3.0*alpha[4]*fl[9]+1.732050807568877*alpha[8]*fl[8]+3.0*alpha[3]*fl[7]+1.732050807568877*alpha[6]*fl[6]+3.0*alpha[1]*fl[5]+1.732050807568877*alpha[4]*fl[4]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_y; 
  incr[3] = 0.125*(1.732050807568877*alpha[8]*fl[15]+1.732050807568877*alpha[4]*fl[14]+alpha[8]*fl[13]+1.732050807568877*alpha[1]*fl[11]+alpha[4]*fl[10]+1.732050807568877*alpha[0]*fl[7]+alpha[1]*fl[6]+1.732050807568877*fl[5]*alpha[6]+fl[1]*alpha[6]+alpha[0]*fl[3]+1.732050807568877*fl[2]*alpha[3]+fl[0]*alpha[3])*dfac_y; 
  incr[4] = 0.125*(1.732050807568877*alpha[6]*fl[15]+1.732050807568877*alpha[3]*fl[14]+alpha[6]*fl[13]+1.732050807568877*alpha[1]*fl[12]+alpha[3]*fl[10]+1.732050807568877*alpha[0]*fl[9]+alpha[1]*fl[8]+1.732050807568877*fl[5]*alpha[8]+fl[1]*alpha[8]+alpha[0]*fl[4]+1.732050807568877*fl[2]*alpha[4]+fl[0]*alpha[4])*dfac_y; 
  incr[5] = -0.125*(3.0*alpha[4]*fl[12]+3.0*alpha[3]*fl[11]+3.0*alpha[8]*fl[9]+1.732050807568877*alpha[4]*fl[8]+1.732050807568877*fl[4]*alpha[8]+3.0*alpha[6]*fl[7]+1.732050807568877*alpha[3]*fl[6]+1.732050807568877*fl[3]*alpha[6]+3.0*alpha[0]*fl[5]+3.0*alpha[1]*fl[2]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_y; 
  incr[6] = 0.125*(1.732050807568877*alpha[4]*fl[15]+1.732050807568877*alpha[8]*fl[14]+alpha[4]*fl[13]+1.732050807568877*alpha[0]*fl[11]+alpha[8]*fl[10]+1.732050807568877*alpha[1]*fl[7]+alpha[0]*fl[6]+1.732050807568877*fl[2]*alpha[6]+fl[0]*alpha[6]+1.732050807568877*alpha[3]*fl[5]+alpha[1]*fl[3]+fl[1]*alpha[3])*dfac_y; 
  incr[7] = -0.125*(3.0*alpha[8]*fl[15]+3.0*alpha[4]*fl[14]+1.732050807568877*alpha[8]*fl[13]+3.0*alpha[1]*fl[11]+1.732050807568877*alpha[4]*fl[10]+3.0*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[6]+3.0*fl[5]*alpha[6]+1.732050807568877*fl[1]*alpha[6]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_y; 
  incr[8] = 0.125*(1.732050807568877*alpha[3]*fl[15]+1.732050807568877*alpha[6]*fl[14]+alpha[3]*fl[13]+1.732050807568877*alpha[0]*fl[12]+alpha[6]*fl[10]+1.732050807568877*alpha[1]*fl[9]+alpha[0]*fl[8]+1.732050807568877*fl[2]*alpha[8]+fl[0]*alpha[8]+1.732050807568877*alpha[4]*fl[5]+alpha[1]*fl[4]+fl[1]*alpha[4])*dfac_y; 
  incr[9] = -0.125*(3.0*alpha[6]*fl[15]+3.0*alpha[3]*fl[14]+1.732050807568877*alpha[6]*fl[13]+3.0*alpha[1]*fl[12]+1.732050807568877*alpha[3]*fl[10]+3.0*alpha[0]*fl[9]+1.732050807568877*alpha[1]*fl[8]+3.0*fl[5]*alpha[8]+1.732050807568877*fl[1]*alpha[8]+1.732050807568877*alpha[0]*fl[4]+3.0*fl[2]*alpha[4]+1.732050807568877*fl[0]*alpha[4])*dfac_y; 
  incr[10] = 0.125*(1.732050807568877*alpha[1]*fl[15]+1.732050807568877*alpha[0]*fl[14]+alpha[1]*fl[13]+1.732050807568877*alpha[6]*fl[12]+1.732050807568877*alpha[8]*fl[11]+alpha[0]*fl[10]+1.732050807568877*alpha[3]*fl[9]+alpha[6]*fl[8]+fl[6]*alpha[8]+1.732050807568877*alpha[4]*fl[7]+alpha[3]*fl[4]+fl[3]*alpha[4])*dfac_y; 
  incr[11] = -0.125*(3.0*alpha[4]*fl[15]+3.0*alpha[8]*fl[14]+1.732050807568877*alpha[4]*fl[13]+3.0*alpha[0]*fl[11]+1.732050807568877*alpha[8]*fl[10]+3.0*alpha[1]*fl[7]+1.732050807568877*alpha[0]*fl[6]+3.0*fl[2]*alpha[6]+1.732050807568877*fl[0]*alpha[6]+3.0*alpha[3]*fl[5]+1.732050807568877*alpha[1]*fl[3]+1.732050807568877*fl[1]*alpha[3])*dfac_y; 
  incr[12] = -0.125*(3.0*alpha[3]*fl[15]+3.0*alpha[6]*fl[14]+1.732050807568877*alpha[3]*fl[13]+3.0*alpha[0]*fl[12]+1.732050807568877*alpha[6]*fl[10]+3.0*alpha[1]*fl[9]+1.732050807568877*alpha[0]*fl[8]+3.0*fl[2]*alpha[8]+1.732050807568877*fl[0]*alpha[8]+3.0*alpha[4]*fl[5]+1.732050807568877*alpha[1]*fl[4]+1.732050807568877*fl[1]*alpha[4])*dfac_y; 
  incr[13] = 0.125*(1.732050807568877*alpha[0]*fl[15]+1.732050807568877*alpha[1]*fl[14]+alpha[0]*fl[13]+1.732050807568877*alpha[3]*fl[12]+1.732050807568877*alpha[4]*fl[11]+alpha[1]*fl[10]+1.732050807568877*alpha[6]*fl[9]+alpha[3]*fl[8]+1.732050807568877*fl[7]*alpha[8]+fl[3]*alpha[8]+alpha[4]*fl[6]+fl[4]*alpha[6])*dfac_y; 
  incr[14] = -0.125*(3.0*alpha[1]*fl[15]+3.0*alpha[0]*fl[14]+1.732050807568877*alpha[1]*fl[13]+3.0*alpha[6]*fl[12]+3.0*alpha[8]*fl[11]+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[3]*fl[9]+1.732050807568877*alpha[6]*fl[8]+1.732050807568877*fl[6]*alpha[8]+3.0*alpha[4]*fl[7]+1.732050807568877*alpha[3]*fl[4]+1.732050807568877*fl[3]*alpha[4])*dfac_y; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+3.0*alpha[1]*fl[14]+1.732050807568877*alpha[0]*fl[13]+3.0*alpha[3]*fl[12]+3.0*alpha[4]*fl[11]+1.732050807568877*alpha[1]*fl[10]+3.0*alpha[6]*fl[9]+1.732050807568877*alpha[3]*fl[8]+3.0*fl[7]*alpha[8]+1.732050807568877*fl[3]*alpha[8]+1.732050807568877*alpha[4]*fl[6]+1.732050807568877*fl[4]*alpha[6])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[8]*fr[12]+1.732050807568877*alpha[6]*fr[11]+1.732050807568877*alpha[4]*fr[9]-1.0*alpha[8]*fr[8]+1.732050807568877*alpha[3]*fr[7]-1.0*alpha[6]*fr[6]+1.732050807568877*alpha[1]*fr[5]-1.0*alpha[4]*fr[4]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_y; 
  incr[1] = -0.125*(1.732050807568877*alpha[4]*fr[12]+1.732050807568877*alpha[3]*fr[11]+1.732050807568877*alpha[8]*fr[9]-1.0*alpha[4]*fr[8]-1.0*fr[4]*alpha[8]+1.732050807568877*alpha[6]*fr[7]-1.0*alpha[3]*fr[6]-1.0*fr[3]*alpha[6]+1.732050807568877*alpha[0]*fr[5]+1.732050807568877*alpha[1]*fr[2]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_y; 
  incr[2] = 0.125*(3.0*alpha[8]*fr[12]+3.0*alpha[6]*fr[11]+3.0*alpha[4]*fr[9]-1.732050807568877*alpha[8]*fr[8]+3.0*alpha[3]*fr[7]-1.732050807568877*alpha[6]*fr[6]+3.0*alpha[1]*fr[5]-1.732050807568877*alpha[4]*fr[4]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_y; 
  incr[3] = -0.125*(1.732050807568877*alpha[8]*fr[15]+1.732050807568877*alpha[4]*fr[14]-1.0*alpha[8]*fr[13]+1.732050807568877*alpha[1]*fr[11]-1.0*alpha[4]*fr[10]+1.732050807568877*alpha[0]*fr[7]-1.0*alpha[1]*fr[6]+1.732050807568877*fr[5]*alpha[6]-1.0*fr[1]*alpha[6]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_y; 
  incr[4] = -0.125*(1.732050807568877*alpha[6]*fr[15]+1.732050807568877*alpha[3]*fr[14]-1.0*alpha[6]*fr[13]+1.732050807568877*alpha[1]*fr[12]-1.0*alpha[3]*fr[10]+1.732050807568877*alpha[0]*fr[9]-1.0*alpha[1]*fr[8]+1.732050807568877*fr[5]*alpha[8]-1.0*fr[1]*alpha[8]-1.0*alpha[0]*fr[4]+1.732050807568877*fr[2]*alpha[4]-1.0*fr[0]*alpha[4])*dfac_y; 
  incr[5] = 0.125*(3.0*alpha[4]*fr[12]+3.0*alpha[3]*fr[11]+3.0*alpha[8]*fr[9]-1.732050807568877*alpha[4]*fr[8]-1.732050807568877*fr[4]*alpha[8]+3.0*alpha[6]*fr[7]-1.732050807568877*alpha[3]*fr[6]-1.732050807568877*fr[3]*alpha[6]+3.0*alpha[0]*fr[5]+3.0*alpha[1]*fr[2]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_y; 
  incr[6] = -0.125*(1.732050807568877*alpha[4]*fr[15]+1.732050807568877*alpha[8]*fr[14]-1.0*alpha[4]*fr[13]+1.732050807568877*alpha[0]*fr[11]-1.0*alpha[8]*fr[10]+1.732050807568877*alpha[1]*fr[7]-1.0*alpha[0]*fr[6]+1.732050807568877*fr[2]*alpha[6]-1.0*fr[0]*alpha[6]+1.732050807568877*alpha[3]*fr[5]-1.0*alpha[1]*fr[3]-1.0*fr[1]*alpha[3])*dfac_y; 
  incr[7] = 0.125*(3.0*alpha[8]*fr[15]+3.0*alpha[4]*fr[14]-1.732050807568877*alpha[8]*fr[13]+3.0*alpha[1]*fr[11]-1.732050807568877*alpha[4]*fr[10]+3.0*alpha[0]*fr[7]-1.732050807568877*alpha[1]*fr[6]+3.0*fr[5]*alpha[6]-1.732050807568877*fr[1]*alpha[6]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_y; 
  incr[8] = -0.125*(1.732050807568877*alpha[3]*fr[15]+1.732050807568877*alpha[6]*fr[14]-1.0*alpha[3]*fr[13]+1.732050807568877*alpha[0]*fr[12]-1.0*alpha[6]*fr[10]+1.732050807568877*alpha[1]*fr[9]-1.0*alpha[0]*fr[8]+1.732050807568877*fr[2]*alpha[8]-1.0*fr[0]*alpha[8]+1.732050807568877*alpha[4]*fr[5]-1.0*alpha[1]*fr[4]-1.0*fr[1]*alpha[4])*dfac_y; 
  incr[9] = 0.125*(3.0*alpha[6]*fr[15]+3.0*alpha[3]*fr[14]-1.732050807568877*alpha[6]*fr[13]+3.0*alpha[1]*fr[12]-1.732050807568877*alpha[3]*fr[10]+3.0*alpha[0]*fr[9]-1.732050807568877*alpha[1]*fr[8]+3.0*fr[5]*alpha[8]-1.732050807568877*fr[1]*alpha[8]-1.732050807568877*alpha[0]*fr[4]+3.0*fr[2]*alpha[4]-1.732050807568877*fr[0]*alpha[4])*dfac_y; 
  incr[10] = -0.125*(1.732050807568877*alpha[1]*fr[15]+1.732050807568877*alpha[0]*fr[14]-1.0*alpha[1]*fr[13]+1.732050807568877*alpha[6]*fr[12]+1.732050807568877*alpha[8]*fr[11]-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[3]*fr[9]-1.0*alpha[6]*fr[8]-1.0*fr[6]*alpha[8]+1.732050807568877*alpha[4]*fr[7]-1.0*alpha[3]*fr[4]-1.0*fr[3]*alpha[4])*dfac_y; 
  incr[11] = 0.125*(3.0*alpha[4]*fr[15]+3.0*alpha[8]*fr[14]-1.732050807568877*alpha[4]*fr[13]+3.0*alpha[0]*fr[11]-1.732050807568877*alpha[8]*fr[10]+3.0*alpha[1]*fr[7]-1.732050807568877*alpha[0]*fr[6]+3.0*fr[2]*alpha[6]-1.732050807568877*fr[0]*alpha[6]+3.0*alpha[3]*fr[5]-1.732050807568877*alpha[1]*fr[3]-1.732050807568877*fr[1]*alpha[3])*dfac_y; 
  incr[12] = 0.125*(3.0*alpha[3]*fr[15]+3.0*alpha[6]*fr[14]-1.732050807568877*alpha[3]*fr[13]+3.0*alpha[0]*fr[12]-1.732050807568877*alpha[6]*fr[10]+3.0*alpha[1]*fr[9]-1.732050807568877*alpha[0]*fr[8]+3.0*fr[2]*alpha[8]-1.732050807568877*fr[0]*alpha[8]+3.0*alpha[4]*fr[5]-1.732050807568877*alpha[1]*fr[4]-1.732050807568877*fr[1]*alpha[4])*dfac_y; 
  incr[13] = -0.125*(1.732050807568877*alpha[0]*fr[15]+1.732050807568877*alpha[1]*fr[14]-1.0*alpha[0]*fr[13]+1.732050807568877*alpha[3]*fr[12]+1.732050807568877*alpha[4]*fr[11]-1.0*alpha[1]*fr[10]+1.732050807568877*alpha[6]*fr[9]-1.0*alpha[3]*fr[8]+1.732050807568877*fr[7]*alpha[8]-1.0*fr[3]*alpha[8]-1.0*alpha[4]*fr[6]-1.0*fr[4]*alpha[6])*dfac_y; 
  incr[14] = 0.125*(3.0*alpha[1]*fr[15]+3.0*alpha[0]*fr[14]-1.732050807568877*alpha[1]*fr[13]+3.0*alpha[6]*fr[12]+3.0*alpha[8]*fr[11]-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[3]*fr[9]-1.732050807568877*alpha[6]*fr[8]-1.732050807568877*fr[6]*alpha[8]+3.0*alpha[4]*fr[7]-1.732050807568877*alpha[3]*fr[4]-1.732050807568877*fr[3]*alpha[4])*dfac_y; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]+3.0*alpha[1]*fr[14]-1.732050807568877*alpha[0]*fr[13]+3.0*alpha[3]*fr[12]+3.0*alpha[4]*fr[11]-1.732050807568877*alpha[1]*fr[10]+3.0*alpha[6]*fr[9]-1.732050807568877*alpha[3]*fr[8]+3.0*fr[7]*alpha[8]-1.732050807568877*fr[3]*alpha[8]-1.732050807568877*alpha[4]*fr[6]-1.732050807568877*fr[4]*alpha[6])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf2x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.5*((1.732050807568877*BdriftX[0]*Bmag[1]*dfac_v*dfac_x*wm+((1.732050807568877*BdriftY[1]*Phi[3]+1.732050807568877*BdriftY[0]*Phi[2])*dfac_v*dfac_y+1.732050807568877*BdriftX[0]*Phi[1]*dfac_v*dfac_x)*q_)*wv-1.732050807568877*BdriftX[0]*Bmag[1]*dfac_x*wm+(((-1.732050807568877*BdriftY[1]*Phi[3])-1.732050807568877*BdriftY[0]*Phi[2])*dfac_y-1.732050807568877*BdriftX[0]*Phi[1]*dfac_x)*q_))/(dfac_v*q_); 

  double alpha[16]; 
  alpha[0] = (-(1.732050807568877*BdriftX[0]*Bmag[1]*dfac_x*wm*wv)/q_)-1.732050807568877*BdriftY[1]*Phi[3]*dfac_y*wv-1.732050807568877*BdriftY[0]*Phi[2]*dfac_y*wv-1.732050807568877*BdriftX[0]*Phi[1]*dfac_x*wv+(1.732050807568877*BdriftX[0]*Bmag[1]*dfac_x*wm)/(dfac_v*q_)+(1.732050807568877*BdriftY[1]*Phi[3]*dfac_y)/dfac_v+(1.732050807568877*BdriftY[0]*Phi[2]*dfac_y)/dfac_v+(1.732050807568877*BdriftX[0]*Phi[1]*dfac_x)/dfac_v; 
  alpha[1] = (-(1.732050807568877*BdriftX[1]*Bmag[1]*dfac_x*wm*wv)/q_)-1.732050807568877*BdriftY[0]*Phi[3]*dfac_y*wv-1.732050807568877*BdriftY[1]*Phi[2]*dfac_y*wv-1.732050807568877*BdriftX[1]*Phi[1]*dfac_x*wv+(1.732050807568877*BdriftX[1]*Bmag[1]*dfac_x*wm)/(dfac_v*q_)+(1.732050807568877*BdriftY[0]*Phi[3]*dfac_y)/dfac_v+(1.732050807568877*BdriftY[1]*Phi[2]*dfac_y)/dfac_v+(1.732050807568877*BdriftX[1]*Phi[1]*dfac_x)/dfac_v; 
  alpha[2] = (1.732050807568877*BdriftX[0]*Phi[3]*dfac_x)/dfac_v-1.732050807568877*BdriftX[0]*Phi[3]*dfac_x*wv; 
  alpha[4] = (1.0*BdriftX[0]*Bmag[1]*dfac_x)/(dfac_m*dfac_v*q_)-(1.0*BdriftX[0]*Bmag[1]*dfac_x*wv)/(dfac_m*q_); 
  alpha[5] = (1.732050807568877*BdriftX[1]*Phi[3]*dfac_x)/dfac_v-1.732050807568877*BdriftX[1]*Phi[3]*dfac_x*wv; 
  alpha[8] = (1.0*BdriftX[1]*Bmag[1]*dfac_x)/(dfac_m*dfac_v*q_)-(1.0*BdriftX[1]*Bmag[1]*dfac_x*wv)/(dfac_m*q_); 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[8]*fl[13]+1.732050807568877*alpha[5]*fl[11]+1.732050807568877*alpha[4]*fl[10]+alpha[8]*fl[8]+1.732050807568877*alpha[2]*fl[7]+1.732050807568877*alpha[1]*fl[6]+alpha[5]*fl[5]+alpha[4]*fl[4]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.125*(1.732050807568877*alpha[4]*fl[13]+1.732050807568877*alpha[2]*fl[11]+1.732050807568877*alpha[8]*fl[10]+alpha[4]*fl[8]+fl[4]*alpha[8]+1.732050807568877*alpha[5]*fl[7]+1.732050807568877*alpha[0]*fl[6]+alpha[2]*fl[5]+fl[2]*alpha[5]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.125*(1.732050807568877*alpha[8]*fl[15]+1.732050807568877*alpha[4]*fl[14]+alpha[8]*fl[12]+1.732050807568877*alpha[1]*fl[11]+alpha[4]*fl[9]+1.732050807568877*alpha[0]*fl[7]+1.732050807568877*alpha[5]*fl[6]+alpha[1]*fl[5]+fl[1]*alpha[5]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.125*(3.0*alpha[8]*fl[13]+3.0*alpha[5]*fl[11]+3.0*alpha[4]*fl[10]+1.732050807568877*alpha[8]*fl[8]+3.0*alpha[2]*fl[7]+3.0*alpha[1]*fl[6]+1.732050807568877*alpha[5]*fl[5]+1.732050807568877*alpha[4]*fl[4]+3.0*alpha[0]*fl[3]+1.732050807568877*alpha[2]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[4] = 0.125*(1.732050807568877*alpha[5]*fl[15]+1.732050807568877*alpha[2]*fl[14]+1.732050807568877*alpha[1]*fl[13]+alpha[5]*fl[12]+1.732050807568877*alpha[0]*fl[10]+alpha[2]*fl[9]+alpha[1]*fl[8]+1.732050807568877*fl[6]*alpha[8]+fl[1]*alpha[8]+alpha[0]*fl[4]+1.732050807568877*fl[3]*alpha[4]+fl[0]*alpha[4])*dfac_v; 
  incr[5] = 0.125*(1.732050807568877*alpha[4]*fl[15]+1.732050807568877*alpha[8]*fl[14]+alpha[4]*fl[12]+1.732050807568877*alpha[0]*fl[11]+alpha[8]*fl[9]+1.732050807568877*alpha[1]*fl[7]+1.732050807568877*alpha[2]*fl[6]+alpha[0]*fl[5]+1.732050807568877*fl[3]*alpha[5]+fl[0]*alpha[5]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.125*(3.0*alpha[4]*fl[13]+3.0*alpha[2]*fl[11]+3.0*alpha[8]*fl[10]+1.732050807568877*alpha[4]*fl[8]+1.732050807568877*fl[4]*alpha[8]+3.0*alpha[5]*fl[7]+3.0*alpha[0]*fl[6]+1.732050807568877*alpha[2]*fl[5]+1.732050807568877*fl[2]*alpha[5]+3.0*alpha[1]*fl[3]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 
  incr[7] = -0.125*(3.0*alpha[8]*fl[15]+3.0*alpha[4]*fl[14]+1.732050807568877*alpha[8]*fl[12]+3.0*alpha[1]*fl[11]+1.732050807568877*alpha[4]*fl[9]+3.0*alpha[0]*fl[7]+3.0*alpha[5]*fl[6]+1.732050807568877*alpha[1]*fl[5]+1.732050807568877*fl[1]*alpha[5]+3.0*alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+1.732050807568877*fl[0]*alpha[2])*dfac_v; 
  incr[8] = 0.125*(1.732050807568877*alpha[2]*fl[15]+1.732050807568877*alpha[5]*fl[14]+1.732050807568877*alpha[0]*fl[13]+alpha[2]*fl[12]+1.732050807568877*alpha[1]*fl[10]+alpha[5]*fl[9]+alpha[0]*fl[8]+1.732050807568877*fl[3]*alpha[8]+fl[0]*alpha[8]+1.732050807568877*alpha[4]*fl[6]+alpha[1]*fl[4]+fl[1]*alpha[4])*dfac_v; 
  incr[9] = 0.125*(1.732050807568877*alpha[1]*fl[15]+1.732050807568877*alpha[0]*fl[14]+1.732050807568877*alpha[5]*fl[13]+alpha[1]*fl[12]+1.732050807568877*alpha[8]*fl[11]+1.732050807568877*alpha[2]*fl[10]+alpha[0]*fl[9]+alpha[5]*fl[8]+fl[5]*alpha[8]+1.732050807568877*alpha[4]*fl[7]+alpha[2]*fl[4]+fl[2]*alpha[4])*dfac_v; 
  incr[10] = -0.125*(3.0*alpha[5]*fl[15]+3.0*alpha[2]*fl[14]+3.0*alpha[1]*fl[13]+1.732050807568877*alpha[5]*fl[12]+3.0*alpha[0]*fl[10]+1.732050807568877*alpha[2]*fl[9]+1.732050807568877*alpha[1]*fl[8]+3.0*fl[6]*alpha[8]+1.732050807568877*fl[1]*alpha[8]+1.732050807568877*alpha[0]*fl[4]+3.0*fl[3]*alpha[4]+1.732050807568877*fl[0]*alpha[4])*dfac_v; 
  incr[11] = -0.125*(3.0*alpha[4]*fl[15]+3.0*alpha[8]*fl[14]+1.732050807568877*alpha[4]*fl[12]+3.0*alpha[0]*fl[11]+1.732050807568877*alpha[8]*fl[9]+3.0*alpha[1]*fl[7]+3.0*alpha[2]*fl[6]+1.732050807568877*alpha[0]*fl[5]+3.0*fl[3]*alpha[5]+1.732050807568877*fl[0]*alpha[5]+1.732050807568877*alpha[1]*fl[2]+1.732050807568877*fl[1]*alpha[2])*dfac_v; 
  incr[12] = 0.125*(1.732050807568877*alpha[0]*fl[15]+1.732050807568877*alpha[1]*fl[14]+1.732050807568877*alpha[2]*fl[13]+alpha[0]*fl[12]+1.732050807568877*alpha[4]*fl[11]+1.732050807568877*alpha[5]*fl[10]+alpha[1]*fl[9]+alpha[2]*fl[8]+1.732050807568877*fl[7]*alpha[8]+fl[2]*alpha[8]+alpha[4]*fl[5]+fl[4]*alpha[5])*dfac_v; 
  incr[13] = -0.125*(3.0*alpha[2]*fl[15]+3.0*alpha[5]*fl[14]+3.0*alpha[0]*fl[13]+1.732050807568877*alpha[2]*fl[12]+3.0*alpha[1]*fl[10]+1.732050807568877*alpha[5]*fl[9]+1.732050807568877*alpha[0]*fl[8]+3.0*fl[3]*alpha[8]+1.732050807568877*fl[0]*alpha[8]+3.0*alpha[4]*fl[6]+1.732050807568877*alpha[1]*fl[4]+1.732050807568877*fl[1]*alpha[4])*dfac_v; 
  incr[14] = -0.125*(3.0*alpha[1]*fl[15]+3.0*alpha[0]*fl[14]+3.0*alpha[5]*fl[13]+1.732050807568877*alpha[1]*fl[12]+3.0*alpha[8]*fl[11]+3.0*alpha[2]*fl[10]+1.732050807568877*alpha[0]*fl[9]+1.732050807568877*alpha[5]*fl[8]+1.732050807568877*fl[5]*alpha[8]+3.0*alpha[4]*fl[7]+1.732050807568877*alpha[2]*fl[4]+1.732050807568877*fl[2]*alpha[4])*dfac_v; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+3.0*alpha[1]*fl[14]+3.0*alpha[2]*fl[13]+1.732050807568877*alpha[0]*fl[12]+3.0*alpha[4]*fl[11]+3.0*alpha[5]*fl[10]+1.732050807568877*alpha[1]*fl[9]+1.732050807568877*alpha[2]*fl[8]+3.0*fl[7]*alpha[8]+1.732050807568877*fl[2]*alpha[8]+1.732050807568877*alpha[4]*fl[5]+1.732050807568877*fl[4]*alpha[5])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[8]*fr[13]+1.732050807568877*alpha[5]*fr[11]+1.732050807568877*alpha[4]*fr[10]-1.0*alpha[8]*fr[8]+1.732050807568877*alpha[2]*fr[7]+1.732050807568877*alpha[1]*fr[6]-1.0*alpha[5]*fr[5]-1.0*alpha[4]*fr[4]+1.732050807568877*alpha[0]*fr[3]-1.0*alpha[2]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.125*(1.732050807568877*alpha[4]*fr[13]+1.732050807568877*alpha[2]*fr[11]+1.732050807568877*alpha[8]*fr[10]-1.0*alpha[4]*fr[8]-1.0*fr[4]*alpha[8]+1.732050807568877*alpha[5]*fr[7]+1.732050807568877*alpha[0]*fr[6]-1.0*alpha[2]*fr[5]-1.0*fr[2]*alpha[5]+1.732050807568877*alpha[1]*fr[3]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = -0.125*(1.732050807568877*alpha[8]*fr[15]+1.732050807568877*alpha[4]*fr[14]-1.0*alpha[8]*fr[12]+1.732050807568877*alpha[1]*fr[11]-1.0*alpha[4]*fr[9]+1.732050807568877*alpha[0]*fr[7]+1.732050807568877*alpha[5]*fr[6]-1.0*alpha[1]*fr[5]-1.0*fr[1]*alpha[5]+1.732050807568877*alpha[2]*fr[3]-1.0*alpha[0]*fr[2]-1.0*fr[0]*alpha[2])*dfac_v; 
  incr[3] = 0.125*(3.0*alpha[8]*fr[13]+3.0*alpha[5]*fr[11]+3.0*alpha[4]*fr[10]-1.732050807568877*alpha[8]*fr[8]+3.0*alpha[2]*fr[7]+3.0*alpha[1]*fr[6]-1.732050807568877*alpha[5]*fr[5]-1.732050807568877*alpha[4]*fr[4]+3.0*alpha[0]*fr[3]-1.732050807568877*alpha[2]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[4] = -0.125*(1.732050807568877*alpha[5]*fr[15]+1.732050807568877*alpha[2]*fr[14]+1.732050807568877*alpha[1]*fr[13]-1.0*alpha[5]*fr[12]+1.732050807568877*alpha[0]*fr[10]-1.0*alpha[2]*fr[9]-1.0*alpha[1]*fr[8]+1.732050807568877*fr[6]*alpha[8]-1.0*fr[1]*alpha[8]-1.0*alpha[0]*fr[4]+1.732050807568877*fr[3]*alpha[4]-1.0*fr[0]*alpha[4])*dfac_v; 
  incr[5] = -0.125*(1.732050807568877*alpha[4]*fr[15]+1.732050807568877*alpha[8]*fr[14]-1.0*alpha[4]*fr[12]+1.732050807568877*alpha[0]*fr[11]-1.0*alpha[8]*fr[9]+1.732050807568877*alpha[1]*fr[7]+1.732050807568877*alpha[2]*fr[6]-1.0*alpha[0]*fr[5]+1.732050807568877*fr[3]*alpha[5]-1.0*fr[0]*alpha[5]-1.0*alpha[1]*fr[2]-1.0*fr[1]*alpha[2])*dfac_v; 
  incr[6] = 0.125*(3.0*alpha[4]*fr[13]+3.0*alpha[2]*fr[11]+3.0*alpha[8]*fr[10]-1.732050807568877*alpha[4]*fr[8]-1.732050807568877*fr[4]*alpha[8]+3.0*alpha[5]*fr[7]+3.0*alpha[0]*fr[6]-1.732050807568877*alpha[2]*fr[5]-1.732050807568877*fr[2]*alpha[5]+3.0*alpha[1]*fr[3]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 
  incr[7] = 0.125*(3.0*alpha[8]*fr[15]+3.0*alpha[4]*fr[14]-1.732050807568877*alpha[8]*fr[12]+3.0*alpha[1]*fr[11]-1.732050807568877*alpha[4]*fr[9]+3.0*alpha[0]*fr[7]+3.0*alpha[5]*fr[6]-1.732050807568877*alpha[1]*fr[5]-1.732050807568877*fr[1]*alpha[5]+3.0*alpha[2]*fr[3]-1.732050807568877*alpha[0]*fr[2]-1.732050807568877*fr[0]*alpha[2])*dfac_v; 
  incr[8] = -0.125*(1.732050807568877*alpha[2]*fr[15]+1.732050807568877*alpha[5]*fr[14]+1.732050807568877*alpha[0]*fr[13]-1.0*alpha[2]*fr[12]+1.732050807568877*alpha[1]*fr[10]-1.0*alpha[5]*fr[9]-1.0*alpha[0]*fr[8]+1.732050807568877*fr[3]*alpha[8]-1.0*fr[0]*alpha[8]+1.732050807568877*alpha[4]*fr[6]-1.0*alpha[1]*fr[4]-1.0*fr[1]*alpha[4])*dfac_v; 
  incr[9] = -0.125*(1.732050807568877*alpha[1]*fr[15]+1.732050807568877*alpha[0]*fr[14]+1.732050807568877*alpha[5]*fr[13]-1.0*alpha[1]*fr[12]+1.732050807568877*alpha[8]*fr[11]+1.732050807568877*alpha[2]*fr[10]-1.0*alpha[0]*fr[9]-1.0*alpha[5]*fr[8]-1.0*fr[5]*alpha[8]+1.732050807568877*alpha[4]*fr[7]-1.0*alpha[2]*fr[4]-1.0*fr[2]*alpha[4])*dfac_v; 
  incr[10] = 0.125*(3.0*alpha[5]*fr[15]+3.0*alpha[2]*fr[14]+3.0*alpha[1]*fr[13]-1.732050807568877*alpha[5]*fr[12]+3.0*alpha[0]*fr[10]-1.732050807568877*alpha[2]*fr[9]-1.732050807568877*alpha[1]*fr[8]+3.0*fr[6]*alpha[8]-1.732050807568877*fr[1]*alpha[8]-1.732050807568877*alpha[0]*fr[4]+3.0*fr[3]*alpha[4]-1.732050807568877*fr[0]*alpha[4])*dfac_v; 
  incr[11] = 0.125*(3.0*alpha[4]*fr[15]+3.0*alpha[8]*fr[14]-1.732050807568877*alpha[4]*fr[12]+3.0*alpha[0]*fr[11]-1.732050807568877*alpha[8]*fr[9]+3.0*alpha[1]*fr[7]+3.0*alpha[2]*fr[6]-1.732050807568877*alpha[0]*fr[5]+3.0*fr[3]*alpha[5]-1.732050807568877*fr[0]*alpha[5]-1.732050807568877*alpha[1]*fr[2]-1.732050807568877*fr[1]*alpha[2])*dfac_v; 
  incr[12] = -0.125*(1.732050807568877*alpha[0]*fr[15]+1.732050807568877*alpha[1]*fr[14]+1.732050807568877*alpha[2]*fr[13]-1.0*alpha[0]*fr[12]+1.732050807568877*alpha[4]*fr[11]+1.732050807568877*alpha[5]*fr[10]-1.0*alpha[1]*fr[9]-1.0*alpha[2]*fr[8]+1.732050807568877*fr[7]*alpha[8]-1.0*fr[2]*alpha[8]-1.0*alpha[4]*fr[5]-1.0*fr[4]*alpha[5])*dfac_v; 
  incr[13] = 0.125*(3.0*alpha[2]*fr[15]+3.0*alpha[5]*fr[14]+3.0*alpha[0]*fr[13]-1.732050807568877*alpha[2]*fr[12]+3.0*alpha[1]*fr[10]-1.732050807568877*alpha[5]*fr[9]-1.732050807568877*alpha[0]*fr[8]+3.0*fr[3]*alpha[8]-1.732050807568877*fr[0]*alpha[8]+3.0*alpha[4]*fr[6]-1.732050807568877*alpha[1]*fr[4]-1.732050807568877*fr[1]*alpha[4])*dfac_v; 
  incr[14] = 0.125*(3.0*alpha[1]*fr[15]+3.0*alpha[0]*fr[14]+3.0*alpha[5]*fr[13]-1.732050807568877*alpha[1]*fr[12]+3.0*alpha[8]*fr[11]+3.0*alpha[2]*fr[10]-1.732050807568877*alpha[0]*fr[9]-1.732050807568877*alpha[5]*fr[8]-1.732050807568877*fr[5]*alpha[8]+3.0*alpha[4]*fr[7]-1.732050807568877*alpha[2]*fr[4]-1.732050807568877*fr[2]*alpha[4])*dfac_v; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]+3.0*alpha[1]*fr[14]+3.0*alpha[2]*fr[13]-1.732050807568877*alpha[0]*fr[12]+3.0*alpha[4]*fr[11]+3.0*alpha[5]*fr[10]-1.732050807568877*alpha[1]*fr[9]-1.732050807568877*alpha[2]*fr[8]+3.0*fr[7]*alpha[8]-1.732050807568877*fr[2]*alpha[8]-1.732050807568877*alpha[4]*fr[5]-1.732050807568877*fr[4]*alpha[5])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.5*((3.464101615137754*BdriftX[1]-2.0*BdriftX[0])*m_*wv2+((6.0*Apar[3]-3.464101615137754*Apar[2])*dfac_y+(1.732050807568877*Apar[0]-3.0*Apar[1])*BdriftX[1]+1.732050807568877*BdriftX[0]*Apar[1]-1.0*Apar[0]*BdriftX[0])*q_*wv+((5.196152422706631*BmagInv[1]-3.0*BmagInv[0])*Phi[3]+(1.732050807568877*BmagInv[0]-3.0*BmagInv[1])*Phi[2])*dfac_y*q_))/q_; 

  double alpha[16]; 
  alpha[0] = (-(3.464101615137754*BdriftX[1]*m_*wv2)/q_)+(2.0*BdriftX[0]*m_*wv2)/q_-6.0*Apar[3]*dfac_y*wv+3.464101615137754*Apar[2]*dfac_y*wv+3.0*Apar[1]*BdriftX[1]*wv-1.732050807568877*Apar[0]*BdriftX[1]*wv-1.732050807568877*BdriftX[0]*Apar[1]*wv+1.0*Apar[0]*BdriftX[0]*wv-5.196152422706631*BmagInv[1]*Phi[3]*dfac_y+3.0*BmagInv[0]*Phi[3]*dfac_y+3.0*BmagInv[1]*Phi[2]*dfac_y-1.732050807568877*BmagInv[0]*Phi[2]*dfac_y; 
  alpha[2] = 3.0*BdriftX[1]*Apar[3]*wv-1.732050807568878*BdriftX[0]*Apar[3]*wv-1.732050807568878*BdriftX[1]*Apar[2]*wv+1.0*BdriftX[0]*Apar[2]*wv; 
  alpha[3] = (1.154700538379251*BdriftX[0]*m_*wv)/(dfac_v*q_)-(2.0*BdriftX[1]*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[3]*fl[6]+1.732050807568877*alpha[2]*fl[5]+alpha[3]*fl[3]+alpha[2]*fl[2]+1.732050807568877*alpha[0]*fl[1]+alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.125*(3.0*alpha[3]*fl[6]+3.0*alpha[2]*fl[5]+1.732050807568877*alpha[3]*fl[3]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.125*(1.732050807568877*alpha[3]*fl[11]+alpha[3]*fl[7]+1.732050807568877*alpha[0]*fl[5]+alpha[0]*fl[2]+1.732050807568877*fl[1]*alpha[2]+fl[0]*alpha[2])*dfac_x; 
  incr[3] = 0.125*(1.732050807568877*alpha[2]*fl[11]+alpha[2]*fl[7]+1.732050807568877*alpha[0]*fl[6]+alpha[0]*fl[3]+1.732050807568877*fl[1]*alpha[3]+fl[0]*alpha[3])*dfac_x; 
  incr[4] = 0.125*(1.732050807568877*alpha[3]*fl[13]+1.732050807568877*alpha[2]*fl[12]+alpha[3]*fl[10]+alpha[2]*fl[9]+1.732050807568877*alpha[0]*fl[8]+alpha[0]*fl[4])*dfac_x; 
  incr[5] = -0.125*(3.0*alpha[3]*fl[11]+1.732050807568877*alpha[3]*fl[7]+3.0*alpha[0]*fl[5]+1.732050807568877*alpha[0]*fl[2]+3.0*fl[1]*alpha[2]+1.732050807568877*fl[0]*alpha[2])*dfac_x; 
  incr[6] = -0.125*(3.0*alpha[2]*fl[11]+1.732050807568877*alpha[2]*fl[7]+3.0*alpha[0]*fl[6]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[1]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_x; 
  incr[7] = 0.125*(1.732050807568877*alpha[0]*fl[11]+alpha[0]*fl[7]+1.732050807568877*alpha[2]*fl[6]+1.732050807568877*alpha[3]*fl[5]+alpha[2]*fl[3]+fl[2]*alpha[3])*dfac_x; 
  incr[8] = -0.125*(3.0*alpha[3]*fl[13]+3.0*alpha[2]*fl[12]+1.732050807568877*alpha[3]*fl[10]+1.732050807568877*alpha[2]*fl[9]+3.0*alpha[0]*fl[8]+1.732050807568877*alpha[0]*fl[4])*dfac_x; 
  incr[9] = 0.125*(1.732050807568877*alpha[3]*fl[15]+alpha[3]*fl[14]+1.732050807568877*alpha[0]*fl[12]+alpha[0]*fl[9]+1.732050807568877*alpha[2]*fl[8]+alpha[2]*fl[4])*dfac_x; 
  incr[10] = 0.125*(1.732050807568877*alpha[2]*fl[15]+alpha[2]*fl[14]+1.732050807568877*alpha[0]*fl[13]+alpha[0]*fl[10]+1.732050807568877*alpha[3]*fl[8]+alpha[3]*fl[4])*dfac_x; 
  incr[11] = -0.125*(3.0*alpha[0]*fl[11]+1.732050807568877*alpha[0]*fl[7]+3.0*alpha[2]*fl[6]+3.0*alpha[3]*fl[5]+1.732050807568877*alpha[2]*fl[3]+1.732050807568877*fl[2]*alpha[3])*dfac_x; 
  incr[12] = -0.125*(3.0*alpha[3]*fl[15]+1.732050807568877*alpha[3]*fl[14]+3.0*alpha[0]*fl[12]+1.732050807568877*alpha[0]*fl[9]+3.0*alpha[2]*fl[8]+1.732050807568877*alpha[2]*fl[4])*dfac_x; 
  incr[13] = -0.125*(3.0*alpha[2]*fl[15]+1.732050807568877*alpha[2]*fl[14]+3.0*alpha[0]*fl[13]+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[3]*fl[8]+1.732050807568877*alpha[3]*fl[4])*dfac_x; 
  incr[14] = 0.125*(1.732050807568877*alpha[0]*fl[15]+alpha[0]*fl[14]+1.732050807568877*alpha[2]*fl[13]+1.732050807568877*alpha[3]*fl[12]+alpha[2]*fl[10]+alpha[3]*fl[9])*dfac_x; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+1.732050807568877*alpha[0]*fl[14]+3.0*alpha[2]*fl[13]+3.0*alpha[3]*fl[12]+1.732050807568877*alpha[2]*fl[10]+1.732050807568877*alpha[3]*fl[9])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[3]*fr[6]+1.732050807568877*alpha[2]*fr[5]-1.0*alpha[3]*fr[3]-1.0*alpha[2]*fr[2]+1.732050807568877*alpha[0]*fr[1]-1.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.125*(3.0*alpha[3]*fr[6]+3.0*alpha[2]*fr[5]-1.732050807568877*alpha[3]*fr[3]-1.732050807568877*alpha[2]*fr[2]+3.0*alpha[0]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.125*(1.732050807568877*alpha[3]*fr[11]-1.0*alpha[3]*fr[7]+1.732050807568877*alpha[0]*fr[5]-1.0*alpha[0]*fr[2]+1.732050807568877*fr[1]*alpha[2]-1.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = -0.125*(1.732050807568877*alpha[2]*fr[11]-1.0*alpha[2]*fr[7]+1.732050807568877*alpha[0]*fr[6]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[1]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_x; 
  incr[4] = -0.125*(1.732050807568877*alpha[3]*fr[13]+1.732050807568877*alpha[2]*fr[12]-1.0*alpha[3]*fr[10]-1.0*alpha[2]*fr[9]+1.732050807568877*alpha[0]*fr[8]-1.0*alpha[0]*fr[4])*dfac_x; 
  incr[5] = 0.125*(3.0*alpha[3]*fr[11]-1.732050807568877*alpha[3]*fr[7]+3.0*alpha[0]*fr[5]-1.732050807568877*alpha[0]*fr[2]+3.0*fr[1]*alpha[2]-1.732050807568877*fr[0]*alpha[2])*dfac_x; 
  incr[6] = 0.125*(3.0*alpha[2]*fr[11]-1.732050807568877*alpha[2]*fr[7]+3.0*alpha[0]*fr[6]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[1]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_x; 
  incr[7] = -0.125*(1.732050807568877*alpha[0]*fr[11]-1.0*alpha[0]*fr[7]+1.732050807568877*alpha[2]*fr[6]+1.732050807568877*alpha[3]*fr[5]-1.0*alpha[2]*fr[3]-1.0*fr[2]*alpha[3])*dfac_x; 
  incr[8] = 0.125*(3.0*alpha[3]*fr[13]+3.0*alpha[2]*fr[12]-1.732050807568877*alpha[3]*fr[10]-1.732050807568877*alpha[2]*fr[9]+3.0*alpha[0]*fr[8]-1.732050807568877*alpha[0]*fr[4])*dfac_x; 
  incr[9] = -0.125*(1.732050807568877*alpha[3]*fr[15]-1.0*alpha[3]*fr[14]+1.732050807568877*alpha[0]*fr[12]-1.0*alpha[0]*fr[9]+1.732050807568877*alpha[2]*fr[8]-1.0*alpha[2]*fr[4])*dfac_x; 
  incr[10] = -0.125*(1.732050807568877*alpha[2]*fr[15]-1.0*alpha[2]*fr[14]+1.732050807568877*alpha[0]*fr[13]-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[3]*fr[8]-1.0*alpha[3]*fr[4])*dfac_x; 
  incr[11] = 0.125*(3.0*alpha[0]*fr[11]-1.732050807568877*alpha[0]*fr[7]+3.0*alpha[2]*fr[6]+3.0*alpha[3]*fr[5]-1.732050807568877*alpha[2]*fr[3]-1.732050807568877*fr[2]*alpha[3])*dfac_x; 
  incr[12] = 0.125*(3.0*alpha[3]*fr[15]-1.732050807568877*alpha[3]*fr[14]+3.0*alpha[0]*fr[12]-1.732050807568877*alpha[0]*fr[9]+3.0*alpha[2]*fr[8]-1.732050807568877*alpha[2]*fr[4])*dfac_x; 
  incr[13] = 0.125*(3.0*alpha[2]*fr[15]-1.732050807568877*alpha[2]*fr[14]+3.0*alpha[0]*fr[13]-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[3]*fr[8]-1.732050807568877*alpha[3]*fr[4])*dfac_x; 
  incr[14] = -0.125*(1.732050807568877*alpha[0]*fr[15]-1.0*alpha[0]*fr[14]+1.732050807568877*alpha[2]*fr[13]+1.732050807568877*alpha[3]*fr[12]-1.0*alpha[2]*fr[10]-1.0*alpha[3]*fr[9])*dfac_x; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]-1.732050807568877*alpha[0]*fr[14]+3.0*alpha[2]*fr[13]+3.0*alpha[3]*fr[12]-1.732050807568877*alpha[2]*fr[10]-1.732050807568877*alpha[3]*fr[9])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
  outl[8] += incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += -1.0*incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_Y_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.5*(2.0*BdriftY[0]*m_*wv2+((6.0*Apar[3]-3.464101615137754*Apar[1])*dfac_x-1.732050807568877*BdriftY[1]*Apar[3]-1.732050807568877*BdriftY[0]*Apar[2]+Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*q_*wv+1.732050807568877*BmagInv[0]*Bmag[1]*dfac_x*wm+(1.732050807568877*BmagInv[0]*Phi[1]-3.0*BmagInv[0]*Phi[3])*dfac_x*q_))/q_; 

  double alpha[16]; 
  alpha[0] = (2.0*BdriftY[0]*m_*wv2)/q_+6.0*Apar[3]*dfac_x*wv-3.464101615137754*Apar[1]*dfac_x*wv-1.732050807568877*BdriftY[1]*Apar[3]*wv-1.732050807568877*BdriftY[0]*Apar[2]*wv+Apar[1]*BdriftY[1]*wv+Apar[0]*BdriftY[0]*wv+(1.732050807568877*BmagInv[0]*Bmag[1]*dfac_x*wm)/q_-3.0*BmagInv[0]*Phi[3]*dfac_x+1.732050807568877*BmagInv[0]*Phi[1]*dfac_x; 
  alpha[1] = (2.0*BdriftY[1]*m_*wv2)/q_-1.732050807568877*BdriftY[0]*Apar[3]*wv-1.732050807568877*BdriftY[1]*Apar[2]*wv+Apar[0]*BdriftY[1]*wv+BdriftY[0]*Apar[1]*wv+(1.732050807568877*Bmag[1]*BmagInv[1]*dfac_x*wm)/q_-3.0*BmagInv[1]*Phi[3]*dfac_x+1.732050807568877*BmagInv[1]*Phi[1]*dfac_x; 
  alpha[3] = (1.154700538379252*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alpha[4] = (BmagInv[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alpha[6] = (1.154700538379252*BdriftY[1]*m_*wv)/(dfac_v*q_); 
  alpha[8] = (Bmag[1]*BmagInv[1]*dfac_x)/(dfac_m*q_); 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[8]*fl[12]+1.732050807568877*alpha[6]*fl[11]+1.732050807568877*alpha[4]*fl[9]+alpha[8]*fl[8]+1.732050807568877*alpha[3]*fl[7]+alpha[6]*fl[6]+1.732050807568877*alpha[1]*fl[5]+alpha[4]*fl[4]+alpha[3]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_y; 
  incr[1] = 0.125*(1.732050807568877*alpha[4]*fl[12]+1.732050807568877*alpha[3]*fl[11]+1.732050807568877*alpha[8]*fl[9]+alpha[4]*fl[8]+fl[4]*alpha[8]+1.732050807568877*alpha[6]*fl[7]+alpha[3]*fl[6]+fl[3]*alpha[6]+1.732050807568877*alpha[0]*fl[5]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_y; 
  incr[2] = -0.125*(3.0*alpha[8]*fl[12]+3.0*alpha[6]*fl[11]+3.0*alpha[4]*fl[9]+1.732050807568877*alpha[8]*fl[8]+3.0*alpha[3]*fl[7]+1.732050807568877*alpha[6]*fl[6]+3.0*alpha[1]*fl[5]+1.732050807568877*alpha[4]*fl[4]+1.732050807568877*alpha[3]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_y; 
  incr[3] = 0.125*(1.732050807568877*alpha[8]*fl[15]+1.732050807568877*alpha[4]*fl[14]+alpha[8]*fl[13]+1.732050807568877*alpha[1]*fl[11]+alpha[4]*fl[10]+1.732050807568877*alpha[0]*fl[7]+alpha[1]*fl[6]+1.732050807568877*fl[5]*alpha[6]+fl[1]*alpha[6]+alpha[0]*fl[3]+1.732050807568877*fl[2]*alpha[3]+fl[0]*alpha[3])*dfac_y; 
  incr[4] = 0.125*(1.732050807568877*alpha[6]*fl[15]+1.732050807568877*alpha[3]*fl[14]+alpha[6]*fl[13]+1.732050807568877*alpha[1]*fl[12]+alpha[3]*fl[10]+1.732050807568877*alpha[0]*fl[9]+alpha[1]*fl[8]+1.732050807568877*fl[5]*alpha[8]+fl[1]*alpha[8]+alpha[0]*fl[4]+1.732050807568877*fl[2]*alpha[4]+fl[0]*alpha[4])*dfac_y; 
  incr[5] = -0.125*(3.0*alpha[4]*fl[12]+3.0*alpha[3]*fl[11]+3.0*alpha[8]*fl[9]+1.732050807568877*alpha[4]*fl[8]+1.732050807568877*fl[4]*alpha[8]+3.0*alpha[6]*fl[7]+1.732050807568877*alpha[3]*fl[6]+1.732050807568877*fl[3]*alpha[6]+3.0*alpha[0]*fl[5]+3.0*alpha[1]*fl[2]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_y; 
  incr[6] = 0.125*(1.732050807568877*alpha[4]*fl[15]+1.732050807568877*alpha[8]*fl[14]+alpha[4]*fl[13]+1.732050807568877*alpha[0]*fl[11]+alpha[8]*fl[10]+1.732050807568877*alpha[1]*fl[7]+alpha[0]*fl[6]+1.732050807568877*fl[2]*alpha[6]+fl[0]*alpha[6]+1.732050807568877*alpha[3]*fl[5]+alpha[1]*fl[3]+fl[1]*alpha[3])*dfac_y; 
  incr[7] = -0.125*(3.0*alpha[8]*fl[15]+3.0*alpha[4]*fl[14]+1.732050807568877*alpha[8]*fl[13]+3.0*alpha[1]*fl[11]+1.732050807568877*alpha[4]*fl[10]+3.0*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[6]+3.0*fl[5]*alpha[6]+1.732050807568877*fl[1]*alpha[6]+1.732050807568877*alpha[0]*fl[3]+3.0*fl[2]*alpha[3]+1.732050807568877*fl[0]*alpha[3])*dfac_y; 
  incr[8] = 0.125*(1.732050807568877*alpha[3]*fl[15]+1.732050807568877*alpha[6]*fl[14]+alpha[3]*fl[13]+1.732050807568877*alpha[0]*fl[12]+alpha[6]*fl[10]+1.732050807568877*alpha[1]*fl[9]+alpha[0]*fl[8]+1.732050807568877*fl[2]*alpha[8]+fl[0]*alpha[8]+1.732050807568877*alpha[4]*fl[5]+alpha[1]*fl[4]+fl[1]*alpha[4])*dfac_y; 
  incr[9] = -0.125*(3.0*alpha[6]*fl[15]+3.0*alpha[3]*fl[14]+1.732050807568877*alpha[6]*fl[13]+3.0*alpha[1]*fl[12]+1.732050807568877*alpha[3]*fl[10]+3.0*alpha[0]*fl[9]+1.732050807568877*alpha[1]*fl[8]+3.0*fl[5]*alpha[8]+1.732050807568877*fl[1]*alpha[8]+1.732050807568877*alpha[0]*fl[4]+3.0*fl[2]*alpha[4]+1.732050807568877*fl[0]*alpha[4])*dfac_y; 
  incr[10] = 0.125*(1.732050807568877*alpha[1]*fl[15]+1.732050807568877*alpha[0]*fl[14]+alpha[1]*fl[13]+1.732050807568877*alpha[6]*fl[12]+1.732050807568877*alpha[8]*fl[11]+alpha[0]*fl[10]+1.732050807568877*alpha[3]*fl[9]+alpha[6]*fl[8]+fl[6]*alpha[8]+1.732050807568877*alpha[4]*fl[7]+alpha[3]*fl[4]+fl[3]*alpha[4])*dfac_y; 
  incr[11] = -0.125*(3.0*alpha[4]*fl[15]+3.0*alpha[8]*fl[14]+1.732050807568877*alpha[4]*fl[13]+3.0*alpha[0]*fl[11]+1.732050807568877*alpha[8]*fl[10]+3.0*alpha[1]*fl[7]+1.732050807568877*alpha[0]*fl[6]+3.0*fl[2]*alpha[6]+1.732050807568877*fl[0]*alpha[6]+3.0*alpha[3]*fl[5]+1.732050807568877*alpha[1]*fl[3]+1.732050807568877*fl[1]*alpha[3])*dfac_y; 
  incr[12] = -0.125*(3.0*alpha[3]*fl[15]+3.0*alpha[6]*fl[14]+1.732050807568877*alpha[3]*fl[13]+3.0*alpha[0]*fl[12]+1.732050807568877*alpha[6]*fl[10]+3.0*alpha[1]*fl[9]+1.732050807568877*alpha[0]*fl[8]+3.0*fl[2]*alpha[8]+1.732050807568877*fl[0]*alpha[8]+3.0*alpha[4]*fl[5]+1.732050807568877*alpha[1]*fl[4]+1.732050807568877*fl[1]*alpha[4])*dfac_y; 
  incr[13] = 0.125*(1.732050807568877*alpha[0]*fl[15]+1.732050807568877*alpha[1]*fl[14]+alpha[0]*fl[13]+1.732050807568877*alpha[3]*fl[12]+1.732050807568877*alpha[4]*fl[11]+alpha[1]*fl[10]+1.732050807568877*alpha[6]*fl[9]+alpha[3]*fl[8]+1.732050807568877*fl[7]*alpha[8]+fl[3]*alpha[8]+alpha[4]*fl[6]+fl[4]*alpha[6])*dfac_y; 
  incr[14] = -0.125*(3.0*alpha[1]*fl[15]+3.0*alpha[0]*fl[14]+1.732050807568877*alpha[1]*fl[13]+3.0*alpha[6]*fl[12]+3.0*alpha[8]*fl[11]+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[3]*fl[9]+1.732050807568877*alpha[6]*fl[8]+1.732050807568877*fl[6]*alpha[8]+3.0*alpha[4]*fl[7]+1.732050807568877*alpha[3]*fl[4]+1.732050807568877*fl[3]*alpha[4])*dfac_y; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+3.0*alpha[1]*fl[14]+1.732050807568877*alpha[0]*fl[13]+3.0*alpha[3]*fl[12]+3.0*alpha[4]*fl[11]+1.732050807568877*alpha[1]*fl[10]+3.0*alpha[6]*fl[9]+1.732050807568877*alpha[3]*fl[8]+3.0*fl[7]*alpha[8]+1.732050807568877*fl[3]*alpha[8]+1.732050807568877*alpha[4]*fl[6]+1.732050807568877*fl[4]*alpha[6])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[8]*fr[12]+1.732050807568877*alpha[6]*fr[11]+1.732050807568877*alpha[4]*fr[9]-1.0*alpha[8]*fr[8]+1.732050807568877*alpha[3]*fr[7]-1.0*alpha[6]*fr[6]+1.732050807568877*alpha[1]*fr[5]-1.0*alpha[4]*fr[4]-1.0*alpha[3]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_y; 
  incr[1] = -0.125*(1.732050807568877*alpha[4]*fr[12]+1.732050807568877*alpha[3]*fr[11]+1.732050807568877*alpha[8]*fr[9]-1.0*alpha[4]*fr[8]-1.0*fr[4]*alpha[8]+1.732050807568877*alpha[6]*fr[7]-1.0*alpha[3]*fr[6]-1.0*fr[3]*alpha[6]+1.732050807568877*alpha[0]*fr[5]+1.732050807568877*alpha[1]*fr[2]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_y; 
  incr[2] = 0.125*(3.0*alpha[8]*fr[12]+3.0*alpha[6]*fr[11]+3.0*alpha[4]*fr[9]-1.732050807568877*alpha[8]*fr[8]+3.0*alpha[3]*fr[7]-1.732050807568877*alpha[6]*fr[6]+3.0*alpha[1]*fr[5]-1.732050807568877*alpha[4]*fr[4]-1.732050807568877*alpha[3]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_y; 
  incr[3] = -0.125*(1.732050807568877*alpha[8]*fr[15]+1.732050807568877*alpha[4]*fr[14]-1.0*alpha[8]*fr[13]+1.732050807568877*alpha[1]*fr[11]-1.0*alpha[4]*fr[10]+1.732050807568877*alpha[0]*fr[7]-1.0*alpha[1]*fr[6]+1.732050807568877*fr[5]*alpha[6]-1.0*fr[1]*alpha[6]-1.0*alpha[0]*fr[3]+1.732050807568877*fr[2]*alpha[3]-1.0*fr[0]*alpha[3])*dfac_y; 
  incr[4] = -0.125*(1.732050807568877*alpha[6]*fr[15]+1.732050807568877*alpha[3]*fr[14]-1.0*alpha[6]*fr[13]+1.732050807568877*alpha[1]*fr[12]-1.0*alpha[3]*fr[10]+1.732050807568877*alpha[0]*fr[9]-1.0*alpha[1]*fr[8]+1.732050807568877*fr[5]*alpha[8]-1.0*fr[1]*alpha[8]-1.0*alpha[0]*fr[4]+1.732050807568877*fr[2]*alpha[4]-1.0*fr[0]*alpha[4])*dfac_y; 
  incr[5] = 0.125*(3.0*alpha[4]*fr[12]+3.0*alpha[3]*fr[11]+3.0*alpha[8]*fr[9]-1.732050807568877*alpha[4]*fr[8]-1.732050807568877*fr[4]*alpha[8]+3.0*alpha[6]*fr[7]-1.732050807568877*alpha[3]*fr[6]-1.732050807568877*fr[3]*alpha[6]+3.0*alpha[0]*fr[5]+3.0*alpha[1]*fr[2]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_y; 
  incr[6] = -0.125*(1.732050807568877*alpha[4]*fr[15]+1.732050807568877*alpha[8]*fr[14]-1.0*alpha[4]*fr[13]+1.732050807568877*alpha[0]*fr[11]-1.0*alpha[8]*fr[10]+1.732050807568877*alpha[1]*fr[7]-1.0*alpha[0]*fr[6]+1.732050807568877*fr[2]*alpha[6]-1.0*fr[0]*alpha[6]+1.732050807568877*alpha[3]*fr[5]-1.0*alpha[1]*fr[3]-1.0*fr[1]*alpha[3])*dfac_y; 
  incr[7] = 0.125*(3.0*alpha[8]*fr[15]+3.0*alpha[4]*fr[14]-1.732050807568877*alpha[8]*fr[13]+3.0*alpha[1]*fr[11]-1.732050807568877*alpha[4]*fr[10]+3.0*alpha[0]*fr[7]-1.732050807568877*alpha[1]*fr[6]+3.0*fr[5]*alpha[6]-1.732050807568877*fr[1]*alpha[6]-1.732050807568877*alpha[0]*fr[3]+3.0*fr[2]*alpha[3]-1.732050807568877*fr[0]*alpha[3])*dfac_y; 
  incr[8] = -0.125*(1.732050807568877*alpha[3]*fr[15]+1.732050807568877*alpha[6]*fr[14]-1.0*alpha[3]*fr[13]+1.732050807568877*alpha[0]*fr[12]-1.0*alpha[6]*fr[10]+1.732050807568877*alpha[1]*fr[9]-1.0*alpha[0]*fr[8]+1.732050807568877*fr[2]*alpha[8]-1.0*fr[0]*alpha[8]+1.732050807568877*alpha[4]*fr[5]-1.0*alpha[1]*fr[4]-1.0*fr[1]*alpha[4])*dfac_y; 
  incr[9] = 0.125*(3.0*alpha[6]*fr[15]+3.0*alpha[3]*fr[14]-1.732050807568877*alpha[6]*fr[13]+3.0*alpha[1]*fr[12]-1.732050807568877*alpha[3]*fr[10]+3.0*alpha[0]*fr[9]-1.732050807568877*alpha[1]*fr[8]+3.0*fr[5]*alpha[8]-1.732050807568877*fr[1]*alpha[8]-1.732050807568877*alpha[0]*fr[4]+3.0*fr[2]*alpha[4]-1.732050807568877*fr[0]*alpha[4])*dfac_y; 
  incr[10] = -0.125*(1.732050807568877*alpha[1]*fr[15]+1.732050807568877*alpha[0]*fr[14]-1.0*alpha[1]*fr[13]+1.732050807568877*alpha[6]*fr[12]+1.732050807568877*alpha[8]*fr[11]-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[3]*fr[9]-1.0*alpha[6]*fr[8]-1.0*fr[6]*alpha[8]+1.732050807568877*alpha[4]*fr[7]-1.0*alpha[3]*fr[4]-1.0*fr[3]*alpha[4])*dfac_y; 
  incr[11] = 0.125*(3.0*alpha[4]*fr[15]+3.0*alpha[8]*fr[14]-1.732050807568877*alpha[4]*fr[13]+3.0*alpha[0]*fr[11]-1.732050807568877*alpha[8]*fr[10]+3.0*alpha[1]*fr[7]-1.732050807568877*alpha[0]*fr[6]+3.0*fr[2]*alpha[6]-1.732050807568877*fr[0]*alpha[6]+3.0*alpha[3]*fr[5]-1.732050807568877*alpha[1]*fr[3]-1.732050807568877*fr[1]*alpha[3])*dfac_y; 
  incr[12] = 0.125*(3.0*alpha[3]*fr[15]+3.0*alpha[6]*fr[14]-1.732050807568877*alpha[3]*fr[13]+3.0*alpha[0]*fr[12]-1.732050807568877*alpha[6]*fr[10]+3.0*alpha[1]*fr[9]-1.732050807568877*alpha[0]*fr[8]+3.0*fr[2]*alpha[8]-1.732050807568877*fr[0]*alpha[8]+3.0*alpha[4]*fr[5]-1.732050807568877*alpha[1]*fr[4]-1.732050807568877*fr[1]*alpha[4])*dfac_y; 
  incr[13] = -0.125*(1.732050807568877*alpha[0]*fr[15]+1.732050807568877*alpha[1]*fr[14]-1.0*alpha[0]*fr[13]+1.732050807568877*alpha[3]*fr[12]+1.732050807568877*alpha[4]*fr[11]-1.0*alpha[1]*fr[10]+1.732050807568877*alpha[6]*fr[9]-1.0*alpha[3]*fr[8]+1.732050807568877*fr[7]*alpha[8]-1.0*fr[3]*alpha[8]-1.0*alpha[4]*fr[6]-1.0*fr[4]*alpha[6])*dfac_y; 
  incr[14] = 0.125*(3.0*alpha[1]*fr[15]+3.0*alpha[0]*fr[14]-1.732050807568877*alpha[1]*fr[13]+3.0*alpha[6]*fr[12]+3.0*alpha[8]*fr[11]-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[3]*fr[9]-1.732050807568877*alpha[6]*fr[8]-1.732050807568877*fr[6]*alpha[8]+3.0*alpha[4]*fr[7]-1.732050807568877*alpha[3]*fr[4]-1.732050807568877*fr[3]*alpha[4])*dfac_y; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]+3.0*alpha[1]*fr[14]-1.732050807568877*alpha[0]*fr[13]+3.0*alpha[3]*fr[12]+3.0*alpha[4]*fr[11]-1.732050807568877*alpha[1]*fr[10]+3.0*alpha[6]*fr[9]-1.732050807568877*alpha[3]*fr[8]+3.0*fr[7]*alpha[8]-1.732050807568877*fr[3]*alpha[8]-1.732050807568877*alpha[4]*fr[6]-1.732050807568877*fr[4]*alpha[6])*dfac_y; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += -1.0*incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += incr[9]; 
  outl[10] += -1.0*incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += incr[12]; 
  outl[13] += -1.0*incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf2x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double dfac_v = 2.0/dxv[2]; 
  double dfac_m = 2.0/dxv[3]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double wv = w[2]; 
  double wm = w[3]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*((3.464101615137754*BdriftX[0]*Bmag[1]*dfac_v*dfac_x*m_*wm+((3.464101615137754*BdriftY[1]*Phi[3]+3.464101615137754*BdriftY[0]*Phi[2])*dfac_v*dfac_y+3.464101615137754*BdriftX[0]*Phi[1]*dfac_v*dfac_x)*m_*q_)*wv+((6.0*Bmag[1]*Apar[2]*dfac_v*dfac_x*dfac_y+(1.732050807568877*Apar[1]*BdriftX[1]+1.732050807568877*Apar[0]*BdriftX[0])*Bmag[1]*dfac_v*dfac_x)*q_-3.464101615137754*BdriftX[0]*Bmag[1]*dfac_x*m_)*wm+(((-3.464101615137754*BdriftY[1]*Phi[3])-3.464101615137754*BdriftY[0]*Phi[2])*dfac_y-3.464101615137754*BdriftX[0]*Phi[1]*dfac_x)*m_*q_+(((6.0*Phi[1]*Apar[2]-6.0*Apar[1]*Phi[2])*dfac_v*dfac_x+((1.732050807568877*Apar[0]*BdriftY[1]+1.732050807568877*BdriftY[0]*Apar[1])*Phi[3]+(1.732050807568877*Apar[1]*BdriftY[1]+1.732050807568877*Apar[0]*BdriftY[0])*Phi[2])*dfac_v)*dfac_y+((1.732050807568877*BdriftX[1]*Apar[3]+1.732050807568877*BdriftX[0]*Apar[2])*Phi[3]+(1.732050807568877*Apar[1]*BdriftX[1]+1.732050807568877*Apar[0]*BdriftX[0])*Phi[1])*dfac_v*dfac_x+4.0*dApardt[0]*dfac_v)*q2))/(dfac_v*m_*q_); 

  double alpha[16]; 
  alpha[0] = (-(1.732050807568877*BdriftX[0]*Bmag[1]*dfac_x*wm*wv)/q_)-1.732050807568877*BdriftY[1]*Phi[3]*dfac_y*wv-1.732050807568877*BdriftY[0]*Phi[2]*dfac_y*wv-1.732050807568877*BdriftX[0]*Phi[1]*dfac_x*wv+(1.732050807568877*BdriftX[0]*Bmag[1]*dfac_x*wm)/(dfac_v*q_)-(3.0*Bmag[1]*Apar[2]*dfac_x*dfac_y*wm)/m_-(0.8660254037844386*Apar[1]*BdriftX[1]*Bmag[1]*dfac_x*wm)/m_-(0.8660254037844386*Apar[0]*BdriftX[0]*Bmag[1]*dfac_x*wm)/m_+(3.0*Apar[1]*Phi[2]*dfac_x*dfac_y*q_)/m_-(3.0*Phi[1]*Apar[2]*dfac_x*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[1]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[1]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*Apar[1]*BdriftY[1]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[0]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*BdriftX[1]*Apar[3]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[2]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*Apar[1]*BdriftX[1]*Phi[1]*dfac_x*q_)/m_-(0.8660254037844386*Apar[0]*BdriftX[0]*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[0]*q_)/m_+(1.732050807568877*BdriftY[1]*Phi[3]*dfac_y)/dfac_v+(1.732050807568877*BdriftY[0]*Phi[2]*dfac_y)/dfac_v+(1.732050807568877*BdriftX[0]*Phi[1]*dfac_x)/dfac_v; 
  alpha[1] = (-(1.732050807568877*BdriftX[1]*Bmag[1]*dfac_x*wm*wv)/q_)-1.732050807568877*BdriftY[0]*Phi[3]*dfac_y*wv-1.732050807568877*BdriftY[1]*Phi[2]*dfac_y*wv-1.732050807568877*BdriftX[1]*Phi[1]*dfac_x*wv+(1.732050807568877*BdriftX[1]*Bmag[1]*dfac_x*wm)/(dfac_v*q_)-(3.0*Bmag[1]*Apar[3]*dfac_x*dfac_y*wm)/m_-(0.8660254037844386*Apar[0]*BdriftX[1]*Bmag[1]*dfac_x*wm)/m_-(0.8660254037844386*BdriftX[0]*Apar[1]*Bmag[1]*dfac_x*wm)/m_+(3.0*Apar[1]*Phi[3]*dfac_x*dfac_y*q_)/m_-(3.0*Phi[1]*Apar[3]*dfac_x*dfac_y*q_)/m_-(1.55884572681199*Apar[1]*BdriftY[1]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[0]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[1]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[1]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[3]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[1]*Apar[2]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*Apar[0]*BdriftX[1]*Phi[1]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[1]*Phi[1]*dfac_x*q_)/m_-(2.0*dApardt[1]*q_)/m_+(1.732050807568877*BdriftY[0]*Phi[3]*dfac_y)/dfac_v+(1.732050807568877*BdriftY[1]*Phi[2]*dfac_y)/dfac_v+(1.732050807568877*BdriftX[1]*Phi[1]*dfac_x)/dfac_v; 
  alpha[2] = (-1.732050807568877*BdriftX[0]*Phi[3]*dfac_x*wv)-(0.8660254037844386*BdriftX[1]*Bmag[1]*Apar[3]*dfac_x*wm)/m_-(0.8660254037844386*BdriftX[0]*Bmag[1]*Apar[2]*dfac_x*wm)/m_-(3.0*Apar[2]*Phi[3]*dfac_x*dfac_y*q_)/m_+(3.0*Phi[2]*Apar[3]*dfac_x*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[3]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[1]*Apar[2]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[1]*Phi[2]*Apar[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[2]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*Apar[1]*BdriftX[1]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*Apar[0]*BdriftX[0]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[1]*Phi[1]*Apar[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Phi[1]*Apar[2]*dfac_x*q_)/m_-(2.0*dApardt[2]*q_)/m_+(1.732050807568877*BdriftX[0]*Phi[3]*dfac_x)/dfac_v; 
  alpha[4] = (-(1.0*BdriftX[0]*Bmag[1]*dfac_x*wv)/(dfac_m*q_))+(1.0*BdriftX[0]*Bmag[1]*dfac_x)/(dfac_m*dfac_v*q_)-(1.732050807568877*Bmag[1]*Apar[2]*dfac_x*dfac_y)/(dfac_m*m_)-(0.5*Apar[1]*BdriftX[1]*Bmag[1]*dfac_x)/(dfac_m*m_)-(0.5*Apar[0]*BdriftX[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[5] = (-1.732050807568877*BdriftX[1]*Phi[3]*dfac_x*wv)-(0.8660254037844386*BdriftX[0]*Bmag[1]*Apar[3]*dfac_x*wm)/m_-(0.8660254037844386*BdriftX[1]*Bmag[1]*Apar[2]*dfac_x*wm)/m_-(1.55884572681199*BdriftY[1]*Apar[3]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[2]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Phi[2]*Apar[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[1]*Apar[2]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftX[1]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[1]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Phi[1]*Apar[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[1]*Phi[1]*Apar[2]*dfac_x*q_)/m_-(2.0*dApardt[3]*q_)/m_+(1.732050807568877*BdriftX[1]*Phi[3]*dfac_x)/dfac_v; 
  alpha[8] = (-(1.0*BdriftX[1]*Bmag[1]*dfac_x*wv)/(dfac_m*q_))+(1.0*BdriftX[1]*Bmag[1]*dfac_x)/(dfac_m*dfac_v*q_)-(1.732050807568877*Bmag[1]*Apar[3]*dfac_x*dfac_y)/(dfac_m*m_)-(0.5*Apar[0]*BdriftX[1]*Bmag[1]*dfac_x)/(dfac_m*m_)-(0.5*BdriftX[0]*Apar[1]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alpha[9] = (-(0.5*BdriftX[1]*Bmag[1]*Apar[3]*dfac_x)/(dfac_m*m_))-(0.5*BdriftX[0]*Bmag[1]*Apar[2]*dfac_x)/(dfac_m*m_); 
  alpha[12] = (-(0.5*BdriftX[0]*Bmag[1]*Apar[3]*dfac_x)/(dfac_m*m_))-(0.5*BdriftX[1]*Bmag[1]*Apar[2]*dfac_x)/(dfac_m*m_); 
  if (alpha0>0) { 
  incr[0] = 0.125*(1.732050807568877*alpha[12]*fl[15]+1.732050807568877*alpha[9]*fl[14]+1.732050807568877*alpha[8]*fl[13]+alpha[12]*fl[12]+1.732050807568877*alpha[5]*fl[11]+1.732050807568877*alpha[4]*fl[10]+alpha[9]*fl[9]+alpha[8]*fl[8]+1.732050807568877*alpha[2]*fl[7]+1.732050807568877*alpha[1]*fl[6]+alpha[5]*fl[5]+alpha[4]*fl[4]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.125*(1.732050807568877*alpha[9]*fl[15]+1.732050807568877*alpha[12]*fl[14]+1.732050807568877*alpha[4]*fl[13]+alpha[9]*fl[12]+fl[9]*alpha[12]+1.732050807568877*alpha[2]*fl[11]+1.732050807568877*alpha[8]*fl[10]+alpha[4]*fl[8]+fl[4]*alpha[8]+1.732050807568877*alpha[5]*fl[7]+1.732050807568877*alpha[0]*fl[6]+alpha[2]*fl[5]+fl[2]*alpha[5]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.125*(1.732050807568877*alpha[8]*fl[15]+1.732050807568877*alpha[4]*fl[14]+1.732050807568877*alpha[12]*fl[13]+alpha[8]*fl[12]+fl[8]*alpha[12]+1.732050807568877*alpha[1]*fl[11]+1.732050807568877*alpha[9]*fl[10]+alpha[4]*fl[9]+fl[4]*alpha[9]+1.732050807568877*alpha[0]*fl[7]+1.732050807568877*alpha[5]*fl[6]+alpha[1]*fl[5]+fl[1]*alpha[5]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.125*(3.0*alpha[12]*fl[15]+3.0*alpha[9]*fl[14]+3.0*alpha[8]*fl[13]+1.732050807568877*alpha[12]*fl[12]+3.0*alpha[5]*fl[11]+3.0*alpha[4]*fl[10]+1.732050807568877*alpha[9]*fl[9]+1.732050807568877*alpha[8]*fl[8]+3.0*alpha[2]*fl[7]+3.0*alpha[1]*fl[6]+1.732050807568877*alpha[5]*fl[5]+1.732050807568877*alpha[4]*fl[4]+3.0*alpha[0]*fl[3]+1.732050807568877*alpha[2]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[4] = 0.125*(1.732050807568877*alpha[5]*fl[15]+1.732050807568877*alpha[2]*fl[14]+1.732050807568877*alpha[1]*fl[13]+alpha[5]*fl[12]+1.732050807568877*fl[11]*alpha[12]+fl[5]*alpha[12]+1.732050807568877*alpha[0]*fl[10]+alpha[2]*fl[9]+1.732050807568877*fl[7]*alpha[9]+fl[2]*alpha[9]+alpha[1]*fl[8]+1.732050807568877*fl[6]*alpha[8]+fl[1]*alpha[8]+alpha[0]*fl[4]+1.732050807568877*fl[3]*alpha[4]+fl[0]*alpha[4])*dfac_v; 
  incr[5] = 0.125*(1.732050807568877*alpha[4]*fl[15]+1.732050807568877*alpha[8]*fl[14]+1.732050807568877*alpha[9]*fl[13]+alpha[4]*fl[12]+1.732050807568877*fl[10]*alpha[12]+fl[4]*alpha[12]+1.732050807568877*alpha[0]*fl[11]+alpha[8]*fl[9]+fl[8]*alpha[9]+1.732050807568877*alpha[1]*fl[7]+1.732050807568877*alpha[2]*fl[6]+alpha[0]*fl[5]+1.732050807568877*fl[3]*alpha[5]+fl[0]*alpha[5]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.125*(3.0*alpha[9]*fl[15]+3.0*alpha[12]*fl[14]+3.0*alpha[4]*fl[13]+1.732050807568877*alpha[9]*fl[12]+1.732050807568877*fl[9]*alpha[12]+3.0*alpha[2]*fl[11]+3.0*alpha[8]*fl[10]+1.732050807568877*alpha[4]*fl[8]+1.732050807568877*fl[4]*alpha[8]+3.0*alpha[5]*fl[7]+3.0*alpha[0]*fl[6]+1.732050807568877*alpha[2]*fl[5]+1.732050807568877*fl[2]*alpha[5]+3.0*alpha[1]*fl[3]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 
  incr[7] = -0.125*(3.0*alpha[8]*fl[15]+3.0*alpha[4]*fl[14]+3.0*alpha[12]*fl[13]+1.732050807568877*alpha[8]*fl[12]+1.732050807568877*fl[8]*alpha[12]+3.0*alpha[1]*fl[11]+3.0*alpha[9]*fl[10]+1.732050807568877*alpha[4]*fl[9]+1.732050807568877*fl[4]*alpha[9]+3.0*alpha[0]*fl[7]+3.0*alpha[5]*fl[6]+1.732050807568877*alpha[1]*fl[5]+1.732050807568877*fl[1]*alpha[5]+3.0*alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+1.732050807568877*fl[0]*alpha[2])*dfac_v; 
  incr[8] = 0.125*(1.732050807568877*alpha[2]*fl[15]+1.732050807568877*alpha[5]*fl[14]+1.732050807568877*alpha[0]*fl[13]+alpha[2]*fl[12]+1.732050807568877*fl[7]*alpha[12]+fl[2]*alpha[12]+1.732050807568877*alpha[9]*fl[11]+1.732050807568877*alpha[1]*fl[10]+alpha[5]*fl[9]+fl[5]*alpha[9]+alpha[0]*fl[8]+1.732050807568877*fl[3]*alpha[8]+fl[0]*alpha[8]+1.732050807568877*alpha[4]*fl[6]+alpha[1]*fl[4]+fl[1]*alpha[4])*dfac_v; 
  incr[9] = 0.125*(1.732050807568877*alpha[1]*fl[15]+1.732050807568877*alpha[0]*fl[14]+1.732050807568877*alpha[5]*fl[13]+alpha[1]*fl[12]+1.732050807568877*fl[6]*alpha[12]+fl[1]*alpha[12]+1.732050807568877*alpha[8]*fl[11]+1.732050807568877*alpha[2]*fl[10]+alpha[0]*fl[9]+1.732050807568877*fl[3]*alpha[9]+fl[0]*alpha[9]+alpha[5]*fl[8]+fl[5]*alpha[8]+1.732050807568877*alpha[4]*fl[7]+alpha[2]*fl[4]+fl[2]*alpha[4])*dfac_v; 
  incr[10] = -0.125*(3.0*alpha[5]*fl[15]+3.0*alpha[2]*fl[14]+3.0*alpha[1]*fl[13]+1.732050807568877*alpha[5]*fl[12]+3.0*fl[11]*alpha[12]+1.732050807568877*fl[5]*alpha[12]+3.0*alpha[0]*fl[10]+1.732050807568877*alpha[2]*fl[9]+3.0*fl[7]*alpha[9]+1.732050807568877*fl[2]*alpha[9]+1.732050807568877*alpha[1]*fl[8]+3.0*fl[6]*alpha[8]+1.732050807568877*fl[1]*alpha[8]+1.732050807568877*alpha[0]*fl[4]+3.0*fl[3]*alpha[4]+1.732050807568877*fl[0]*alpha[4])*dfac_v; 
  incr[11] = -0.125*(3.0*alpha[4]*fl[15]+3.0*alpha[8]*fl[14]+3.0*alpha[9]*fl[13]+1.732050807568877*alpha[4]*fl[12]+3.0*fl[10]*alpha[12]+1.732050807568877*fl[4]*alpha[12]+3.0*alpha[0]*fl[11]+1.732050807568877*alpha[8]*fl[9]+1.732050807568877*fl[8]*alpha[9]+3.0*alpha[1]*fl[7]+3.0*alpha[2]*fl[6]+1.732050807568877*alpha[0]*fl[5]+3.0*fl[3]*alpha[5]+1.732050807568877*fl[0]*alpha[5]+1.732050807568877*alpha[1]*fl[2]+1.732050807568877*fl[1]*alpha[2])*dfac_v; 
  incr[12] = 0.125*(1.732050807568877*alpha[0]*fl[15]+1.732050807568877*alpha[1]*fl[14]+1.732050807568877*alpha[2]*fl[13]+alpha[0]*fl[12]+1.732050807568877*fl[3]*alpha[12]+fl[0]*alpha[12]+1.732050807568877*alpha[4]*fl[11]+1.732050807568877*alpha[5]*fl[10]+alpha[1]*fl[9]+1.732050807568877*fl[6]*alpha[9]+fl[1]*alpha[9]+alpha[2]*fl[8]+1.732050807568877*fl[7]*alpha[8]+fl[2]*alpha[8]+alpha[4]*fl[5]+fl[4]*alpha[5])*dfac_v; 
  incr[13] = -0.125*(3.0*alpha[2]*fl[15]+3.0*alpha[5]*fl[14]+3.0*alpha[0]*fl[13]+1.732050807568877*alpha[2]*fl[12]+3.0*fl[7]*alpha[12]+1.732050807568877*fl[2]*alpha[12]+3.0*alpha[9]*fl[11]+3.0*alpha[1]*fl[10]+1.732050807568877*alpha[5]*fl[9]+1.732050807568877*fl[5]*alpha[9]+1.732050807568877*alpha[0]*fl[8]+3.0*fl[3]*alpha[8]+1.732050807568877*fl[0]*alpha[8]+3.0*alpha[4]*fl[6]+1.732050807568877*alpha[1]*fl[4]+1.732050807568877*fl[1]*alpha[4])*dfac_v; 
  incr[14] = -0.125*(3.0*alpha[1]*fl[15]+3.0*alpha[0]*fl[14]+3.0*alpha[5]*fl[13]+1.732050807568877*alpha[1]*fl[12]+3.0*fl[6]*alpha[12]+1.732050807568877*fl[1]*alpha[12]+3.0*alpha[8]*fl[11]+3.0*alpha[2]*fl[10]+1.732050807568877*alpha[0]*fl[9]+3.0*fl[3]*alpha[9]+1.732050807568877*fl[0]*alpha[9]+1.732050807568877*alpha[5]*fl[8]+1.732050807568877*fl[5]*alpha[8]+3.0*alpha[4]*fl[7]+1.732050807568877*alpha[2]*fl[4]+1.732050807568877*fl[2]*alpha[4])*dfac_v; 
  incr[15] = -0.125*(3.0*alpha[0]*fl[15]+3.0*alpha[1]*fl[14]+3.0*alpha[2]*fl[13]+1.732050807568877*alpha[0]*fl[12]+3.0*fl[3]*alpha[12]+1.732050807568877*fl[0]*alpha[12]+3.0*alpha[4]*fl[11]+3.0*alpha[5]*fl[10]+1.732050807568877*alpha[1]*fl[9]+3.0*fl[6]*alpha[9]+1.732050807568877*fl[1]*alpha[9]+1.732050807568877*alpha[2]*fl[8]+3.0*fl[7]*alpha[8]+1.732050807568877*fl[2]*alpha[8]+1.732050807568877*alpha[4]*fl[5]+1.732050807568877*fl[4]*alpha[5])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } else { 
  incr[0] = -0.125*(1.732050807568877*alpha[12]*fr[15]+1.732050807568877*alpha[9]*fr[14]+1.732050807568877*alpha[8]*fr[13]-1.0*alpha[12]*fr[12]+1.732050807568877*alpha[5]*fr[11]+1.732050807568877*alpha[4]*fr[10]-1.0*alpha[9]*fr[9]-1.0*alpha[8]*fr[8]+1.732050807568877*alpha[2]*fr[7]+1.732050807568877*alpha[1]*fr[6]-1.0*alpha[5]*fr[5]-1.0*alpha[4]*fr[4]+1.732050807568877*alpha[0]*fr[3]-1.0*alpha[2]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.125*(1.732050807568877*alpha[9]*fr[15]+1.732050807568877*alpha[12]*fr[14]+1.732050807568877*alpha[4]*fr[13]-1.0*alpha[9]*fr[12]-1.0*fr[9]*alpha[12]+1.732050807568877*alpha[2]*fr[11]+1.732050807568877*alpha[8]*fr[10]-1.0*alpha[4]*fr[8]-1.0*fr[4]*alpha[8]+1.732050807568877*alpha[5]*fr[7]+1.732050807568877*alpha[0]*fr[6]-1.0*alpha[2]*fr[5]-1.0*fr[2]*alpha[5]+1.732050807568877*alpha[1]*fr[3]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = -0.125*(1.732050807568877*alpha[8]*fr[15]+1.732050807568877*alpha[4]*fr[14]+1.732050807568877*alpha[12]*fr[13]-1.0*alpha[8]*fr[12]-1.0*fr[8]*alpha[12]+1.732050807568877*alpha[1]*fr[11]+1.732050807568877*alpha[9]*fr[10]-1.0*alpha[4]*fr[9]-1.0*fr[4]*alpha[9]+1.732050807568877*alpha[0]*fr[7]+1.732050807568877*alpha[5]*fr[6]-1.0*alpha[1]*fr[5]-1.0*fr[1]*alpha[5]+1.732050807568877*alpha[2]*fr[3]-1.0*alpha[0]*fr[2]-1.0*fr[0]*alpha[2])*dfac_v; 
  incr[3] = 0.125*(3.0*alpha[12]*fr[15]+3.0*alpha[9]*fr[14]+3.0*alpha[8]*fr[13]-1.732050807568877*alpha[12]*fr[12]+3.0*alpha[5]*fr[11]+3.0*alpha[4]*fr[10]-1.732050807568877*alpha[9]*fr[9]-1.732050807568877*alpha[8]*fr[8]+3.0*alpha[2]*fr[7]+3.0*alpha[1]*fr[6]-1.732050807568877*alpha[5]*fr[5]-1.732050807568877*alpha[4]*fr[4]+3.0*alpha[0]*fr[3]-1.732050807568877*alpha[2]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[4] = -0.125*(1.732050807568877*alpha[5]*fr[15]+1.732050807568877*alpha[2]*fr[14]+1.732050807568877*alpha[1]*fr[13]-1.0*alpha[5]*fr[12]+1.732050807568877*fr[11]*alpha[12]-1.0*fr[5]*alpha[12]+1.732050807568877*alpha[0]*fr[10]-1.0*alpha[2]*fr[9]+1.732050807568877*fr[7]*alpha[9]-1.0*fr[2]*alpha[9]-1.0*alpha[1]*fr[8]+1.732050807568877*fr[6]*alpha[8]-1.0*fr[1]*alpha[8]-1.0*alpha[0]*fr[4]+1.732050807568877*fr[3]*alpha[4]-1.0*fr[0]*alpha[4])*dfac_v; 
  incr[5] = -0.125*(1.732050807568877*alpha[4]*fr[15]+1.732050807568877*alpha[8]*fr[14]+1.732050807568877*alpha[9]*fr[13]-1.0*alpha[4]*fr[12]+1.732050807568877*fr[10]*alpha[12]-1.0*fr[4]*alpha[12]+1.732050807568877*alpha[0]*fr[11]-1.0*alpha[8]*fr[9]-1.0*fr[8]*alpha[9]+1.732050807568877*alpha[1]*fr[7]+1.732050807568877*alpha[2]*fr[6]-1.0*alpha[0]*fr[5]+1.732050807568877*fr[3]*alpha[5]-1.0*fr[0]*alpha[5]-1.0*alpha[1]*fr[2]-1.0*fr[1]*alpha[2])*dfac_v; 
  incr[6] = 0.125*(3.0*alpha[9]*fr[15]+3.0*alpha[12]*fr[14]+3.0*alpha[4]*fr[13]-1.732050807568877*alpha[9]*fr[12]-1.732050807568877*fr[9]*alpha[12]+3.0*alpha[2]*fr[11]+3.0*alpha[8]*fr[10]-1.732050807568877*alpha[4]*fr[8]-1.732050807568877*fr[4]*alpha[8]+3.0*alpha[5]*fr[7]+3.0*alpha[0]*fr[6]-1.732050807568877*alpha[2]*fr[5]-1.732050807568877*fr[2]*alpha[5]+3.0*alpha[1]*fr[3]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 
  incr[7] = 0.125*(3.0*alpha[8]*fr[15]+3.0*alpha[4]*fr[14]+3.0*alpha[12]*fr[13]-1.732050807568877*alpha[8]*fr[12]-1.732050807568877*fr[8]*alpha[12]+3.0*alpha[1]*fr[11]+3.0*alpha[9]*fr[10]-1.732050807568877*alpha[4]*fr[9]-1.732050807568877*fr[4]*alpha[9]+3.0*alpha[0]*fr[7]+3.0*alpha[5]*fr[6]-1.732050807568877*alpha[1]*fr[5]-1.732050807568877*fr[1]*alpha[5]+3.0*alpha[2]*fr[3]-1.732050807568877*alpha[0]*fr[2]-1.732050807568877*fr[0]*alpha[2])*dfac_v; 
  incr[8] = -0.125*(1.732050807568877*alpha[2]*fr[15]+1.732050807568877*alpha[5]*fr[14]+1.732050807568877*alpha[0]*fr[13]-1.0*alpha[2]*fr[12]+1.732050807568877*fr[7]*alpha[12]-1.0*fr[2]*alpha[12]+1.732050807568877*alpha[9]*fr[11]+1.732050807568877*alpha[1]*fr[10]-1.0*alpha[5]*fr[9]-1.0*fr[5]*alpha[9]-1.0*alpha[0]*fr[8]+1.732050807568877*fr[3]*alpha[8]-1.0*fr[0]*alpha[8]+1.732050807568877*alpha[4]*fr[6]-1.0*alpha[1]*fr[4]-1.0*fr[1]*alpha[4])*dfac_v; 
  incr[9] = -0.125*(1.732050807568877*alpha[1]*fr[15]+1.732050807568877*alpha[0]*fr[14]+1.732050807568877*alpha[5]*fr[13]-1.0*alpha[1]*fr[12]+1.732050807568877*fr[6]*alpha[12]-1.0*fr[1]*alpha[12]+1.732050807568877*alpha[8]*fr[11]+1.732050807568877*alpha[2]*fr[10]-1.0*alpha[0]*fr[9]+1.732050807568877*fr[3]*alpha[9]-1.0*fr[0]*alpha[9]-1.0*alpha[5]*fr[8]-1.0*fr[5]*alpha[8]+1.732050807568877*alpha[4]*fr[7]-1.0*alpha[2]*fr[4]-1.0*fr[2]*alpha[4])*dfac_v; 
  incr[10] = 0.125*(3.0*alpha[5]*fr[15]+3.0*alpha[2]*fr[14]+3.0*alpha[1]*fr[13]-1.732050807568877*alpha[5]*fr[12]+3.0*fr[11]*alpha[12]-1.732050807568877*fr[5]*alpha[12]+3.0*alpha[0]*fr[10]-1.732050807568877*alpha[2]*fr[9]+3.0*fr[7]*alpha[9]-1.732050807568877*fr[2]*alpha[9]-1.732050807568877*alpha[1]*fr[8]+3.0*fr[6]*alpha[8]-1.732050807568877*fr[1]*alpha[8]-1.732050807568877*alpha[0]*fr[4]+3.0*fr[3]*alpha[4]-1.732050807568877*fr[0]*alpha[4])*dfac_v; 
  incr[11] = 0.125*(3.0*alpha[4]*fr[15]+3.0*alpha[8]*fr[14]+3.0*alpha[9]*fr[13]-1.732050807568877*alpha[4]*fr[12]+3.0*fr[10]*alpha[12]-1.732050807568877*fr[4]*alpha[12]+3.0*alpha[0]*fr[11]-1.732050807568877*alpha[8]*fr[9]-1.732050807568877*fr[8]*alpha[9]+3.0*alpha[1]*fr[7]+3.0*alpha[2]*fr[6]-1.732050807568877*alpha[0]*fr[5]+3.0*fr[3]*alpha[5]-1.732050807568877*fr[0]*alpha[5]-1.732050807568877*alpha[1]*fr[2]-1.732050807568877*fr[1]*alpha[2])*dfac_v; 
  incr[12] = -0.125*(1.732050807568877*alpha[0]*fr[15]+1.732050807568877*alpha[1]*fr[14]+1.732050807568877*alpha[2]*fr[13]-1.0*alpha[0]*fr[12]+1.732050807568877*fr[3]*alpha[12]-1.0*fr[0]*alpha[12]+1.732050807568877*alpha[4]*fr[11]+1.732050807568877*alpha[5]*fr[10]-1.0*alpha[1]*fr[9]+1.732050807568877*fr[6]*alpha[9]-1.0*fr[1]*alpha[9]-1.0*alpha[2]*fr[8]+1.732050807568877*fr[7]*alpha[8]-1.0*fr[2]*alpha[8]-1.0*alpha[4]*fr[5]-1.0*fr[4]*alpha[5])*dfac_v; 
  incr[13] = 0.125*(3.0*alpha[2]*fr[15]+3.0*alpha[5]*fr[14]+3.0*alpha[0]*fr[13]-1.732050807568877*alpha[2]*fr[12]+3.0*fr[7]*alpha[12]-1.732050807568877*fr[2]*alpha[12]+3.0*alpha[9]*fr[11]+3.0*alpha[1]*fr[10]-1.732050807568877*alpha[5]*fr[9]-1.732050807568877*fr[5]*alpha[9]-1.732050807568877*alpha[0]*fr[8]+3.0*fr[3]*alpha[8]-1.732050807568877*fr[0]*alpha[8]+3.0*alpha[4]*fr[6]-1.732050807568877*alpha[1]*fr[4]-1.732050807568877*fr[1]*alpha[4])*dfac_v; 
  incr[14] = 0.125*(3.0*alpha[1]*fr[15]+3.0*alpha[0]*fr[14]+3.0*alpha[5]*fr[13]-1.732050807568877*alpha[1]*fr[12]+3.0*fr[6]*alpha[12]-1.732050807568877*fr[1]*alpha[12]+3.0*alpha[8]*fr[11]+3.0*alpha[2]*fr[10]-1.732050807568877*alpha[0]*fr[9]+3.0*fr[3]*alpha[9]-1.732050807568877*fr[0]*alpha[9]-1.732050807568877*alpha[5]*fr[8]-1.732050807568877*fr[5]*alpha[8]+3.0*alpha[4]*fr[7]-1.732050807568877*alpha[2]*fr[4]-1.732050807568877*fr[2]*alpha[4])*dfac_v; 
  incr[15] = 0.125*(3.0*alpha[0]*fr[15]+3.0*alpha[1]*fr[14]+3.0*alpha[2]*fr[13]-1.732050807568877*alpha[0]*fr[12]+3.0*fr[3]*alpha[12]-1.732050807568877*fr[0]*alpha[12]+3.0*alpha[4]*fr[11]+3.0*alpha[5]*fr[10]-1.732050807568877*alpha[1]*fr[9]+3.0*fr[6]*alpha[9]-1.732050807568877*fr[1]*alpha[9]-1.732050807568877*alpha[2]*fr[8]+3.0*fr[7]*alpha[8]-1.732050807568877*fr[2]*alpha[8]-1.732050807568877*alpha[4]*fr[5]-1.732050807568877*fr[4]*alpha[5])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 
  outr[4] += incr[4]; 
  outr[5] += incr[5]; 
  outr[6] += incr[6]; 
  outr[7] += incr[7]; 
  outr[8] += incr[8]; 
  outr[9] += incr[9]; 
  outr[10] += incr[10]; 
  outr[11] += incr[11]; 
  outr[12] += incr[12]; 
  outr[13] += incr[13]; 
  outr[14] += incr[14]; 
  outr[15] += incr[15]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += incr[7]; 
  outl[8] += -1.0*incr[8]; 
  outl[9] += -1.0*incr[9]; 
  outl[10] += incr[10]; 
  outl[11] += incr[11]; 
  outl[12] += -1.0*incr[12]; 
  outl[13] += incr[13]; 
  outl[14] += incr[14]; 
  outl[15] += incr[15]; 
  } 
  return std::abs(alpha0); 
} 
