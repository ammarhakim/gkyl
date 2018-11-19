#include <GyrokineticModDecl.h> 
double EmGyrokineticSurf2x2vSer_X_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -0.0625*geoZ[0]*dfac_y*(BmagInv[0]*(3.0*Apar[3]-1.732050807568877*Apar[2])*wv+2.0*(1.732050807568877*Phi[2]-3.0*Phi[3])); 

  double alpha[8]; 
  alpha[0] = -0.25*geoZ[0]*dfac_y*(BmagInv[0]*(4.242640687119286*Apar[3]-2.449489742783178*Apar[2])*wv-8.485281374238571*Phi[3]+4.898979485566357*Phi[2]); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[5]+fl[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[6]+fl[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[8]+fl[4])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[2])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(3.0*fl[6]+1.732050807568877*fl[3])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[11]+fl[7])*dfac_x; 
  incr[8] = -0.1767766952966368*alpha[0]*(3.0*fl[8]+1.732050807568877*fl[4])*dfac_x; 
  incr[9] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[12]+fl[9])*dfac_x; 
  incr[10] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[13]+fl[10])*dfac_x; 
  incr[11] = -0.1767766952966368*alpha[0]*(3.0*fl[11]+1.732050807568877*fl[7])*dfac_x; 
  incr[12] = -0.1767766952966368*alpha[0]*(3.0*fl[12]+1.732050807568877*fl[9])*dfac_x; 
  incr[13] = -0.1767766952966368*alpha[0]*(3.0*fl[13]+1.732050807568877*fl[10])*dfac_x; 
  incr[14] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[15]+fl[14])*dfac_x; 
  incr[15] = -0.1767766952966368*alpha[0]*(3.0*fl[15]+1.732050807568877*fl[14])*dfac_x; 

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
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[6]-1.0*fr[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[8]-1.0*fr[4])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[2])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(3.0*fr[6]-1.732050807568877*fr[3])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[11]-1.0*fr[7])*dfac_x; 
  incr[8] = 0.1767766952966368*alpha[0]*(3.0*fr[8]-1.732050807568877*fr[4])*dfac_x; 
  incr[9] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[12]-1.0*fr[9])*dfac_x; 
  incr[10] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[13]-1.0*fr[10])*dfac_x; 
  incr[11] = 0.1767766952966368*alpha[0]*(3.0*fr[11]-1.732050807568877*fr[7])*dfac_x; 
  incr[12] = 0.1767766952966368*alpha[0]*(3.0*fr[12]-1.732050807568877*fr[9])*dfac_x; 
  incr[13] = 0.1767766952966368*alpha[0]*(3.0*fr[13]-1.732050807568877*fr[10])*dfac_x; 
  incr[14] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[15]-1.0*fr[14])*dfac_x; 
  incr[15] = 0.1767766952966368*alpha[0]*(3.0*fr[15]-1.732050807568877*fr[14])*dfac_x; 

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
double EmGyrokineticSurf2x2vSer_Y_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = 0.0625*geoZ[0]*dfac_x*(BmagInv[0]*(3.0*Apar[3]-1.732050807568877*Apar[1])*wv+2.0*(1.732050807568877*Phi[1]-3.0*Phi[3])); 

  double alpha[8]; 
  alpha[0] = 0.25*geoZ[0]*dfac_x*(BmagInv[0]*(4.242640687119286*Apar[3]-2.449489742783178*Apar[1])*wv-8.485281374238571*Phi[3]+4.898979485566357*Phi[1]); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[2]+fl[0])*dfac_y; 
  incr[1] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[5]+fl[1])*dfac_y; 
  incr[2] = -0.1767766952966368*alpha[0]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_y; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[7]+fl[3])*dfac_y; 
  incr[4] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[9]+fl[4])*dfac_y; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[1])*dfac_y; 
  incr[6] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[11]+fl[6])*dfac_y; 
  incr[7] = -0.1767766952966368*alpha[0]*(3.0*fl[7]+1.732050807568877*fl[3])*dfac_y; 
  incr[8] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[12]+fl[8])*dfac_y; 
  incr[9] = -0.1767766952966368*alpha[0]*(3.0*fl[9]+1.732050807568877*fl[4])*dfac_y; 
  incr[10] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[14]+fl[10])*dfac_y; 
  incr[11] = -0.1767766952966368*alpha[0]*(3.0*fl[11]+1.732050807568877*fl[6])*dfac_y; 
  incr[12] = -0.1767766952966368*alpha[0]*(3.0*fl[12]+1.732050807568877*fl[8])*dfac_y; 
  incr[13] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[15]+fl[13])*dfac_y; 
  incr[14] = -0.1767766952966368*alpha[0]*(3.0*fl[14]+1.732050807568877*fl[10])*dfac_y; 
  incr[15] = -0.1767766952966368*alpha[0]*(3.0*fl[15]+1.732050807568877*fl[13])*dfac_y; 

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
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_y; 
  incr[1] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[1])*dfac_y; 
  incr[2] = 0.1767766952966368*alpha[0]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_y; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[7]-1.0*fr[3])*dfac_y; 
  incr[4] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[9]-1.0*fr[4])*dfac_y; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[1])*dfac_y; 
  incr[6] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[11]-1.0*fr[6])*dfac_y; 
  incr[7] = 0.1767766952966368*alpha[0]*(3.0*fr[7]-1.732050807568877*fr[3])*dfac_y; 
  incr[8] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[12]-1.0*fr[8])*dfac_y; 
  incr[9] = 0.1767766952966368*alpha[0]*(3.0*fr[9]-1.732050807568877*fr[4])*dfac_y; 
  incr[10] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[14]-1.0*fr[10])*dfac_y; 
  incr[11] = 0.1767766952966368*alpha[0]*(3.0*fr[11]-1.732050807568877*fr[6])*dfac_y; 
  incr[12] = 0.1767766952966368*alpha[0]*(3.0*fr[12]-1.732050807568877*fr[8])*dfac_y; 
  incr[13] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[15]-1.0*fr[13])*dfac_y; 
  incr[14] = 0.1767766952966368*alpha[0]*(3.0*fr[14]-1.732050807568877*fr[10])*dfac_y; 
  incr[15] = 0.1767766952966368*alpha[0]*(3.0*fr[15]-1.732050807568877*fr[13])*dfac_y; 

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
double EmGyrokineticSurf2x2vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.03125*(3.0*BmagInv[0]*geoZ[0]*(Apar[1]*Phi[2]-1.0*Phi[1]*Apar[2])*dfac_x*dfac_y-8.0*dApardt[0])*q_)/m_; 

  double alpha[8]; 
  alpha[0] = (0.1767766952966368*(BmagInv[0]*geoZ[0]*(3.0*Apar[1]*Phi[2]-3.0*Phi[1]*Apar[2])*dfac_x*dfac_y-8.0*dApardt[0])*q_)/m_; 
  alpha[1] = (0.1767766952966368*(BmagInv[0]*geoZ[0]*(3.0*Apar[1]*Phi[3]-3.0*Phi[1]*Apar[3])*dfac_x*dfac_y-8.0*dApardt[1])*q_)/m_; 
  alpha[2] = -(0.1767766952966368*(BmagInv[0]*geoZ[0]*(3.0*Apar[2]*Phi[3]-3.0*Phi[2]*Apar[3])*dfac_x*dfac_y+8.0*dApardt[2])*q_)/m_; 
  alpha[4] = -(1.414213562373095*dApardt[3]*q_)/m_; 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[11]+alpha[2]*fl[7]+alpha[1]*fl[6])+alpha[4]*fl[5]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[11]+alpha[4]*fl[7]+alpha[0]*fl[6])+alpha[2]*fl[5]+fl[2]*alpha[4]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[11]+alpha[0]*fl[7]+alpha[4]*fl[6])+alpha[1]*fl[5]+fl[1]*alpha[4]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.1767766952966368*(3.0*(alpha[4]*fl[11]+alpha[2]*fl[7]+alpha[1]*fl[6])+1.732050807568877*alpha[4]*fl[5]+3.0*alpha[0]*fl[3]+1.732050807568877*(alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+alpha[4]*fl[12]+1.732050807568877*alpha[0]*fl[10]+alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4])*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[11]+alpha[1]*fl[7]+alpha[2]*fl[6])+alpha[0]*fl[5]+(1.732050807568877*fl[3]+fl[0])*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*(alpha[2]*fl[11]+alpha[4]*fl[7]+alpha[0]*fl[6])+1.732050807568877*(alpha[2]*fl[5]+fl[2]*alpha[4])+3.0*alpha[1]*fl[3]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[1]*fl[11]+alpha[0]*fl[7]+alpha[4]*fl[6])+1.732050807568877*(alpha[1]*fl[5]+fl[1]*alpha[4])+3.0*alpha[2]*fl[3]+1.732050807568877*(alpha[0]*fl[2]+fl[0]*alpha[2]))*dfac_v; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+alpha[2]*fl[12]+1.732050807568877*alpha[1]*fl[10]+alpha[4]*fl[9]+alpha[0]*fl[8]+alpha[1]*fl[4])*dfac_v; 
  incr[9] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+alpha[1]*fl[12]+1.732050807568877*alpha[2]*fl[10]+alpha[0]*fl[9]+alpha[4]*fl[8]+alpha[2]*fl[4])*dfac_v; 
  incr[10] = -0.1767766952966368*(3.0*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+1.732050807568877*alpha[4]*fl[12]+3.0*alpha[0]*fl[10]+1.732050807568877*(alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[0]*fl[4]))*dfac_v; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[0]*fl[11]+alpha[1]*fl[7]+alpha[2]*fl[6])+1.732050807568877*alpha[0]*fl[5]+3.0*fl[3]*alpha[4]+1.732050807568877*(fl[0]*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2]))*dfac_v; 
  incr[12] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+alpha[0]*fl[12]+1.732050807568877*alpha[4]*fl[10]+alpha[1]*fl[9]+alpha[2]*fl[8]+alpha[4]*fl[4])*dfac_v; 
  incr[13] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+1.732050807568877*alpha[2]*fl[12]+3.0*alpha[1]*fl[10]+1.732050807568877*(alpha[4]*fl[9]+alpha[0]*fl[8]+alpha[1]*fl[4]))*dfac_v; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+1.732050807568877*alpha[1]*fl[12]+3.0*alpha[2]*fl[10]+1.732050807568877*(alpha[0]*fl[9]+alpha[4]*fl[8]+alpha[2]*fl[4]))*dfac_v; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+1.732050807568877*alpha[0]*fl[12]+3.0*alpha[4]*fl[10]+1.732050807568877*(alpha[1]*fl[9]+alpha[2]*fl[8]+alpha[4]*fl[4]))*dfac_v; 

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
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[11]+alpha[2]*fr[7]+alpha[1]*fr[6])-1.0*alpha[4]*fr[5]+1.732050807568877*alpha[0]*fr[3]-1.0*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[11]+alpha[4]*fr[7]+alpha[0]*fr[6])-1.0*(alpha[2]*fr[5]+fr[2]*alpha[4])+1.732050807568877*alpha[1]*fr[3]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[11]+alpha[0]*fr[7]+alpha[4]*fr[6])-1.0*(alpha[1]*fr[5]+fr[1]*alpha[4])+1.732050807568877*alpha[2]*fr[3]-1.0*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[4]*fr[11]+alpha[2]*fr[7]+alpha[1]*fr[6])-1.732050807568877*alpha[4]*fr[5]+3.0*alpha[0]*fr[3]-1.732050807568877*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.0*alpha[4]*fr[12]+1.732050807568877*alpha[0]*fr[10]-1.0*(alpha[2]*fr[9]+alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[11]+alpha[1]*fr[7]+alpha[2]*fr[6])-1.0*alpha[0]*fr[5]+1.732050807568877*fr[3]*alpha[4]-1.0*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[2]*fr[11]+alpha[4]*fr[7]+alpha[0]*fr[6])-1.732050807568877*(alpha[2]*fr[5]+fr[2]*alpha[4])+3.0*alpha[1]*fr[3]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[1]*fr[11]+alpha[0]*fr[7]+alpha[4]*fr[6])-1.732050807568877*(alpha[1]*fr[5]+fr[1]*alpha[4])+3.0*alpha[2]*fr[3]-1.732050807568877*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.0*alpha[2]*fr[12]+1.732050807568877*alpha[1]*fr[10]-1.0*(alpha[4]*fr[9]+alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.0*alpha[1]*fr[12]+1.732050807568877*alpha[2]*fr[10]-1.0*(alpha[0]*fr[9]+alpha[4]*fr[8]+alpha[2]*fr[4]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.732050807568877*alpha[4]*fr[12]+3.0*alpha[0]*fr[10]-1.732050807568877*(alpha[2]*fr[9]+alpha[1]*fr[8]+alpha[0]*fr[4]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[0]*fr[11]+alpha[1]*fr[7]+alpha[2]*fr[6])-1.732050807568877*alpha[0]*fr[5]+3.0*fr[3]*alpha[4]-1.732050807568877*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.0*alpha[0]*fr[12]+1.732050807568877*alpha[4]*fr[10]-1.0*(alpha[1]*fr[9]+alpha[2]*fr[8]+alpha[4]*fr[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.732050807568877*alpha[2]*fr[12]+3.0*alpha[1]*fr[10]-1.732050807568877*(alpha[4]*fr[9]+alpha[0]*fr[8]+alpha[1]*fr[4]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.732050807568877*alpha[1]*fr[12]+3.0*alpha[2]*fr[10]-1.732050807568877*(alpha[0]*fr[9]+alpha[4]*fr[8]+alpha[2]*fr[4]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.732050807568877*alpha[0]*fr[12]+3.0*alpha[4]*fr[10]-1.732050807568877*(alpha[1]*fr[9]+alpha[2]*fr[8]+alpha[4]*fr[4]))*dfac_v; 

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
double EmGyrokineticSurf2x2vSer_X_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -0.0625*dfac_y*((((9.0*BmagInv[1]-5.196152422706631*BmagInv[0])*geoZ[1]+geoZ[0]*(3.0*BmagInv[0]-5.196152422706631*BmagInv[1]))*Apar[3]+((3.0*BmagInv[0]-5.196152422706631*BmagInv[1])*geoZ[1]+geoZ[0]*(3.0*BmagInv[1]-1.732050807568877*BmagInv[0]))*Apar[2])*wv+2.0*((5.196152422706631*geoZ[1]-3.0*geoZ[0])*Phi[3]+(1.732050807568877*geoZ[0]-3.0*geoZ[1])*Phi[2])); 

  double alpha[8]; 
  alpha[0] = -0.25*dfac_y*((((12.72792206135786*BmagInv[1]-7.348469228349534*BmagInv[0])*geoZ[1]+geoZ[0]*(4.242640687119286*BmagInv[0]-7.348469228349534*BmagInv[1]))*Apar[3]+((4.242640687119286*BmagInv[0]-7.348469228349534*BmagInv[1])*geoZ[1]+geoZ[0]*(4.242640687119286*BmagInv[1]-2.449489742783178*BmagInv[0]))*Apar[2])*wv+(14.69693845669907*geoZ[1]-8.485281374238571*geoZ[0])*Phi[3]+(4.898979485566357*geoZ[0]-8.485281374238571*geoZ[1])*Phi[2]); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.1767766952966368*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[5]+fl[2])*dfac_x; 
  incr[3] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[6]+fl[3])*dfac_x; 
  incr[4] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[8]+fl[4])*dfac_x; 
  incr[5] = -0.1767766952966368*alpha[0]*(3.0*fl[5]+1.732050807568877*fl[2])*dfac_x; 
  incr[6] = -0.1767766952966368*alpha[0]*(3.0*fl[6]+1.732050807568877*fl[3])*dfac_x; 
  incr[7] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[11]+fl[7])*dfac_x; 
  incr[8] = -0.1767766952966368*alpha[0]*(3.0*fl[8]+1.732050807568877*fl[4])*dfac_x; 
  incr[9] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[12]+fl[9])*dfac_x; 
  incr[10] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[13]+fl[10])*dfac_x; 
  incr[11] = -0.1767766952966368*alpha[0]*(3.0*fl[11]+1.732050807568877*fl[7])*dfac_x; 
  incr[12] = -0.1767766952966368*alpha[0]*(3.0*fl[12]+1.732050807568877*fl[9])*dfac_x; 
  incr[13] = -0.1767766952966368*alpha[0]*(3.0*fl[13]+1.732050807568877*fl[10])*dfac_x; 
  incr[14] = 0.1767766952966368*alpha[0]*(1.732050807568877*fl[15]+fl[14])*dfac_x; 
  incr[15] = -0.1767766952966368*alpha[0]*(3.0*fl[15]+1.732050807568877*fl[14])*dfac_x; 

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
  incr[0] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.1767766952966368*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[5]-1.0*fr[2])*dfac_x; 
  incr[3] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[6]-1.0*fr[3])*dfac_x; 
  incr[4] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[8]-1.0*fr[4])*dfac_x; 
  incr[5] = 0.1767766952966368*alpha[0]*(3.0*fr[5]-1.732050807568877*fr[2])*dfac_x; 
  incr[6] = 0.1767766952966368*alpha[0]*(3.0*fr[6]-1.732050807568877*fr[3])*dfac_x; 
  incr[7] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[11]-1.0*fr[7])*dfac_x; 
  incr[8] = 0.1767766952966368*alpha[0]*(3.0*fr[8]-1.732050807568877*fr[4])*dfac_x; 
  incr[9] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[12]-1.0*fr[9])*dfac_x; 
  incr[10] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[13]-1.0*fr[10])*dfac_x; 
  incr[11] = 0.1767766952966368*alpha[0]*(3.0*fr[11]-1.732050807568877*fr[7])*dfac_x; 
  incr[12] = 0.1767766952966368*alpha[0]*(3.0*fr[12]-1.732050807568877*fr[9])*dfac_x; 
  incr[13] = 0.1767766952966368*alpha[0]*(3.0*fr[13]-1.732050807568877*fr[10])*dfac_x; 
  incr[14] = -0.1767766952966368*alpha[0]*(1.732050807568877*fr[15]-1.0*fr[14])*dfac_x; 
  incr[15] = 0.1767766952966368*alpha[0]*(3.0*fr[15]-1.732050807568877*fr[14])*dfac_x; 

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
double EmGyrokineticSurf2x2vSer_Y_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = (0.03125*dfac_x*(wv*(3.464101615137754*Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*m_*wv+(3.0*(((2.0*BmagInv[1]-1.0*BmagInv[0]*Bmag[1])*geoZ[1]+geoZ[0]*(2.0*BmagInv[0]-1.0*Bmag[1]*BmagInv[1]))*Apar[3]-1.0*Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2])+1.732050807568877*(((Apar[0]*Bmag[1]-2.0*Apar[1])*BmagInv[1]+BmagInv[0]*Apar[1]*Bmag[1])*geoZ[1]+geoZ[0]*(Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(Apar[0]*Bmag[1]-2.0*Apar[1]))))*q_)+4.0*geoZ[0]*(1.732050807568877*Bmag[1]*wm+(1.732050807568877*Phi[1]-3.0*Phi[3])*q_)))/q_; 

  double alpha[8]; 
  alpha[0] = (0.1767766952966368*dfac_x*(3.464101615137754*Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*m_*wv2+(((6.0*BmagInv[1]-3.0*BmagInv[0]*Bmag[1])*geoZ[1]+geoZ[0]*(6.0*BmagInv[0]-3.0*Bmag[1]*BmagInv[1]))*Apar[3]-3.0*Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2]+((1.732050807568877*Apar[0]*Bmag[1]-3.464101615137754*Apar[1])*BmagInv[1]+1.732050807568877*BmagInv[0]*Apar[1]*Bmag[1])*geoZ[1]+geoZ[0]*(1.732050807568877*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(1.732050807568877*Apar[0]*Bmag[1]-3.464101615137754*Apar[1])))*q_*wv+geoZ[0]*(6.928203230275509*Bmag[1]*wm+(6.928203230275509*Phi[1]-12.0*Phi[3])*q_)))/q_; 
  alpha[1] = (0.03535533905932736*dfac_x*(17.32050807568877*Bmag[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*m_*wv2+(((30.0*BmagInv[0]-27.0*Bmag[1]*BmagInv[1])*geoZ[1]+geoZ[0]*(30.0*BmagInv[1]-15.0*BmagInv[0]*Bmag[1]))*Apar[3]-15.0*Bmag[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2]+(15.58845726811989*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(8.660254037844386*Apar[0]*Bmag[1]-17.32050807568877*Apar[1]))*geoZ[1]+geoZ[0]*((8.660254037844386*Apar[0]*Bmag[1]-17.32050807568877*Apar[1])*BmagInv[1]+8.660254037844386*BmagInv[0]*Apar[1]*Bmag[1]))*q_*wv+geoZ[1]*(34.64101615137754*Bmag[1]*wm+(34.64101615137754*Phi[1]-60.0*Phi[3])*q_)))/q_; 
  alpha[2] = (0.3535533905932737*Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*dfac_x*m_*wv)/(dfac_v*q_); 
  alpha[3] = (0.7071067811865475*geoZ[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alpha[4] = (0.3535533905932737*Bmag[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*dfac_x*m_*wv)/(dfac_v*q_); 
  alpha[5] = (0.7071067811865475*Bmag[1]*geoZ[1]*dfac_x)/(dfac_m*q_); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[12]+alpha[4]*fl[11]+alpha[3]*fl[9])+alpha[5]*fl[8]+1.732050807568877*alpha[2]*fl[7]+alpha[4]*fl[6]+1.732050807568877*alpha[1]*fl[5]+alpha[3]*fl[4]+alpha[2]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_y; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[12]+alpha[2]*fl[11]+alpha[5]*fl[9])+alpha[3]*fl[8]+1.732050807568877*alpha[4]*fl[7]+alpha[2]*fl[6]+1.732050807568877*alpha[0]*fl[5]+fl[4]*alpha[5]+fl[3]*alpha[4]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_y; 
  incr[2] = -0.1767766952966368*(3.0*(alpha[5]*fl[12]+alpha[4]*fl[11]+alpha[3]*fl[9])+1.732050807568877*alpha[5]*fl[8]+3.0*alpha[2]*fl[7]+1.732050807568877*alpha[4]*fl[6]+3.0*alpha[1]*fl[5]+1.732050807568877*(alpha[3]*fl[4]+alpha[2]*fl[3])+3.0*alpha[0]*fl[2]+1.732050807568877*(alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_y; 
  incr[3] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[15]+alpha[3]*fl[14])+alpha[5]*fl[13]+1.732050807568877*alpha[1]*fl[11]+alpha[3]*fl[10]+1.732050807568877*alpha[0]*fl[7]+alpha[1]*fl[6]+alpha[4]*(1.732050807568877*fl[5]+fl[1])+alpha[0]*fl[3]+alpha[2]*(1.732050807568877*fl[2]+fl[0]))*dfac_y; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[15]+alpha[2]*fl[14])+alpha[4]*fl[13]+1.732050807568877*alpha[1]*fl[12]+alpha[2]*fl[10]+1.732050807568877*alpha[0]*fl[9]+alpha[1]*fl[8]+alpha[5]*(1.732050807568877*fl[5]+fl[1])+alpha[0]*fl[4]+(1.732050807568877*fl[2]+fl[0])*alpha[3])*dfac_y; 
  incr[5] = -0.1767766952966368*(3.0*(alpha[3]*fl[12]+alpha[2]*fl[11]+alpha[5]*fl[9])+1.732050807568877*alpha[3]*fl[8]+3.0*alpha[4]*fl[7]+1.732050807568877*alpha[2]*fl[6]+3.0*alpha[0]*fl[5]+1.732050807568877*(fl[4]*alpha[5]+fl[3]*alpha[4])+3.0*alpha[1]*fl[2]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_y; 
  incr[6] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[15]+alpha[5]*fl[14])+alpha[3]*fl[13]+1.732050807568877*alpha[0]*fl[11]+alpha[5]*fl[10]+1.732050807568877*alpha[1]*fl[7]+alpha[0]*fl[6]+1.732050807568877*alpha[2]*fl[5]+(1.732050807568877*fl[2]+fl[0])*alpha[4]+alpha[1]*fl[3]+fl[1]*alpha[2])*dfac_y; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[5]*fl[15]+alpha[3]*fl[14])+1.732050807568877*alpha[5]*fl[13]+3.0*alpha[1]*fl[11]+1.732050807568877*alpha[3]*fl[10]+3.0*alpha[0]*fl[7]+1.732050807568877*alpha[1]*fl[6]+3.0*alpha[4]*fl[5]+1.732050807568877*(fl[1]*alpha[4]+alpha[0]*fl[3])+alpha[2]*(3.0*fl[2]+1.732050807568877*fl[0]))*dfac_y; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[4]*fl[14])+alpha[2]*fl[13]+1.732050807568877*alpha[0]*fl[12]+alpha[4]*fl[10]+1.732050807568877*alpha[1]*fl[9]+alpha[0]*fl[8]+1.732050807568877*alpha[3]*fl[5]+(1.732050807568877*fl[2]+fl[0])*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3])*dfac_y; 
  incr[9] = -0.1767766952966368*(3.0*(alpha[4]*fl[15]+alpha[2]*fl[14])+1.732050807568877*alpha[4]*fl[13]+3.0*alpha[1]*fl[12]+1.732050807568877*alpha[2]*fl[10]+3.0*alpha[0]*fl[9]+1.732050807568877*alpha[1]*fl[8]+3.0*alpha[5]*fl[5]+1.732050807568877*(fl[1]*alpha[5]+alpha[0]*fl[4])+(3.0*fl[2]+1.732050807568877*fl[0])*alpha[3])*dfac_y; 
  incr[10] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14])+alpha[1]*fl[13]+1.732050807568877*(alpha[4]*fl[12]+alpha[5]*fl[11])+alpha[0]*fl[10]+1.732050807568877*alpha[2]*fl[9]+alpha[4]*fl[8]+1.732050807568877*alpha[3]*fl[7]+alpha[5]*fl[6]+alpha[2]*fl[4]+alpha[3]*fl[3])*dfac_y; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[3]*fl[15]+alpha[5]*fl[14])+1.732050807568877*alpha[3]*fl[13]+3.0*alpha[0]*fl[11]+1.732050807568877*alpha[5]*fl[10]+3.0*alpha[1]*fl[7]+1.732050807568877*alpha[0]*fl[6]+3.0*(alpha[2]*fl[5]+fl[2]*alpha[4])+1.732050807568877*(fl[0]*alpha[4]+alpha[1]*fl[3]+fl[1]*alpha[2]))*dfac_y; 
  incr[12] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[4]*fl[14])+1.732050807568877*alpha[2]*fl[13]+3.0*alpha[0]*fl[12]+1.732050807568877*alpha[4]*fl[10]+3.0*alpha[1]*fl[9]+1.732050807568877*alpha[0]*fl[8]+3.0*(alpha[3]*fl[5]+fl[2]*alpha[5])+1.732050807568877*(fl[0]*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3]))*dfac_y; 
  incr[13] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14])+alpha[0]*fl[13]+1.732050807568877*(alpha[2]*fl[12]+alpha[3]*fl[11])+alpha[1]*fl[10]+1.732050807568877*alpha[4]*fl[9]+alpha[2]*fl[8]+1.732050807568877*alpha[5]*fl[7]+alpha[3]*fl[6]+fl[3]*alpha[5]+alpha[4]*fl[4])*dfac_y; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14])+1.732050807568877*alpha[1]*fl[13]+3.0*(alpha[4]*fl[12]+alpha[5]*fl[11])+1.732050807568877*alpha[0]*fl[10]+3.0*alpha[2]*fl[9]+1.732050807568877*alpha[4]*fl[8]+3.0*alpha[3]*fl[7]+1.732050807568877*(alpha[5]*fl[6]+alpha[2]*fl[4]+alpha[3]*fl[3]))*dfac_y; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14])+1.732050807568877*alpha[0]*fl[13]+3.0*(alpha[2]*fl[12]+alpha[3]*fl[11])+1.732050807568877*alpha[1]*fl[10]+3.0*alpha[4]*fl[9]+1.732050807568877*alpha[2]*fl[8]+3.0*alpha[5]*fl[7]+1.732050807568877*(alpha[3]*fl[6]+fl[3]*alpha[5]+alpha[4]*fl[4]))*dfac_y; 

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
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[12]+alpha[4]*fr[11]+alpha[3]*fr[9])-1.0*alpha[5]*fr[8]+1.732050807568877*alpha[2]*fr[7]-1.0*alpha[4]*fr[6]+1.732050807568877*alpha[1]*fr[5]-1.0*(alpha[3]*fr[4]+alpha[2]*fr[3])+1.732050807568877*alpha[0]*fr[2]-1.0*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_y; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[12]+alpha[2]*fr[11]+alpha[5]*fr[9])-1.0*alpha[3]*fr[8]+1.732050807568877*alpha[4]*fr[7]-1.0*alpha[2]*fr[6]+1.732050807568877*alpha[0]*fr[5]-1.0*(fr[4]*alpha[5]+fr[3]*alpha[4])+1.732050807568877*alpha[1]*fr[2]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_y; 
  incr[2] = 0.1767766952966368*(3.0*(alpha[5]*fr[12]+alpha[4]*fr[11]+alpha[3]*fr[9])-1.732050807568877*alpha[5]*fr[8]+3.0*alpha[2]*fr[7]-1.732050807568877*alpha[4]*fr[6]+3.0*alpha[1]*fr[5]-1.732050807568877*(alpha[3]*fr[4]+alpha[2]*fr[3])+3.0*alpha[0]*fr[2]-1.732050807568877*(alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_y; 
  incr[3] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[15]+alpha[3]*fr[14])-1.0*alpha[5]*fr[13]+1.732050807568877*alpha[1]*fr[11]-1.0*alpha[3]*fr[10]+1.732050807568877*alpha[0]*fr[7]-1.0*alpha[1]*fr[6]+1.732050807568877*alpha[4]*fr[5]-1.0*(fr[1]*alpha[4]+alpha[0]*fr[3])+alpha[2]*(1.732050807568877*fr[2]-1.0*fr[0]))*dfac_y; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[15]+alpha[2]*fr[14])-1.0*alpha[4]*fr[13]+1.732050807568877*alpha[1]*fr[12]-1.0*alpha[2]*fr[10]+1.732050807568877*alpha[0]*fr[9]-1.0*alpha[1]*fr[8]+1.732050807568877*alpha[5]*fr[5]-1.0*(fr[1]*alpha[5]+alpha[0]*fr[4])+(1.732050807568877*fr[2]-1.0*fr[0])*alpha[3])*dfac_y; 
  incr[5] = 0.1767766952966368*(3.0*(alpha[3]*fr[12]+alpha[2]*fr[11]+alpha[5]*fr[9])-1.732050807568877*alpha[3]*fr[8]+3.0*alpha[4]*fr[7]-1.732050807568877*alpha[2]*fr[6]+3.0*alpha[0]*fr[5]-1.732050807568877*(fr[4]*alpha[5]+fr[3]*alpha[4])+3.0*alpha[1]*fr[2]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_y; 
  incr[6] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[15]+alpha[5]*fr[14])-1.0*alpha[3]*fr[13]+1.732050807568877*alpha[0]*fr[11]-1.0*alpha[5]*fr[10]+1.732050807568877*alpha[1]*fr[7]-1.0*alpha[0]*fr[6]+1.732050807568877*(alpha[2]*fr[5]+fr[2]*alpha[4])-1.0*(fr[0]*alpha[4]+alpha[1]*fr[3]+fr[1]*alpha[2]))*dfac_y; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[5]*fr[15]+alpha[3]*fr[14])-1.732050807568877*alpha[5]*fr[13]+3.0*alpha[1]*fr[11]-1.732050807568877*alpha[3]*fr[10]+3.0*alpha[0]*fr[7]-1.732050807568877*alpha[1]*fr[6]+3.0*alpha[4]*fr[5]-1.732050807568877*(fr[1]*alpha[4]+alpha[0]*fr[3])+alpha[2]*(3.0*fr[2]-1.732050807568877*fr[0]))*dfac_y; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[4]*fr[14])-1.0*alpha[2]*fr[13]+1.732050807568877*alpha[0]*fr[12]-1.0*alpha[4]*fr[10]+1.732050807568877*alpha[1]*fr[9]-1.0*alpha[0]*fr[8]+1.732050807568877*(alpha[3]*fr[5]+fr[2]*alpha[5])-1.0*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_y; 
  incr[9] = 0.1767766952966368*(3.0*(alpha[4]*fr[15]+alpha[2]*fr[14])-1.732050807568877*alpha[4]*fr[13]+3.0*alpha[1]*fr[12]-1.732050807568877*alpha[2]*fr[10]+3.0*alpha[0]*fr[9]-1.732050807568877*alpha[1]*fr[8]+3.0*alpha[5]*fr[5]-1.732050807568877*(fr[1]*alpha[5]+alpha[0]*fr[4])+(3.0*fr[2]-1.732050807568877*fr[0])*alpha[3])*dfac_y; 
  incr[10] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14])-1.0*alpha[1]*fr[13]+1.732050807568877*(alpha[4]*fr[12]+alpha[5]*fr[11])-1.0*alpha[0]*fr[10]+1.732050807568877*alpha[2]*fr[9]-1.0*alpha[4]*fr[8]+1.732050807568877*alpha[3]*fr[7]-1.0*(alpha[5]*fr[6]+alpha[2]*fr[4]+alpha[3]*fr[3]))*dfac_y; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[3]*fr[15]+alpha[5]*fr[14])-1.732050807568877*alpha[3]*fr[13]+3.0*alpha[0]*fr[11]-1.732050807568877*alpha[5]*fr[10]+3.0*alpha[1]*fr[7]-1.732050807568877*alpha[0]*fr[6]+3.0*(alpha[2]*fr[5]+fr[2]*alpha[4])-1.732050807568877*(fr[0]*alpha[4]+alpha[1]*fr[3]+fr[1]*alpha[2]))*dfac_y; 
  incr[12] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[4]*fr[14])-1.732050807568877*alpha[2]*fr[13]+3.0*alpha[0]*fr[12]-1.732050807568877*alpha[4]*fr[10]+3.0*alpha[1]*fr[9]-1.732050807568877*alpha[0]*fr[8]+3.0*(alpha[3]*fr[5]+fr[2]*alpha[5])-1.732050807568877*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_y; 
  incr[13] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14])-1.0*alpha[0]*fr[13]+1.732050807568877*(alpha[2]*fr[12]+alpha[3]*fr[11])-1.0*alpha[1]*fr[10]+1.732050807568877*alpha[4]*fr[9]-1.0*alpha[2]*fr[8]+1.732050807568877*alpha[5]*fr[7]-1.0*(alpha[3]*fr[6]+fr[3]*alpha[5]+alpha[4]*fr[4]))*dfac_y; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14])-1.732050807568877*alpha[1]*fr[13]+3.0*(alpha[4]*fr[12]+alpha[5]*fr[11])-1.732050807568877*alpha[0]*fr[10]+3.0*alpha[2]*fr[9]-1.732050807568877*alpha[4]*fr[8]+3.0*alpha[3]*fr[7]-1.732050807568877*(alpha[5]*fr[6]+alpha[2]*fr[4]+alpha[3]*fr[3]))*dfac_y; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14])-1.732050807568877*alpha[0]*fr[13]+3.0*(alpha[2]*fr[12]+alpha[3]*fr[11])-1.732050807568877*alpha[1]*fr[10]+3.0*alpha[4]*fr[9]-1.732050807568877*alpha[2]*fr[8]+3.0*alpha[5]*fr[7]-1.732050807568877*(alpha[3]*fr[6]+fr[3]*alpha[5]+alpha[4]*fr[4]))*dfac_y; 

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
double EmGyrokineticSurf2x2vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
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
  double alpha0 = -(0.003125*(dfac_v*(30.0*Bmag[1]*dfac_x*dfac_y*(((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Phi[2])*m_*wv+((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2])*wm)+(3.0*(((9.0*Apar[1]*Bmag[1]*BmagInv[1]+5.0*BmagInv[0]*(Apar[0]*Bmag[1]-2.0*Apar[1]))*geoZ[1]+5.0*geoZ[0]*((Apar[0]*Bmag[1]-2.0*Apar[1])*BmagInv[1]+BmagInv[0]*Apar[1]*Bmag[1]))*Phi[3]+5.0*(2.0*Phi[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3]+(((Apar[0]*Bmag[1]-2.0*Apar[1])*BmagInv[1]+BmagInv[0]*Apar[1]*Bmag[1])*geoZ[1]+geoZ[0]*(Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(Apar[0]*Bmag[1]-2.0*Apar[1])))*Phi[2]+2.0*Phi[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2]))*dfac_x*dfac_y+80.0*dApardt[0])*q_)-30.0*Bmag[1]*((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Phi[2])*dfac_x*dfac_y*m_))/(dfac_v*m_); 

  double alpha[8]; 
  alpha[0] = -(0.01767766952966368*(dfac_v*(30.0*Bmag[1]*dfac_x*dfac_y*(((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Phi[2])*m_*wv+((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2])*wm)+((((27.0*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(15.0*Apar[0]*Bmag[1]-30.0*Apar[1]))*geoZ[1]+geoZ[0]*((15.0*Apar[0]*Bmag[1]-30.0*Apar[1])*BmagInv[1]+15.0*BmagInv[0]*Apar[1]*Bmag[1]))*Phi[3]+30.0*Phi[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3]+(((15.0*Apar[0]*Bmag[1]-30.0*Apar[1])*BmagInv[1]+15.0*BmagInv[0]*Apar[1]*Bmag[1])*geoZ[1]+geoZ[0]*(15.0*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(15.0*Apar[0]*Bmag[1]-30.0*Apar[1])))*Phi[2]+30.0*Phi[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2])*dfac_x*dfac_y+80.0*dApardt[0])*q_)-30.0*Bmag[1]*((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Phi[2])*dfac_x*dfac_y*m_))/(dfac_v*m_); 
  alpha[1] = -(0.01767766952966368*(dfac_v*(Bmag[1]*dfac_x*dfac_y*(((54.0*BmagInv[1]*geoZ[1]+30.0*BmagInv[0]*geoZ[0])*Phi[3]+30.0*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[2])*m_*wv+((54.0*BmagInv[1]*geoZ[1]+30.0*BmagInv[0]*geoZ[0])*Apar[3]+30.0*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2])*wm)+(((((27.0*Apar[0]*Bmag[1]-54.0*Apar[1])*BmagInv[1]+27.0*BmagInv[0]*Apar[1]*Bmag[1])*geoZ[1]+geoZ[0]*(27.0*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(15.0*Apar[0]*Bmag[1]-30.0*Apar[1])))*Phi[3]+Phi[1]*(54.0*BmagInv[1]*geoZ[1]+30.0*BmagInv[0]*geoZ[0])*Apar[3]+((27.0*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(15.0*Apar[0]*Bmag[1]-30.0*Apar[1]))*geoZ[1]+geoZ[0]*((15.0*Apar[0]*Bmag[1]-30.0*Apar[1])*BmagInv[1]+15.0*BmagInv[0]*Apar[1]*Bmag[1]))*Phi[2]+30.0*Phi[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2])*dfac_x*dfac_y+80.0*dApardt[1])*q_)+Bmag[1]*(((-54.0*BmagInv[1]*geoZ[1])-30.0*BmagInv[0]*geoZ[0])*Phi[3]-30.0*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[2])*dfac_x*dfac_y*m_))/(dfac_v*m_); 
  alpha[2] = -(0.01767766952966368*(((Bmag[1]*(27.0*BmagInv[1]*geoZ[1]+15.0*BmagInv[0]*geoZ[0])*Apar[3]+((30.0*BmagInv[1]+15.0*BmagInv[0]*Bmag[1])*geoZ[1]+geoZ[0]*(15.0*Bmag[1]*BmagInv[1]+30.0*BmagInv[0]))*Apar[2])*Phi[3]+Phi[2]*(((15.0*BmagInv[0]*Bmag[1]-30.0*BmagInv[1])*geoZ[1]+geoZ[0]*(15.0*Bmag[1]*BmagInv[1]-30.0*BmagInv[0]))*Apar[3]+15.0*Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2]))*dfac_x*dfac_y+80.0*dApardt[2])*q_)/m_; 
  alpha[3] = -(0.3061862178478971*Bmag[1]*((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2])*dfac_x*dfac_y)/(dfac_m*m_); 
  alpha[4] = -(0.01767766952966368*(((27.0*Bmag[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3]+((27.0*Bmag[1]*BmagInv[1]+30.0*BmagInv[0])*geoZ[1]+geoZ[0]*(30.0*BmagInv[1]+15.0*BmagInv[0]*Bmag[1]))*Apar[2])*Phi[3]+Phi[2]*(((27.0*Bmag[1]*BmagInv[1]-30.0*BmagInv[0])*geoZ[1]+geoZ[0]*(15.0*BmagInv[0]*Bmag[1]-30.0*BmagInv[1]))*Apar[3]+15.0*Bmag[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2]))*dfac_x*dfac_y+80.0*dApardt[3])*q_)/m_; 
  alpha[5] = -(0.03535533905932736*Bmag[1]*((15.58845726811989*BmagInv[1]*geoZ[1]+8.660254037844386*BmagInv[0]*geoZ[0])*Apar[3]+8.660254037844386*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2])*dfac_x*dfac_y)/(dfac_m*m_); 
  if (alpha0>0) { 
  incr[0] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[13]+alpha[4]*fl[11]+alpha[3]*fl[10])+alpha[5]*fl[8]+1.732050807568877*(alpha[2]*fl[7]+alpha[1]*fl[6])+alpha[4]*fl[5]+alpha[3]*fl[4]+1.732050807568877*alpha[0]*fl[3]+alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[13]+alpha[2]*fl[11]+alpha[5]*fl[10])+alpha[3]*fl[8]+1.732050807568877*(alpha[4]*fl[7]+alpha[0]*fl[6])+alpha[2]*fl[5]+fl[4]*alpha[5]+fl[2]*alpha[4]+1.732050807568877*alpha[1]*fl[3]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = 0.1767766952966368*(1.732050807568877*(alpha[5]*fl[15]+alpha[3]*fl[14])+alpha[5]*fl[12]+1.732050807568877*alpha[1]*fl[11]+alpha[3]*fl[9]+1.732050807568877*(alpha[0]*fl[7]+alpha[4]*fl[6])+alpha[1]*fl[5]+fl[1]*alpha[4]+1.732050807568877*alpha[2]*fl[3]+alpha[0]*fl[2]+fl[0]*alpha[2])*dfac_v; 
  incr[3] = -0.1767766952966368*(3.0*(alpha[5]*fl[13]+alpha[4]*fl[11]+alpha[3]*fl[10])+1.732050807568877*alpha[5]*fl[8]+3.0*(alpha[2]*fl[7]+alpha[1]*fl[6])+1.732050807568877*(alpha[4]*fl[5]+alpha[3]*fl[4])+3.0*alpha[0]*fl[3]+1.732050807568877*(alpha[2]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0]))*dfac_v; 
  incr[4] = 0.1767766952966368*(1.732050807568877*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+alpha[4]*fl[12]+1.732050807568877*alpha[0]*fl[10]+alpha[2]*fl[9]+alpha[1]*fl[8]+alpha[5]*(1.732050807568877*fl[6]+fl[1])+alpha[0]*fl[4]+alpha[3]*(1.732050807568877*fl[3]+fl[0]))*dfac_v; 
  incr[5] = 0.1767766952966368*(1.732050807568877*(alpha[3]*fl[15]+alpha[5]*fl[14])+alpha[3]*fl[12]+1.732050807568877*alpha[0]*fl[11]+alpha[5]*fl[9]+1.732050807568877*(alpha[1]*fl[7]+alpha[2]*fl[6])+alpha[0]*fl[5]+(1.732050807568877*fl[3]+fl[0])*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2])*dfac_v; 
  incr[6] = -0.1767766952966368*(3.0*(alpha[3]*fl[13]+alpha[2]*fl[11]+alpha[5]*fl[10])+1.732050807568877*alpha[3]*fl[8]+3.0*(alpha[4]*fl[7]+alpha[0]*fl[6])+1.732050807568877*(alpha[2]*fl[5]+fl[4]*alpha[5]+fl[2]*alpha[4])+3.0*alpha[1]*fl[3]+1.732050807568877*(alpha[0]*fl[1]+fl[0]*alpha[1]))*dfac_v; 
  incr[7] = -0.1767766952966368*(3.0*(alpha[5]*fl[15]+alpha[3]*fl[14])+1.732050807568877*alpha[5]*fl[12]+3.0*alpha[1]*fl[11]+1.732050807568877*alpha[3]*fl[9]+3.0*(alpha[0]*fl[7]+alpha[4]*fl[6])+1.732050807568877*(alpha[1]*fl[5]+fl[1]*alpha[4])+3.0*alpha[2]*fl[3]+1.732050807568877*(alpha[0]*fl[2]+fl[0]*alpha[2]))*dfac_v; 
  incr[8] = 0.1767766952966368*(1.732050807568877*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+alpha[2]*fl[12]+1.732050807568877*alpha[1]*fl[10]+alpha[4]*fl[9]+alpha[0]*fl[8]+1.732050807568877*alpha[3]*fl[6]+(1.732050807568877*fl[3]+fl[0])*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3])*dfac_v; 
  incr[9] = 0.1767766952966368*(1.732050807568877*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+alpha[1]*fl[12]+1.732050807568877*(alpha[5]*fl[11]+alpha[2]*fl[10])+alpha[0]*fl[9]+alpha[4]*fl[8]+1.732050807568877*alpha[3]*fl[7]+alpha[5]*fl[5]+alpha[2]*fl[4]+fl[2]*alpha[3])*dfac_v; 
  incr[10] = -0.1767766952966368*(3.0*(alpha[4]*fl[15]+alpha[2]*fl[14]+alpha[1]*fl[13])+1.732050807568877*alpha[4]*fl[12]+3.0*alpha[0]*fl[10]+1.732050807568877*(alpha[2]*fl[9]+alpha[1]*fl[8])+3.0*alpha[5]*fl[6]+1.732050807568877*(fl[1]*alpha[5]+alpha[0]*fl[4])+alpha[3]*(3.0*fl[3]+1.732050807568877*fl[0]))*dfac_v; 
  incr[11] = -0.1767766952966368*(3.0*(alpha[3]*fl[15]+alpha[5]*fl[14])+1.732050807568877*alpha[3]*fl[12]+3.0*alpha[0]*fl[11]+1.732050807568877*alpha[5]*fl[9]+3.0*(alpha[1]*fl[7]+alpha[2]*fl[6])+1.732050807568877*alpha[0]*fl[5]+3.0*fl[3]*alpha[4]+1.732050807568877*(fl[0]*alpha[4]+alpha[1]*fl[2]+fl[1]*alpha[2]))*dfac_v; 
  incr[12] = 0.1767766952966368*(1.732050807568877*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+alpha[0]*fl[12]+1.732050807568877*(alpha[3]*fl[11]+alpha[4]*fl[10])+alpha[1]*fl[9]+alpha[2]*fl[8]+1.732050807568877*alpha[5]*fl[7]+alpha[3]*fl[5]+fl[2]*alpha[5]+alpha[4]*fl[4])*dfac_v; 
  incr[13] = -0.1767766952966368*(3.0*(alpha[2]*fl[15]+alpha[4]*fl[14]+alpha[0]*fl[13])+1.732050807568877*alpha[2]*fl[12]+3.0*alpha[1]*fl[10]+1.732050807568877*(alpha[4]*fl[9]+alpha[0]*fl[8])+3.0*(alpha[3]*fl[6]+fl[3]*alpha[5])+1.732050807568877*(fl[0]*alpha[5]+alpha[1]*fl[4]+fl[1]*alpha[3]))*dfac_v; 
  incr[14] = -0.1767766952966368*(3.0*(alpha[1]*fl[15]+alpha[0]*fl[14]+alpha[4]*fl[13])+1.732050807568877*alpha[1]*fl[12]+3.0*(alpha[5]*fl[11]+alpha[2]*fl[10])+1.732050807568877*(alpha[0]*fl[9]+alpha[4]*fl[8])+3.0*alpha[3]*fl[7]+1.732050807568877*(alpha[5]*fl[5]+alpha[2]*fl[4]+fl[2]*alpha[3]))*dfac_v; 
  incr[15] = -0.1767766952966368*(3.0*(alpha[0]*fl[15]+alpha[1]*fl[14]+alpha[2]*fl[13])+1.732050807568877*alpha[0]*fl[12]+3.0*(alpha[3]*fl[11]+alpha[4]*fl[10])+1.732050807568877*(alpha[1]*fl[9]+alpha[2]*fl[8])+3.0*alpha[5]*fl[7]+1.732050807568877*(alpha[3]*fl[5]+fl[2]*alpha[5]+alpha[4]*fl[4]))*dfac_v; 

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
  incr[0] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[13]+alpha[4]*fr[11]+alpha[3]*fr[10])-1.0*alpha[5]*fr[8]+1.732050807568877*(alpha[2]*fr[7]+alpha[1]*fr[6])-1.0*(alpha[4]*fr[5]+alpha[3]*fr[4])+1.732050807568877*alpha[0]*fr[3]-1.0*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[1] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[13]+alpha[2]*fr[11]+alpha[5]*fr[10])-1.0*alpha[3]*fr[8]+1.732050807568877*(alpha[4]*fr[7]+alpha[0]*fr[6])-1.0*(alpha[2]*fr[5]+fr[4]*alpha[5]+fr[2]*alpha[4])+1.732050807568877*alpha[1]*fr[3]-1.0*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[2] = -0.1767766952966368*(1.732050807568877*(alpha[5]*fr[15]+alpha[3]*fr[14])-1.0*alpha[5]*fr[12]+1.732050807568877*alpha[1]*fr[11]-1.0*alpha[3]*fr[9]+1.732050807568877*(alpha[0]*fr[7]+alpha[4]*fr[6])-1.0*(alpha[1]*fr[5]+fr[1]*alpha[4])+1.732050807568877*alpha[2]*fr[3]-1.0*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[3] = 0.1767766952966368*(3.0*(alpha[5]*fr[13]+alpha[4]*fr[11]+alpha[3]*fr[10])-1.732050807568877*alpha[5]*fr[8]+3.0*(alpha[2]*fr[7]+alpha[1]*fr[6])-1.732050807568877*(alpha[4]*fr[5]+alpha[3]*fr[4])+3.0*alpha[0]*fr[3]-1.732050807568877*(alpha[2]*fr[2]+alpha[1]*fr[1]+alpha[0]*fr[0]))*dfac_v; 
  incr[4] = -0.1767766952966368*(1.732050807568877*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.0*alpha[4]*fr[12]+1.732050807568877*alpha[0]*fr[10]-1.0*(alpha[2]*fr[9]+alpha[1]*fr[8])+1.732050807568877*alpha[5]*fr[6]-1.0*(fr[1]*alpha[5]+alpha[0]*fr[4])+alpha[3]*(1.732050807568877*fr[3]-1.0*fr[0]))*dfac_v; 
  incr[5] = -0.1767766952966368*(1.732050807568877*(alpha[3]*fr[15]+alpha[5]*fr[14])-1.0*alpha[3]*fr[12]+1.732050807568877*alpha[0]*fr[11]-1.0*alpha[5]*fr[9]+1.732050807568877*(alpha[1]*fr[7]+alpha[2]*fr[6])-1.0*alpha[0]*fr[5]+1.732050807568877*fr[3]*alpha[4]-1.0*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[6] = 0.1767766952966368*(3.0*(alpha[3]*fr[13]+alpha[2]*fr[11]+alpha[5]*fr[10])-1.732050807568877*alpha[3]*fr[8]+3.0*(alpha[4]*fr[7]+alpha[0]*fr[6])-1.732050807568877*(alpha[2]*fr[5]+fr[4]*alpha[5]+fr[2]*alpha[4])+3.0*alpha[1]*fr[3]-1.732050807568877*(alpha[0]*fr[1]+fr[0]*alpha[1]))*dfac_v; 
  incr[7] = 0.1767766952966368*(3.0*(alpha[5]*fr[15]+alpha[3]*fr[14])-1.732050807568877*alpha[5]*fr[12]+3.0*alpha[1]*fr[11]-1.732050807568877*alpha[3]*fr[9]+3.0*(alpha[0]*fr[7]+alpha[4]*fr[6])-1.732050807568877*(alpha[1]*fr[5]+fr[1]*alpha[4])+3.0*alpha[2]*fr[3]-1.732050807568877*(alpha[0]*fr[2]+fr[0]*alpha[2]))*dfac_v; 
  incr[8] = -0.1767766952966368*(1.732050807568877*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.0*alpha[2]*fr[12]+1.732050807568877*alpha[1]*fr[10]-1.0*(alpha[4]*fr[9]+alpha[0]*fr[8])+1.732050807568877*(alpha[3]*fr[6]+fr[3]*alpha[5])-1.0*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_v; 
  incr[9] = -0.1767766952966368*(1.732050807568877*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.0*alpha[1]*fr[12]+1.732050807568877*(alpha[5]*fr[11]+alpha[2]*fr[10])-1.0*(alpha[0]*fr[9]+alpha[4]*fr[8])+1.732050807568877*alpha[3]*fr[7]-1.0*(alpha[5]*fr[5]+alpha[2]*fr[4]+fr[2]*alpha[3]))*dfac_v; 
  incr[10] = 0.1767766952966368*(3.0*(alpha[4]*fr[15]+alpha[2]*fr[14]+alpha[1]*fr[13])-1.732050807568877*alpha[4]*fr[12]+3.0*alpha[0]*fr[10]-1.732050807568877*(alpha[2]*fr[9]+alpha[1]*fr[8])+3.0*alpha[5]*fr[6]-1.732050807568877*(fr[1]*alpha[5]+alpha[0]*fr[4])+alpha[3]*(3.0*fr[3]-1.732050807568877*fr[0]))*dfac_v; 
  incr[11] = 0.1767766952966368*(3.0*(alpha[3]*fr[15]+alpha[5]*fr[14])-1.732050807568877*alpha[3]*fr[12]+3.0*alpha[0]*fr[11]-1.732050807568877*alpha[5]*fr[9]+3.0*(alpha[1]*fr[7]+alpha[2]*fr[6])-1.732050807568877*alpha[0]*fr[5]+3.0*fr[3]*alpha[4]-1.732050807568877*(fr[0]*alpha[4]+alpha[1]*fr[2]+fr[1]*alpha[2]))*dfac_v; 
  incr[12] = -0.1767766952966368*(1.732050807568877*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.0*alpha[0]*fr[12]+1.732050807568877*(alpha[3]*fr[11]+alpha[4]*fr[10])-1.0*(alpha[1]*fr[9]+alpha[2]*fr[8])+1.732050807568877*alpha[5]*fr[7]-1.0*(alpha[3]*fr[5]+fr[2]*alpha[5]+alpha[4]*fr[4]))*dfac_v; 
  incr[13] = 0.1767766952966368*(3.0*(alpha[2]*fr[15]+alpha[4]*fr[14]+alpha[0]*fr[13])-1.732050807568877*alpha[2]*fr[12]+3.0*alpha[1]*fr[10]-1.732050807568877*(alpha[4]*fr[9]+alpha[0]*fr[8])+3.0*(alpha[3]*fr[6]+fr[3]*alpha[5])-1.732050807568877*(fr[0]*alpha[5]+alpha[1]*fr[4]+fr[1]*alpha[3]))*dfac_v; 
  incr[14] = 0.1767766952966368*(3.0*(alpha[1]*fr[15]+alpha[0]*fr[14]+alpha[4]*fr[13])-1.732050807568877*alpha[1]*fr[12]+3.0*(alpha[5]*fr[11]+alpha[2]*fr[10])-1.732050807568877*(alpha[0]*fr[9]+alpha[4]*fr[8])+3.0*alpha[3]*fr[7]-1.732050807568877*(alpha[5]*fr[5]+alpha[2]*fr[4]+fr[2]*alpha[3]))*dfac_v; 
  incr[15] = 0.1767766952966368*(3.0*(alpha[0]*fr[15]+alpha[1]*fr[14]+alpha[2]*fr[13])-1.732050807568877*alpha[0]*fr[12]+3.0*(alpha[3]*fr[11]+alpha[4]*fr[10])-1.732050807568877*(alpha[1]*fr[9]+alpha[2]*fr[8])+3.0*alpha[5]*fr[7]-1.732050807568877*(alpha[3]*fr[5]+fr[2]*alpha[5]+alpha[4]*fr[4]))*dfac_v; 

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
