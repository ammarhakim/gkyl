#include <GyrokineticModDecl.h> 
double GyrokineticSurf2x2vSer_X_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[16]; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.5*((5.196152422706631*BmagInv[1]-3.0*BmagInv[0])*Phi[3]+(1.732050807568877*BmagInv[0]-3.0*BmagInv[1])*Phi[2])*dfac_y;

  if (alpha0>0) { 
  incr[0] = -0.125*(9.0*BmagInv[1]*fl[1]*Phi[3]-5.196152422706631*BmagInv[0]*fl[1]*Phi[3]+5.196152422706631*fl[0]*BmagInv[1]*Phi[3]-3.0*BmagInv[0]*fl[0]*Phi[3]-5.196152422706631*BmagInv[1]*fl[1]*Phi[2]+3.0*BmagInv[0]*fl[1]*Phi[2]-3.0*fl[0]*BmagInv[1]*Phi[2]+1.732050807568877*BmagInv[0]*fl[0]*Phi[2])*dfac_x*dfac_y; 
  incr[1] = 0.375*(5.196152422706631*BmagInv[1]*fl[1]*Phi[3]-3.0*BmagInv[0]*fl[1]*Phi[3]+3.0*fl[0]*BmagInv[1]*Phi[3]-1.732050807568877*BmagInv[0]*fl[0]*Phi[3]-3.0*BmagInv[1]*fl[1]*Phi[2]+1.732050807568877*BmagInv[0]*fl[1]*Phi[2]-1.732050807568877*fl[0]*BmagInv[1]*Phi[2]+BmagInv[0]*fl[0]*Phi[2])*dfac_x*dfac_y; 
  incr[2] = -0.125*(9.0*BmagInv[1]*Phi[3]*fl[5]-5.196152422706631*BmagInv[0]*Phi[3]*fl[5]-5.196152422706631*BmagInv[1]*Phi[2]*fl[5]+3.0*BmagInv[0]*Phi[2]*fl[5]+5.196152422706631*BmagInv[1]*fl[2]*Phi[3]-3.0*BmagInv[0]*fl[2]*Phi[3]-3.0*BmagInv[1]*Phi[2]*fl[2]+1.732050807568877*BmagInv[0]*Phi[2]*fl[2])*dfac_x*dfac_y; 
  incr[3] = -0.125*(9.0*BmagInv[1]*Phi[3]*fl[6]-5.196152422706631*BmagInv[0]*Phi[3]*fl[6]-5.196152422706631*BmagInv[1]*Phi[2]*fl[6]+3.0*BmagInv[0]*Phi[2]*fl[6]+5.196152422706631*BmagInv[1]*Phi[3]*fl[3]-3.0*BmagInv[0]*Phi[3]*fl[3]-3.0*BmagInv[1]*Phi[2]*fl[3]+1.732050807568877*BmagInv[0]*Phi[2]*fl[3])*dfac_x*dfac_y; 
  incr[4] = -0.125*(9.0*BmagInv[1]*Phi[3]*fl[8]-5.196152422706631*BmagInv[0]*Phi[3]*fl[8]-5.196152422706631*BmagInv[1]*Phi[2]*fl[8]+3.0*BmagInv[0]*Phi[2]*fl[8]+5.196152422706631*BmagInv[1]*Phi[3]*fl[4]-3.0*BmagInv[0]*Phi[3]*fl[4]-3.0*BmagInv[1]*Phi[2]*fl[4]+1.732050807568877*BmagInv[0]*Phi[2]*fl[4])*dfac_x*dfac_y; 
  incr[5] = 0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fl[5]-3.0*BmagInv[0]*Phi[3]*fl[5]-3.0*BmagInv[1]*Phi[2]*fl[5]+1.732050807568877*BmagInv[0]*Phi[2]*fl[5]+3.0*BmagInv[1]*fl[2]*Phi[3]-1.732050807568877*BmagInv[0]*fl[2]*Phi[3]-1.732050807568877*BmagInv[1]*Phi[2]*fl[2]+BmagInv[0]*Phi[2]*fl[2])*dfac_x*dfac_y; 
  incr[6] = 0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fl[6]-3.0*BmagInv[0]*Phi[3]*fl[6]-3.0*BmagInv[1]*Phi[2]*fl[6]+1.732050807568877*BmagInv[0]*Phi[2]*fl[6]+3.0*BmagInv[1]*Phi[3]*fl[3]-1.732050807568877*BmagInv[0]*Phi[3]*fl[3]-1.732050807568877*BmagInv[1]*Phi[2]*fl[3]+BmagInv[0]*Phi[2]*fl[3])*dfac_x*dfac_y; 
  incr[7] = -0.125*(9.0*BmagInv[1]*Phi[3]*fl[11]-5.196152422706631*BmagInv[0]*Phi[3]*fl[11]-5.196152422706631*BmagInv[1]*Phi[2]*fl[11]+3.0*BmagInv[0]*Phi[2]*fl[11]+5.196152422706631*BmagInv[1]*Phi[3]*fl[7]-3.0*BmagInv[0]*Phi[3]*fl[7]-3.0*BmagInv[1]*Phi[2]*fl[7]+1.732050807568877*BmagInv[0]*Phi[2]*fl[7])*dfac_x*dfac_y; 
  incr[8] = 0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fl[8]-3.0*BmagInv[0]*Phi[3]*fl[8]-3.0*BmagInv[1]*Phi[2]*fl[8]+1.732050807568877*BmagInv[0]*Phi[2]*fl[8]+3.0*BmagInv[1]*Phi[3]*fl[4]-1.732050807568877*BmagInv[0]*Phi[3]*fl[4]-1.732050807568877*BmagInv[1]*Phi[2]*fl[4]+BmagInv[0]*Phi[2]*fl[4])*dfac_x*dfac_y; 
  incr[9] = -0.125*(9.0*BmagInv[1]*Phi[3]*fl[12]-5.196152422706631*BmagInv[0]*Phi[3]*fl[12]-5.196152422706631*BmagInv[1]*Phi[2]*fl[12]+3.0*BmagInv[0]*Phi[2]*fl[12]+5.196152422706631*BmagInv[1]*Phi[3]*fl[9]-3.0*BmagInv[0]*Phi[3]*fl[9]-3.0*BmagInv[1]*Phi[2]*fl[9]+1.732050807568877*BmagInv[0]*Phi[2]*fl[9])*dfac_x*dfac_y; 
  incr[10] = -0.125*(9.0*BmagInv[1]*Phi[3]*fl[13]-5.196152422706631*BmagInv[0]*Phi[3]*fl[13]-5.196152422706631*BmagInv[1]*Phi[2]*fl[13]+3.0*BmagInv[0]*Phi[2]*fl[13]+5.196152422706631*BmagInv[1]*Phi[3]*fl[10]-3.0*BmagInv[0]*Phi[3]*fl[10]-3.0*BmagInv[1]*Phi[2]*fl[10]+1.732050807568877*BmagInv[0]*Phi[2]*fl[10])*dfac_x*dfac_y; 
  incr[11] = 0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fl[11]-3.0*BmagInv[0]*Phi[3]*fl[11]-3.0*BmagInv[1]*Phi[2]*fl[11]+1.732050807568877*BmagInv[0]*Phi[2]*fl[11]+3.0*BmagInv[1]*Phi[3]*fl[7]-1.732050807568877*BmagInv[0]*Phi[3]*fl[7]-1.732050807568877*BmagInv[1]*Phi[2]*fl[7]+BmagInv[0]*Phi[2]*fl[7])*dfac_x*dfac_y; 
  incr[12] = 0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fl[12]-3.0*BmagInv[0]*Phi[3]*fl[12]-3.0*BmagInv[1]*Phi[2]*fl[12]+1.732050807568877*BmagInv[0]*Phi[2]*fl[12]+3.0*BmagInv[1]*Phi[3]*fl[9]-1.732050807568877*BmagInv[0]*Phi[3]*fl[9]-1.732050807568877*BmagInv[1]*Phi[2]*fl[9]+BmagInv[0]*Phi[2]*fl[9])*dfac_x*dfac_y; 
  incr[13] = 0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fl[13]-3.0*BmagInv[0]*Phi[3]*fl[13]-3.0*BmagInv[1]*Phi[2]*fl[13]+1.732050807568877*BmagInv[0]*Phi[2]*fl[13]+3.0*BmagInv[1]*Phi[3]*fl[10]-1.732050807568877*BmagInv[0]*Phi[3]*fl[10]-1.732050807568877*BmagInv[1]*Phi[2]*fl[10]+BmagInv[0]*Phi[2]*fl[10])*dfac_x*dfac_y; 
  incr[14] = -0.125*(9.0*BmagInv[1]*Phi[3]*fl[15]-5.196152422706631*BmagInv[0]*Phi[3]*fl[15]-5.196152422706631*BmagInv[1]*Phi[2]*fl[15]+3.0*BmagInv[0]*Phi[2]*fl[15]+5.196152422706631*BmagInv[1]*Phi[3]*fl[14]-3.0*BmagInv[0]*Phi[3]*fl[14]-3.0*BmagInv[1]*Phi[2]*fl[14]+1.732050807568877*BmagInv[0]*Phi[2]*fl[14])*dfac_x*dfac_y; 
  incr[15] = 0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fl[15]-3.0*BmagInv[0]*Phi[3]*fl[15]-3.0*BmagInv[1]*Phi[2]*fl[15]+1.732050807568877*BmagInv[0]*Phi[2]*fl[15]+3.0*BmagInv[1]*Phi[3]*fl[14]-1.732050807568877*BmagInv[0]*Phi[3]*fl[14]-1.732050807568877*BmagInv[1]*Phi[2]*fl[14]+BmagInv[0]*Phi[2]*fl[14])*dfac_x*dfac_y; 

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
  incr[0] = 0.125*(9.0*BmagInv[1]*fr[1]*Phi[3]-5.196152422706631*BmagInv[0]*fr[1]*Phi[3]-5.196152422706631*fr[0]*BmagInv[1]*Phi[3]+3.0*BmagInv[0]*fr[0]*Phi[3]-5.196152422706631*BmagInv[1]*fr[1]*Phi[2]+3.0*BmagInv[0]*fr[1]*Phi[2]+3.0*fr[0]*BmagInv[1]*Phi[2]-1.732050807568877*BmagInv[0]*fr[0]*Phi[2])*dfac_x*dfac_y; 
  incr[1] = -0.375*(5.196152422706631*BmagInv[1]*fr[1]*Phi[3]-3.0*BmagInv[0]*fr[1]*Phi[3]-3.0*fr[0]*BmagInv[1]*Phi[3]+1.732050807568877*BmagInv[0]*fr[0]*Phi[3]-3.0*BmagInv[1]*fr[1]*Phi[2]+1.732050807568877*BmagInv[0]*fr[1]*Phi[2]+1.732050807568877*fr[0]*BmagInv[1]*Phi[2]-1.0*BmagInv[0]*fr[0]*Phi[2])*dfac_x*dfac_y; 
  incr[2] = 0.125*(9.0*BmagInv[1]*Phi[3]*fr[5]-5.196152422706631*BmagInv[0]*Phi[3]*fr[5]-5.196152422706631*BmagInv[1]*Phi[2]*fr[5]+3.0*BmagInv[0]*Phi[2]*fr[5]-5.196152422706631*BmagInv[1]*fr[2]*Phi[3]+3.0*BmagInv[0]*fr[2]*Phi[3]+3.0*BmagInv[1]*Phi[2]*fr[2]-1.732050807568877*BmagInv[0]*Phi[2]*fr[2])*dfac_x*dfac_y; 
  incr[3] = 0.125*(9.0*BmagInv[1]*Phi[3]*fr[6]-5.196152422706631*BmagInv[0]*Phi[3]*fr[6]-5.196152422706631*BmagInv[1]*Phi[2]*fr[6]+3.0*BmagInv[0]*Phi[2]*fr[6]-5.196152422706631*BmagInv[1]*Phi[3]*fr[3]+3.0*BmagInv[0]*Phi[3]*fr[3]+3.0*BmagInv[1]*Phi[2]*fr[3]-1.732050807568877*BmagInv[0]*Phi[2]*fr[3])*dfac_x*dfac_y; 
  incr[4] = 0.125*(9.0*BmagInv[1]*Phi[3]*fr[8]-5.196152422706631*BmagInv[0]*Phi[3]*fr[8]-5.196152422706631*BmagInv[1]*Phi[2]*fr[8]+3.0*BmagInv[0]*Phi[2]*fr[8]-5.196152422706631*BmagInv[1]*Phi[3]*fr[4]+3.0*BmagInv[0]*Phi[3]*fr[4]+3.0*BmagInv[1]*Phi[2]*fr[4]-1.732050807568877*BmagInv[0]*Phi[2]*fr[4])*dfac_x*dfac_y; 
  incr[5] = -0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fr[5]-3.0*BmagInv[0]*Phi[3]*fr[5]-3.0*BmagInv[1]*Phi[2]*fr[5]+1.732050807568877*BmagInv[0]*Phi[2]*fr[5]-3.0*BmagInv[1]*fr[2]*Phi[3]+1.732050807568877*BmagInv[0]*fr[2]*Phi[3]+1.732050807568877*BmagInv[1]*Phi[2]*fr[2]-1.0*BmagInv[0]*Phi[2]*fr[2])*dfac_x*dfac_y; 
  incr[6] = -0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fr[6]-3.0*BmagInv[0]*Phi[3]*fr[6]-3.0*BmagInv[1]*Phi[2]*fr[6]+1.732050807568877*BmagInv[0]*Phi[2]*fr[6]-3.0*BmagInv[1]*Phi[3]*fr[3]+1.732050807568877*BmagInv[0]*Phi[3]*fr[3]+1.732050807568877*BmagInv[1]*Phi[2]*fr[3]-1.0*BmagInv[0]*Phi[2]*fr[3])*dfac_x*dfac_y; 
  incr[7] = 0.125*(9.0*BmagInv[1]*Phi[3]*fr[11]-5.196152422706631*BmagInv[0]*Phi[3]*fr[11]-5.196152422706631*BmagInv[1]*Phi[2]*fr[11]+3.0*BmagInv[0]*Phi[2]*fr[11]-5.196152422706631*BmagInv[1]*Phi[3]*fr[7]+3.0*BmagInv[0]*Phi[3]*fr[7]+3.0*BmagInv[1]*Phi[2]*fr[7]-1.732050807568877*BmagInv[0]*Phi[2]*fr[7])*dfac_x*dfac_y; 
  incr[8] = -0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fr[8]-3.0*BmagInv[0]*Phi[3]*fr[8]-3.0*BmagInv[1]*Phi[2]*fr[8]+1.732050807568877*BmagInv[0]*Phi[2]*fr[8]-3.0*BmagInv[1]*Phi[3]*fr[4]+1.732050807568877*BmagInv[0]*Phi[3]*fr[4]+1.732050807568877*BmagInv[1]*Phi[2]*fr[4]-1.0*BmagInv[0]*Phi[2]*fr[4])*dfac_x*dfac_y; 
  incr[9] = 0.125*(9.0*BmagInv[1]*Phi[3]*fr[12]-5.196152422706631*BmagInv[0]*Phi[3]*fr[12]-5.196152422706631*BmagInv[1]*Phi[2]*fr[12]+3.0*BmagInv[0]*Phi[2]*fr[12]-5.196152422706631*BmagInv[1]*Phi[3]*fr[9]+3.0*BmagInv[0]*Phi[3]*fr[9]+3.0*BmagInv[1]*Phi[2]*fr[9]-1.732050807568877*BmagInv[0]*Phi[2]*fr[9])*dfac_x*dfac_y; 
  incr[10] = 0.125*(9.0*BmagInv[1]*Phi[3]*fr[13]-5.196152422706631*BmagInv[0]*Phi[3]*fr[13]-5.196152422706631*BmagInv[1]*Phi[2]*fr[13]+3.0*BmagInv[0]*Phi[2]*fr[13]-5.196152422706631*BmagInv[1]*Phi[3]*fr[10]+3.0*BmagInv[0]*Phi[3]*fr[10]+3.0*BmagInv[1]*Phi[2]*fr[10]-1.732050807568877*BmagInv[0]*Phi[2]*fr[10])*dfac_x*dfac_y; 
  incr[11] = -0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fr[11]-3.0*BmagInv[0]*Phi[3]*fr[11]-3.0*BmagInv[1]*Phi[2]*fr[11]+1.732050807568877*BmagInv[0]*Phi[2]*fr[11]-3.0*BmagInv[1]*Phi[3]*fr[7]+1.732050807568877*BmagInv[0]*Phi[3]*fr[7]+1.732050807568877*BmagInv[1]*Phi[2]*fr[7]-1.0*BmagInv[0]*Phi[2]*fr[7])*dfac_x*dfac_y; 
  incr[12] = -0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fr[12]-3.0*BmagInv[0]*Phi[3]*fr[12]-3.0*BmagInv[1]*Phi[2]*fr[12]+1.732050807568877*BmagInv[0]*Phi[2]*fr[12]-3.0*BmagInv[1]*Phi[3]*fr[9]+1.732050807568877*BmagInv[0]*Phi[3]*fr[9]+1.732050807568877*BmagInv[1]*Phi[2]*fr[9]-1.0*BmagInv[0]*Phi[2]*fr[9])*dfac_x*dfac_y; 
  incr[13] = -0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fr[13]-3.0*BmagInv[0]*Phi[3]*fr[13]-3.0*BmagInv[1]*Phi[2]*fr[13]+1.732050807568877*BmagInv[0]*Phi[2]*fr[13]-3.0*BmagInv[1]*Phi[3]*fr[10]+1.732050807568877*BmagInv[0]*Phi[3]*fr[10]+1.732050807568877*BmagInv[1]*Phi[2]*fr[10]-1.0*BmagInv[0]*Phi[2]*fr[10])*dfac_x*dfac_y; 
  incr[14] = 0.125*(9.0*BmagInv[1]*Phi[3]*fr[15]-5.196152422706631*BmagInv[0]*Phi[3]*fr[15]-5.196152422706631*BmagInv[1]*Phi[2]*fr[15]+3.0*BmagInv[0]*Phi[2]*fr[15]-5.196152422706631*BmagInv[1]*Phi[3]*fr[14]+3.0*BmagInv[0]*Phi[3]*fr[14]+3.0*BmagInv[1]*Phi[2]*fr[14]-1.732050807568877*BmagInv[0]*Phi[2]*fr[14])*dfac_x*dfac_y; 
  incr[15] = -0.375*(5.196152422706631*BmagInv[1]*Phi[3]*fr[15]-3.0*BmagInv[0]*Phi[3]*fr[15]-3.0*BmagInv[1]*Phi[2]*fr[15]+1.732050807568877*BmagInv[0]*Phi[2]*fr[15]-3.0*BmagInv[1]*Phi[3]*fr[14]+1.732050807568877*BmagInv[0]*Phi[3]*fr[14]+1.732050807568877*BmagInv[1]*Phi[2]*fr[14]-1.0*BmagInv[0]*Phi[2]*fr[14])*dfac_x*dfac_y; 

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
double GyrokineticSurf2x2vSer_Y_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[16]; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.5*(2.0*BcurvY[0]*m_*wv2+1.732050807568877*BmagInv[0]*Bmag[1]*dfac_x*wm+(1.732050807568877*BmagInv[0]*Phi[1]-3.0*BmagInv[0]*Phi[3])*dfac_x*q_))/q_;

  if (alpha0>0) { 
  incr[0] = (0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[1]*fl[5]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[0]*fl[2]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[1]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[0]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[11]*dfac_m*m_*wv+6.0*BcurvY[0]*fl[7]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[6]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[3]*dfac_m*m_*wv+9.0*Bmag[1]*BmagInv[1]*fl[5]*dfac_m*dfac_v*dfac_x*wm+9.0*BmagInv[0]*Bmag[1]*fl[2]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[1]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*fl[0]*Bmag[1]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[1]*Phi[3]*fl[5]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fl[5]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[0]*fl[2]*Phi[3]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*fl[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*fl[0]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fl[2]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[1]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*fl[0]*Phi[1]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*Bmag[1]*BmagInv[1]*fl[12]*dfac_v*dfac_x+5.196152422706631*BmagInv[0]*Bmag[1]*fl[9]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[8]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[4]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[1] = (0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[0]*fl[5]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[1]*fl[2]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[1]*dfac_m*dfac_v*m_*wv2+6.0*fl[0]*BcurvY[1]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[11]*dfac_m*m_*wv+6.0*BcurvY[1]*fl[7]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[6]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[3]*dfac_m*m_*wv+9.0*BmagInv[0]*Bmag[1]*fl[5]*dfac_m*dfac_v*dfac_x*wm+9.0*Bmag[1]*BmagInv[1]*fl[2]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[1]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*fl[0]*Bmag[1]*BmagInv[1]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[0]*Phi[3]*fl[5]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fl[5]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[1]*fl[2]*Phi[3]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*fl[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_-9.0*fl[0]*BmagInv[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fl[2]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[1]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*fl[0]*BmagInv[1]*Phi[1]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Bmag[1]*fl[12]*dfac_v*dfac_x+5.196152422706631*Bmag[1]*BmagInv[1]*fl[9]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[8]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[4]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[2] = -(0.125*dfac_y*(6.0*BcurvY[1]*fl[5]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[2]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[1]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[0]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[11]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[7]*dfac_m*m_*wv+2.0*BcurvY[1]*fl[6]*dfac_m*m_*wv+2.0*BcurvY[0]*fl[3]*dfac_m*m_*wv+5.196152422706631*Bmag[1]*BmagInv[1]*fl[5]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[2]*dfac_m*dfac_v*dfac_x*wm+3.0*Bmag[1]*BmagInv[1]*fl[1]*dfac_m*dfac_v*dfac_x*wm+3.0*BmagInv[0]*fl[0]*Bmag[1]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[1]*Phi[3]*fl[5]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[5]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*fl[2]*Phi[3]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*fl[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*fl[0]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[2]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[1]*Phi[1]*fl[1]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*fl[0]*Phi[1]*dfac_m*dfac_v*dfac_x*q_+3.0*Bmag[1]*BmagInv[1]*fl[12]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[9]*dfac_v*dfac_x+1.732050807568877*Bmag[1]*BmagInv[1]*fl[8]*dfac_v*dfac_x+1.732050807568877*BmagInv[0]*Bmag[1]*fl[4]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[3] = (0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[1]*fl[11]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[0]*fl[7]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[6]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[3]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[5]*dfac_m*m_*wv+6.0*BcurvY[0]*fl[2]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[1]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[0]*dfac_m*m_*wv+9.0*Bmag[1]*BmagInv[1]*fl[11]*dfac_m*dfac_v*dfac_x*wm+9.0*BmagInv[0]*Bmag[1]*fl[7]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[6]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[3]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[1]*Phi[3]*fl[11]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fl[11]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[0]*Phi[3]*fl[7]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fl[7]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fl[6]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[6]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fl[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*Bmag[1]*BmagInv[1]*fl[15]*dfac_v*dfac_x+5.196152422706631*BmagInv[0]*Bmag[1]*fl[14]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[13]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[10]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[4] = (0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[1]*fl[12]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[0]*fl[9]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[8]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[4]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[15]*dfac_m*m_*wv+6.0*BcurvY[0]*fl[14]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[13]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[10]*dfac_m*m_*wv+9.0*Bmag[1]*BmagInv[1]*fl[12]*dfac_m*dfac_v*dfac_x*wm+9.0*BmagInv[0]*Bmag[1]*fl[9]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[8]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[4]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[1]*Phi[3]*fl[12]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fl[12]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[0]*Phi[3]*fl[9]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fl[9]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fl[8]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[8]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fl[4]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[4]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*Bmag[1]*BmagInv[1]*fl[5]*dfac_v*dfac_x+5.196152422706631*BmagInv[0]*Bmag[1]*fl[2]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[1]*dfac_v*dfac_x+3.0*BmagInv[0]*fl[0]*Bmag[1]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[5] = -(0.125*dfac_y*(6.0*BcurvY[0]*fl[5]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[2]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[1]*dfac_m*dfac_v*m_*wv2+3.464101615137754*fl[0]*BcurvY[1]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[11]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[7]*dfac_m*m_*wv+2.0*BcurvY[0]*fl[6]*dfac_m*m_*wv+2.0*BcurvY[1]*fl[3]*dfac_m*m_*wv+5.196152422706631*BmagInv[0]*Bmag[1]*fl[5]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[2]*dfac_m*dfac_v*dfac_x*wm+3.0*BmagInv[0]*Bmag[1]*fl[1]*dfac_m*dfac_v*dfac_x*wm+3.0*fl[0]*Bmag[1]*BmagInv[1]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[0]*Phi[3]*fl[5]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[5]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*fl[2]*Phi[3]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*fl[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*fl[0]*BmagInv[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[2]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Phi[1]*fl[1]*dfac_m*dfac_v*dfac_x*q_+3.0*fl[0]*BmagInv[1]*Phi[1]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Bmag[1]*fl[12]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[9]*dfac_v*dfac_x+1.732050807568877*BmagInv[0]*Bmag[1]*fl[8]*dfac_v*dfac_x+1.732050807568877*Bmag[1]*BmagInv[1]*fl[4]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[6] = (0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[0]*fl[11]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[1]*fl[7]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[6]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[3]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[5]*dfac_m*m_*wv+6.0*BcurvY[1]*fl[2]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[1]*dfac_m*m_*wv+3.464101615137754*fl[0]*BcurvY[1]*dfac_m*m_*wv+9.0*BmagInv[0]*Bmag[1]*fl[11]*dfac_m*dfac_v*dfac_x*wm+9.0*Bmag[1]*BmagInv[1]*fl[7]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[6]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[3]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[0]*Phi[3]*fl[11]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fl[11]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[1]*Phi[3]*fl[7]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fl[7]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fl[6]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[6]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fl[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Bmag[1]*fl[15]*dfac_v*dfac_x+5.196152422706631*Bmag[1]*BmagInv[1]*fl[14]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[13]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[10]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[7] = -(0.125*dfac_y*(6.0*BcurvY[1]*fl[11]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[7]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[6]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[3]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[5]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[2]*dfac_m*m_*wv+2.0*BcurvY[1]*fl[1]*dfac_m*m_*wv+2.0*BcurvY[0]*fl[0]*dfac_m*m_*wv+5.196152422706631*Bmag[1]*BmagInv[1]*fl[11]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[7]*dfac_m*dfac_v*dfac_x*wm+3.0*Bmag[1]*BmagInv[1]*fl[6]*dfac_m*dfac_v*dfac_x*wm+3.0*BmagInv[0]*Bmag[1]*fl[3]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[1]*Phi[3]*fl[11]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[11]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fl[7]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[7]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[3]*fl[6]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[1]*Phi[1]*fl[6]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[3]*fl[3]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Phi[1]*fl[3]*dfac_m*dfac_v*dfac_x*q_+3.0*Bmag[1]*BmagInv[1]*fl[15]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[14]*dfac_v*dfac_x+1.732050807568877*Bmag[1]*BmagInv[1]*fl[13]*dfac_v*dfac_x+1.732050807568877*BmagInv[0]*Bmag[1]*fl[10]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[8] = (0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[0]*fl[12]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[1]*fl[9]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[8]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[4]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[15]*dfac_m*m_*wv+6.0*BcurvY[1]*fl[14]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[13]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[10]*dfac_m*m_*wv+9.0*BmagInv[0]*Bmag[1]*fl[12]*dfac_m*dfac_v*dfac_x*wm+9.0*Bmag[1]*BmagInv[1]*fl[9]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[8]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[4]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[0]*Phi[3]*fl[12]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fl[12]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[1]*Phi[3]*fl[9]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fl[9]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fl[8]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[8]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fl[4]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[4]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Bmag[1]*fl[5]*dfac_v*dfac_x+5.196152422706631*Bmag[1]*BmagInv[1]*fl[2]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[1]*dfac_v*dfac_x+3.0*fl[0]*Bmag[1]*BmagInv[1]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[9] = -(0.125*dfac_y*(6.0*BcurvY[1]*fl[12]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[9]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[8]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[4]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[15]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[14]*dfac_m*m_*wv+2.0*BcurvY[1]*fl[13]*dfac_m*m_*wv+2.0*BcurvY[0]*fl[10]*dfac_m*m_*wv+5.196152422706631*Bmag[1]*BmagInv[1]*fl[12]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[9]*dfac_m*dfac_v*dfac_x*wm+3.0*Bmag[1]*BmagInv[1]*fl[8]*dfac_m*dfac_v*dfac_x*wm+3.0*BmagInv[0]*Bmag[1]*fl[4]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[1]*Phi[3]*fl[12]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[12]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fl[9]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[9]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[3]*fl[8]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[1]*Phi[1]*fl[8]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[3]*fl[4]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Phi[1]*fl[4]*dfac_m*dfac_v*dfac_x*q_+3.0*Bmag[1]*BmagInv[1]*fl[5]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[2]*dfac_v*dfac_x+1.732050807568877*Bmag[1]*BmagInv[1]*fl[1]*dfac_v*dfac_x+1.732050807568877*BmagInv[0]*fl[0]*Bmag[1]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[10] = (0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[1]*fl[15]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[0]*fl[14]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[13]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[10]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[12]*dfac_m*m_*wv+6.0*BcurvY[0]*fl[9]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[8]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[4]*dfac_m*m_*wv+9.0*Bmag[1]*BmagInv[1]*fl[15]*dfac_m*dfac_v*dfac_x*wm+9.0*BmagInv[0]*Bmag[1]*fl[14]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[13]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[10]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[1]*Phi[3]*fl[15]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fl[15]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[0]*Phi[3]*fl[14]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fl[14]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fl[13]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[13]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fl[10]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[10]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*Bmag[1]*BmagInv[1]*fl[11]*dfac_v*dfac_x+5.196152422706631*BmagInv[0]*Bmag[1]*fl[7]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[6]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[3]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[11] = -(0.125*dfac_y*(6.0*BcurvY[0]*fl[11]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[7]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[6]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[3]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[5]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[2]*dfac_m*m_*wv+2.0*BcurvY[0]*fl[1]*dfac_m*m_*wv+2.0*fl[0]*BcurvY[1]*dfac_m*m_*wv+5.196152422706631*BmagInv[0]*Bmag[1]*fl[11]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[7]*dfac_m*dfac_v*dfac_x*wm+3.0*BmagInv[0]*Bmag[1]*fl[6]*dfac_m*dfac_v*dfac_x*wm+3.0*Bmag[1]*BmagInv[1]*fl[3]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[0]*Phi[3]*fl[11]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[11]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fl[7]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[7]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[3]*fl[6]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Phi[1]*fl[6]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[3]*fl[3]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[1]*Phi[1]*fl[3]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Bmag[1]*fl[15]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[14]*dfac_v*dfac_x+1.732050807568877*BmagInv[0]*Bmag[1]*fl[13]*dfac_v*dfac_x+1.732050807568877*Bmag[1]*BmagInv[1]*fl[10]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[12] = -(0.125*dfac_y*(6.0*BcurvY[0]*fl[12]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[9]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[8]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[4]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[15]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[14]*dfac_m*m_*wv+2.0*BcurvY[0]*fl[13]*dfac_m*m_*wv+2.0*BcurvY[1]*fl[10]*dfac_m*m_*wv+5.196152422706631*BmagInv[0]*Bmag[1]*fl[12]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[9]*dfac_m*dfac_v*dfac_x*wm+3.0*BmagInv[0]*Bmag[1]*fl[8]*dfac_m*dfac_v*dfac_x*wm+3.0*Bmag[1]*BmagInv[1]*fl[4]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[0]*Phi[3]*fl[12]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[12]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fl[9]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[9]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[3]*fl[8]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Phi[1]*fl[8]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[3]*fl[4]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[1]*Phi[1]*fl[4]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Bmag[1]*fl[5]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[2]*dfac_v*dfac_x+1.732050807568877*BmagInv[0]*Bmag[1]*fl[1]*dfac_v*dfac_x+1.732050807568877*fl[0]*Bmag[1]*BmagInv[1]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[13] = (0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[0]*fl[15]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[1]*fl[14]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[13]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[10]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[12]*dfac_m*m_*wv+6.0*BcurvY[1]*fl[9]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[8]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[4]*dfac_m*m_*wv+9.0*BmagInv[0]*Bmag[1]*fl[15]*dfac_m*dfac_v*dfac_x*wm+9.0*Bmag[1]*BmagInv[1]*fl[14]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[13]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[10]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[0]*Phi[3]*fl[15]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fl[15]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[1]*Phi[3]*fl[14]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fl[14]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fl[13]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[13]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fl[10]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[10]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Bmag[1]*fl[11]*dfac_v*dfac_x+5.196152422706631*Bmag[1]*BmagInv[1]*fl[7]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[6]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[3]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[14] = -(0.125*dfac_y*(6.0*BcurvY[1]*fl[15]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fl[14]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[13]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[10]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[12]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fl[9]*dfac_m*m_*wv+2.0*BcurvY[1]*fl[8]*dfac_m*m_*wv+2.0*BcurvY[0]*fl[4]*dfac_m*m_*wv+5.196152422706631*Bmag[1]*BmagInv[1]*fl[15]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fl[14]*dfac_m*dfac_v*dfac_x*wm+3.0*Bmag[1]*BmagInv[1]*fl[13]*dfac_m*dfac_v*dfac_x*wm+3.0*BmagInv[0]*Bmag[1]*fl[10]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[1]*Phi[3]*fl[15]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[15]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fl[14]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[14]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[3]*fl[13]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[1]*Phi[1]*fl[13]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[3]*fl[10]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Phi[1]*fl[10]*dfac_m*dfac_v*dfac_x*q_+3.0*Bmag[1]*BmagInv[1]*fl[11]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fl[7]*dfac_v*dfac_x+1.732050807568877*Bmag[1]*BmagInv[1]*fl[6]*dfac_v*dfac_x+1.732050807568877*BmagInv[0]*Bmag[1]*fl[3]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[15] = -(0.125*dfac_y*(6.0*BcurvY[0]*fl[15]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fl[14]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[13]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fl[10]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fl[12]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fl[9]*dfac_m*m_*wv+2.0*BcurvY[0]*fl[8]*dfac_m*m_*wv+2.0*BcurvY[1]*fl[4]*dfac_m*m_*wv+5.196152422706631*BmagInv[0]*Bmag[1]*fl[15]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fl[14]*dfac_m*dfac_v*dfac_x*wm+3.0*BmagInv[0]*Bmag[1]*fl[13]*dfac_m*dfac_v*dfac_x*wm+3.0*Bmag[1]*BmagInv[1]*fl[10]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[0]*Phi[3]*fl[15]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fl[15]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fl[14]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fl[14]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[3]*fl[13]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Phi[1]*fl[13]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[3]*fl[10]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[1]*Phi[1]*fl[10]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Bmag[1]*fl[11]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fl[7]*dfac_v*dfac_x+1.732050807568877*BmagInv[0]*Bmag[1]*fl[6]*dfac_v*dfac_x+1.732050807568877*Bmag[1]*BmagInv[1]*fl[3]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 

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
  incr[0] = -(0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[1]*fr[5]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[0]*fr[2]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[1]*fr[1]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[0]*fr[0]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fr[11]*dfac_m*m_*wv+6.0*BcurvY[0]*fr[7]*dfac_m*m_*wv-3.464101615137754*BcurvY[1]*fr[6]*dfac_m*m_*wv-3.464101615137754*BcurvY[0]*fr[3]*dfac_m*m_*wv+9.0*Bmag[1]*BmagInv[1]*fr[5]*dfac_m*dfac_v*dfac_x*wm+9.0*BmagInv[0]*Bmag[1]*fr[2]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*Bmag[1]*BmagInv[1]*fr[1]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*BmagInv[0]*fr[0]*Bmag[1]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[1]*Phi[3]*fr[5]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fr[5]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[0]*fr[2]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*fr[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*fr[0]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fr[2]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[1]*fr[1]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*fr[0]*Phi[1]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*Bmag[1]*BmagInv[1]*fr[12]*dfac_v*dfac_x+5.196152422706631*BmagInv[0]*Bmag[1]*fr[9]*dfac_v*dfac_x-3.0*Bmag[1]*BmagInv[1]*fr[8]*dfac_v*dfac_x-3.0*BmagInv[0]*Bmag[1]*fr[4]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[1] = -(0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[0]*fr[5]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[1]*fr[2]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[0]*fr[1]*dfac_m*dfac_v*m_*wv2-6.0*fr[0]*BcurvY[1]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fr[11]*dfac_m*m_*wv+6.0*BcurvY[1]*fr[7]*dfac_m*m_*wv-3.464101615137754*BcurvY[0]*fr[6]*dfac_m*m_*wv-3.464101615137754*BcurvY[1]*fr[3]*dfac_m*m_*wv+9.0*BmagInv[0]*Bmag[1]*fr[5]*dfac_m*dfac_v*dfac_x*wm+9.0*Bmag[1]*BmagInv[1]*fr[2]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*BmagInv[0]*Bmag[1]*fr[1]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*fr[0]*Bmag[1]*BmagInv[1]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[0]*Phi[3]*fr[5]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fr[5]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[1]*fr[2]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*fr[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+9.0*fr[0]*BmagInv[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fr[2]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[1]*fr[1]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*fr[0]*BmagInv[1]*Phi[1]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Bmag[1]*fr[12]*dfac_v*dfac_x+5.196152422706631*Bmag[1]*BmagInv[1]*fr[9]*dfac_v*dfac_x-3.0*BmagInv[0]*Bmag[1]*fr[8]*dfac_v*dfac_x-3.0*Bmag[1]*BmagInv[1]*fr[4]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[2] = (0.125*dfac_y*(6.0*BcurvY[1]*fr[5]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fr[2]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[1]*fr[1]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[0]*fr[0]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fr[11]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fr[7]*dfac_m*m_*wv-2.0*BcurvY[1]*fr[6]*dfac_m*m_*wv-2.0*BcurvY[0]*fr[3]*dfac_m*m_*wv+5.196152422706631*Bmag[1]*BmagInv[1]*fr[5]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fr[2]*dfac_m*dfac_v*dfac_x*wm-3.0*Bmag[1]*BmagInv[1]*fr[1]*dfac_m*dfac_v*dfac_x*wm-3.0*BmagInv[0]*fr[0]*Bmag[1]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[1]*Phi[3]*fr[5]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fr[5]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*fr[2]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*fr[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*fr[0]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fr[2]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[1]*Phi[1]*fr[1]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[0]*fr[0]*Phi[1]*dfac_m*dfac_v*dfac_x*q_+3.0*Bmag[1]*BmagInv[1]*fr[12]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fr[9]*dfac_v*dfac_x-1.732050807568877*Bmag[1]*BmagInv[1]*fr[8]*dfac_v*dfac_x-1.732050807568877*BmagInv[0]*Bmag[1]*fr[4]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[3] = -(0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[1]*fr[11]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[0]*fr[7]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[1]*fr[6]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[0]*fr[3]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fr[5]*dfac_m*m_*wv+6.0*BcurvY[0]*fr[2]*dfac_m*m_*wv-3.464101615137754*BcurvY[1]*fr[1]*dfac_m*m_*wv-3.464101615137754*BcurvY[0]*fr[0]*dfac_m*m_*wv+9.0*Bmag[1]*BmagInv[1]*fr[11]*dfac_m*dfac_v*dfac_x*wm+9.0*BmagInv[0]*Bmag[1]*fr[7]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*Bmag[1]*BmagInv[1]*fr[6]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*BmagInv[0]*Bmag[1]*fr[3]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[1]*Phi[3]*fr[11]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fr[11]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[0]*Phi[3]*fr[7]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fr[7]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[3]*fr[6]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[1]*fr[6]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[3]*fr[3]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[1]*fr[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*Bmag[1]*BmagInv[1]*fr[15]*dfac_v*dfac_x+5.196152422706631*BmagInv[0]*Bmag[1]*fr[14]*dfac_v*dfac_x-3.0*Bmag[1]*BmagInv[1]*fr[13]*dfac_v*dfac_x-3.0*BmagInv[0]*Bmag[1]*fr[10]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[4] = -(0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[1]*fr[12]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[0]*fr[9]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[1]*fr[8]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[0]*fr[4]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fr[15]*dfac_m*m_*wv+6.0*BcurvY[0]*fr[14]*dfac_m*m_*wv-3.464101615137754*BcurvY[1]*fr[13]*dfac_m*m_*wv-3.464101615137754*BcurvY[0]*fr[10]*dfac_m*m_*wv+9.0*Bmag[1]*BmagInv[1]*fr[12]*dfac_m*dfac_v*dfac_x*wm+9.0*BmagInv[0]*Bmag[1]*fr[9]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*Bmag[1]*BmagInv[1]*fr[8]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*BmagInv[0]*Bmag[1]*fr[4]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[1]*Phi[3]*fr[12]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fr[12]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[0]*Phi[3]*fr[9]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fr[9]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[3]*fr[8]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[1]*fr[8]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[3]*fr[4]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[1]*fr[4]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*Bmag[1]*BmagInv[1]*fr[5]*dfac_v*dfac_x+5.196152422706631*BmagInv[0]*Bmag[1]*fr[2]*dfac_v*dfac_x-3.0*Bmag[1]*BmagInv[1]*fr[1]*dfac_v*dfac_x-3.0*BmagInv[0]*fr[0]*Bmag[1]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[5] = (0.125*dfac_y*(6.0*BcurvY[0]*fr[5]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fr[2]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[0]*fr[1]*dfac_m*dfac_v*m_*wv2-3.464101615137754*fr[0]*BcurvY[1]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fr[11]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fr[7]*dfac_m*m_*wv-2.0*BcurvY[0]*fr[6]*dfac_m*m_*wv-2.0*BcurvY[1]*fr[3]*dfac_m*m_*wv+5.196152422706631*BmagInv[0]*Bmag[1]*fr[5]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fr[2]*dfac_m*dfac_v*dfac_x*wm-3.0*BmagInv[0]*Bmag[1]*fr[1]*dfac_m*dfac_v*dfac_x*wm-3.0*fr[0]*Bmag[1]*BmagInv[1]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[0]*Phi[3]*fr[5]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fr[5]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*fr[2]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*fr[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*fr[0]*BmagInv[1]*Phi[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fr[2]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[0]*Phi[1]*fr[1]*dfac_m*dfac_v*dfac_x*q_-3.0*fr[0]*BmagInv[1]*Phi[1]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Bmag[1]*fr[12]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fr[9]*dfac_v*dfac_x-1.732050807568877*BmagInv[0]*Bmag[1]*fr[8]*dfac_v*dfac_x-1.732050807568877*Bmag[1]*BmagInv[1]*fr[4]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[6] = -(0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[0]*fr[11]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[1]*fr[7]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[0]*fr[6]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[1]*fr[3]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fr[5]*dfac_m*m_*wv+6.0*BcurvY[1]*fr[2]*dfac_m*m_*wv-3.464101615137754*BcurvY[0]*fr[1]*dfac_m*m_*wv-3.464101615137754*fr[0]*BcurvY[1]*dfac_m*m_*wv+9.0*BmagInv[0]*Bmag[1]*fr[11]*dfac_m*dfac_v*dfac_x*wm+9.0*Bmag[1]*BmagInv[1]*fr[7]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*BmagInv[0]*Bmag[1]*fr[6]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*Bmag[1]*BmagInv[1]*fr[3]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[0]*Phi[3]*fr[11]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fr[11]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[1]*Phi[3]*fr[7]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fr[7]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[3]*fr[6]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[1]*fr[6]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[3]*fr[3]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[1]*fr[3]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Bmag[1]*fr[15]*dfac_v*dfac_x+5.196152422706631*Bmag[1]*BmagInv[1]*fr[14]*dfac_v*dfac_x-3.0*BmagInv[0]*Bmag[1]*fr[13]*dfac_v*dfac_x-3.0*Bmag[1]*BmagInv[1]*fr[10]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[7] = (0.125*dfac_y*(6.0*BcurvY[1]*fr[11]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fr[7]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[1]*fr[6]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[0]*fr[3]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fr[5]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fr[2]*dfac_m*m_*wv-2.0*BcurvY[1]*fr[1]*dfac_m*m_*wv-2.0*BcurvY[0]*fr[0]*dfac_m*m_*wv+5.196152422706631*Bmag[1]*BmagInv[1]*fr[11]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fr[7]*dfac_m*dfac_v*dfac_x*wm-3.0*Bmag[1]*BmagInv[1]*fr[6]*dfac_m*dfac_v*dfac_x*wm-3.0*BmagInv[0]*Bmag[1]*fr[3]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[1]*Phi[3]*fr[11]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fr[11]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fr[7]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fr[7]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[3]*fr[6]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[1]*Phi[1]*fr[6]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[3]*fr[3]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[0]*Phi[1]*fr[3]*dfac_m*dfac_v*dfac_x*q_+3.0*Bmag[1]*BmagInv[1]*fr[15]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fr[14]*dfac_v*dfac_x-1.732050807568877*Bmag[1]*BmagInv[1]*fr[13]*dfac_v*dfac_x-1.732050807568877*BmagInv[0]*Bmag[1]*fr[10]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[8] = -(0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[0]*fr[12]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[1]*fr[9]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[0]*fr[8]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[1]*fr[4]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fr[15]*dfac_m*m_*wv+6.0*BcurvY[1]*fr[14]*dfac_m*m_*wv-3.464101615137754*BcurvY[0]*fr[13]*dfac_m*m_*wv-3.464101615137754*BcurvY[1]*fr[10]*dfac_m*m_*wv+9.0*BmagInv[0]*Bmag[1]*fr[12]*dfac_m*dfac_v*dfac_x*wm+9.0*Bmag[1]*BmagInv[1]*fr[9]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*BmagInv[0]*Bmag[1]*fr[8]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*Bmag[1]*BmagInv[1]*fr[4]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[0]*Phi[3]*fr[12]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fr[12]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[1]*Phi[3]*fr[9]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fr[9]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[3]*fr[8]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[1]*fr[8]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[3]*fr[4]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[1]*fr[4]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Bmag[1]*fr[5]*dfac_v*dfac_x+5.196152422706631*Bmag[1]*BmagInv[1]*fr[2]*dfac_v*dfac_x-3.0*BmagInv[0]*Bmag[1]*fr[1]*dfac_v*dfac_x-3.0*fr[0]*Bmag[1]*BmagInv[1]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[9] = (0.125*dfac_y*(6.0*BcurvY[1]*fr[12]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fr[9]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[1]*fr[8]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[0]*fr[4]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fr[15]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fr[14]*dfac_m*m_*wv-2.0*BcurvY[1]*fr[13]*dfac_m*m_*wv-2.0*BcurvY[0]*fr[10]*dfac_m*m_*wv+5.196152422706631*Bmag[1]*BmagInv[1]*fr[12]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fr[9]*dfac_m*dfac_v*dfac_x*wm-3.0*Bmag[1]*BmagInv[1]*fr[8]*dfac_m*dfac_v*dfac_x*wm-3.0*BmagInv[0]*Bmag[1]*fr[4]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[1]*Phi[3]*fr[12]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fr[12]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fr[9]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fr[9]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[3]*fr[8]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[1]*Phi[1]*fr[8]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[3]*fr[4]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[0]*Phi[1]*fr[4]*dfac_m*dfac_v*dfac_x*q_+3.0*Bmag[1]*BmagInv[1]*fr[5]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fr[2]*dfac_v*dfac_x-1.732050807568877*Bmag[1]*BmagInv[1]*fr[1]*dfac_v*dfac_x-1.732050807568877*BmagInv[0]*fr[0]*Bmag[1]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[10] = -(0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[1]*fr[15]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[0]*fr[14]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[1]*fr[13]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[0]*fr[10]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fr[12]*dfac_m*m_*wv+6.0*BcurvY[0]*fr[9]*dfac_m*m_*wv-3.464101615137754*BcurvY[1]*fr[8]*dfac_m*m_*wv-3.464101615137754*BcurvY[0]*fr[4]*dfac_m*m_*wv+9.0*Bmag[1]*BmagInv[1]*fr[15]*dfac_m*dfac_v*dfac_x*wm+9.0*BmagInv[0]*Bmag[1]*fr[14]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*Bmag[1]*BmagInv[1]*fr[13]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*BmagInv[0]*Bmag[1]*fr[10]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[1]*Phi[3]*fr[15]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fr[15]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[0]*Phi[3]*fr[14]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fr[14]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[3]*fr[13]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[1]*fr[13]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[3]*fr[10]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[1]*fr[10]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*Bmag[1]*BmagInv[1]*fr[11]*dfac_v*dfac_x+5.196152422706631*BmagInv[0]*Bmag[1]*fr[7]*dfac_v*dfac_x-3.0*Bmag[1]*BmagInv[1]*fr[6]*dfac_v*dfac_x-3.0*BmagInv[0]*Bmag[1]*fr[3]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[11] = (0.125*dfac_y*(6.0*BcurvY[0]*fr[11]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fr[7]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[0]*fr[6]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[1]*fr[3]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fr[5]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fr[2]*dfac_m*m_*wv-2.0*BcurvY[0]*fr[1]*dfac_m*m_*wv-2.0*fr[0]*BcurvY[1]*dfac_m*m_*wv+5.196152422706631*BmagInv[0]*Bmag[1]*fr[11]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fr[7]*dfac_m*dfac_v*dfac_x*wm-3.0*BmagInv[0]*Bmag[1]*fr[6]*dfac_m*dfac_v*dfac_x*wm-3.0*Bmag[1]*BmagInv[1]*fr[3]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[0]*Phi[3]*fr[11]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fr[11]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fr[7]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fr[7]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[3]*fr[6]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[0]*Phi[1]*fr[6]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[3]*fr[3]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[1]*Phi[1]*fr[3]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Bmag[1]*fr[15]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fr[14]*dfac_v*dfac_x-1.732050807568877*BmagInv[0]*Bmag[1]*fr[13]*dfac_v*dfac_x-1.732050807568877*Bmag[1]*BmagInv[1]*fr[10]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[12] = (0.125*dfac_y*(6.0*BcurvY[0]*fr[12]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fr[9]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[0]*fr[8]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[1]*fr[4]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fr[15]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fr[14]*dfac_m*m_*wv-2.0*BcurvY[0]*fr[13]*dfac_m*m_*wv-2.0*BcurvY[1]*fr[10]*dfac_m*m_*wv+5.196152422706631*BmagInv[0]*Bmag[1]*fr[12]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fr[9]*dfac_m*dfac_v*dfac_x*wm-3.0*BmagInv[0]*Bmag[1]*fr[8]*dfac_m*dfac_v*dfac_x*wm-3.0*Bmag[1]*BmagInv[1]*fr[4]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[0]*Phi[3]*fr[12]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fr[12]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fr[9]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fr[9]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[3]*fr[8]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[0]*Phi[1]*fr[8]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[3]*fr[4]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[1]*Phi[1]*fr[4]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Bmag[1]*fr[5]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fr[2]*dfac_v*dfac_x-1.732050807568877*BmagInv[0]*Bmag[1]*fr[1]*dfac_v*dfac_x-1.732050807568877*fr[0]*Bmag[1]*BmagInv[1]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[13] = -(0.04166666666666666*dfac_y*(10.39230484541326*BcurvY[0]*fr[15]*dfac_m*dfac_v*m_*wv2+10.39230484541326*BcurvY[1]*fr[14]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[0]*fr[13]*dfac_m*dfac_v*m_*wv2-6.0*BcurvY[1]*fr[10]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fr[12]*dfac_m*m_*wv+6.0*BcurvY[1]*fr[9]*dfac_m*m_*wv-3.464101615137754*BcurvY[0]*fr[8]*dfac_m*m_*wv-3.464101615137754*BcurvY[1]*fr[4]*dfac_m*m_*wv+9.0*BmagInv[0]*Bmag[1]*fr[15]*dfac_m*dfac_v*dfac_x*wm+9.0*Bmag[1]*BmagInv[1]*fr[14]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*BmagInv[0]*Bmag[1]*fr[13]*dfac_m*dfac_v*dfac_x*wm-5.196152422706631*Bmag[1]*BmagInv[1]*fr[10]*dfac_m*dfac_v*dfac_x*wm-15.58845726811989*BmagInv[0]*Phi[3]*fr[15]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[1]*fr[15]*dfac_m*dfac_v*dfac_x*q_-15.58845726811989*BmagInv[1]*Phi[3]*fr[14]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[1]*fr[14]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[0]*Phi[3]*fr[13]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[0]*Phi[1]*fr[13]*dfac_m*dfac_v*dfac_x*q_+9.0*BmagInv[1]*Phi[3]*fr[10]*dfac_m*dfac_v*dfac_x*q_-5.196152422706631*BmagInv[1]*Phi[1]*fr[10]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Bmag[1]*fr[11]*dfac_v*dfac_x+5.196152422706631*Bmag[1]*BmagInv[1]*fr[7]*dfac_v*dfac_x-3.0*BmagInv[0]*Bmag[1]*fr[6]*dfac_v*dfac_x-3.0*Bmag[1]*BmagInv[1]*fr[3]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[14] = (0.125*dfac_y*(6.0*BcurvY[1]*fr[15]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[0]*fr[14]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[1]*fr[13]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[0]*fr[10]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[1]*fr[12]*dfac_m*m_*wv+3.464101615137754*BcurvY[0]*fr[9]*dfac_m*m_*wv-2.0*BcurvY[1]*fr[8]*dfac_m*m_*wv-2.0*BcurvY[0]*fr[4]*dfac_m*m_*wv+5.196152422706631*Bmag[1]*BmagInv[1]*fr[15]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*BmagInv[0]*Bmag[1]*fr[14]*dfac_m*dfac_v*dfac_x*wm-3.0*Bmag[1]*BmagInv[1]*fr[13]*dfac_m*dfac_v*dfac_x*wm-3.0*BmagInv[0]*Bmag[1]*fr[10]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[1]*Phi[3]*fr[15]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fr[15]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[0]*Phi[3]*fr[14]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fr[14]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[3]*fr[13]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[1]*Phi[1]*fr[13]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[3]*fr[10]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[0]*Phi[1]*fr[10]*dfac_m*dfac_v*dfac_x*q_+3.0*Bmag[1]*BmagInv[1]*fr[11]*dfac_v*dfac_x+3.0*BmagInv[0]*Bmag[1]*fr[7]*dfac_v*dfac_x-1.732050807568877*Bmag[1]*BmagInv[1]*fr[6]*dfac_v*dfac_x-1.732050807568877*BmagInv[0]*Bmag[1]*fr[3]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 
  incr[15] = (0.125*dfac_y*(6.0*BcurvY[0]*fr[15]*dfac_m*dfac_v*m_*wv2+6.0*BcurvY[1]*fr[14]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[0]*fr[13]*dfac_m*dfac_v*m_*wv2-3.464101615137754*BcurvY[1]*fr[10]*dfac_m*dfac_v*m_*wv2+3.464101615137754*BcurvY[0]*fr[12]*dfac_m*m_*wv+3.464101615137754*BcurvY[1]*fr[9]*dfac_m*m_*wv-2.0*BcurvY[0]*fr[8]*dfac_m*m_*wv-2.0*BcurvY[1]*fr[4]*dfac_m*m_*wv+5.196152422706631*BmagInv[0]*Bmag[1]*fr[15]*dfac_m*dfac_v*dfac_x*wm+5.196152422706631*Bmag[1]*BmagInv[1]*fr[14]*dfac_m*dfac_v*dfac_x*wm-3.0*BmagInv[0]*Bmag[1]*fr[13]*dfac_m*dfac_v*dfac_x*wm-3.0*Bmag[1]*BmagInv[1]*fr[10]*dfac_m*dfac_v*dfac_x*wm-9.0*BmagInv[0]*Phi[3]*fr[15]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[1]*fr[15]*dfac_m*dfac_v*dfac_x*q_-9.0*BmagInv[1]*Phi[3]*fr[14]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[1]*fr[14]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[0]*Phi[3]*fr[13]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[0]*Phi[1]*fr[13]*dfac_m*dfac_v*dfac_x*q_+5.196152422706631*BmagInv[1]*Phi[3]*fr[10]*dfac_m*dfac_v*dfac_x*q_-3.0*BmagInv[1]*Phi[1]*fr[10]*dfac_m*dfac_v*dfac_x*q_+3.0*BmagInv[0]*Bmag[1]*fr[11]*dfac_v*dfac_x+3.0*Bmag[1]*BmagInv[1]*fr[7]*dfac_v*dfac_x-1.732050807568877*BmagInv[0]*Bmag[1]*fr[6]*dfac_v*dfac_x-1.732050807568877*Bmag[1]*BmagInv[1]*fr[3]*dfac_v*dfac_x))/(dfac_m*dfac_v*q_); 

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
double GyrokineticSurf2x2vSer_Vpar_P1(const double q_, const double m_, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
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
  double Ghat[16]; 
  double incr[16]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.5*((1.732050807568877*BcurvY[1]*Phi[3]+1.732050807568877*BcurvY[0]*Phi[2])*dfac_v*dfac_y*wv+((-1.732050807568877*BcurvY[1]*Phi[3])-1.732050807568877*BcurvY[0]*Phi[2])*dfac_y))/dfac_v;

  if (alpha0>0) { 
  incr[0] = -0.125*(3.0*BcurvY[0]*Phi[3]*fl[6]+3.0*BcurvY[1]*Phi[2]*fl[6]+3.0*BcurvY[1]*Phi[3]*fl[3]+3.0*BcurvY[0]*Phi[2]*fl[3]+1.732050807568877*BcurvY[0]*fl[1]*Phi[3]+1.732050807568877*fl[0]*BcurvY[1]*Phi[3]+1.732050807568877*BcurvY[1]*fl[1]*Phi[2]+1.732050807568877*BcurvY[0]*fl[0]*Phi[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[1] = -0.025*(27.0*BcurvY[1]*Phi[3]*fl[6]+15.0*BcurvY[0]*Phi[2]*fl[6]+15.0*BcurvY[0]*Phi[3]*fl[3]+15.0*BcurvY[1]*Phi[2]*fl[3]+15.58845726811989*BcurvY[1]*fl[1]*Phi[3]+8.660254037844386*BcurvY[0]*fl[0]*Phi[3]+8.660254037844386*BcurvY[0]*fl[1]*Phi[2]+8.660254037844386*fl[0]*BcurvY[1]*Phi[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[2] = -0.125*(3.0*BcurvY[0]*Phi[3]*fl[11]+3.0*BcurvY[1]*Phi[2]*fl[11]+3.0*BcurvY[1]*Phi[3]*fl[7]+3.0*BcurvY[0]*Phi[2]*fl[7]+1.732050807568877*BcurvY[0]*Phi[3]*fl[5]+1.732050807568877*BcurvY[1]*Phi[2]*fl[5]+1.732050807568877*BcurvY[1]*fl[2]*Phi[3]+1.732050807568877*BcurvY[0]*Phi[2]*fl[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[3] = 0.375*(1.732050807568877*BcurvY[0]*Phi[3]*fl[6]+1.732050807568877*BcurvY[1]*Phi[2]*fl[6]+1.732050807568877*BcurvY[1]*Phi[3]*fl[3]+1.732050807568877*BcurvY[0]*Phi[2]*fl[3]+BcurvY[0]*fl[1]*Phi[3]+fl[0]*BcurvY[1]*Phi[3]+BcurvY[1]*fl[1]*Phi[2]+BcurvY[0]*fl[0]*Phi[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[4] = -0.125*(3.0*BcurvY[0]*Phi[3]*fl[13]+3.0*BcurvY[1]*Phi[2]*fl[13]+3.0*BcurvY[1]*Phi[3]*fl[10]+3.0*BcurvY[0]*Phi[2]*fl[10]+1.732050807568877*BcurvY[0]*Phi[3]*fl[8]+1.732050807568877*BcurvY[1]*Phi[2]*fl[8]+1.732050807568877*BcurvY[1]*Phi[3]*fl[4]+1.732050807568877*BcurvY[0]*Phi[2]*fl[4])*dfac_y*(dfac_v*wv-1.0); 
  incr[5] = -0.025*(27.0*BcurvY[1]*Phi[3]*fl[11]+15.0*BcurvY[0]*Phi[2]*fl[11]+15.0*BcurvY[0]*Phi[3]*fl[7]+15.0*BcurvY[1]*Phi[2]*fl[7]+15.58845726811989*BcurvY[1]*Phi[3]*fl[5]+8.660254037844386*BcurvY[0]*Phi[2]*fl[5]+8.660254037844386*BcurvY[0]*fl[2]*Phi[3]+8.660254037844386*BcurvY[1]*Phi[2]*fl[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[6] = 0.075*(15.58845726811989*BcurvY[1]*Phi[3]*fl[6]+8.660254037844386*BcurvY[0]*Phi[2]*fl[6]+8.660254037844386*BcurvY[0]*Phi[3]*fl[3]+8.660254037844386*BcurvY[1]*Phi[2]*fl[3]+9.0*BcurvY[1]*fl[1]*Phi[3]+5.0*BcurvY[0]*fl[0]*Phi[3]+5.0*BcurvY[0]*fl[1]*Phi[2]+5.0*fl[0]*BcurvY[1]*Phi[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[7] = 0.375*(1.732050807568877*BcurvY[0]*Phi[3]*fl[11]+1.732050807568877*BcurvY[1]*Phi[2]*fl[11]+1.732050807568877*BcurvY[1]*Phi[3]*fl[7]+1.732050807568877*BcurvY[0]*Phi[2]*fl[7]+BcurvY[0]*Phi[3]*fl[5]+BcurvY[1]*Phi[2]*fl[5]+BcurvY[1]*fl[2]*Phi[3]+BcurvY[0]*Phi[2]*fl[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[8] = -0.025*(27.0*BcurvY[1]*Phi[3]*fl[13]+15.0*BcurvY[0]*Phi[2]*fl[13]+15.0*BcurvY[0]*Phi[3]*fl[10]+15.0*BcurvY[1]*Phi[2]*fl[10]+15.58845726811989*BcurvY[1]*Phi[3]*fl[8]+8.660254037844386*BcurvY[0]*Phi[2]*fl[8]+8.660254037844386*BcurvY[0]*Phi[3]*fl[4]+8.660254037844386*BcurvY[1]*Phi[2]*fl[4])*dfac_y*(dfac_v*wv-1.0); 
  incr[9] = -0.125*(3.0*BcurvY[0]*Phi[3]*fl[15]+3.0*BcurvY[1]*Phi[2]*fl[15]+3.0*BcurvY[1]*Phi[3]*fl[14]+3.0*BcurvY[0]*Phi[2]*fl[14]+1.732050807568877*BcurvY[0]*Phi[3]*fl[12]+1.732050807568877*BcurvY[1]*Phi[2]*fl[12]+1.732050807568877*BcurvY[1]*Phi[3]*fl[9]+1.732050807568877*BcurvY[0]*Phi[2]*fl[9])*dfac_y*(dfac_v*wv-1.0); 
  incr[10] = 0.375*(1.732050807568877*BcurvY[0]*Phi[3]*fl[13]+1.732050807568877*BcurvY[1]*Phi[2]*fl[13]+1.732050807568877*BcurvY[1]*Phi[3]*fl[10]+1.732050807568877*BcurvY[0]*Phi[2]*fl[10]+BcurvY[0]*Phi[3]*fl[8]+BcurvY[1]*Phi[2]*fl[8]+BcurvY[1]*Phi[3]*fl[4]+BcurvY[0]*Phi[2]*fl[4])*dfac_y*(dfac_v*wv-1.0); 
  incr[11] = 0.075*(15.58845726811989*BcurvY[1]*Phi[3]*fl[11]+8.660254037844386*BcurvY[0]*Phi[2]*fl[11]+8.660254037844386*BcurvY[0]*Phi[3]*fl[7]+8.660254037844386*BcurvY[1]*Phi[2]*fl[7]+9.0*BcurvY[1]*Phi[3]*fl[5]+5.0*BcurvY[0]*Phi[2]*fl[5]+5.0*BcurvY[0]*fl[2]*Phi[3]+5.0*BcurvY[1]*Phi[2]*fl[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[12] = -0.025*(27.0*BcurvY[1]*Phi[3]*fl[15]+15.0*BcurvY[0]*Phi[2]*fl[15]+15.0*BcurvY[0]*Phi[3]*fl[14]+15.0*BcurvY[1]*Phi[2]*fl[14]+15.58845726811989*BcurvY[1]*Phi[3]*fl[12]+8.660254037844386*BcurvY[0]*Phi[2]*fl[12]+8.660254037844386*BcurvY[0]*Phi[3]*fl[9]+8.660254037844386*BcurvY[1]*Phi[2]*fl[9])*dfac_y*(dfac_v*wv-1.0); 
  incr[13] = 0.075*(15.58845726811989*BcurvY[1]*Phi[3]*fl[13]+8.660254037844386*BcurvY[0]*Phi[2]*fl[13]+8.660254037844386*BcurvY[0]*Phi[3]*fl[10]+8.660254037844386*BcurvY[1]*Phi[2]*fl[10]+9.0*BcurvY[1]*Phi[3]*fl[8]+5.0*BcurvY[0]*Phi[2]*fl[8]+5.0*BcurvY[0]*Phi[3]*fl[4]+5.0*BcurvY[1]*Phi[2]*fl[4])*dfac_y*(dfac_v*wv-1.0); 
  incr[14] = 0.375*(1.732050807568877*BcurvY[0]*Phi[3]*fl[15]+1.732050807568877*BcurvY[1]*Phi[2]*fl[15]+1.732050807568877*BcurvY[1]*Phi[3]*fl[14]+1.732050807568877*BcurvY[0]*Phi[2]*fl[14]+BcurvY[0]*Phi[3]*fl[12]+BcurvY[1]*Phi[2]*fl[12]+BcurvY[1]*Phi[3]*fl[9]+BcurvY[0]*Phi[2]*fl[9])*dfac_y*(dfac_v*wv-1.0); 
  incr[15] = 0.075*(15.58845726811989*BcurvY[1]*Phi[3]*fl[15]+8.660254037844386*BcurvY[0]*Phi[2]*fl[15]+8.660254037844386*BcurvY[0]*Phi[3]*fl[14]+8.660254037844386*BcurvY[1]*Phi[2]*fl[14]+9.0*BcurvY[1]*Phi[3]*fl[12]+5.0*BcurvY[0]*Phi[2]*fl[12]+5.0*BcurvY[0]*Phi[3]*fl[9]+5.0*BcurvY[1]*Phi[2]*fl[9])*dfac_y*(dfac_v*wv-1.0); 

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
  incr[0] = 0.125*(3.0*BcurvY[0]*Phi[3]*fr[6]+3.0*BcurvY[1]*Phi[2]*fr[6]+3.0*BcurvY[1]*Phi[3]*fr[3]+3.0*BcurvY[0]*Phi[2]*fr[3]-1.732050807568877*BcurvY[0]*fr[1]*Phi[3]-1.732050807568877*fr[0]*BcurvY[1]*Phi[3]-1.732050807568877*BcurvY[1]*fr[1]*Phi[2]-1.732050807568877*BcurvY[0]*fr[0]*Phi[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[1] = 0.025*(27.0*BcurvY[1]*Phi[3]*fr[6]+15.0*BcurvY[0]*Phi[2]*fr[6]+15.0*BcurvY[0]*Phi[3]*fr[3]+15.0*BcurvY[1]*Phi[2]*fr[3]-15.58845726811989*BcurvY[1]*fr[1]*Phi[3]-8.660254037844386*BcurvY[0]*fr[0]*Phi[3]-8.660254037844386*BcurvY[0]*fr[1]*Phi[2]-8.660254037844386*fr[0]*BcurvY[1]*Phi[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[2] = 0.125*(3.0*BcurvY[0]*Phi[3]*fr[11]+3.0*BcurvY[1]*Phi[2]*fr[11]+3.0*BcurvY[1]*Phi[3]*fr[7]+3.0*BcurvY[0]*Phi[2]*fr[7]-1.732050807568877*BcurvY[0]*Phi[3]*fr[5]-1.732050807568877*BcurvY[1]*Phi[2]*fr[5]-1.732050807568877*BcurvY[1]*fr[2]*Phi[3]-1.732050807568877*BcurvY[0]*Phi[2]*fr[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[3] = -0.375*(1.732050807568877*BcurvY[0]*Phi[3]*fr[6]+1.732050807568877*BcurvY[1]*Phi[2]*fr[6]+1.732050807568877*BcurvY[1]*Phi[3]*fr[3]+1.732050807568877*BcurvY[0]*Phi[2]*fr[3]-1.0*BcurvY[0]*fr[1]*Phi[3]-1.0*fr[0]*BcurvY[1]*Phi[3]-1.0*BcurvY[1]*fr[1]*Phi[2]-1.0*BcurvY[0]*fr[0]*Phi[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[4] = 0.125*(3.0*BcurvY[0]*Phi[3]*fr[13]+3.0*BcurvY[1]*Phi[2]*fr[13]+3.0*BcurvY[1]*Phi[3]*fr[10]+3.0*BcurvY[0]*Phi[2]*fr[10]-1.732050807568877*BcurvY[0]*Phi[3]*fr[8]-1.732050807568877*BcurvY[1]*Phi[2]*fr[8]-1.732050807568877*BcurvY[1]*Phi[3]*fr[4]-1.732050807568877*BcurvY[0]*Phi[2]*fr[4])*dfac_y*(dfac_v*wv-1.0); 
  incr[5] = 0.025*(27.0*BcurvY[1]*Phi[3]*fr[11]+15.0*BcurvY[0]*Phi[2]*fr[11]+15.0*BcurvY[0]*Phi[3]*fr[7]+15.0*BcurvY[1]*Phi[2]*fr[7]-15.58845726811989*BcurvY[1]*Phi[3]*fr[5]-8.660254037844386*BcurvY[0]*Phi[2]*fr[5]-8.660254037844386*BcurvY[0]*fr[2]*Phi[3]-8.660254037844386*BcurvY[1]*Phi[2]*fr[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[6] = -0.075*(15.58845726811989*BcurvY[1]*Phi[3]*fr[6]+8.660254037844386*BcurvY[0]*Phi[2]*fr[6]+8.660254037844386*BcurvY[0]*Phi[3]*fr[3]+8.660254037844386*BcurvY[1]*Phi[2]*fr[3]-9.0*BcurvY[1]*fr[1]*Phi[3]-5.0*BcurvY[0]*fr[0]*Phi[3]-5.0*BcurvY[0]*fr[1]*Phi[2]-5.0*fr[0]*BcurvY[1]*Phi[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[7] = -0.375*(1.732050807568877*BcurvY[0]*Phi[3]*fr[11]+1.732050807568877*BcurvY[1]*Phi[2]*fr[11]+1.732050807568877*BcurvY[1]*Phi[3]*fr[7]+1.732050807568877*BcurvY[0]*Phi[2]*fr[7]-1.0*BcurvY[0]*Phi[3]*fr[5]-1.0*BcurvY[1]*Phi[2]*fr[5]-1.0*BcurvY[1]*fr[2]*Phi[3]-1.0*BcurvY[0]*Phi[2]*fr[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[8] = 0.025*(27.0*BcurvY[1]*Phi[3]*fr[13]+15.0*BcurvY[0]*Phi[2]*fr[13]+15.0*BcurvY[0]*Phi[3]*fr[10]+15.0*BcurvY[1]*Phi[2]*fr[10]-15.58845726811989*BcurvY[1]*Phi[3]*fr[8]-8.660254037844386*BcurvY[0]*Phi[2]*fr[8]-8.660254037844386*BcurvY[0]*Phi[3]*fr[4]-8.660254037844386*BcurvY[1]*Phi[2]*fr[4])*dfac_y*(dfac_v*wv-1.0); 
  incr[9] = 0.125*(3.0*BcurvY[0]*Phi[3]*fr[15]+3.0*BcurvY[1]*Phi[2]*fr[15]+3.0*BcurvY[1]*Phi[3]*fr[14]+3.0*BcurvY[0]*Phi[2]*fr[14]-1.732050807568877*BcurvY[0]*Phi[3]*fr[12]-1.732050807568877*BcurvY[1]*Phi[2]*fr[12]-1.732050807568877*BcurvY[1]*Phi[3]*fr[9]-1.732050807568877*BcurvY[0]*Phi[2]*fr[9])*dfac_y*(dfac_v*wv-1.0); 
  incr[10] = -0.375*(1.732050807568877*BcurvY[0]*Phi[3]*fr[13]+1.732050807568877*BcurvY[1]*Phi[2]*fr[13]+1.732050807568877*BcurvY[1]*Phi[3]*fr[10]+1.732050807568877*BcurvY[0]*Phi[2]*fr[10]-1.0*BcurvY[0]*Phi[3]*fr[8]-1.0*BcurvY[1]*Phi[2]*fr[8]-1.0*BcurvY[1]*Phi[3]*fr[4]-1.0*BcurvY[0]*Phi[2]*fr[4])*dfac_y*(dfac_v*wv-1.0); 
  incr[11] = -0.075*(15.58845726811989*BcurvY[1]*Phi[3]*fr[11]+8.660254037844386*BcurvY[0]*Phi[2]*fr[11]+8.660254037844386*BcurvY[0]*Phi[3]*fr[7]+8.660254037844386*BcurvY[1]*Phi[2]*fr[7]-9.0*BcurvY[1]*Phi[3]*fr[5]-5.0*BcurvY[0]*Phi[2]*fr[5]-5.0*BcurvY[0]*fr[2]*Phi[3]-5.0*BcurvY[1]*Phi[2]*fr[2])*dfac_y*(dfac_v*wv-1.0); 
  incr[12] = 0.025*(27.0*BcurvY[1]*Phi[3]*fr[15]+15.0*BcurvY[0]*Phi[2]*fr[15]+15.0*BcurvY[0]*Phi[3]*fr[14]+15.0*BcurvY[1]*Phi[2]*fr[14]-15.58845726811989*BcurvY[1]*Phi[3]*fr[12]-8.660254037844386*BcurvY[0]*Phi[2]*fr[12]-8.660254037844386*BcurvY[0]*Phi[3]*fr[9]-8.660254037844386*BcurvY[1]*Phi[2]*fr[9])*dfac_y*(dfac_v*wv-1.0); 
  incr[13] = -0.075*(15.58845726811989*BcurvY[1]*Phi[3]*fr[13]+8.660254037844386*BcurvY[0]*Phi[2]*fr[13]+8.660254037844386*BcurvY[0]*Phi[3]*fr[10]+8.660254037844386*BcurvY[1]*Phi[2]*fr[10]-9.0*BcurvY[1]*Phi[3]*fr[8]-5.0*BcurvY[0]*Phi[2]*fr[8]-5.0*BcurvY[0]*Phi[3]*fr[4]-5.0*BcurvY[1]*Phi[2]*fr[4])*dfac_y*(dfac_v*wv-1.0); 
  incr[14] = -0.375*(1.732050807568877*BcurvY[0]*Phi[3]*fr[15]+1.732050807568877*BcurvY[1]*Phi[2]*fr[15]+1.732050807568877*BcurvY[1]*Phi[3]*fr[14]+1.732050807568877*BcurvY[0]*Phi[2]*fr[14]-1.0*BcurvY[0]*Phi[3]*fr[12]-1.0*BcurvY[1]*Phi[2]*fr[12]-1.0*BcurvY[1]*Phi[3]*fr[9]-1.0*BcurvY[0]*Phi[2]*fr[9])*dfac_y*(dfac_v*wv-1.0); 
  incr[15] = -0.075*(15.58845726811989*BcurvY[1]*Phi[3]*fr[15]+8.660254037844386*BcurvY[0]*Phi[2]*fr[15]+8.660254037844386*BcurvY[0]*Phi[3]*fr[14]+8.660254037844386*BcurvY[1]*Phi[2]*fr[14]-9.0*BcurvY[1]*Phi[3]*fr[12]-5.0*BcurvY[0]*Phi[2]*fr[12]-5.0*BcurvY[0]*Phi[3]*fr[9]-5.0*BcurvY[1]*Phi[2]*fr[9])*dfac_y*(dfac_v*wv-1.0); 

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
