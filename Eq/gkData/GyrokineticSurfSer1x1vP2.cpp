#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x1vSer_X_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.5*wv; 

  double alpha[8]; 
  alpha[0] = 2.0*wv; 
  alpha[2] = 1.154700538379252/dfac_v; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  double favg[8]; 

  favg[0] = 1.0*fr[0]+fl[0]; 
  favg[1] = -1.0*fr[1]+fl[1]; 
  favg[2] = 1.0*fr[2]+fl[2]; 
  favg[3] = -1.0*fr[3]+fl[3]; 
  favg[4] = 1.0*fr[4]+fl[4]; 
  favg[5] = 1.0*fr[5]+fl[5]; 
  favg[6] = 1.0*fr[6]+fl[6]; 
  favg[7] = -1.0*fr[7]+fl[7]; 
  double fjump[8]; 

  fjump[0] = amax*(1.0*fr[0]-fl[0]); 
  fjump[1] = amax*(-1.0*fr[1]-fl[1]); 
  fjump[2] = amax*(1.0*fr[2]-fl[2]); 
  fjump[3] = amax*(-1.0*fr[3]-fl[3]); 
  fjump[4] = amax*(1.0*fr[4]-fl[4]); 
  fjump[5] = amax*(1.0*fr[5]-fl[5]); 
  fjump[6] = amax*(1.0*fr[6]-fl[6]); 
  fjump[7] = amax*(-1.0*fr[7]-fl[7]); 
  Ghat[0] = alpha[2]*(0.5590169943749476*favg[6]+0.4330127018922193*favg[3]+0.25*favg[2])-1.118033988749895*fjump[4]+alpha[0]*(0.5590169943749475*favg[4]+0.4330127018922193*favg[1]+0.25*favg[0])-0.8660254037844386*fjump[1]-0.5*fjump[0]; 
  Ghat[2] = alpha[2]*(0.3872983346207417*favg[7]+0.223606797749979*favg[5]+0.5590169943749475*favg[4]+0.4330127018922193*favg[1]+0.25*favg[0])-1.118033988749895*fjump[6]+alpha[0]*(0.5590169943749476*favg[6]+0.4330127018922193*favg[3]+0.25*favg[2])-0.8660254037844386*fjump[3]-0.5*fjump[2]; 
  Ghat[5] = (-0.8660254037844387*fjump[7])+alpha[0]*(0.4330127018922194*favg[7]+0.25*favg[5])+alpha[2]*(0.5000000000000001*favg[6]+0.3872983346207416*favg[3]+0.223606797749979*favg[2])-0.5*fjump[5]; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[4] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[5] = 0.5*Ghat[5]*dfac_x; 
  incr[6] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[7] = -0.8660254037844387*Ghat[5]*dfac_x; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.6123724356957944*Phi[1]*dfac_x*q_)/m_; 

  double alpha[8]; 
  alpha[0] = -(2.449489742783178*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(5.477225575051662*Phi[2]*dfac_x*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  double favg[8]; 

  favg[0] = 1.0*fr[0]+fl[0]; 
  favg[1] = 1.0*fr[1]+fl[1]; 
  favg[2] = -1.0*fr[2]+fl[2]; 
  favg[3] = -1.0*fr[3]+fl[3]; 
  favg[4] = 1.0*fr[4]+fl[4]; 
  favg[5] = 1.0*fr[5]+fl[5]; 
  favg[6] = -1.0*fr[6]+fl[6]; 
  favg[7] = 1.0*fr[7]+fl[7]; 
  double fjump[8]; 

  fjump[0] = amax*(1.0*fr[0]-fl[0]); 
  fjump[1] = amax*(1.0*fr[1]-fl[1]); 
  fjump[2] = amax*(-1.0*fr[2]-fl[2]); 
  fjump[3] = amax*(-1.0*fr[3]-fl[3]); 
  fjump[4] = amax*(1.0*fr[4]-fl[4]); 
  fjump[5] = amax*(1.0*fr[5]-fl[5]); 
  fjump[6] = amax*(-1.0*fr[6]-fl[6]); 
  fjump[7] = amax*(1.0*fr[7]-fl[7]); 
  Ghat[0] = alpha[1]*(0.5590169943749476*favg[7]+0.4330127018922193*favg[3]+0.25*favg[1])-1.118033988749895*fjump[5]+alpha[0]*(0.5590169943749475*favg[5]+0.4330127018922193*favg[2]+0.25*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (-1.118033988749895*fjump[7])+alpha[0]*(0.5590169943749476*favg[7]+0.4330127018922193*favg[3]+0.25*favg[1])+alpha[1]*(0.3872983346207417*favg[6]+0.5590169943749475*favg[5]+0.223606797749979*favg[4]+0.4330127018922193*favg[2]+0.25*favg[0])-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = alpha[1]*(0.5000000000000001*favg[7]+0.3872983346207416*favg[3]+0.223606797749979*favg[1])-0.8660254037844387*fjump[6]+alpha[0]*(0.4330127018922194*favg[6]+0.25*favg[4])-0.5*fjump[4]; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[4] = 0.5*Ghat[4]*dfac_v; 
  incr[5] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[6] = -0.8660254037844387*Ghat[4]*dfac_v; 
  incr[7] = 1.118033988749895*Ghat[1]*dfac_v; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_X_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.5*wv; 

  double alpha[8]; 
  alpha[0] = 2.0*wv; 
  alpha[2] = 1.154700538379252/dfac_v; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  double favg[8]; 

  favg[0] = 1.0*fr[0]+fl[0]; 
  favg[1] = -1.0*fr[1]+fl[1]; 
  favg[2] = 1.0*fr[2]+fl[2]; 
  favg[3] = -1.0*fr[3]+fl[3]; 
  favg[4] = 1.0*fr[4]+fl[4]; 
  favg[5] = 1.0*fr[5]+fl[5]; 
  favg[6] = 1.0*fr[6]+fl[6]; 
  favg[7] = -1.0*fr[7]+fl[7]; 
  double fjump[8]; 

  fjump[0] = amax*(1.0*fr[0]-fl[0]); 
  fjump[1] = amax*(-1.0*fr[1]-fl[1]); 
  fjump[2] = amax*(1.0*fr[2]-fl[2]); 
  fjump[3] = amax*(-1.0*fr[3]-fl[3]); 
  fjump[4] = amax*(1.0*fr[4]-fl[4]); 
  fjump[5] = amax*(1.0*fr[5]-fl[5]); 
  fjump[6] = amax*(1.0*fr[6]-fl[6]); 
  fjump[7] = amax*(-1.0*fr[7]-fl[7]); 
  Ghat[0] = alpha[2]*(0.5590169943749476*favg[6]+0.4330127018922193*favg[3]+0.25*favg[2])-1.118033988749895*fjump[4]+alpha[0]*(0.5590169943749475*favg[4]+0.4330127018922193*favg[1]+0.25*favg[0])-0.8660254037844386*fjump[1]-0.5*fjump[0]; 
  Ghat[2] = alpha[2]*(0.3872983346207417*favg[7]+0.223606797749979*favg[5]+0.5590169943749475*favg[4]+0.4330127018922193*favg[1]+0.25*favg[0])-1.118033988749895*fjump[6]+alpha[0]*(0.5590169943749476*favg[6]+0.4330127018922193*favg[3]+0.25*favg[2])-0.8660254037844386*fjump[3]-0.5*fjump[2]; 
  Ghat[5] = (-0.8660254037844387*fjump[7])+alpha[0]*(0.4330127018922194*favg[7]+0.25*favg[5])+alpha[2]*(0.5000000000000001*favg[6]+0.3872983346207416*favg[3]+0.223606797749979*favg[2])-0.5*fjump[5]; 

  incr[0] = 0.5*Ghat[0]*dfac_x; 
  incr[1] = -0.8660254037844386*Ghat[0]*dfac_x; 
  incr[2] = 0.5*Ghat[2]*dfac_x; 
  incr[3] = -0.8660254037844386*Ghat[2]*dfac_x; 
  incr[4] = 1.118033988749895*Ghat[0]*dfac_x; 
  incr[5] = 0.5*Ghat[5]*dfac_x; 
  incr[6] = 1.118033988749895*Ghat[2]*dfac_x; 
  incr[7] = -0.8660254037844387*Ghat[5]*dfac_x; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += -1.0*incr[6]; 
  outl[7] += incr[7]; 
return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double cfl, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.6123724356957944*Phi[1]*dfac_x*q_)/m_; 

  double alpha[8]; 
  alpha[0] = -(2.449489742783178*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(5.477225575051662*Phi[2]*dfac_x*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  double Ghat[8]; 
  double favg[8]; 

  favg[0] = 1.0*fr[0]+fl[0]; 
  favg[1] = 1.0*fr[1]+fl[1]; 
  favg[2] = -1.0*fr[2]+fl[2]; 
  favg[3] = -1.0*fr[3]+fl[3]; 
  favg[4] = 1.0*fr[4]+fl[4]; 
  favg[5] = 1.0*fr[5]+fl[5]; 
  favg[6] = -1.0*fr[6]+fl[6]; 
  favg[7] = 1.0*fr[7]+fl[7]; 
  double fjump[8]; 

  fjump[0] = amax*(1.0*fr[0]-fl[0]); 
  fjump[1] = amax*(1.0*fr[1]-fl[1]); 
  fjump[2] = amax*(-1.0*fr[2]-fl[2]); 
  fjump[3] = amax*(-1.0*fr[3]-fl[3]); 
  fjump[4] = amax*(1.0*fr[4]-fl[4]); 
  fjump[5] = amax*(1.0*fr[5]-fl[5]); 
  fjump[6] = amax*(-1.0*fr[6]-fl[6]); 
  fjump[7] = amax*(1.0*fr[7]-fl[7]); 
  Ghat[0] = alpha[1]*(0.5590169943749476*favg[7]+0.4330127018922193*favg[3]+0.25*favg[1])-1.118033988749895*fjump[5]+alpha[0]*(0.5590169943749475*favg[5]+0.4330127018922193*favg[2]+0.25*favg[0])-0.8660254037844386*fjump[2]-0.5*fjump[0]; 
  Ghat[1] = (-1.118033988749895*fjump[7])+alpha[0]*(0.5590169943749476*favg[7]+0.4330127018922193*favg[3]+0.25*favg[1])+alpha[1]*(0.3872983346207417*favg[6]+0.5590169943749475*favg[5]+0.223606797749979*favg[4]+0.4330127018922193*favg[2]+0.25*favg[0])-0.8660254037844386*fjump[3]-0.5*fjump[1]; 
  Ghat[4] = alpha[1]*(0.5000000000000001*favg[7]+0.3872983346207416*favg[3]+0.223606797749979*favg[1])-0.8660254037844387*fjump[6]+alpha[0]*(0.4330127018922194*favg[6]+0.25*favg[4])-0.5*fjump[4]; 

  incr[0] = 0.5*Ghat[0]*dfac_v; 
  incr[1] = 0.5*Ghat[1]*dfac_v; 
  incr[2] = -0.8660254037844386*Ghat[0]*dfac_v; 
  incr[3] = -0.8660254037844386*Ghat[1]*dfac_v; 
  incr[4] = 0.5*Ghat[4]*dfac_v; 
  incr[5] = 1.118033988749895*Ghat[0]*dfac_v; 
  incr[6] = -0.8660254037844387*Ghat[4]*dfac_v; 
  incr[7] = 1.118033988749895*Ghat[1]*dfac_v; 

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
  outl[3] += incr[3]; 
  outl[4] += -1.0*incr[4]; 
  outl[5] += -1.0*incr[5]; 
  outl[6] += incr[6]; 
  outl[7] += -1.0*incr[7]; 
return std::abs(alpha0); 
} 
