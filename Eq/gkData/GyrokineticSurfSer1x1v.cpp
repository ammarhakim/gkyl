#include <GyrokineticModDecl.h> 
double GyrokineticSurf1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[4]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.7071067811865475*Gradpar[0]*wv; 

  double alpha[4]; 
  alpha[0] = 1.414213562373095*Gradpar[0]*wv; 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fl[3]+fl[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(3.0*fl[3]+1.732050807568877*fl[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fr[3]-1.0*fr[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(3.0*fr[3]-1.732050807568877*fr[2])*dfac_x; 

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
double GyrokineticSurf1x1vSer_X_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.7071067811865475*Gradpar[0]*wv; 

  double alpha[8]; 
  alpha[0] = 1.414213562373095*Gradpar[0]*wv; 
  alpha[2] = (0.8164965809277261*Gradpar[0])/dfac_v; 
  if (alpha0>0) { 
  incr[0] = 0.08333333333333333*(6.708203932499369*alpha[2]*fl[6]+6.708203932499369*alpha[0]*fl[4]+5.196152422706631*alpha[2]*fl[3]+3.0*alpha[2]*fl[2]+5.196152422706631*alpha[0]*fl[1]+3.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.25*(3.872983346207417*alpha[2]*fl[6]+3.872983346207417*alpha[0]*fl[4]+3.0*alpha[2]*fl[3]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01666666666666667*(23.2379000772445*alpha[2]*fl[7]+33.54101966249684*alpha[0]*fl[6]+13.41640786499874*alpha[2]*fl[5]+33.54101966249685*alpha[2]*fl[4]+25.98076211353316*alpha[0]*fl[3]+15.0*alpha[0]*fl[2]+25.98076211353316*fl[1]*alpha[2]+15.0*fl[0]*alpha[2])*dfac_x; 
  incr[3] = -0.05*(13.41640786499874*alpha[2]*fl[7]+19.36491673103708*alpha[0]*fl[6]+7.745966692414834*alpha[2]*fl[5]+19.36491673103709*alpha[2]*fl[4]+15.0*alpha[0]*fl[3]+8.660254037844386*alpha[0]*fl[2]+15.0*fl[1]*alpha[2]+8.660254037844386*fl[0]*alpha[2])*dfac_x; 
  incr[4] = 0.08333333333333333*(15.0*alpha[2]*fl[6]+15.0*alpha[0]*fl[4]+11.61895003862225*alpha[2]*fl[3]+6.708203932499369*alpha[2]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[5] = 0.01666666666666667*(25.98076211353316*alpha[0]*fl[7]+30.0*alpha[2]*fl[6]+15.0*alpha[0]*fl[5]+23.2379000772445*alpha[2]*fl[3]+13.41640786499874*alpha[2]*fl[2])*dfac_x; 
  incr[6] = 0.01666666666666667*(51.96152422706631*alpha[2]*fl[7]+75.0*alpha[0]*fl[6]+30.0*alpha[2]*fl[5]+75.00000000000001*alpha[2]*fl[4]+58.09475019311126*alpha[0]*fl[3]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*fl[1]*alpha[2]+33.54101966249684*fl[0]*alpha[2])*dfac_x; 
  incr[7] = -0.05*(15.0*alpha[0]*fl[7]+17.32050807568877*alpha[2]*fl[6]+8.660254037844387*alpha[0]*fl[5]+13.41640786499874*alpha[2]*fl[3]+7.745966692414834*alpha[2]*fl[2])*dfac_x; 

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
  } else { 
  incr[0] = 0.08333333333333333*(6.708203932499369*alpha[2]*fr[6]+6.708203932499369*alpha[0]*fr[4]-5.196152422706631*alpha[2]*fr[3]+3.0*alpha[2]*fr[2]-5.196152422706631*alpha[0]*fr[1]+3.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = -0.25*(3.872983346207417*alpha[2]*fr[6]+3.872983346207417*alpha[0]*fr[4]-3.0*alpha[2]*fr[3]+1.732050807568877*alpha[2]*fr[2]-3.0*alpha[0]*fr[1]+1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01666666666666667*(23.2379000772445*alpha[2]*fr[7]-33.54101966249684*alpha[0]*fr[6]-13.41640786499874*alpha[2]*fr[5]-33.54101966249685*alpha[2]*fr[4]+25.98076211353316*alpha[0]*fr[3]-15.0*alpha[0]*fr[2]+25.98076211353316*fr[1]*alpha[2]-15.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05*(13.41640786499874*alpha[2]*fr[7]-19.36491673103708*alpha[0]*fr[6]-7.745966692414834*alpha[2]*fr[5]-19.36491673103709*alpha[2]*fr[4]+15.0*alpha[0]*fr[3]-8.660254037844386*alpha[0]*fr[2]+15.0*fr[1]*alpha[2]-8.660254037844386*fr[0]*alpha[2])*dfac_x; 
  incr[4] = 0.08333333333333333*(15.0*alpha[2]*fr[6]+15.0*alpha[0]*fr[4]-11.61895003862225*alpha[2]*fr[3]+6.708203932499369*alpha[2]*fr[2]-11.61895003862225*alpha[0]*fr[1]+6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[5] = -0.01666666666666667*(25.98076211353316*alpha[0]*fr[7]-30.0*alpha[2]*fr[6]-15.0*alpha[0]*fr[5]+23.2379000772445*alpha[2]*fr[3]-13.41640786499874*alpha[2]*fr[2])*dfac_x; 
  incr[6] = -0.01666666666666667*(51.96152422706631*alpha[2]*fr[7]-75.0*alpha[0]*fr[6]-30.0*alpha[2]*fr[5]-75.00000000000001*alpha[2]*fr[4]+58.09475019311126*alpha[0]*fr[3]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*fr[1]*alpha[2]-33.54101966249684*fr[0]*alpha[2])*dfac_x; 
  incr[7] = 0.05*(15.0*alpha[0]*fr[7]-17.32050807568877*alpha[2]*fr[6]-8.660254037844387*alpha[0]*fr[5]+13.41640786499874*alpha[2]*fr[3]-7.745966692414834*alpha[2]*fr[2])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[4]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.8660254037844386*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[4]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[2]+fl[0])*dfac_v; 
  incr[1] = 0.25*alpha[0]*(1.732050807568877*fl[3]+fl[1])*dfac_v; 
  incr[2] = -0.25*alpha[0]*(3.0*fl[2]+1.732050807568877*fl[0])*dfac_v; 
  incr[3] = -0.25*alpha[0]*(3.0*fl[3]+1.732050807568877*fl[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[2]-1.0*fr[0])*dfac_v; 
  incr[1] = -0.25*alpha[0]*(1.732050807568877*fr[3]-1.0*fr[1])*dfac_v; 
  incr[2] = 0.25*alpha[0]*(3.0*fr[2]-1.732050807568877*fr[0])*dfac_v; 
  incr[3] = 0.25*alpha[0]*(3.0*fr[3]-1.732050807568877*fr[1])*dfac_v; 

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
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.8660254037844386*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 

  double alpha[8]; 
  alpha[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alpha[1] = -(3.872983346207417*Gradpar[0]*Phi[2]*dfac_x*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-1.118033988749895*fr[5]*amax)+1.118033988749895*fl[5]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax+0.5590169943749476*alpha[1]*fr[7]+0.5590169943749476*alpha[1]*fl[7]+0.5590169943749475*alpha[0]*fr[5]+0.5590169943749475*alpha[0]*fl[5]-0.4330127018922193*alpha[1]*fr[3]+0.4330127018922193*alpha[1]*fl[3]-0.4330127018922193*alpha[0]*fr[2]+0.4330127018922193*alpha[0]*fl[2]+0.25*alpha[1]*fr[1]+0.25*alpha[1]*fl[1]+0.25*alpha[0]*fr[0]+0.25*alpha[0]*fl[0]; 
  Ghat[1] = (-1.118033988749895*fr[7]*amax)+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax+0.5590169943749476*alpha[0]*fr[7]+0.5590169943749476*alpha[0]*fl[7]-0.3872983346207417*alpha[1]*fr[6]+0.3872983346207417*alpha[1]*fl[6]+0.5590169943749475*alpha[1]*fr[5]+0.5590169943749475*alpha[1]*fl[5]+0.223606797749979*alpha[1]*fr[4]+0.223606797749979*alpha[1]*fl[4]-0.4330127018922193*alpha[0]*fr[3]+0.4330127018922193*alpha[0]*fl[3]-0.4330127018922193*alpha[1]*fr[2]+0.4330127018922193*alpha[1]*fl[2]+0.25*alpha[0]*fr[1]+0.25*alpha[0]*fl[1]+0.25*fr[0]*alpha[1]+0.25*fl[0]*alpha[1]; 
  Ghat[4] = 0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax+0.5000000000000001*alpha[1]*fr[7]+0.5000000000000001*alpha[1]*fl[7]-0.4330127018922193*alpha[0]*fr[6]+0.4330127018922193*alpha[0]*fl[6]+0.25*alpha[0]*fr[4]+0.25*alpha[0]*fl[4]-0.3872983346207416*alpha[1]*fr[3]+0.3872983346207416*alpha[1]*fl[3]+0.223606797749979*alpha[1]*fr[1]+0.223606797749979*alpha[1]*fl[1]; 

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
double EmGyrokineticSurf1x1vSer_X_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[4]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = 0.3535533905932737*(2.449489742783178*geoY[0]*Apar[1]*dfac_x+2.0*Gradpar[0])*wv; 

  double alpha[4]; 
  alpha[0] = 1.732050807568877*geoY[0]*Apar[1]*dfac_x*wv+1.414213562373095*Gradpar[0]*wv; 
  if (alpha0>0) { 
  incr[0] = 0.25*alpha[0]*(1.732050807568877*fl[1]+fl[0])*dfac_x; 
  incr[1] = -0.25*alpha[0]*(3.0*fl[1]+1.732050807568877*fl[0])*dfac_x; 
  incr[2] = 0.25*alpha[0]*(1.732050807568877*fl[3]+fl[2])*dfac_x; 
  incr[3] = -0.25*alpha[0]*(3.0*fl[3]+1.732050807568877*fl[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.25*alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])*dfac_x; 
  incr[1] = 0.25*alpha[0]*(3.0*fr[1]-1.732050807568877*fr[0])*dfac_x; 
  incr[2] = -0.25*alpha[0]*(1.732050807568877*fr[3]-1.0*fr[2])*dfac_x; 
  incr[3] = 0.25*alpha[0]*(3.0*fr[3]-1.732050807568877*fr[2])*dfac_x; 

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
double EmGyrokineticSurf1x1vSer_X_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -0.3535533905932737*((9.48683298050514*geoY[0]*Apar[2]-2.449489742783178*geoY[0]*Apar[1])*dfac_x-2.0*Gradpar[0])*wv; 

  double alpha[8]; 
  alpha[0] = (-6.708203932499369*geoY[0]*Apar[2]*dfac_x*wv)+1.732050807568877*geoY[0]*Apar[1]*dfac_x*wv+1.414213562373095*Gradpar[0]*wv; 
  alpha[2] = (-(3.872983346207417*geoY[0]*Apar[2]*dfac_x)/dfac_v)+(1.0*geoY[0]*Apar[1]*dfac_x)/dfac_v+(0.8164965809277259*Gradpar[0])/dfac_v; 
  if (alpha0>0) { 
  incr[0] = 0.08333333333333333*(6.708203932499369*alpha[2]*fl[6]+6.708203932499369*alpha[0]*fl[4]+5.196152422706631*alpha[2]*fl[3]+3.0*alpha[2]*fl[2]+5.196152422706631*alpha[0]*fl[1]+3.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.25*(3.872983346207417*alpha[2]*fl[6]+3.872983346207417*alpha[0]*fl[4]+3.0*alpha[2]*fl[3]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01666666666666667*(23.2379000772445*alpha[2]*fl[7]+33.54101966249684*alpha[0]*fl[6]+13.41640786499874*alpha[2]*fl[5]+33.54101966249685*alpha[2]*fl[4]+25.98076211353316*alpha[0]*fl[3]+15.0*alpha[0]*fl[2]+25.98076211353316*fl[1]*alpha[2]+15.0*fl[0]*alpha[2])*dfac_x; 
  incr[3] = -0.05*(13.41640786499874*alpha[2]*fl[7]+19.36491673103708*alpha[0]*fl[6]+7.745966692414834*alpha[2]*fl[5]+19.36491673103709*alpha[2]*fl[4]+15.0*alpha[0]*fl[3]+8.660254037844386*alpha[0]*fl[2]+15.0*fl[1]*alpha[2]+8.660254037844386*fl[0]*alpha[2])*dfac_x; 
  incr[4] = 0.08333333333333333*(15.0*alpha[2]*fl[6]+15.0*alpha[0]*fl[4]+11.61895003862225*alpha[2]*fl[3]+6.708203932499369*alpha[2]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[5] = 0.01666666666666667*(25.98076211353316*alpha[0]*fl[7]+30.0*alpha[2]*fl[6]+15.0*alpha[0]*fl[5]+23.2379000772445*alpha[2]*fl[3]+13.41640786499874*alpha[2]*fl[2])*dfac_x; 
  incr[6] = 0.01666666666666667*(51.96152422706631*alpha[2]*fl[7]+75.0*alpha[0]*fl[6]+30.0*alpha[2]*fl[5]+75.00000000000001*alpha[2]*fl[4]+58.09475019311126*alpha[0]*fl[3]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*fl[1]*alpha[2]+33.54101966249684*fl[0]*alpha[2])*dfac_x; 
  incr[7] = -0.05*(15.0*alpha[0]*fl[7]+17.32050807568877*alpha[2]*fl[6]+8.660254037844387*alpha[0]*fl[5]+13.41640786499874*alpha[2]*fl[3]+7.745966692414834*alpha[2]*fl[2])*dfac_x; 

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
  } else { 
  incr[0] = 0.08333333333333333*(6.708203932499369*alpha[2]*fr[6]+6.708203932499369*alpha[0]*fr[4]-5.196152422706631*alpha[2]*fr[3]+3.0*alpha[2]*fr[2]-5.196152422706631*alpha[0]*fr[1]+3.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = -0.25*(3.872983346207417*alpha[2]*fr[6]+3.872983346207417*alpha[0]*fr[4]-3.0*alpha[2]*fr[3]+1.732050807568877*alpha[2]*fr[2]-3.0*alpha[0]*fr[1]+1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01666666666666667*(23.2379000772445*alpha[2]*fr[7]-33.54101966249684*alpha[0]*fr[6]-13.41640786499874*alpha[2]*fr[5]-33.54101966249685*alpha[2]*fr[4]+25.98076211353316*alpha[0]*fr[3]-15.0*alpha[0]*fr[2]+25.98076211353316*fr[1]*alpha[2]-15.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05*(13.41640786499874*alpha[2]*fr[7]-19.36491673103708*alpha[0]*fr[6]-7.745966692414834*alpha[2]*fr[5]-19.36491673103709*alpha[2]*fr[4]+15.0*alpha[0]*fr[3]-8.660254037844386*alpha[0]*fr[2]+15.0*fr[1]*alpha[2]-8.660254037844386*fr[0]*alpha[2])*dfac_x; 
  incr[4] = 0.08333333333333333*(15.0*alpha[2]*fr[6]+15.0*alpha[0]*fr[4]-11.61895003862225*alpha[2]*fr[3]+6.708203932499369*alpha[2]*fr[2]-11.61895003862225*alpha[0]*fr[1]+6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[5] = -0.01666666666666667*(25.98076211353316*alpha[0]*fr[7]-30.0*alpha[2]*fr[6]-15.0*alpha[0]*fr[5]+23.2379000772445*alpha[2]*fr[3]-13.41640786499874*alpha[2]*fr[2])*dfac_x; 
  incr[6] = -0.01666666666666667*(51.96152422706631*alpha[2]*fr[7]-75.0*alpha[0]*fr[6]-30.0*alpha[2]*fr[5]-75.00000000000001*alpha[2]*fr[4]+58.09475019311126*alpha[0]*fr[3]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*fr[1]*alpha[2]-33.54101966249684*fr[0]*alpha[2])*dfac_x; 
  incr[7] = 0.05*(15.0*alpha[0]*fr[7]-17.32050807568877*alpha[2]*fr[6]-8.660254037844387*alpha[0]*fr[5]+13.41640786499874*alpha[2]*fr[3]-7.745966692414834*alpha[2]*fr[2])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x1vSer_Vpar_P1_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[4]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(4.242640687119286*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x+3.464101615137754*Gradpar[0]*Phi[1]*dfac_x+2.828427124746191*dApardt[0])*q_)/m_; 

  double alpha[4]; 
  alpha[0] = (-(2.121320343559643*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x*q_)/m_)-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*alpha[1]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*alpha[0]*fl[3]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*alpha[1]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[3] = -0.25*(3.0*alpha[0]*fl[3]+3.0*alpha[1]*fl[2]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*alpha[1]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*alpha[0]*fr[3]+1.732050807568877*alpha[1]*fr[2]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(3.0*alpha[1]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[3] = 0.25*(3.0*alpha[0]*fr[3]+3.0*alpha[1]*fr[2]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 

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
double EmGyrokineticSurf1x1vSer_Vpar_P2_Bvars_0(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*((21.21320343559643*geoY[0]*Apar[2]*Phi[2]+4.242640687119286*geoY[0]*Apar[1]*Phi[1])*dfac_x*dfac_x+3.464101615137754*Gradpar[0]*Phi[1]*dfac_x+2.828427124746191*dApardt[0])*q_)/m_; 

  double alpha[8]; 
  alpha[0] = (-(10.60660171779821*geoY[0]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_)-(2.121320343559643*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha[1] = (-(4.74341649025257*geoY[0]*Apar[1]*Phi[2]*dfac_x*dfac_x*q_)/m_)-(4.74341649025257*geoY[0]*Phi[1]*Apar[2]*dfac_x*dfac_x*q_)/m_-(3.872983346207417*Gradpar[0]*Phi[2]*dfac_x*q_)/m_-(1.414213562373095*dApardt[1]*q_)/m_; 
  alpha[4] = (-(9.48683298050514*geoY[0]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_)-(1.414213562373095*dApardt[2]*q_)/m_; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-1.118033988749895*fr[5]*amax)+1.118033988749895*fl[5]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax+0.5590169943749476*alpha[1]*fr[7]+0.5590169943749476*alpha[1]*fl[7]-0.4330127018922194*alpha[4]*fr[6]+0.4330127018922194*alpha[4]*fl[6]+0.5590169943749475*alpha[0]*fr[5]+0.5590169943749475*alpha[0]*fl[5]+0.25*alpha[4]*fr[4]+0.25*alpha[4]*fl[4]-0.4330127018922193*alpha[1]*fr[3]+0.4330127018922193*alpha[1]*fl[3]-0.4330127018922193*alpha[0]*fr[2]+0.4330127018922193*alpha[0]*fl[2]+0.25*alpha[1]*fr[1]+0.25*alpha[1]*fl[1]+0.25*alpha[0]*fr[0]+0.25*alpha[0]*fl[0]; 
  Ghat[1] = (-1.118033988749895*fr[7]*amax)+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax+0.5000000000000001*alpha[4]*fr[7]+0.5590169943749476*alpha[0]*fr[7]+0.5000000000000001*alpha[4]*fl[7]+0.5590169943749476*alpha[0]*fl[7]-0.3872983346207417*alpha[1]*fr[6]+0.3872983346207417*alpha[1]*fl[6]+0.5590169943749475*alpha[1]*fr[5]+0.5590169943749475*alpha[1]*fl[5]+0.223606797749979*alpha[1]*fr[4]+0.223606797749979*alpha[1]*fl[4]-0.3872983346207416*fr[3]*alpha[4]+0.3872983346207416*fl[3]*alpha[4]+0.223606797749979*fr[1]*alpha[4]+0.223606797749979*fl[1]*alpha[4]-0.4330127018922193*alpha[0]*fr[3]+0.4330127018922193*alpha[0]*fl[3]-0.4330127018922193*alpha[1]*fr[2]+0.4330127018922193*alpha[1]*fl[2]+0.25*alpha[0]*fr[1]+0.25*alpha[0]*fl[1]+0.25*fr[0]*alpha[1]+0.25*fl[0]*alpha[1]; 
  Ghat[4] = 0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax+0.5000000000000001*alpha[1]*fr[7]+0.5000000000000001*alpha[1]*fl[7]-0.276641667586244*alpha[4]*fr[6]-0.4330127018922193*alpha[0]*fr[6]+0.276641667586244*alpha[4]*fl[6]+0.4330127018922193*alpha[0]*fl[6]+0.5590169943749475*alpha[4]*fr[5]+0.5590169943749475*alpha[4]*fl[5]+0.159719141249985*alpha[4]*fr[4]+0.25*alpha[0]*fr[4]+0.159719141249985*alpha[4]*fl[4]+0.25*alpha[0]*fl[4]-0.4330127018922193*fr[2]*alpha[4]+0.4330127018922193*fl[2]*alpha[4]+0.25*fr[0]*alpha[4]+0.25*fl[0]*alpha[4]-0.3872983346207416*alpha[1]*fr[3]+0.3872983346207416*alpha[1]*fl[3]+0.223606797749979*alpha[1]*fr[1]+0.223606797749979*alpha[1]*fl[1]; 

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
double GyrokineticSurf1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[4]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(((7.348469228349534*Bmag[1]*BmagInv[1]-4.242640687119286*BmagInv[0]*Bmag[1])*geoY[1]-4.242640687119286*geoY[0]*Bmag[1]*BmagInv[1]+2.449489742783178*BmagInv[0]*geoY[0]*Bmag[1])*dfac_x*m_*wv2+(4.898979485566357*Gradpar[1]-2.828427124746191*Gradpar[0])*q_*wv))/q_; 

  double alpha[4]; 
  alpha[0] = (-(3.674234614174767*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv2)/q_)+(2.121320343559643*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv2)/q_+(2.121320343559643*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv2)/q_-(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv2)/q_-2.449489742783178*Gradpar[1]*wv+1.414213562373095*Gradpar[0]*wv; 
  alpha[2] = (-(2.121320343559642*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_))+(1.224744871391589*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(1.224744871391589*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv)/(dfac_v*q_)-(0.7071067811865475*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*alpha[2]*fl[3]+alpha[2]*fl[2]+1.732050807568877*alpha[0]*fl[1]+alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.25*(3.0*alpha[2]*fl[3]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.25*(1.732050807568877*alpha[0]*fl[3]+alpha[0]*fl[2]+1.732050807568877*fl[1]*alpha[2]+fl[0]*alpha[2])*dfac_x; 
  incr[3] = -0.25*(3.0*alpha[0]*fl[3]+1.732050807568877*alpha[0]*fl[2]+3.0*fl[1]*alpha[2]+1.732050807568877*fl[0]*alpha[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*alpha[2]*fr[3]-1.0*alpha[2]*fr[2]+1.732050807568877*alpha[0]*fr[1]-1.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.25*(3.0*alpha[2]*fr[3]-1.732050807568877*alpha[2]*fr[2]+3.0*alpha[0]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.25*(1.732050807568877*alpha[0]*fr[3]-1.0*alpha[0]*fr[2]+1.732050807568877*fr[1]*alpha[2]-1.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = 0.25*(3.0*alpha[0]*fr[3]-1.732050807568877*alpha[0]*fr[2]+3.0*fr[1]*alpha[2]-1.732050807568877*fr[0]*alpha[2])*dfac_x; 

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
double GyrokineticSurf1x1vSer_X_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.08333333333333333*((((142.3024947075771*Bmag[2]-36.74234614174767*Bmag[1])*BmagInv[2]+(63.63961030678928*BmagInv[0]-110.227038425243*BmagInv[1])*Bmag[2]+28.46049894151542*Bmag[1]*BmagInv[1]-16.43167672515498*BmagInv[0]*Bmag[1])*geoY[2]+((63.63961030678928*geoY[0]-110.227038425243*geoY[1])*Bmag[2]+28.46049894151542*Bmag[1]*geoY[1]-16.43167672515498*geoY[0]*Bmag[1])*BmagInv[2]+((85.38149682454625*BmagInv[1]-49.29503017546494*BmagInv[0])*geoY[1]-49.29503017546494*geoY[0]*BmagInv[1]+28.46049894151542*BmagInv[0]*geoY[0])*Bmag[2]+(12.72792206135786*BmagInv[0]*Bmag[1]-22.0454076850486*Bmag[1]*BmagInv[1])*geoY[1]+12.72792206135786*geoY[0]*Bmag[1]*BmagInv[1]-7.348469228349534*BmagInv[0]*geoY[0]*Bmag[1])*dfac_v2*dfac_x*m_*wv2+(18.97366596101028*Gradpar[2]-14.69693845669907*Gradpar[1]+8.485281374238571*Gradpar[0])*dfac_v2*q_*wv+(((47.43416490252571*Bmag[2]-12.24744871391589*Bmag[1])*BmagInv[2]+(21.21320343559643*BmagInv[0]-36.74234614174767*BmagInv[1])*Bmag[2]+9.48683298050514*Bmag[1]*BmagInv[1]-5.477225575051662*BmagInv[0]*Bmag[1])*geoY[2]+((21.21320343559643*geoY[0]-36.74234614174767*geoY[1])*Bmag[2]+9.48683298050514*Bmag[1]*geoY[1]-5.477225575051662*geoY[0]*Bmag[1])*BmagInv[2]+((28.46049894151542*BmagInv[1]-16.43167672515498*BmagInv[0])*geoY[1]-16.43167672515498*geoY[0]*BmagInv[1]+9.48683298050514*BmagInv[0]*geoY[0])*Bmag[2]+(4.242640687119286*BmagInv[0]*Bmag[1]-7.348469228349534*Bmag[1]*BmagInv[1])*geoY[1]+4.242640687119286*geoY[0]*Bmag[1]*BmagInv[1]-2.449489742783178*BmagInv[0]*geoY[0]*Bmag[1])*dfac_x*m_))/(dfac_v2*q_); 

  double alpha[8]; 
  alpha[0] = (23.71708245126285*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_*wv2)/q_-(6.123724356957945*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_*wv2)/q_-(18.37117307087383*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_*wv2)/q_+(10.60660171779821*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_*wv2)/q_+(4.743416490252569*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_*wv2)/q_-(2.738612787525831*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_*wv2)/q_-(18.37117307087383*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv2)/q_+(10.60660171779821*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv2)/q_+(4.743416490252569*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_*wv2)/q_-(2.738612787525831*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_*wv2)/q_+(14.23024947075771*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_*wv2)/q_-(8.215838362577491*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_*wv2)/q_-(8.215838362577491*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_*wv2)/q_+(4.743416490252569*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_*wv2)/q_-(3.674234614174767*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv2)/q_+(2.121320343559643*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv2)/q_+(2.121320343559643*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv2)/q_-(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv2)/q_+3.16227766016838*Gradpar[2]*wv-2.449489742783178*Gradpar[1]*wv+1.414213562373095*Gradpar[0]*wv+(7.90569415042095*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(2.041241452319315*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(6.123724356957944*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(3.535533905932738*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(1.58113883008419*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(0.912870929175277*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(6.123724356957944*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(3.535533905932738*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(1.58113883008419*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)-(0.912870929175277*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(4.743416490252569*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(2.73861278752583*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(2.73861278752583*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)+(1.58113883008419*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(1.224744871391589*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.7071067811865476*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.7071067811865476*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_)/(dfac_v2*q_)-(0.408248290463863*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_)/(dfac_v2*q_); 
  alpha[2] = (27.38612787525831*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(7.071067811865476*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(21.21320343559643*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)+(12.24744871391589*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)+(5.477225575051662*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(3.16227766016838*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(21.21320343559643*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(12.24744871391589*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(5.477225575051662*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)-(3.16227766016838*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(16.43167672515498*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(9.486832980505138*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(9.486832980505138*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)+(5.477225575051662*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(4.242640687119286*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(2.449489742783178*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(2.449489742783178*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv)/(dfac_v*q_)-(1.414213562373095*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv)/(dfac_v*q_)+(1.825741858350554*Gradpar[2])/dfac_v-(1.414213562373095*Gradpar[1])/dfac_v+(0.816496580927726*Gradpar[0])/dfac_v; 
  alpha[5] = (7.071067811865476*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(1.825741858350554*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(5.477225575051662*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(3.16227766016838*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(1.414213562373095*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(0.8164965809277259*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(5.477225575051662*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(3.16227766016838*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(1.414213562373095*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)-(0.8164965809277259*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(4.242640687119285*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(2.449489742783178*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(2.449489742783178*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)+(1.414213562373095*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(1.095445115010332*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.6324555320336759*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.6324555320336759*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_)/(dfac_v2*q_)-(0.3651483716701108*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_)/(dfac_v2*q_); 
  if (alpha0>0) { 
  incr[0] = 0.01666666666666667*(25.98076211353316*alpha[5]*fl[7]+33.54101966249684*alpha[2]*fl[6]+15.0*alpha[5]*fl[5]+33.54101966249685*alpha[0]*fl[4]+25.98076211353316*alpha[2]*fl[3]+15.0*alpha[2]*fl[2]+25.98076211353316*alpha[0]*fl[1]+15.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.05*(15.0*alpha[5]*fl[7]+19.36491673103708*alpha[2]*fl[6]+8.660254037844386*alpha[5]*fl[5]+19.36491673103709*alpha[0]*fl[4]+15.0*alpha[2]*fl[3]+8.660254037844386*alpha[2]*fl[2]+15.0*alpha[0]*fl[1]+8.660254037844386*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01666666666666667*(23.2379000772445*alpha[2]*fl[7]+30.0*alpha[5]*fl[6]+33.54101966249684*alpha[0]*fl[6]+13.41640786499874*alpha[2]*fl[5]+23.2379000772445*fl[3]*alpha[5]+13.41640786499874*fl[2]*alpha[5]+33.54101966249685*alpha[2]*fl[4]+25.98076211353316*alpha[0]*fl[3]+15.0*alpha[0]*fl[2]+25.98076211353316*fl[1]*alpha[2]+15.0*fl[0]*alpha[2])*dfac_x; 
  incr[3] = -0.05*(13.41640786499874*alpha[2]*fl[7]+17.32050807568877*alpha[5]*fl[6]+19.36491673103708*alpha[0]*fl[6]+7.745966692414834*alpha[2]*fl[5]+13.41640786499874*fl[3]*alpha[5]+7.745966692414834*fl[2]*alpha[5]+19.36491673103709*alpha[2]*fl[4]+15.0*alpha[0]*fl[3]+8.660254037844386*alpha[0]*fl[2]+15.0*fl[1]*alpha[2]+8.660254037844386*fl[0]*alpha[2])*dfac_x; 
  incr[4] = 0.08333333333333333*(11.61895003862225*alpha[5]*fl[7]+15.0*alpha[2]*fl[6]+6.708203932499369*alpha[5]*fl[5]+15.0*alpha[0]*fl[4]+11.61895003862225*alpha[2]*fl[3]+6.708203932499369*alpha[2]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[5] = 0.002380952380952381*(116.1895003862225*alpha[5]*fl[7]+181.8653347947321*alpha[0]*fl[7]+210.0*alpha[2]*fl[6]+67.0820393249937*alpha[5]*fl[5]+105.0*alpha[0]*fl[5]+234.787137637478*fl[4]*alpha[5]+181.8653347947321*fl[1]*alpha[5]+105.0*fl[0]*alpha[5]+162.6653005407115*alpha[2]*fl[3]+93.91485505499116*alpha[2]*fl[2])*dfac_x; 
  incr[6] = 0.01666666666666667*(51.96152422706631*alpha[2]*fl[7]+67.0820393249937*alpha[5]*fl[6]+75.0*alpha[0]*fl[6]+30.0*alpha[2]*fl[5]+51.96152422706632*fl[3]*alpha[5]+30.0*fl[2]*alpha[5]+75.00000000000001*alpha[2]*fl[4]+58.09475019311126*alpha[0]*fl[3]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*fl[1]*alpha[2]+33.54101966249684*fl[0]*alpha[2])*dfac_x; 
  incr[7] = -0.007142857142857143*(67.0820393249937*alpha[5]*fl[7]+105.0*alpha[0]*fl[7]+121.2435565298214*alpha[2]*fl[6]+38.72983346207417*alpha[5]*fl[5]+60.62177826491071*alpha[0]*fl[5]+135.5544171172596*fl[4]*alpha[5]+105.0*fl[1]*alpha[5]+60.62177826491071*fl[0]*alpha[5]+93.91485505499116*alpha[2]*fl[3]+54.22176684690384*alpha[2]*fl[2])*dfac_x; 

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
  } else { 
  incr[0] = -0.01666666666666667*(25.98076211353316*alpha[5]*fr[7]-33.54101966249684*alpha[2]*fr[6]-15.0*alpha[5]*fr[5]-33.54101966249685*alpha[0]*fr[4]+25.98076211353316*alpha[2]*fr[3]-15.0*alpha[2]*fr[2]+25.98076211353316*alpha[0]*fr[1]-15.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.05*(15.0*alpha[5]*fr[7]-19.36491673103708*alpha[2]*fr[6]-8.660254037844386*alpha[5]*fr[5]-19.36491673103709*alpha[0]*fr[4]+15.0*alpha[2]*fr[3]-8.660254037844386*alpha[2]*fr[2]+15.0*alpha[0]*fr[1]-8.660254037844386*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01666666666666667*(23.2379000772445*alpha[2]*fr[7]-30.0*alpha[5]*fr[6]-33.54101966249684*alpha[0]*fr[6]-13.41640786499874*alpha[2]*fr[5]+23.2379000772445*fr[3]*alpha[5]-13.41640786499874*fr[2]*alpha[5]-33.54101966249685*alpha[2]*fr[4]+25.98076211353316*alpha[0]*fr[3]-15.0*alpha[0]*fr[2]+25.98076211353316*fr[1]*alpha[2]-15.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05*(13.41640786499874*alpha[2]*fr[7]-17.32050807568877*alpha[5]*fr[6]-19.36491673103708*alpha[0]*fr[6]-7.745966692414834*alpha[2]*fr[5]+13.41640786499874*fr[3]*alpha[5]-7.745966692414834*fr[2]*alpha[5]-19.36491673103709*alpha[2]*fr[4]+15.0*alpha[0]*fr[3]-8.660254037844386*alpha[0]*fr[2]+15.0*fr[1]*alpha[2]-8.660254037844386*fr[0]*alpha[2])*dfac_x; 
  incr[4] = -0.08333333333333333*(11.61895003862225*alpha[5]*fr[7]-15.0*alpha[2]*fr[6]-6.708203932499369*alpha[5]*fr[5]-15.0*alpha[0]*fr[4]+11.61895003862225*alpha[2]*fr[3]-6.708203932499369*alpha[2]*fr[2]+11.61895003862225*alpha[0]*fr[1]-6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[5] = -0.002380952380952381*(116.1895003862225*alpha[5]*fr[7]+181.8653347947321*alpha[0]*fr[7]-210.0*alpha[2]*fr[6]-67.0820393249937*alpha[5]*fr[5]-105.0*alpha[0]*fr[5]-234.787137637478*fr[4]*alpha[5]+181.8653347947321*fr[1]*alpha[5]-105.0*fr[0]*alpha[5]+162.6653005407115*alpha[2]*fr[3]-93.91485505499116*alpha[2]*fr[2])*dfac_x; 
  incr[6] = -0.01666666666666667*(51.96152422706631*alpha[2]*fr[7]-67.0820393249937*alpha[5]*fr[6]-75.0*alpha[0]*fr[6]-30.0*alpha[2]*fr[5]+51.96152422706632*fr[3]*alpha[5]-30.0*fr[2]*alpha[5]-75.00000000000001*alpha[2]*fr[4]+58.09475019311126*alpha[0]*fr[3]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*fr[1]*alpha[2]-33.54101966249684*fr[0]*alpha[2])*dfac_x; 
  incr[7] = 0.007142857142857143*(67.0820393249937*alpha[5]*fr[7]+105.0*alpha[0]*fr[7]-121.2435565298214*alpha[2]*fr[6]-38.72983346207417*alpha[5]*fr[5]-60.62177826491071*alpha[0]*fr[5]-135.5544171172596*fr[4]*alpha[5]+105.0*fr[1]*alpha[5]-60.62177826491071*fr[0]*alpha[5]+93.91485505499116*alpha[2]*fr[3]-54.22176684690384*alpha[2]*fr[2])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double GyrokineticSurf1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[4]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.25*((3.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]+3.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_v*dfac_x*dfac_x*m_*wv-3.464101615137754*Gradpar[0]*Phi[1]*dfac_v*dfac_x*q_+((-3.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1])-3.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_x*dfac_x*m_))/(dfac_v*m_); 

  double alpha[4]; 
  alpha[0] = 1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*wv-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[1] = 1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*wv-(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x*q_)/m_-(1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*alpha[1]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*alpha[0]*fl[3]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*alpha[1]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[3] = -0.25*(3.0*alpha[0]*fl[3]+3.0*alpha[1]*fl[2]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*alpha[1]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*alpha[0]*fr[3]+1.732050807568877*alpha[1]*fr[2]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(3.0*alpha[1]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[3] = 0.25*(3.0*alpha[0]*fr[3]+3.0*alpha[1]*fr[2]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 

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
double GyrokineticSurf1x1vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.03571428571428571*((((165.0*Bmag[2]*BmagInv[2]+93.91485505499116*BmagInv[0]*Bmag[2]+42.0*Bmag[1]*BmagInv[1])*Phi[2]+21.0*Bmag[1]*Phi[1]*BmagInv[2]+42.0*BmagInv[1]*Phi[1]*Bmag[2])*geoY[2]+((93.91485505499116*geoY[0]*Bmag[2]+42.0*Bmag[1]*geoY[1])*BmagInv[2]+(189.0*BmagInv[1]*geoY[1]+105.0*BmagInv[0]*geoY[0])*Bmag[2]+46.95742752749558*BmagInv[0]*Bmag[1]*geoY[1]+46.95742752749558*geoY[0]*Bmag[1]*BmagInv[1])*Phi[2]+42.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]+(46.95742752749558*BmagInv[0]*Phi[1]*geoY[1]+46.95742752749558*geoY[0]*BmagInv[1]*Phi[1])*Bmag[2]+21.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]+21.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_v*dfac_x*dfac_x*m_*wv+((-54.22176684690384*Gradpar[1]*Phi[2])-24.24871130596428*Gradpar[0]*Phi[1])*dfac_v*dfac_x*q_+((((-165.0*Bmag[2]*BmagInv[2])-93.91485505499116*BmagInv[0]*Bmag[2]-42.0*Bmag[1]*BmagInv[1])*Phi[2]-21.0*Bmag[1]*Phi[1]*BmagInv[2]-42.0*BmagInv[1]*Phi[1]*Bmag[2])*geoY[2]+(((-93.91485505499116*geoY[0]*Bmag[2])-42.0*Bmag[1]*geoY[1])*BmagInv[2]+((-189.0*BmagInv[1]*geoY[1])-105.0*BmagInv[0]*geoY[0])*Bmag[2]-46.95742752749558*BmagInv[0]*Bmag[1]*geoY[1]-46.95742752749558*geoY[0]*Bmag[1]*BmagInv[1])*Phi[2]-42.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]+((-46.95742752749558*BmagInv[0]*Phi[1]*geoY[1])-46.95742752749558*geoY[0]*BmagInv[1]*Phi[1])*Bmag[2]-21.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]-21.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_x*dfac_x*m_))/(dfac_v*m_); 

  double alpha[8]; 
  alpha[0] = 11.78571428571428*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+6.708203932499369*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+3.0*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+1.5*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+3.0*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+6.708203932499369*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+3.0*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+13.5*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+7.5*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+3.354101966249684*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+3.354101966249684*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*wv+3.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+3.354101966249684*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+3.354101966249684*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*wv-(3.872983346207417*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(11.78571428571428*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(6.708203932499369*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(3.0*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(1.5*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(3.0*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(6.708203932499369*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(13.5*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(7.5*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249684*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249684*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249684*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249684*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[1] = 5.270731661249505*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+17.24966725499838*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+3.0*BmagInv[0]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+5.270731661249505*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+3.0*BmagInv[0]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+1.341640786499874*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*wv+17.24966725499838*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+3.0*geoY[0]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+13.5*BmagInv[0]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+13.5*geoY[0]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+6.037383539249432*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+3.354101966249686*BmagInv[0]*geoY[0]*Bmag[1]*Phi[2]*dfac_x*dfac_x*wv+3.0*geoY[0]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+1.341640786499874*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wv+6.037383539249432*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+3.354101966249686*BmagInv[0]*geoY[0]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*wv-(3.464101615137754*Gradpar[2]*Phi[2]*dfac_x*q_)/m_-(3.872983346207417*Gradpar[0]*Phi[2]*dfac_x*q_)/m_-(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x*q_)/m_-(5.270731661249505*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(17.24966725499838*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(3.0*BmagInv[0]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(5.270731661249505*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(3.0*BmagInv[0]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(1.341640786499874*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x)/dfac_v-(17.24966725499838*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*geoY[0]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(13.5*BmagInv[0]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(13.5*geoY[0]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(6.037383539249432*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249686*BmagInv[0]*geoY[0]*Bmag[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*geoY[0]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(1.341640786499874*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(6.037383539249432*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249686*BmagInv[0]*geoY[0]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[4] = 19.1662969499982*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+11.78571428571428*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+5.270731661249505*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+0.95831484749991*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+5.270731661249505*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*wv+11.78571428571428*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+5.270731661249505*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+17.24966725499838*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+6.708203932499371*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+3.0*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+3.0*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*wv+5.270731661249505*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+1.5*geoY[0]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x*wv+3.0*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+3.0*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+1.341640786499874*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv-(3.464101615137754*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(1.732050807568877*Phi[1]*Gradpar[2]*dfac_x*q_)/m_-(19.1662969499982*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(11.78571428571428*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(5.270731661249505*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(0.95831484749991*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(5.270731661249505*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x)/dfac_v-(11.78571428571428*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(5.270731661249505*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(17.24966725499838*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(6.708203932499371*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(5.270731661249505*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(1.5*geoY[0]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(3.0*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(3.0*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(1.341640786499874*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-1.118033988749895*fr[5]*amax)+1.118033988749895*fl[5]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax+0.5590169943749476*alpha[1]*fr[7]+0.5590169943749476*alpha[1]*fl[7]-0.4330127018922194*alpha[4]*fr[6]+0.4330127018922194*alpha[4]*fl[6]+0.5590169943749475*alpha[0]*fr[5]+0.5590169943749475*alpha[0]*fl[5]+0.25*alpha[4]*fr[4]+0.25*alpha[4]*fl[4]-0.4330127018922193*alpha[1]*fr[3]+0.4330127018922193*alpha[1]*fl[3]-0.4330127018922193*alpha[0]*fr[2]+0.4330127018922193*alpha[0]*fl[2]+0.25*alpha[1]*fr[1]+0.25*alpha[1]*fl[1]+0.25*alpha[0]*fr[0]+0.25*alpha[0]*fl[0]; 
  Ghat[1] = (-1.118033988749895*fr[7]*amax)+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax+0.5000000000000001*alpha[4]*fr[7]+0.5590169943749476*alpha[0]*fr[7]+0.5000000000000001*alpha[4]*fl[7]+0.5590169943749476*alpha[0]*fl[7]-0.3872983346207417*alpha[1]*fr[6]+0.3872983346207417*alpha[1]*fl[6]+0.5590169943749475*alpha[1]*fr[5]+0.5590169943749475*alpha[1]*fl[5]+0.223606797749979*alpha[1]*fr[4]+0.223606797749979*alpha[1]*fl[4]-0.3872983346207416*fr[3]*alpha[4]+0.3872983346207416*fl[3]*alpha[4]+0.223606797749979*fr[1]*alpha[4]+0.223606797749979*fl[1]*alpha[4]-0.4330127018922193*alpha[0]*fr[3]+0.4330127018922193*alpha[0]*fl[3]-0.4330127018922193*alpha[1]*fr[2]+0.4330127018922193*alpha[1]*fl[2]+0.25*alpha[0]*fr[1]+0.25*alpha[0]*fl[1]+0.25*fr[0]*alpha[1]+0.25*fl[0]*alpha[1]; 
  Ghat[4] = 0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax+0.5000000000000001*alpha[1]*fr[7]+0.5000000000000001*alpha[1]*fl[7]-0.276641667586244*alpha[4]*fr[6]-0.4330127018922193*alpha[0]*fr[6]+0.276641667586244*alpha[4]*fl[6]+0.4330127018922193*alpha[0]*fl[6]+0.5590169943749475*alpha[4]*fr[5]+0.5590169943749475*alpha[4]*fl[5]+0.159719141249985*alpha[4]*fr[4]+0.25*alpha[0]*fr[4]+0.159719141249985*alpha[4]*fl[4]+0.25*alpha[0]*fl[4]-0.4330127018922193*fr[2]*alpha[4]+0.4330127018922193*fl[2]*alpha[4]+0.25*fr[0]*alpha[4]+0.25*fl[0]*alpha[4]-0.3872983346207416*alpha[1]*fr[3]+0.3872983346207416*alpha[1]*fl[3]+0.223606797749979*alpha[1]*fr[1]+0.223606797749979*alpha[1]*fl[1]; 

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
double EmGyrokineticSurf1x1vSer_X_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[4]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = -(0.25*(((7.348469228349534*Bmag[1]*BmagInv[1]-4.242640687119286*BmagInv[0]*Bmag[1])*geoY[1]-4.242640687119286*geoY[0]*Bmag[1]*BmagInv[1]+2.449489742783178*BmagInv[0]*geoY[0]*Bmag[1])*dfac_x*m_*wv2+((((5.196152422706631*Apar[0]-9.0*Apar[1])*Bmag[1]*BmagInv[1]+(5.196152422706631*BmagInv[0]*Apar[1]-3.0*Apar[0]*BmagInv[0])*Bmag[1]+6.0*Apar[1])*geoY[1]+(5.196152422706631*geoY[0]*Apar[1]-3.0*Apar[0]*geoY[0])*Bmag[1]*BmagInv[1]+(1.732050807568877*Apar[0]*BmagInv[0]*geoY[0]-3.0*BmagInv[0]*geoY[0]*Apar[1])*Bmag[1]-3.464101615137754*geoY[0]*Apar[1])*dfac_x+4.898979485566357*Gradpar[1]-2.828427124746191*Gradpar[0])*q_*wv))/q_; 

  double alpha[4]; 
  alpha[0] = (-(3.674234614174767*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv2)/q_)+(2.121320343559643*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv2)/q_+(2.121320343559643*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv2)/q_-(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv2)/q_+4.5*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*wv-2.598076211353316*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*wv-2.598076211353316*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*dfac_x*wv+1.5*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*wv-3.0*Apar[1]*geoY[1]*dfac_x*wv-2.598076211353316*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*dfac_x*wv+1.5*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*wv+1.5*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*dfac_x*wv-0.8660254037844386*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*wv+1.732050807568877*geoY[0]*Apar[1]*dfac_x*wv-2.449489742783178*Gradpar[1]*wv+1.414213562373095*Gradpar[0]*wv; 
  alpha[2] = (-(2.121320343559642*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_))+(1.224744871391589*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(1.224744871391589*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv)/(dfac_v*q_)-(0.7071067811865475*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv)/(dfac_v*q_); 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*alpha[2]*fl[3]+alpha[2]*fl[2]+1.732050807568877*alpha[0]*fl[1]+alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.25*(3.0*alpha[2]*fl[3]+1.732050807568877*alpha[2]*fl[2]+3.0*alpha[0]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.25*(1.732050807568877*alpha[0]*fl[3]+alpha[0]*fl[2]+1.732050807568877*fl[1]*alpha[2]+fl[0]*alpha[2])*dfac_x; 
  incr[3] = -0.25*(3.0*alpha[0]*fl[3]+1.732050807568877*alpha[0]*fl[2]+3.0*fl[1]*alpha[2]+1.732050807568877*fl[0]*alpha[2])*dfac_x; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += incr[1]; 
  outl[2] += -1.0*incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*alpha[2]*fr[3]-1.0*alpha[2]*fr[2]+1.732050807568877*alpha[0]*fr[1]-1.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.25*(3.0*alpha[2]*fr[3]-1.732050807568877*alpha[2]*fr[2]+3.0*alpha[0]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.25*(1.732050807568877*alpha[0]*fr[3]-1.0*alpha[0]*fr[2]+1.732050807568877*fr[1]*alpha[2]-1.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = 0.25*(3.0*alpha[0]*fr[3]-1.732050807568877*alpha[0]*fr[2]+3.0*fr[1]*alpha[2]-1.732050807568877*fr[0]*alpha[2])*dfac_x; 

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
double EmGyrokineticSurf1x1vSer_X_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.08333333333333333*((((142.3024947075771*Bmag[2]-36.74234614174767*Bmag[1])*BmagInv[2]+(63.63961030678928*BmagInv[0]-110.227038425243*BmagInv[1])*Bmag[2]+28.46049894151542*Bmag[1]*BmagInv[1]-16.43167672515498*BmagInv[0]*Bmag[1])*geoY[2]+((63.63961030678928*geoY[0]-110.227038425243*geoY[1])*Bmag[2]+28.46049894151542*Bmag[1]*geoY[1]-16.43167672515498*geoY[0]*Bmag[1])*BmagInv[2]+((85.38149682454625*BmagInv[1]-49.29503017546494*BmagInv[0])*geoY[1]-49.29503017546494*geoY[0]*BmagInv[1]+28.46049894151542*BmagInv[0]*geoY[0])*Bmag[2]+(12.72792206135786*BmagInv[0]*Bmag[1]-22.0454076850486*Bmag[1]*BmagInv[1])*geoY[1]+12.72792206135786*geoY[0]*Bmag[1]*BmagInv[1]-7.348469228349534*BmagInv[0]*geoY[0]*Bmag[1])*dfac_v2*dfac_x*m_*wv2+(((((225.0*Apar[2]-174.2842505793337*Apar[1]+100.6230589874906*Apar[0])*Bmag[2]-58.09475019311126*Bmag[1]*Apar[2]+(45.0*Apar[1]-25.98076211353316*Apar[0])*Bmag[1])*BmagInv[2]+((100.6230589874906*BmagInv[0]-174.2842505793337*BmagInv[1])*Apar[2]+(135.0*Apar[1]-77.94228634059945*Apar[0])*BmagInv[1]-77.94228634059945*BmagInv[0]*Apar[1]+45.0*Apar[0]*BmagInv[0])*Bmag[2]+(45.0*Bmag[1]*BmagInv[1]-25.98076211353316*BmagInv[0]*Bmag[1]-90.0)*Apar[2]+(20.12461179749811*Apar[0]-34.85685011586674*Apar[1])*Bmag[1]*BmagInv[1]+(20.12461179749811*BmagInv[0]*Apar[1]-11.61895003862225*Apar[0]*BmagInv[0])*Bmag[1]+23.2379000772445*Apar[1])*geoY[2]+(((100.6230589874906*geoY[0]-174.2842505793337*geoY[1])*Apar[2]+(135.0*Apar[1]-77.94228634059945*Apar[0])*geoY[1]-77.94228634059945*geoY[0]*Apar[1]+45.0*Apar[0]*geoY[0])*Bmag[2]+(45.0*Bmag[1]*geoY[1]-25.98076211353316*geoY[0]*Bmag[1])*Apar[2]+(20.12461179749811*Apar[0]-34.85685011586674*Apar[1])*Bmag[1]*geoY[1]+(20.12461179749811*geoY[0]*Apar[1]-11.61895003862225*Apar[0]*geoY[0])*Bmag[1])*BmagInv[2]+(((135.0*BmagInv[1]-77.94228634059945*BmagInv[0])*geoY[1]-77.94228634059945*geoY[0]*BmagInv[1]+45.0*BmagInv[0]*geoY[0])*Apar[2]+((60.37383539249433*Apar[0]-104.5705503476002*Apar[1])*BmagInv[1]+60.37383539249433*BmagInv[0]*Apar[1]-34.85685011586674*Apar[0]*BmagInv[0])*geoY[1]+(60.37383539249433*geoY[0]*Apar[1]-34.85685011586674*Apar[0]*geoY[0])*BmagInv[1]-34.85685011586674*BmagInv[0]*geoY[0]*Apar[1]+20.12461179749811*Apar[0]*BmagInv[0]*geoY[0])*Bmag[2]+(((-34.85685011586674*Bmag[1]*BmagInv[1])+20.12461179749811*BmagInv[0]*Bmag[1]+69.71370023173348)*geoY[1]+20.12461179749811*geoY[0]*Bmag[1]*BmagInv[1]-11.61895003862225*BmagInv[0]*geoY[0]*Bmag[1]-40.24922359499622*geoY[0])*Apar[2]+((27.0*Apar[1]-15.58845726811989*Apar[0])*Bmag[1]*BmagInv[1]+(9.0*Apar[0]*BmagInv[0]-15.58845726811989*BmagInv[0]*Apar[1])*Bmag[1]-18.0*Apar[1])*geoY[1]+(9.0*Apar[0]*geoY[0]-15.58845726811989*geoY[0]*Apar[1])*Bmag[1]*BmagInv[1]+(9.0*BmagInv[0]*geoY[0]*Apar[1]-5.196152422706631*Apar[0]*BmagInv[0]*geoY[0])*Bmag[1]+10.39230484541326*geoY[0]*Apar[1])*dfac_v2*dfac_x+(18.97366596101028*Gradpar[2]-14.69693845669907*Gradpar[1]+8.485281374238571*Gradpar[0])*dfac_v2)*q_*wv+(((47.43416490252571*Bmag[2]-12.24744871391589*Bmag[1])*BmagInv[2]+(21.21320343559643*BmagInv[0]-36.74234614174767*BmagInv[1])*Bmag[2]+9.48683298050514*Bmag[1]*BmagInv[1]-5.477225575051662*BmagInv[0]*Bmag[1])*geoY[2]+((21.21320343559643*geoY[0]-36.74234614174767*geoY[1])*Bmag[2]+9.48683298050514*Bmag[1]*geoY[1]-5.477225575051662*geoY[0]*Bmag[1])*BmagInv[2]+((28.46049894151542*BmagInv[1]-16.43167672515498*BmagInv[0])*geoY[1]-16.43167672515498*geoY[0]*BmagInv[1]+9.48683298050514*BmagInv[0]*geoY[0])*Bmag[2]+(4.242640687119286*BmagInv[0]*Bmag[1]-7.348469228349534*Bmag[1]*BmagInv[1])*geoY[1]+4.242640687119286*geoY[0]*Bmag[1]*BmagInv[1]-2.449489742783178*BmagInv[0]*geoY[0]*Bmag[1])*dfac_x*m_))/(dfac_v2*q_); 

  double alpha[8]; 
  alpha[0] = (23.71708245126285*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_*wv2)/q_-(6.123724356957945*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_*wv2)/q_-(18.37117307087383*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_*wv2)/q_+(10.60660171779821*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_*wv2)/q_+(4.743416490252569*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_*wv2)/q_-(2.738612787525831*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_*wv2)/q_-(18.37117307087383*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv2)/q_+(10.60660171779821*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv2)/q_+(4.743416490252569*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_*wv2)/q_-(2.738612787525831*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_*wv2)/q_+(14.23024947075771*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_*wv2)/q_-(8.21583836257749*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_*wv2)/q_-(8.21583836257749*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_*wv2)/q_+(4.743416490252569*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_*wv2)/q_-(3.674234614174767*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv2)/q_+(2.121320343559643*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv2)/q_+(2.121320343559643*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv2)/q_-(1.224744871391589*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv2)/q_+37.5*Apar[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*wv-29.04737509655563*Apar[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*wv+16.77050983124843*Apar[0]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*wv-9.682458365518544*Bmag[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x*wv+7.5*Apar[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*wv-4.330127018922193*Apar[0]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*wv-29.04737509655563*BmagInv[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*wv+16.77050983124843*BmagInv[0]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*wv+22.5*Apar[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*wv-12.99038105676658*Apar[0]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*wv-12.99038105676658*BmagInv[0]*Apar[1]*Bmag[2]*geoY[2]*dfac_x*wv+7.5*Apar[0]*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*wv+7.5*Bmag[1]*BmagInv[1]*Apar[2]*geoY[2]*dfac_x*wv-4.330127018922193*BmagInv[0]*Bmag[1]*Apar[2]*geoY[2]*dfac_x*wv-15.0*Apar[2]*geoY[2]*dfac_x*wv-5.809475019311124*Apar[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*wv+3.354101966249685*Apar[0]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*wv+3.354101966249685*BmagInv[0]*Apar[1]*Bmag[1]*geoY[2]*dfac_x*wv-1.936491673103709*Apar[0]*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*wv+3.872983346207417*Apar[1]*geoY[2]*dfac_x*wv-29.04737509655563*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*wv+16.77050983124843*geoY[0]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*wv+22.5*Apar[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*wv-12.99038105676658*Apar[0]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*wv-12.99038105676658*geoY[0]*Apar[1]*Bmag[2]*BmagInv[2]*dfac_x*wv+7.5*Apar[0]*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*wv+7.5*Bmag[1]*geoY[1]*Apar[2]*BmagInv[2]*dfac_x*wv-4.330127018922193*geoY[0]*Bmag[1]*Apar[2]*BmagInv[2]*dfac_x*wv-5.809475019311124*Apar[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*wv+3.354101966249685*Apar[0]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*wv+3.354101966249685*geoY[0]*Apar[1]*Bmag[1]*BmagInv[2]*dfac_x*wv-1.936491673103709*Apar[0]*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*wv+22.5*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*wv-12.99038105676658*BmagInv[0]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*wv-12.99038105676658*geoY[0]*BmagInv[1]*Apar[2]*Bmag[2]*dfac_x*wv+7.5*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*dfac_x*wv-17.42842505793337*Apar[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*wv+10.06230589874905*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*wv+10.06230589874905*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*dfac_x*wv-5.809475019311124*Apar[0]*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*wv+10.06230589874905*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*dfac_x*wv-5.809475019311124*Apar[0]*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*wv-5.809475019311124*BmagInv[0]*geoY[0]*Apar[1]*Bmag[2]*dfac_x*wv+3.354101966249685*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*wv-5.809475019311124*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*dfac_x*wv+3.354101966249685*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*dfac_x*wv+11.61895003862225*geoY[1]*Apar[2]*dfac_x*wv+3.354101966249685*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*dfac_x*wv-1.936491673103709*BmagInv[0]*geoY[0]*Bmag[1]*Apar[2]*dfac_x*wv-6.708203932499369*geoY[0]*Apar[2]*dfac_x*wv+4.5*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*wv-2.598076211353315*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*wv-2.598076211353315*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*dfac_x*wv+1.5*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*wv-3.0*Apar[1]*geoY[1]*dfac_x*wv-2.598076211353315*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*dfac_x*wv+1.5*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*wv+1.5*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*dfac_x*wv-0.8660254037844386*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*wv+1.732050807568877*geoY[0]*Apar[1]*dfac_x*wv+3.16227766016838*Gradpar[2]*wv-2.449489742783178*Gradpar[1]*wv+1.414213562373095*Gradpar[0]*wv+(7.90569415042095*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(2.041241452319315*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(6.123724356957944*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(3.535533905932738*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(1.58113883008419*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(0.912870929175277*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(6.123724356957944*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(3.535533905932738*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(1.58113883008419*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)-(0.912870929175277*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(4.743416490252569*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(2.73861278752583*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(2.73861278752583*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)+(1.58113883008419*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(1.224744871391589*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.7071067811865476*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.7071067811865476*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_)/(dfac_v2*q_)-(0.408248290463863*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_)/(dfac_v2*q_); 
  alpha[2] = (27.38612787525831*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(7.071067811865477*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(21.21320343559643*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)+(12.24744871391589*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)+(5.477225575051662*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(3.16227766016838*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_*wv)/(dfac_v*q_)-(21.21320343559643*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(12.24744871391589*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(5.477225575051662*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)-(3.16227766016838*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_*wv)/(dfac_v*q_)+(16.43167672515499*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(9.486832980505138*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(9.486832980505138*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)+(5.477225575051662*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_*wv)/(dfac_v*q_)-(4.242640687119286*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(2.449489742783178*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_*wv)/(dfac_v*q_)+(2.449489742783178*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_*wv)/(dfac_v*q_)-(1.414213562373095*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_*wv)/(dfac_v*q_)+(21.65063509461097*Apar[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v-(16.77050983124843*Apar[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v+(9.682458365518544*Apar[0]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v-(5.590169943749475*Bmag[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v+(4.330127018922193*Apar[1]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v-(2.5*Apar[0]*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x)/dfac_v-(16.77050983124843*BmagInv[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x)/dfac_v+(9.682458365518544*BmagInv[0]*Apar[2]*Bmag[2]*geoY[2]*dfac_x)/dfac_v+(12.99038105676658*Apar[1]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x)/dfac_v-(7.5*Apar[0]*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x)/dfac_v-(7.5*BmagInv[0]*Apar[1]*Bmag[2]*geoY[2]*dfac_x)/dfac_v+(4.330127018922193*Apar[0]*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x)/dfac_v+(4.330127018922193*Bmag[1]*BmagInv[1]*Apar[2]*geoY[2]*dfac_x)/dfac_v-(2.5*BmagInv[0]*Bmag[1]*Apar[2]*geoY[2]*dfac_x)/dfac_v-(8.660254037844386*Apar[2]*geoY[2]*dfac_x)/dfac_v-(3.354101966249685*Apar[1]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x)/dfac_v+(1.936491673103709*Apar[0]*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x)/dfac_v+(1.936491673103709*BmagInv[0]*Apar[1]*Bmag[1]*geoY[2]*dfac_x)/dfac_v-(1.118033988749895*Apar[0]*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x)/dfac_v+(2.23606797749979*Apar[1]*geoY[2]*dfac_x)/dfac_v-(16.77050983124843*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v+(9.682458365518544*geoY[0]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v+(12.99038105676658*Apar[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v-(7.5*Apar[0]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v-(7.5*geoY[0]*Apar[1]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v+(4.330127018922193*Apar[0]*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x)/dfac_v+(4.330127018922193*Bmag[1]*geoY[1]*Apar[2]*BmagInv[2]*dfac_x)/dfac_v-(2.5*geoY[0]*Bmag[1]*Apar[2]*BmagInv[2]*dfac_x)/dfac_v-(3.354101966249685*Apar[1]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x)/dfac_v+(1.936491673103709*Apar[0]*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x)/dfac_v+(1.936491673103709*geoY[0]*Apar[1]*Bmag[1]*BmagInv[2]*dfac_x)/dfac_v-(1.118033988749895*Apar[0]*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x)/dfac_v+(12.99038105676658*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x)/dfac_v-(7.5*BmagInv[0]*geoY[1]*Apar[2]*Bmag[2]*dfac_x)/dfac_v-(7.5*geoY[0]*BmagInv[1]*Apar[2]*Bmag[2]*dfac_x)/dfac_v+(4.330127018922193*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*dfac_x)/dfac_v-(10.06230589874905*Apar[1]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x)/dfac_v+(5.809475019311124*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x)/dfac_v+(5.809475019311124*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*dfac_x)/dfac_v-(3.354101966249685*Apar[0]*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x)/dfac_v+(5.809475019311124*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*dfac_x)/dfac_v-(3.354101966249685*Apar[0]*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x)/dfac_v-(3.354101966249685*BmagInv[0]*geoY[0]*Apar[1]*Bmag[2]*dfac_x)/dfac_v+(1.936491673103709*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x)/dfac_v-(3.354101966249685*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*dfac_x)/dfac_v+(1.936491673103709*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*dfac_x)/dfac_v+(6.708203932499369*geoY[1]*Apar[2]*dfac_x)/dfac_v+(1.936491673103709*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*dfac_x)/dfac_v-(1.118033988749895*BmagInv[0]*geoY[0]*Bmag[1]*Apar[2]*dfac_x)/dfac_v-(3.872983346207417*geoY[0]*Apar[2]*dfac_x)/dfac_v+(2.598076211353315*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x)/dfac_v-(1.5*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x)/dfac_v-(1.5*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*dfac_x)/dfac_v+(0.8660254037844386*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x)/dfac_v-(1.732050807568877*Apar[1]*geoY[1]*dfac_x)/dfac_v-(1.5*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*dfac_x)/dfac_v+(0.8660254037844386*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x)/dfac_v+(0.8660254037844386*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*dfac_x)/dfac_v-(0.5*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x)/dfac_v+(1.0*geoY[0]*Apar[1]*dfac_x)/dfac_v+(1.825741858350554*Gradpar[2])/dfac_v-(1.414213562373095*Gradpar[1])/dfac_v+(0.816496580927726*Gradpar[0])/dfac_v; 
  alpha[5] = (7.071067811865476*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(1.825741858350554*Bmag[1]*BmagInv[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(5.477225575051662*BmagInv[1]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(3.16227766016838*BmagInv[0]*Bmag[2]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)+(1.414213562373095*Bmag[1]*BmagInv[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(0.8164965809277259*BmagInv[0]*Bmag[1]*geoY[2]*dfac_x*m_)/(dfac_v2*q_)-(5.477225575051662*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(3.16227766016838*geoY[0]*Bmag[2]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(1.414213562373095*Bmag[1]*geoY[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)-(0.8164965809277259*geoY[0]*Bmag[1]*BmagInv[2]*dfac_x*m_)/(dfac_v2*q_)+(4.242640687119285*BmagInv[1]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(2.449489742783178*BmagInv[0]*geoY[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(2.449489742783178*geoY[0]*BmagInv[1]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)+(1.414213562373095*BmagInv[0]*geoY[0]*Bmag[2]*dfac_x*m_)/(dfac_v2*q_)-(1.095445115010332*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.6324555320336759*BmagInv[0]*Bmag[1]*geoY[1]*dfac_x*m_)/(dfac_v2*q_)+(0.6324555320336759*geoY[0]*Bmag[1]*BmagInv[1]*dfac_x*m_)/(dfac_v2*q_)-(0.3651483716701108*BmagInv[0]*geoY[0]*Bmag[1]*dfac_x*m_)/(dfac_v2*q_); 
  if (alpha0>0) { 
  incr[0] = 0.01666666666666667*(25.98076211353316*alpha[5]*fl[7]+33.54101966249684*alpha[2]*fl[6]+15.0*alpha[5]*fl[5]+33.54101966249685*alpha[0]*fl[4]+25.98076211353316*alpha[2]*fl[3]+15.0*alpha[2]*fl[2]+25.98076211353316*alpha[0]*fl[1]+15.0*alpha[0]*fl[0])*dfac_x; 
  incr[1] = -0.05*(15.0*alpha[5]*fl[7]+19.36491673103708*alpha[2]*fl[6]+8.660254037844386*alpha[5]*fl[5]+19.36491673103709*alpha[0]*fl[4]+15.0*alpha[2]*fl[3]+8.660254037844386*alpha[2]*fl[2]+15.0*alpha[0]*fl[1]+8.660254037844386*alpha[0]*fl[0])*dfac_x; 
  incr[2] = 0.01666666666666667*(23.2379000772445*alpha[2]*fl[7]+30.0*alpha[5]*fl[6]+33.54101966249684*alpha[0]*fl[6]+13.41640786499874*alpha[2]*fl[5]+23.2379000772445*fl[3]*alpha[5]+13.41640786499874*fl[2]*alpha[5]+33.54101966249685*alpha[2]*fl[4]+25.98076211353316*alpha[0]*fl[3]+15.0*alpha[0]*fl[2]+25.98076211353316*fl[1]*alpha[2]+15.0*fl[0]*alpha[2])*dfac_x; 
  incr[3] = -0.05*(13.41640786499874*alpha[2]*fl[7]+17.32050807568877*alpha[5]*fl[6]+19.36491673103708*alpha[0]*fl[6]+7.745966692414834*alpha[2]*fl[5]+13.41640786499874*fl[3]*alpha[5]+7.745966692414834*fl[2]*alpha[5]+19.36491673103709*alpha[2]*fl[4]+15.0*alpha[0]*fl[3]+8.660254037844386*alpha[0]*fl[2]+15.0*fl[1]*alpha[2]+8.660254037844386*fl[0]*alpha[2])*dfac_x; 
  incr[4] = 0.08333333333333333*(11.61895003862225*alpha[5]*fl[7]+15.0*alpha[2]*fl[6]+6.708203932499369*alpha[5]*fl[5]+15.0*alpha[0]*fl[4]+11.61895003862225*alpha[2]*fl[3]+6.708203932499369*alpha[2]*fl[2]+11.61895003862225*alpha[0]*fl[1]+6.708203932499369*alpha[0]*fl[0])*dfac_x; 
  incr[5] = 0.002380952380952381*(116.1895003862225*alpha[5]*fl[7]+181.8653347947321*alpha[0]*fl[7]+210.0*alpha[2]*fl[6]+67.0820393249937*alpha[5]*fl[5]+105.0*alpha[0]*fl[5]+234.787137637478*fl[4]*alpha[5]+181.8653347947321*fl[1]*alpha[5]+105.0*fl[0]*alpha[5]+162.6653005407115*alpha[2]*fl[3]+93.91485505499116*alpha[2]*fl[2])*dfac_x; 
  incr[6] = 0.01666666666666667*(51.96152422706631*alpha[2]*fl[7]+67.0820393249937*alpha[5]*fl[6]+75.0*alpha[0]*fl[6]+30.0*alpha[2]*fl[5]+51.96152422706632*fl[3]*alpha[5]+30.0*fl[2]*alpha[5]+75.00000000000001*alpha[2]*fl[4]+58.09475019311126*alpha[0]*fl[3]+33.54101966249684*alpha[0]*fl[2]+58.09475019311126*fl[1]*alpha[2]+33.54101966249684*fl[0]*alpha[2])*dfac_x; 
  incr[7] = -0.007142857142857143*(67.0820393249937*alpha[5]*fl[7]+105.0*alpha[0]*fl[7]+121.2435565298214*alpha[2]*fl[6]+38.72983346207417*alpha[5]*fl[5]+60.62177826491071*alpha[0]*fl[5]+135.5544171172596*fl[4]*alpha[5]+105.0*fl[1]*alpha[5]+60.62177826491071*fl[0]*alpha[5]+93.91485505499116*alpha[2]*fl[3]+54.22176684690384*alpha[2]*fl[2])*dfac_x; 

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
  } else { 
  incr[0] = -0.01666666666666667*(25.98076211353316*alpha[5]*fr[7]-33.54101966249684*alpha[2]*fr[6]-15.0*alpha[5]*fr[5]-33.54101966249685*alpha[0]*fr[4]+25.98076211353316*alpha[2]*fr[3]-15.0*alpha[2]*fr[2]+25.98076211353316*alpha[0]*fr[1]-15.0*alpha[0]*fr[0])*dfac_x; 
  incr[1] = 0.05*(15.0*alpha[5]*fr[7]-19.36491673103708*alpha[2]*fr[6]-8.660254037844386*alpha[5]*fr[5]-19.36491673103709*alpha[0]*fr[4]+15.0*alpha[2]*fr[3]-8.660254037844386*alpha[2]*fr[2]+15.0*alpha[0]*fr[1]-8.660254037844386*alpha[0]*fr[0])*dfac_x; 
  incr[2] = -0.01666666666666667*(23.2379000772445*alpha[2]*fr[7]-30.0*alpha[5]*fr[6]-33.54101966249684*alpha[0]*fr[6]-13.41640786499874*alpha[2]*fr[5]+23.2379000772445*fr[3]*alpha[5]-13.41640786499874*fr[2]*alpha[5]-33.54101966249685*alpha[2]*fr[4]+25.98076211353316*alpha[0]*fr[3]-15.0*alpha[0]*fr[2]+25.98076211353316*fr[1]*alpha[2]-15.0*fr[0]*alpha[2])*dfac_x; 
  incr[3] = 0.05*(13.41640786499874*alpha[2]*fr[7]-17.32050807568877*alpha[5]*fr[6]-19.36491673103708*alpha[0]*fr[6]-7.745966692414834*alpha[2]*fr[5]+13.41640786499874*fr[3]*alpha[5]-7.745966692414834*fr[2]*alpha[5]-19.36491673103709*alpha[2]*fr[4]+15.0*alpha[0]*fr[3]-8.660254037844386*alpha[0]*fr[2]+15.0*fr[1]*alpha[2]-8.660254037844386*fr[0]*alpha[2])*dfac_x; 
  incr[4] = -0.08333333333333333*(11.61895003862225*alpha[5]*fr[7]-15.0*alpha[2]*fr[6]-6.708203932499369*alpha[5]*fr[5]-15.0*alpha[0]*fr[4]+11.61895003862225*alpha[2]*fr[3]-6.708203932499369*alpha[2]*fr[2]+11.61895003862225*alpha[0]*fr[1]-6.708203932499369*alpha[0]*fr[0])*dfac_x; 
  incr[5] = -0.002380952380952381*(116.1895003862225*alpha[5]*fr[7]+181.8653347947321*alpha[0]*fr[7]-210.0*alpha[2]*fr[6]-67.0820393249937*alpha[5]*fr[5]-105.0*alpha[0]*fr[5]-234.787137637478*fr[4]*alpha[5]+181.8653347947321*fr[1]*alpha[5]-105.0*fr[0]*alpha[5]+162.6653005407115*alpha[2]*fr[3]-93.91485505499116*alpha[2]*fr[2])*dfac_x; 
  incr[6] = -0.01666666666666667*(51.96152422706631*alpha[2]*fr[7]-67.0820393249937*alpha[5]*fr[6]-75.0*alpha[0]*fr[6]-30.0*alpha[2]*fr[5]+51.96152422706632*fr[3]*alpha[5]-30.0*fr[2]*alpha[5]-75.00000000000001*alpha[2]*fr[4]+58.09475019311126*alpha[0]*fr[3]-33.54101966249684*alpha[0]*fr[2]+58.09475019311126*fr[1]*alpha[2]-33.54101966249684*fr[0]*alpha[2])*dfac_x; 
  incr[7] = 0.007142857142857143*(67.0820393249937*alpha[5]*fr[7]+105.0*alpha[0]*fr[7]-121.2435565298214*alpha[2]*fr[6]-38.72983346207417*alpha[5]*fr[5]-60.62177826491071*alpha[0]*fr[5]-135.5544171172596*fr[4]*alpha[5]+105.0*fr[1]*alpha[5]-60.62177826491071*fr[0]*alpha[5]+93.91485505499116*alpha[2]*fr[3]-54.22176684690384*alpha[2]*fr[2])*dfac_x; 

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
  } 
  return std::abs(alpha0); 
} 
double EmGyrokineticSurf1x1vSer_Vpar_P1_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[4]; 
  double incr[4]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.125*((6.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]+6.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_v*dfac_x*dfac_x*m_*wv+(((4.242640687119286*Apar[0]*Bmag[1]*BmagInv[1]+4.242640687119286*BmagInv[0]*Apar[1]*Bmag[1])*Phi[1]*geoY[1]+(4.242640687119286*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]+4.242640687119286*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]-8.485281374238571*geoY[0]*Apar[1])*Phi[1])*dfac_v*dfac_x*dfac_x-6.928203230275509*Gradpar[0]*Phi[1]*dfac_v*dfac_x-5.656854249492382*dApardt[0]*dfac_v)*q_+((-6.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1])-6.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_x*dfac_x*m_))/(dfac_v*m_); 

  double alpha[4]; 
  alpha[0] = 1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*wv+(1.060660171779821*Apar[0]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*BmagInv[0]*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(2.121320343559643*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(1.414213562373096*dApardt[0]*q_)/m_-(1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[1] = 1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*wv+(1.909188309203679*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_-(2.121320343559643*Apar[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x*q_)/m_-(1.414213562373096*dApardt[1]*q_)/m_-(1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  if (alpha0>0) { 
  incr[0] = 0.25*(1.732050807568877*alpha[1]*fl[3]+1.732050807568877*alpha[0]*fl[2]+alpha[1]*fl[1]+alpha[0]*fl[0])*dfac_v; 
  incr[1] = 0.25*(1.732050807568877*alpha[0]*fl[3]+1.732050807568877*alpha[1]*fl[2]+alpha[0]*fl[1]+fl[0]*alpha[1])*dfac_v; 
  incr[2] = -0.25*(3.0*alpha[1]*fl[3]+3.0*alpha[0]*fl[2]+1.732050807568877*alpha[1]*fl[1]+1.732050807568877*alpha[0]*fl[0])*dfac_v; 
  incr[3] = -0.25*(3.0*alpha[0]*fl[3]+3.0*alpha[1]*fl[2]+1.732050807568877*alpha[0]*fl[1]+1.732050807568877*fl[0]*alpha[1])*dfac_v; 

  outr[0] += incr[0]; 
  outr[1] += incr[1]; 
  outr[2] += incr[2]; 
  outr[3] += incr[3]; 

  outl[0] += -1.0*incr[0]; 
  outl[1] += -1.0*incr[1]; 
  outl[2] += incr[2]; 
  outl[3] += incr[3]; 
  } else { 
  incr[0] = -0.25*(1.732050807568877*alpha[1]*fr[3]+1.732050807568877*alpha[0]*fr[2]-1.0*alpha[1]*fr[1]-1.0*alpha[0]*fr[0])*dfac_v; 
  incr[1] = -0.25*(1.732050807568877*alpha[0]*fr[3]+1.732050807568877*alpha[1]*fr[2]-1.0*alpha[0]*fr[1]-1.0*fr[0]*alpha[1])*dfac_v; 
  incr[2] = 0.25*(3.0*alpha[1]*fr[3]+3.0*alpha[0]*fr[2]-1.732050807568877*alpha[1]*fr[1]-1.732050807568877*alpha[0]*fr[0])*dfac_v; 
  incr[3] = 0.25*(3.0*alpha[0]*fr[3]+3.0*alpha[1]*fr[2]-1.732050807568877*alpha[0]*fr[1]-1.732050807568877*fr[0]*alpha[1])*dfac_v; 

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
double EmGyrokineticSurf1x1vSer_Vpar_P2_Bvars_1(const double q_, const double m_, const double dt, const double *w, const double *dxv, const double amax_in, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *dApardt, const double *fl, const double *fr, double *outl, double *outr) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double Ghat[8]; 
  double incr[8]; 
  // surface-averaged phase velocity in this direction 
  double alpha0 = (0.003571428571428571*((((1650.0*Bmag[2]*BmagInv[2]+939.1485505499119*BmagInv[0]*Bmag[2]+420.0*Bmag[1]*BmagInv[1])*Phi[2]+210.0*Bmag[1]*Phi[1]*BmagInv[2]+420.0*BmagInv[1]*Phi[1]*Bmag[2])*geoY[2]+((939.1485505499119*geoY[0]*Bmag[2]+420.0*Bmag[1]*geoY[1])*BmagInv[2]+(1890.0*BmagInv[1]*geoY[1]+1050.0*BmagInv[0]*geoY[0])*Bmag[2]+469.5742752749559*BmagInv[0]*Bmag[1]*geoY[1]+469.5742752749559*geoY[0]*Bmag[1]*BmagInv[1])*Phi[2]+420.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]+(469.5742752749559*BmagInv[0]*Phi[1]*geoY[1]+469.5742752749559*geoY[0]*BmagInv[1]*Phi[1])*Bmag[2]+210.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]+210.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_v*dfac_x*dfac_x*m_*wv+((((((1897.366596101029*Apar[2]+1166.726188957804*Apar[0])*Bmag[2]+521.7758139277827*Apar[1]*Bmag[1])*BmagInv[2]+(1166.726188957804*BmagInv[0]*Apar[2]+1707.629936490926*Apar[1]*BmagInv[1]+664.0783086353599*Apar[0]*BmagInv[0])*Bmag[2]+(521.7758139277827*Bmag[1]*BmagInv[1]-1328.15661727072)*Apar[2]+296.98484809835*Apar[0]*Bmag[1]*BmagInv[1]+296.98484809835*BmagInv[0]*Apar[1]*Bmag[1])*Phi[2]+(521.7758139277827*Apar[1]*Phi[1]*Bmag[2]+94.86832980505142*Bmag[1]*Phi[1]*Apar[2]+148.492424049175*Apar[0]*Bmag[1]*Phi[1])*BmagInv[2]+(521.7758139277827*BmagInv[1]*Phi[1]*Apar[2]+(296.98484809835*Apar[0]*BmagInv[1]+296.98484809835*BmagInv[0]*Apar[1])*Phi[1])*Bmag[2]+148.492424049175*BmagInv[0]*Bmag[1]*Phi[1]*Apar[2]+132.815661727072*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1])*geoY[2]+(((1166.726188957804*geoY[0]*Apar[2]+1707.629936490926*Apar[1]*geoY[1]+664.0783086353599*Apar[0]*geoY[0])*Bmag[2]+521.7758139277827*Bmag[1]*geoY[1]*Apar[2]+296.98484809835*Apar[0]*Bmag[1]*geoY[1]+296.98484809835*geoY[0]*Apar[1]*Bmag[1])*BmagInv[2]+((1707.629936490926*BmagInv[1]*geoY[1]+664.0783086353599*BmagInv[0]*geoY[0])*Apar[2]+(1336.431816442575*Apar[0]*BmagInv[1]+1336.431816442575*BmagInv[0]*Apar[1])*geoY[1]+1336.431816442575*geoY[0]*Apar[1]*BmagInv[1]+742.462120245875*Apar[0]*BmagInv[0]*geoY[0])*Bmag[2]+(296.98484809835*BmagInv[0]*Bmag[1]*geoY[1]+296.98484809835*geoY[0]*Bmag[1]*BmagInv[1]-1484.92424049175*geoY[0])*Apar[2]+(597.6704777718237*Apar[1]*Bmag[1]*BmagInv[1]+332.0391543176799*Apar[0]*BmagInv[0]*Bmag[1]-664.0783086353599*Apar[1])*geoY[1]+332.0391543176799*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]+332.0391543176799*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1])*Phi[2]+((521.7758139277827*Phi[1]*geoY[1]*Apar[2]+296.98484809835*Apar[0]*Phi[1]*geoY[1]+296.98484809835*geoY[0]*Apar[1]*Phi[1])*Bmag[2]+148.492424049175*geoY[0]*Bmag[1]*Phi[1]*Apar[2]+132.815661727072*Apar[1]*Bmag[1]*Phi[1]*geoY[1])*BmagInv[2]+((296.98484809835*BmagInv[0]*Phi[1]*geoY[1]+296.98484809835*geoY[0]*BmagInv[1]*Phi[1])*Apar[2]+(597.6704777718237*Apar[1]*BmagInv[1]+332.0391543176799*Apar[0]*BmagInv[0])*Phi[1]*geoY[1]+(332.0391543176799*Apar[0]*geoY[0]*BmagInv[1]+332.0391543176799*BmagInv[0]*geoY[0]*Apar[1])*Phi[1])*Bmag[2]+(132.815661727072*Bmag[1]*BmagInv[1]-664.0783086353599)*Phi[1]*geoY[1]*Apar[2]+(148.492424049175*Apar[0]*Bmag[1]*BmagInv[1]+148.492424049175*BmagInv[0]*Apar[1]*Bmag[1])*Phi[1]*geoY[1]+(148.492424049175*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]+148.492424049175*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]-296.98484809835*geoY[0]*Apar[1])*Phi[1])*dfac_v*dfac_x*dfac_x+((-542.2176684690385*Gradpar[1]*Phi[2])-242.4871130596428*Gradpar[0]*Phi[1])*dfac_v*dfac_x-197.9898987322334*dApardt[0]*dfac_v)*q_+((((-1650.0*Bmag[2]*BmagInv[2])-939.1485505499119*BmagInv[0]*Bmag[2]-420.0*Bmag[1]*BmagInv[1])*Phi[2]-210.0*Bmag[1]*Phi[1]*BmagInv[2]-420.0*BmagInv[1]*Phi[1]*Bmag[2])*geoY[2]+(((-939.1485505499119*geoY[0]*Bmag[2])-420.0*Bmag[1]*geoY[1])*BmagInv[2]+((-1890.0*BmagInv[1]*geoY[1])-1050.0*BmagInv[0]*geoY[0])*Bmag[2]-469.5742752749559*BmagInv[0]*Bmag[1]*geoY[1]-469.5742752749559*geoY[0]*Bmag[1]*BmagInv[1])*Phi[2]-420.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]+((-469.5742752749559*BmagInv[0]*Phi[1]*geoY[1])-469.5742752749559*geoY[0]*BmagInv[1]*Phi[1])*Bmag[2]-210.0*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]-210.0*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1])*dfac_x*dfac_x*m_))/(dfac_v*m_); 

  double alpha[8]; 
  alpha[0] = 11.78571428571428*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+6.708203932499371*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+3.0*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+1.5*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+3.0*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+6.708203932499371*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+3.0*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+13.5*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+7.5*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+3.354101966249686*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+3.354101966249686*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*wv+3.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+3.354101966249686*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+3.354101966249686*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*wv+(13.55261854357878*Apar[2]*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(8.33375849255574*Apar[0]*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Apar[1]*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(8.33375849255574*BmagInv[0]*Apar[2]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*Apar[1]*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(4.743416490252571*Apar[0]*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Bmag[1]*BmagInv[1]*Apar[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_-(9.486832980505143*Apar[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*BmagInv[0]*Apar[1]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Apar[1]*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(0.6776309271789387*Bmag[1]*Phi[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*BmagInv[1]*Phi[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*BmagInv[0]*Apar[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*BmagInv[0]*Bmag[1]*Phi[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(0.9486832980505141*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_+(8.33375849255574*geoY[0]*Apar[2]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*Apar[1]*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(4.743416490252571*Apar[0]*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Bmag[1]*geoY[1]*Apar[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*geoY[0]*Apar[1]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(4.743416490252571*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(9.545941546018392*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(9.545941546018392*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(9.545941546018392*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.303300858899107*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_-(10.60660171779822*geoY[0]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(4.269074841227312*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.371708245126285*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_-(4.743416490252571*Apar[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.371708245126285*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.371708245126285*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Phi[1]*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*geoY[0]*Apar[1]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*geoY[0]*Bmag[1]*Phi[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(0.9486832980505141*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*BmagInv[0]*Phi[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*geoY[0]*BmagInv[1]*Phi[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(4.269074841227312*Apar[1]*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(2.371708245126285*Apar[0]*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(2.371708245126285*Apar[0]*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(2.371708245126285*BmagInv[0]*geoY[0]*Apar[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(0.9486832980505141*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*q_)/m_-(4.743416490252571*Phi[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*BmagInv[0]*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(2.121320343559643*geoY[0]*Apar[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(3.872983346207417*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_-(1.414213562373095*dApardt[0]*q_)/m_-(11.78571428571428*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(6.708203932499369*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(3.0*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(1.5*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(3.0*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(6.708203932499369*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(13.5*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(7.5*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249685*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249685*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249685*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249685*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(1.5*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[1] = 5.270731661249505*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+17.24966725499838*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+3.0*BmagInv[0]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+5.270731661249505*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+3.0*BmagInv[0]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+1.341640786499874*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*wv+17.24966725499838*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+3.0*geoY[0]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+13.5*BmagInv[0]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+13.5*geoY[0]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+6.037383539249432*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+3.354101966249686*BmagInv[0]*geoY[0]*Bmag[1]*Phi[2]*dfac_x*dfac_x*wv+3.0*geoY[0]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+1.341640786499874*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*wv+6.037383539249432*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+3.354101966249686*BmagInv[0]*geoY[0]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*wv+(20.45558902718227*Apar[1]*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(6.060915267313267*Bmag[1]*Apar[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Apar[0]*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(20.45558902718227*BmagInv[1]*Apar[2]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*Apar[0]*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*BmagInv[0]*Apar[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*BmagInv[0]*Bmag[1]*Apar[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.45482374058194*Apar[1]*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*BmagInv[0]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_-(4.242640687119287*Apar[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(6.060915267313267*Phi[1]*Apar[2]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Apar[0]*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.666751698511148*Apar[1]*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*BmagInv[0]*Phi[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(5.45482374058194*Apar[1]*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*BmagInv[0]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.666751698511148*Bmag[1]*BmagInv[1]*Phi[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*q_)/m_-(4.242640687119287*Phi[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(0.9486832980505141*Apar[0]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_+(0.9486832980505141*BmagInv[0]*Apar[1]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_+(20.45558902718227*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*Apar[0]*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*geoY[0]*Apar[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*geoY[0]*Bmag[1]*Apar[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.45482374058194*Apar[1]*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*geoY[0]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*BmagInv[0]*geoY[1]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*geoY[0]*BmagInv[1]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(20.45558902718227*Apar[1]*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(9.545941546018392*Apar[0]*BmagInv[0]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(9.545941546018392*Apar[0]*geoY[0]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(9.545941546018392*BmagInv[0]*geoY[0]*Apar[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.45482374058194*Bmag[1]*BmagInv[1]*geoY[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_-(19.09188309203678*geoY[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*BmagInv[0]*geoY[0]*Bmag[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(4.269074841227312*Apar[0]*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(4.269074841227312*BmagInv[0]*Apar[1]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(4.269074841227312*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.371708245126285*Apar[0]*BmagInv[0]*geoY[0]*Bmag[1]*Phi[2]*dfac_x*dfac_x*q_)/m_-(4.743416490252571*geoY[0]*Apar[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*geoY[0]*Phi[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(5.45482374058194*Apar[1]*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*geoY[0]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(1.666751698511148*Bmag[1]*Phi[1]*geoY[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(0.9486832980505141*Apar[0]*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(0.9486832980505141*geoY[0]*Apar[1]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(5.45482374058194*BmagInv[1]*Phi[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*BmagInv[0]*geoY[0]*Phi[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(4.269074841227312*Apar[0]*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(4.269074841227312*BmagInv[0]*Apar[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(4.269074841227312*geoY[0]*Apar[1]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(2.371708245126285*Apar[0]*BmagInv[0]*geoY[0]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(0.9486832980505141*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*q_)/m_+(0.9486832980505141*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*Apar[2]*dfac_x*dfac_x*q_)/m_-(4.743416490252571*geoY[0]*Phi[1]*Apar[2]*dfac_x*dfac_x*q_)/m_+(1.909188309203678*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_-(2.121320343559643*Apar[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(3.464101615137754*Gradpar[2]*Phi[2]*dfac_x*q_)/m_-(3.872983346207417*Gradpar[0]*Phi[2]*dfac_x*q_)/m_-(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x*q_)/m_-(1.414213562373095*dApardt[1]*q_)/m_-(5.270731661249505*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(17.24966725499838*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(3.0*BmagInv[0]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(5.270731661249505*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(3.0*BmagInv[0]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(1.341640786499874*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x)/dfac_v-(17.24966725499838*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*geoY[0]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(13.5*BmagInv[0]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(13.5*geoY[0]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(6.037383539249432*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249686*BmagInv[0]*geoY[0]*Bmag[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*geoY[0]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(1.341640786499874*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(6.037383539249432*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(3.354101966249686*BmagInv[0]*geoY[0]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v-(1.5*geoY[0]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x)/dfac_v; 
  alpha[4] = 19.16629694999821*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+11.78571428571429*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+5.270731661249505*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*wv+0.9583148474999101*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*wv+5.270731661249505*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*wv+1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*wv+11.78571428571429*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+5.270731661249505*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*wv+17.24966725499838*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+6.708203932499371*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*wv+3.0*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*wv+3.0*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*wv+5.270731661249505*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*wv+1.5*geoY[0]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x*wv+3.0*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*wv+3.0*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*wv+1.341640786499874*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*wv+(25.13902355192434*Apar[2]*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(13.55261854357878*Apar[0]*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(6.060915267313267*Apar[1]*Bmag[1]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(13.55261854357878*BmagInv[0]*Apar[2]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(20.45558902718227*Apar[1]*BmagInv[1]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(8.33375849255574*Apar[0]*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(6.060915267313267*Bmag[1]*BmagInv[1]*Apar[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_-(16.66751698511148*Apar[2]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Apar[0]*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*BmagInv[0]*Apar[1]*Bmag[1]*Phi[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(6.060915267313267*Apar[1]*Phi[1]*Bmag[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(2.272843225242474*Bmag[1]*Phi[1]*Apar[2]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(0.6776309271789387*Apar[0]*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(6.060915267313267*BmagInv[1]*Phi[1]*Apar[2]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Apar[0]*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*BmagInv[0]*Apar[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(0.6776309271789387*BmagInv[0]*Bmag[1]*Phi[1]*Apar[2]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.666751698511148*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*BmagInv[0]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_-(2.121320343559643*Apar[1]*Phi[1]*geoY[2]*dfac_x*dfac_x*q_)/m_+(13.55261854357878*geoY[0]*Apar[2]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(20.45558902718227*Apar[1]*geoY[1]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(8.33375849255574*Apar[0]*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(6.060915267313267*Bmag[1]*geoY[1]*Apar[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Apar[0]*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*geoY[0]*Apar[1]*Bmag[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(20.45558902718227*BmagInv[1]*geoY[1]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(8.33375849255574*BmagInv[0]*geoY[0]*Apar[2]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*Apar[0]*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*BmagInv[0]*Apar[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(12.1973566892209*geoY[0]*Apar[1]*BmagInv[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(4.743416490252571*Apar[0]*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*BmagInv[0]*Bmag[1]*geoY[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*geoY[0]*Bmag[1]*BmagInv[1]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_-(9.486832980505143*geoY[0]*Apar[2]*Phi[2]*dfac_x*dfac_x*q_)/m_+(5.45482374058194*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_-(4.242640687119287*Apar[1]*geoY[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*BmagInv[0]*geoY[0]*Apar[1]*Bmag[1]*Phi[2]*dfac_x*dfac_x*q_)/m_+(6.060915267313267*Phi[1]*geoY[1]*Apar[2]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*Apar[0]*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*geoY[0]*Apar[1]*Phi[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(0.6776309271789387*geoY[0]*Bmag[1]*Phi[1]*Apar[2]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(1.666751698511148*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*Apar[0]*geoY[0]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*BmagInv[0]*Phi[1]*geoY[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(3.726970099484163*geoY[0]*BmagInv[1]*Phi[1]*Apar[2]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(5.45482374058194*Apar[1]*BmagInv[1]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*Apar[0]*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(2.121320343559643*BmagInv[0]*geoY[0]*Apar[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x*q_)/m_+(1.666751698511148*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*q_)/m_-(4.242640687119287*Phi[1]*geoY[1]*Apar[2]*dfac_x*dfac_x*q_)/m_+(1.060660171779821*BmagInv[0]*geoY[0]*Bmag[1]*Phi[1]*Apar[2]*dfac_x*dfac_x*q_)/m_+(0.9486832980505142*Apar[0]*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(0.9486832980505142*BmagInv[0]*Apar[1]*Bmag[1]*Phi[1]*geoY[1]*dfac_x*dfac_x*q_)/m_+(0.9486832980505142*geoY[0]*Apar[1]*Bmag[1]*BmagInv[1]*Phi[1]*dfac_x*dfac_x*q_)/m_-(3.464101615137754*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(1.732050807568877*Phi[1]*Gradpar[2]*dfac_x*q_)/m_-(1.414213562373096*dApardt[2]*q_)/m_-(19.1662969499982*Bmag[2]*BmagInv[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(11.78571428571429*BmagInv[0]*Bmag[2]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(5.270731661249505*Bmag[1]*BmagInv[1]*Phi[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(0.9583148474999101*Bmag[1]*Phi[1]*BmagInv[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(5.270731661249505*BmagInv[1]*Phi[1]*Bmag[2]*geoY[2]*dfac_x*dfac_x)/dfac_v-(1.5*BmagInv[0]*Bmag[1]*Phi[1]*geoY[2]*dfac_x*dfac_x)/dfac_v-(11.78571428571429*geoY[0]*Bmag[2]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(5.270731661249505*Bmag[1]*geoY[1]*BmagInv[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(17.24966725499838*BmagInv[1]*geoY[1]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(6.708203932499371*BmagInv[0]*geoY[0]*Bmag[2]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*BmagInv[0]*Bmag[1]*geoY[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(3.0*geoY[0]*Bmag[1]*BmagInv[1]*Phi[2]*dfac_x*dfac_x)/dfac_v-(5.270731661249505*Phi[1]*geoY[1]*Bmag[2]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(1.5*geoY[0]*Bmag[1]*Phi[1]*BmagInv[2]*dfac_x*dfac_x)/dfac_v-(3.0*BmagInv[0]*Phi[1]*geoY[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(3.0*geoY[0]*BmagInv[1]*Phi[1]*Bmag[2]*dfac_x*dfac_x)/dfac_v-(1.341640786499874*Bmag[1]*BmagInv[1]*Phi[1]*geoY[1]*dfac_x*dfac_x)/dfac_v; 
  double amax = 0.0; 
  bool upwind = false; 
  if(upwind) 
    amax = std::abs(alpha0); 
  else 
    amax = amax_in; 

  Ghat[0] = (-1.118033988749895*fr[5]*amax)+1.118033988749895*fl[5]*amax+0.8660254037844386*fr[2]*amax+0.8660254037844386*fl[2]*amax-0.5*fr[0]*amax+0.5*fl[0]*amax+0.5590169943749476*alpha[1]*fr[7]+0.5590169943749476*alpha[1]*fl[7]-0.4330127018922194*alpha[4]*fr[6]+0.4330127018922194*alpha[4]*fl[6]+0.5590169943749475*alpha[0]*fr[5]+0.5590169943749475*alpha[0]*fl[5]+0.25*alpha[4]*fr[4]+0.25*alpha[4]*fl[4]-0.4330127018922193*alpha[1]*fr[3]+0.4330127018922193*alpha[1]*fl[3]-0.4330127018922193*alpha[0]*fr[2]+0.4330127018922193*alpha[0]*fl[2]+0.25*alpha[1]*fr[1]+0.25*alpha[1]*fl[1]+0.25*alpha[0]*fr[0]+0.25*alpha[0]*fl[0]; 
  Ghat[1] = (-1.118033988749895*fr[7]*amax)+1.118033988749895*fl[7]*amax+0.8660254037844386*fr[3]*amax+0.8660254037844386*fl[3]*amax-0.5*fr[1]*amax+0.5*fl[1]*amax+0.5000000000000001*alpha[4]*fr[7]+0.5590169943749476*alpha[0]*fr[7]+0.5000000000000001*alpha[4]*fl[7]+0.5590169943749476*alpha[0]*fl[7]-0.3872983346207417*alpha[1]*fr[6]+0.3872983346207417*alpha[1]*fl[6]+0.5590169943749475*alpha[1]*fr[5]+0.5590169943749475*alpha[1]*fl[5]+0.223606797749979*alpha[1]*fr[4]+0.223606797749979*alpha[1]*fl[4]-0.3872983346207416*fr[3]*alpha[4]+0.3872983346207416*fl[3]*alpha[4]+0.223606797749979*fr[1]*alpha[4]+0.223606797749979*fl[1]*alpha[4]-0.4330127018922193*alpha[0]*fr[3]+0.4330127018922193*alpha[0]*fl[3]-0.4330127018922193*alpha[1]*fr[2]+0.4330127018922193*alpha[1]*fl[2]+0.25*alpha[0]*fr[1]+0.25*alpha[0]*fl[1]+0.25*fr[0]*alpha[1]+0.25*fl[0]*alpha[1]; 
  Ghat[4] = 0.8660254037844386*fr[6]*amax+0.8660254037844386*fl[6]*amax-0.5*fr[4]*amax+0.5*fl[4]*amax+0.5000000000000001*alpha[1]*fr[7]+0.5000000000000001*alpha[1]*fl[7]-0.276641667586244*alpha[4]*fr[6]-0.4330127018922193*alpha[0]*fr[6]+0.276641667586244*alpha[4]*fl[6]+0.4330127018922193*alpha[0]*fl[6]+0.5590169943749475*alpha[4]*fr[5]+0.5590169943749475*alpha[4]*fl[5]+0.159719141249985*alpha[4]*fr[4]+0.25*alpha[0]*fr[4]+0.159719141249985*alpha[4]*fl[4]+0.25*alpha[0]*fl[4]-0.4330127018922193*fr[2]*alpha[4]+0.4330127018922193*fl[2]*alpha[4]+0.25*fr[0]*alpha[4]+0.25*fl[0]*alpha[4]-0.3872983346207416*alpha[1]*fr[3]+0.3872983346207416*alpha[1]*fl[3]+0.223606797749979*alpha[1]*fr[1]+0.223606797749979*alpha[1]*fl[1]; 

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
