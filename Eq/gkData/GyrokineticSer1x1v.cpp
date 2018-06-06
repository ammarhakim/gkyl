#include <GyrokineticModDecl.h> 
double GyrokineticVol1x1vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  out[1] += 1.732050807568877*f[0]*dfac_x*wv; 
  out[2] += -(2.121320343559642*f[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  out[3] += (0.7071067811865475*(2.449489742783178*f[2]*dfac_x*m_*wv-3.0*Phi[1]*f[1]*dfac_v*dfac_x*q_))/m_; 
  double cflFreq = 0.0; 
  cflFreq += fabs(0.7071067811865475*Bmag[0]*wv)*dxInv; 
  cflFreq += fabs(-(0.8660254037844386*Bmag[0]*Phi[1]*dfac_x*q_)/m_)*dvInv; 
  return cflFreq; 
} 
double GyrokineticVol1x1vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  out[1] += (1.732050807568877*f[0]*dfac_v*dfac_x*wv+f[2]*dfac_x)/dfac_v; 
  out[2] += -(2.121320343559642*(2.23606797749979*f[1]*Phi[2]+f[0]*Phi[1])*dfac_v*dfac_x*q_)/m_; 
  out[3] += (0.1*(17.32050807568877*f[2]*dfac_v*dfac_x*m_*wv-21.21320343559643*(2.0*Phi[2]*f[4]+2.23606797749979*f[0]*Phi[2]+Phi[1]*f[1])*dfac_v2*dfac_x*q_+2.0*(4.47213595499958*f[5]+5.0*f[0])*dfac_x*m_))/(dfac_v*m_); 
  out[4] += (3.872983346207417*f[1]*dfac_v*dfac_x*wv+2.23606797749979*f[3]*dfac_x)/dfac_v; 
  out[5] += -(2.121320343559642*(5.0*Phi[2]*f[3]+2.23606797749979*Phi[1]*f[2])*dfac_v*dfac_x*q_)/m_; 
  out[6] += (0.03333333333333333*(116.1895003862225*f[3]*dfac_v*dfac_x*m_*wv-63.63961030678928*(Phi[1]*f[4]+2.0*f[1]*Phi[2])*dfac_v2*dfac_x*q_+10.0*(6.0*f[7]+6.708203932499369*f[1])*dfac_x*m_))/(dfac_v*m_); 
  out[7] += (0.03333333333333333*(51.96152422706632*f[5]*dfac_v*dfac_x*m_*wv-21.21320343559643*(13.41640786499874*Phi[2]*f[6]+6.708203932499369*Phi[1]*f[3]+15.0*Phi[2]*f[2])*dfac_v2*dfac_x*q_+26.83281572999747*f[2]*dfac_x*m_))/(dfac_v*m_); 
  double cflFreq = 0.0; 
  cflFreq += fabs(0.7071067811865475*Bmag[0]*wv)*dxInv; 
  cflFreq += fabs(-(0.8660254037844386*Bmag[0]*Phi[1]*dfac_x*q_)/m_)*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x1vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  out[1] += 1.732050807568877*f[0]*dfac_x*wv; 
  out[2] += -(2.121320343559642*f[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  out[3] += (0.7071067811865475*(2.449489742783178*f[2]*dfac_x*m_*wv-3.0*Phi[1]*f[1]*dfac_v*dfac_x*q_))/m_; 
  double cflFreq = 0.0; 
  cflFreq += fabs(0.7071067811865475*Bmag[0]*wv)*dxInv; 
  cflFreq += fabs(-(0.8660254037844386*Bmag[0]*Phi[1]*dfac_x*q_)/m_)*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x1vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  out[1] += (1.732050807568877*f[0]*dfac_v*dfac_x*wv+f[2]*dfac_x)/dfac_v; 
  out[2] += -(2.121320343559642*(2.23606797749979*f[1]*Phi[2]+f[0]*Phi[1])*dfac_v*dfac_x*q_)/m_; 
  out[3] += (0.1*(17.32050807568877*f[2]*dfac_v*dfac_x*m_*wv-21.21320343559643*(2.0*Phi[2]*f[4]+2.23606797749979*f[0]*Phi[2]+Phi[1]*f[1])*dfac_v2*dfac_x*q_+2.0*(4.47213595499958*f[5]+5.0*f[0])*dfac_x*m_))/(dfac_v*m_); 
  out[4] += (3.872983346207417*f[1]*dfac_v*dfac_x*wv+2.23606797749979*f[3]*dfac_x)/dfac_v; 
  out[5] += -(2.121320343559642*(5.0*Phi[2]*f[3]+2.23606797749979*Phi[1]*f[2])*dfac_v*dfac_x*q_)/m_; 
  out[6] += (0.03333333333333333*(116.1895003862225*f[3]*dfac_v*dfac_x*m_*wv-63.63961030678928*(Phi[1]*f[4]+2.0*f[1]*Phi[2])*dfac_v2*dfac_x*q_+10.0*(6.0*f[7]+6.708203932499369*f[1])*dfac_x*m_))/(dfac_v*m_); 
  out[7] += (0.03333333333333333*(51.96152422706632*f[5]*dfac_v*dfac_x*m_*wv-21.21320343559643*(13.41640786499874*Phi[2]*f[6]+6.708203932499369*Phi[1]*f[3]+15.0*Phi[2]*f[2])*dfac_v2*dfac_x*q_+26.83281572999747*f[2]*dfac_x*m_))/(dfac_v*m_); 
  double cflFreq = 0.0; 
  cflFreq += fabs(0.7071067811865475*Bmag[0]*wv)*dxInv; 
  cflFreq += fabs(-(0.8660254037844386*Bmag[0]*Phi[1]*dfac_x*q_)/m_)*dvInv; 
  return cflFreq; 
} 
double GyrokineticVol1x1vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  out[1] += 1.732050807568877*f[0]*dfac_x*wv; 
  out[2] += -(2.121320343559642*f[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  out[3] += (0.7071067811865475*(2.449489742783178*f[2]*dfac_x*m_*wv-3.0*Phi[1]*f[1]*dfac_v*dfac_x*q_))/m_; 
  double cflFreq = 0.0; 
  cflFreq += fabs(0.7071067811865475*Bmag[0]*wv)*dxInv; 
  cflFreq += fabs(-(0.8660254037844386*Bmag[0]*Phi[1]*dfac_x*q_)/m_)*dvInv; 
  return cflFreq; 
} 
double GyrokineticVol1x1vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  out[1] += (1.732050807568877*f[0]*dfac_v*dfac_x*wv+f[2]*dfac_x)/dfac_v; 
  out[2] += -(2.121320343559642*(2.23606797749979*f[1]*Phi[2]+f[0]*Phi[1])*dfac_v*dfac_x*q_)/m_; 
  out[3] += (0.1*(17.32050807568877*f[2]*dfac_v*dfac_x*m_*wv-21.21320343559643*(2.0*Phi[2]*f[4]+2.23606797749979*f[0]*Phi[2]+Phi[1]*f[1])*dfac_v2*dfac_x*q_+2.0*(4.47213595499958*f[5]+5.0*f[0])*dfac_x*m_))/(dfac_v*m_); 
  out[4] += (3.872983346207417*f[1]*dfac_v*dfac_x*wv+2.23606797749979*f[3]*dfac_x)/dfac_v; 
  out[5] += -(2.121320343559642*(5.0*Phi[2]*f[3]+2.23606797749979*Phi[1]*f[2])*dfac_v*dfac_x*q_)/m_; 
  out[6] += (0.03333333333333333*(116.1895003862225*f[3]*dfac_v*dfac_x*m_*wv-63.63961030678928*(Phi[1]*f[4]+2.0*f[1]*Phi[2])*dfac_v2*dfac_x*q_+10.0*(6.0*f[7]+6.708203932499369*f[1])*dfac_x*m_))/(dfac_v*m_); 
  out[7] += (0.03333333333333333*(51.96152422706632*f[5]*dfac_v*dfac_x*m_*wv-21.21320343559643*(13.41640786499874*Phi[2]*f[6]+6.708203932499369*Phi[1]*f[3]+15.0*Phi[2]*f[2])*dfac_v2*dfac_x*q_+26.83281572999747*f[2]*dfac_x*m_))/(dfac_v*m_); 
  double cflFreq = 0.0; 
  cflFreq += fabs((0.7071067811865475*Bmag[0]-0.7905694150420947*Bmag[2])*wv)*dxInv; 
  cflFreq += fabs(-(1.224744871391589*Phi[1]*(0.7071067811865475*Bmag[0]-0.7905694150420947*Bmag[2])*dfac_x*q_)/m_)*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x1vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  out[1] += 1.732050807568877*f[0]*dfac_x*wv; 
  out[2] += -(2.121320343559642*f[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  out[3] += (0.7071067811865475*(2.449489742783178*f[2]*dfac_x*m_*wv-3.0*Phi[1]*f[1]*dfac_v*dfac_x*q_))/m_; 
  double cflFreq = 0.0; 
  cflFreq += fabs(0.7071067811865475*Bmag[0]*wv)*dxInv; 
  cflFreq += fabs(-(0.8660254037844386*Bmag[0]*Phi[1]*dfac_x*q_)/m_)*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x1vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  out[1] += (1.732050807568877*f[0]*dfac_v*dfac_x*wv+f[2]*dfac_x)/dfac_v; 
  out[2] += -(2.121320343559642*(2.23606797749979*f[1]*Phi[2]+f[0]*Phi[1])*dfac_v*dfac_x*q_)/m_; 
  out[3] += (0.1*(17.32050807568877*f[2]*dfac_v*dfac_x*m_*wv-21.21320343559643*(2.0*Phi[2]*f[4]+2.23606797749979*f[0]*Phi[2]+Phi[1]*f[1])*dfac_v2*dfac_x*q_+2.0*(4.47213595499958*f[5]+5.0*f[0])*dfac_x*m_))/(dfac_v*m_); 
  out[4] += (3.872983346207417*f[1]*dfac_v*dfac_x*wv+2.23606797749979*f[3]*dfac_x)/dfac_v; 
  out[5] += -(2.121320343559642*(5.0*Phi[2]*f[3]+2.23606797749979*Phi[1]*f[2])*dfac_v*dfac_x*q_)/m_; 
  out[6] += (0.03333333333333333*(116.1895003862225*f[3]*dfac_v*dfac_x*m_*wv-63.63961030678928*(Phi[1]*f[4]+2.0*f[1]*Phi[2])*dfac_v2*dfac_x*q_+10.0*(6.0*f[7]+6.708203932499369*f[1])*dfac_x*m_))/(dfac_v*m_); 
  out[7] += (0.03333333333333333*(51.96152422706632*f[5]*dfac_v*dfac_x*m_*wv-21.21320343559643*(13.41640786499874*Phi[2]*f[6]+6.708203932499369*Phi[1]*f[3]+15.0*Phi[2]*f[2])*dfac_v2*dfac_x*q_+26.83281572999747*f[2]*dfac_x*m_))/(dfac_v*m_); 
  double cflFreq = 0.0; 
  cflFreq += fabs((0.7071067811865475*Bmag[0]-0.7905694150420947*Bmag[2])*wv)*dxInv; 
  cflFreq += fabs(-(1.224744871391589*Phi[1]*(0.7071067811865475*Bmag[0]-0.7905694150420947*Bmag[2])*dfac_x*q_)/m_)*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticStep2Vol1x1vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[1]; 
  double dfac_v = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[3] += -(0.2449489742783178*(4.47213595499958*dApardt[1]*f[4]+5.0*dApardt[0]*f[1]+5.0*f[0]*dApardt[1])*dfac_v*q_)/m_; 
  return fabs(-(0.5*Bmag[0]*dApardt[0]*q_)/m_)*dvInv; 
} 
double EmGyrokineticStep2Vol1x1vSerP2(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[1]; 
  double dfac_v = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[2]*f[4]+dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[3] += -(0.2449489742783178*(4.47213595499958*dApardt[1]*f[4]+4.47213595499958*f[1]*dApardt[2]+5.0*dApardt[0]*f[1]+5.0*f[0]*dApardt[1])*dfac_v*q_)/m_; 
  out[5] += -(0.7071067811865475*(3.872983346207417*dApardt[2]*f[6]+3.872983346207417*dApardt[1]*f[3]+3.872983346207417*dApardt[0]*f[2])*dfac_v*q_)/m_; 
  out[6] += -(0.07824607964359516*(10.0*dApardt[2]*f[4]+15.65247584249853*dApardt[0]*f[4]+15.65247584249853*f[0]*dApardt[2]+14.0*dApardt[1]*f[1])*dfac_v*q_)/m_; 
  out[7] += -(0.1414213562373095*(17.32050807568877*dApardt[1]*f[6]+17.32050807568877*dApardt[2]*f[3]+19.36491673103708*dApardt[0]*f[3]+19.36491673103708*dApardt[1]*f[2])*dfac_v*q_)/m_; 
  return fabs((0.125*((-5.0*Bmag[2]*dApardt[2]*q_)+4.47213595499958*Bmag[0]*dApardt[2]*q_+4.47213595499958*dApardt[0]*Bmag[2]*q_-4.0*Bmag[0]*dApardt[0]*q_))/m_)*dvInv; 
} 
