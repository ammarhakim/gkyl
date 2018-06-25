#include <GyrokineticModDecl.h> 
double GyrokineticVol1x1vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphax[4]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*wv; 
  out[1] += 0.8660254037844386*alphax[0]*f[0]*dfac_x; 
  out[3] += 0.8660254037844386*alphax[0]*f[2]*dfac_x; 
  cflFreq += fabs(0.5*alphax[0])*dxInv; 
  double alphav[4]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  out[2] += 0.8660254037844386*alphav[0]*f[0]*dfac_v; 
  out[3] += 0.8660254037844386*alphav[0]*f[1]*dfac_v; 
  cflFreq += fabs(0.5*alphav[0])*dvInv; 
  return cflFreq; 
} 
double GyrokineticVol1x1vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphax[8]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*wv; 
  alphax[2] = (0.8164965809277261*Gradpar[0])/dfac_v; 
  out[1] += 0.8660254037844386*(alphax[2]*f[2]+alphax[0]*f[0])*dfac_x; 
  out[3] += 0.1732050807568877*(4.47213595499958*alphax[2]*f[5]+5.0*alphax[0]*f[2]+5.0*f[0]*alphax[2])*dfac_x; 
  out[4] += 1.936491673103709*(alphax[2]*f[3]+alphax[0]*f[1])*dfac_x; 
  out[6] += 0.5*(3.464101615137754*alphax[2]*f[7]+3.872983346207417*alphax[0]*f[3]+3.872983346207417*f[1]*alphax[2])*dfac_x; 
  out[7] += 0.3872983346207417*(2.23606797749979*alphax[0]*f[5]+2.0*alphax[2]*f[2])*dfac_x; 
  cflFreq += fabs(0.5*alphax[0])*dxInv; 
  double alphav[8]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = -(3.872983346207417*Gradpar[0]*Phi[2]*dfac_x*q_)/m_; 
  out[2] += 0.8660254037844386*(alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[3] += 0.1732050807568877*(4.47213595499958*alphav[1]*f[4]+5.0*alphav[0]*f[1]+5.0*f[0]*alphav[1])*dfac_v; 
  out[5] += 1.936491673103709*(alphav[1]*f[3]+alphav[0]*f[2])*dfac_v; 
  out[6] += 0.3872983346207417*(2.23606797749979*alphav[0]*f[4]+2.0*alphav[1]*f[1])*dfac_v; 
  out[7] += 0.5*(3.464101615137754*alphav[1]*f[6]+3.872983346207417*alphav[0]*f[3]+3.872983346207417*alphav[1]*f[2])*dfac_v; 
  cflFreq += fabs(0.5*alphav[0])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x1vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphax[4]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*wv; 
  out[1] += 0.8660254037844386*alphax[0]*f[0]*dfac_x; 
  out[3] += 0.8660254037844386*alphax[0]*f[2]*dfac_x; 
  cflFreq += fabs(0.5*alphax[0])*dxInv; 
  double alphav[4]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  out[2] += 0.8660254037844386*alphav[0]*f[0]*dfac_v; 
  out[3] += 0.8660254037844386*alphav[0]*f[1]*dfac_v; 
  cflFreq += fabs(0.5*alphav[0])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x1vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphax[8]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*wv; 
  alphax[2] = (0.8164965809277261*Gradpar[0])/dfac_v; 
  out[1] += 0.8660254037844386*(alphax[2]*f[2]+alphax[0]*f[0])*dfac_x; 
  out[3] += 0.1732050807568877*(4.47213595499958*alphax[2]*f[5]+5.0*alphax[0]*f[2]+5.0*f[0]*alphax[2])*dfac_x; 
  out[4] += 1.936491673103709*(alphax[2]*f[3]+alphax[0]*f[1])*dfac_x; 
  out[6] += 0.5*(3.464101615137754*alphax[2]*f[7]+3.872983346207417*alphax[0]*f[3]+3.872983346207417*f[1]*alphax[2])*dfac_x; 
  out[7] += 0.3872983346207417*(2.23606797749979*alphax[0]*f[5]+2.0*alphax[2]*f[2])*dfac_x; 
  cflFreq += fabs(0.5*alphax[0])*dxInv; 
  double alphav[8]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = -(3.872983346207417*Gradpar[0]*Phi[2]*dfac_x*q_)/m_; 
  out[2] += 0.8660254037844386*(alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[3] += 0.1732050807568877*(4.47213595499958*alphav[1]*f[4]+5.0*alphav[0]*f[1]+5.0*f[0]*alphav[1])*dfac_v; 
  out[5] += 1.936491673103709*(alphav[1]*f[3]+alphav[0]*f[2])*dfac_v; 
  out[6] += 0.3872983346207417*(2.23606797749979*alphav[0]*f[4]+2.0*alphav[1]*f[1])*dfac_v; 
  out[7] += 0.5*(3.464101615137754*alphav[1]*f[6]+3.872983346207417*alphav[0]*f[3]+3.872983346207417*alphav[1]*f[2])*dfac_v; 
  cflFreq += fabs(0.5*alphav[0])*dvInv; 
  return cflFreq; 
} 
double GyrokineticVol1x1vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphax[4]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*wv; 
  alphax[1] = 1.414213562373095*Gradpar[1]*wv; 
  out[1] += 0.8660254037844386*(alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[3] += 0.8660254037844386*(alphax[1]*f[3]+alphax[0]*f[2])*dfac_x; 
  cflFreq += fabs(0.5*alphax[0])*dxInv; 
  double alphav[4]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = -(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x*q_)/m_; 
  out[2] += 0.8660254037844386*(alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[3] += 0.8660254037844386*(alphav[0]*f[1]+f[0]*alphav[1])*dfac_v; 
  cflFreq += fabs(0.5*alphav[0])*dvInv; 
  return cflFreq; 
} 
double GyrokineticVol1x1vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphax[8]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*wv; 
  alphax[1] = 1.414213562373095*Gradpar[1]*wv; 
  alphax[2] = (0.8164965809277261*Gradpar[0])/dfac_v; 
  alphax[3] = (0.8164965809277261*Gradpar[1])/dfac_v; 
  alphax[4] = 1.414213562373095*Gradpar[2]*wv; 
  alphax[6] = (0.816496580927726*Gradpar[2])/dfac_v; 
  out[1] += 0.8660254037844386*(alphax[6]*f[6]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[3] += 0.1*(7.745966692414834*alphax[3]*f[7]+8.660254037844387*alphax[4]*f[6]+8.660254037844387*f[4]*alphax[6]+7.745966692414834*alphax[2]*f[5]+8.660254037844386*alphax[1]*f[3]+8.660254037844386*f[1]*alphax[3]+8.660254037844386*alphax[0]*f[2]+8.660254037844386*f[0]*alphax[2])*dfac_x; 
  out[4] += 0.1*(17.32050807568877*alphax[3]*f[6]+17.32050807568877*f[3]*alphax[6]+17.32050807568877*alphax[1]*f[4]+17.32050807568877*f[1]*alphax[4]+19.36491673103709*alphax[2]*f[3]+19.36491673103709*f[2]*alphax[3]+19.36491673103709*alphax[0]*f[1]+19.36491673103709*f[0]*alphax[1])*dfac_x; 
  out[6] += 0.1*(15.49193338482967*alphax[6]*f[7]+17.32050807568877*alphax[2]*f[7]+17.32050807568877*alphax[1]*f[6]+17.32050807568877*f[1]*alphax[6]+17.32050807568877*alphax[3]*f[5]+17.32050807568877*alphax[3]*f[4]+17.32050807568877*f[3]*alphax[4]+19.36491673103708*alphax[0]*f[3]+19.36491673103708*f[0]*alphax[3]+19.36491673103708*alphax[1]*f[2]+19.36491673103708*f[1]*alphax[2])*dfac_x; 
  out[7] += 0.1*(8.660254037844386*alphax[1]*f[7]+7.745966692414834*alphax[6]*f[6]+8.660254037844387*alphax[0]*f[5]+7.745966692414834*alphax[3]*f[3]+7.745966692414834*alphax[2]*f[2])*dfac_x; 
  cflFreq += fabs(0.5*alphax[0]-0.5590169943749475*alphax[4])*dxInv; 
  double alphav[8]; 
  alphav[0] = (-(3.872983346207417*Gradpar[1]*Phi[2]*dfac_x*q_)/m_)-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = (-(3.464101615137754*Gradpar[2]*Phi[2]*dfac_x*q_)/m_)-(3.872983346207417*Gradpar[0]*Phi[2]*dfac_x*q_)/m_-(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x*q_)/m_; 
  alphav[4] = (-(3.464101615137754*Gradpar[1]*Phi[2]*dfac_x*q_)/m_)-(1.732050807568877*Phi[1]*Gradpar[2]*dfac_x*q_)/m_; 
  out[2] += 0.8660254037844386*(alphav[4]*f[4]+alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[3] += 0.1732050807568877*(4.47213595499958*alphav[1]*f[4]+4.47213595499958*f[1]*alphav[4]+5.0*alphav[0]*f[1]+5.0*f[0]*alphav[1])*dfac_v; 
  out[5] += 0.5*(3.872983346207417*alphav[4]*f[6]+3.872983346207417*alphav[1]*f[3]+3.872983346207417*alphav[0]*f[2])*dfac_v; 
  out[6] += 0.05532833351724881*(10.0*alphav[4]*f[4]+15.65247584249853*alphav[0]*f[4]+15.65247584249853*f[0]*alphav[4]+14.0*alphav[1]*f[1])*dfac_v; 
  out[7] += 0.1*(17.32050807568877*alphav[1]*f[6]+17.32050807568877*f[3]*alphav[4]+19.36491673103708*alphav[0]*f[3]+19.36491673103708*alphav[1]*f[2])*dfac_v; 
  cflFreq += fabs(0.5*alphav[0]-0.5590169943749475*alphav[4])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x1vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphax[4]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*wv; 
  alphax[1] = 1.414213562373095*Gradpar[1]*wv; 
  out[1] += 0.8660254037844386*(alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[3] += 0.8660254037844386*(alphax[1]*f[3]+alphax[0]*f[2])*dfac_x; 
  cflFreq += fabs(0.5*alphax[0])*dxInv; 
  double alphav[4]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = -(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x*q_)/m_; 
  out[2] += 0.8660254037844386*(alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[3] += 0.8660254037844386*(alphav[0]*f[1]+f[0]*alphav[1])*dfac_v; 
  cflFreq += fabs(0.5*alphav[0])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x1vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphax[8]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*wv; 
  alphax[1] = 1.414213562373095*Gradpar[1]*wv; 
  alphax[2] = (0.8164965809277261*Gradpar[0])/dfac_v; 
  alphax[3] = (0.8164965809277261*Gradpar[1])/dfac_v; 
  alphax[4] = 1.414213562373095*Gradpar[2]*wv; 
  alphax[6] = (0.816496580927726*Gradpar[2])/dfac_v; 
  out[1] += 0.8660254037844386*(alphax[6]*f[6]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[3] += 0.1*(7.745966692414834*alphax[3]*f[7]+8.660254037844387*alphax[4]*f[6]+8.660254037844387*f[4]*alphax[6]+7.745966692414834*alphax[2]*f[5]+8.660254037844386*alphax[1]*f[3]+8.660254037844386*f[1]*alphax[3]+8.660254037844386*alphax[0]*f[2]+8.660254037844386*f[0]*alphax[2])*dfac_x; 
  out[4] += 0.1*(17.32050807568877*alphax[3]*f[6]+17.32050807568877*f[3]*alphax[6]+17.32050807568877*alphax[1]*f[4]+17.32050807568877*f[1]*alphax[4]+19.36491673103709*alphax[2]*f[3]+19.36491673103709*f[2]*alphax[3]+19.36491673103709*alphax[0]*f[1]+19.36491673103709*f[0]*alphax[1])*dfac_x; 
  out[6] += 0.1*(15.49193338482967*alphax[6]*f[7]+17.32050807568877*alphax[2]*f[7]+17.32050807568877*alphax[1]*f[6]+17.32050807568877*f[1]*alphax[6]+17.32050807568877*alphax[3]*f[5]+17.32050807568877*alphax[3]*f[4]+17.32050807568877*f[3]*alphax[4]+19.36491673103708*alphax[0]*f[3]+19.36491673103708*f[0]*alphax[3]+19.36491673103708*alphax[1]*f[2]+19.36491673103708*f[1]*alphax[2])*dfac_x; 
  out[7] += 0.1*(8.660254037844386*alphax[1]*f[7]+7.745966692414834*alphax[6]*f[6]+8.660254037844387*alphax[0]*f[5]+7.745966692414834*alphax[3]*f[3]+7.745966692414834*alphax[2]*f[2])*dfac_x; 
  cflFreq += fabs(0.5*alphax[0]-0.5590169943749475*alphax[4])*dxInv; 
  double alphav[8]; 
  alphav[0] = (-(3.872983346207417*Gradpar[1]*Phi[2]*dfac_x*q_)/m_)-(1.732050807568877*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = (-(3.464101615137754*Gradpar[2]*Phi[2]*dfac_x*q_)/m_)-(3.872983346207417*Gradpar[0]*Phi[2]*dfac_x*q_)/m_-(1.732050807568877*Gradpar[1]*Phi[1]*dfac_x*q_)/m_; 
  alphav[4] = (-(3.464101615137754*Gradpar[1]*Phi[2]*dfac_x*q_)/m_)-(1.732050807568877*Phi[1]*Gradpar[2]*dfac_x*q_)/m_; 
  out[2] += 0.8660254037844386*(alphav[4]*f[4]+alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[3] += 0.1732050807568877*(4.47213595499958*alphav[1]*f[4]+4.47213595499958*f[1]*alphav[4]+5.0*alphav[0]*f[1]+5.0*f[0]*alphav[1])*dfac_v; 
  out[5] += 0.5*(3.872983346207417*alphav[4]*f[6]+3.872983346207417*alphav[1]*f[3]+3.872983346207417*alphav[0]*f[2])*dfac_v; 
  out[6] += 0.05532833351724881*(10.0*alphav[4]*f[4]+15.65247584249853*alphav[0]*f[4]+15.65247584249853*f[0]*alphav[4]+14.0*alphav[1]*f[1])*dfac_v; 
  out[7] += 0.1*(17.32050807568877*alphav[1]*f[6]+17.32050807568877*f[3]*alphav[4]+19.36491673103708*alphav[0]*f[3]+19.36491673103708*alphav[1]*f[2])*dfac_v; 
  cflFreq += fabs(0.5*alphav[0]-0.5590169943749475*alphav[4])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticStep2Vol1x1vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[1]; 
  double dfac_v = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[3] += -(0.2449489742783178*(4.47213595499958*dApardt[1]*f[4]+5.0*dApardt[0]*f[1]+5.0*f[0]*dApardt[1])*dfac_v*q_)/m_; 
  return fabs(-(0.7071067811865475*dApardt[0]*q_)/m_)*dvInv; 
} 
double EmGyrokineticStep2Vol1x1vSerP2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[1]; 
  double dfac_v = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[2]*f[4]+dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[3] += -(0.2449489742783178*(4.47213595499958*dApardt[1]*f[4]+4.47213595499958*f[1]*dApardt[2]+5.0*dApardt[0]*f[1]+5.0*f[0]*dApardt[1])*dfac_v*q_)/m_; 
  out[5] += -(0.7071067811865475*(3.872983346207417*dApardt[2]*f[6]+3.872983346207417*dApardt[1]*f[3]+3.872983346207417*dApardt[0]*f[2])*dfac_v*q_)/m_; 
  out[6] += -(0.07824607964359516*(10.0*dApardt[2]*f[4]+15.65247584249853*dApardt[0]*f[4]+15.65247584249853*f[0]*dApardt[2]+14.0*dApardt[1]*f[1])*dfac_v*q_)/m_; 
  out[7] += -(0.1414213562373095*(17.32050807568877*dApardt[1]*f[6]+17.32050807568877*dApardt[2]*f[3]+19.36491673103708*dApardt[0]*f[3]+19.36491673103708*dApardt[1]*f[2])*dfac_v*q_)/m_; 
  return fabs((0.3535533905932737*(2.23606797749979*dApardt[2]*q_-2.0*dApardt[0]*q_))/m_)*dvInv; 
} 
