#include <GyrokineticModDecl.h> 
double GyrokineticVol1x2vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dmInv = 1.0/dxv[2]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[8]; 
  alphax[0] = 2.0*Gradpar[0]*wv; 
  out[1] += 0.6123724356957944*alphax[0]*f[0]*dfac_x; 
  out[4] += 0.6123724356957944*alphax[0]*f[2]*dfac_x; 
  out[5] += 0.6123724356957944*alphax[0]*f[3]*dfac_x; 
  out[7] += 0.6123724356957944*alphax[0]*f[6]*dfac_x; 
  cflFreq += fabs(0.3535533905932737*alphax[0])*dxInv; 
  double alphav[8]; 
  alphav[0] = -(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  out[2] += 0.6123724356957944*alphav[0]*f[0]*dfac_v; 
  out[4] += 0.6123724356957944*alphav[0]*f[1]*dfac_v; 
  out[6] += 0.6123724356957944*alphav[0]*f[3]*dfac_v; 
  out[7] += 0.6123724356957944*alphav[0]*f[5]*dfac_v; 
  cflFreq += fabs(0.3535533905932737*alphav[0])*dvInv; 
  return cflFreq; 
} 
double GyrokineticVol1x2vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dmInv = 1.0/dxv[2]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[20]; 
  alphax[0] = 2.0*Gradpar[0]*wv; 
  alphax[2] = (1.154700538379252*Gradpar[0])/dfac_v; 
  out[1] += 0.6123724356957944*(alphax[2]*f[2]+alphax[0]*f[0])*dfac_x; 
  out[4] += 0.1224744871391589*(4.47213595499958*alphax[2]*f[8]+5.0*alphax[0]*f[2]+5.0*f[0]*alphax[2])*dfac_x; 
  out[5] += 0.6123724356957944*(alphax[2]*f[6]+alphax[0]*f[3])*dfac_x; 
  out[7] += 1.369306393762915*(alphax[2]*f[4]+alphax[0]*f[1])*dfac_x; 
  out[10] += 0.07071067811865474*(7.745966692414834*alphax[2]*f[14]+8.660254037844386*alphax[0]*f[6]+8.660254037844386*alphax[2]*f[3])*dfac_x; 
  out[11] += 0.3535533905932737*(3.464101615137754*alphax[2]*f[12]+3.872983346207417*alphax[0]*f[4]+3.872983346207417*f[1]*alphax[2])*dfac_x; 
  out[12] += 0.273861278752583*(2.23606797749979*alphax[0]*f[8]+2.0*alphax[2]*f[2])*dfac_x; 
  out[13] += 1.369306393762915*(alphax[2]*f[10]+alphax[0]*f[5])*dfac_x; 
  out[15] += 0.07071067811865474*(8.660254037844386*alphax[2]*f[16]+8.660254037844387*alphax[0]*f[9])*dfac_x; 
  out[17] += 0.6123724356957944*(2.0*alphax[2]*f[18]+2.23606797749979*alphax[0]*f[10]+2.23606797749979*alphax[2]*f[5])*dfac_x; 
  out[18] += 0.1581138830084189*(3.872983346207417*alphax[0]*f[14]+3.464101615137754*alphax[2]*f[6])*dfac_x; 
  out[19] += 0.07071067811865474*(8.660254037844387*alphax[0]*f[16]+8.660254037844386*alphax[2]*f[9])*dfac_x; 
  cflFreq += fabs(0.3535533905932737*alphax[0])*dxInv; 
  double alphav[20]; 
  alphav[0] = -(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = -(5.477225575051662*Gradpar[0]*Phi[2]*dfac_x*q_)/m_; 
  out[2] += 0.6123724356957944*(alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[4] += 0.1224744871391589*(4.47213595499958*alphav[1]*f[7]+5.0*alphav[0]*f[1]+5.0*f[0]*alphav[1])*dfac_v; 
  out[6] += 0.6123724356957944*(alphav[1]*f[5]+alphav[0]*f[3])*dfac_v; 
  out[8] += 1.369306393762915*(alphav[1]*f[4]+alphav[0]*f[2])*dfac_v; 
  out[10] += 0.07071067811865474*(7.745966692414834*alphav[1]*f[13]+8.660254037844386*alphav[0]*f[5]+8.660254037844386*alphav[1]*f[3])*dfac_v; 
  out[11] += 0.273861278752583*(2.23606797749979*alphav[0]*f[7]+2.0*alphav[1]*f[1])*dfac_v; 
  out[12] += 0.3535533905932737*(3.464101615137754*alphav[1]*f[11]+3.872983346207417*alphav[0]*f[4]+3.872983346207417*alphav[1]*f[2])*dfac_v; 
  out[14] += 1.369306393762915*(alphav[1]*f[10]+alphav[0]*f[6])*dfac_v; 
  out[16] += 0.07071067811865474*(8.660254037844386*alphav[1]*f[15]+8.660254037844387*alphav[0]*f[9])*dfac_v; 
  out[17] += 0.1581138830084189*(3.872983346207417*alphav[0]*f[13]+3.464101615137754*alphav[1]*f[5])*dfac_v; 
  out[18] += 0.6123724356957944*(2.0*alphav[1]*f[17]+2.23606797749979*alphav[0]*f[10]+2.23606797749979*alphav[1]*f[6])*dfac_v; 
  out[19] += 0.07071067811865474*(8.660254037844387*alphav[0]*f[15]+8.660254037844386*alphav[1]*f[9])*dfac_v; 
  cflFreq += fabs(0.3535533905932737*alphav[0])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x2vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dmInv = 1.0/dxv[2]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[8]; 
  alphax[0] = 2.0*Gradpar[0]*wv; 
  out[1] += 0.6123724356957944*alphax[0]*f[0]*dfac_x; 
  out[4] += 0.6123724356957944*alphax[0]*f[2]*dfac_x; 
  out[5] += 0.6123724356957944*alphax[0]*f[3]*dfac_x; 
  out[7] += 0.6123724356957944*alphax[0]*f[6]*dfac_x; 
  cflFreq += fabs(0.3535533905932737*alphax[0])*dxInv; 
  double alphav[8]; 
  alphav[0] = -(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  out[2] += 0.6123724356957944*alphav[0]*f[0]*dfac_v; 
  out[4] += 0.6123724356957944*alphav[0]*f[1]*dfac_v; 
  out[6] += 0.6123724356957944*alphav[0]*f[3]*dfac_v; 
  out[7] += 0.6123724356957944*alphav[0]*f[5]*dfac_v; 
  cflFreq += fabs(0.3535533905932737*alphav[0])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x2vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dmInv = 1.0/dxv[2]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[20]; 
  alphax[0] = 2.0*Gradpar[0]*wv; 
  alphax[2] = (1.154700538379252*Gradpar[0])/dfac_v; 
  out[1] += 0.6123724356957944*(alphax[2]*f[2]+alphax[0]*f[0])*dfac_x; 
  out[4] += 0.1224744871391589*(4.47213595499958*alphax[2]*f[8]+5.0*alphax[0]*f[2]+5.0*f[0]*alphax[2])*dfac_x; 
  out[5] += 0.6123724356957944*(alphax[2]*f[6]+alphax[0]*f[3])*dfac_x; 
  out[7] += 1.369306393762915*(alphax[2]*f[4]+alphax[0]*f[1])*dfac_x; 
  out[10] += 0.07071067811865474*(7.745966692414834*alphax[2]*f[14]+8.660254037844386*alphax[0]*f[6]+8.660254037844386*alphax[2]*f[3])*dfac_x; 
  out[11] += 0.3535533905932737*(3.464101615137754*alphax[2]*f[12]+3.872983346207417*alphax[0]*f[4]+3.872983346207417*f[1]*alphax[2])*dfac_x; 
  out[12] += 0.273861278752583*(2.23606797749979*alphax[0]*f[8]+2.0*alphax[2]*f[2])*dfac_x; 
  out[13] += 1.369306393762915*(alphax[2]*f[10]+alphax[0]*f[5])*dfac_x; 
  out[15] += 0.07071067811865474*(8.660254037844386*alphax[2]*f[16]+8.660254037844387*alphax[0]*f[9])*dfac_x; 
  out[17] += 0.6123724356957944*(2.0*alphax[2]*f[18]+2.23606797749979*alphax[0]*f[10]+2.23606797749979*alphax[2]*f[5])*dfac_x; 
  out[18] += 0.1581138830084189*(3.872983346207417*alphax[0]*f[14]+3.464101615137754*alphax[2]*f[6])*dfac_x; 
  out[19] += 0.07071067811865474*(8.660254037844387*alphax[0]*f[16]+8.660254037844386*alphax[2]*f[9])*dfac_x; 
  cflFreq += fabs(0.3535533905932737*alphax[0])*dxInv; 
  double alphav[20]; 
  alphav[0] = -(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = -(5.477225575051662*Gradpar[0]*Phi[2]*dfac_x*q_)/m_; 
  out[2] += 0.6123724356957944*(alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[4] += 0.1224744871391589*(4.47213595499958*alphav[1]*f[7]+5.0*alphav[0]*f[1]+5.0*f[0]*alphav[1])*dfac_v; 
  out[6] += 0.6123724356957944*(alphav[1]*f[5]+alphav[0]*f[3])*dfac_v; 
  out[8] += 1.369306393762915*(alphav[1]*f[4]+alphav[0]*f[2])*dfac_v; 
  out[10] += 0.07071067811865474*(7.745966692414834*alphav[1]*f[13]+8.660254037844386*alphav[0]*f[5]+8.660254037844386*alphav[1]*f[3])*dfac_v; 
  out[11] += 0.273861278752583*(2.23606797749979*alphav[0]*f[7]+2.0*alphav[1]*f[1])*dfac_v; 
  out[12] += 0.3535533905932737*(3.464101615137754*alphav[1]*f[11]+3.872983346207417*alphav[0]*f[4]+3.872983346207417*alphav[1]*f[2])*dfac_v; 
  out[14] += 1.369306393762915*(alphav[1]*f[10]+alphav[0]*f[6])*dfac_v; 
  out[16] += 0.07071067811865474*(8.660254037844386*alphav[1]*f[15]+8.660254037844387*alphav[0]*f[9])*dfac_v; 
  out[17] += 0.1581138830084189*(3.872983346207417*alphav[0]*f[13]+3.464101615137754*alphav[1]*f[5])*dfac_v; 
  out[18] += 0.6123724356957944*(2.0*alphav[1]*f[17]+2.23606797749979*alphav[0]*f[10]+2.23606797749979*alphav[1]*f[6])*dfac_v; 
  out[19] += 0.07071067811865474*(8.660254037844387*alphav[0]*f[15]+8.660254037844386*alphav[1]*f[9])*dfac_v; 
  cflFreq += fabs(0.3535533905932737*alphav[0])*dvInv; 
  return cflFreq; 
} 
double GyrokineticVol1x2vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dmInv = 1.0/dxv[2]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[8]; 
  alphax[0] = 2.0*Gradpar[0]*wv; 
  alphax[1] = 2.0*Gradpar[1]*wv; 
  out[1] += 0.6123724356957944*(alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[4] += 0.6123724356957944*(alphax[1]*f[4]+alphax[0]*f[2])*dfac_x; 
  out[5] += 0.6123724356957944*(alphax[1]*f[5]+alphax[0]*f[3])*dfac_x; 
  out[7] += 0.6123724356957944*(alphax[1]*f[7]+alphax[0]*f[6])*dfac_x; 
  cflFreq += fabs(0.3535533905932737*alphax[0])*dxInv; 
  double alphav[8]; 
  alphav[0] = (-(2.449489742783178*Gradpar[0]*Bmag[1]*dfac_x*wm)/m_)-(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = (-(2.449489742783178*Bmag[1]*Gradpar[1]*dfac_x*wm)/m_)-(2.449489742783178*Gradpar[1]*Phi[1]*dfac_x*q_)/m_; 
  alphav[3] = -(1.414213562373095*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alphav[5] = -(1.414213562373095*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  out[2] += 0.6123724356957944*(alphav[5]*f[5]+alphav[3]*f[3]+alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[4] += 0.6123724356957944*(alphav[3]*f[5]+f[3]*alphav[5]+alphav[0]*f[1]+f[0]*alphav[1])*dfac_v; 
  out[6] += 0.6123724356957944*(alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[3]+f[0]*alphav[3])*dfac_v; 
  out[7] += 0.6123724356957944*(alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[3]+f[1]*alphav[3])*dfac_v; 
  cflFreq += fabs(0.3535533905932737*alphav[0])*dvInv; 
  return cflFreq; 
} 
double GyrokineticVol1x2vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dmInv = 1.0/dxv[2]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[20]; 
  alphax[0] = 2.0*Gradpar[0]*wv; 
  alphax[1] = 2.0*Gradpar[1]*wv; 
  alphax[2] = (1.154700538379252*Gradpar[0])/dfac_v; 
  alphax[4] = (1.154700538379252*Gradpar[1])/dfac_v; 
  alphax[7] = 2.0*Gradpar[2]*wv; 
  alphax[11] = (1.154700538379251*Gradpar[2])/dfac_v; 
  out[1] += 0.6123724356957944*(alphax[11]*f[11]+alphax[7]*f[7]+alphax[4]*f[4]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[4] += 0.07071067811865474*(7.745966692414834*alphax[4]*f[12]+8.660254037844387*alphax[7]*f[11]+8.660254037844387*f[7]*alphax[11]+7.745966692414834*alphax[2]*f[8]+8.660254037844386*alphax[1]*f[4]+8.660254037844386*f[1]*alphax[4]+8.660254037844386*alphax[0]*f[2]+8.660254037844386*f[0]*alphax[2])*dfac_x; 
  out[5] += 0.07071067811865474*(8.660254037844387*alphax[11]*f[17]+8.660254037844387*alphax[7]*f[13]+8.660254037844386*alphax[4]*f[10]+8.660254037844386*alphax[2]*f[6]+8.660254037844386*alphax[1]*f[5]+8.660254037844386*alphax[0]*f[3])*dfac_x; 
  out[7] += 0.07071067811865474*(17.32050807568877*alphax[4]*f[11]+17.32050807568877*f[4]*alphax[11]+17.32050807568877*alphax[1]*f[7]+17.32050807568877*f[1]*alphax[7]+19.36491673103709*alphax[2]*f[4]+19.36491673103709*f[2]*alphax[4]+19.36491673103709*alphax[0]*f[1]+19.36491673103709*f[0]*alphax[1])*dfac_x; 
  out[10] += 0.07071067811865474*(7.745966692414834*alphax[4]*f[18]+8.660254037844386*alphax[7]*f[17]+7.745966692414834*alphax[2]*f[14]+8.660254037844386*alphax[11]*f[13]+8.660254037844386*alphax[1]*f[10]+8.660254037844386*alphax[0]*f[6]+8.660254037844386*alphax[4]*f[5]+8.660254037844386*alphax[2]*f[3])*dfac_x; 
  out[11] += 0.07071067811865474*(15.49193338482967*alphax[11]*f[12]+17.32050807568877*alphax[2]*f[12]+17.32050807568877*alphax[1]*f[11]+17.32050807568877*f[1]*alphax[11]+17.32050807568877*alphax[4]*f[8]+17.32050807568877*alphax[4]*f[7]+17.32050807568877*f[4]*alphax[7]+19.36491673103708*alphax[0]*f[4]+19.36491673103708*f[0]*alphax[4]+19.36491673103708*alphax[1]*f[2]+19.36491673103708*f[1]*alphax[2])*dfac_x; 
  out[12] += 0.07071067811865474*(8.660254037844386*alphax[1]*f[12]+7.745966692414834*alphax[11]*f[11]+8.660254037844387*alphax[0]*f[8]+7.745966692414834*alphax[4]*f[4]+7.745966692414834*alphax[2]*f[2])*dfac_x; 
  out[13] += 0.07071067811865474*(17.32050807568877*alphax[4]*f[17]+17.32050807568877*alphax[1]*f[13]+17.32050807568877*f[10]*alphax[11]+19.36491673103708*alphax[2]*f[10]+17.32050807568877*f[5]*alphax[7]+19.36491673103708*alphax[4]*f[6]+19.36491673103708*alphax[0]*f[5]+19.36491673103708*alphax[1]*f[3])*dfac_x; 
  out[15] += 0.07071067811865474*(8.660254037844387*alphax[4]*f[19]+8.660254037844386*alphax[2]*f[16]+8.660254037844386*alphax[1]*f[15]+8.660254037844387*alphax[0]*f[9])*dfac_x; 
  out[17] += 0.07071067811865474*(15.49193338482967*alphax[11]*f[18]+17.32050807568877*alphax[2]*f[18]+17.32050807568877*alphax[1]*f[17]+17.32050807568877*alphax[4]*f[14]+17.32050807568877*alphax[4]*f[13]+17.32050807568877*f[5]*alphax[11]+17.32050807568877*alphax[7]*f[10]+19.36491673103709*alphax[0]*f[10]+19.36491673103709*alphax[1]*f[6]+19.36491673103709*alphax[2]*f[5]+19.36491673103709*f[3]*alphax[4])*dfac_x; 
  out[18] += 0.07071067811865474*(8.660254037844386*alphax[1]*f[18]+7.745966692414834*alphax[11]*f[17]+8.660254037844387*alphax[0]*f[14]+7.745966692414834*alphax[4]*f[10]+7.745966692414834*alphax[2]*f[6])*dfac_x; 
  out[19] += 0.07071067811865474*(8.660254037844386*alphax[1]*f[19]+8.660254037844387*alphax[0]*f[16]+8.660254037844387*alphax[4]*f[15]+8.660254037844386*alphax[2]*f[9])*dfac_x; 
  cflFreq += fabs(0.3535533905932737*alphax[0]-0.3952847075210473*alphax[7])*dxInv; 
  double alphav[20]; 
  alphav[0] = (-(5.477225575051662*Gradpar[1]*Bmag[2]*dfac_x*wm)/m_)-(2.449489742783178*Gradpar[0]*Bmag[1]*dfac_x*wm)/m_-(5.477225575051662*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = (-(4.898979485566357*Bmag[2]*Gradpar[2]*dfac_x*wm)/m_)-(5.477225575051662*Gradpar[0]*Bmag[2]*dfac_x*wm)/m_-(2.449489742783178*Bmag[1]*Gradpar[1]*dfac_x*wm)/m_-(4.898979485566357*Gradpar[2]*Phi[2]*dfac_x*q_)/m_-(5.477225575051662*Gradpar[0]*Phi[2]*dfac_x*q_)/m_-(2.449489742783178*Gradpar[1]*Phi[1]*dfac_x*q_)/m_; 
  alphav[3] = (-(3.16227766016838*Gradpar[1]*Bmag[2]*dfac_x)/(dfac_m*m_))-(1.414213562373095*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alphav[5] = (-(2.82842712474619*Bmag[2]*Gradpar[2]*dfac_x)/(dfac_m*m_))-(3.16227766016838*Gradpar[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(1.414213562373095*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  alphav[7] = (-(2.449489742783178*Bmag[1]*Gradpar[2]*dfac_x*wm)/m_)-(4.898979485566357*Gradpar[1]*Bmag[2]*dfac_x*wm)/m_-(4.898979485566357*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(2.449489742783178*Phi[1]*Gradpar[2]*dfac_x*q_)/m_; 
  alphav[13] = (-(1.414213562373095*Bmag[1]*Gradpar[2]*dfac_x)/(dfac_m*m_))-(2.82842712474619*Gradpar[1]*Bmag[2]*dfac_x)/(dfac_m*m_); 
  out[2] += 0.6123724356957944*(alphav[13]*f[13]+alphav[7]*f[7]+alphav[5]*f[5]+alphav[3]*f[3]+alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[4] += 0.07071067811865474*(7.745966692414834*alphav[5]*f[13]+7.745966692414834*f[5]*alphav[13]+7.745966692414834*alphav[1]*f[7]+7.745966692414834*f[1]*alphav[7]+8.660254037844386*alphav[3]*f[5]+8.660254037844386*f[3]*alphav[5]+8.660254037844386*alphav[0]*f[1]+8.660254037844386*f[0]*alphav[1])*dfac_v; 
  out[6] += 0.07071067811865474*(7.745966692414834*alphav[5]*f[15]+8.660254037844387*alphav[7]*f[13]+8.660254037844387*f[7]*alphav[13]+7.745966692414834*alphav[3]*f[9]+8.660254037844386*alphav[1]*f[5]+8.660254037844386*f[1]*alphav[5]+8.660254037844386*alphav[0]*f[3]+8.660254037844386*f[0]*alphav[3])*dfac_v; 
  out[8] += 0.3535533905932737*(3.872983346207417*alphav[13]*f[17]+3.872983346207417*alphav[7]*f[11]+3.872983346207417*alphav[5]*f[10]+3.872983346207417*alphav[3]*f[6]+3.872983346207417*alphav[1]*f[4]+3.872983346207417*alphav[0]*f[2])*dfac_v; 
  out[10] += 0.07071067811865474*(6.928203230275509*alphav[13]*f[15]+7.745966692414834*alphav[3]*f[15]+7.745966692414834*alphav[1]*f[13]+7.745966692414834*f[1]*alphav[13]+7.745966692414834*alphav[5]*f[9]+7.745966692414834*alphav[5]*f[7]+7.745966692414834*f[5]*alphav[7]+8.660254037844386*alphav[0]*f[5]+8.660254037844386*f[0]*alphav[5]+8.660254037844386*alphav[1]*f[3]+8.660254037844386*f[1]*alphav[3])*dfac_v; 
  out[11] += 0.01010152544552211*(38.72983346207417*alphav[13]*f[13]+60.6217782649107*alphav[3]*f[13]+60.6217782649107*f[3]*alphav[13]+38.72983346207417*alphav[7]*f[7]+60.62177826491071*alphav[0]*f[7]+60.62177826491071*f[0]*alphav[7]+54.22176684690384*alphav[5]*f[5]+54.22176684690384*alphav[1]*f[1])*dfac_v; 
  out[12] += 0.07071067811865474*(17.32050807568877*alphav[5]*f[17]+17.32050807568877*f[10]*alphav[13]+17.32050807568877*alphav[1]*f[11]+19.36491673103708*alphav[3]*f[10]+17.32050807568877*f[4]*alphav[7]+19.36491673103708*alphav[5]*f[6]+19.36491673103708*alphav[0]*f[4]+19.36491673103708*alphav[1]*f[2])*dfac_v; 
  out[14] += 0.07071067811865474*(17.32050807568877*alphav[5]*f[19]+19.36491673103708*alphav[7]*f[17]+17.32050807568877*alphav[3]*f[16]+19.36491673103708*f[11]*alphav[13]+19.36491673103708*alphav[1]*f[10]+19.36491673103708*alphav[0]*f[6]+19.36491673103708*f[4]*alphav[5]+19.36491673103708*f[2]*alphav[3])*dfac_v; 
  out[16] += 0.07071067811865474*(8.660254037844386*alphav[1]*f[15]+7.745966692414834*alphav[13]*f[13]+8.660254037844387*alphav[0]*f[9]+7.745966692414834*alphav[5]*f[5]+7.745966692414834*alphav[3]*f[3])*dfac_v; 
  out[17] += 0.002020305089104421*(242.4871130596428*alphav[5]*f[15]+193.6491673103708*alphav[7]*f[13]+303.1088913245536*alphav[0]*f[13]+271.1088342345192*f[9]*alphav[13]+193.6491673103708*f[7]*alphav[13]+303.1088913245536*f[0]*alphav[13]+303.1088913245535*alphav[3]*f[7]+303.1088913245535*f[3]*alphav[7]+271.1088342345192*alphav[1]*f[5]+271.1088342345192*f[1]*alphav[5])*dfac_v; 
  out[18] += 0.07071067811865474*(15.49193338482967*alphav[13]*f[19]+17.32050807568877*alphav[3]*f[19]+17.32050807568877*alphav[1]*f[17]+17.32050807568877*alphav[5]*f[16]+17.32050807568877*f[4]*alphav[13]+17.32050807568877*alphav[5]*f[11]+17.32050807568877*alphav[7]*f[10]+19.36491673103709*alphav[0]*f[10]+19.36491673103709*alphav[1]*f[6]+19.36491673103709*f[2]*alphav[5]+19.36491673103709*alphav[3]*f[4])*dfac_v; 
  out[19] += 0.01414213562373095*(38.72983346207417*alphav[7]*f[15]+43.30127018922195*alphav[0]*f[15]+34.64101615137755*alphav[5]*f[13]+34.64101615137755*f[5]*alphav[13]+43.30127018922193*alphav[1]*f[9]+38.72983346207418*alphav[3]*f[5]+38.72983346207418*f[3]*alphav[5])*dfac_v; 
  cflFreq += fabs(0.3535533905932737*alphav[0]-0.3952847075210473*alphav[7])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x2vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dmInv = 1.0/dxv[2]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[8]; 
  alphax[0] = 2.0*Gradpar[0]*wv; 
  alphax[1] = 2.0*Gradpar[1]*wv; 
  out[1] += 0.6123724356957944*(alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[4] += 0.6123724356957944*(alphax[1]*f[4]+alphax[0]*f[2])*dfac_x; 
  out[5] += 0.6123724356957944*(alphax[1]*f[5]+alphax[0]*f[3])*dfac_x; 
  out[7] += 0.6123724356957944*(alphax[1]*f[7]+alphax[0]*f[6])*dfac_x; 
  cflFreq += fabs(0.3535533905932737*alphax[0])*dxInv; 
  double alphav[8]; 
  alphav[0] = (-(2.449489742783178*Gradpar[0]*Bmag[1]*dfac_x*wm)/m_)-(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = (-(2.449489742783178*Bmag[1]*Gradpar[1]*dfac_x*wm)/m_)-(2.449489742783178*Gradpar[1]*Phi[1]*dfac_x*q_)/m_; 
  alphav[3] = -(1.414213562373095*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alphav[5] = -(1.414213562373095*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  out[2] += 0.6123724356957944*(alphav[5]*f[5]+alphav[3]*f[3]+alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[4] += 0.6123724356957944*(alphav[3]*f[5]+f[3]*alphav[5]+alphav[0]*f[1]+f[0]*alphav[1])*dfac_v; 
  out[6] += 0.6123724356957944*(alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[3]+f[0]*alphav[3])*dfac_v; 
  out[7] += 0.6123724356957944*(alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[3]+f[1]*alphav[3])*dfac_v; 
  cflFreq += fabs(0.3535533905932737*alphav[0])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol1x2vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dvInv = 1.0/dxv[1]; 
  double dmInv = 1.0/dxv[2]; 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[20]; 
  alphax[0] = 2.0*Gradpar[0]*wv; 
  alphax[1] = 2.0*Gradpar[1]*wv; 
  alphax[2] = (1.154700538379252*Gradpar[0])/dfac_v; 
  alphax[4] = (1.154700538379252*Gradpar[1])/dfac_v; 
  alphax[7] = 2.0*Gradpar[2]*wv; 
  alphax[11] = (1.154700538379251*Gradpar[2])/dfac_v; 
  out[1] += 0.6123724356957944*(alphax[11]*f[11]+alphax[7]*f[7]+alphax[4]*f[4]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[4] += 0.07071067811865474*(7.745966692414834*alphax[4]*f[12]+8.660254037844387*alphax[7]*f[11]+8.660254037844387*f[7]*alphax[11]+7.745966692414834*alphax[2]*f[8]+8.660254037844386*alphax[1]*f[4]+8.660254037844386*f[1]*alphax[4]+8.660254037844386*alphax[0]*f[2]+8.660254037844386*f[0]*alphax[2])*dfac_x; 
  out[5] += 0.07071067811865474*(8.660254037844387*alphax[11]*f[17]+8.660254037844387*alphax[7]*f[13]+8.660254037844386*alphax[4]*f[10]+8.660254037844386*alphax[2]*f[6]+8.660254037844386*alphax[1]*f[5]+8.660254037844386*alphax[0]*f[3])*dfac_x; 
  out[7] += 0.07071067811865474*(17.32050807568877*alphax[4]*f[11]+17.32050807568877*f[4]*alphax[11]+17.32050807568877*alphax[1]*f[7]+17.32050807568877*f[1]*alphax[7]+19.36491673103709*alphax[2]*f[4]+19.36491673103709*f[2]*alphax[4]+19.36491673103709*alphax[0]*f[1]+19.36491673103709*f[0]*alphax[1])*dfac_x; 
  out[10] += 0.07071067811865474*(7.745966692414834*alphax[4]*f[18]+8.660254037844386*alphax[7]*f[17]+7.745966692414834*alphax[2]*f[14]+8.660254037844386*alphax[11]*f[13]+8.660254037844386*alphax[1]*f[10]+8.660254037844386*alphax[0]*f[6]+8.660254037844386*alphax[4]*f[5]+8.660254037844386*alphax[2]*f[3])*dfac_x; 
  out[11] += 0.07071067811865474*(15.49193338482967*alphax[11]*f[12]+17.32050807568877*alphax[2]*f[12]+17.32050807568877*alphax[1]*f[11]+17.32050807568877*f[1]*alphax[11]+17.32050807568877*alphax[4]*f[8]+17.32050807568877*alphax[4]*f[7]+17.32050807568877*f[4]*alphax[7]+19.36491673103708*alphax[0]*f[4]+19.36491673103708*f[0]*alphax[4]+19.36491673103708*alphax[1]*f[2]+19.36491673103708*f[1]*alphax[2])*dfac_x; 
  out[12] += 0.07071067811865474*(8.660254037844386*alphax[1]*f[12]+7.745966692414834*alphax[11]*f[11]+8.660254037844387*alphax[0]*f[8]+7.745966692414834*alphax[4]*f[4]+7.745966692414834*alphax[2]*f[2])*dfac_x; 
  out[13] += 0.07071067811865474*(17.32050807568877*alphax[4]*f[17]+17.32050807568877*alphax[1]*f[13]+17.32050807568877*f[10]*alphax[11]+19.36491673103708*alphax[2]*f[10]+17.32050807568877*f[5]*alphax[7]+19.36491673103708*alphax[4]*f[6]+19.36491673103708*alphax[0]*f[5]+19.36491673103708*alphax[1]*f[3])*dfac_x; 
  out[15] += 0.07071067811865474*(8.660254037844387*alphax[4]*f[19]+8.660254037844386*alphax[2]*f[16]+8.660254037844386*alphax[1]*f[15]+8.660254037844387*alphax[0]*f[9])*dfac_x; 
  out[17] += 0.07071067811865474*(15.49193338482967*alphax[11]*f[18]+17.32050807568877*alphax[2]*f[18]+17.32050807568877*alphax[1]*f[17]+17.32050807568877*alphax[4]*f[14]+17.32050807568877*alphax[4]*f[13]+17.32050807568877*f[5]*alphax[11]+17.32050807568877*alphax[7]*f[10]+19.36491673103709*alphax[0]*f[10]+19.36491673103709*alphax[1]*f[6]+19.36491673103709*alphax[2]*f[5]+19.36491673103709*f[3]*alphax[4])*dfac_x; 
  out[18] += 0.07071067811865474*(8.660254037844386*alphax[1]*f[18]+7.745966692414834*alphax[11]*f[17]+8.660254037844387*alphax[0]*f[14]+7.745966692414834*alphax[4]*f[10]+7.745966692414834*alphax[2]*f[6])*dfac_x; 
  out[19] += 0.07071067811865474*(8.660254037844386*alphax[1]*f[19]+8.660254037844387*alphax[0]*f[16]+8.660254037844387*alphax[4]*f[15]+8.660254037844386*alphax[2]*f[9])*dfac_x; 
  cflFreq += fabs(0.3535533905932737*alphax[0]-0.3952847075210473*alphax[7])*dxInv; 
  double alphav[20]; 
  alphav[0] = (-(5.477225575051662*Gradpar[1]*Bmag[2]*dfac_x*wm)/m_)-(2.449489742783178*Gradpar[0]*Bmag[1]*dfac_x*wm)/m_-(5.477225575051662*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(2.449489742783178*Gradpar[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = (-(4.898979485566357*Bmag[2]*Gradpar[2]*dfac_x*wm)/m_)-(5.477225575051662*Gradpar[0]*Bmag[2]*dfac_x*wm)/m_-(2.449489742783178*Bmag[1]*Gradpar[1]*dfac_x*wm)/m_-(4.898979485566357*Gradpar[2]*Phi[2]*dfac_x*q_)/m_-(5.477225575051662*Gradpar[0]*Phi[2]*dfac_x*q_)/m_-(2.449489742783178*Gradpar[1]*Phi[1]*dfac_x*q_)/m_; 
  alphav[3] = (-(3.16227766016838*Gradpar[1]*Bmag[2]*dfac_x)/(dfac_m*m_))-(1.414213562373095*Gradpar[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alphav[5] = (-(2.82842712474619*Bmag[2]*Gradpar[2]*dfac_x)/(dfac_m*m_))-(3.16227766016838*Gradpar[0]*Bmag[2]*dfac_x)/(dfac_m*m_)-(1.414213562373095*Bmag[1]*Gradpar[1]*dfac_x)/(dfac_m*m_); 
  alphav[7] = (-(2.449489742783178*Bmag[1]*Gradpar[2]*dfac_x*wm)/m_)-(4.898979485566357*Gradpar[1]*Bmag[2]*dfac_x*wm)/m_-(4.898979485566357*Gradpar[1]*Phi[2]*dfac_x*q_)/m_-(2.449489742783178*Phi[1]*Gradpar[2]*dfac_x*q_)/m_; 
  alphav[13] = (-(1.414213562373095*Bmag[1]*Gradpar[2]*dfac_x)/(dfac_m*m_))-(2.82842712474619*Gradpar[1]*Bmag[2]*dfac_x)/(dfac_m*m_); 
  out[2] += 0.6123724356957944*(alphav[13]*f[13]+alphav[7]*f[7]+alphav[5]*f[5]+alphav[3]*f[3]+alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[4] += 0.07071067811865474*(7.745966692414834*alphav[5]*f[13]+7.745966692414834*f[5]*alphav[13]+7.745966692414834*alphav[1]*f[7]+7.745966692414834*f[1]*alphav[7]+8.660254037844386*alphav[3]*f[5]+8.660254037844386*f[3]*alphav[5]+8.660254037844386*alphav[0]*f[1]+8.660254037844386*f[0]*alphav[1])*dfac_v; 
  out[6] += 0.07071067811865474*(7.745966692414834*alphav[5]*f[15]+8.660254037844387*alphav[7]*f[13]+8.660254037844387*f[7]*alphav[13]+7.745966692414834*alphav[3]*f[9]+8.660254037844386*alphav[1]*f[5]+8.660254037844386*f[1]*alphav[5]+8.660254037844386*alphav[0]*f[3]+8.660254037844386*f[0]*alphav[3])*dfac_v; 
  out[8] += 0.3535533905932737*(3.872983346207417*alphav[13]*f[17]+3.872983346207417*alphav[7]*f[11]+3.872983346207417*alphav[5]*f[10]+3.872983346207417*alphav[3]*f[6]+3.872983346207417*alphav[1]*f[4]+3.872983346207417*alphav[0]*f[2])*dfac_v; 
  out[10] += 0.07071067811865474*(6.928203230275509*alphav[13]*f[15]+7.745966692414834*alphav[3]*f[15]+7.745966692414834*alphav[1]*f[13]+7.745966692414834*f[1]*alphav[13]+7.745966692414834*alphav[5]*f[9]+7.745966692414834*alphav[5]*f[7]+7.745966692414834*f[5]*alphav[7]+8.660254037844386*alphav[0]*f[5]+8.660254037844386*f[0]*alphav[5]+8.660254037844386*alphav[1]*f[3]+8.660254037844386*f[1]*alphav[3])*dfac_v; 
  out[11] += 0.01010152544552211*(38.72983346207417*alphav[13]*f[13]+60.6217782649107*alphav[3]*f[13]+60.6217782649107*f[3]*alphav[13]+38.72983346207417*alphav[7]*f[7]+60.62177826491071*alphav[0]*f[7]+60.62177826491071*f[0]*alphav[7]+54.22176684690384*alphav[5]*f[5]+54.22176684690384*alphav[1]*f[1])*dfac_v; 
  out[12] += 0.07071067811865474*(17.32050807568877*alphav[5]*f[17]+17.32050807568877*f[10]*alphav[13]+17.32050807568877*alphav[1]*f[11]+19.36491673103708*alphav[3]*f[10]+17.32050807568877*f[4]*alphav[7]+19.36491673103708*alphav[5]*f[6]+19.36491673103708*alphav[0]*f[4]+19.36491673103708*alphav[1]*f[2])*dfac_v; 
  out[14] += 0.07071067811865474*(17.32050807568877*alphav[5]*f[19]+19.36491673103708*alphav[7]*f[17]+17.32050807568877*alphav[3]*f[16]+19.36491673103708*f[11]*alphav[13]+19.36491673103708*alphav[1]*f[10]+19.36491673103708*alphav[0]*f[6]+19.36491673103708*f[4]*alphav[5]+19.36491673103708*f[2]*alphav[3])*dfac_v; 
  out[16] += 0.07071067811865474*(8.660254037844386*alphav[1]*f[15]+7.745966692414834*alphav[13]*f[13]+8.660254037844387*alphav[0]*f[9]+7.745966692414834*alphav[5]*f[5]+7.745966692414834*alphav[3]*f[3])*dfac_v; 
  out[17] += 0.002020305089104421*(242.4871130596428*alphav[5]*f[15]+193.6491673103708*alphav[7]*f[13]+303.1088913245536*alphav[0]*f[13]+271.1088342345192*f[9]*alphav[13]+193.6491673103708*f[7]*alphav[13]+303.1088913245536*f[0]*alphav[13]+303.1088913245535*alphav[3]*f[7]+303.1088913245535*f[3]*alphav[7]+271.1088342345192*alphav[1]*f[5]+271.1088342345192*f[1]*alphav[5])*dfac_v; 
  out[18] += 0.07071067811865474*(15.49193338482967*alphav[13]*f[19]+17.32050807568877*alphav[3]*f[19]+17.32050807568877*alphav[1]*f[17]+17.32050807568877*alphav[5]*f[16]+17.32050807568877*f[4]*alphav[13]+17.32050807568877*alphav[5]*f[11]+17.32050807568877*alphav[7]*f[10]+19.36491673103709*alphav[0]*f[10]+19.36491673103709*alphav[1]*f[6]+19.36491673103709*f[2]*alphav[5]+19.36491673103709*alphav[3]*f[4])*dfac_v; 
  out[19] += 0.01414213562373095*(38.72983346207417*alphav[7]*f[15]+43.30127018922195*alphav[0]*f[15]+34.64101615137755*alphav[5]*f[13]+34.64101615137755*f[5]*alphav[13]+43.30127018922193*alphav[1]*f[9]+38.72983346207418*alphav[3]*f[5]+38.72983346207418*f[3]*alphav[5])*dfac_v; 
  cflFreq += fabs(0.3535533905932737*alphav[0]-0.3952847075210473*alphav[7])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticStep2Vol1x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[1]; 
  double dfac_v = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[4] += -(0.2449489742783178*(4.47213595499958*dApardt[1]*f[7]+5.0*dApardt[0]*f[1]+5.0*f[0]*dApardt[1])*dfac_v*q_)/m_; 
  out[6] += -(1.224744871391589*(dApardt[1]*f[5]+dApardt[0]*f[3])*dfac_v*q_)/m_; 
  out[7] += -(0.1414213562373095*(7.745966692414834*dApardt[1]*f[13]+8.660254037844386*dApardt[0]*f[5]+8.660254037844386*dApardt[1]*f[3])*dfac_v*q_)/m_; 
  return fabs(-(0.7071067811865475*dApardt[0]*q_)/m_)*dvInv; 
} 
double EmGyrokineticStep2Vol1x2vSerP2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[1]; 
  double dfac_v = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[2]*f[7]+dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[4] += -(0.2449489742783178*(4.47213595499958*dApardt[1]*f[7]+4.47213595499958*f[1]*dApardt[2]+5.0*dApardt[0]*f[1]+5.0*f[0]*dApardt[1])*dfac_v*q_)/m_; 
  out[6] += -(0.1414213562373095*(8.660254037844387*dApardt[2]*f[13]+8.660254037844386*dApardt[1]*f[5]+8.660254037844386*dApardt[0]*f[3])*dfac_v*q_)/m_; 
  out[8] += -(0.7071067811865475*(3.872983346207417*dApardt[2]*f[11]+3.872983346207417*dApardt[1]*f[4]+3.872983346207417*dApardt[0]*f[2])*dfac_v*q_)/m_; 
  out[10] += -(0.1414213562373095*(7.745966692414834*dApardt[1]*f[13]+7.745966692414834*dApardt[2]*f[5]+8.660254037844386*dApardt[0]*f[5]+8.660254037844386*dApardt[1]*f[3])*dfac_v*q_)/m_; 
  out[11] += -(0.07824607964359516*(10.0*dApardt[2]*f[7]+15.65247584249853*dApardt[0]*f[7]+15.65247584249853*f[0]*dApardt[2]+14.0*dApardt[1]*f[1])*dfac_v*q_)/m_; 
  out[12] += -(0.1414213562373095*(17.32050807568877*dApardt[1]*f[11]+17.32050807568877*dApardt[2]*f[4]+19.36491673103708*dApardt[0]*f[4]+19.36491673103708*dApardt[1]*f[2])*dfac_v*q_)/m_; 
  out[14] += -(2.738612787525831*(dApardt[2]*f[17]+dApardt[1]*f[10]+dApardt[0]*f[6])*dfac_v*q_)/m_; 
  out[16] += -(0.1414213562373095*(8.660254037844386*dApardt[1]*f[15]+8.660254037844387*dApardt[0]*f[9])*dfac_v*q_)/m_; 
  out[17] += -(0.02020305089104421*(38.72983346207417*dApardt[2]*f[13]+60.62177826491071*dApardt[0]*f[13]+54.22176684690384*dApardt[1]*f[5]+60.6217782649107*dApardt[2]*f[3])*dfac_v*q_)/m_; 
  out[18] += -(1.224744871391589*(2.0*dApardt[1]*f[17]+2.0*dApardt[2]*f[10]+2.23606797749979*dApardt[0]*f[10]+2.23606797749979*dApardt[1]*f[6])*dfac_v*q_)/m_; 
  out[19] += -(0.1414213562373095*(7.745966692414834*dApardt[2]*f[15]+8.660254037844387*dApardt[0]*f[15]+8.660254037844386*dApardt[1]*f[9])*dfac_v*q_)/m_; 
  return fabs((0.3535533905932737*(2.23606797749979*dApardt[2]*q_-2.0*dApardt[0]*q_))/m_)*dvInv; 
} 
