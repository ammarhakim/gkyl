#include <GyrokineticModDecl.h> 
double GyrokineticVol2x2vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dyInv = 1.0/dxv[1]; 
  double dvInv = 1.0/dxv[2]; 
  double dmInv = 1.0/dxv[3]; 
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
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[16]; 
  alphax[0] = (2.0*BdriftX[0]*m_*wv2)/q_-1.732050807568877*BmagInv[0]*Phi[2]*dfac_y; 
  alphax[1] = -1.732050807568877*BmagInv[0]*Phi[3]*dfac_y; 
  alphax[3] = (1.154700538379252*BdriftX[0]*m_*wv)/(dfac_v*q_); 
  out[1] += 0.4330127018922193*(alphax[3]*f[3]+alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[5] += 0.4330127018922193*(alphax[3]*f[7]+alphax[1]*f[5]+alphax[0]*f[2])*dfac_x; 
  out[6] += 0.4330127018922193*(alphax[1]*f[6]+alphax[0]*f[3]+f[0]*alphax[3])*dfac_x; 
  out[8] += 0.4330127018922193*(alphax[3]*f[10]+alphax[1]*f[8]+alphax[0]*f[4])*dfac_x; 
  out[11] += 0.4330127018922193*(alphax[1]*f[11]+alphax[0]*f[7]+f[2]*alphax[3])*dfac_x; 
  out[12] += 0.4330127018922193*(alphax[3]*f[14]+alphax[1]*f[12]+alphax[0]*f[9])*dfac_x; 
  out[13] += 0.4330127018922193*(alphax[1]*f[13]+alphax[0]*f[10]+alphax[3]*f[4])*dfac_x; 
  out[15] += 0.4330127018922193*(alphax[1]*f[15]+alphax[0]*f[14]+alphax[3]*f[9])*dfac_x; 
  cflFreq += fabs(0.25*alphax[0])*dxInv; 
  double alphay[16]; 
  alphay[0] = (2.0*BdriftY[0]*m_*wv2)/q_+1.732050807568877*BmagInv[0]*Phi[1]*dfac_x; 
  alphay[2] = 1.732050807568877*BmagInv[0]*Phi[3]*dfac_x; 
  alphay[3] = (1.154700538379252*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  out[2] += 0.4330127018922193*(alphay[3]*f[3]+alphay[2]*f[2]+alphay[0]*f[0])*dfac_y; 
  out[5] += 0.4330127018922193*(alphay[3]*f[6]+alphay[2]*f[5]+alphay[0]*f[1])*dfac_y; 
  out[7] += 0.4330127018922193*(alphay[2]*f[7]+alphay[0]*f[3]+f[0]*alphay[3])*dfac_y; 
  out[9] += 0.4330127018922193*(alphay[3]*f[10]+alphay[2]*f[9]+alphay[0]*f[4])*dfac_y; 
  out[11] += 0.4330127018922193*(alphay[2]*f[11]+alphay[0]*f[6]+f[1]*alphay[3])*dfac_y; 
  out[12] += 0.4330127018922193*(alphay[3]*f[13]+alphay[2]*f[12]+alphay[0]*f[8])*dfac_y; 
  out[14] += 0.4330127018922193*(alphay[2]*f[14]+alphay[0]*f[10]+alphay[3]*f[4])*dfac_y; 
  out[15] += 0.4330127018922193*(alphay[2]*f[15]+alphay[0]*f[13]+alphay[3]*f[8])*dfac_y; 
  cflFreq += fabs(0.25*alphay[0])*dyInv; 
  double alphav[16]; 
  alphav[0] = (-1.732050807568877*BdriftY[0]*Phi[2]*dfac_y*wv)-1.732050807568877*BdriftX[0]*Phi[1]*dfac_x*wv; 
  alphav[1] = -1.732050807568877*BdriftY[0]*Phi[3]*dfac_y*wv; 
  alphav[2] = -1.732050807568877*BdriftX[0]*Phi[3]*dfac_x*wv; 
  alphav[3] = (-(1.0*BdriftY[0]*Phi[2]*dfac_y)/dfac_v)-(1.0*BdriftX[0]*Phi[1]*dfac_x)/dfac_v; 
  alphav[6] = -(1.0*BdriftY[0]*Phi[3]*dfac_y)/dfac_v; 
  alphav[7] = -(1.0*BdriftX[0]*Phi[3]*dfac_x)/dfac_v; 
  out[3] += 0.4330127018922193*(alphav[7]*f[7]+alphav[6]*f[6]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[6] += 0.4330127018922193*(alphav[7]*f[11]+alphav[3]*f[6]+f[3]*alphav[6]+alphav[2]*f[5]+alphav[0]*f[1]+f[0]*alphav[1])*dfac_v; 
  out[7] += 0.4330127018922193*(alphav[6]*f[11]+alphav[3]*f[7]+f[3]*alphav[7]+alphav[1]*f[5]+alphav[0]*f[2]+f[0]*alphav[2])*dfac_v; 
  out[10] += 0.4330127018922193*(alphav[7]*f[14]+alphav[6]*f[13]+alphav[3]*f[10]+alphav[2]*f[9]+alphav[1]*f[8]+alphav[0]*f[4])*dfac_v; 
  out[11] += 0.4330127018922193*(alphav[3]*f[11]+alphav[6]*f[7]+f[6]*alphav[7]+alphav[0]*f[5]+alphav[1]*f[2]+f[1]*alphav[2])*dfac_v; 
  out[13] += 0.4330127018922193*(alphav[7]*f[15]+alphav[3]*f[13]+alphav[2]*f[12]+alphav[6]*f[10]+alphav[0]*f[8]+alphav[1]*f[4])*dfac_v; 
  out[14] += 0.4330127018922193*(alphav[6]*f[15]+alphav[3]*f[14]+alphav[1]*f[12]+alphav[7]*f[10]+alphav[0]*f[9]+alphav[2]*f[4])*dfac_v; 
  out[15] += 0.4330127018922193*(alphav[3]*f[15]+alphav[6]*f[14]+alphav[7]*f[13]+alphav[0]*f[12]+alphav[1]*f[9]+alphav[2]*f[8])*dfac_v; 
  cflFreq += fabs(0.25*alphav[0])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol2x2vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dyInv = 1.0/dxv[1]; 
  double dvInv = 1.0/dxv[2]; 
  double dmInv = 1.0/dxv[3]; 
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
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[16]; 
  alphax[0] = (2.0*BdriftX[0]*m_*wv2)/q_+3.464101615137754*Apar[2]*dfac_y*wv+Apar[0]*BdriftX[0]*wv-1.732050807568877*BmagInv[0]*Phi[2]*dfac_y; 
  alphax[1] = 3.464101615137754*Apar[3]*dfac_y*wv+BdriftX[0]*Apar[1]*wv-1.732050807568877*BmagInv[0]*Phi[3]*dfac_y; 
  alphax[2] = BdriftX[0]*Apar[2]*wv; 
  alphax[3] = (1.154700538379252*BdriftX[0]*m_*wv)/(dfac_v*q_); 
  alphax[5] = BdriftX[0]*Apar[3]*wv; 
  out[1] += 0.4330127018922193*(alphax[5]*f[5]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[5] += 0.4330127018922193*(alphax[3]*f[7]+alphax[1]*f[5]+f[1]*alphax[5]+alphax[0]*f[2]+f[0]*alphax[2])*dfac_x; 
  out[6] += 0.4330127018922193*(alphax[5]*f[11]+alphax[2]*f[7]+alphax[1]*f[6]+alphax[0]*f[3]+f[0]*alphax[3])*dfac_x; 
  out[8] += 0.4330127018922193*(alphax[5]*f[12]+alphax[3]*f[10]+alphax[2]*f[9]+alphax[1]*f[8]+alphax[0]*f[4])*dfac_x; 
  out[11] += 0.4330127018922193*(alphax[1]*f[11]+alphax[0]*f[7]+alphax[5]*f[6]+alphax[2]*f[3]+f[2]*alphax[3])*dfac_x; 
  out[12] += 0.4330127018922193*(alphax[3]*f[14]+alphax[1]*f[12]+alphax[0]*f[9]+alphax[5]*f[8]+alphax[2]*f[4])*dfac_x; 
  out[13] += 0.4330127018922193*(alphax[5]*f[15]+alphax[2]*f[14]+alphax[1]*f[13]+alphax[0]*f[10]+alphax[3]*f[4])*dfac_x; 
  out[15] += 0.4330127018922193*(alphax[1]*f[15]+alphax[0]*f[14]+alphax[5]*f[13]+alphax[2]*f[10]+alphax[3]*f[9])*dfac_x; 
  cflFreq += fabs(0.25*alphax[0])*dxInv; 
  double alphay[16]; 
  alphay[0] = (2.0*BdriftY[0]*m_*wv2)/q_-3.464101615137754*Apar[1]*dfac_x*wv+Apar[0]*BdriftY[0]*wv+1.732050807568877*BmagInv[0]*Phi[1]*dfac_x; 
  alphay[1] = BdriftY[0]*Apar[1]*wv; 
  alphay[2] = (-3.464101615137755*Apar[3]*dfac_x*wv)+1.0*BdriftY[0]*Apar[2]*wv+1.732050807568878*BmagInv[0]*Phi[3]*dfac_x; 
  alphay[3] = (1.154700538379252*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alphay[5] = BdriftY[0]*Apar[3]*wv; 
  out[2] += 0.4330127018922193*(alphay[5]*f[5]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0])*dfac_y; 
  out[5] += 0.4330127018922193*(alphay[3]*f[6]+alphay[2]*f[5]+f[2]*alphay[5]+alphay[0]*f[1]+f[0]*alphay[1])*dfac_y; 
  out[7] += 0.4330127018922193*(alphay[5]*f[11]+alphay[2]*f[7]+alphay[1]*f[6]+alphay[0]*f[3]+f[0]*alphay[3])*dfac_y; 
  out[9] += 0.4330127018922193*(alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+alphay[0]*f[4])*dfac_y; 
  out[11] += 0.4330127018922193*(alphay[2]*f[11]+alphay[5]*f[7]+alphay[0]*f[6]+alphay[1]*f[3]+f[1]*alphay[3])*dfac_y; 
  out[12] += 0.4330127018922193*(alphay[3]*f[13]+alphay[2]*f[12]+alphay[5]*f[9]+alphay[0]*f[8]+alphay[1]*f[4])*dfac_y; 
  out[14] += 0.4330127018922193*(alphay[5]*f[15]+alphay[2]*f[14]+alphay[1]*f[13]+alphay[0]*f[10]+alphay[3]*f[4])*dfac_y; 
  out[15] += 0.4330127018922193*(alphay[2]*f[15]+alphay[5]*f[14]+alphay[0]*f[13]+alphay[1]*f[10]+alphay[3]*f[8])*dfac_y; 
  cflFreq += fabs(0.25*alphay[0])*dyInv; 
  double alphav[16]; 
  alphav[0] = (-1.732050807568877*BdriftY[0]*Phi[2]*dfac_y*wv)-1.732050807568877*BdriftX[0]*Phi[1]*dfac_x*wv+(3.0*Apar[1]*Phi[2]*dfac_x*dfac_y*q_)/m_-(3.0*Phi[1]*Apar[2]*dfac_x*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[1]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[0]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[2]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*Apar[0]*BdriftX[0]*Phi[1]*dfac_x*q_)/m_; 
  alphav[1] = (-1.732050807568877*BdriftY[0]*Phi[3]*dfac_y*wv)+(3.0*Apar[1]*Phi[3]*dfac_x*dfac_y*q_)/m_-(3.0*Phi[1]*Apar[3]*dfac_x*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[0]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[1]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[3]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[1]*Phi[1]*dfac_x*q_)/m_; 
  alphav[2] = (-1.732050807568877*BdriftX[0]*Phi[3]*dfac_x*wv)-(3.0*Apar[2]*Phi[3]*dfac_x*dfac_y*q_)/m_+(3.0*Phi[2]*Apar[3]*dfac_x*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[3]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[2]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftX[0]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Phi[1]*Apar[2]*dfac_x*q_)/m_; 
  alphav[3] = (-(1.0*BdriftY[0]*Phi[2]*dfac_y)/dfac_v)-(1.0*BdriftX[0]*Phi[1]*dfac_x)/dfac_v; 
  alphav[5] = (-(0.8660254037844386*BdriftY[0]*Apar[2]*Phi[3]*dfac_y*q_)/m_)-(0.8660254037844386*BdriftY[0]*Phi[2]*Apar[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[1]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Phi[1]*Apar[3]*dfac_x*q_)/m_; 
  alphav[6] = -(1.0*BdriftY[0]*Phi[3]*dfac_y)/dfac_v; 
  alphav[7] = -(1.0*BdriftX[0]*Phi[3]*dfac_x)/dfac_v; 
  out[3] += 0.4330127018922193*(alphav[7]*f[7]+alphav[6]*f[6]+alphav[5]*f[5]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[6] += 0.4330127018922193*(alphav[7]*f[11]+alphav[3]*f[6]+f[3]*alphav[6]+alphav[2]*f[5]+f[2]*alphav[5]+alphav[0]*f[1]+f[0]*alphav[1])*dfac_v; 
  out[7] += 0.4330127018922193*(alphav[6]*f[11]+alphav[3]*f[7]+f[3]*alphav[7]+alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[2]+f[0]*alphav[2])*dfac_v; 
  out[10] += 0.4330127018922193*(alphav[7]*f[14]+alphav[6]*f[13]+alphav[5]*f[12]+alphav[3]*f[10]+alphav[2]*f[9]+alphav[1]*f[8]+alphav[0]*f[4])*dfac_v; 
  out[11] += 0.4330127018922193*(alphav[3]*f[11]+alphav[6]*f[7]+f[6]*alphav[7]+alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[2]+f[1]*alphav[2])*dfac_v; 
  out[13] += 0.4330127018922193*(alphav[7]*f[15]+alphav[3]*f[13]+alphav[2]*f[12]+alphav[6]*f[10]+alphav[5]*f[9]+alphav[0]*f[8]+alphav[1]*f[4])*dfac_v; 
  out[14] += 0.4330127018922193*(alphav[6]*f[15]+alphav[3]*f[14]+alphav[1]*f[12]+alphav[7]*f[10]+alphav[0]*f[9]+alphav[5]*f[8]+alphav[2]*f[4])*dfac_v; 
  out[15] += 0.4330127018922193*(alphav[3]*f[15]+alphav[6]*f[14]+alphav[7]*f[13]+alphav[0]*f[12]+alphav[1]*f[9]+alphav[2]*f[8]+f[4]*alphav[5])*dfac_v; 
  cflFreq += fabs(0.25*alphav[0])*dvInv; 
  return cflFreq; 
} 
double GyrokineticVol2x2vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dyInv = 1.0/dxv[1]; 
  double dvInv = 1.0/dxv[2]; 
  double dmInv = 1.0/dxv[3]; 
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
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[16]; 
  alphax[0] = (2.0*BdriftX[0]*m_*wv2)/q_-1.732050807568877*BmagInv[1]*Phi[3]*dfac_y-1.732050807568877*BmagInv[0]*Phi[2]*dfac_y; 
  alphax[1] = (2.0*BdriftX[1]*m_*wv2)/q_-1.732050807568877*BmagInv[0]*Phi[3]*dfac_y-1.732050807568877*BmagInv[1]*Phi[2]*dfac_y; 
  alphax[3] = (1.154700538379252*BdriftX[0]*m_*wv)/(dfac_v*q_); 
  alphax[6] = (1.154700538379252*BdriftX[1]*m_*wv)/(dfac_v*q_); 
  out[1] += 0.4330127018922193*(alphax[6]*f[6]+alphax[3]*f[3]+alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[5] += 0.4330127018922193*(alphax[6]*f[11]+alphax[3]*f[7]+alphax[1]*f[5]+alphax[0]*f[2])*dfac_x; 
  out[6] += 0.4330127018922193*(alphax[1]*f[6]+f[1]*alphax[6]+alphax[0]*f[3]+f[0]*alphax[3])*dfac_x; 
  out[8] += 0.4330127018922193*(alphax[6]*f[13]+alphax[3]*f[10]+alphax[1]*f[8]+alphax[0]*f[4])*dfac_x; 
  out[11] += 0.4330127018922193*(alphax[1]*f[11]+alphax[0]*f[7]+f[5]*alphax[6]+f[2]*alphax[3])*dfac_x; 
  out[12] += 0.4330127018922193*(alphax[6]*f[15]+alphax[3]*f[14]+alphax[1]*f[12]+alphax[0]*f[9])*dfac_x; 
  out[13] += 0.4330127018922193*(alphax[1]*f[13]+alphax[0]*f[10]+alphax[6]*f[8]+alphax[3]*f[4])*dfac_x; 
  out[15] += 0.4330127018922193*(alphax[1]*f[15]+alphax[0]*f[14]+alphax[6]*f[12]+alphax[3]*f[9])*dfac_x; 
  cflFreq += fabs(0.25*alphax[0])*dxInv; 
  double alphay[16]; 
  alphay[0] = (2.0*BdriftY[0]*m_*wv2)/q_+(1.732050807568877*BmagInv[0]*Bmag[1]*dfac_x*wm)/q_+1.732050807568877*BmagInv[0]*Phi[1]*dfac_x; 
  alphay[1] = (2.0*BdriftY[1]*m_*wv2)/q_+(1.732050807568877*Bmag[1]*BmagInv[1]*dfac_x*wm)/q_+1.732050807568877*BmagInv[1]*Phi[1]*dfac_x; 
  alphay[2] = 1.732050807568877*BmagInv[0]*Phi[3]*dfac_x; 
  alphay[3] = (1.154700538379252*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alphay[4] = (BmagInv[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alphay[5] = 1.732050807568877*BmagInv[1]*Phi[3]*dfac_x; 
  alphay[6] = (1.154700538379252*BdriftY[1]*m_*wv)/(dfac_v*q_); 
  alphay[8] = (Bmag[1]*BmagInv[1]*dfac_x)/(dfac_m*q_); 
  out[2] += 0.4330127018922193*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0])*dfac_y; 
  out[5] += 0.4330127018922193*(alphay[4]*f[8]+f[4]*alphay[8]+alphay[3]*f[6]+f[3]*alphay[6]+alphay[2]*f[5]+f[2]*alphay[5]+alphay[0]*f[1]+f[0]*alphay[1])*dfac_y; 
  out[7] += 0.4330127018922193*(alphay[8]*f[13]+alphay[5]*f[11]+alphay[4]*f[10]+alphay[2]*f[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphay[0]*f[3]+f[0]*alphay[3])*dfac_y; 
  out[9] += 0.4330127018922193*(alphay[6]*f[13]+alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4])*dfac_y; 
  out[11] += 0.4330127018922193*(alphay[4]*f[13]+alphay[2]*f[11]+alphay[8]*f[10]+alphay[5]*f[7]+alphay[0]*f[6]+f[0]*alphay[6]+alphay[1]*f[3]+f[1]*alphay[3])*dfac_y; 
  out[12] += 0.4330127018922193*(alphay[3]*f[13]+alphay[2]*f[12]+alphay[6]*f[10]+alphay[5]*f[9]+alphay[0]*f[8]+f[0]*alphay[8]+alphay[1]*f[4]+f[1]*alphay[4])*dfac_y; 
  out[14] += 0.4330127018922193*(alphay[5]*f[15]+alphay[2]*f[14]+alphay[1]*f[13]+alphay[0]*f[10]+alphay[6]*f[8]+f[6]*alphay[8]+alphay[3]*f[4]+f[3]*alphay[4])*dfac_y; 
  out[15] += 0.4330127018922193*(alphay[2]*f[15]+alphay[5]*f[14]+alphay[0]*f[13]+alphay[1]*f[10]+alphay[3]*f[8]+f[3]*alphay[8]+alphay[4]*f[6]+f[4]*alphay[6])*dfac_y; 
  cflFreq += fabs(0.25*alphay[0])*dyInv; 
  double alphav[16]; 
  alphav[0] = (-(1.732050807568877*BdriftX[0]*Bmag[1]*dfac_x*wm*wv)/q_)-1.732050807568877*BdriftY[1]*Phi[3]*dfac_y*wv-1.732050807568877*BdriftY[0]*Phi[2]*dfac_y*wv-1.732050807568877*BdriftX[0]*Phi[1]*dfac_x*wv; 
  alphav[1] = (-(1.732050807568877*BdriftX[1]*Bmag[1]*dfac_x*wm*wv)/q_)-1.732050807568877*BdriftY[0]*Phi[3]*dfac_y*wv-1.732050807568877*BdriftY[1]*Phi[2]*dfac_y*wv-1.732050807568877*BdriftX[1]*Phi[1]*dfac_x*wv; 
  alphav[2] = -1.732050807568877*BdriftX[0]*Phi[3]*dfac_x*wv; 
  alphav[3] = (-(1.0*BdriftX[0]*Bmag[1]*dfac_x*wm)/(dfac_v*q_))-(1.0*BdriftY[1]*Phi[3]*dfac_y)/dfac_v-(1.0*BdriftY[0]*Phi[2]*dfac_y)/dfac_v-(1.0*BdriftX[0]*Phi[1]*dfac_x)/dfac_v; 
  alphav[4] = -(1.0*BdriftX[0]*Bmag[1]*dfac_x*wv)/(dfac_m*q_); 
  alphav[5] = -1.732050807568877*BdriftX[1]*Phi[3]*dfac_x*wv; 
  alphav[6] = (-(1.0*BdriftX[1]*Bmag[1]*dfac_x*wm)/(dfac_v*q_))-(1.0*BdriftY[0]*Phi[3]*dfac_y)/dfac_v-(1.0*BdriftY[1]*Phi[2]*dfac_y)/dfac_v-(1.0*BdriftX[1]*Phi[1]*dfac_x)/dfac_v; 
  alphav[7] = -(1.0*BdriftX[0]*Phi[3]*dfac_x)/dfac_v; 
  alphav[8] = -(1.0*BdriftX[1]*Bmag[1]*dfac_x*wv)/(dfac_m*q_); 
  alphav[10] = -(0.5773502691896258*BdriftX[0]*Bmag[1]*dfac_x)/(dfac_m*dfac_v*q_); 
  alphav[11] = -(1.0*BdriftX[1]*Phi[3]*dfac_x)/dfac_v; 
  alphav[13] = -(0.5773502691896258*BdriftX[1]*Bmag[1]*dfac_x)/(dfac_m*dfac_v*q_); 
  out[3] += 0.4330127018922193*(alphav[13]*f[13]+alphav[11]*f[11]+alphav[10]*f[10]+alphav[8]*f[8]+alphav[7]*f[7]+alphav[6]*f[6]+alphav[5]*f[5]+alphav[4]*f[4]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[6] += 0.4330127018922193*(alphav[10]*f[13]+f[10]*alphav[13]+alphav[7]*f[11]+f[7]*alphav[11]+alphav[4]*f[8]+f[4]*alphav[8]+alphav[3]*f[6]+f[3]*alphav[6]+alphav[2]*f[5]+f[2]*alphav[5]+alphav[0]*f[1]+f[0]*alphav[1])*dfac_v; 
  out[7] += 0.4330127018922193*(alphav[13]*f[15]+alphav[10]*f[14]+alphav[8]*f[12]+alphav[6]*f[11]+f[6]*alphav[11]+alphav[4]*f[9]+alphav[3]*f[7]+f[3]*alphav[7]+alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[2]+f[0]*alphav[2])*dfac_v; 
  out[10] += 0.4330127018922193*(alphav[11]*f[15]+alphav[7]*f[14]+alphav[6]*f[13]+f[6]*alphav[13]+alphav[5]*f[12]+alphav[3]*f[10]+f[3]*alphav[10]+alphav[2]*f[9]+alphav[1]*f[8]+f[1]*alphav[8]+alphav[0]*f[4]+f[0]*alphav[4])*dfac_v; 
  out[11] += 0.4330127018922193*(alphav[10]*f[15]+alphav[13]*f[14]+alphav[4]*f[12]+alphav[3]*f[11]+f[3]*alphav[11]+alphav[8]*f[9]+alphav[6]*f[7]+f[6]*alphav[7]+alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[2]+f[1]*alphav[2])*dfac_v; 
  out[13] += 0.4330127018922193*(alphav[7]*f[15]+alphav[11]*f[14]+alphav[3]*f[13]+f[3]*alphav[13]+alphav[2]*f[12]+alphav[6]*f[10]+f[6]*alphav[10]+alphav[5]*f[9]+alphav[0]*f[8]+f[0]*alphav[8]+alphav[1]*f[4]+f[1]*alphav[4])*dfac_v; 
  out[14] += 0.4330127018922193*(alphav[6]*f[15]+alphav[3]*f[14]+alphav[11]*f[13]+f[11]*alphav[13]+alphav[1]*f[12]+alphav[7]*f[10]+f[7]*alphav[10]+alphav[0]*f[9]+alphav[5]*f[8]+f[5]*alphav[8]+alphav[2]*f[4]+f[2]*alphav[4])*dfac_v; 
  out[15] += 0.4330127018922193*(alphav[3]*f[15]+alphav[6]*f[14]+alphav[7]*f[13]+f[7]*alphav[13]+alphav[0]*f[12]+alphav[10]*f[11]+f[10]*alphav[11]+alphav[1]*f[9]+alphav[2]*f[8]+f[2]*alphav[8]+alphav[4]*f[5]+f[4]*alphav[5])*dfac_v; 
  cflFreq += fabs(0.25*alphav[0])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticVol2x2vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dxInv = 1.0/dxv[0]; 
  double dyInv = 1.0/dxv[1]; 
  double dvInv = 1.0/dxv[2]; 
  double dmInv = 1.0/dxv[3]; 
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
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphax[16]; 
  alphax[0] = (2.0*BdriftX[0]*m_*wv2)/q_+3.464101615137754*Apar[2]*dfac_y*wv+Apar[1]*BdriftX[1]*wv+Apar[0]*BdriftX[0]*wv-1.732050807568877*BmagInv[1]*Phi[3]*dfac_y-1.732050807568877*BmagInv[0]*Phi[2]*dfac_y; 
  alphax[1] = (2.0*BdriftX[1]*m_*wv2)/q_+3.464101615137754*Apar[3]*dfac_y*wv+Apar[0]*BdriftX[1]*wv+BdriftX[0]*Apar[1]*wv-1.732050807568877*BmagInv[0]*Phi[3]*dfac_y-1.732050807568877*BmagInv[1]*Phi[2]*dfac_y; 
  alphax[2] = BdriftX[1]*Apar[3]*wv+BdriftX[0]*Apar[2]*wv; 
  alphax[3] = (1.154700538379252*BdriftX[0]*m_*wv)/(dfac_v*q_); 
  alphax[5] = BdriftX[0]*Apar[3]*wv+BdriftX[1]*Apar[2]*wv; 
  alphax[6] = (1.154700538379252*BdriftX[1]*m_*wv)/(dfac_v*q_); 
  out[1] += 0.4330127018922193*(alphax[6]*f[6]+alphax[5]*f[5]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0])*dfac_x; 
  out[5] += 0.4330127018922193*(alphax[6]*f[11]+alphax[3]*f[7]+alphax[1]*f[5]+f[1]*alphax[5]+alphax[0]*f[2]+f[0]*alphax[2])*dfac_x; 
  out[6] += 0.4330127018922193*(alphax[5]*f[11]+alphax[2]*f[7]+alphax[1]*f[6]+f[1]*alphax[6]+alphax[0]*f[3]+f[0]*alphax[3])*dfac_x; 
  out[8] += 0.4330127018922193*(alphax[6]*f[13]+alphax[5]*f[12]+alphax[3]*f[10]+alphax[2]*f[9]+alphax[1]*f[8]+alphax[0]*f[4])*dfac_x; 
  out[11] += 0.4330127018922193*(alphax[1]*f[11]+alphax[0]*f[7]+alphax[5]*f[6]+f[5]*alphax[6]+alphax[2]*f[3]+f[2]*alphax[3])*dfac_x; 
  out[12] += 0.4330127018922193*(alphax[6]*f[15]+alphax[3]*f[14]+alphax[1]*f[12]+alphax[0]*f[9]+alphax[5]*f[8]+alphax[2]*f[4])*dfac_x; 
  out[13] += 0.4330127018922193*(alphax[5]*f[15]+alphax[2]*f[14]+alphax[1]*f[13]+alphax[0]*f[10]+alphax[6]*f[8]+alphax[3]*f[4])*dfac_x; 
  out[15] += 0.4330127018922193*(alphax[1]*f[15]+alphax[0]*f[14]+alphax[5]*f[13]+alphax[6]*f[12]+alphax[2]*f[10]+alphax[3]*f[9])*dfac_x; 
  cflFreq += fabs(0.25*alphax[0])*dxInv; 
  double alphay[16]; 
  alphay[0] = (2.0*BdriftY[0]*m_*wv2)/q_-3.464101615137754*Apar[1]*dfac_x*wv+Apar[1]*BdriftY[1]*wv+Apar[0]*BdriftY[0]*wv+(1.732050807568877*BmagInv[0]*Bmag[1]*dfac_x*wm)/q_+1.732050807568877*BmagInv[0]*Phi[1]*dfac_x; 
  alphay[1] = (2.0*BdriftY[1]*m_*wv2)/q_+Apar[0]*BdriftY[1]*wv+BdriftY[0]*Apar[1]*wv+(1.732050807568877*Bmag[1]*BmagInv[1]*dfac_x*wm)/q_+1.732050807568877*BmagInv[1]*Phi[1]*dfac_x; 
  alphay[2] = (-3.464101615137755*Apar[3]*dfac_x*wv)+1.0*BdriftY[1]*Apar[3]*wv+1.0*BdriftY[0]*Apar[2]*wv+1.732050807568878*BmagInv[0]*Phi[3]*dfac_x; 
  alphay[3] = (1.154700538379252*BdriftY[0]*m_*wv)/(dfac_v*q_); 
  alphay[4] = (BmagInv[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alphay[5] = BdriftY[0]*Apar[3]*wv+BdriftY[1]*Apar[2]*wv+1.732050807568877*BmagInv[1]*Phi[3]*dfac_x; 
  alphay[6] = (1.154700538379252*BdriftY[1]*m_*wv)/(dfac_v*q_); 
  alphay[8] = (Bmag[1]*BmagInv[1]*dfac_x)/(dfac_m*q_); 
  out[2] += 0.4330127018922193*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0])*dfac_y; 
  out[5] += 0.4330127018922193*(alphay[4]*f[8]+f[4]*alphay[8]+alphay[3]*f[6]+f[3]*alphay[6]+alphay[2]*f[5]+f[2]*alphay[5]+alphay[0]*f[1]+f[0]*alphay[1])*dfac_y; 
  out[7] += 0.4330127018922193*(alphay[8]*f[13]+alphay[5]*f[11]+alphay[4]*f[10]+alphay[2]*f[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphay[0]*f[3]+f[0]*alphay[3])*dfac_y; 
  out[9] += 0.4330127018922193*(alphay[6]*f[13]+alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4])*dfac_y; 
  out[11] += 0.4330127018922193*(alphay[4]*f[13]+alphay[2]*f[11]+alphay[8]*f[10]+alphay[5]*f[7]+alphay[0]*f[6]+f[0]*alphay[6]+alphay[1]*f[3]+f[1]*alphay[3])*dfac_y; 
  out[12] += 0.4330127018922193*(alphay[3]*f[13]+alphay[2]*f[12]+alphay[6]*f[10]+alphay[5]*f[9]+alphay[0]*f[8]+f[0]*alphay[8]+alphay[1]*f[4]+f[1]*alphay[4])*dfac_y; 
  out[14] += 0.4330127018922193*(alphay[5]*f[15]+alphay[2]*f[14]+alphay[1]*f[13]+alphay[0]*f[10]+alphay[6]*f[8]+f[6]*alphay[8]+alphay[3]*f[4]+f[3]*alphay[4])*dfac_y; 
  out[15] += 0.4330127018922193*(alphay[2]*f[15]+alphay[5]*f[14]+alphay[0]*f[13]+alphay[1]*f[10]+alphay[3]*f[8]+f[3]*alphay[8]+alphay[4]*f[6]+f[4]*alphay[6])*dfac_y; 
  cflFreq += fabs(0.25*alphay[0])*dyInv; 
  double alphav[16]; 
  alphav[0] = (-(1.732050807568877*BdriftX[0]*Bmag[1]*dfac_x*wm*wv)/q_)-1.732050807568877*BdriftY[1]*Phi[3]*dfac_y*wv-1.732050807568877*BdriftY[0]*Phi[2]*dfac_y*wv-1.732050807568877*BdriftX[0]*Phi[1]*dfac_x*wv-(3.0*Bmag[1]*Apar[2]*dfac_x*dfac_y*wm)/m_-(0.8660254037844386*Apar[1]*BdriftX[1]*Bmag[1]*dfac_x*wm)/m_-(0.8660254037844386*Apar[0]*BdriftX[0]*Bmag[1]*dfac_x*wm)/m_+(3.0*Apar[1]*Phi[2]*dfac_x*dfac_y*q_)/m_-(3.0*Phi[1]*Apar[2]*dfac_x*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[1]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[1]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*Apar[1]*BdriftY[1]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[0]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*BdriftX[1]*Apar[3]*Phi[3]*dfac_x*q2)/(m_*q_)-(0.8660254037844386*BdriftX[0]*Apar[2]*Phi[3]*dfac_x*q2)/(m_*q_)-(0.8660254037844386*Apar[1]*BdriftX[1]*Phi[1]*dfac_x*q2)/(m_*q_)-(0.8660254037844386*Apar[0]*BdriftX[0]*Phi[1]*dfac_x*q2)/(m_*q_); 
  alphav[1] = (-(1.732050807568877*BdriftX[1]*Bmag[1]*dfac_x*wm*wv)/q_)-1.732050807568877*BdriftY[0]*Phi[3]*dfac_y*wv-1.732050807568877*BdriftY[1]*Phi[2]*dfac_y*wv-1.732050807568877*BdriftX[1]*Phi[1]*dfac_x*wv-(3.0*Bmag[1]*Apar[3]*dfac_x*dfac_y*wm)/m_-(0.8660254037844386*Apar[0]*BdriftX[1]*Bmag[1]*dfac_x*wm)/m_-(0.8660254037844386*BdriftX[0]*Apar[1]*Bmag[1]*dfac_x*wm)/m_+(3.0*Apar[1]*Phi[3]*dfac_x*dfac_y*q_)/m_-(3.0*Phi[1]*Apar[3]*dfac_x*dfac_y*q_)/m_-(1.55884572681199*Apar[1]*BdriftY[1]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[0]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftY[1]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[1]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[3]*Phi[3]*dfac_x*q2)/(m_*q_)-(0.8660254037844386*BdriftX[1]*Apar[2]*Phi[3]*dfac_x*q2)/(m_*q_)-(0.8660254037844386*Apar[0]*BdriftX[1]*Phi[1]*dfac_x*q2)/(m_*q_)-(0.8660254037844386*BdriftX[0]*Apar[1]*Phi[1]*dfac_x*q2)/(m_*q_); 
  alphav[2] = (-1.732050807568877*BdriftX[0]*Phi[3]*dfac_x*wv)-(0.8660254037844386*BdriftX[1]*Bmag[1]*Apar[3]*dfac_x*wm)/m_-(0.8660254037844386*BdriftX[0]*Bmag[1]*Apar[2]*dfac_x*wm)/m_-(3.0*Apar[2]*Phi[3]*dfac_x*dfac_y*q_)/m_+(3.0*Phi[2]*Apar[3]*dfac_x*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[3]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[1]*Apar[2]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[1]*Phi[2]*Apar[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[2]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*Apar[1]*BdriftX[1]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*Apar[0]*BdriftX[0]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[1]*Phi[1]*Apar[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Phi[1]*Apar[2]*dfac_x*q_)/m_; 
  alphav[3] = (-(1.0*BdriftX[0]*Bmag[1]*dfac_x*wm)/(dfac_v*q_))-(1.0*BdriftY[1]*Phi[3]*dfac_y)/dfac_v-(1.0*BdriftY[0]*Phi[2]*dfac_y)/dfac_v-(1.0*BdriftX[0]*Phi[1]*dfac_x)/dfac_v; 
  alphav[4] = (-(1.0*BdriftX[0]*Bmag[1]*dfac_x*wv)/(dfac_m*q_))-(1.732050807568877*Bmag[1]*Apar[2]*dfac_x*dfac_y)/(dfac_m*m_)-(0.5*Apar[1]*BdriftX[1]*Bmag[1]*dfac_x)/(dfac_m*m_)-(0.5*Apar[0]*BdriftX[0]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alphav[5] = (-1.732050807568877*BdriftX[1]*Phi[3]*dfac_x*wv)-(0.8660254037844386*BdriftX[0]*Bmag[1]*Apar[3]*dfac_x*wm)/m_-(0.8660254037844386*BdriftX[1]*Bmag[1]*Apar[2]*dfac_x*wm)/m_-(1.55884572681199*BdriftY[1]*Apar[3]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Apar[2]*Phi[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[0]*Phi[2]*Apar[3]*dfac_y*q_)/m_-(0.8660254037844386*BdriftY[1]*Apar[2]*Phi[2]*dfac_y*q_)/m_-(0.8660254037844386*Apar[0]*BdriftX[1]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Apar[1]*Phi[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[0]*Phi[1]*Apar[3]*dfac_x*q_)/m_-(0.8660254037844386*BdriftX[1]*Phi[1]*Apar[2]*dfac_x*q_)/m_; 
  alphav[6] = (-(1.0*BdriftX[1]*Bmag[1]*dfac_x*wm)/(dfac_v*q_))-(1.0*BdriftY[0]*Phi[3]*dfac_y)/dfac_v-(1.0*BdriftY[1]*Phi[2]*dfac_y)/dfac_v-(1.0*BdriftX[1]*Phi[1]*dfac_x)/dfac_v; 
  alphav[7] = -(1.0*BdriftX[0]*Phi[3]*dfac_x)/dfac_v; 
  alphav[8] = (-(1.0*BdriftX[1]*Bmag[1]*dfac_x*wv)/(dfac_m*q_))-(1.732050807568877*Bmag[1]*Apar[3]*dfac_x*dfac_y)/(dfac_m*m_)-(0.5*Apar[0]*BdriftX[1]*Bmag[1]*dfac_x)/(dfac_m*m_)-(0.5*BdriftX[0]*Apar[1]*Bmag[1]*dfac_x)/(dfac_m*m_); 
  alphav[9] = (-(0.5*BdriftX[1]*Bmag[1]*Apar[3]*dfac_x)/(dfac_m*m_))-(0.5*BdriftX[0]*Bmag[1]*Apar[2]*dfac_x)/(dfac_m*m_); 
  alphav[10] = -(0.5773502691896258*BdriftX[0]*Bmag[1]*dfac_x)/(dfac_m*dfac_v*q_); 
  alphav[11] = -(1.0*BdriftX[1]*Phi[3]*dfac_x)/dfac_v; 
  alphav[12] = (-(0.5*BdriftX[0]*Bmag[1]*Apar[3]*dfac_x)/(dfac_m*m_))-(0.5*BdriftX[1]*Bmag[1]*Apar[2]*dfac_x)/(dfac_m*m_); 
  alphav[13] = -(0.5773502691896258*BdriftX[1]*Bmag[1]*dfac_x)/(dfac_m*dfac_v*q_); 
  out[3] += 0.4330127018922193*(alphav[13]*f[13]+alphav[12]*f[12]+alphav[11]*f[11]+alphav[10]*f[10]+alphav[9]*f[9]+alphav[8]*f[8]+alphav[7]*f[7]+alphav[6]*f[6]+alphav[5]*f[5]+alphav[4]*f[4]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0])*dfac_v; 
  out[6] += 0.4330127018922193*(alphav[10]*f[13]+f[10]*alphav[13]+alphav[9]*f[12]+f[9]*alphav[12]+alphav[7]*f[11]+f[7]*alphav[11]+alphav[4]*f[8]+f[4]*alphav[8]+alphav[3]*f[6]+f[3]*alphav[6]+alphav[2]*f[5]+f[2]*alphav[5]+alphav[0]*f[1]+f[0]*alphav[1])*dfac_v; 
  out[7] += 0.4330127018922193*(alphav[13]*f[15]+alphav[10]*f[14]+alphav[8]*f[12]+f[8]*alphav[12]+alphav[6]*f[11]+f[6]*alphav[11]+alphav[4]*f[9]+f[4]*alphav[9]+alphav[3]*f[7]+f[3]*alphav[7]+alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[2]+f[0]*alphav[2])*dfac_v; 
  out[10] += 0.4330127018922193*(alphav[11]*f[15]+alphav[7]*f[14]+alphav[6]*f[13]+f[6]*alphav[13]+alphav[5]*f[12]+f[5]*alphav[12]+alphav[3]*f[10]+f[3]*alphav[10]+alphav[2]*f[9]+f[2]*alphav[9]+alphav[1]*f[8]+f[1]*alphav[8]+alphav[0]*f[4]+f[0]*alphav[4])*dfac_v; 
  out[11] += 0.4330127018922193*(alphav[10]*f[15]+alphav[13]*f[14]+alphav[4]*f[12]+f[4]*alphav[12]+alphav[3]*f[11]+f[3]*alphav[11]+alphav[8]*f[9]+f[8]*alphav[9]+alphav[6]*f[7]+f[6]*alphav[7]+alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[2]+f[1]*alphav[2])*dfac_v; 
  out[13] += 0.4330127018922193*(alphav[7]*f[15]+alphav[11]*f[14]+alphav[3]*f[13]+f[3]*alphav[13]+alphav[2]*f[12]+f[2]*alphav[12]+alphav[6]*f[10]+f[6]*alphav[10]+alphav[5]*f[9]+f[5]*alphav[9]+alphav[0]*f[8]+f[0]*alphav[8]+alphav[1]*f[4]+f[1]*alphav[4])*dfac_v; 
  out[14] += 0.4330127018922193*(alphav[6]*f[15]+alphav[3]*f[14]+alphav[11]*f[13]+f[11]*alphav[13]+alphav[1]*f[12]+f[1]*alphav[12]+alphav[7]*f[10]+f[7]*alphav[10]+alphav[0]*f[9]+f[0]*alphav[9]+alphav[5]*f[8]+f[5]*alphav[8]+alphav[2]*f[4]+f[2]*alphav[4])*dfac_v; 
  out[15] += 0.4330127018922193*(alphav[3]*f[15]+alphav[6]*f[14]+alphav[7]*f[13]+f[7]*alphav[13]+alphav[0]*f[12]+f[0]*alphav[12]+alphav[10]*f[11]+f[10]*alphav[11]+alphav[1]*f[9]+f[1]*alphav[9]+alphav[2]*f[8]+f[2]*alphav[8]+alphav[4]*f[5]+f[4]*alphav[5])*dfac_v; 
  cflFreq += fabs(0.25*alphav[0])*dvInv; 
  return cflFreq; 
} 
double EmGyrokineticStep2Vol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[2]; 
  double dfac_v = 2.0/dxv[2]; 
  out[3] += -(0.8660254037844386*(dApardt[3]*f[5]+dApardt[2]*f[2]+dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[6] += -(0.8660254037844386*(dApardt[2]*f[5]+f[2]*dApardt[3]+dApardt[0]*f[1]+f[0]*dApardt[1])*dfac_v*q_)/m_; 
  out[7] += -(0.8660254037844386*(dApardt[1]*f[5]+f[1]*dApardt[3]+dApardt[0]*f[2]+f[0]*dApardt[2])*dfac_v*q_)/m_; 
  out[10] += -(0.8660254037844386*(dApardt[3]*f[12]+dApardt[2]*f[9]+dApardt[1]*f[8]+dApardt[0]*f[4])*dfac_v*q_)/m_; 
  out[11] += -(0.8660254037844386*(dApardt[0]*f[5]+f[0]*dApardt[3]+dApardt[1]*f[2]+f[1]*dApardt[2])*dfac_v*q_)/m_; 
  out[13] += -(0.8660254037844386*(dApardt[2]*f[12]+dApardt[3]*f[9]+dApardt[0]*f[8]+dApardt[1]*f[4])*dfac_v*q_)/m_; 
  out[14] += -(0.8660254037844386*(dApardt[1]*f[12]+dApardt[0]*f[9]+dApardt[3]*f[8]+dApardt[2]*f[4])*dfac_v*q_)/m_; 
  out[15] += -(0.8660254037844386*(dApardt[0]*f[12]+dApardt[1]*f[9]+dApardt[2]*f[8]+dApardt[3]*f[4])*dfac_v*q_)/m_; 
  return fabs(-(0.25*Bmag[0]*dApardt[0]*q_)/m_)*dvInv; 
} 
