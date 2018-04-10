#include <GyrokineticModDecl.h> 
double GyrokineticVol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *BcurvY, const double *Phi, const double *f, double *out) 
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
  out[1] += -0.75*(BmagInv[0]*f[1]*Phi[3]+f[0]*BmagInv[1]*Phi[3]+BmagInv[1]*f[1]*Phi[2]+BmagInv[0]*f[0]*Phi[2])*dfac_x*dfac_y; 
  out[2] += (0.25*(3.464101615137754*(BcurvY[1]*f[1]+BcurvY[0]*f[0])*dfac_m*dfac_v*dfac_y*m_*wv2+2.0*(BcurvY[1]*f[6]+BcurvY[0]*f[3])*dfac_m*dfac_y*m_*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*(Bmag[1]*BmagInv[1]*f[1]*wm+BmagInv[0]*f[0]*Bmag[1]*wm+BmagInv[1]*Phi[3]*f[5]*q_+BmagInv[0]*f[2]*Phi[3]*q_+BmagInv[1]*Phi[1]*f[1]*q_+BmagInv[0]*f[0]*Phi[1]*q_)+1.732050807568877*Bmag[1]*(BmagInv[1]*f[8]+BmagInv[0]*f[4])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[3] += -0.25*(3.0*(BcurvY[0]*f[1]*Phi[3]+f[0]*BcurvY[1]*Phi[3]+BcurvY[1]*f[1]*Phi[2]+BcurvY[0]*f[0]*Phi[2])*dfac_v*dfac_y*wv+1.732050807568877*(BcurvY[0]*Phi[3]*f[6]+BcurvY[1]*Phi[2]*f[6]+BcurvY[1]*Phi[3]*f[3]+BcurvY[0]*Phi[2]*f[3])*dfac_y); 
  out[5] += (0.25*(3.464101615137754*(BcurvY[0]*f[1]+f[0]*BcurvY[1])*dfac_m*dfac_v*dfac_y*m_*wv2+2.0*(BcurvY[0]*f[6]+BcurvY[1]*f[3])*dfac_m*dfac_y*m_*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*(BmagInv[0]*Bmag[1]*f[1]*wm+f[0]*Bmag[1]*BmagInv[1]*wm-1.0*BmagInv[1]*Phi[2]*f[5]*q_-1.0*BmagInv[0]*Phi[2]*f[2]*q_+BmagInv[0]*Phi[1]*f[1]*q_+f[0]*BmagInv[1]*Phi[1]*q_)+1.732050807568877*Bmag[1]*(BmagInv[0]*f[8]+BmagInv[1]*f[4])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[6] += -0.05*(3.0*(9.0*BcurvY[1]*f[1]*Phi[3]+5.0*BcurvY[0]*f[0]*Phi[3]+5.0*BcurvY[0]*f[1]*Phi[2]+5.0*f[0]*BcurvY[1]*Phi[2])*dfac_v*dfac_y*wv+15.0*(BmagInv[0]*Phi[3]*f[6]+BmagInv[1]*Phi[2]*f[6]+BmagInv[1]*Phi[3]*f[3]+BmagInv[0]*Phi[2]*f[3])*dfac_x*dfac_y+1.732050807568877*(9.0*BcurvY[1]*Phi[3]*f[6]+5.0*BcurvY[0]*Phi[2]*f[6]+5.0*BcurvY[0]*Phi[3]*f[3]+5.0*BcurvY[1]*Phi[2]*f[3])*dfac_y); 
  out[7] += (0.25*(1.732050807568877*dfac_m*dfac_v*dfac_y*(2.0*BcurvY[1]*f[6]*m_*wv2+2.0*BcurvY[0]*f[3]*m_*wv2-1.0*BcurvY[0]*Phi[3]*f[11]*q_-1.0*BcurvY[1]*Phi[2]*f[11]*q_-1.0*BcurvY[1]*Phi[3]*f[7]*q_-1.0*BcurvY[0]*Phi[2]*f[7]*q_)-3.0*(BcurvY[0]*Phi[3]*f[5]+BcurvY[1]*Phi[2]*f[5]+BcurvY[1]*f[2]*Phi[3]+BcurvY[0]*Phi[2]*f[2])*dfac_m*dfac_v2*dfac_y*q_*wv+2.0*(BcurvY[1]*f[1]+BcurvY[0]*f[0])*dfac_m*dfac_y*m_*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*(Bmag[1]*BmagInv[1]*f[6]*wm+BmagInv[0]*Bmag[1]*f[3]*wm+BmagInv[1]*Phi[3]*f[11]*q_+BmagInv[0]*Phi[3]*f[7]*q_+BmagInv[1]*Phi[1]*f[6]*q_+BmagInv[0]*Phi[1]*f[3]*q_)+1.732050807568877*Bmag[1]*(BmagInv[1]*f[13]+BmagInv[0]*f[10])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[8] += -0.75*(BmagInv[0]*Phi[3]*f[8]+BmagInv[1]*Phi[2]*f[8]+BmagInv[1]*Phi[3]*f[4]+BmagInv[0]*Phi[2]*f[4])*dfac_x*dfac_y; 
  out[9] += (0.25*(3.464101615137754*(BcurvY[1]*f[8]+BcurvY[0]*f[4])*dfac_m*dfac_v*dfac_y*m_*wv2+2.0*(BcurvY[1]*f[13]+BcurvY[0]*f[10])*dfac_m*dfac_y*m_*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*(Bmag[1]*BmagInv[1]*f[8]*wm+BmagInv[0]*Bmag[1]*f[4]*wm+BmagInv[1]*Phi[3]*f[12]*q_+BmagInv[0]*Phi[3]*f[9]*q_+BmagInv[1]*Phi[1]*f[8]*q_+BmagInv[0]*Phi[1]*f[4]*q_)+1.732050807568877*Bmag[1]*(BmagInv[1]*f[1]+BmagInv[0]*f[0])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[10] += -0.25*(3.0*(BcurvY[0]*Phi[3]*f[8]+BcurvY[1]*Phi[2]*f[8]+BcurvY[1]*Phi[3]*f[4]+BcurvY[0]*Phi[2]*f[4])*dfac_v*dfac_y*wv+1.732050807568877*(BcurvY[0]*Phi[3]*f[13]+BcurvY[1]*Phi[2]*f[13]+BcurvY[1]*Phi[3]*f[10]+BcurvY[0]*Phi[2]*f[10])*dfac_y); 
  out[11] += (0.05*(1.732050807568877*dfac_m*dfac_v*dfac_y*(10.0*BcurvY[0]*f[6]*m_*wv2+10.0*BcurvY[1]*f[3]*m_*wv2-9.0*BcurvY[1]*Phi[3]*f[11]*q_-5.0*BcurvY[0]*Phi[2]*f[11]*q_-5.0*BcurvY[0]*Phi[3]*f[7]*q_-5.0*BcurvY[1]*Phi[2]*f[7]*q_)-3.0*(9.0*BcurvY[1]*Phi[3]*f[5]+5.0*BcurvY[0]*Phi[2]*f[5]+5.0*BcurvY[0]*f[2]*Phi[3]+5.0*BcurvY[1]*Phi[2]*f[2])*dfac_m*dfac_v2*dfac_y*q_*wv+10.0*(BcurvY[0]*f[1]+f[0]*BcurvY[1])*dfac_m*dfac_y*m_*wv+15.0*dfac_m*dfac_v*dfac_x*dfac_y*(BmagInv[0]*Bmag[1]*f[6]*wm+Bmag[1]*BmagInv[1]*f[3]*wm-1.0*BmagInv[1]*Phi[2]*f[11]*q_-1.0*BmagInv[0]*Phi[2]*f[7]*q_+BmagInv[0]*Phi[1]*f[6]*q_+BmagInv[1]*Phi[1]*f[3]*q_)+8.660254037844386*Bmag[1]*(BmagInv[0]*f[13]+BmagInv[1]*f[10])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[12] += (0.25*(3.464101615137754*(BcurvY[0]*f[8]+BcurvY[1]*f[4])*dfac_m*dfac_v*dfac_y*m_*wv2+2.0*(BcurvY[0]*f[13]+BcurvY[1]*f[10])*dfac_m*dfac_y*m_*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*(BmagInv[0]*Bmag[1]*f[8]*wm+Bmag[1]*BmagInv[1]*f[4]*wm-1.0*BmagInv[1]*Phi[2]*f[12]*q_-1.0*BmagInv[0]*Phi[2]*f[9]*q_+BmagInv[0]*Phi[1]*f[8]*q_+BmagInv[1]*Phi[1]*f[4]*q_)+1.732050807568877*Bmag[1]*(BmagInv[0]*f[1]+f[0]*BmagInv[1])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[13] += -0.05*(3.0*(9.0*BcurvY[1]*Phi[3]*f[8]+5.0*BcurvY[0]*Phi[2]*f[8]+5.0*BcurvY[0]*Phi[3]*f[4]+5.0*BcurvY[1]*Phi[2]*f[4])*dfac_v*dfac_y*wv+15.0*(BmagInv[0]*Phi[3]*f[13]+BmagInv[1]*Phi[2]*f[13]+BmagInv[1]*Phi[3]*f[10]+BmagInv[0]*Phi[2]*f[10])*dfac_x*dfac_y+1.732050807568877*(9.0*BcurvY[1]*Phi[3]*f[13]+5.0*BcurvY[0]*Phi[2]*f[13]+5.0*BcurvY[0]*Phi[3]*f[10]+5.0*BcurvY[1]*Phi[2]*f[10])*dfac_y); 
  out[14] += (0.25*(1.732050807568877*dfac_m*dfac_v*dfac_y*(2.0*BcurvY[1]*f[13]*m_*wv2+2.0*BcurvY[0]*f[10]*m_*wv2-1.0*BcurvY[0]*Phi[3]*f[15]*q_-1.0*BcurvY[1]*Phi[2]*f[15]*q_-1.0*BcurvY[1]*Phi[3]*f[14]*q_-1.0*BcurvY[0]*Phi[2]*f[14]*q_)-3.0*(BcurvY[0]*Phi[3]*f[12]+BcurvY[1]*Phi[2]*f[12]+BcurvY[1]*Phi[3]*f[9]+BcurvY[0]*Phi[2]*f[9])*dfac_m*dfac_v2*dfac_y*q_*wv+2.0*(BcurvY[1]*f[8]+BcurvY[0]*f[4])*dfac_m*dfac_y*m_*wv+3.0*dfac_m*dfac_v*dfac_x*dfac_y*(Bmag[1]*BmagInv[1]*f[13]*wm+BmagInv[0]*Bmag[1]*f[10]*wm+BmagInv[1]*Phi[3]*f[15]*q_+BmagInv[0]*Phi[3]*f[14]*q_+BmagInv[1]*Phi[1]*f[13]*q_+BmagInv[0]*Phi[1]*f[10]*q_)+1.732050807568877*Bmag[1]*(BmagInv[1]*f[6]+BmagInv[0]*f[3])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  out[15] += (0.05*(1.732050807568877*dfac_m*dfac_v*dfac_y*(10.0*BcurvY[0]*f[13]*m_*wv2+10.0*BcurvY[1]*f[10]*m_*wv2-9.0*BcurvY[1]*Phi[3]*f[15]*q_-5.0*BcurvY[0]*Phi[2]*f[15]*q_-5.0*BcurvY[0]*Phi[3]*f[14]*q_-5.0*BcurvY[1]*Phi[2]*f[14]*q_)-3.0*(9.0*BcurvY[1]*Phi[3]*f[12]+5.0*BcurvY[0]*Phi[2]*f[12]+5.0*BcurvY[0]*Phi[3]*f[9]+5.0*BcurvY[1]*Phi[2]*f[9])*dfac_m*dfac_v2*dfac_y*q_*wv+10.0*(BcurvY[0]*f[8]+BcurvY[1]*f[4])*dfac_m*dfac_y*m_*wv+15.0*dfac_m*dfac_v*dfac_x*dfac_y*(BmagInv[0]*Bmag[1]*f[13]*wm+Bmag[1]*BmagInv[1]*f[10]*wm-1.0*BmagInv[1]*Phi[2]*f[15]*q_-1.0*BmagInv[0]*Phi[2]*f[14]*q_+BmagInv[0]*Phi[1]*f[13]*q_+BmagInv[1]*Phi[1]*f[10]*q_)+8.660254037844386*Bmag[1]*(BmagInv[0]*f[6]+BmagInv[1]*f[3])*dfac_v*dfac_x*dfac_y))/(dfac_m*dfac_v*q_); 
  double cflFreq = 0.0; 
  cflFreq += fabs(-0.2165063509461096*Bmag[0]*BmagInv[0]*Phi[2]*dfac_y)*dxInv; 
  cflFreq += fabs((0.125*(2.0*BcurvY[0]*Bmag[0]*m_*wv2+1.732050807568877*Bmag[0]*BmagInv[0]*Bmag[1]*dfac_x*wm+1.732050807568877*Bmag[0]*BmagInv[0]*Phi[1]*dfac_x*q_))/q_)*dyInv; 
  cflFreq += fabs(-0.2165063509461096*BcurvY[0]*Bmag[0]*Phi[2]*dfac_y*wv)*dvInv; 
  return cflFreq; 
} 
