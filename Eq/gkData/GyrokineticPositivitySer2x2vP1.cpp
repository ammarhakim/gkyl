#include <GyrokineticModDecl.h> 
double GyrokineticVolPositivity2x2vSer_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double *f, double *out, double *cflRateByDir) 
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
  double m2 = m_*m_; 
  double cflRate = 0.0; 
  cflRateByDir[0] = 0.; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaCtrl;
  double alphax[16]; 
  alphax[0] = dfac_x*((2.0*BdriftX[0]*m_*wv2)/q_-1.732050807568877*BmagInv[0]*Phi[2]*dfac_y); 
  alphax[1] = -1.732050807568877*BmagInv[0]*Phi[3]*dfac_x*dfac_y; 
  alphax[3] = (1.154700538379252*BdriftX[0]*dfac_x*m_*wv)/(dfac_v*q_); 
  cflRateByDir[1] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*((-0.25*alphax[3])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphax[3])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphax[3])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphax[3])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.25*alphax[3])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphax[3])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphax[3])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphax[3])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[1];
  out[1] += 0.4330127018922193*(alphax[3]*f[3]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[5] += 0.4330127018922193*(alphax[3]*f[7]+alphax[1]*f[5]+alphax[0]*f[2]); 
  out[6] += 0.4330127018922193*(alphax[1]*f[6]+alphax[0]*f[3]+f[0]*alphax[3]); 
  out[8] += 0.4330127018922193*(alphax[3]*f[10]+alphax[1]*f[8]+alphax[0]*f[4]); 
  out[11] += 0.4330127018922193*(alphax[1]*f[11]+alphax[0]*f[7]+f[2]*alphax[3]); 
  out[12] += 0.4330127018922193*(alphax[3]*f[14]+alphax[1]*f[12]+alphax[0]*f[9]); 
  out[13] += 0.4330127018922193*(alphax[1]*f[13]+alphax[0]*f[10]+alphax[3]*f[4]); 
  out[15] += 0.4330127018922193*(alphax[1]*f[15]+alphax[0]*f[14]+alphax[3]*f[9]); 
  double alphay[16]; 
  alphay[0] = dfac_y*((2.0*BdriftY[0]*m_*wv2)/q_+1.732050807568877*BmagInv[0]*Phi[1]*dfac_x); 
  alphay[2] = 1.732050807568877*BmagInv[0]*Phi[3]*dfac_x*dfac_y; 
  alphay[3] = (1.154700538379252*BdriftY[0]*dfac_y*m_*wv)/(dfac_v*q_); 
  cflRateByDir[2] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*((-0.25*alphay[3])-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[3])-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[3])-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[3])-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.25*alphay[3])+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[3])+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[3])+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[3])+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[2];
  out[2] += 0.4330127018922193*(alphay[3]*f[3]+alphay[2]*f[2]+alphay[0]*f[0]); 
  out[5] += 0.4330127018922193*(alphay[3]*f[6]+alphay[2]*f[5]+alphay[0]*f[1]); 
  out[7] += 0.4330127018922193*(alphay[2]*f[7]+alphay[0]*f[3]+f[0]*alphay[3]); 
  out[9] += 0.4330127018922193*(alphay[3]*f[10]+alphay[2]*f[9]+alphay[0]*f[4]); 
  out[11] += 0.4330127018922193*(alphay[2]*f[11]+alphay[0]*f[6]+f[1]*alphay[3]); 
  out[12] += 0.4330127018922193*(alphay[3]*f[13]+alphay[2]*f[12]+alphay[0]*f[8]); 
  out[14] += 0.4330127018922193*(alphay[2]*f[14]+alphay[0]*f[10]+alphay[3]*f[4]); 
  out[15] += 0.4330127018922193*(alphay[2]*f[15]+alphay[0]*f[13]+alphay[3]*f[8]); 
  double alphav[16]; 
  alphav[0] = -1.732050807568877*dfac_v*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*wv; 
  alphav[1] = -1.732050807568877*BdriftY[0]*Phi[3]*dfac_v*dfac_y*wv; 
  alphav[2] = -1.732050807568877*BdriftX[0]*Phi[3]*dfac_v*dfac_x*wv; 
  alphav[3] = -1.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x); 
  alphav[6] = -1.0*BdriftY[0]*Phi[3]*dfac_y; 
  alphav[7] = -1.0*BdriftX[0]*Phi[3]*dfac_x; 
  cflRateByDir[3] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphav[3]-1.0*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphav[3]+alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.4330127018922193*(alphav[7]+alphav[6])-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphav[7]-0.4330127018922193*(alphav[6]+alphav[3])-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphav[7])+0.4330127018922193*alphav[6]-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphav[2]+alphav[1]+alphav[0])-0.4330127018922193*(alphav[7]+alphav[6]+alphav[3])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(alphav[7]+alphav[6])-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphav[7]-0.4330127018922193*(alphav[6]+alphav[3])-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphav[7])+0.4330127018922193*alphav[6]-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphav[2]+alphav[1]+alphav[0])-0.4330127018922193*(alphav[7]+alphav[6]+alphav[3])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.4330127018922193*(alphav[7]+alphav[6]))+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphav[7])+0.4330127018922193*(alphav[6]+alphav[3])-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphav[7]+alphav[6]+alphav[3])+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(alphav[7]+alphav[6]))+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphav[7])+0.4330127018922193*(alphav[6]+alphav[3])-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphav[7]+alphav[6]+alphav[3])+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[3];
  out[3] += 0.4330127018922193*(alphav[7]*f[7]+alphav[6]*f[6]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[6] += 0.4330127018922193*(alphav[7]*f[11]+alphav[3]*f[6]+f[3]*alphav[6]+alphav[2]*f[5]+alphav[0]*f[1]+f[0]*alphav[1]); 
  out[7] += 0.4330127018922193*(alphav[6]*f[11]+alphav[3]*f[7]+f[3]*alphav[7]+alphav[1]*f[5]+alphav[0]*f[2]+f[0]*alphav[2]); 
  out[10] += 0.4330127018922193*(alphav[7]*f[14]+alphav[6]*f[13]+alphav[3]*f[10]+alphav[2]*f[9]+alphav[1]*f[8]+alphav[0]*f[4]); 
  out[11] += 0.4330127018922193*(alphav[3]*f[11]+alphav[6]*f[7]+f[6]*alphav[7]+alphav[0]*f[5]+alphav[1]*f[2]+f[1]*alphav[2]); 
  out[13] += 0.4330127018922193*(alphav[7]*f[15]+alphav[3]*f[13]+alphav[2]*f[12]+alphav[6]*f[10]+alphav[0]*f[8]+alphav[1]*f[4]); 
  out[14] += 0.4330127018922193*(alphav[6]*f[15]+alphav[3]*f[14]+alphav[1]*f[12]+alphav[7]*f[10]+alphav[0]*f[9]+alphav[2]*f[4]); 
  out[15] += 0.4330127018922193*(alphav[3]*f[15]+alphav[6]*f[14]+alphav[7]*f[13]+alphav[0]*f[12]+alphav[1]*f[9]+alphav[2]*f[8]); 
  return cflRate; 
} 
double GyrokineticVolPositivity2x2vSer_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double *f, double *out, double *cflRateByDir) 
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
  double m2 = m_*m_; 
  double cflRate = 0.0; 
  cflRateByDir[0] = 0.; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaCtrl;
  double alphax[16]; 
  alphax[0] = dfac_x*((2.0*BdriftX[0]*m_*wv2)/q_-1.732050807568877*(BmagInv[1]*Phi[3]+BmagInv[0]*Phi[2])*dfac_y); 
  alphax[1] = dfac_x*((2.0*BdriftX[1]*m_*wv2)/q_-1.732050807568877*(BmagInv[0]*Phi[3]+BmagInv[1]*Phi[2])*dfac_y); 
  alphax[3] = (1.154700538379252*BdriftX[0]*dfac_x*m_*wv)/(dfac_v*q_); 
  alphax[6] = (1.154700538379252*BdriftX[1]*dfac_x*m_*wv)/(dfac_v*q_); 
  cflRateByDir[1] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.4330127018922193*alphax[6]-0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[6]-0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[6])+0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[6])+0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[6]-0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[6]-0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[6])+0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[6])+0.25*alphax[3]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.4330127018922193*alphax[6])-0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[6])-0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[6]+0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[6]+0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[6])-0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[6])-0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[6]+0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[6]+0.25*alphax[3]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[1];
  out[1] += 0.4330127018922193*(alphax[6]*f[6]+alphax[3]*f[3]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[5] += 0.4330127018922193*(alphax[6]*f[11]+alphax[3]*f[7]+alphax[1]*f[5]+alphax[0]*f[2]); 
  out[6] += 0.4330127018922193*(alphax[1]*f[6]+f[1]*alphax[6]+alphax[0]*f[3]+f[0]*alphax[3]); 
  out[8] += 0.4330127018922193*(alphax[6]*f[13]+alphax[3]*f[10]+alphax[1]*f[8]+alphax[0]*f[4]); 
  out[11] += 0.4330127018922193*(alphax[1]*f[11]+alphax[0]*f[7]+f[5]*alphax[6]+f[2]*alphax[3]); 
  out[12] += 0.4330127018922193*(alphax[6]*f[15]+alphax[3]*f[14]+alphax[1]*f[12]+alphax[0]*f[9]); 
  out[13] += 0.4330127018922193*(alphax[1]*f[13]+alphax[0]*f[10]+alphax[6]*f[8]+alphax[3]*f[4]); 
  out[15] += 0.4330127018922193*(alphax[1]*f[15]+alphax[0]*f[14]+alphax[6]*f[12]+alphax[3]*f[9]); 
  double alphay[16]; 
  alphay[0] = dfac_y*((2.0*BdriftY[0]*m_*wv2)/q_+1.732050807568877*BmagInv[0]*dfac_x*((Bmag[1]*wm)/q_+Phi[1])); 
  alphay[1] = dfac_y*((2.0*BdriftY[1]*m_*wv2)/q_+1.732050807568877*BmagInv[1]*dfac_x*((Bmag[1]*wm)/q_+Phi[1])); 
  alphay[2] = 1.732050807568877*BmagInv[0]*Phi[3]*dfac_x*dfac_y; 
  alphay[3] = (1.154700538379252*BdriftY[0]*dfac_y*m_*wv)/(dfac_v*q_); 
  alphay[4] = (BmagInv[0]*Bmag[1]*dfac_x*dfac_y)/(dfac_m*q_); 
  alphay[5] = 1.732050807568877*BmagInv[1]*Phi[3]*dfac_x*dfac_y; 
  alphay[6] = (1.154700538379252*BdriftY[1]*dfac_y*m_*wv)/(dfac_v*q_); 
  alphay[8] = (Bmag[1]*BmagInv[1]*dfac_x*dfac_y)/(dfac_m*q_); 
  cflRateByDir[2] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[2];
  out[2] += 0.4330127018922193*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[5] += 0.4330127018922193*(alphay[4]*f[8]+f[4]*alphay[8]+alphay[3]*f[6]+f[3]*alphay[6]+alphay[2]*f[5]+f[2]*alphay[5]+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[7] += 0.4330127018922193*(alphay[8]*f[13]+alphay[5]*f[11]+alphay[4]*f[10]+alphay[2]*f[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphay[0]*f[3]+f[0]*alphay[3]); 
  out[9] += 0.4330127018922193*(alphay[6]*f[13]+alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4]); 
  out[11] += 0.4330127018922193*(alphay[4]*f[13]+alphay[2]*f[11]+alphay[8]*f[10]+alphay[5]*f[7]+alphay[0]*f[6]+f[0]*alphay[6]+alphay[1]*f[3]+f[1]*alphay[3]); 
  out[12] += 0.4330127018922193*(alphay[3]*f[13]+alphay[2]*f[12]+alphay[6]*f[10]+alphay[5]*f[9]+alphay[0]*f[8]+f[0]*alphay[8]+alphay[1]*f[4]+f[1]*alphay[4]); 
  out[14] += 0.4330127018922193*(alphay[5]*f[15]+alphay[2]*f[14]+alphay[1]*f[13]+alphay[0]*f[10]+alphay[6]*f[8]+f[6]*alphay[8]+alphay[3]*f[4]+f[3]*alphay[4]); 
  out[15] += 0.4330127018922193*(alphay[2]*f[15]+alphay[5]*f[14]+alphay[0]*f[13]+alphay[1]*f[10]+alphay[3]*f[8]+f[3]*alphay[8]+alphay[4]*f[6]+f[4]*alphay[6]); 
  double alphav[16]; 
  alphav[0] = -1.732050807568877*dfac_v*((BdriftX[0]*Bmag[1]*dfac_x*wm)/q_+(BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*wv; 
  alphav[1] = -1.732050807568877*dfac_v*((BdriftX[1]*Bmag[1]*dfac_x*wm)/q_+(BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*wv; 
  alphav[2] = -1.732050807568877*BdriftX[0]*Phi[3]*dfac_v*dfac_x*wv; 
  alphav[3] = -1.0*((BdriftX[0]*Bmag[1]*dfac_x*wm)/q_+(BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x); 
  alphav[4] = -(1.0*BdriftX[0]*Bmag[1]*dfac_v*dfac_x*wv)/(dfac_m*q_); 
  alphav[5] = -1.732050807568877*BdriftX[1]*Phi[3]*dfac_v*dfac_x*wv; 
  alphav[6] = -1.0*((BdriftX[1]*Bmag[1]*dfac_x*wm)/q_+(BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x); 
  alphav[7] = -1.0*BdriftX[0]*Phi[3]*dfac_x; 
  alphav[8] = -(1.0*BdriftX[1]*Bmag[1]*dfac_v*dfac_x*wv)/(dfac_m*q_); 
  alphav[10] = -(0.5773502691896258*BdriftX[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alphav[11] = -1.0*BdriftX[1]*Phi[3]*dfac_x; 
  alphav[13] = -(0.5773502691896258*BdriftX[1]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  cflRateByDir[3] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphav[3]-1.0*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphav[3]+alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*((-0.4330127018922193*(alphav[13]+alphav[11]))+0.4330127018922193*alphav[10]+0.25*alphav[8]+0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]-0.25*alphav[4]-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(alphav[13]+alphav[11]+alphav[10])-0.25*alphav[8]+0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphav[13])+0.4330127018922193*(alphav[11]+alphav[10])+0.25*alphav[8]-0.4330127018922193*alphav[7]+0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphav[13]-0.4330127018922193*alphav[11]+0.4330127018922193*alphav[10]-0.25*alphav[8]-0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]-0.25*alphav[4]-0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphav[13]-0.4330127018922193*(alphav[11]+alphav[10])-0.25*alphav[8]+0.4330127018922193*(alphav[7]+alphav[6])+0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphav[13])+0.4330127018922193*alphav[11]-0.4330127018922193*alphav[10]+0.25*alphav[8]+0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]-0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(alphav[13]+alphav[11])-0.4330127018922193*alphav[10]-0.25*alphav[8]-0.4330127018922193*alphav[7]+0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(alphav[13]+alphav[11]+alphav[10]))+0.25*alphav[8]-0.4330127018922193*(alphav[7]+alphav[6])+0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.4330127018922193*(alphav[13]+alphav[11])-0.4330127018922193*alphav[10]+0.25*alphav[8]-0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]-0.25*alphav[4]+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(alphav[13]+alphav[11]+alphav[10]))-0.25*alphav[8]-0.4330127018922193*alphav[7]+0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphav[13]-0.4330127018922193*(alphav[11]+alphav[10])+0.25*alphav[8]+0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphav[13])+0.4330127018922193*alphav[11]-0.4330127018922193*alphav[10]-0.25*alphav[8]+0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]-0.25*alphav[4]+0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphav[13])+0.4330127018922193*(alphav[11]+alphav[10])-0.25*alphav[8]-0.4330127018922193*(alphav[7]+alphav[6])+0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphav[13]-0.4330127018922193*alphav[11]+0.4330127018922193*alphav[10]+0.25*alphav[8]-0.4330127018922193*alphav[7]+0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]+0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(alphav[13]+alphav[11]))+0.4330127018922193*alphav[10]-0.25*alphav[8]+0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphav[13]+alphav[11]+alphav[10])+0.25*alphav[8]+0.4330127018922193*(alphav[7]+alphav[6])+0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[3];
  out[3] += 0.4330127018922193*(alphav[13]*f[13]+alphav[11]*f[11]+alphav[10]*f[10]+alphav[8]*f[8]+alphav[7]*f[7]+alphav[6]*f[6]+alphav[5]*f[5]+alphav[4]*f[4]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[6] += 0.4330127018922193*(alphav[10]*f[13]+f[10]*alphav[13]+alphav[7]*f[11]+f[7]*alphav[11]+alphav[4]*f[8]+f[4]*alphav[8]+alphav[3]*f[6]+f[3]*alphav[6]+alphav[2]*f[5]+f[2]*alphav[5]+alphav[0]*f[1]+f[0]*alphav[1]); 
  out[7] += 0.4330127018922193*(alphav[13]*f[15]+alphav[10]*f[14]+alphav[8]*f[12]+alphav[6]*f[11]+f[6]*alphav[11]+alphav[4]*f[9]+alphav[3]*f[7]+f[3]*alphav[7]+alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[2]+f[0]*alphav[2]); 
  out[10] += 0.4330127018922193*(alphav[11]*f[15]+alphav[7]*f[14]+alphav[6]*f[13]+f[6]*alphav[13]+alphav[5]*f[12]+alphav[3]*f[10]+f[3]*alphav[10]+alphav[2]*f[9]+alphav[1]*f[8]+f[1]*alphav[8]+alphav[0]*f[4]+f[0]*alphav[4]); 
  out[11] += 0.4330127018922193*(alphav[10]*f[15]+alphav[13]*f[14]+alphav[4]*f[12]+alphav[3]*f[11]+f[3]*alphav[11]+alphav[8]*f[9]+alphav[6]*f[7]+f[6]*alphav[7]+alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[2]+f[1]*alphav[2]); 
  out[13] += 0.4330127018922193*(alphav[7]*f[15]+alphav[11]*f[14]+alphav[3]*f[13]+f[3]*alphav[13]+alphav[2]*f[12]+alphav[6]*f[10]+f[6]*alphav[10]+alphav[5]*f[9]+alphav[0]*f[8]+f[0]*alphav[8]+alphav[1]*f[4]+f[1]*alphav[4]); 
  out[14] += 0.4330127018922193*(alphav[6]*f[15]+alphav[3]*f[14]+alphav[11]*f[13]+f[11]*alphav[13]+alphav[1]*f[12]+alphav[7]*f[10]+f[7]*alphav[10]+alphav[0]*f[9]+alphav[5]*f[8]+f[5]*alphav[8]+alphav[2]*f[4]+f[2]*alphav[4]); 
  out[15] += 0.4330127018922193*(alphav[3]*f[15]+alphav[6]*f[14]+alphav[7]*f[13]+f[7]*alphav[13]+alphav[0]*f[12]+alphav[10]*f[11]+f[10]*alphav[11]+alphav[1]*f[9]+alphav[2]*f[8]+f[2]*alphav[8]+alphav[4]*f[5]+f[4]*alphav[5]); 
  return cflRate; 
} 
