#include <GyrokineticModDecl.h> 
double EmGyrokineticVolPositivity2x2vSer_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, 
                        const double *f, double *outX, double *outV, double *cflRateByDir) 
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
  alphax[0] = dfac_x*((2.0*BdriftX[0]*m_*wv2)/q_+(1.732050807568877*BmagInv[0]*Apar[2]*dfac_y+Apar[0]*BdriftX[0])*wv-1.732050807568877*BmagInv[0]*Phi[2]*dfac_y); 
  alphax[1] = dfac_x*((1.732050807568877*BmagInv[0]*Apar[3]*dfac_y+BdriftX[0]*Apar[1])*wv-1.732050807568877*BmagInv[0]*Phi[3]*dfac_y); 
  alphax[2] = BdriftX[0]*Apar[2]*dfac_x*wv; 
  alphax[3] = (1.154700538379252*BdriftX[0]*dfac_x*m_*wv)/(dfac_v*q_); 
  alphax[5] = BdriftX[0]*Apar[3]*dfac_x*wv; 
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
  alphaL = 0.0625*(0.4330127018922193*alphax[5]-0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[5])-0.25*alphax[3]+0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[5])+0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[5]-0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[5])-0.25*alphax[3]+0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[5])+0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.4330127018922193*alphax[5])-0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[5])+0.25*alphax[3]-0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[5]+0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[5])-0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[5])+0.25*alphax[3]-0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[5]+0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[1];
  outX[1] += 0.4330127018922193*(alphax[5]*f[5]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  outX[5] += 0.4330127018922193*(alphax[3]*f[7]+alphax[1]*f[5]+f[1]*alphax[5]+alphax[0]*f[2]+f[0]*alphax[2]); 
  outX[6] += 0.4330127018922193*(alphax[5]*f[11]+alphax[2]*f[7]+alphax[1]*f[6]+alphax[0]*f[3]+f[0]*alphax[3]); 
  outX[8] += 0.4330127018922193*(alphax[5]*f[12]+alphax[3]*f[10]+alphax[2]*f[9]+alphax[1]*f[8]+alphax[0]*f[4]); 
  outX[11] += 0.4330127018922193*(alphax[1]*f[11]+alphax[0]*f[7]+alphax[5]*f[6]+alphax[2]*f[3]+f[2]*alphax[3]); 
  outX[12] += 0.4330127018922193*(alphax[3]*f[14]+alphax[1]*f[12]+alphax[0]*f[9]+alphax[5]*f[8]+alphax[2]*f[4]); 
  outX[13] += 0.4330127018922193*(alphax[5]*f[15]+alphax[2]*f[14]+alphax[1]*f[13]+alphax[0]*f[10]+alphax[3]*f[4]); 
  outX[15] += 0.4330127018922193*(alphax[1]*f[15]+alphax[0]*f[14]+alphax[5]*f[13]+alphax[2]*f[10]+alphax[3]*f[9]); 
  double alphay[16]; 
  alphay[0] = dfac_y*((2.0*BdriftY[0]*m_*wv2)/q_+(Apar[0]*BdriftY[0]-1.732050807568877*BmagInv[0]*Apar[1]*dfac_x)*wv+1.732050807568877*BmagInv[0]*Phi[1]*dfac_x); 
  alphay[1] = BdriftY[0]*Apar[1]*dfac_y*wv; 
  alphay[2] = -1.0*dfac_y*(1.732050807568877*BmagInv[0]*Apar[3]*dfac_x*wv-1.0*(BdriftY[0]*Apar[2]*wv+1.732050807568877*BmagInv[0]*Phi[3]*dfac_x)); 
  alphay[3] = (1.154700538379252*BdriftY[0]*dfac_y*m_*wv)/(dfac_v*q_); 
  alphay[5] = BdriftY[0]*Apar[3]*dfac_y*wv; 
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
  alphaL = 0.0625*(0.4330127018922193*alphay[5]-0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphay[5])-0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphay[5]+0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphay[5])+0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphay[5]-0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphay[5])-0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphay[5]+0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphay[5])+0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.4330127018922193*alphay[5])-0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[5]-0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphay[5])+0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[5]+0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphay[5])-0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[5]-0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphay[5])+0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[5]+0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflRateByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[2];
  outX[2] += 0.4330127018922193*(alphay[5]*f[5]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  outX[5] += 0.4330127018922193*(alphay[3]*f[6]+alphay[2]*f[5]+f[2]*alphay[5]+alphay[0]*f[1]+f[0]*alphay[1]); 
  outX[7] += 0.4330127018922193*(alphay[5]*f[11]+alphay[2]*f[7]+alphay[1]*f[6]+alphay[0]*f[3]+f[0]*alphay[3]); 
  outX[9] += 0.4330127018922193*(alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+alphay[0]*f[4]); 
  outX[11] += 0.4330127018922193*(alphay[2]*f[11]+alphay[5]*f[7]+alphay[0]*f[6]+alphay[1]*f[3]+f[1]*alphay[3]); 
  outX[12] += 0.4330127018922193*(alphay[3]*f[13]+alphay[2]*f[12]+alphay[5]*f[9]+alphay[0]*f[8]+alphay[1]*f[4]); 
  outX[14] += 0.4330127018922193*(alphay[5]*f[15]+alphay[2]*f[14]+alphay[1]*f[13]+alphay[0]*f[10]+alphay[3]*f[4]); 
  outX[15] += 0.4330127018922193*(alphay[2]*f[15]+alphay[5]*f[14]+alphay[0]*f[13]+alphay[1]*f[10]+alphay[3]*f[8]); 
  double alphav[16]; 
  alphav[0] = -1.0*dfac_v*(1.732050807568877*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*wv+((0.8660254037844386*(BdriftY[0]*(Apar[1]*Phi[3]+Apar[0]*Phi[2])*dfac_y+BdriftX[0]*(Apar[2]*Phi[3]+Apar[0]*Phi[1])*dfac_x)-1.0*BmagInv[0]*(1.5*Apar[1]*Phi[2]-1.5*Phi[1]*Apar[2])*dfac_x*dfac_y)*q_)/m_); 
  alphav[1] = -1.0*dfac_v*(1.732050807568877*BdriftY[0]*Phi[3]*dfac_y*wv+((0.8660254037844386*(BdriftY[0]*(Apar[0]*Phi[3]+Apar[1]*Phi[2])*dfac_y+BdriftX[0]*(Apar[3]*Phi[3]+Apar[1]*Phi[1])*dfac_x)-1.0*BmagInv[0]*(1.5*Apar[1]*Phi[3]-1.5*Phi[1]*Apar[3])*dfac_x*dfac_y)*q_)/m_); 
  alphav[2] = -1.0*dfac_v*(1.732050807568877*BdriftX[0]*Phi[3]*dfac_x*wv+((0.8660254037844386*(BdriftY[0]*(Apar[3]*Phi[3]+Apar[2]*Phi[2])*dfac_y+BdriftX[0]*(Apar[0]*Phi[3]+Phi[1]*Apar[2])*dfac_x)-1.0*BmagInv[0]*(1.5*Phi[2]*Apar[3]-1.5*Apar[2]*Phi[3])*dfac_x*dfac_y)*q_)/m_); 
  alphav[3] = -1.0*(BdriftY[0]*Phi[2]*dfac_y+BdriftX[0]*Phi[1]*dfac_x); 
  alphav[5] = -(0.8660254037844386*dfac_v*(BdriftY[0]*(Apar[2]*Phi[3]+Phi[2]*Apar[3])*dfac_y+BdriftX[0]*(Apar[1]*Phi[3]+Phi[1]*Apar[3])*dfac_x)*q_)/m_; 
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
  alphaL = 0.0625*(0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*alphav[5]-0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphav[7])+0.4330127018922193*alphav[6]-0.25*alphav[5]-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(alphav[7]+alphav[6]))+0.25*alphav[5]-0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*alphav[5]-0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphav[7])+0.4330127018922193*alphav[6]-0.25*alphav[5]-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(alphav[7]+alphav[6]))+0.25*alphav[5]-0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.4330127018922193*(alphav[7]+alphav[6]))+0.25*alphav[5]+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphav[7])+0.4330127018922193*alphav[6]-0.25*alphav[5]+0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*alphav[5]+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]+0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(alphav[7]+alphav[6]))+0.25*alphav[5]+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphav[7])+0.4330127018922193*alphav[6]-0.25*alphav[5]+0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*alphav[5]+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]+0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[3];
  outV[3] += 0.4330127018922193*(alphav[7]*f[7]+alphav[6]*f[6]+alphav[5]*f[5]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  outV[6] += 0.4330127018922193*(alphav[7]*f[11]+alphav[3]*f[6]+f[3]*alphav[6]+alphav[2]*f[5]+f[2]*alphav[5]+alphav[0]*f[1]+f[0]*alphav[1]); 
  outV[7] += 0.4330127018922193*(alphav[6]*f[11]+alphav[3]*f[7]+f[3]*alphav[7]+alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[2]+f[0]*alphav[2]); 
  outV[10] += 0.4330127018922193*(alphav[7]*f[14]+alphav[6]*f[13]+alphav[5]*f[12]+alphav[3]*f[10]+alphav[2]*f[9]+alphav[1]*f[8]+alphav[0]*f[4]); 
  outV[11] += 0.4330127018922193*(alphav[3]*f[11]+alphav[6]*f[7]+f[6]*alphav[7]+alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[2]+f[1]*alphav[2]); 
  outV[13] += 0.4330127018922193*(alphav[7]*f[15]+alphav[3]*f[13]+alphav[2]*f[12]+alphav[6]*f[10]+alphav[5]*f[9]+alphav[0]*f[8]+alphav[1]*f[4]); 
  outV[14] += 0.4330127018922193*(alphav[6]*f[15]+alphav[3]*f[14]+alphav[1]*f[12]+alphav[7]*f[10]+alphav[0]*f[9]+alphav[5]*f[8]+alphav[2]*f[4]); 
  outV[15] += 0.4330127018922193*(alphav[3]*f[15]+alphav[6]*f[14]+alphav[7]*f[13]+alphav[0]*f[12]+alphav[1]*f[9]+alphav[2]*f[8]+f[4]*alphav[5]); 
  return cflRate; 
} 
double EmGyrokineticStep2VolPositivity2x2vSer_P1(const double q_, const double m_, const double *w, const double *dxv, 
                           const double *ohmMod, const double *dApardt, 
                           const double *f, double *outVol, double *outSurf, double *cflRateByDir) 
{ 
  double dvInv = 1.0/dxv[2]; 
  double dfac_v = 2.0/dxv[2]; 
  outVol[3] += -(0.8660254037844386*(dApardt[3]*f[5]+dApardt[2]*f[2]+dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  outVol[6] += -(0.8660254037844386*(dApardt[2]*f[5]+f[2]*dApardt[3]+dApardt[0]*f[1]+f[0]*dApardt[1])*dfac_v*q_)/m_; 
  outVol[7] += -(0.8660254037844386*(dApardt[1]*f[5]+f[1]*dApardt[3]+dApardt[0]*f[2]+f[0]*dApardt[2])*dfac_v*q_)/m_; 
  outVol[10] += -(0.8660254037844386*(dApardt[3]*f[12]+dApardt[2]*f[9]+dApardt[1]*f[8]+dApardt[0]*f[4])*dfac_v*q_)/m_; 
  outVol[11] += -(0.8660254037844386*(dApardt[0]*f[5]+f[0]*dApardt[3]+dApardt[1]*f[2]+f[1]*dApardt[2])*dfac_v*q_)/m_; 
  outVol[13] += -(0.8660254037844386*(dApardt[2]*f[12]+dApardt[3]*f[9]+dApardt[0]*f[8]+dApardt[1]*f[4])*dfac_v*q_)/m_; 
  outVol[14] += -(0.8660254037844386*(dApardt[1]*f[12]+dApardt[0]*f[9]+dApardt[3]*f[8]+dApardt[2]*f[4])*dfac_v*q_)/m_; 
  outVol[15] += -(0.8660254037844386*(dApardt[0]*f[12]+dApardt[1]*f[9]+dApardt[2]*f[8]+dApardt[3]*f[4])*dfac_v*q_)/m_; 
  outSurf[0] += -(0.5*(dApardt[3]*ohmMod[5]+dApardt[2]*ohmMod[2]+dApardt[1]*ohmMod[1]+dApardt[0]*ohmMod[0])*q_)/m_; 
  outSurf[1] += -(0.5*(dApardt[2]*ohmMod[5]+ohmMod[2]*dApardt[3]+dApardt[0]*ohmMod[1]+ohmMod[0]*dApardt[1])*q_)/m_; 
  outSurf[2] += -(0.5*(dApardt[1]*ohmMod[5]+ohmMod[1]*dApardt[3]+dApardt[0]*ohmMod[2]+ohmMod[0]*dApardt[2])*q_)/m_; 
  outSurf[3] += -(0.5*(dApardt[3]*ohmMod[11]+dApardt[2]*ohmMod[7]+dApardt[1]*ohmMod[6]+dApardt[0]*ohmMod[3])*q_)/m_; 
  outSurf[4] += -(0.5*(dApardt[3]*ohmMod[12]+dApardt[2]*ohmMod[9]+dApardt[1]*ohmMod[8]+dApardt[0]*ohmMod[4])*q_)/m_; 
  outSurf[5] += -(0.5*(dApardt[0]*ohmMod[5]+ohmMod[0]*dApardt[3]+dApardt[1]*ohmMod[2]+ohmMod[1]*dApardt[2])*q_)/m_; 
  outSurf[6] += -(0.5*(dApardt[2]*ohmMod[11]+dApardt[3]*ohmMod[7]+dApardt[0]*ohmMod[6]+dApardt[1]*ohmMod[3])*q_)/m_; 
  outSurf[7] += -(0.5*(dApardt[1]*ohmMod[11]+dApardt[0]*ohmMod[7]+dApardt[3]*ohmMod[6]+dApardt[2]*ohmMod[3])*q_)/m_; 
  outSurf[8] += -(0.5*(dApardt[2]*ohmMod[12]+dApardt[3]*ohmMod[9]+dApardt[0]*ohmMod[8]+dApardt[1]*ohmMod[4])*q_)/m_; 
  outSurf[9] += -(0.5*(dApardt[1]*ohmMod[12]+dApardt[0]*ohmMod[9]+dApardt[3]*ohmMod[8]+dApardt[2]*ohmMod[4])*q_)/m_; 
  outSurf[10] += -(0.5*(dApardt[3]*ohmMod[15]+dApardt[2]*ohmMod[14]+dApardt[1]*ohmMod[13]+dApardt[0]*ohmMod[10])*q_)/m_; 
  outSurf[11] += -(0.5*(dApardt[0]*ohmMod[11]+dApardt[1]*ohmMod[7]+dApardt[2]*ohmMod[6]+dApardt[3]*ohmMod[3])*q_)/m_; 
  outSurf[12] += -(0.5*(dApardt[0]*ohmMod[12]+dApardt[1]*ohmMod[9]+dApardt[2]*ohmMod[8]+dApardt[3]*ohmMod[4])*q_)/m_; 
  outSurf[13] += -(0.5*(dApardt[2]*ohmMod[15]+dApardt[3]*ohmMod[14]+dApardt[0]*ohmMod[13]+dApardt[1]*ohmMod[10])*q_)/m_; 
  outSurf[14] += -(0.5*(dApardt[1]*ohmMod[15]+dApardt[0]*ohmMod[14]+dApardt[3]*ohmMod[13]+dApardt[2]*ohmMod[10])*q_)/m_; 
  outSurf[15] += -(0.5*(dApardt[0]*ohmMod[15]+dApardt[1]*ohmMod[14]+dApardt[2]*ohmMod[13]+dApardt[3]*ohmMod[10])*q_)/m_; 
  double cflRate = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -(0.25*dApardt[0]*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = -(0.25*dApardt[0]*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[3];return cflRate; 
} 
double EmGyrokineticVolPositivity2x2vSer_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, 
                        const double *f, double *outX, double *outV, double *cflRateByDir) 
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
  alphax[0] = dfac_x*((2.0*BdriftX[0]*m_*wv2)/q_+(1.732050807568877*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_y+Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*wv-1.732050807568877*(BmagInv[1]*Phi[3]+BmagInv[0]*Phi[2])*dfac_y); 
  alphax[1] = dfac_x*((2.0*BdriftX[1]*m_*wv2)/q_+(1.732050807568877*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y+Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*wv-1.732050807568877*(BmagInv[0]*Phi[3]+BmagInv[1]*Phi[2])*dfac_y); 
  alphax[2] = (BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*dfac_x*wv; 
  alphax[3] = (1.154700538379252*BdriftX[0]*dfac_x*m_*wv)/(dfac_v*q_); 
  alphax[5] = (BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*dfac_x*wv; 
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
  alphaL = 0.0625*(0.4330127018922193*(alphax[6]+alphax[5])-0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[6])+0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(alphax[6]+alphax[5]))+0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(alphax[6]+alphax[5])-0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[6])+0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(alphax[6]+alphax[5]))+0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*((-0.4330127018922193*(alphax[6]+alphax[5]))-0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[6])+0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphax[6]+alphax[5])+0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(alphax[6]+alphax[5]))-0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[6])+0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphax[6]+alphax[5])+0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflRateByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[1];
  outX[1] += 0.4330127018922193*(alphax[6]*f[6]+alphax[5]*f[5]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  outX[5] += 0.4330127018922193*(alphax[6]*f[11]+alphax[3]*f[7]+alphax[1]*f[5]+f[1]*alphax[5]+alphax[0]*f[2]+f[0]*alphax[2]); 
  outX[6] += 0.4330127018922193*(alphax[5]*f[11]+alphax[2]*f[7]+alphax[1]*f[6]+f[1]*alphax[6]+alphax[0]*f[3]+f[0]*alphax[3]); 
  outX[8] += 0.4330127018922193*(alphax[6]*f[13]+alphax[5]*f[12]+alphax[3]*f[10]+alphax[2]*f[9]+alphax[1]*f[8]+alphax[0]*f[4]); 
  outX[11] += 0.4330127018922193*(alphax[1]*f[11]+alphax[0]*f[7]+alphax[5]*f[6]+f[5]*alphax[6]+alphax[2]*f[3]+f[2]*alphax[3]); 
  outX[12] += 0.4330127018922193*(alphax[6]*f[15]+alphax[3]*f[14]+alphax[1]*f[12]+alphax[0]*f[9]+alphax[5]*f[8]+alphax[2]*f[4]); 
  outX[13] += 0.4330127018922193*(alphax[5]*f[15]+alphax[2]*f[14]+alphax[1]*f[13]+alphax[0]*f[10]+alphax[6]*f[8]+alphax[3]*f[4]); 
  outX[15] += 0.4330127018922193*(alphax[1]*f[15]+alphax[0]*f[14]+alphax[5]*f[13]+alphax[6]*f[12]+alphax[2]*f[10]+alphax[3]*f[9]); 
  double alphay[16]; 
  alphay[0] = dfac_y*((2.0*BdriftY[0]*m_*wv2)/q_+(Apar[1]*(BdriftY[1]-1.732050807568877*BmagInv[0]*dfac_x)+Apar[0]*BdriftY[0])*wv+1.732050807568877*BmagInv[0]*dfac_x*((Bmag[1]*wm)/q_+Phi[1])); 
  alphay[1] = dfac_y*((2.0*BdriftY[1]*m_*wv2)/q_+((-1.732050807568877*Apar[1]*BmagInv[1]*dfac_x)+Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*wv+1.732050807568877*BmagInv[1]*dfac_x*((Bmag[1]*wm)/q_+Phi[1])); 
  alphay[2] = -1.0*dfac_y*(1.732050807568877*BmagInv[0]*Apar[3]*dfac_x*wv-1.0*((BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2])*wv+1.732050807568877*BmagInv[0]*Phi[3]*dfac_x)); 
  alphay[3] = (1.154700538379252*BdriftY[0]*dfac_y*m_*wv)/(dfac_v*q_); 
  alphay[4] = (BmagInv[0]*Bmag[1]*dfac_x*dfac_y)/(dfac_m*q_); 
  alphay[5] = -1.0*dfac_y*(1.732050807568877*BmagInv[1]*Apar[3]*dfac_x*wv-1.0*((BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2])*wv+1.732050807568877*BmagInv[1]*Phi[3]*dfac_x)); 
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
  outX[2] += 0.4330127018922193*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  outX[5] += 0.4330127018922193*(alphay[4]*f[8]+f[4]*alphay[8]+alphay[3]*f[6]+f[3]*alphay[6]+alphay[2]*f[5]+f[2]*alphay[5]+alphay[0]*f[1]+f[0]*alphay[1]); 
  outX[7] += 0.4330127018922193*(alphay[8]*f[13]+alphay[5]*f[11]+alphay[4]*f[10]+alphay[2]*f[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphay[0]*f[3]+f[0]*alphay[3]); 
  outX[9] += 0.4330127018922193*(alphay[6]*f[13]+alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4]); 
  outX[11] += 0.4330127018922193*(alphay[4]*f[13]+alphay[2]*f[11]+alphay[8]*f[10]+alphay[5]*f[7]+alphay[0]*f[6]+f[0]*alphay[6]+alphay[1]*f[3]+f[1]*alphay[3]); 
  outX[12] += 0.4330127018922193*(alphay[3]*f[13]+alphay[2]*f[12]+alphay[6]*f[10]+alphay[5]*f[9]+alphay[0]*f[8]+f[0]*alphay[8]+alphay[1]*f[4]+f[1]*alphay[4]); 
  outX[14] += 0.4330127018922193*(alphay[5]*f[15]+alphay[2]*f[14]+alphay[1]*f[13]+alphay[0]*f[10]+alphay[6]*f[8]+f[6]*alphay[8]+alphay[3]*f[4]+f[3]*alphay[4]); 
  outX[15] += 0.4330127018922193*(alphay[2]*f[15]+alphay[5]*f[14]+alphay[0]*f[13]+alphay[1]*f[10]+alphay[3]*f[8]+f[3]*alphay[8]+alphay[4]*f[6]+f[4]*alphay[6]); 
  double alphav[16]; 
  alphav[0] = -1.0*dfac_v*(1.732050807568877*((BdriftX[0]*Bmag[1]*dfac_x*wm)/q_+(BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x)*wv+(Bmag[1]*dfac_x*(1.5*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_y+0.8660254037844386*(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0]))*wm+((0.8660254037844386*(((Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[3]+(Apar[1]*BdriftY[1]+Apar[0]*BdriftY[0])*Phi[2])*dfac_y+((BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*Phi[3]+(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[1])*dfac_x)-1.0*(BmagInv[1]*(1.5*Apar[1]*Phi[3]-1.5*Phi[1]*Apar[3])+BmagInv[0]*(1.5*Apar[1]*Phi[2]-1.5*Phi[1]*Apar[2]))*dfac_x*dfac_y)*q2)/q_)/m_); 
  alphav[1] = -1.0*dfac_v*(1.732050807568877*((BdriftX[1]*Bmag[1]*dfac_x*wm)/q_+(BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x)*wv+(Bmag[1]*dfac_x*(1.5*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y+0.8660254037844386*(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1]))*wm+((1.732050807568877*((0.5*(Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])*Phi[2]-1.0*((-0.9*Apar[1]*BdriftY[1])-0.5*Apar[0]*BdriftY[0])*Phi[3])*dfac_y+0.5*((BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*Phi[3]+(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*Phi[1])*dfac_x)-1.0*(BmagInv[0]*(1.5*Apar[1]*Phi[3]-1.5*Phi[1]*Apar[3])+BmagInv[1]*(1.5*Apar[1]*Phi[2]-1.5*Phi[1]*Apar[2]))*dfac_x*dfac_y)*q2)/q_)/m_); 
  alphav[2] = -1.0*dfac_v*(1.732050807568877*BdriftX[0]*Phi[3]*dfac_x*wv+(0.8660254037844386*Bmag[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*dfac_x*wm+(0.8660254037844386*(((BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2])*Phi[3]+Phi[2]*(BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2]))*dfac_y+((Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])*Phi[3]+Phi[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2]))*dfac_x)-1.0*BmagInv[0]*(1.5*Phi[2]*Apar[3]-1.5*Apar[2]*Phi[3])*dfac_x*dfac_y)*q_)/m_); 
  alphav[3] = -1.0*((BdriftX[0]*Bmag[1]*dfac_x*wm)/q_+(BdriftY[1]*Phi[3]+BdriftY[0]*Phi[2])*dfac_y+BdriftX[0]*Phi[1]*dfac_x); 
  alphav[4] = -(1.0*Bmag[1]*dfac_v*dfac_x*((BdriftX[0]*wv)/q_+(0.8660254037844386*(BmagInv[1]*Apar[3]+BmagInv[0]*Apar[2])*dfac_y+0.5*(Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0]))/m_))/dfac_m; 
  alphav[5] = -1.0*dfac_v*(1.732050807568877*BdriftX[1]*Phi[3]*dfac_x*wv+(0.8660254037844386*Bmag[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*dfac_x*wm+(1.732050807568877*((0.5*Phi[2]*(BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2])-1.0*((-0.9*BdriftY[1]*Apar[3])-0.5*BdriftY[0]*Apar[2])*Phi[3])*dfac_y+0.5*((Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])*Phi[3]+Phi[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2]))*dfac_x)-1.0*BmagInv[1]*(1.5*Phi[2]*Apar[3]-1.5*Apar[2]*Phi[3])*dfac_x*dfac_y)*q_)/m_); 
  alphav[6] = -1.0*((BdriftX[1]*Bmag[1]*dfac_x*wm)/q_+(BdriftY[0]*Phi[3]+BdriftY[1]*Phi[2])*dfac_y+BdriftX[1]*Phi[1]*dfac_x); 
  alphav[7] = -1.0*BdriftX[0]*Phi[3]*dfac_x; 
  alphav[8] = -(1.0*Bmag[1]*dfac_v*dfac_x*((BdriftX[1]*wv)/q_+(0.8660254037844386*(BmagInv[0]*Apar[3]+BmagInv[1]*Apar[2])*dfac_y+0.5*(Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1]))/m_))/dfac_m; 
  alphav[9] = -(0.5*Bmag[1]*(BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2])*dfac_v*dfac_x)/(dfac_m*m_); 
  alphav[10] = -(0.5773502691896258*BdriftX[0]*Bmag[1]*dfac_x)/(dfac_m*q_); 
  alphav[11] = -1.0*BdriftX[1]*Phi[3]*dfac_x; 
  alphav[12] = -(0.5*Bmag[1]*(BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2])*dfac_v*dfac_x)/(dfac_m*m_); 
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
  alphaL = 0.0625*((-0.4330127018922193*alphav[13])-0.25*alphav[12]-0.4330127018922193*alphav[11]+0.4330127018922193*alphav[10]+0.25*(alphav[9]+alphav[8])+0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]-0.25*alphav[4]-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphav[13]+0.25*alphav[12]+0.4330127018922193*(alphav[11]+alphav[10])+0.25*alphav[9]-0.25*alphav[8]+0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphav[13])+0.25*alphav[12]+0.4330127018922193*(alphav[11]+alphav[10])-0.25*alphav[9]+0.25*alphav[8]-0.4330127018922193*alphav[7]+0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphav[13]-0.25*alphav[12]-0.4330127018922193*alphav[11]+0.4330127018922193*alphav[10]-0.25*(alphav[9]+alphav[8])-0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]-0.25*alphav[4]-0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphav[13]+0.25*alphav[12]-0.4330127018922193*(alphav[11]+alphav[10])-0.25*(alphav[9]+alphav[8])+0.4330127018922193*(alphav[7]+alphav[6])+0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphav[13])-0.25*alphav[12]+0.4330127018922193*alphav[11]-0.4330127018922193*alphav[10]-0.25*alphav[9]+0.25*alphav[8]+0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]-0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphav[13]-0.25*alphav[12]+0.4330127018922193*alphav[11]-0.4330127018922193*alphav[10]+0.25*alphav[9]-0.25*alphav[8]-0.4330127018922193*alphav[7]+0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphav[13])+0.25*alphav[12]-0.4330127018922193*(alphav[11]+alphav[10])+0.25*(alphav[9]+alphav[8])-0.4330127018922193*(alphav[7]+alphav[6])+0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.4330127018922193*alphav[13]-0.25*alphav[12]+0.4330127018922193*alphav[11]-0.4330127018922193*alphav[10]+0.25*(alphav[9]+alphav[8])-0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]-0.25*alphav[4]+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphav[13])+0.25*alphav[12]-0.4330127018922193*(alphav[11]+alphav[10])+0.25*alphav[9]-0.25*alphav[8]-0.4330127018922193*alphav[7]+0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphav[13]+0.25*alphav[12]-0.4330127018922193*(alphav[11]+alphav[10])-0.25*alphav[9]+0.25*alphav[8]+0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphav[13])-0.25*alphav[12]+0.4330127018922193*alphav[11]-0.4330127018922193*alphav[10]-0.25*(alphav[9]+alphav[8])+0.4330127018922193*(alphav[7]+alphav[6])+0.25*alphav[5]-0.25*alphav[4]+0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphav[13])+0.25*alphav[12]+0.4330127018922193*(alphav[11]+alphav[10])-0.25*(alphav[9]+alphav[8])-0.4330127018922193*(alphav[7]+alphav[6])+0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphav[13]-0.25*alphav[12]-0.4330127018922193*alphav[11]+0.4330127018922193*alphav[10]-0.25*alphav[9]+0.25*alphav[8]-0.4330127018922193*alphav[7]+0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]+0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphav[13])-0.25*alphav[12]-0.4330127018922193*alphav[11]+0.4330127018922193*alphav[10]+0.25*alphav[9]-0.25*alphav[8]+0.4330127018922193*alphav[7]-0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphav[13]+0.25*alphav[12]+0.4330127018922193*(alphav[11]+alphav[10])+0.25*(alphav[9]+alphav[8])+0.4330127018922193*(alphav[7]+alphav[6])+0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  cflRateByDir[3] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#endif 
  cflRateByDir[0] += cflRateByDir[3];
  outV[3] += 0.4330127018922193*(alphav[13]*f[13]+alphav[12]*f[12]+alphav[11]*f[11]+alphav[10]*f[10]+alphav[9]*f[9]+alphav[8]*f[8]+alphav[7]*f[7]+alphav[6]*f[6]+alphav[5]*f[5]+alphav[4]*f[4]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  outV[6] += 0.4330127018922193*(alphav[10]*f[13]+f[10]*alphav[13]+alphav[9]*f[12]+f[9]*alphav[12]+alphav[7]*f[11]+f[7]*alphav[11]+alphav[4]*f[8]+f[4]*alphav[8]+alphav[3]*f[6]+f[3]*alphav[6]+alphav[2]*f[5]+f[2]*alphav[5]+alphav[0]*f[1]+f[0]*alphav[1]); 
  outV[7] += 0.4330127018922193*(alphav[13]*f[15]+alphav[10]*f[14]+alphav[8]*f[12]+f[8]*alphav[12]+alphav[6]*f[11]+f[6]*alphav[11]+alphav[4]*f[9]+f[4]*alphav[9]+alphav[3]*f[7]+f[3]*alphav[7]+alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[2]+f[0]*alphav[2]); 
  outV[10] += 0.4330127018922193*(alphav[11]*f[15]+alphav[7]*f[14]+alphav[6]*f[13]+f[6]*alphav[13]+alphav[5]*f[12]+f[5]*alphav[12]+alphav[3]*f[10]+f[3]*alphav[10]+alphav[2]*f[9]+f[2]*alphav[9]+alphav[1]*f[8]+f[1]*alphav[8]+alphav[0]*f[4]+f[0]*alphav[4]); 
  outV[11] += 0.4330127018922193*(alphav[10]*f[15]+alphav[13]*f[14]+alphav[4]*f[12]+f[4]*alphav[12]+alphav[3]*f[11]+f[3]*alphav[11]+alphav[8]*f[9]+f[8]*alphav[9]+alphav[6]*f[7]+f[6]*alphav[7]+alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[2]+f[1]*alphav[2]); 
  outV[13] += 0.4330127018922193*(alphav[7]*f[15]+alphav[11]*f[14]+alphav[3]*f[13]+f[3]*alphav[13]+alphav[2]*f[12]+f[2]*alphav[12]+alphav[6]*f[10]+f[6]*alphav[10]+alphav[5]*f[9]+f[5]*alphav[9]+alphav[0]*f[8]+f[0]*alphav[8]+alphav[1]*f[4]+f[1]*alphav[4]); 
  outV[14] += 0.4330127018922193*(alphav[6]*f[15]+alphav[3]*f[14]+alphav[11]*f[13]+f[11]*alphav[13]+alphav[1]*f[12]+f[1]*alphav[12]+alphav[7]*f[10]+f[7]*alphav[10]+alphav[0]*f[9]+f[0]*alphav[9]+alphav[5]*f[8]+f[5]*alphav[8]+alphav[2]*f[4]+f[2]*alphav[4]); 
  outV[15] += 0.4330127018922193*(alphav[3]*f[15]+alphav[6]*f[14]+alphav[7]*f[13]+f[7]*alphav[13]+alphav[0]*f[12]+f[0]*alphav[12]+alphav[10]*f[11]+f[10]*alphav[11]+alphav[1]*f[9]+f[1]*alphav[9]+alphav[2]*f[8]+f[2]*alphav[8]+alphav[4]*f[5]+f[4]*alphav[5]); 
  return cflRate; 
} 
