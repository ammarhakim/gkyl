#include <GyrokineticModDecl.h> 
double GyrokineticVol2x0vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, double *cflRateCtrl, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_y = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wy = w[1]; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaCtrl;
  double alphax[4]; 
  alphax[0] = -0.8660254037844386*BmagInv[0]*Phi[2]*dfac_x*dfac_y; 
  alphax[1] = -0.8660254037844386*BmagInv[0]*Phi[3]*dfac_x*dfac_y; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.25*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.5*alphax[0]-0.8660254037844386*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.5*alphax[0]-0.8660254037844386*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  alphaCtrl = 0.5*alphax[0]-0.8660254037844386*alphax[1]; 
  cflRateCtrl[0] = -0.5*(alphaCtrl-std::abs(alphaCtrl)); 
  alphaCtrl = 0.8660254037844386*alphax[1]+0.5*alphax[0]; 
  cflRateCtrl[1] = 0.5*(alphaCtrl+std::abs(alphaCtrl)); 
  alphaCtrl = 0.5*alphax[0]-0.8660254037844386*alphax[1]; 
  cflRateCtrl[2] = -0.5*(alphaCtrl-std::abs(alphaCtrl)); 
  alphaCtrl = 0.8660254037844386*alphax[1]+0.5*alphax[0]; 
  cflRateCtrl[3] = 0.5*(alphaCtrl+std::abs(alphaCtrl)); 
  double alphay[4]; 
  alphay[0] = 0.8660254037844386*BmagInv[0]*Phi[1]*dfac_x*dfac_y; 
  alphay[2] = 0.8660254037844386*BmagInv[0]*Phi[3]*dfac_x*dfac_y; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.25*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.5*alphay[0]-0.8660254037844386*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.5*alphay[0]-0.8660254037844386*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(0.8660254037844386*alphay[2]+0.5*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.8660254037844386*alphay[2]+0.5*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  alphaCtrl = 0.5*alphay[0]-0.8660254037844386*alphay[2]; 
  cflRateCtrl[0] += -0.5*(alphaCtrl-std::abs(alphaCtrl)); 
  alphaCtrl = 0.5*alphay[0]-0.8660254037844386*alphay[2]; 
  cflRateCtrl[1] += -0.5*(alphaCtrl-std::abs(alphaCtrl)); 
  alphaCtrl = 0.8660254037844386*alphay[2]+0.5*alphay[0]; 
  cflRateCtrl[2] += 0.5*(alphaCtrl+std::abs(alphaCtrl)); 
  alphaCtrl = 0.8660254037844386*alphay[2]+0.5*alphay[0]; 
  cflRateCtrl[3] += 0.5*(alphaCtrl+std::abs(alphaCtrl)); 
  out[1] += 0.8660254037844386*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphay[2]*f[2]+alphay[0]*f[0]); 
  out[3] += 0.8660254037844386*((alphay[2]+alphax[1])*f[3]+alphax[0]*f[2]+alphay[0]*f[1]); 
  return cflFreq; 
} 
