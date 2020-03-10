#include <GyrokineticModDecl.h> 
double GyrokineticVolPositivity1x1vSer_P1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double *f, double *out, double *positivityWeightByDir) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflRate = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaCtrl;
  double alphax[4]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*dfac_x*wv; 
  positivityWeightByDir[1] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*alphax[0]; 
  positivityWeightByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*alphax[0]; 
  positivityWeightByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*alphax[0]; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); //std::abs(alphaL); 
  alphaL = 0.125*alphax[0]; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); //std::abs(alphaL); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*alphax[0]; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); //std::abs(alphaL); 
  alphaR = 0.125*alphax[0]; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); //std::abs(alphaL); 
#endif 

  out[1] += 0.8660254037844386*alphax[0]*f[0]; 
  out[3] += 0.8660254037844386*alphax[0]*f[2]; 
  double alphav[4]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  positivityWeightByDir[2] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*alphav[0]; 
  positivityWeightByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*alphav[0]; 
  positivityWeightByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*alphav[0]; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); //std::abs(alphaL); 
  alphaL = 0.125*alphav[0]; 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); //std::abs(alphaL); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*alphav[0]; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); //std::abs(alphaL); 
  alphaR = 0.125*alphav[0]; 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); //std::abs(alphaL); 
#endif 

  out[2] += 0.8660254037844386*alphav[0]*f[0]; 
  out[3] += 0.8660254037844386*alphav[0]*f[1]; 
  return cflRate; 
} 
double GyrokineticVolPositivity1x1vSer_P1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, 
                        const double *Bmag, const double *BmagInv, const double *Gradpar, 
                        const double *BdriftX, const double *BdriftY, const double *Phi, 
                        const double *f, double *out, double *positivityWeightByDir) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflRate = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaCtrl;
  double alphax[4]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*dfac_x*wv; 
  alphax[1] = 1.414213562373095*Gradpar[1]*dfac_x*wv; 
  positivityWeightByDir[1] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.25*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  positivityWeightByDir[1] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(1.732050807568877*alphax[1]+alphax[0]); 
  positivityWeightByDir[1] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.5*alphax[0]-0.8660254037844386*alphax[1]); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); //std::abs(alphaL); 
  alphaL = 0.25*(0.5*alphax[0]-0.8660254037844386*alphax[1]); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[1] += -0.5*(alphaL-std::abs(alphaL)); //std::abs(alphaL); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); //std::abs(alphaL); 
  alphaR = 0.25*(0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[1] += 0.5*(alphaR+std::abs(alphaR)); //std::abs(alphaL); 
#endif 

  out[1] += 0.8660254037844386*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphax[1]*f[3]+alphax[0]*f[2]); 
  double alphav[4]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  alphav[1] = -(1.732050807568877*Gradpar[1]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  positivityWeightByDir[2] = 0.; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*alphav[0]; 
  positivityWeightByDir[2] += std::abs(alphaL); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*alphav[0]; 
  positivityWeightByDir[2] += std::abs(alphaR); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.5*alphav[0]-0.5*alphav[1]); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); //std::abs(alphaL); 
  alphaL = 0.125*(alphav[1]+alphav[0]); 
  cflRate += -0.5*(alphaL-std::abs(alphaL)); 
  positivityWeightByDir[2] += -0.5*(alphaL-std::abs(alphaL)); //std::abs(alphaL); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(0.5*alphav[0]-0.5*alphav[1]); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); //std::abs(alphaL); 
  alphaR = 0.125*(alphav[1]+alphav[0]); 
  cflRate += 0.5*(alphaR+std::abs(alphaR)); 
  positivityWeightByDir[2] += 0.5*(alphaR+std::abs(alphaR)); //std::abs(alphaL); 
#endif 

  out[2] += 0.8660254037844386*(alphav[1]*f[1]+alphav[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphav[0]*f[1]+f[0]*alphav[1]); 
  return cflRate; 
} 
