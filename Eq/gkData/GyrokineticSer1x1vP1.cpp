#include <GyrokineticModDecl.h> 
double GyrokineticVol1x1vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphaM = 0.0; 
  double alphaP = 0.0; 
  double alphax[4]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*dfac_x*wv; 
  alphaM = 0.25*alphax[0]; 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.25*alphax[0]; 
  if(alphaP>0) cflFreq += alphaP; 
  double alphav[4]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  alphaM = 0.25*alphav[0]; 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.25*alphav[0]; 
  if(alphaP>0) cflFreq += alphaP; 
  out[1] += 0.8660254037844386*alphax[0]*f[0]; 
  out[2] += 0.8660254037844386*alphav[0]*f[0]; 
  out[3] += 0.8660254037844386*(alphax[0]*f[2]+alphav[0]*f[1]); 
  return cflFreq; 
} 
double GyrokineticVol1x1vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphaM = 0.0; 
  double alphaP = 0.0; 
  double alphax[4]; 
  alphax[0] = 1.414213562373095*Gradpar[0]*dfac_x*wv; 
  alphax[1] = 1.414213562373095*Gradpar[1]*dfac_x*wv; 
  alphaM = -0.25*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.25*(1.732050807568877*alphax[1]+alphax[0]); 
  if(alphaP>0) cflFreq += alphaP; 
  double alphav[4]; 
  alphav[0] = -(1.732050807568877*Gradpar[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  alphav[1] = -(1.732050807568877*Gradpar[1]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
  alphaM = 0.25*alphav[0]; 
  if(alphaM<0) cflFreq += -alphaM; 
  alphaP = 0.25*alphav[0]; 
  if(alphaP>0) cflFreq += alphaP; 
  out[1] += 0.8660254037844386*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphav[1]*f[1]+alphav[0]*f[0]); 
  out[3] += 0.8660254037844386*(alphax[1]*f[3]+alphax[0]*f[2]+alphav[0]*f[1]+f[0]*alphav[1]); 
  return cflFreq; 
} 
