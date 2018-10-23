#include <GyrokineticModDecl.h> 
double EmGyrokineticVol1x2vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
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
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = 2.0*Gradpar[0]*dfac_x*wv; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0441941738241592*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0441941738241592*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0441941738241592*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0441941738241592*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0441941738241592*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0441941738241592*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0441941738241592*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0441941738241592*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[8]; 
  alphav[0] = -(2.449489742783178*Gradpar[0]*Phi[1]*dfac_v*dfac_x*q_)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0441941738241592*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0441941738241592*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0441941738241592*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0441941738241592*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0441941738241592*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0441941738241592*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0441941738241592*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0441941738241592*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.6123724356957944*alphax[0]*f[0]; 
  out[2] += 0.6123724356957944*alphav[0]*f[0]; 
  out[4] += 0.6123724356957944*(alphax[0]*f[2]+alphav[0]*f[1]); 
  out[5] += 0.6123724356957944*alphax[0]*f[3]; 
  out[6] += 0.6123724356957944*alphav[0]*f[3]; 
  out[7] += 0.6123724356957944*(alphax[0]*f[6]+alphav[0]*f[5]); 
  return cflFreq; 
} 
double EmGyrokineticStep2Vol1x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[1]; 
  double dfac_v = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[4] += -(1.224744871391589*(dApardt[0]*f[1]+f[0]*dApardt[1])*dfac_v*q_)/m_; 
  out[6] += -(1.224744871391589*(dApardt[1]*f[5]+dApardt[0]*f[3])*dfac_v*q_)/m_; 
  out[7] += -(1.224744871391589*(dApardt[0]*f[5]+dApardt[1]*f[3])*dfac_v*q_)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -(0.3535533905932737*dApardt[0]*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = -(0.3535533905932737*dApardt[0]*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = -(0.125*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.08838834764831843*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.125*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.08838834764831843*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = -(0.125*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.08838834764831843*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.125*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.08838834764831843*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
#endif 
  return cflFreq; 
} 
double EmGyrokineticVol1x2vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *BdriftX, const double *BdriftY, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
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
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = 2.0*Gradpar[0]*dfac_x*wv; 
  alphax[1] = 2.0*Gradpar[1]*dfac_x*wv; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(2.449489742783178*alphax[1]-1.414213562373095*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.3535533905932737*alphax[0]-0.6123724356957944*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*alphax[0]-0.6123724356957944*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*alphax[0]-0.6123724356957944*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*alphax[0]-0.6123724356957944*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[8]; 
  alphav[0] = -(2.449489742783178*Gradpar[0]*dfac_v*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 
  alphav[1] = -(2.449489742783178*Gradpar[1]*dfac_v*dfac_x*(Bmag[1]*wm+Phi[1]*q_))/m_; 
  alphav[3] = -(1.414213562373095*Gradpar[0]*Bmag[1]*dfac_v*dfac_x)/(dfac_m*m_); 
  alphav[5] = -(1.414213562373095*Bmag[1]*Gradpar[1]*dfac_v*dfac_x)/(dfac_m*m_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.3535533905932737*alphav[5]-0.3535533905932737*(alphav[3]+alphav[1])+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*(alphav[1]+alphav[0])-0.3535533905932737*(alphav[5]+alphav[3])); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.3535533905932737*alphav[5])+0.3535533905932737*alphav[3]-0.3535533905932737*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.04419417382415921*(alphav[5]+alphav[3]+alphav[1]+alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(0.3535533905932737*alphav[5]-0.3535533905932737*(alphav[3]+alphav[1])+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.3535533905932737*(alphav[1]+alphav[0])-0.3535533905932737*(alphav[5]+alphav[3])); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*((-0.3535533905932737*alphav[5])+0.3535533905932737*alphav[3]-0.3535533905932737*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.04419417382415921*(alphav[5]+alphav[3]+alphav[1]+alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.6123724356957944*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphav[5]*f[5]+alphav[3]*f[3]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphav[3]*f[5]+f[3]*alphav[5]+alphax[1]*f[4]+alphax[0]*f[2]+alphav[0]*f[1]+f[0]*alphav[1]); 
  out[5] += 0.6123724356957944*(alphax[1]*f[5]+alphax[0]*f[3]); 
  out[6] += 0.6123724356957944*(alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[3]+f[0]*alphav[3]); 
  out[7] += 0.6123724356957944*(alphax[1]*f[7]+alphax[0]*f[6]+alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[3]+f[1]*alphav[3]); 
  return cflFreq; 
} 
