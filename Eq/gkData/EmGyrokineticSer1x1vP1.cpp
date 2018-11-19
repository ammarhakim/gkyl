#include <GyrokineticModDecl.h> 
double EmGyrokineticVol1x1vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[4]; 
  alphax[0] = (1.224744871391589*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2+1.414213562373095*Gradpar[0]*dfac_x)*wv; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*alphax[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*alphax[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[4]; 
  alphav[0] = (Phi[1]*dfac_v*((-1.5*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2)-1.732050807568877*Gradpar[0]*dfac_x)*q_)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.8660254037844386*alphax[0]*f[0]; 
  out[2] += 0.8660254037844386*alphav[0]*f[0]; 
  out[3] += 0.8660254037844386*(alphax[0]*f[2]+alphav[0]*f[1]); 
  return cflFreq; 
} 
double EmGyrokineticStep2Vol1x1vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  double dvInv = 1.0/dxv[1]; 
  double dfac_v = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[1]*f[1]+dApardt[0]*f[0])*dfac_v*q_)/m_; 
  out[3] += -(1.224744871391589*(dApardt[0]*f[1]+f[0]*dApardt[1])*dfac_v*q_)/m_; 
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
  alphaL = -(0.25*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.1767766952966369*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = -(0.25*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.1767766952966369*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
#endif 
  return cflFreq; 
} 
double EmGyrokineticVol1x1vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[4]; 
  alphax[0] = -1.0*((1.224744871391589*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_*wv2)/q_-1.0*(1.732050807568877*((Apar[1]*(0.7071067811865475*BmagInv[1]-0.5*BmagInv[0]*Bmag[1])-0.5*Apar[0]*Bmag[1]*BmagInv[1])*geoY[1]+geoY[0]*(BmagInv[0]*(0.7071067811865475*Apar[1]-0.5*Apar[0]*Bmag[1])-0.5*Apar[1]*Bmag[1]*BmagInv[1]))*dfac_x2+1.414213562373095*Gradpar[0]*dfac_x)*wv); 
  alphax[1] = -1.0*((1.224744871391589*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_*wv2)/q_-1.0*(1.732050807568877*((BmagInv[0]*(0.7071067811865475*Apar[1]-0.5*Apar[0]*Bmag[1])-0.9*Apar[1]*Bmag[1]*BmagInv[1])*geoY[1]+geoY[0]*(Apar[1]*(0.7071067811865475*BmagInv[1]-0.5*BmagInv[0]*Bmag[1])-0.5*Apar[0]*Bmag[1]*BmagInv[1]))*dfac_x2+1.414213562373095*Gradpar[1]*dfac_x)*wv); 
  alphax[2] = -(0.7071067811865475*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_*wv)/(dfac_v*q_); 
  alphax[3] = -(0.7071067811865475*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_*wv)/(dfac_v*q_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.25*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(1.732050807568877*alphax[1]+alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.8660254037844386*alphax[3]-0.5*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.25*((-0.8660254037844386*alphax[3])+0.5*alphax[2]-0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*((-0.8660254037844386*alphax[3])-0.5*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.25*(0.8660254037844386*alphax[3]+0.5*alphax[2]+0.8660254037844386*alphax[1]+0.5*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[4]; 
  alphav[0] = Phi[1]*dfac_v*(1.5*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*wv+((((1.060660171779821*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(1.060660171779821*BmagInv[0]*Bmag[1]-1.5*BmagInv[1]))*geoY[1]+geoY[0]*(1.060660171779821*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(1.060660171779821*Apar[0]*Bmag[1]-1.5*Apar[1])))*dfac_x2-1.732050807568877*Gradpar[0]*dfac_x)*q_)/m_); 
  alphav[1] = Phi[1]*dfac_v*(1.5*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*wv+((((1.909188309203678*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(1.060660171779821*Apar[0]*Bmag[1]-1.5*Apar[1]))*geoY[1]+geoY[0]*(1.060660171779821*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(1.060660171779821*BmagInv[0]*Bmag[1]-1.5*BmagInv[1])))*dfac_x2-1.732050807568877*Gradpar[1]*dfac_x)*q_)/m_); 
  alphav[2] = 0.8660254037844386*Bmag[1]*Phi[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2; 
  alphav[3] = 0.8660254037844386*Bmag[1]*Phi[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.25*(1.732050807568877*alphav[2]-1.0*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(1.732050807568877*alphav[2]+alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.8660254037844386*alphav[3]-0.8660254037844386*alphav[2]-0.5*alphav[1]+0.5*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.25*(0.5*(alphav[1]+alphav[0])-0.8660254037844386*(alphav[3]+alphav[2])); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*((-0.8660254037844386*alphav[3])+0.8660254037844386*alphav[2]-0.5*alphav[1]+0.5*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.25*(0.8660254037844386*(alphav[3]+alphav[2])+0.5*(alphav[1]+alphav[0])); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.8660254037844386*(alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[3] += 0.8660254037844386*((alphav[2]+alphax[1])*f[3]+f[1]*alphax[3]+f[2]*(alphav[3]+alphax[0])+f[0]*alphax[2]+alphav[0]*f[1]+f[0]*alphav[1]); 
  return cflFreq; 
} 
