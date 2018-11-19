#include <GyrokineticModDecl.h> 
double EmGyrokineticVol2x2vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*geoZ[0]*dfac_x*dfac_y*(0.5*BmagInv[0]*Apar[2]*wv-1.0*Phi[2]); 
  alphax[1] = 1.732050807568877*geoZ[0]*dfac_x*dfac_y*(0.5*BmagInv[0]*Apar[3]*wv-1.0*Phi[3]); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphay[16]; 
  alphay[0] = 1.732050807568877*geoZ[0]*dfac_x*dfac_y*(Phi[1]-0.5*BmagInv[0]*Apar[1]*wv); 
  alphay[2] = 1.732050807568877*geoZ[0]*dfac_x*dfac_y*(Phi[3]-0.5*BmagInv[0]*Apar[3]*wv); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[16]; 
  alphav[0] = (BmagInv[0]*geoZ[0]*(0.75*Apar[1]*Phi[2]-0.75*Phi[1]*Apar[2])*dfac_v*dfac_x*dfac_y*q_)/m_; 
  alphav[1] = (BmagInv[0]*geoZ[0]*(0.75*Apar[1]*Phi[3]-0.75*Phi[1]*Apar[3])*dfac_v*dfac_x*dfac_y*q_)/m_; 
  alphav[2] = (BmagInv[0]*geoZ[0]*(0.75*Phi[2]*Apar[3]-0.75*Apar[2]*Phi[3])*dfac_v*dfac_x*dfac_y*q_)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.125*alphav[0]; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*alphav[0]; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*alphav[0]-0.25*(alphav[2]+alphav[1])); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*(alphav[1]+alphav[0])-0.25*alphav[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.015625*(alphav[2]+alphav[1]+alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphav[0]-0.25*(alphav[2]+alphav[1])); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*(alphav[1]+alphav[0])-0.25*alphav[2]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.015625*(alphav[2]+alphav[1]+alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.25*alphav[0]-0.25*(alphav[2]+alphav[1])); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*(alphav[1]+alphav[0])-0.25*alphav[2]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.015625*(alphav[2]+alphav[1]+alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*alphav[0]-0.25*(alphav[2]+alphav[1])); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*(alphav[1]+alphav[0])-0.25*alphav[2]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.015625*(alphav[2]+alphav[1]+alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.4330127018922193*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[2]*f[2]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[5] += 0.4330127018922193*((alphay[2]+alphax[1])*f[5]+alphax[0]*f[2]+alphay[0]*f[1]); 
  out[6] += 0.4330127018922193*(alphax[1]*f[6]+alphav[2]*f[5]+alphax[0]*f[3]+alphav[0]*f[1]+f[0]*alphav[1]); 
  out[7] += 0.4330127018922193*(alphay[2]*f[7]+alphav[1]*f[5]+alphay[0]*f[3]+alphav[0]*f[2]+f[0]*alphav[2]); 
  out[8] += 0.4330127018922193*(alphax[1]*f[8]+alphax[0]*f[4]); 
  out[9] += 0.4330127018922193*(alphay[2]*f[9]+alphay[0]*f[4]); 
  out[10] += 0.4330127018922193*(alphav[2]*f[9]+alphav[1]*f[8]+alphav[0]*f[4]); 
  out[11] += 0.4330127018922193*((alphay[2]+alphax[1])*f[11]+alphax[0]*f[7]+alphay[0]*f[6]+alphav[0]*f[5]+alphav[1]*f[2]+f[1]*alphav[2]); 
  out[12] += 0.4330127018922193*((alphay[2]+alphax[1])*f[12]+alphax[0]*f[9]+alphay[0]*f[8]); 
  out[13] += 0.4330127018922193*(alphax[1]*f[13]+alphav[2]*f[12]+alphax[0]*f[10]+alphav[0]*f[8]+alphav[1]*f[4]); 
  out[14] += 0.4330127018922193*(alphay[2]*f[14]+alphav[1]*f[12]+alphay[0]*f[10]+alphav[0]*f[9]+alphav[2]*f[4]); 
  out[15] += 0.4330127018922193*((alphay[2]+alphax[1])*f[15]+alphax[0]*f[14]+alphay[0]*f[13]+alphav[0]*f[12]+alphav[1]*f[9]+alphav[2]*f[8]); 
  return cflFreq; 
} 
double EmGyrokineticStep2Vol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -(0.25*dApardt[0]*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = -(0.25*dApardt[0]*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*dfac_v*q_)/m_; 
  if(alphaR>0) cflFreq += alphaR; 
#endif 
  return cflFreq; 
} 
double EmGyrokineticVol2x2vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double *f, double *out) 
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
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*dfac_x*dfac_y*(0.5*((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2])*wv-1.0*(geoZ[1]*Phi[3]+geoZ[0]*Phi[2])); 
  alphax[1] = 1.732050807568877*dfac_x*dfac_y*((0.5*(BmagInv[0]*geoZ[0]*Apar[3]+(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2])+0.9*BmagInv[1]*geoZ[1]*Apar[3])*wv-1.0*(geoZ[0]*Phi[3]+geoZ[1]*Phi[2])); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphay[16]; 
  alphay[0] = 1.732050807568877*dfac_x*dfac_y*((0.5*Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*m_*wv2)/q_+geoZ[0]*((0.25*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(0.25*Apar[0]*Bmag[1]-0.5*Apar[1]))*wv+(Bmag[1]*wm)/q_+Phi[1])+(0.25*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(0.25*BmagInv[0]*Bmag[1]-0.5*BmagInv[1]))*geoZ[1]*wv); 
  alphay[1] = 1.732050807568877*dfac_x*dfac_y*((0.5*Bmag[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*m_*wv2)/q_+((0.45*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(0.25*Apar[0]*Bmag[1]-0.5*Apar[1]))*geoZ[1]+geoZ[0]*(0.25*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(0.25*BmagInv[0]*Bmag[1]-0.5*BmagInv[1])))*wv+geoZ[1]*((Bmag[1]*wm)/q_+Phi[1])); 
  alphay[2] = 1.732050807568877*dfac_x*dfac_y*(geoZ[0]*(0.25*BmagInv[0]*Bmag[1]*Apar[2]*wv+Phi[3])+((0.25*BmagInv[0]*Bmag[1]*geoZ[1]+geoZ[0]*(0.25*Bmag[1]*BmagInv[1]-0.5*BmagInv[0]))*Apar[3]+BmagInv[1]*geoZ[1]*(0.25*Bmag[1]*Apar[2]-0.5*Apar[3]))*wv); 
  alphay[3] = (0.5*Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*dfac_x*dfac_y*m_*wv)/(dfac_v*q_); 
  alphay[4] = (geoZ[0]*Bmag[1]*dfac_x*dfac_y)/(dfac_m*q_); 
  alphay[5] = 1.732050807568877*dfac_x*dfac_y*((0.25*Bmag[1]*(BmagInv[0]*geoZ[0]*Apar[3]+(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2])+(0.45*Bmag[1]*BmagInv[1]*geoZ[1]-0.5*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1]))*Apar[3])*wv+geoZ[1]*Phi[3]); 
  alphay[6] = (0.5*Bmag[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*dfac_x*dfac_y*m_*wv)/(dfac_v*q_); 
  alphay[8] = (Bmag[1]*geoZ[1]*dfac_x*dfac_y)/(dfac_m*q_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[16]; 
  alphav[0] = dfac_v*dfac_x*dfac_y*(((((BmagInv[0]*(0.75*Apar[1]-0.375*Apar[0]*Bmag[1])-0.675*Apar[1]*Bmag[1]*BmagInv[1])*geoZ[1]+geoZ[0]*(Apar[1]*(0.75*BmagInv[1]-0.375*BmagInv[0]*Bmag[1])-0.375*Apar[0]*Bmag[1]*BmagInv[1]))*Phi[3]-0.75*Phi[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3]+((Apar[1]*(0.75*BmagInv[1]-0.375*BmagInv[0]*Bmag[1])-0.375*Apar[0]*Bmag[1]*BmagInv[1])*geoZ[1]+geoZ[0]*(BmagInv[0]*(0.75*Apar[1]-0.375*Apar[0]*Bmag[1])-0.375*Apar[1]*Bmag[1]*BmagInv[1]))*Phi[2]-0.75*Phi[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2])*q_-0.75*Bmag[1]*((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2])*wm)/m_-0.75*Bmag[1]*((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Phi[2])*wv); 
  alphav[1] = dfac_v*dfac_x*dfac_y*(Bmag[1]*((-0.75*(BmagInv[0]*geoZ[0]*Phi[3]+(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[2]))-1.35*BmagInv[1]*geoZ[1]*Phi[3])*wv+(Bmag[1]*((-0.75*(BmagInv[0]*geoZ[0]*Apar[3]+(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2]))-1.35*BmagInv[1]*geoZ[1]*Apar[3])*wm+(((Apar[1]*(1.35*BmagInv[1]-0.675*BmagInv[0]*Bmag[1])-0.675*Apar[0]*Bmag[1]*BmagInv[1])*geoZ[1]+geoZ[0]*(BmagInv[0]*(0.75*Apar[1]-0.375*Apar[0]*Bmag[1])-0.675*Apar[1]*Bmag[1]*BmagInv[1]))*Phi[3]+Phi[1]*((-1.35*BmagInv[1]*geoZ[1])-0.75*BmagInv[0]*geoZ[0])*Apar[3]+((BmagInv[0]*(0.75*Apar[1]-0.375*Apar[0]*Bmag[1])-0.675*Apar[1]*Bmag[1]*BmagInv[1])*geoZ[1]+geoZ[0]*(Apar[1]*(0.75*BmagInv[1]-0.375*BmagInv[0]*Bmag[1])-0.375*Apar[0]*Bmag[1]*BmagInv[1]))*Phi[2]-0.75*Phi[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2])*q_)/m_); 
  alphav[2] = (((Bmag[1]*((-0.675*BmagInv[1]*geoZ[1])-0.375*BmagInv[0]*geoZ[0])*Apar[3]+(((-0.75*BmagInv[1])-0.375*BmagInv[0]*Bmag[1])*geoZ[1]+geoZ[0]*((-0.375*Bmag[1]*BmagInv[1])-0.75*BmagInv[0]))*Apar[2])*Phi[3]+Phi[2]*(((0.75*BmagInv[1]-0.375*BmagInv[0]*Bmag[1])*geoZ[1]+geoZ[0]*(0.75*BmagInv[0]-0.375*Bmag[1]*BmagInv[1]))*Apar[3]-0.375*Bmag[1]*(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2]))*dfac_v*dfac_x*dfac_y*q_)/m_; 
  alphav[3] = -0.4330127018922193*Bmag[1]*((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Phi[2])*dfac_x*dfac_y; 
  alphav[4] = -(0.4330127018922193*Bmag[1]*((BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3]+(BmagInv[1]*geoZ[1]+BmagInv[0]*geoZ[0])*Apar[2])*dfac_v*dfac_x*dfac_y)/(dfac_m*m_); 
  alphav[5] = ((((((-0.675*Bmag[1]*BmagInv[1])-0.75*BmagInv[0])*geoZ[1]+geoZ[0]*((-0.75*BmagInv[1])-0.375*BmagInv[0]*Bmag[1]))*Apar[2]-0.675*Bmag[1]*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[3])*Phi[3]+Phi[2]*((0.75*(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])-0.675*Bmag[1]*BmagInv[1]*geoZ[1])*Apar[3]-0.375*Bmag[1]*(BmagInv[0]*geoZ[0]*Apar[3]+(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2])))*dfac_v*dfac_x*dfac_y*q_)/m_; 
  alphav[6] = 1.732050807568877*Bmag[1]*((-0.25*(BmagInv[0]*geoZ[0]*Phi[3]+(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Phi[2]))-0.45*BmagInv[1]*geoZ[1]*Phi[3])*dfac_x*dfac_y; 
  alphav[8] = (1.732050807568877*Bmag[1]*((-0.25*(BmagInv[0]*geoZ[0]*Apar[3]+(BmagInv[0]*geoZ[1]+geoZ[0]*BmagInv[1])*Apar[2]))-0.45*BmagInv[1]*geoZ[1]*Apar[3])*dfac_v*dfac_x*dfac_y)/(dfac_m*m_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(1.732050807568877*alphav[3]-1.0*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(1.732050807568877*alphav[3]+alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0625*(0.25*alphav[8]+0.4330127018922193*alphav[6]+0.25*alphav[5]-0.25*alphav[4]-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*((-0.25*alphav[8])-0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphav[8]+0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*((-0.25*alphav[8])-0.4330127018922193*alphav[6]+0.25*alphav[5]-0.25*alphav[4]-0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*((-0.25*alphav[8])+0.4330127018922193*alphav[6]+0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphav[8]-0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]-0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*((-0.25*alphav[8])+0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]-0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.0625*(0.25*alphav[8]-0.4330127018922193*alphav[6]+0.25*(alphav[5]+alphav[4])-0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0625*(0.25*alphav[8]-0.4330127018922193*alphav[6]+0.25*alphav[5]-0.25*alphav[4]+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*((-0.25*alphav[8])+0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*alphav[8]-0.4330127018922193*alphav[6]-0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*((-0.25*alphav[8])+0.4330127018922193*alphav[6]+0.25*alphav[5]-0.25*alphav[4]+0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*((-0.25*alphav[8])-0.4330127018922193*alphav[6]+0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]-0.25*(alphav[2]+alphav[1])+0.25*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*alphav[8]+0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]+0.4330127018922193*alphav[3]-0.25*alphav[2]+0.25*(alphav[1]+alphav[0])); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*((-0.25*alphav[8])-0.4330127018922193*alphav[6]-0.25*alphav[5]+0.25*alphav[4]+0.4330127018922193*alphav[3]+0.25*alphav[2]-0.25*alphav[1]+0.25*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.0625*(0.25*alphav[8]+0.4330127018922193*alphav[6]+0.25*(alphav[5]+alphav[4])+0.4330127018922193*alphav[3]+0.25*(alphav[2]+alphav[1]+alphav[0])); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.4330127018922193*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphav[8]*f[8]+alphav[6]*f[6]+alphav[5]*f[5]+alphav[4]*f[4]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[5] += 0.4330127018922193*(alphay[4]*f[8]+f[4]*alphay[8]+alphay[3]*f[6]+f[3]*alphay[6]+(alphay[2]+alphax[1])*f[5]+f[2]*(alphay[5]+alphax[0])+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[6] += 0.4330127018922193*(alphav[4]*f[8]+f[4]*alphav[8]+(alphav[3]+alphax[1])*f[6]+f[3]*alphav[6]+alphav[2]*f[5]+f[2]*alphav[5]+alphax[0]*f[3]+alphav[0]*f[1]+f[0]*alphav[1]); 
  out[7] += 0.4330127018922193*(alphay[8]*f[13]+alphav[8]*f[12]+(alphav[6]+alphay[5])*f[11]+alphay[4]*f[10]+alphav[4]*f[9]+(alphav[3]+alphay[2])*f[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphav[1]*f[5]+f[1]*alphav[5]+alphay[0]*f[3]+f[0]*alphay[3]+alphav[0]*f[2]+f[0]*alphav[2]); 
  out[8] += 0.4330127018922193*(alphax[1]*f[8]+alphax[0]*f[4]); 
  out[9] += 0.4330127018922193*(alphay[6]*f[13]+alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4]); 
  out[10] += 0.4330127018922193*(alphav[6]*f[13]+alphav[5]*f[12]+alphav[3]*f[10]+alphav[2]*f[9]+alphav[1]*f[8]+f[1]*alphav[8]+alphav[0]*f[4]+f[0]*alphav[4]); 
  out[11] += 0.4330127018922193*(alphay[4]*f[13]+alphav[4]*f[12]+(alphav[3]+alphay[2]+alphax[1])*f[11]+alphay[8]*f[10]+alphav[8]*f[9]+(alphav[6]+alphay[5]+alphax[0])*f[7]+alphay[0]*f[6]+f[0]*alphay[6]+alphav[0]*f[5]+f[0]*alphav[5]+alphay[1]*f[3]+f[1]*alphay[3]+alphav[1]*f[2]+f[1]*alphav[2]); 
  out[12] += 0.4330127018922193*(alphay[3]*f[13]+(alphay[2]+alphax[1])*f[12]+alphay[6]*f[10]+(alphay[5]+alphax[0])*f[9]+alphay[0]*f[8]+f[0]*alphay[8]+alphay[1]*f[4]+f[1]*alphay[4]); 
  out[13] += 0.4330127018922193*((alphav[3]+alphax[1])*f[13]+alphav[2]*f[12]+(alphav[6]+alphax[0])*f[10]+alphav[5]*f[9]+alphav[0]*f[8]+f[0]*alphav[8]+alphav[1]*f[4]+f[1]*alphav[4]); 
  out[14] += 0.4330127018922193*((alphav[6]+alphay[5])*f[15]+(alphav[3]+alphay[2])*f[14]+alphay[1]*f[13]+alphav[1]*f[12]+alphay[0]*f[10]+alphav[0]*f[9]+(alphay[6]+alphav[5])*f[8]+f[6]*alphay[8]+f[5]*alphav[8]+(alphay[3]+alphav[2])*f[4]+f[3]*alphay[4]+f[2]*alphav[4]); 
  out[15] += 0.4330127018922193*((alphav[3]+alphay[2]+alphax[1])*f[15]+(alphav[6]+alphay[5]+alphax[0])*f[14]+alphay[0]*f[13]+alphav[0]*f[12]+alphay[1]*f[10]+alphav[1]*f[9]+(alphay[3]+alphav[2])*f[8]+f[3]*alphay[8]+f[2]*alphav[8]+alphay[4]*f[6]+f[4]*alphay[6]+alphav[4]*f[5]+f[4]*alphav[5]); 
  return cflFreq; 
} 
