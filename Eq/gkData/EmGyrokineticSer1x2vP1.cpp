#include <GyrokineticModDecl.h> 
double EmGyrokineticVol1x2vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double* dApardt, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = (1.732050807568877*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2+2.0*Gradpar[0]*dfac_x)*wv; 
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
  alphav[0] = -(1.0*Phi[1]*dfac_v*(2.121320343559642*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2+2.449489742783178*Gradpar[0]*dfac_x)*q_)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -(0.125*(2.828427124746191*dApardt[0]*dfac_v*q_-1.414213562373095*alphav[0]*m_))/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = -(0.125*(2.828427124746191*dApardt[0]*dfac_v*q_-1.414213562373095*alphav[0]*m_))/m_; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.3535533905932737*alphav[0]-(1.0*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*alphav[0]-(0.7071067811865475*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*alphav[0]-(1.0*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.3535533905932737*alphav[0]-(0.7071067811865475*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(0.3535533905932737*alphav[0]-(1.0*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.3535533905932737*alphav[0]-(0.7071067811865475*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.3535533905932737*alphav[0]-(1.0*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.3535533905932737*alphav[0]-(0.7071067811865475*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_); 
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
  return 0.0; 
} 
double EmGyrokineticVol1x2vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *Apar, const double* dApardt, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_x = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wx = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double dfac_x2 = dfac_x*dfac_x; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = -1.0*((1.732050807568877*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_*wv2)/q_+(1.732050807568877*(0.7071067811865475*Bmag[1]*(BmagInv[0]*Apar[1]*geoY[1]+geoY[0]*(Apar[1]*BmagInv[1]+Apar[0]*BmagInv[0]))+(0.7071067811865475*Apar[0]*Bmag[1]-1.0*Apar[1])*BmagInv[1]*geoY[1])*dfac_x2-1.0*(1.732050807568877*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2+2.0*Gradpar[0]*dfac_x))*wv); 
  alphax[1] = -1.0*((1.732050807568877*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_*wv2)/q_+(1.732050807568877*(BmagInv[0]*(0.7071067811865475*Apar[0]*Bmag[1]-1.0*Apar[1])*geoY[1]+geoY[0]*(0.7071067811865475*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(0.7071067811865475*BmagInv[0]*Bmag[1]-1.0*BmagInv[1])))*dfac_x2-1.0*(2.0*Gradpar[1]*dfac_x-2.204540768504859*Apar[1]*Bmag[1]*BmagInv[1]*geoY[1]*dfac_x2))*wv); 
  alphax[2] = -(1.0*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*m_*wv)/(dfac_v*q_); 
  alphax[4] = -(1.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*m_*wv)/(dfac_v*q_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(2.449489742783178*alphax[1]-1.414213562373095*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.6123724356957944*alphax[4]-0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.6123724356957944*alphax[4])+0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*(0.6123724356957944*alphax[4]-0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-0.6123724356957944*alphax[4])+0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*((-0.6123724356957944*alphax[4])-0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.6123724356957944*alphax[4]+0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*((-0.6123724356957944*alphax[4])-0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*(0.6123724356957944*alphax[4]+0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  double alphav[8]; 
  alphav[0] = dfac_v*(2.121320343559642*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*((Bmag[1]*wm)/q_+Phi[1])*wv+(((1.5*Bmag[1]*(BmagInv[0]*Apar[1]*geoY[1]+geoY[0]*(Apar[1]*BmagInv[1]+Apar[0]*BmagInv[0]))+(1.5*Apar[0]*Bmag[1]-2.121320343559642*Apar[1])*BmagInv[1]*geoY[1])*dfac_x2-1.0*(2.121320343559642*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2+2.449489742783178*Gradpar[0]*dfac_x))*(Bmag[1]*wm+(Phi[1]*q2)/q_))/m_); 
  alphav[1] = dfac_v*(2.121320343559642*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*((Bmag[1]*wm)/q_+Phi[1])*wv+((((2.7*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(1.5*Apar[0]*Bmag[1]-2.121320343559642*Apar[1]))*geoY[1]+geoY[0]*(1.5*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(1.5*BmagInv[0]*Bmag[1]-2.121320343559642*BmagInv[1])))*dfac_x2-2.449489742783178*Gradpar[1]*dfac_x)*(Bmag[1]*wm+(Phi[1]*q2)/q_))/m_); 
  alphav[2] = 1.224744871391589*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*((Bmag[1]*wm)/q_+Phi[1]); 
  alphav[3] = (Bmag[1]*dfac_v*((1.224744871391589*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2*wv)/q_+(1.732050807568877*(0.5*Bmag[1]*(BmagInv[0]*Apar[1]*geoY[1]+geoY[0]*(Apar[1]*BmagInv[1]+Apar[0]*BmagInv[0]))+(0.5*Apar[0]*Bmag[1]-0.7071067811865475*Apar[1])*BmagInv[1]*geoY[1])*dfac_x2-1.0*(1.224744871391589*BmagInv[0]*geoY[0]*Apar[1]*dfac_x2+1.414213562373095*Gradpar[0]*dfac_x))/m_))/dfac_m; 
  alphav[4] = 1.224744871391589*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*((Bmag[1]*wm)/q_+Phi[1]); 
  alphav[5] = (Bmag[1]*dfac_v*((1.224744871391589*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2*wv)/q_+(1.732050807568877*((0.9*Apar[1]*Bmag[1]*BmagInv[1]+BmagInv[0]*(0.5*Apar[0]*Bmag[1]-0.7071067811865475*Apar[1]))*geoY[1]+geoY[0]*(0.5*Apar[0]*Bmag[1]*BmagInv[1]+Apar[1]*(0.5*BmagInv[0]*Bmag[1]-0.7071067811865475*BmagInv[1])))*dfac_x2-1.414213562373095*Gradpar[1]*dfac_x)/m_))/dfac_m; 
  alphav[6] = (0.7071067811865475*Bmag[1]*Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0])*dfac_x2)/(dfac_m*q_); 
  alphav[7] = (0.7071067811865475*Bmag[1]*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*dfac_x2)/(dfac_m*q_); 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -(0.125*(2.828427124746191*dApardt[0]*dfac_v*q_+(2.449489742783178*alphav[2]-1.414213562373095*alphav[0])*m_))/m_; 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate surface-averaged alpha on right 
  alphaR = -(0.125*(2.828427124746191*dApardt[0]*dfac_v*q_+((-2.449489742783178*alphav[2])-1.414213562373095*alphav[0])*m_))/m_; 
  if(alphaR>0) cflFreq += alphaR; 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*((-(1.0*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_)-0.6123724356957944*alphav[7]+0.6123724356957944*alphav[6]+0.3535533905932737*alphav[5]+0.6123724356957944*alphav[4]-0.3535533905932737*alphav[3]-0.6123724356957944*alphav[2]-0.3535533905932737*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-(0.7071067811865475*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_)+0.6123724356957944*(alphav[7]+alphav[6])-0.3535533905932737*alphav[5]-0.6123724356957944*alphav[4]-0.3535533905932737*alphav[3]-0.6123724356957944*alphav[2]+0.3535533905932737*(alphav[1]+alphav[0])); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-(1.0*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_)+0.6123724356957944*alphav[7]-0.6123724356957944*alphav[6]-0.3535533905932737*alphav[5]+0.6123724356957944*alphav[4]+0.3535533905932737*alphav[3]-0.6123724356957944*alphav[2]-0.3535533905932737*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaL<0) cflFreq += -alphaL; 
  alphaL = 0.125*((-(0.7071067811865475*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_)-0.6123724356957944*(alphav[7]+alphav[6])+0.3535533905932737*alphav[5]-0.6123724356957944*alphav[4]+0.3535533905932737*alphav[3]-0.6123724356957944*alphav[2]+0.3535533905932737*(alphav[1]+alphav[0])); 
  if(alphaL<0) cflFreq += -alphaL; 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*((-(1.0*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_)+0.6123724356957944*alphav[7]-0.6123724356957944*alphav[6]+0.3535533905932737*alphav[5]-0.6123724356957944*alphav[4]-0.3535533905932737*alphav[3]+0.6123724356957944*alphav[2]-0.3535533905932737*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*((-(0.7071067811865475*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_)-0.6123724356957944*(alphav[7]+alphav[6])-0.3535533905932737*alphav[5]+0.6123724356957944*alphav[4]-0.3535533905932737*alphav[3]+0.6123724356957944*alphav[2]+0.3535533905932737*(alphav[1]+alphav[0])); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*((-(1.0*(0.7071067811865475*dApardt[0]-0.7071067811865475*dApardt[1])*dfac_v*q_)/m_)-0.6123724356957944*alphav[7]+0.6123724356957944*alphav[6]-0.3535533905932737*alphav[5]-0.6123724356957944*alphav[4]+0.3535533905932737*alphav[3]+0.6123724356957944*alphav[2]-0.3535533905932737*alphav[1]+0.3535533905932737*alphav[0]); 
  if(alphaR>0) cflFreq += alphaR; 
  alphaR = 0.125*((-(0.7071067811865475*(dApardt[1]+dApardt[0])*dfac_v*q_)/m_)+0.6123724356957944*(alphav[7]+alphav[6])+0.3535533905932737*alphav[5]+0.6123724356957944*alphav[4]+0.3535533905932737*alphav[3]+0.6123724356957944*alphav[2]+0.3535533905932737*(alphav[1]+alphav[0])); 
  if(alphaR>0) cflFreq += alphaR; 
#endif 

  out[1] += 0.6123724356957944*(alphax[4]*f[4]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphav[7]*f[7]+alphav[6]*f[6]+alphav[5]*f[5]+alphav[4]*f[4]+alphav[3]*f[3]+alphav[2]*f[2]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphav[6]*f[7]+f[6]*alphav[7]+alphav[3]*f[5]+f[3]*alphav[5]+(alphav[2]+alphax[1])*f[4]+f[1]*alphax[4]+f[2]*(alphav[4]+alphax[0])+f[0]*alphax[2]+alphav[0]*f[1]+f[0]*alphav[1]); 
  out[5] += 0.6123724356957944*(alphax[4]*f[7]+alphax[2]*f[6]+alphax[1]*f[5]+alphax[0]*f[3]); 
  out[6] += 0.6123724356957944*(alphav[4]*f[7]+f[4]*alphav[7]+alphav[2]*f[6]+f[2]*alphav[6]+alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[3]+f[0]*alphav[3]); 
  out[7] += 0.6123724356957944*((alphav[2]+alphax[1])*f[7]+f[2]*alphav[7]+(alphav[4]+alphax[0])*f[6]+f[4]*alphav[6]+(alphax[4]+alphav[0])*f[5]+f[0]*alphav[5]+(alphax[2]+alphav[1])*f[3]+f[1]*alphav[3]); 
  return cflFreq; 
} 
