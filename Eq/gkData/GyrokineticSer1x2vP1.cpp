#include <GyrokineticModDecl.h> 
double GyrokineticVol1x2vSerP1_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *jacobTotInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_z = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wz = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double hamil[8]; 
  hamil[0] = (0.2357022603955158*(3.0*dfac_v2*(2.0*m_*wv2+2.828427124746191*(Bmag[0]*wm+Phi[0]*q_))+2.0*m_))/dfac_v2; 
  hamil[1] = 2.0*Phi[1]*q_; 
  hamil[2] = (1.632993161855453*m_*wv)/dfac_v; 
  hamil[3] = (1.154700538379252*Bmag[0])/dfac_m; 
  double BstarX_by_Bmag[8]; 
  double BstarY_by_Bmag[8]; 
  double BstarZ_by_Bmag[8]; 
  BstarZ_by_Bmag[0] = 1.414213562373095*Gradpar[0]*jacobTotInv[0]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaz[8]; 
  alphaz[0] = (0.6123724356957944*BstarZ_by_Bmag[0]*hamil[2]*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0441941738241592*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0441941738241592*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphav[8]; 
  alphav[0] = -(0.6123724356957944*BstarZ_by_Bmag[0]*hamil[1]*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.0441941738241592*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.0441941738241592*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*alphaz[0]*f[0]; 
  out[2] += 0.6123724356957944*alphav[0]*f[0]; 
  out[4] += 0.6123724356957944*(alphaz[0]*f[2]+alphav[0]*f[1]); 
  out[5] += 0.6123724356957944*alphaz[0]*f[3]; 
  out[6] += 0.6123724356957944*alphav[0]*f[3]; 
  out[7] += 0.6123724356957944*(alphaz[0]*f[6]+alphav[0]*f[5]); 
  return cflFreq; 
} 
double GyrokineticVol1x2vSerP1_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *jacobTotInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_z = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double dfac_m = 2.0/dxv[2]; 
  double wz = w[0]; 
  double wv = w[1]; 
  double wm = w[2]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double hamil[8]; 
  hamil[0] = (0.2357022603955158*(3.0*dfac_v2*(2.0*m_*wv2+2.828427124746191*(Bmag[0]*wm+Phi[0]*q_))+2.0*m_))/dfac_v2; 
  hamil[1] = 2.0*(Bmag[1]*wm+Phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wv)/dfac_v; 
  hamil[3] = (1.154700538379252*Bmag[0])/dfac_m; 
  hamil[5] = (1.154700538379252*Bmag[1])/dfac_m; 
  double BstarX_by_Bmag[8]; 
  double BstarY_by_Bmag[8]; 
  double BstarZ_by_Bmag[8]; 
  BstarZ_by_Bmag[0] = 1.414213562373095*(Gradpar[1]*jacobTotInv[1]+Gradpar[0]*jacobTotInv[0]); 
  BstarZ_by_Bmag[1] = 1.414213562373095*(Gradpar[0]*jacobTotInv[1]+jacobTotInv[0]*Gradpar[1]); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaz[8]; 
  alphaz[0] = (0.6123724356957944*BstarZ_by_Bmag[0]*hamil[2]*dfac_v*dfac_z)/m_; 
  alphaz[1] = (0.6123724356957944*BstarZ_by_Bmag[1]*hamil[2]*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = -0.125*(2.449489742783178*alphaz[1]-1.414213562373095*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(2.449489742783178*alphaz[1]+1.414213562373095*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.3535533905932737*alphaz[0]-0.6123724356957944*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphaz[0]-0.6123724356957944*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphaz[0]-0.6123724356957944*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphaz[0]-0.6123724356957944*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphav[8]; 
  alphav[0] = -(0.6123724356957944*BstarZ_by_Bmag[0]*hamil[1]*dfac_v*dfac_z)/m_; 
  alphav[1] = -(0.6123724356957944*BstarZ_by_Bmag[1]*hamil[1]*dfac_v*dfac_z)/m_; 
  alphav[3] = -(0.6123724356957944*BstarZ_by_Bmag[0]*hamil[5]*dfac_v*dfac_z)/m_; 
  alphav[5] = -(0.6123724356957944*BstarZ_by_Bmag[1]*hamil[5]*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.3535533905932737*alphav[5]-0.3535533905932737*(alphav[3]+alphav[1])+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*(alphav[1]+alphav[0])-0.3535533905932737*(alphav[5]+alphav[3])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.3535533905932737*alphav[5])+0.3535533905932737*alphav[3]-0.3535533905932737*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.04419417382415921*(alphav[5]+alphav[3]+alphav[1]+alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(0.3535533905932737*alphav[5]-0.3535533905932737*(alphav[3]+alphav[1])+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*(alphav[1]+alphav[0])-0.3535533905932737*(alphav[5]+alphav[3])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.3535533905932737*alphav[5])+0.3535533905932737*alphav[3]-0.3535533905932737*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.04419417382415921*(alphav[5]+alphav[3]+alphav[1]+alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*(alphaz[1]*f[1]+alphaz[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphav[5]*f[5]+alphav[3]*f[3]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphav[3]*f[5]+f[3]*alphav[5]+alphaz[1]*f[4]+alphaz[0]*f[2]+alphav[0]*f[1]+f[0]*alphav[1]); 
  out[5] += 0.6123724356957944*(alphaz[1]*f[5]+alphaz[0]*f[3]); 
  out[6] += 0.6123724356957944*(alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[3]+f[0]*alphav[3]); 
  out[7] += 0.6123724356957944*(alphaz[1]*f[7]+alphaz[0]*f[6]+alphav[0]*f[5]+f[0]*alphav[5]+alphav[1]*f[3]+f[1]*alphav[3]); 
  return cflFreq; 
} 
