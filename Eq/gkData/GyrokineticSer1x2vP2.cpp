#include <GyrokineticModDecl.h> 
double GyrokineticVol1x2vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *f, double *out) 
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
  double hamil[20]; 
  hamil[0] = (0.2357022603955158*(3.0*dfac_v2*(2.0*m_*wv2+2.828427124746191*(Bmag[0]*wm+Phi[0]*q_))+2.0*m_))/dfac_v2; 
  hamil[1] = 2.0*Phi[1]*q_; 
  hamil[2] = (1.632993161855453*m_*wv)/dfac_v; 
  hamil[3] = (1.154700538379252*Bmag[0])/dfac_m; 
  hamil[7] = 2.0*Phi[2]*q_; 
  hamil[8] = (0.421637021355784*m_)/dfac_v2; 
  double BstarX_by_Bmag[20]; 
  double BstarY_by_Bmag[20]; 
  double BstarZ_by_Bmag[20]; 
  BstarZ_by_Bmag[0] = 2.0*Gradpar[0]; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaz[20]; 
  alphaz[0] = (0.6123724356957944*BstarZ_by_Bmag[0]*hamil[2]*dfac_v*dfac_z)/m_; 
  alphaz[2] = (1.369306393762915*BstarZ_by_Bmag[0]*hamil[8]*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.3535533905932737*alphaz[0]-0.4743416490252568*alphaz[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphaz[2]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphaz[0]-0.4743416490252568*alphaz[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphaz[2]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphaz[0]-0.4743416490252568*alphaz[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphaz[2]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(0.3535533905932737*alphaz[0]-0.4743416490252568*alphaz[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphaz[2]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alphaz[0]-0.4743416490252568*alphaz[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphaz[2]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alphaz[0]-0.4743416490252568*alphaz[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphaz[2]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphav[20]; 
  alphav[0] = -(0.6123724356957944*BstarZ_by_Bmag[0]*hamil[1]*dfac_v*dfac_z)/m_; 
  alphav[1] = -(1.369306393762915*BstarZ_by_Bmag[0]*hamil[7]*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alphav[0]-0.4743416490252568*alphav[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*(alphaz[2]*f[2]+alphaz[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphav[1]*f[1]+alphav[0]*f[0]); 
  out[4] += 0.1224744871391589*(4.47213595499958*(alphaz[2]*f[8]+alphav[1]*f[7])+5.0*(alphaz[0]*f[2]+f[0]*alphaz[2]+alphav[0]*f[1]+f[0]*alphav[1])); 
  out[5] += 0.6123724356957944*(alphaz[2]*f[6]+alphaz[0]*f[3]); 
  out[6] += 0.6123724356957944*(alphav[1]*f[5]+alphav[0]*f[3]); 
  out[7] += 1.369306393762915*(alphaz[2]*f[4]+alphaz[0]*f[1]); 
  out[8] += 1.369306393762915*(alphav[1]*f[4]+alphav[0]*f[2]); 
  out[10] += 0.07071067811865474*(7.745966692414834*(alphaz[2]*f[14]+alphav[1]*f[13])+8.660254037844386*(alphaz[0]*f[6]+alphav[0]*f[5]+(alphaz[2]+alphav[1])*f[3])); 
  out[11] += 0.07071067811865474*(17.32050807568877*alphaz[2]*f[12]+8.660254037844387*alphav[0]*f[7]+19.36491673103708*alphaz[0]*f[4]+f[1]*(19.36491673103708*alphaz[2]+7.745966692414834*alphav[1])); 
  out[12] += 0.07071067811865474*(17.32050807568877*alphav[1]*f[11]+8.660254037844387*alphaz[0]*f[8]+19.36491673103708*alphav[0]*f[4]+(7.745966692414834*alphaz[2]+19.36491673103708*alphav[1])*f[2]); 
  out[13] += 1.369306393762915*(alphaz[2]*f[10]+alphaz[0]*f[5]); 
  out[14] += 1.369306393762915*(alphav[1]*f[10]+alphav[0]*f[6]); 
  out[15] += 0.07071067811865474*(8.660254037844386*alphaz[2]*f[16]+8.660254037844387*alphaz[0]*f[9]); 
  out[16] += 0.07071067811865474*(8.660254037844386*alphav[1]*f[15]+8.660254037844387*alphav[0]*f[9]); 
  out[17] += 0.07071067811865474*(17.32050807568877*alphaz[2]*f[18]+8.660254037844387*alphav[0]*f[13]+19.36491673103709*alphaz[0]*f[10]+(19.36491673103709*alphaz[2]+7.745966692414834*alphav[1])*f[5]); 
  out[18] += 0.07071067811865474*(17.32050807568877*alphav[1]*f[17]+8.660254037844387*alphaz[0]*f[14]+19.36491673103709*alphav[0]*f[10]+(7.745966692414834*alphaz[2]+19.36491673103709*alphav[1])*f[6]); 
  out[19] += 0.07071067811865474*(8.660254037844387*(alphaz[0]*f[16]+alphav[0]*f[15])+8.660254037844386*(alphaz[2]+alphav[1])*f[9]); 
  return cflFreq; 
} 
double GyrokineticVol1x2vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *f, double *out) 
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
  double hamil[20]; 
  hamil[0] = (0.2357022603955158*(3.0*dfac_v2*(2.0*m_*wv2+2.828427124746191*(Bmag[0]*wm+Phi[0]*q_))+2.0*m_))/dfac_v2; 
  hamil[1] = 2.0*(Bmag[1]*wm+Phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wv)/dfac_v; 
  hamil[3] = (1.154700538379252*Bmag[0])/dfac_m; 
  hamil[5] = (1.154700538379252*Bmag[1])/dfac_m; 
  hamil[7] = 2.0*(Bmag[2]*wm+Phi[2]*q_); 
  hamil[8] = (0.421637021355784*m_)/dfac_v2; 
  hamil[13] = (1.154700538379251*Bmag[2])/dfac_m; 
  double BstarX_by_Bmag[20]; 
  double BstarY_by_Bmag[20]; 
  double BstarZ_by_Bmag[20]; 
  BstarX_by_Bmag[0] = (1.732050807568877*((Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2])*geoY[2]+Bmag[2]*(2.0*geoY[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))*dfac_z*m_*wv)/q_; 
  BstarX_by_Bmag[1] = (0.04948716593053935*((122.9837387624885*Bmag[2]*BmagInv[2]+14.0*(5.0*BmagInv[0]*Bmag[2]+2.23606797749979*Bmag[1]*BmagInv[1]))*geoY[2]+7.0*(2.0*(5.0*geoY[0]*Bmag[2]+2.23606797749979*Bmag[1]*geoY[1])*BmagInv[2]+2.23606797749979*(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+5.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_z*m_*wv)/q_; 
  BstarX_by_Bmag[2] = (((Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2])*geoY[2]+Bmag[2]*(2.0*geoY[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))*dfac_z*m_)/(dfac_v*q_); 
  BstarX_by_Bmag[4] = (0.02857142857142857*((122.9837387624885*Bmag[2]*BmagInv[2]+14.0*(5.0*BmagInv[0]*Bmag[2]+2.23606797749979*Bmag[1]*BmagInv[1]))*geoY[2]+7.0*(2.0*(5.0*geoY[0]*Bmag[2]+2.23606797749979*Bmag[1]*geoY[1])*BmagInv[2]+2.23606797749979*(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+5.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_z*m_)/(dfac_v*q_); 
  BstarX_by_Bmag[7] = (0.04948716593053935*((11.18033988749895*(2.0*Bmag[1]*BmagInv[2]+11.0*BmagInv[1]*Bmag[2])+35.0*BmagInv[0]*Bmag[1])*geoY[2]+(122.9837387624885*geoY[1]*Bmag[2]+35.0*geoY[0]*Bmag[1])*BmagInv[2]+14.0*(5.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+2.23606797749979*Bmag[1]*BmagInv[1]*geoY[1]))*dfac_z*m_*wv)/q_; 
  BstarX_by_Bmag[11] = (0.063887656499994*((5.0*(2.0*Bmag[1]*BmagInv[2]+11.0*BmagInv[1]*Bmag[2])+15.65247584249853*BmagInv[0]*Bmag[1])*geoY[2]+(55.0*geoY[1]*Bmag[2]+15.65247584249853*geoY[0]*Bmag[1])*BmagInv[2]+14.0*(2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*BmagInv[1]*geoY[1]))*dfac_z*m_)/(dfac_v*q_); 
  BstarY_by_Bmag[0] = -(1.732050807568877*((Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2])*geoX[2]+Bmag[2]*(2.0*geoX[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoX[1]+BmagInv[0]*geoX[0]))*dfac_z*m_*wv)/q_; 
  BstarY_by_Bmag[1] = -(0.04948716593053935*((122.9837387624885*Bmag[2]*BmagInv[2]+14.0*(5.0*BmagInv[0]*Bmag[2]+2.23606797749979*Bmag[1]*BmagInv[1]))*geoX[2]+7.0*(2.0*(5.0*geoX[0]*Bmag[2]+2.23606797749979*Bmag[1]*geoX[1])*BmagInv[2]+2.23606797749979*(9.0*BmagInv[1]*geoX[1]+5.0*BmagInv[0]*geoX[0])*Bmag[2]+5.0*Bmag[1]*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1])))*dfac_z*m_*wv)/q_; 
  BstarY_by_Bmag[2] = -(1.0*((Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2])*geoX[2]+Bmag[2]*(2.0*geoX[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoX[1]+BmagInv[0]*geoX[0]))*dfac_z*m_)/(dfac_v*q_); 
  BstarY_by_Bmag[4] = -(0.02857142857142857*((122.9837387624885*Bmag[2]*BmagInv[2]+14.0*(5.0*BmagInv[0]*Bmag[2]+2.23606797749979*Bmag[1]*BmagInv[1]))*geoX[2]+7.0*(2.0*(5.0*geoX[0]*Bmag[2]+2.23606797749979*Bmag[1]*geoX[1])*BmagInv[2]+2.23606797749979*(9.0*BmagInv[1]*geoX[1]+5.0*BmagInv[0]*geoX[0])*Bmag[2]+5.0*Bmag[1]*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1])))*dfac_z*m_)/(dfac_v*q_); 
  BstarY_by_Bmag[7] = -(0.04948716593053935*((11.18033988749895*(2.0*Bmag[1]*BmagInv[2]+11.0*BmagInv[1]*Bmag[2])+35.0*BmagInv[0]*Bmag[1])*geoX[2]+(122.9837387624885*geoX[1]*Bmag[2]+35.0*geoX[0]*Bmag[1])*BmagInv[2]+14.0*(5.0*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1])*Bmag[2]+2.23606797749979*Bmag[1]*BmagInv[1]*geoX[1]))*dfac_z*m_*wv)/q_; 
  BstarY_by_Bmag[11] = -(0.063887656499994*((5.0*(2.0*Bmag[1]*BmagInv[2]+11.0*BmagInv[1]*Bmag[2])+15.65247584249853*BmagInv[0]*Bmag[1])*geoX[2]+(55.0*geoX[1]*Bmag[2]+15.65247584249853*geoX[0]*Bmag[1])*BmagInv[2]+14.0*(2.23606797749979*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1])*Bmag[2]+Bmag[1]*BmagInv[1]*geoX[1]))*dfac_z*m_)/(dfac_v*q_); 
  BstarZ_by_Bmag[0] = 2.0*Gradpar[0]; 
  BstarZ_by_Bmag[1] = 2.0*Gradpar[1]; 
  BstarZ_by_Bmag[7] = 2.0*Gradpar[2]; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaz[20]; 
  alphaz[0] = (0.6123724356957944*BstarZ_by_Bmag[0]*hamil[2]*dfac_v*dfac_z)/m_; 
  alphaz[1] = (0.6123724356957944*BstarZ_by_Bmag[1]*hamil[2]*dfac_v*dfac_z)/m_; 
  alphaz[2] = (1.369306393762915*BstarZ_by_Bmag[0]*hamil[8]*dfac_v*dfac_z)/m_; 
  alphaz[4] = (1.369306393762915*BstarZ_by_Bmag[1]*hamil[8]*dfac_v*dfac_z)/m_; 
  alphaz[7] = (0.6123724356957944*hamil[2]*BstarZ_by_Bmag[7]*dfac_v*dfac_z)/m_; 
  alphaz[11] = (1.369306393762915*BstarZ_by_Bmag[7]*hamil[8]*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.125*(3.16227766016838*alphaz[7]-2.449489742783178*alphaz[1]+1.414213562373095*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.125*(3.16227766016838*alphaz[7]+2.449489742783178*alphaz[1]+1.414213562373095*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*((-1.060660171779821*alphaz[11])+0.7905694150420947*alphaz[7]+0.8215838362577489*alphaz[4]-0.4743416490252568*alphaz[2]-0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.7905694150420947*alphaz[7]-0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.060660171779821*alphaz[11]+0.7905694150420947*alphaz[7]-0.8215838362577489*alphaz[4]+0.4743416490252568*alphaz[2]-0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-1.060660171779821*alphaz[11])+0.7905694150420947*alphaz[7]+0.8215838362577489*alphaz[4]-0.4743416490252568*alphaz[2]-0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.060660171779821*alphaz[11]+0.7905694150420947*alphaz[7]-0.8215838362577489*alphaz[4]+0.4743416490252568*alphaz[2]-0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-1.060660171779821*alphaz[11])+0.7905694150420947*alphaz[7]+0.8215838362577489*alphaz[4]-0.4743416490252568*alphaz[2]-0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.7905694150420947*alphaz[7]-0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.060660171779821*alphaz[11]+0.7905694150420947*alphaz[7]-0.8215838362577489*alphaz[4]+0.4743416490252568*alphaz[2]-0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*((-1.060660171779821*alphaz[11])+0.7905694150420947*alphaz[7]-0.8215838362577489*alphaz[4]-0.4743416490252568*alphaz[2]+0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.7905694150420947*alphaz[7]+0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(1.060660171779821*alphaz[11]+0.7905694150420947*alphaz[7]+0.8215838362577489*alphaz[4]+0.4743416490252568*alphaz[2]+0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-1.060660171779821*alphaz[11])+0.7905694150420947*alphaz[7]-0.8215838362577489*alphaz[4]-0.4743416490252568*alphaz[2]+0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(1.060660171779821*alphaz[11]+0.7905694150420947*alphaz[7]+0.8215838362577489*alphaz[4]+0.4743416490252568*alphaz[2]+0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-1.060660171779821*alphaz[11])+0.7905694150420947*alphaz[7]-0.8215838362577489*alphaz[4]-0.4743416490252568*alphaz[2]+0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.7905694150420947*alphaz[7]+0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(1.060660171779821*alphaz[11]+0.7905694150420947*alphaz[7]+0.8215838362577489*alphaz[4]+0.4743416490252568*alphaz[2]+0.6123724356957944*alphaz[1]+0.3535533905932737*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphav[20]; 
  alphav[0] = -(0.6123724356957944*(2.23606797749979*BstarZ_by_Bmag[1]*hamil[7]+BstarZ_by_Bmag[0]*hamil[1])*dfac_v*dfac_z)/m_; 
  alphav[1] = -(0.6123724356957944*((2.0*BstarZ_by_Bmag[7]+2.23606797749979*BstarZ_by_Bmag[0])*hamil[7]+BstarZ_by_Bmag[1]*hamil[1])*dfac_v*dfac_z)/m_; 
  alphav[3] = -(0.3535533905932737*(3.872983346207417*BstarZ_by_Bmag[1]*hamil[13]+1.732050807568877*BstarZ_by_Bmag[0]*hamil[5])*dfac_v*dfac_z)/m_; 
  alphav[5] = -(0.3535533905932737*(3.872983346207417*(0.8944271909999159*BstarZ_by_Bmag[7]+BstarZ_by_Bmag[0])*hamil[13]+1.732050807568877*BstarZ_by_Bmag[1]*hamil[5])*dfac_v*dfac_z)/m_; 
  alphav[7] = -(0.6123724356957944*(2.0*BstarZ_by_Bmag[1]*hamil[7]+hamil[1]*BstarZ_by_Bmag[7])*dfac_v*dfac_z)/m_; 
  alphav[13] = -(0.3535533905932737*(3.464101615137754*BstarZ_by_Bmag[1]*hamil[13]+1.732050807568877*hamil[5]*BstarZ_by_Bmag[7])*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.1767766952966368*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.1767766952966368*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.125*((-0.4242640687119285*alphav[13])+0.3162277660168379*alphav[7]+0.6363961030678926*alphav[5]-0.4743416490252568*(alphav[3]+alphav[1])+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.5303300858899104*alphav[13]-0.3952847075210473*alphav[7]-0.4743416490252568*alphav[3]+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.4242640687119285*alphav[13])+0.3162277660168379*alphav[7]-0.6363961030678926*alphav[5]-0.4743416490252568*alphav[3]+0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3162277660168379*alphav[7]-0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3162277660168379*alphav[7]+0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4242640687119285*alphav[13]+0.3162277660168379*alphav[7]-0.6363961030678926*alphav[5]+0.4743416490252568*alphav[3]-0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.5303300858899104*alphav[13])-0.3952847075210473*alphav[7]+0.4743416490252568*alphav[3]+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4242640687119285*alphav[13]+0.3162277660168379*alphav[7]+0.6363961030678926*alphav[5]+0.4743416490252568*(alphav[3]+alphav[1])+0.3535533905932737*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.125*((-0.4242640687119285*alphav[13])+0.3162277660168379*alphav[7]+0.6363961030678926*alphav[5]-0.4743416490252568*(alphav[3]+alphav[1])+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.5303300858899104*alphav[13]-0.3952847075210473*alphav[7]-0.4743416490252568*alphav[3]+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.4242640687119285*alphav[13])+0.3162277660168379*alphav[7]-0.6363961030678926*alphav[5]-0.4743416490252568*alphav[3]+0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3162277660168379*alphav[7]-0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3162277660168379*alphav[7]+0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4242640687119285*alphav[13]+0.3162277660168379*alphav[7]-0.6363961030678926*alphav[5]+0.4743416490252568*alphav[3]-0.4743416490252568*alphav[1]+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.5303300858899104*alphav[13])-0.3952847075210473*alphav[7]+0.4743416490252568*alphav[3]+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4242640687119285*alphav[13]+0.3162277660168379*alphav[7]+0.6363961030678926*alphav[5]+0.4743416490252568*(alphav[3]+alphav[1])+0.3535533905932737*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*(alphaz[11]*f[11]+alphaz[7]*f[7]+alphaz[4]*f[4]+alphaz[2]*f[2]+alphaz[1]*f[1]+alphaz[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphav[13]*f[13]+alphav[7]*f[7]+alphav[5]*f[5]+alphav[3]*f[3]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[4] += 0.07071067811865474*(7.745966692414834*(alphav[5]*f[13]+f[5]*alphav[13]+alphaz[4]*f[12])+8.660254037844387*(alphaz[7]*f[11]+f[7]*alphaz[11])+7.745966692414834*(alphaz[2]*f[8]+alphav[1]*f[7]+f[1]*alphav[7])+8.660254037844386*(alphav[3]*f[5]+f[3]*alphav[5]+alphaz[1]*f[4]+f[1]*alphaz[4]+alphaz[0]*f[2]+f[0]*alphaz[2]+alphav[0]*f[1]+f[0]*alphav[1])); 
  out[5] += 0.07071067811865474*(8.660254037844387*(alphaz[11]*f[17]+alphaz[7]*f[13])+8.660254037844386*(alphaz[4]*f[10]+alphaz[2]*f[6]+alphaz[1]*f[5]+alphaz[0]*f[3])); 
  out[6] += 0.07071067811865474*(7.745966692414834*alphav[5]*f[15]+8.660254037844387*(alphav[7]*f[13]+f[7]*alphav[13])+7.745966692414834*alphav[3]*f[9]+8.660254037844386*(alphav[1]*f[5]+f[1]*alphav[5]+alphav[0]*f[3]+f[0]*alphav[3])); 
  out[7] += 0.07071067811865474*(17.32050807568877*(alphaz[4]*f[11]+f[4]*alphaz[11])+17.32050807568877*(alphaz[1]*f[7]+f[1]*alphaz[7])+19.36491673103709*(alphaz[2]*f[4]+f[2]*alphaz[4]+alphaz[0]*f[1]+f[0]*alphaz[1])); 
  out[8] += 1.369306393762915*(alphav[13]*f[17]+alphav[7]*f[11]+alphav[5]*f[10]+alphav[3]*f[6]+alphav[1]*f[4]+alphav[0]*f[2]); 
  out[10] += 0.07071067811865474*(7.745966692414834*alphaz[4]*f[18]+8.660254037844386*alphaz[7]*f[17]+6.928203230275509*alphav[13]*f[15]+7.745966692414834*(alphav[3]*f[15]+alphaz[2]*f[14])+8.660254037844386*alphaz[11]*f[13]+7.745966692414834*(alphav[1]*f[13]+f[1]*alphav[13])+8.660254037844386*alphaz[1]*f[10]+7.745966692414834*(alphav[5]*(f[9]+f[7])+f[5]*alphav[7])+8.660254037844386*(alphaz[0]*f[6]+(alphaz[4]+alphav[0])*f[5]+f[0]*alphav[5]+(alphaz[2]+alphav[1])*f[3]+f[1]*alphav[3])); 
  out[11] += 0.01010152544552211*(38.72983346207417*alphav[13]*f[13]+60.6217782649107*(alphav[3]*f[13]+f[3]*alphav[13])+108.4435336938077*alphaz[11]*f[12]+121.2435565298214*(alphaz[2]*f[12]+alphaz[1]*f[11]+f[1]*alphaz[11])+121.2435565298214*alphaz[4]*f[8]+(38.72983346207417*alphav[7]+121.2435565298214*alphaz[4]+60.62177826491071*alphav[0])*f[7]+121.2435565298214*f[4]*alphaz[7]+60.62177826491071*f[0]*alphav[7]+54.22176684690384*alphav[5]*f[5]+135.5544171172596*(alphaz[0]*f[4]+f[0]*alphaz[4]+alphaz[1]*f[2])+f[1]*(135.5544171172596*alphaz[2]+54.22176684690384*alphav[1])); 
  out[12] += 0.07071067811865474*(17.32050807568877*alphav[5]*f[17]+17.32050807568877*f[10]*alphav[13]+8.660254037844386*alphaz[1]*f[12]+(7.745966692414834*alphaz[11]+17.32050807568877*alphav[1])*f[11]+19.36491673103708*alphav[3]*f[10]+8.660254037844387*alphaz[0]*f[8]+17.32050807568877*f[4]*alphav[7]+19.36491673103708*alphav[5]*f[6]+(7.745966692414834*alphaz[4]+19.36491673103708*alphav[0])*f[4]+(7.745966692414834*alphaz[2]+19.36491673103708*alphav[1])*f[2]); 
  out[13] += 0.07071067811865474*(17.32050807568877*alphaz[4]*f[17]+17.32050807568877*alphaz[1]*f[13]+f[10]*(17.32050807568877*alphaz[11]+19.36491673103708*alphaz[2])+17.32050807568877*f[5]*alphaz[7]+19.36491673103708*(alphaz[4]*f[6]+alphaz[0]*f[5]+alphaz[1]*f[3])); 
  out[14] += 0.07071067811865474*(17.32050807568877*alphav[5]*f[19]+19.36491673103708*alphav[7]*f[17]+17.32050807568877*alphav[3]*f[16]+19.36491673103708*(f[11]*alphav[13]+alphav[1]*f[10]+alphav[0]*f[6]+f[4]*alphav[5]+f[2]*alphav[3])); 
  out[15] += 0.07071067811865474*(8.660254037844387*alphaz[4]*f[19]+8.660254037844386*(alphaz[2]*f[16]+alphaz[1]*f[15])+8.660254037844387*alphaz[0]*f[9]); 
  out[16] += 0.07071067811865474*(8.660254037844386*alphav[1]*f[15]+7.745966692414834*alphav[13]*f[13]+8.660254037844387*alphav[0]*f[9]+7.745966692414834*(alphav[5]*f[5]+alphav[3]*f[3])); 
  out[17] += 0.002020305089104421*(542.2176684690384*alphaz[11]*f[18]+606.217782649107*(alphaz[2]*f[18]+alphaz[1]*f[17])+242.4871130596428*alphav[5]*f[15]+606.2177826491072*alphaz[4]*f[14]+(193.6491673103708*alphav[7]+606.2177826491072*alphaz[4]+303.1088913245536*alphav[0])*f[13]+(271.1088342345192*f[9]+193.6491673103708*f[7]+303.1088913245536*f[0])*alphav[13]+606.2177826491072*f[5]*alphaz[11]+(606.217782649107*alphaz[7]+677.7720855862981*alphaz[0])*f[10]+303.1088913245535*(alphav[3]*f[7]+f[3]*alphav[7])+677.7720855862981*(alphaz[1]*f[6]+alphaz[2]*f[5])+271.1088342345192*(alphav[1]*f[5]+f[1]*alphav[5])+677.7720855862981*f[3]*alphaz[4]); 
  out[18] += 0.07071067811865474*((15.49193338482967*alphav[13]+17.32050807568877*alphav[3])*f[19]+8.660254037844386*alphaz[1]*f[18]+(7.745966692414834*alphaz[11]+17.32050807568877*alphav[1])*f[17]+17.32050807568877*alphav[5]*f[16]+8.660254037844387*alphaz[0]*f[14]+17.32050807568877*(f[4]*alphav[13]+alphav[5]*f[11])+(17.32050807568877*alphav[7]+7.745966692414834*alphaz[4]+19.36491673103709*alphav[0])*f[10]+7.745966692414834*alphaz[2]*f[6]+19.36491673103709*(alphav[1]*f[6]+f[2]*alphav[5]+alphav[3]*f[4])); 
  out[19] += 0.01414213562373095*(43.30127018922193*alphaz[1]*f[19]+43.30127018922195*alphaz[0]*f[16]+(38.72983346207417*alphav[7]+43.30127018922195*(alphaz[4]+alphav[0]))*f[15]+34.64101615137755*(alphav[5]*f[13]+f[5]*alphav[13])+43.30127018922193*(alphaz[2]+alphav[1])*f[9]+38.72983346207418*(alphav[3]*f[5]+f[3]*alphav[5])); 
  return cflFreq; 
} 
