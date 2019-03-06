#include <GyrokineticModDecl.h> 
double GyrokineticVol1x1vSerP2_Bvars_0(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_z = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wz = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double hamil[8]; 
  hamil[0] = (0.3333333333333333*(3.0*dfac_v2*(m_*wv2+1.414213562373095*Phi[0]*q_)+m_))/dfac_v2; 
  hamil[1] = 1.414213562373095*Phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wv)/dfac_v; 
  hamil[4] = 1.414213562373095*Phi[2]*q_; 
  hamil[5] = (0.2981423969999719*m_)/dfac_v2; 
  double BstarX_by_Bmag[8]; 
  double BstarY_by_Bmag[8]; 
  double BstarZ_by_Bmag[8]; 
  BstarZ_by_Bmag[0] = 1.414213562373095*Gradpar[0]; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaz[8]; 
  alphaz[0] = (0.8660254037844386*BstarZ_by_Bmag[0]*hamil[2]*dfac_v*dfac_z)/m_; 
  alphaz[2] = (1.936491673103709*BstarZ_by_Bmag[0]*hamil[5]*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.5*alphaz[0]-0.6708203932499369*alphaz[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.6708203932499369*alphaz[2]+0.5*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(0.5*alphaz[0]-0.6708203932499369*alphaz[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.6708203932499369*alphaz[2]+0.5*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphav[8]; 
  alphav[0] = -(0.8660254037844386*BstarZ_by_Bmag[0]*hamil[1]*dfac_v*dfac_z)/m_; 
  alphav[1] = -(1.936491673103709*BstarZ_by_Bmag[0]*hamil[4]*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.5*alphav[0]-0.6708203932499369*alphav[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(0.5*alphav[0]-0.6708203932499369*alphav[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.8660254037844386*(alphaz[2]*f[2]+alphaz[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphav[1]*f[1]+alphav[0]*f[0]); 
  out[3] += 0.1732050807568877*(4.47213595499958*(alphaz[2]*f[5]+alphav[1]*f[4])+5.0*(alphaz[0]*f[2]+f[0]*alphaz[2]+alphav[0]*f[1]+f[0]*alphav[1])); 
  out[4] += 1.936491673103709*(alphaz[2]*f[3]+alphaz[0]*f[1]); 
  out[5] += 1.936491673103709*(alphav[1]*f[3]+alphav[0]*f[2]); 
  out[6] += 0.1*(17.32050807568877*alphaz[2]*f[7]+8.660254037844387*alphav[0]*f[4]+19.36491673103708*alphaz[0]*f[3]+f[1]*(19.36491673103708*alphaz[2]+7.745966692414834*alphav[1])); 
  out[7] += 0.1*(17.32050807568877*alphav[1]*f[6]+8.660254037844387*alphaz[0]*f[5]+19.36491673103708*alphav[0]*f[3]+(7.745966692414834*alphaz[2]+19.36491673103708*alphav[1])*f[2]); 
  return cflFreq; 
} 
double GyrokineticVol1x1vSerP2_Bvars_1(const double q_, const double m_, const double *w, const double *dxv, const double *Bmag, const double *BmagInv, const double *Gradpar, const double *geoX, const double *geoY, const double *geoZ, const double *Phi, const double *f, double *out) 
{ 
// w[NDIM]: Cell-center coordinates. dxv[NDIM]: Cell spacing. H/f: Input Hamiltonian/distribution function. out: Incremented output 
  double dfac_z = 2.0/dxv[0]; 
  double dfac_v = 2.0/dxv[1]; 
  double wz = w[0]; 
  double wv = w[1]; 
  double wv2 = wv*wv; 
  double dfac_v2 = dfac_v*dfac_v; 
  double q2 = q_*q_; 
  double m2 = m_*m_; 
  double hamil[8]; 
  hamil[0] = (0.3333333333333333*(3.0*dfac_v2*(m_*wv2+1.414213562373095*Phi[0]*q_)+m_))/dfac_v2; 
  hamil[1] = 1.414213562373095*Phi[1]*q_; 
  hamil[2] = (1.154700538379252*m_*wv)/dfac_v; 
  hamil[4] = 1.414213562373095*Phi[2]*q_; 
  hamil[5] = (0.2981423969999719*m_)/dfac_v2; 
  double BstarX_by_Bmag[8]; 
  double BstarY_by_Bmag[8]; 
  double BstarZ_by_Bmag[8]; 
  BstarX_by_Bmag[0] = (1.224744871391589*((Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2])*geoY[2]+Bmag[2]*(2.0*geoY[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))*dfac_z*m_*wv)/q_; 
  BstarX_by_Bmag[1] = (0.03499271061118826*((122.9837387624885*Bmag[2]*BmagInv[2]+7.0*(10.0*BmagInv[0]*Bmag[2]+4.47213595499958*Bmag[1]*BmagInv[1]))*geoY[2]+7.0*(2.0*(5.0*geoY[0]*Bmag[2]+2.23606797749979*Bmag[1]*geoY[1])*BmagInv[2]+2.23606797749979*(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+5.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_z*m_*wv)/q_; 
  BstarX_by_Bmag[2] = (0.7071067811865475*((Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2])*geoY[2]+Bmag[2]*(2.0*geoY[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoY[1]+BmagInv[0]*geoY[0]))*dfac_z*m_)/(dfac_v*q_); 
  BstarX_by_Bmag[3] = (0.02020305089104421*((122.9837387624885*Bmag[2]*BmagInv[2]+7.0*(10.0*BmagInv[0]*Bmag[2]+4.47213595499958*Bmag[1]*BmagInv[1]))*geoY[2]+7.0*(2.0*(5.0*geoY[0]*Bmag[2]+2.23606797749979*Bmag[1]*geoY[1])*BmagInv[2]+2.23606797749979*(9.0*BmagInv[1]*geoY[1]+5.0*BmagInv[0]*geoY[0])*Bmag[2]+5.0*Bmag[1]*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])))*dfac_z*m_)/(dfac_v*q_); 
  BstarX_by_Bmag[4] = (0.03499271061118826*((11.18033988749895*(2.0*Bmag[1]*BmagInv[2]+11.0*BmagInv[1]*Bmag[2])+35.0*BmagInv[0]*Bmag[1])*geoY[2]+(122.9837387624885*geoY[1]*Bmag[2]+35.0*geoY[0]*Bmag[1])*BmagInv[2]+7.0*(10.0*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+4.47213595499958*Bmag[1]*BmagInv[1]*geoY[1]))*dfac_z*m_*wv)/q_; 
  BstarX_by_Bmag[6] = (0.04517539514526257*((5.0*(2.0*Bmag[1]*BmagInv[2]+11.0*BmagInv[1]*Bmag[2])+15.65247584249853*BmagInv[0]*Bmag[1])*geoY[2]+(55.0*geoY[1]*Bmag[2]+15.65247584249853*geoY[0]*Bmag[1])*BmagInv[2]+14.0*(2.23606797749979*(BmagInv[0]*geoY[1]+geoY[0]*BmagInv[1])*Bmag[2]+Bmag[1]*BmagInv[1]*geoY[1]))*dfac_z*m_)/(dfac_v*q_); 
  BstarY_by_Bmag[0] = -(1.224744871391589*((Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2])*geoX[2]+Bmag[2]*(2.0*geoX[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoX[1]+BmagInv[0]*geoX[0]))*dfac_z*m_*wv)/q_; 
  BstarY_by_Bmag[1] = -(0.03499271061118826*((122.9837387624885*Bmag[2]*BmagInv[2]+7.0*(10.0*BmagInv[0]*Bmag[2]+4.47213595499958*Bmag[1]*BmagInv[1]))*geoX[2]+7.0*(2.0*(5.0*geoX[0]*Bmag[2]+2.23606797749979*Bmag[1]*geoX[1])*BmagInv[2]+2.23606797749979*(9.0*BmagInv[1]*geoX[1]+5.0*BmagInv[0]*geoX[0])*Bmag[2]+5.0*Bmag[1]*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1])))*dfac_z*m_*wv)/q_; 
  BstarY_by_Bmag[2] = -(0.7071067811865475*((Bmag[1]*BmagInv[2]+2.0*BmagInv[1]*Bmag[2])*geoX[2]+Bmag[2]*(2.0*geoX[1]*BmagInv[2]+2.23606797749979*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1]))+Bmag[1]*(BmagInv[1]*geoX[1]+BmagInv[0]*geoX[0]))*dfac_z*m_)/(dfac_v*q_); 
  BstarY_by_Bmag[3] = -(0.02020305089104421*((122.9837387624885*Bmag[2]*BmagInv[2]+7.0*(10.0*BmagInv[0]*Bmag[2]+4.47213595499958*Bmag[1]*BmagInv[1]))*geoX[2]+7.0*(2.0*(5.0*geoX[0]*Bmag[2]+2.23606797749979*Bmag[1]*geoX[1])*BmagInv[2]+2.23606797749979*(9.0*BmagInv[1]*geoX[1]+5.0*BmagInv[0]*geoX[0])*Bmag[2]+5.0*Bmag[1]*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1])))*dfac_z*m_)/(dfac_v*q_); 
  BstarY_by_Bmag[4] = -(0.03499271061118826*((11.18033988749895*(2.0*Bmag[1]*BmagInv[2]+11.0*BmagInv[1]*Bmag[2])+35.0*BmagInv[0]*Bmag[1])*geoX[2]+(122.9837387624885*geoX[1]*Bmag[2]+35.0*geoX[0]*Bmag[1])*BmagInv[2]+7.0*(10.0*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1])*Bmag[2]+4.47213595499958*Bmag[1]*BmagInv[1]*geoX[1]))*dfac_z*m_*wv)/q_; 
  BstarY_by_Bmag[6] = -(0.04517539514526257*((5.0*(2.0*Bmag[1]*BmagInv[2]+11.0*BmagInv[1]*Bmag[2])+15.65247584249853*BmagInv[0]*Bmag[1])*geoX[2]+(55.0*geoX[1]*Bmag[2]+15.65247584249853*geoX[0]*Bmag[1])*BmagInv[2]+14.0*(2.23606797749979*(BmagInv[0]*geoX[1]+geoX[0]*BmagInv[1])*Bmag[2]+Bmag[1]*BmagInv[1]*geoX[1]))*dfac_z*m_)/(dfac_v*q_); 
  BstarZ_by_Bmag[0] = 1.414213562373095*Gradpar[0]; 
  BstarZ_by_Bmag[1] = 1.414213562373095*Gradpar[1]; 
  BstarZ_by_Bmag[4] = 1.414213562373095*Gradpar[2]; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphaz[8]; 
  alphaz[0] = (0.8660254037844386*BstarZ_by_Bmag[0]*hamil[2]*dfac_v*dfac_z)/m_; 
  alphaz[1] = (0.8660254037844386*BstarZ_by_Bmag[1]*hamil[2]*dfac_v*dfac_z)/m_; 
  alphaz[2] = (1.936491673103709*BstarZ_by_Bmag[0]*hamil[5]*dfac_v*dfac_z)/m_; 
  alphaz[3] = (1.936491673103709*BstarZ_by_Bmag[1]*hamil[5]*dfac_v*dfac_z)/m_; 
  alphaz[4] = (0.8660254037844386*hamil[2]*BstarZ_by_Bmag[4]*dfac_v*dfac_z)/m_; 
  alphaz[6] = (1.936491673103709*BstarZ_by_Bmag[4]*hamil[5]*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*(2.23606797749979*alphaz[4]-1.732050807568877*alphaz[1]+alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*(2.23606797749979*alphaz[4]+1.732050807568877*alphaz[1]+alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*((-1.5*alphaz[6])+1.118033988749895*alphaz[4]+1.161895003862225*alphaz[3]-0.6708203932499369*alphaz[2]-0.8660254037844386*alphaz[1]+0.5*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(1.118033988749895*alphaz[4]-0.8660254037844386*alphaz[1]+0.5*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(1.5*alphaz[6]+1.118033988749895*alphaz[4]-1.161895003862225*alphaz[3]+0.6708203932499369*alphaz[2]-0.8660254037844386*alphaz[1]+0.5*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*((-1.5*alphaz[6])+1.118033988749895*alphaz[4]-1.161895003862225*alphaz[3]-0.6708203932499369*alphaz[2]+0.8660254037844386*alphaz[1]+0.5*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(1.118033988749895*alphaz[4]+0.8660254037844386*alphaz[1]+0.5*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(1.5*alphaz[6]+1.118033988749895*alphaz[4]+1.161895003862225*alphaz[3]+0.6708203932499369*alphaz[2]+0.8660254037844386*alphaz[1]+0.5*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphav[8]; 
  alphav[0] = -(0.8660254037844386*(2.23606797749979*BstarZ_by_Bmag[1]*hamil[4]+BstarZ_by_Bmag[0]*hamil[1])*dfac_v*dfac_z)/m_; 
  alphav[1] = -(1.732050807568877*(0.5*(2.23606797749979*BstarZ_by_Bmag[0]*hamil[4]+BstarZ_by_Bmag[1]*hamil[1])+BstarZ_by_Bmag[4]*hamil[4])*dfac_v*dfac_z)/m_; 
  alphav[4] = (1.732050807568877*((-1.0*BstarZ_by_Bmag[1]*hamil[4])-0.5*hamil[1]*BstarZ_by_Bmag[4])*dfac_v*dfac_z)/m_; 
#if cflType == SURFAVG 
  // evaluate surface-averaged alpha on left 
  alphaL = 0.25*alphav[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate surface-averaged alpha on right 
  alphaR = 0.25*alphav[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // evaluate alpha at left surface quadrature points 
  alphaL = 0.25*(0.4472135954999579*alphav[4]-0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.5*alphav[0]-0.5590169943749475*alphav[4]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.25*(0.4472135954999579*alphav[4]+0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // evaluate alpha at right surface quadrature points 
  alphaR = 0.25*(0.4472135954999579*alphav[4]-0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.5*alphav[0]-0.5590169943749475*alphav[4]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.25*(0.4472135954999579*alphav[4]+0.6708203932499369*alphav[1]+0.5*alphav[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.8660254037844386*(alphaz[6]*f[6]+alphaz[4]*f[4]+alphaz[3]*f[3]+alphaz[2]*f[2]+alphaz[1]*f[1]+alphaz[0]*f[0]); 
  out[2] += 0.8660254037844386*(alphav[4]*f[4]+alphav[1]*f[1]+alphav[0]*f[0]); 
  out[3] += 0.1*(7.745966692414834*alphaz[3]*f[7]+8.660254037844387*(alphaz[4]*f[6]+f[4]*alphaz[6])+7.745966692414834*(alphaz[2]*f[5]+alphav[1]*f[4]+f[1]*alphav[4])+8.660254037844386*(alphaz[1]*f[3]+f[1]*alphaz[3]+alphaz[0]*f[2]+f[0]*alphaz[2]+alphav[0]*f[1]+f[0]*alphav[1])); 
  out[4] += 0.1*(17.32050807568877*(alphaz[3]*f[6]+f[3]*alphaz[6])+17.32050807568877*(alphaz[1]*f[4]+f[1]*alphaz[4])+19.36491673103709*(alphaz[2]*f[3]+f[2]*alphaz[3]+alphaz[0]*f[1]+f[0]*alphaz[1])); 
  out[5] += 1.936491673103709*(alphav[4]*f[6]+alphav[1]*f[3]+alphav[0]*f[2]); 
  out[6] += 0.01428571428571429*(108.4435336938077*alphaz[6]*f[7]+121.2435565298214*(alphaz[2]*f[7]+alphaz[1]*f[6]+f[1]*alphaz[6])+121.2435565298214*alphaz[3]*f[5]+(38.72983346207417*alphav[4]+121.2435565298214*alphaz[3]+60.62177826491071*alphav[0])*f[4]+121.2435565298214*f[3]*alphaz[4]+60.62177826491071*f[0]*alphav[4]+135.5544171172596*(alphaz[0]*f[3]+f[0]*alphaz[3]+alphaz[1]*f[2])+f[1]*(135.5544171172596*alphaz[2]+54.22176684690384*alphav[1])); 
  out[7] += 0.1*(8.660254037844386*alphaz[1]*f[7]+(7.745966692414834*alphaz[6]+17.32050807568877*alphav[1])*f[6]+8.660254037844387*alphaz[0]*f[5]+f[3]*(17.32050807568877*alphav[4]+7.745966692414834*alphaz[3]+19.36491673103708*alphav[0])+(7.745966692414834*alphaz[2]+19.36491673103708*alphav[1])*f[2]); 
  return cflFreq; 
} 
