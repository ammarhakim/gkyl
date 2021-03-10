#include <GyrokineticModDecl.h> 
double EmGyrokineticSimpleHelicalVol1x2vSerP2_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential .
  // f: Distribution function.
  // out: output increment.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wvpar = w[1];
  double rdvpar2 = 2.0/dxv[1];
  double wmu = w[2];
  double rdmu2 = 2.0/dxv[2];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[2]*w[2];
  double rdmu2Sq = rdmu2*rdmu2;

  double hamil[20]; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+2.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*phi[1]*q_; 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[7] = 2.0*phi[2]*q_; 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 

  double BstarZdBmag[20]; 
  BstarZdBmag[0] = 2.0*cmag[0]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[20]; 
  alphax[0] = (0.6123724356957944*BstarZdBmag[0]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[2] = (1.369306393762915*BstarZdBmag[0]*hamil[8]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.1767766952966368*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.1767766952966368*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alphax[0]-0.4743416490252568*alphax[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphax[2]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[20]; 
  alphavpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamil[1]*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(1.369306393762915*BstarZdBmag[0]*hamil[7]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.1767766952966368*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.1767766952966368*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.3535533905932737*alphavpar[0]-0.4743416490252568*alphavpar[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphavpar[0]-0.4743416490252568*alphavpar[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alphavpar[0]-0.4743416490252568*alphavpar[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.3535533905932737*alphavpar[0]-0.4743416490252568*alphavpar[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alphavpar[0]-0.4743416490252568*alphavpar[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alphavpar[0]-0.4743416490252568*alphavpar[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*(alphax[2]*f[2]+alphax[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[4] += 0.1224744871391589*(4.47213595499958*(alphax[2]*f[8]+alphavpar[1]*f[7])+5.0*(alphax[0]*f[2]+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1])); 
  out[5] += 0.6123724356957944*(alphax[2]*f[6]+alphax[0]*f[3]); 
  out[6] += 0.6123724356957944*(alphavpar[1]*f[5]+alphavpar[0]*f[3]); 
  out[7] += 1.369306393762915*(alphax[2]*f[4]+alphax[0]*f[1]); 
  out[8] += 1.369306393762915*(alphavpar[1]*f[4]+alphavpar[0]*f[2]); 
  out[10] += 0.07071067811865474*(7.745966692414834*(alphax[2]*f[14]+alphavpar[1]*f[13])+8.660254037844386*(alphax[0]*f[6]+alphavpar[0]*f[5]+(alphax[2]+alphavpar[1])*f[3])); 
  out[11] += 0.07071067811865474*(17.32050807568877*alphax[2]*f[12]+8.660254037844387*alphavpar[0]*f[7]+19.36491673103708*alphax[0]*f[4]+f[1]*(19.36491673103708*alphax[2]+7.745966692414834*alphavpar[1])); 
  out[12] += 0.07071067811865474*(17.32050807568877*alphavpar[1]*f[11]+8.660254037844387*alphax[0]*f[8]+19.36491673103708*alphavpar[0]*f[4]+(7.745966692414834*alphax[2]+19.36491673103708*alphavpar[1])*f[2]); 
  out[13] += 1.369306393762915*(alphax[2]*f[10]+alphax[0]*f[5]); 
  out[14] += 1.369306393762915*(alphavpar[1]*f[10]+alphavpar[0]*f[6]); 
  out[15] += 0.07071067811865474*(8.660254037844386*alphax[2]*f[16]+8.660254037844387*alphax[0]*f[9]); 
  out[16] += 0.07071067811865474*(8.660254037844386*alphavpar[1]*f[15]+8.660254037844387*alphavpar[0]*f[9]); 
  out[17] += 0.07071067811865474*(17.32050807568877*alphax[2]*f[18]+8.660254037844387*alphavpar[0]*f[13]+19.36491673103709*alphax[0]*f[10]+(19.36491673103709*alphax[2]+7.745966692414834*alphavpar[1])*f[5]); 
  out[18] += 0.07071067811865474*(17.32050807568877*alphavpar[1]*f[17]+8.660254037844387*alphax[0]*f[14]+19.36491673103709*alphavpar[0]*f[10]+(7.745966692414834*alphax[2]+19.36491673103709*alphavpar[1])*f[6]); 
  out[19] += 0.07071067811865474*(8.660254037844387*(alphax[0]*f[16]+alphavpar[0]*f[15])+8.660254037844386*(alphax[2]+alphavpar[1])*f[9]); 
  return cflFreq; 
} 
double EmGyrokineticSimpleHelicalVol1x2vSerP2_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential .
  // f: Distribution function.
  // out: output increment.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wvpar = w[1];
  double rdvpar2 = 2.0/dxv[1];
  double wmu = w[2];
  double rdmu2 = 2.0/dxv[2];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wvparSq = w[1]*w[1];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[2]*w[2];
  double rdmu2Sq = rdmu2*rdmu2;

  double hamil[20]; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+2.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[7] = 2.0*(bmag[2]*wmu+phi[2]*q_); 
  hamil[8] = (0.421637021355784*m_)/rdvpar2Sq; 
  hamil[13] = (1.154700538379251*bmag[2])/rdmu2; 

  double BstarZdBmag[20]; 
  BstarZdBmag[0] = 2.0*cmag[0]; 
  BstarZdBmag[1] = 2.0*cmag[1]; 
  BstarZdBmag[7] = 2.0*cmag[2]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[20]; 
  alphax[0] = (0.6123724356957944*BstarZdBmag[0]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.6123724356957944*BstarZdBmag[1]*hamil[2]*rdvpar2*rdx2)/m_; 
  alphax[2] = (1.369306393762915*BstarZdBmag[0]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphax[4] = (1.369306393762915*BstarZdBmag[1]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphax[7] = (0.6123724356957944*hamil[2]*BstarZdBmag[7]*rdvpar2*rdx2)/m_; 
  alphax[11] = (1.369306393762915*BstarZdBmag[7]*hamil[8]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*(3.16227766016838*alphax[7]-2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(3.16227766016838*alphax[7]+2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*((-1.060660171779821*alphax[11])+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.7905694150420947*alphax[7]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.060660171779821*alphax[11]+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-1.060660171779821*alphax[11])+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.060660171779821*alphax[11]+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-1.060660171779821*alphax[11])+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.7905694150420947*alphax[7]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(1.060660171779821*alphax[11]+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*((-1.060660171779821*alphax[11])+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.7905694150420947*alphax[7]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(1.060660171779821*alphax[11]+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-1.060660171779821*alphax[11])+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(1.060660171779821*alphax[11]+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-1.060660171779821*alphax[11])+0.7905694150420947*alphax[7]-0.8215838362577489*alphax[4]-0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.7905694150420947*alphax[7]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(1.060660171779821*alphax[11]+0.7905694150420947*alphax[7]+0.8215838362577489*alphax[4]+0.4743416490252568*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[20]; 
  alphavpar[0] = -(0.6123724356957944*(2.23606797749979*BstarZdBmag[1]*hamil[7]+BstarZdBmag[0]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[1] = -(0.6123724356957944*((2.0*BstarZdBmag[7]+2.23606797749979*BstarZdBmag[0])*hamil[7]+BstarZdBmag[1]*hamil[1])*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(0.3535533905932737*(3.872983346207417*BstarZdBmag[1]*hamil[13]+1.732050807568877*BstarZdBmag[0]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphavpar[5] = -(0.3535533905932737*(3.872983346207417*(0.8944271909999159*BstarZdBmag[7]+BstarZdBmag[0])*hamil[13]+1.732050807568877*BstarZdBmag[1]*hamil[5])*rdvpar2*rdx2)/m_; 
  alphavpar[7] = -(0.6123724356957944*(2.0*BstarZdBmag[1]*hamil[7]+hamil[1]*BstarZdBmag[7])*rdvpar2*rdx2)/m_; 
  alphavpar[13] = -(0.3535533905932737*(3.464101615137754*BstarZdBmag[1]*hamil[13]+1.732050807568877*hamil[5]*BstarZdBmag[7])*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.1767766952966368*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.1767766952966368*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*((-0.4242640687119285*alphavpar[13])+0.3162277660168379*alphavpar[7]+0.6363961030678926*alphavpar[5]-0.4743416490252568*(alphavpar[3]+alphavpar[1])+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.5303300858899104*alphavpar[13]-0.3952847075210473*alphavpar[7]-0.4743416490252568*alphavpar[3]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.4242640687119285*alphavpar[13])+0.3162277660168379*alphavpar[7]-0.6363961030678926*alphavpar[5]-0.4743416490252568*alphavpar[3]+0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3162277660168379*alphavpar[7]-0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3162277660168379*alphavpar[7]+0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4242640687119285*alphavpar[13]+0.3162277660168379*alphavpar[7]-0.6363961030678926*alphavpar[5]+0.4743416490252568*alphavpar[3]-0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.5303300858899104*alphavpar[13])-0.3952847075210473*alphavpar[7]+0.4743416490252568*alphavpar[3]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4242640687119285*alphavpar[13]+0.3162277660168379*alphavpar[7]+0.6363961030678926*alphavpar[5]+0.4743416490252568*(alphavpar[3]+alphavpar[1])+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*((-0.4242640687119285*alphavpar[13])+0.3162277660168379*alphavpar[7]+0.6363961030678926*alphavpar[5]-0.4743416490252568*(alphavpar[3]+alphavpar[1])+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.5303300858899104*alphavpar[13]-0.3952847075210473*alphavpar[7]-0.4743416490252568*alphavpar[3]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.4242640687119285*alphavpar[13])+0.3162277660168379*alphavpar[7]-0.6363961030678926*alphavpar[5]-0.4743416490252568*alphavpar[3]+0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3162277660168379*alphavpar[7]-0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3162277660168379*alphavpar[7]+0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4242640687119285*alphavpar[13]+0.3162277660168379*alphavpar[7]-0.6363961030678926*alphavpar[5]+0.4743416490252568*alphavpar[3]-0.4743416490252568*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.5303300858899104*alphavpar[13])-0.3952847075210473*alphavpar[7]+0.4743416490252568*alphavpar[3]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4242640687119285*alphavpar[13]+0.3162277660168379*alphavpar[7]+0.6363961030678926*alphavpar[5]+0.4743416490252568*(alphavpar[3]+alphavpar[1])+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*(alphax[11]*f[11]+alphax[7]*f[7]+alphax[4]*f[4]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphavpar[13]*f[13]+alphavpar[7]*f[7]+alphavpar[5]*f[5]+alphavpar[3]*f[3]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[4] += 0.07071067811865474*(7.745966692414834*(alphavpar[5]*f[13]+f[5]*alphavpar[13]+alphax[4]*f[12])+8.660254037844387*(alphax[7]*f[11]+f[7]*alphax[11])+7.745966692414834*(alphax[2]*f[8]+alphavpar[1]*f[7]+f[1]*alphavpar[7])+8.660254037844386*(alphavpar[3]*f[5]+f[3]*alphavpar[5]+alphax[1]*f[4]+f[1]*alphax[4]+alphax[0]*f[2]+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1])); 
  out[5] += 0.07071067811865474*(8.660254037844387*(alphax[11]*f[17]+alphax[7]*f[13])+8.660254037844386*(alphax[4]*f[10]+alphax[2]*f[6]+alphax[1]*f[5]+alphax[0]*f[3])); 
  out[6] += 0.07071067811865474*(7.745966692414834*alphavpar[5]*f[15]+8.660254037844387*(alphavpar[7]*f[13]+f[7]*alphavpar[13])+7.745966692414834*alphavpar[3]*f[9]+8.660254037844386*(alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphavpar[0]*f[3]+f[0]*alphavpar[3])); 
  out[7] += 0.07071067811865474*(17.32050807568877*(alphax[4]*f[11]+f[4]*alphax[11])+17.32050807568877*(alphax[1]*f[7]+f[1]*alphax[7])+19.36491673103709*(alphax[2]*f[4]+f[2]*alphax[4]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[8] += 1.369306393762915*(alphavpar[13]*f[17]+alphavpar[7]*f[11]+alphavpar[5]*f[10]+alphavpar[3]*f[6]+alphavpar[1]*f[4]+alphavpar[0]*f[2]); 
  out[10] += 0.07071067811865474*(7.745966692414834*alphax[4]*f[18]+8.660254037844386*alphax[7]*f[17]+6.928203230275509*alphavpar[13]*f[15]+7.745966692414834*(alphavpar[3]*f[15]+alphax[2]*f[14])+8.660254037844386*alphax[11]*f[13]+7.745966692414834*(alphavpar[1]*f[13]+f[1]*alphavpar[13])+8.660254037844386*alphax[1]*f[10]+7.745966692414834*(alphavpar[5]*(f[9]+f[7])+f[5]*alphavpar[7])+8.660254037844386*(alphax[0]*f[6]+(alphax[4]+alphavpar[0])*f[5]+f[0]*alphavpar[5]+(alphax[2]+alphavpar[1])*f[3]+f[1]*alphavpar[3])); 
  out[11] += 0.01010152544552211*(38.72983346207417*alphavpar[13]*f[13]+60.6217782649107*(alphavpar[3]*f[13]+f[3]*alphavpar[13])+108.4435336938077*alphax[11]*f[12]+121.2435565298214*(alphax[2]*f[12]+alphax[1]*f[11]+f[1]*alphax[11])+121.2435565298214*alphax[4]*f[8]+(38.72983346207417*alphavpar[7]+121.2435565298214*alphax[4]+60.62177826491071*alphavpar[0])*f[7]+121.2435565298214*f[4]*alphax[7]+60.62177826491071*f[0]*alphavpar[7]+54.22176684690384*alphavpar[5]*f[5]+135.5544171172596*(alphax[0]*f[4]+f[0]*alphax[4]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphavpar[1])); 
  out[12] += 0.07071067811865474*(17.32050807568877*alphavpar[5]*f[17]+17.32050807568877*f[10]*alphavpar[13]+8.660254037844386*alphax[1]*f[12]+(7.745966692414834*alphax[11]+17.32050807568877*alphavpar[1])*f[11]+19.36491673103708*alphavpar[3]*f[10]+8.660254037844387*alphax[0]*f[8]+17.32050807568877*f[4]*alphavpar[7]+19.36491673103708*alphavpar[5]*f[6]+(7.745966692414834*alphax[4]+19.36491673103708*alphavpar[0])*f[4]+(7.745966692414834*alphax[2]+19.36491673103708*alphavpar[1])*f[2]); 
  out[13] += 0.07071067811865474*(17.32050807568877*alphax[4]*f[17]+17.32050807568877*alphax[1]*f[13]+f[10]*(17.32050807568877*alphax[11]+19.36491673103708*alphax[2])+17.32050807568877*f[5]*alphax[7]+19.36491673103708*(alphax[4]*f[6]+alphax[0]*f[5]+alphax[1]*f[3])); 
  out[14] += 0.07071067811865474*(17.32050807568877*alphavpar[5]*f[19]+19.36491673103708*alphavpar[7]*f[17]+17.32050807568877*alphavpar[3]*f[16]+19.36491673103708*(f[11]*alphavpar[13]+alphavpar[1]*f[10]+alphavpar[0]*f[6]+f[4]*alphavpar[5]+f[2]*alphavpar[3])); 
  out[15] += 0.07071067811865474*(8.660254037844387*alphax[4]*f[19]+8.660254037844386*(alphax[2]*f[16]+alphax[1]*f[15])+8.660254037844387*alphax[0]*f[9]); 
  out[16] += 0.07071067811865474*(8.660254037844386*alphavpar[1]*f[15]+7.745966692414834*alphavpar[13]*f[13]+8.660254037844387*alphavpar[0]*f[9]+7.745966692414834*(alphavpar[5]*f[5]+alphavpar[3]*f[3])); 
  out[17] += 0.002020305089104421*(542.2176684690384*alphax[11]*f[18]+606.217782649107*(alphax[2]*f[18]+alphax[1]*f[17])+242.4871130596428*alphavpar[5]*f[15]+606.2177826491072*alphax[4]*f[14]+(193.6491673103708*alphavpar[7]+606.2177826491072*alphax[4]+303.1088913245536*alphavpar[0])*f[13]+(271.1088342345192*f[9]+193.6491673103708*f[7]+303.1088913245536*f[0])*alphavpar[13]+606.2177826491072*f[5]*alphax[11]+(606.217782649107*alphax[7]+677.7720855862981*alphax[0])*f[10]+303.1088913245535*(alphavpar[3]*f[7]+f[3]*alphavpar[7])+677.7720855862981*(alphax[1]*f[6]+alphax[2]*f[5])+271.1088342345192*(alphavpar[1]*f[5]+f[1]*alphavpar[5])+677.7720855862981*f[3]*alphax[4]); 
  out[18] += 0.07071067811865474*((15.49193338482967*alphavpar[13]+17.32050807568877*alphavpar[3])*f[19]+8.660254037844386*alphax[1]*f[18]+(7.745966692414834*alphax[11]+17.32050807568877*alphavpar[1])*f[17]+17.32050807568877*alphavpar[5]*f[16]+8.660254037844387*alphax[0]*f[14]+17.32050807568877*(f[4]*alphavpar[13]+alphavpar[5]*f[11])+(17.32050807568877*alphavpar[7]+7.745966692414834*alphax[4]+19.36491673103709*alphavpar[0])*f[10]+7.745966692414834*alphax[2]*f[6]+19.36491673103709*(alphavpar[1]*f[6]+f[2]*alphavpar[5]+alphavpar[3]*f[4])); 
  out[19] += 0.01414213562373095*(43.30127018922193*alphax[1]*f[19]+43.30127018922195*alphax[0]*f[16]+(38.72983346207417*alphavpar[7]+43.30127018922195*(alphax[4]+alphavpar[0]))*f[15]+34.64101615137755*(alphavpar[5]*f[13]+f[5]*alphavpar[13])+43.30127018922193*(alphax[2]+alphavpar[1])*f[9]+38.72983346207418*(alphavpar[3]*f[5]+f[3]*alphavpar[5])); 
  return cflFreq; 
} 
double EmGyrokineticSimpleHelicalStep2Vol1x2vSerP2(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // dApardt: time derivative of Apar.
  // f: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[1]; 
  out[2] += -(1.224744871391589*(dApardt[2]*f[7]+dApardt[1]*f[1]+dApardt[0]*f[0])*q_*rdvpar2)/m_; 
  out[4] += -(0.2449489742783178*(4.47213595499958*dApardt[1]*f[7]+4.47213595499958*f[1]*dApardt[2]+5.0*dApardt[0]*f[1]+5.0*f[0]*dApardt[1])*q_*rdvpar2)/m_; 
  out[6] += -(0.1414213562373095*(8.660254037844387*dApardt[2]*f[13]+8.660254037844386*dApardt[1]*f[5]+8.660254037844386*dApardt[0]*f[3])*q_*rdvpar2)/m_; 
  out[8] += -(0.7071067811865475*(3.872983346207417*dApardt[2]*f[11]+3.872983346207417*dApardt[1]*f[4]+3.872983346207417*dApardt[0]*f[2])*q_*rdvpar2)/m_; 
  out[10] += -(0.1414213562373095*(7.745966692414834*dApardt[1]*f[13]+7.745966692414834*dApardt[2]*f[5]+8.660254037844386*dApardt[0]*f[5]+8.660254037844386*dApardt[1]*f[3])*q_*rdvpar2)/m_; 
  out[11] += -(0.07824607964359516*(10.0*dApardt[2]*f[7]+15.65247584249853*dApardt[0]*f[7]+15.65247584249853*f[0]*dApardt[2]+14.0*dApardt[1]*f[1])*q_*rdvpar2)/m_; 
  out[12] += -(0.1414213562373095*(17.32050807568877*dApardt[1]*f[11]+17.32050807568877*dApardt[2]*f[4]+19.36491673103708*dApardt[0]*f[4]+19.36491673103708*dApardt[1]*f[2])*q_*rdvpar2)/m_; 
  out[14] += -(2.738612787525831*(dApardt[2]*f[17]+dApardt[1]*f[10]+dApardt[0]*f[6])*q_*rdvpar2)/m_; 
  out[16] += -(0.1414213562373095*(8.660254037844386*dApardt[1]*f[15]+8.660254037844387*dApardt[0]*f[9])*q_*rdvpar2)/m_; 
  out[17] += -(0.02020305089104421*(38.72983346207417*dApardt[2]*f[13]+60.62177826491071*dApardt[0]*f[13]+54.22176684690384*dApardt[1]*f[5]+60.6217782649107*dApardt[2]*f[3])*q_*rdvpar2)/m_; 
  out[18] += -(1.224744871391589*(2.0*dApardt[1]*f[17]+2.0*dApardt[2]*f[10]+2.23606797749979*dApardt[0]*f[10]+2.23606797749979*dApardt[1]*f[6])*q_*rdvpar2)/m_; 
  out[19] += -(0.1414213562373095*(7.745966692414834*dApardt[2]*f[15]+8.660254037844387*dApardt[0]*f[15]+8.660254037844386*dApardt[1]*f[9])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -(0.3535533905932737*dApardt[0]*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = -(0.3535533905932737*dApardt[0]*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.125*(0.6324555320336759*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336759*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336759*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336759*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336759*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.125*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.125*(0.6324555320336759*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.125*(0.6324555320336759*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336759*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336759*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336759*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336759*dApardt[2]-0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.125*(0.7071067811865475*dApardt[0]-0.7905694150420947*dApardt[2])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.125*(0.6324555320336759*dApardt[2]+0.9486832980505137*dApardt[1]+0.7071067811865475*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

return cflFreq; 
} 
