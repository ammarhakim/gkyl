#include <PenalizedGyrokineticModDecl.h> 
double PenalizedGyrokineticGenGeoVol1x2vSerP1_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *f, double *out) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // penaltyChi: penalization factor multiplying the mirror force.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuB[8]; 
  hamilNoMuB[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2Sq; 
  hamilNoMuB[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuB[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 

  double hamilMuB[8]; 
  hamilMuB[0] = 2.0*bmag[0]*wmu; 
  hamilMuB[3] = (1.154700538379252*bmag[0])/rdmu2; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = (0.6123724356957944*BstarZdBmag[0]*hamilNoMuB[2]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.1767766952966368*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.1767766952966368*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0441941738241592*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0441941738241592*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[8]; 
  alphavpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamilNoMuB[1]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.1767766952966368*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.1767766952966368*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0441941738241592*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0441941738241592*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*alphax[0]*f[0]; 
  out[2] += 0.6123724356957944*alphavpar[0]*f[0]; 
  out[4] += 0.6123724356957944*(alphax[0]*f[2]+alphavpar[0]*f[1]); 
  out[5] += 0.6123724356957944*alphax[0]*f[3]; 
  out[6] += 0.6123724356957944*alphavpar[0]*f[3]; 
  out[7] += 0.6123724356957944*(alphax[0]*f[6]+alphavpar[0]*f[5]); 
  return cflFreq; 
} 
double PenalizedGyrokineticGenGeoVol1x2vSerP1_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *f, double *out) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // penaltyChi: penalization factor multiplying the mirror force.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuB[8]; 
  hamilNoMuB[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2Sq; 
  hamilNoMuB[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuB[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 

  double hamilMuB[8]; 
  hamilMuB[0] = 2.0*bmag[0]*wmu; 
  hamilMuB[1] = 2.0*bmag[1]*wmu; 
  hamilMuB[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamilMuB[5] = (1.154700538379252*bmag[1])/rdmu2; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmag[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2*wvpar+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = (0.6123724356957944*BstarZdBmag[0]*hamilNoMuB[2]*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.6123724356957944*BstarZdBmag[1]*hamilNoMuB[2]*rdvpar2*rdx2)/m_; 
  alphax[2] = (0.6123724356957944*BstarZdBmag[2]*hamilNoMuB[2]*rdvpar2*rdx2)/m_; 
  alphax[4] = (0.6123724356957944*hamilNoMuB[2]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(2.449489742783178*alphax[1]-1.414213562373095*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.6123724356957944*alphax[4]-0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*alphax[4])+0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*alphax[4]-0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*alphax[4])+0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*((-0.6123724356957944*alphax[4])-0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alphax[4]+0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*alphax[4])-0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alphax[4]+0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[8]; 
  alphavpar[0] = (1.732050807568877*(BstarZdBmag[0]*((-0.3535533905932737*hamilNoMuB[1])-0.125*penaltyChi[0]*hamilMuB[1])-0.125*(hamilMuB[5]*(BstarZdBmag[1]*penaltyChi[5]+BstarZdBmag[0]*penaltyChi[3])+BstarZdBmag[1]*hamilMuB[1]*penaltyChi[1]))*rdvpar2*rdx2)/m_; 
  alphavpar[1] = (1.732050807568877*(BstarZdBmag[1]*((-0.3535533905932737*hamilNoMuB[1])-0.125*penaltyChi[0]*hamilMuB[1])-0.125*(hamilMuB[5]*(BstarZdBmag[0]*penaltyChi[5]+BstarZdBmag[1]*penaltyChi[3])+BstarZdBmag[0]*hamilMuB[1]*penaltyChi[1]))*rdvpar2*rdx2)/m_; 
  alphavpar[2] = (1.732050807568877*(((-0.3535533905932737*hamilNoMuB[1])-0.125*penaltyChi[0]*hamilMuB[1])*BstarZdBmag[2]-0.125*(hamilMuB[5]*(BstarZdBmag[4]*penaltyChi[5]+BstarZdBmag[2]*penaltyChi[3])+hamilMuB[1]*penaltyChi[1]*BstarZdBmag[4]))*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(0.2165063509461096*(BstarZdBmag[1]*(hamilMuB[1]*penaltyChi[5]+penaltyChi[1]*hamilMuB[5])+BstarZdBmag[0]*(penaltyChi[0]*hamilMuB[5]+hamilMuB[1]*penaltyChi[3]))*rdvpar2*rdx2)/m_; 
  alphavpar[4] = (1.732050807568877*((-0.125*BstarZdBmag[2]*hamilMuB[5]*penaltyChi[5])+BstarZdBmag[4]*((-0.125*penaltyChi[3]*hamilMuB[5])-0.3535533905932737*hamilNoMuB[1])-0.125*hamilMuB[1]*(penaltyChi[0]*BstarZdBmag[4]+penaltyChi[1]*BstarZdBmag[2]))*rdvpar2*rdx2)/m_; 
  alphavpar[5] = -(0.2165063509461096*(BstarZdBmag[0]*(hamilMuB[1]*penaltyChi[5]+penaltyChi[1]*hamilMuB[5])+BstarZdBmag[1]*(penaltyChi[0]*hamilMuB[5]+hamilMuB[1]*penaltyChi[3]))*rdvpar2*rdx2)/m_; 
  alphavpar[6] = -(0.2165063509461096*(BstarZdBmag[4]*(hamilMuB[1]*penaltyChi[5]+penaltyChi[1]*hamilMuB[5])+BstarZdBmag[2]*(penaltyChi[0]*hamilMuB[5]+hamilMuB[1]*penaltyChi[3]))*rdvpar2*rdx2)/m_; 
  alphavpar[7] = -(0.2165063509461096*(hamilMuB[1]*BstarZdBmag[2]*penaltyChi[5]+(penaltyChi[0]*BstarZdBmag[4]+penaltyChi[1]*BstarZdBmag[2])*hamilMuB[5]+hamilMuB[1]*penaltyChi[3]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(2.449489742783178*alphavpar[2]-1.414213562373095*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.449489742783178*alphavpar[2]+1.414213562373095*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*((-0.6123724356957944*alphavpar[7])+0.6123724356957944*alphavpar[6]+0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*(alphavpar[7]+alphavpar[6])-0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*alphavpar[7]-0.6123724356957944*alphavpar[6]-0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*(alphavpar[7]+alphavpar[6]))+0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.6123724356957944*alphavpar[7]-0.6123724356957944*alphavpar[6]+0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*(alphavpar[7]+alphavpar[6]))-0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*alphavpar[7])+0.6123724356957944*alphavpar[6]-0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*(alphavpar[7]+alphavpar[6])+0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*(alphax[4]*f[4]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphavpar[6]*f[7]+f[6]*alphavpar[7]+alphavpar[3]*f[5]+f[3]*alphavpar[5]+(alphavpar[2]+alphax[1])*f[4]+f[1]*alphax[4]+f[2]*(alphavpar[4]+alphax[0])+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[5] += 0.6123724356957944*(alphax[4]*f[7]+alphax[2]*f[6]+alphax[1]*f[5]+alphax[0]*f[3]); 
  out[6] += 0.6123724356957944*(alphavpar[4]*f[7]+f[4]*alphavpar[7]+alphavpar[2]*f[6]+f[2]*alphavpar[6]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphavpar[0]*f[3]+f[0]*alphavpar[3]); 
  out[7] += 0.6123724356957944*((alphavpar[2]+alphax[1])*f[7]+f[2]*alphavpar[7]+(alphavpar[4]+alphax[0])*f[6]+f[4]*alphavpar[6]+(alphax[4]+alphavpar[0])*f[5]+f[0]*alphavpar[5]+(alphax[2]+alphavpar[1])*f[3]+f[1]*alphavpar[3]); 
  return cflFreq; 
} 
double PenalizedGyrokineticGenGeoVol1x2vSerP1_Bvarsz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *f, double *out) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // penaltyChi: penalization factor multiplying the mirror force.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuB[8]; 
  hamilNoMuB[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2Sq; 
  hamilNoMuB[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuB[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 

  double hamilMuB[8]; 
  hamilMuB[0] = 2.0*bmag[0]*wmu; 
  hamilMuB[3] = (1.154700538379252*bmag[0])/rdmu2; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = (0.6123724356957944*BstarZdBmag[0]*hamilNoMuB[2]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.1767766952966368*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.1767766952966368*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0441941738241592*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphax[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0441941738241592*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphax[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[8]; 
  alphavpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamilNoMuB[1]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.1767766952966368*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.1767766952966368*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0441941738241592*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0441941738241592*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*alphax[0]*f[0]; 
  out[2] += 0.6123724356957944*alphavpar[0]*f[0]; 
  out[4] += 0.6123724356957944*(alphax[0]*f[2]+alphavpar[0]*f[1]); 
  out[5] += 0.6123724356957944*alphax[0]*f[3]; 
  out[6] += 0.6123724356957944*alphavpar[0]*f[3]; 
  out[7] += 0.6123724356957944*(alphax[0]*f[6]+alphavpar[0]*f[5]); 
  return cflFreq; 
} 
double PenalizedGyrokineticGenGeoVol1x2vSerP1_Bvarsxz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *f, double *out) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // penaltyChi: penalization factor multiplying the mirror force.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuB[8]; 
  hamilNoMuB[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2Sq; 
  hamilNoMuB[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuB[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 

  double hamilMuB[8]; 
  hamilMuB[0] = 2.0*bmag[0]*wmu; 
  hamilMuB[1] = 2.0*bmag[1]*wmu; 
  hamilMuB[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamilMuB[5] = (1.154700538379252*bmag[1])/rdmu2; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmag[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2*wvpar+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[8]; 
  alphax[0] = (0.6123724356957944*BstarZdBmag[0]*hamilNoMuB[2]*rdvpar2*rdx2)/m_; 
  alphax[1] = (0.6123724356957944*BstarZdBmag[1]*hamilNoMuB[2]*rdvpar2*rdx2)/m_; 
  alphax[2] = (0.6123724356957944*BstarZdBmag[2]*hamilNoMuB[2]*rdvpar2*rdx2)/m_; 
  alphax[4] = (0.6123724356957944*hamilNoMuB[2]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(2.449489742783178*alphax[1]-1.414213562373095*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.6123724356957944*alphax[4]-0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*alphax[4])+0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*alphax[4]-0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*alphax[4])+0.3535533905932737*alphax[2]-0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*((-0.6123724356957944*alphax[4])-0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alphax[4]+0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*alphax[4])-0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alphax[4]+0.3535533905932737*alphax[2]+0.6123724356957944*alphax[1]+0.3535533905932737*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[8]; 
  alphavpar[0] = (1.732050807568877*(BstarZdBmag[0]*((-0.3535533905932737*hamilNoMuB[1])-0.125*penaltyChi[0]*hamilMuB[1])-0.125*(hamilMuB[5]*(BstarZdBmag[1]*penaltyChi[5]+BstarZdBmag[0]*penaltyChi[3])+BstarZdBmag[1]*hamilMuB[1]*penaltyChi[1]))*rdvpar2*rdx2)/m_; 
  alphavpar[1] = (1.732050807568877*(BstarZdBmag[1]*((-0.3535533905932737*hamilNoMuB[1])-0.125*penaltyChi[0]*hamilMuB[1])-0.125*(hamilMuB[5]*(BstarZdBmag[0]*penaltyChi[5]+BstarZdBmag[1]*penaltyChi[3])+BstarZdBmag[0]*hamilMuB[1]*penaltyChi[1]))*rdvpar2*rdx2)/m_; 
  alphavpar[2] = (1.732050807568877*(((-0.3535533905932737*hamilNoMuB[1])-0.125*penaltyChi[0]*hamilMuB[1])*BstarZdBmag[2]-0.125*(hamilMuB[5]*(BstarZdBmag[4]*penaltyChi[5]+BstarZdBmag[2]*penaltyChi[3])+hamilMuB[1]*penaltyChi[1]*BstarZdBmag[4]))*rdvpar2*rdx2)/m_; 
  alphavpar[3] = -(0.2165063509461096*(BstarZdBmag[1]*(hamilMuB[1]*penaltyChi[5]+penaltyChi[1]*hamilMuB[5])+BstarZdBmag[0]*(penaltyChi[0]*hamilMuB[5]+hamilMuB[1]*penaltyChi[3]))*rdvpar2*rdx2)/m_; 
  alphavpar[4] = (1.732050807568877*((-0.125*BstarZdBmag[2]*hamilMuB[5]*penaltyChi[5])+BstarZdBmag[4]*((-0.125*penaltyChi[3]*hamilMuB[5])-0.3535533905932737*hamilNoMuB[1])-0.125*hamilMuB[1]*(penaltyChi[0]*BstarZdBmag[4]+penaltyChi[1]*BstarZdBmag[2]))*rdvpar2*rdx2)/m_; 
  alphavpar[5] = -(0.2165063509461096*(BstarZdBmag[0]*(hamilMuB[1]*penaltyChi[5]+penaltyChi[1]*hamilMuB[5])+BstarZdBmag[1]*(penaltyChi[0]*hamilMuB[5]+hamilMuB[1]*penaltyChi[3]))*rdvpar2*rdx2)/m_; 
  alphavpar[6] = -(0.2165063509461096*(BstarZdBmag[4]*(hamilMuB[1]*penaltyChi[5]+penaltyChi[1]*hamilMuB[5])+BstarZdBmag[2]*(penaltyChi[0]*hamilMuB[5]+hamilMuB[1]*penaltyChi[3]))*rdvpar2*rdx2)/m_; 
  alphavpar[7] = -(0.2165063509461096*(hamilMuB[1]*BstarZdBmag[2]*penaltyChi[5]+(penaltyChi[0]*BstarZdBmag[4]+penaltyChi[1]*BstarZdBmag[2])*hamilMuB[5]+hamilMuB[1]*penaltyChi[3]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(2.449489742783178*alphavpar[2]-1.414213562373095*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.449489742783178*alphavpar[2]+1.414213562373095*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*((-0.6123724356957944*alphavpar[7])+0.6123724356957944*alphavpar[6]+0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*(alphavpar[7]+alphavpar[6])-0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*alphavpar[7]-0.6123724356957944*alphavpar[6]-0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*(alphavpar[7]+alphavpar[6]))+0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]-0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.6123724356957944*alphavpar[7]-0.6123724356957944*alphavpar[6]+0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*(alphavpar[7]+alphavpar[6]))-0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]-0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*alphavpar[7])+0.6123724356957944*alphavpar[6]-0.3535533905932737*alphavpar[5]-0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]-0.3535533905932737*alphavpar[1]+0.3535533905932737*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*(alphavpar[7]+alphavpar[6])+0.3535533905932737*alphavpar[5]+0.6123724356957944*alphavpar[4]+0.3535533905932737*alphavpar[3]+0.6123724356957944*alphavpar[2]+0.3535533905932737*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.6123724356957944*(alphax[4]*f[4]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.6123724356957944*(alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[4] += 0.6123724356957944*(alphavpar[6]*f[7]+f[6]*alphavpar[7]+alphavpar[3]*f[5]+f[3]*alphavpar[5]+(alphavpar[2]+alphax[1])*f[4]+f[1]*alphax[4]+f[2]*(alphavpar[4]+alphax[0])+f[0]*alphax[2]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[5] += 0.6123724356957944*(alphax[4]*f[7]+alphax[2]*f[6]+alphax[1]*f[5]+alphax[0]*f[3]); 
  out[6] += 0.6123724356957944*(alphavpar[4]*f[7]+f[4]*alphavpar[7]+alphavpar[2]*f[6]+f[2]*alphavpar[6]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphavpar[0]*f[3]+f[0]*alphavpar[3]); 
  out[7] += 0.6123724356957944*((alphavpar[2]+alphax[1])*f[7]+f[2]*alphavpar[7]+(alphavpar[4]+alphax[0])*f[6]+f[4]*alphavpar[6]+(alphax[4]+alphavpar[0])*f[5]+f[0]*alphavpar[5]+(alphax[2]+alphavpar[1])*f[3]+f[1]*alphavpar[3]); 
  return cflFreq; 
} 
