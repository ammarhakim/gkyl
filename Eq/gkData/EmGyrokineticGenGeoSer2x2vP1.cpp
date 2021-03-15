#include <GyrokineticModDecl.h> 
double EmGyrokineticGenGeoVol2x2vSerP1_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential .
  // f: Distribution function.
  // out: output increment.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wvpar = w[2];
  double rdvpar2 = 2.0/dxv[2];
  double wmu = w[3];
  double rdmu2 = 2.0/dxv[3];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wvparSq = w[2]*w[2];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[3]*w[3];
  double rdmu2Sq = rdmu2*rdmu2;

  double hamil[16]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 

  double BstarXdBmag[16]; 
  BstarXdBmag[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2; 
  BstarXdBmag[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2; 

  double BstarYdBmag[16]; 
  BstarYdBmag[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2; 
  BstarYdBmag[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[0]*hamil[3]*rdvpar2)/m_-(0.25*b_z[0]*jacobTotInv[0]*hamil[2]*rdy2)/q_); 
  alphax[1] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[1]*hamil[3]*rdvpar2)/m_-(0.25*b_z[0]*jacobTotInv[0]*hamil[5]*rdy2)/q_); 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[16]; 
  alphay[0] = 0.4330127018922193*((b_z[0]*jacobTotInv[0]*hamil[1]*rdx2)/q_+(BstarYdBmag[0]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[2] = 0.4330127018922193*((b_z[0]*jacobTotInv[0]*hamil[5]*rdx2)/q_+(BstarYdBmag[2]*hamil[3]*rdvpar2)/m_)*rdy2; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[16]; 
  alphavpar[0] = -(0.4330127018922193*rdvpar2*(BstarYdBmag[0]*hamil[2]*rdy2+BstarXdBmag[0]*hamil[1]*rdx2))/m_; 
  alphavpar[1] = -(0.4330127018922193*rdvpar2*(BstarYdBmag[0]*hamil[5]*rdy2+BstarXdBmag[1]*hamil[1]*rdx2))/m_; 
  alphavpar[2] = -(0.4330127018922193*rdvpar2*(BstarYdBmag[2]*hamil[2]*rdy2+BstarXdBmag[0]*hamil[5]*rdx2))/m_; 
  alphavpar[5] = -(0.4330127018922193*hamil[5]*rdvpar2*(BstarYdBmag[2]*rdy2+BstarXdBmag[1]*rdx2))/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphavpar[5]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphavpar[1]+alphavpar[0])-0.25*(alphavpar[5]+alphavpar[2])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[5])+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.015625*(alphavpar[5]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[5]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphavpar[1]+alphavpar[0])-0.25*(alphavpar[5]+alphavpar[2])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[5])+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.015625*(alphavpar[5]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*alphavpar[5]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphavpar[1]+alphavpar[0])-0.25*(alphavpar[5]+alphavpar[2])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[5])+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.015625*(alphavpar[5]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[5]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphavpar[1]+alphavpar[0])-0.25*(alphavpar[5]+alphavpar[2])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[5])+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.015625*(alphavpar[5]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[2]*f[2]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphavpar[5]*f[5]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[5] += 0.4330127018922193*((alphay[2]+alphax[1])*f[5]+alphax[0]*f[2]+alphay[0]*f[1]); 
  out[6] += 0.4330127018922193*(alphax[1]*f[6]+alphavpar[2]*f[5]+f[2]*alphavpar[5]+alphax[0]*f[3]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[7] += 0.4330127018922193*(alphay[2]*f[7]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphay[0]*f[3]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[8] += 0.4330127018922193*(alphax[1]*f[8]+alphax[0]*f[4]); 
  out[9] += 0.4330127018922193*(alphay[2]*f[9]+alphay[0]*f[4]); 
  out[10] += 0.4330127018922193*(alphavpar[5]*f[12]+alphavpar[2]*f[9]+alphavpar[1]*f[8]+alphavpar[0]*f[4]); 
  out[11] += 0.4330127018922193*((alphay[2]+alphax[1])*f[11]+alphax[0]*f[7]+alphay[0]*f[6]+alphavpar[0]*f[5]+f[0]*alphavpar[5]+alphavpar[1]*f[2]+f[1]*alphavpar[2]); 
  out[12] += 0.4330127018922193*((alphay[2]+alphax[1])*f[12]+alphax[0]*f[9]+alphay[0]*f[8]); 
  out[13] += 0.4330127018922193*(alphax[1]*f[13]+alphavpar[2]*f[12]+alphax[0]*f[10]+alphavpar[5]*f[9]+alphavpar[0]*f[8]+alphavpar[1]*f[4]); 
  out[14] += 0.4330127018922193*(alphay[2]*f[14]+alphavpar[1]*f[12]+alphay[0]*f[10]+alphavpar[0]*f[9]+alphavpar[5]*f[8]+alphavpar[2]*f[4]); 
  out[15] += 0.4330127018922193*((alphay[2]+alphax[1])*f[15]+alphax[0]*f[14]+alphay[0]*f[13]+alphavpar[0]*f[12]+alphavpar[1]*f[9]+alphavpar[2]*f[8]+f[4]*alphavpar[5]); 
  return cflFreq; 
} 
double EmGyrokineticGenGeoVol2x2vSerP1_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential .
  // f: Distribution function.
  // out: output increment.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wvpar = w[2];
  double rdvpar2 = 2.0/dxv[2];
  double wmu = w[3];
  double rdmu2 = 2.0/dxv[3];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wvparSq = w[2]*w[2];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[3]*w[3];
  double rdmu2Sq = rdmu2*rdmu2;

  double hamil[16]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 

  double BstarXdBmag[16]; 
  BstarXdBmag[0] = 0.8660254037844386*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])*rdy2; 
  BstarXdBmag[1] = 0.1732050807568877*((9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])*rdy2; 

  double BstarYdBmag[16]; 
  BstarYdBmag[0] = -(0.8660254037844386*rdx2*(2.0*jacobTotInv[0]*b_z[1]*m_*wvpar+(2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*q_))/q_; 
  BstarYdBmag[1] = -(0.8660254037844386*rdx2*(2.0*b_z[1]*jacobTotInv[1]*m_*wvpar+((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])*q_))/q_; 
  BstarYdBmag[2] = -0.8660254037844386*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])*rdx2; 
  BstarYdBmag[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[5] = -0.8660254037844386*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])*rdx2; 
  BstarYdBmag[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[0]*hamil[3]*rdvpar2)/m_-(0.25*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[5]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[2])*rdy2)/q_); 
  alphax[1] = 1.732050807568877*rdx2*((((-0.25*(b_z[0]*jacobTotInv[0]*hamil[5]+(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[2]))-0.45*b_z[1]*jacobTotInv[1]*hamil[5])*rdy2)/q_+(0.25*BstarXdBmag[1]*hamil[3]*rdvpar2)/m_); 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[16]; 
  alphay[0] = 0.4330127018922193*((hamil[1]*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*rdx2)/q_+(BstarYdBmag[0]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[1] = 0.4330127018922193*((hamil[1]*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*rdx2)/q_+(BstarYdBmag[1]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[2] = 0.4330127018922193*(((b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[5]*rdx2)/q_+(BstarYdBmag[2]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[3] = (0.4330127018922193*BstarYdBmag[3]*hamil[3]*rdvpar2*rdy2)/m_; 
  alphay[4] = (0.4330127018922193*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[8]*rdx2*rdy2)/q_; 
  alphay[5] = 0.4330127018922193*(((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[5]*rdx2)/q_+(hamil[3]*BstarYdBmag[5]*rdvpar2)/m_)*rdy2; 
  alphay[6] = (0.4330127018922193*hamil[3]*BstarYdBmag[6]*rdvpar2*rdy2)/m_; 
  alphay[8] = (0.4330127018922193*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[8]*rdx2*rdy2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[16]; 
  alphavpar[0] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[1]*hamil[5]+BstarYdBmag[0]*hamil[2])*rdy2+BstarXdBmag[0]*hamil[1]*rdx2))/m_; 
  alphavpar[1] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[0]*hamil[5]+BstarYdBmag[1]*hamil[2])*rdy2+BstarXdBmag[1]*hamil[1]*rdx2))/m_; 
  alphavpar[2] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[5]*hamil[5]+BstarYdBmag[2]*hamil[2])*rdy2+BstarXdBmag[0]*hamil[5]*rdx2))/m_; 
  alphavpar[3] = -(0.4330127018922193*(hamil[5]*BstarYdBmag[6]+hamil[2]*BstarYdBmag[3])*rdvpar2*rdy2)/m_; 
  alphavpar[4] = -(0.4330127018922193*BstarXdBmag[0]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphavpar[5] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[2]*hamil[5]+hamil[2]*BstarYdBmag[5])*rdy2+BstarXdBmag[1]*hamil[5]*rdx2))/m_; 
  alphavpar[6] = -(0.4330127018922193*(hamil[2]*BstarYdBmag[6]+BstarYdBmag[3]*hamil[5])*rdvpar2*rdy2)/m_; 
  alphavpar[8] = -(0.4330127018922193*BstarXdBmag[1]*hamil[8]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphavpar[3]-1.0*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphavpar[3]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphavpar[8]*f[8]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[5] += 0.4330127018922193*(alphay[4]*f[8]+f[4]*alphay[8]+alphay[3]*f[6]+f[3]*alphay[6]+(alphay[2]+alphax[1])*f[5]+f[2]*(alphay[5]+alphax[0])+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[6] += 0.4330127018922193*(alphavpar[4]*f[8]+f[4]*alphavpar[8]+(alphavpar[3]+alphax[1])*f[6]+f[3]*alphavpar[6]+alphavpar[2]*f[5]+f[2]*alphavpar[5]+alphax[0]*f[3]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[7] += 0.4330127018922193*(alphay[8]*f[13]+alphavpar[8]*f[12]+(alphavpar[6]+alphay[5])*f[11]+alphay[4]*f[10]+alphavpar[4]*f[9]+(alphavpar[3]+alphay[2])*f[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphay[0]*f[3]+f[0]*alphay[3]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[8] += 0.4330127018922193*(alphax[1]*f[8]+alphax[0]*f[4]); 
  out[9] += 0.4330127018922193*(alphay[6]*f[13]+alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4]); 
  out[10] += 0.4330127018922193*(alphavpar[6]*f[13]+alphavpar[5]*f[12]+alphavpar[3]*f[10]+alphavpar[2]*f[9]+alphavpar[1]*f[8]+f[1]*alphavpar[8]+alphavpar[0]*f[4]+f[0]*alphavpar[4]); 
  out[11] += 0.4330127018922193*(alphay[4]*f[13]+alphavpar[4]*f[12]+(alphavpar[3]+alphay[2]+alphax[1])*f[11]+alphay[8]*f[10]+alphavpar[8]*f[9]+(alphavpar[6]+alphay[5]+alphax[0])*f[7]+alphay[0]*f[6]+f[0]*alphay[6]+alphavpar[0]*f[5]+f[0]*alphavpar[5]+alphay[1]*f[3]+f[1]*alphay[3]+alphavpar[1]*f[2]+f[1]*alphavpar[2]); 
  out[12] += 0.4330127018922193*(alphay[3]*f[13]+(alphay[2]+alphax[1])*f[12]+alphay[6]*f[10]+(alphay[5]+alphax[0])*f[9]+alphay[0]*f[8]+f[0]*alphay[8]+alphay[1]*f[4]+f[1]*alphay[4]); 
  out[13] += 0.4330127018922193*((alphavpar[3]+alphax[1])*f[13]+alphavpar[2]*f[12]+(alphavpar[6]+alphax[0])*f[10]+alphavpar[5]*f[9]+alphavpar[0]*f[8]+f[0]*alphavpar[8]+alphavpar[1]*f[4]+f[1]*alphavpar[4]); 
  out[14] += 0.4330127018922193*((alphavpar[6]+alphay[5])*f[15]+(alphavpar[3]+alphay[2])*f[14]+alphay[1]*f[13]+alphavpar[1]*f[12]+alphay[0]*f[10]+alphavpar[0]*f[9]+(alphay[6]+alphavpar[5])*f[8]+f[6]*alphay[8]+f[5]*alphavpar[8]+(alphay[3]+alphavpar[2])*f[4]+f[3]*alphay[4]+f[2]*alphavpar[4]); 
  out[15] += 0.4330127018922193*((alphavpar[3]+alphay[2]+alphax[1])*f[15]+(alphavpar[6]+alphay[5]+alphax[0])*f[14]+alphay[0]*f[13]+alphavpar[0]*f[12]+alphay[1]*f[10]+alphavpar[1]*f[9]+(alphay[3]+alphavpar[2])*f[8]+f[3]*alphay[8]+f[2]*alphavpar[8]+alphay[4]*f[6]+f[4]*alphay[6]+alphavpar[4]*f[5]+f[4]*alphavpar[5]); 
  return cflFreq; 
} 
double EmGyrokineticGenGeoVol2x2vSerP1_Bvarsz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential .
  // f: Distribution function.
  // out: output increment.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wvpar = w[2];
  double rdvpar2 = 2.0/dxv[2];
  double wmu = w[3];
  double rdmu2 = 2.0/dxv[3];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wvparSq = w[2]*w[2];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[3]*w[3];
  double rdmu2Sq = rdmu2*rdmu2;

  double hamil[16]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 

  double BstarXdBmag[16]; 
  BstarXdBmag[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2; 
  BstarXdBmag[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2; 

  double BstarYdBmag[16]; 
  BstarYdBmag[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2; 
  BstarYdBmag[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[0]*hamil[3]*rdvpar2)/m_-(0.25*b_z[0]*jacobTotInv[0]*hamil[2]*rdy2)/q_); 
  alphax[1] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[1]*hamil[3]*rdvpar2)/m_-(0.25*b_z[0]*jacobTotInv[0]*hamil[5]*rdy2)/q_); 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[16]; 
  alphay[0] = 0.4330127018922193*((b_z[0]*jacobTotInv[0]*hamil[1]*rdx2)/q_+(BstarYdBmag[0]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[2] = 0.4330127018922193*((b_z[0]*jacobTotInv[0]*hamil[5]*rdx2)/q_+(BstarYdBmag[2]*hamil[3]*rdvpar2)/m_)*rdy2; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[0]-0.4330127018922193*alphay[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[16]; 
  alphavpar[0] = -(0.4330127018922193*rdvpar2*(BstarYdBmag[0]*hamil[2]*rdy2+BstarXdBmag[0]*hamil[1]*rdx2))/m_; 
  alphavpar[1] = -(0.4330127018922193*rdvpar2*(BstarYdBmag[0]*hamil[5]*rdy2+BstarXdBmag[1]*hamil[1]*rdx2))/m_; 
  alphavpar[2] = -(0.4330127018922193*rdvpar2*(BstarYdBmag[2]*hamil[2]*rdy2+BstarXdBmag[0]*hamil[5]*rdx2))/m_; 
  alphavpar[5] = -(0.4330127018922193*hamil[5]*rdvpar2*(BstarYdBmag[2]*rdy2+BstarXdBmag[1]*rdx2))/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphavpar[5]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphavpar[1]+alphavpar[0])-0.25*(alphavpar[5]+alphavpar[2])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[5])+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.015625*(alphavpar[5]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[5]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphavpar[1]+alphavpar[0])-0.25*(alphavpar[5]+alphavpar[2])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[5])+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.015625*(alphavpar[5]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*alphavpar[5]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphavpar[1]+alphavpar[0])-0.25*(alphavpar[5]+alphavpar[2])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[5])+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.015625*(alphavpar[5]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[5]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphavpar[1]+alphavpar[0])-0.25*(alphavpar[5]+alphavpar[2])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[5])+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.015625*(alphavpar[5]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[2]*f[2]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphavpar[5]*f[5]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[5] += 0.4330127018922193*((alphay[2]+alphax[1])*f[5]+alphax[0]*f[2]+alphay[0]*f[1]); 
  out[6] += 0.4330127018922193*(alphax[1]*f[6]+alphavpar[2]*f[5]+f[2]*alphavpar[5]+alphax[0]*f[3]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[7] += 0.4330127018922193*(alphay[2]*f[7]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphay[0]*f[3]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[8] += 0.4330127018922193*(alphax[1]*f[8]+alphax[0]*f[4]); 
  out[9] += 0.4330127018922193*(alphay[2]*f[9]+alphay[0]*f[4]); 
  out[10] += 0.4330127018922193*(alphavpar[5]*f[12]+alphavpar[2]*f[9]+alphavpar[1]*f[8]+alphavpar[0]*f[4]); 
  out[11] += 0.4330127018922193*((alphay[2]+alphax[1])*f[11]+alphax[0]*f[7]+alphay[0]*f[6]+alphavpar[0]*f[5]+f[0]*alphavpar[5]+alphavpar[1]*f[2]+f[1]*alphavpar[2]); 
  out[12] += 0.4330127018922193*((alphay[2]+alphax[1])*f[12]+alphax[0]*f[9]+alphay[0]*f[8]); 
  out[13] += 0.4330127018922193*(alphax[1]*f[13]+alphavpar[2]*f[12]+alphax[0]*f[10]+alphavpar[5]*f[9]+alphavpar[0]*f[8]+alphavpar[1]*f[4]); 
  out[14] += 0.4330127018922193*(alphay[2]*f[14]+alphavpar[1]*f[12]+alphay[0]*f[10]+alphavpar[0]*f[9]+alphavpar[5]*f[8]+alphavpar[2]*f[4]); 
  out[15] += 0.4330127018922193*((alphay[2]+alphax[1])*f[15]+alphax[0]*f[14]+alphay[0]*f[13]+alphavpar[0]*f[12]+alphavpar[1]*f[9]+alphavpar[2]*f[8]+f[4]*alphavpar[5]); 
  return cflFreq; 
} 
double EmGyrokineticGenGeoVol2x2vSerP1_Bvarsxz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // phi: electrostatic potential .
  // f: Distribution function.
  // out: output increment.

  double wx = w[0];
  double rdx2 = 2.0/dxv[0];
  double wy = w[1];
  double rdy2 = 2.0/dxv[1];
  double wvpar = w[2];
  double rdvpar2 = 2.0/dxv[2];
  double wmu = w[3];
  double rdmu2 = 2.0/dxv[3];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wvparSq = w[2]*w[2];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[3]*w[3];
  double rdmu2Sq = rdmu2*rdmu2;

  double hamil[16]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 

  double BstarXdBmag[16]; 
  BstarXdBmag[0] = 0.8660254037844386*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])*rdy2; 
  BstarXdBmag[1] = 0.1732050807568877*((9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])*rdy2; 

  double BstarYdBmag[16]; 
  BstarYdBmag[0] = -(0.8660254037844386*rdx2*(2.0*jacobTotInv[0]*b_z[1]*m_*wvpar+(2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*q_))/q_; 
  BstarYdBmag[1] = -(0.8660254037844386*rdx2*(2.0*b_z[1]*jacobTotInv[1]*m_*wvpar+((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])*q_))/q_; 
  BstarYdBmag[2] = -0.8660254037844386*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])*rdx2; 
  BstarYdBmag[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[5] = -0.8660254037844386*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])*rdx2; 
  BstarYdBmag[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[0]*hamil[3]*rdvpar2)/m_-(0.25*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[5]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[2])*rdy2)/q_); 
  alphax[1] = 1.732050807568877*rdx2*((((-0.25*(b_z[0]*jacobTotInv[0]*hamil[5]+(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[2]))-0.45*b_z[1]*jacobTotInv[1]*hamil[5])*rdy2)/q_+(0.25*BstarXdBmag[1]*hamil[3]*rdvpar2)/m_); 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphax[0]-0.4330127018922193*alphax[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[16]; 
  alphay[0] = 0.4330127018922193*((hamil[1]*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*rdx2)/q_+(BstarYdBmag[0]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[1] = 0.4330127018922193*((hamil[1]*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*rdx2)/q_+(BstarYdBmag[1]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[2] = 0.4330127018922193*(((b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[5]*rdx2)/q_+(BstarYdBmag[2]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[3] = (0.4330127018922193*BstarYdBmag[3]*hamil[3]*rdvpar2*rdy2)/m_; 
  alphay[4] = (0.4330127018922193*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[8]*rdx2*rdy2)/q_; 
  alphay[5] = 0.4330127018922193*(((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[5]*rdx2)/q_+(hamil[3]*BstarYdBmag[5]*rdvpar2)/m_)*rdy2; 
  alphay[6] = (0.4330127018922193*hamil[3]*BstarYdBmag[6]*rdvpar2*rdy2)/m_; 
  alphay[8] = (0.4330127018922193*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[8]*rdx2*rdy2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])-0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))+0.4330127018922193*alphay[5]-0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]-0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]+0.4330127018922193*alphay[5]-0.25*alphay[4]+0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphay[8])+0.25*alphay[6]-0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphay[8]-0.25*alphay[6]+0.4330127018922193*alphay[5]+0.25*alphay[4]-0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alphay[8]+alphay[6]))-0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alphay[8]+alphay[6])+0.4330127018922193*alphay[5]+0.25*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[16]; 
  alphavpar[0] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[1]*hamil[5]+BstarYdBmag[0]*hamil[2])*rdy2+BstarXdBmag[0]*hamil[1]*rdx2))/m_; 
  alphavpar[1] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[0]*hamil[5]+BstarYdBmag[1]*hamil[2])*rdy2+BstarXdBmag[1]*hamil[1]*rdx2))/m_; 
  alphavpar[2] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[5]*hamil[5]+BstarYdBmag[2]*hamil[2])*rdy2+BstarXdBmag[0]*hamil[5]*rdx2))/m_; 
  alphavpar[3] = -(0.4330127018922193*(hamil[5]*BstarYdBmag[6]+hamil[2]*BstarYdBmag[3])*rdvpar2*rdy2)/m_; 
  alphavpar[4] = -(0.4330127018922193*BstarXdBmag[0]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphavpar[5] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[2]*hamil[5]+hamil[2]*BstarYdBmag[5])*rdy2+BstarXdBmag[1]*hamil[5]*rdx2))/m_; 
  alphavpar[6] = -(0.4330127018922193*(hamil[2]*BstarYdBmag[6]+BstarYdBmag[3]*hamil[5])*rdvpar2*rdy2)/m_; 
  alphavpar[8] = -(0.4330127018922193*BstarXdBmag[1]*hamil[8]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphavpar[3]-1.0*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphavpar[3]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]-0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])+0.4330127018922193*alphavpar[6]+0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alphavpar[8])-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alphavpar[8]+0.4330127018922193*alphavpar[6]+0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphavpar[8]*f[8]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[5] += 0.4330127018922193*(alphay[4]*f[8]+f[4]*alphay[8]+alphay[3]*f[6]+f[3]*alphay[6]+(alphay[2]+alphax[1])*f[5]+f[2]*(alphay[5]+alphax[0])+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[6] += 0.4330127018922193*(alphavpar[4]*f[8]+f[4]*alphavpar[8]+(alphavpar[3]+alphax[1])*f[6]+f[3]*alphavpar[6]+alphavpar[2]*f[5]+f[2]*alphavpar[5]+alphax[0]*f[3]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[7] += 0.4330127018922193*(alphay[8]*f[13]+alphavpar[8]*f[12]+(alphavpar[6]+alphay[5])*f[11]+alphay[4]*f[10]+alphavpar[4]*f[9]+(alphavpar[3]+alphay[2])*f[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphay[0]*f[3]+f[0]*alphay[3]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[8] += 0.4330127018922193*(alphax[1]*f[8]+alphax[0]*f[4]); 
  out[9] += 0.4330127018922193*(alphay[6]*f[13]+alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4]); 
  out[10] += 0.4330127018922193*(alphavpar[6]*f[13]+alphavpar[5]*f[12]+alphavpar[3]*f[10]+alphavpar[2]*f[9]+alphavpar[1]*f[8]+f[1]*alphavpar[8]+alphavpar[0]*f[4]+f[0]*alphavpar[4]); 
  out[11] += 0.4330127018922193*(alphay[4]*f[13]+alphavpar[4]*f[12]+(alphavpar[3]+alphay[2]+alphax[1])*f[11]+alphay[8]*f[10]+alphavpar[8]*f[9]+(alphavpar[6]+alphay[5]+alphax[0])*f[7]+alphay[0]*f[6]+f[0]*alphay[6]+alphavpar[0]*f[5]+f[0]*alphavpar[5]+alphay[1]*f[3]+f[1]*alphay[3]+alphavpar[1]*f[2]+f[1]*alphavpar[2]); 
  out[12] += 0.4330127018922193*(alphay[3]*f[13]+(alphay[2]+alphax[1])*f[12]+alphay[6]*f[10]+(alphay[5]+alphax[0])*f[9]+alphay[0]*f[8]+f[0]*alphay[8]+alphay[1]*f[4]+f[1]*alphay[4]); 
  out[13] += 0.4330127018922193*((alphavpar[3]+alphax[1])*f[13]+alphavpar[2]*f[12]+(alphavpar[6]+alphax[0])*f[10]+alphavpar[5]*f[9]+alphavpar[0]*f[8]+f[0]*alphavpar[8]+alphavpar[1]*f[4]+f[1]*alphavpar[4]); 
  out[14] += 0.4330127018922193*((alphavpar[6]+alphay[5])*f[15]+(alphavpar[3]+alphay[2])*f[14]+alphay[1]*f[13]+alphavpar[1]*f[12]+alphay[0]*f[10]+alphavpar[0]*f[9]+(alphay[6]+alphavpar[5])*f[8]+f[6]*alphay[8]+f[5]*alphavpar[8]+(alphay[3]+alphavpar[2])*f[4]+f[3]*alphay[4]+f[2]*alphavpar[4]); 
  out[15] += 0.4330127018922193*((alphavpar[3]+alphay[2]+alphax[1])*f[15]+(alphavpar[6]+alphay[5]+alphax[0])*f[14]+alphay[0]*f[13]+alphavpar[0]*f[12]+alphay[1]*f[10]+alphavpar[1]*f[9]+(alphay[3]+alphavpar[2])*f[8]+f[3]*alphay[8]+f[2]*alphavpar[8]+alphay[4]*f[6]+f[4]*alphay[6]+alphavpar[4]*f[5]+f[4]*alphavpar[5]); 
  return cflFreq; 
} 
double EmGyrokineticGenGeoStep2Vol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // dApardt: time derivative of Apar.
  // f: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[2]; 
  out[3] += -(0.8660254037844386*(dApardt[3]*f[5]+dApardt[2]*f[2]+dApardt[1]*f[1]+dApardt[0]*f[0])*q_*rdvpar2)/m_; 
  out[6] += -(0.8660254037844386*(dApardt[2]*f[5]+f[2]*dApardt[3]+dApardt[0]*f[1]+f[0]*dApardt[1])*q_*rdvpar2)/m_; 
  out[7] += -(0.8660254037844386*(dApardt[1]*f[5]+f[1]*dApardt[3]+dApardt[0]*f[2]+f[0]*dApardt[2])*q_*rdvpar2)/m_; 
  out[10] += -(0.8660254037844386*(dApardt[3]*f[12]+dApardt[2]*f[9]+dApardt[1]*f[8]+dApardt[0]*f[4])*q_*rdvpar2)/m_; 
  out[11] += -(0.8660254037844386*(dApardt[0]*f[5]+f[0]*dApardt[3]+dApardt[1]*f[2]+f[1]*dApardt[2])*q_*rdvpar2)/m_; 
  out[13] += -(0.8660254037844386*(dApardt[2]*f[12]+dApardt[3]*f[9]+dApardt[0]*f[8]+dApardt[1]*f[4])*q_*rdvpar2)/m_; 
  out[14] += -(0.8660254037844386*(dApardt[1]*f[12]+dApardt[0]*f[9]+dApardt[3]*f[8]+dApardt[2]*f[4])*q_*rdvpar2)/m_; 
  out[15] += -(0.8660254037844386*(dApardt[0]*f[12]+dApardt[1]*f[9]+dApardt[2]*f[8]+dApardt[3]*f[4])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -(0.25*dApardt[0]*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = -(0.25*dApardt[0]*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*(0.5*dApardt[3]-0.5*(dApardt[2]+dApardt[1])+0.5*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*(0.5*(dApardt[1]+dApardt[0])-0.5*(dApardt[3]+dApardt[2]))*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0625*((-0.5*dApardt[3])+0.5*dApardt[2]-0.5*dApardt[1]+0.5*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

return cflFreq; 
} 
