#include <GyrokineticModDecl.h> 
double EmGyrokineticGenGeoVol3x2vSerP1_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f, double *out) 
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
  double wz = w[2];
  double rdz2 = 2.0/dxv[2];
  double wvpar = w[3];
  double rdvpar2 = 2.0/dxv[3];
  double wmu = w[4];
  double rdmu2 = 2.0/dxv[4];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wzSq = w[2]*w[2];
  double rdz2Sq = rdz2*rdz2;
  double wvparSq = w[3]*w[3];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[4]*w[4];
  double rdmu2Sq = rdmu2*rdmu2;

  double hamil[32]; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(4.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+4.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = 2.0*phi[3]*q_; 
  hamil[4] = (3.265986323710906*m_*wvpar)/rdvpar2; 
  hamil[5] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[6] = 2.0*phi[4]*q_; 
  hamil[7] = 2.0*phi[5]*q_; 
  hamil[8] = 2.0*phi[6]*q_; 
  hamil[16] = 2.0*phi[7]*q_; 

  double BstarXdBmag[32]; 
  BstarXdBmag[0] = -0.4330127018922193*jacobTotInv[0]*(b_y[0]*Apar[3]*rdz2-1.0*b_z[0]*Apar[2]*rdy2); 
  BstarXdBmag[1] = -0.4330127018922193*jacobTotInv[0]*(b_y[0]*Apar[5]*rdz2-1.0*b_z[0]*Apar[4]*rdy2); 
  BstarXdBmag[2] = -0.4330127018922193*b_y[0]*jacobTotInv[0]*Apar[6]*rdz2; 
  BstarXdBmag[3] = 0.4330127018922193*b_z[0]*jacobTotInv[0]*Apar[6]*rdy2; 
  BstarXdBmag[6] = -0.4330127018922193*b_y[0]*jacobTotInv[0]*Apar[7]*rdz2; 
  BstarXdBmag[7] = 0.4330127018922193*b_z[0]*jacobTotInv[0]*Apar[7]*rdy2; 

  double BstarYdBmag[32]; 
  BstarYdBmag[0] = 0.4330127018922193*jacobTotInv[0]*(b_x[0]*Apar[3]*rdz2-1.0*b_z[0]*Apar[1]*rdx2); 
  BstarYdBmag[1] = 0.4330127018922193*b_x[0]*jacobTotInv[0]*Apar[5]*rdz2; 
  BstarYdBmag[2] = 0.4330127018922193*jacobTotInv[0]*(b_x[0]*Apar[6]*rdz2-1.0*b_z[0]*Apar[4]*rdx2); 
  BstarYdBmag[3] = -0.4330127018922193*b_z[0]*jacobTotInv[0]*Apar[5]*rdx2; 
  BstarYdBmag[6] = 0.4330127018922193*b_x[0]*jacobTotInv[0]*Apar[7]*rdz2; 
  BstarYdBmag[8] = -0.4330127018922193*b_z[0]*jacobTotInv[0]*Apar[7]*rdx2; 

  double BstarZdBmag[32]; 
  BstarZdBmag[0] = -0.25*jacobTotInv[0]*(1.732050807568877*b_x[0]*Apar[2]*rdy2-1.0*(1.732050807568877*b_y[0]*Apar[1]*rdx2+2.828427124746191*cmag[0])); 
  BstarZdBmag[1] = -0.4330127018922193*b_x[0]*jacobTotInv[0]*Apar[4]*rdy2; 
  BstarZdBmag[2] = 0.4330127018922193*b_y[0]*jacobTotInv[0]*Apar[4]*rdx2; 
  BstarZdBmag[3] = -0.4330127018922193*jacobTotInv[0]*(b_x[0]*Apar[6]*rdy2-1.0*b_y[0]*Apar[5]*rdx2); 
  BstarZdBmag[7] = -0.4330127018922193*b_x[0]*jacobTotInv[0]*Apar[7]*rdy2; 
  BstarZdBmag[8] = 0.4330127018922193*b_y[0]*jacobTotInv[0]*Apar[7]*rdx2; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[32]; 
  alphax[0] = 1.732050807568877*rdx2*((jacobTotInv[0]*(0.125*b_y[0]*hamil[3]*rdz2-0.125*b_z[0]*hamil[2]*rdy2))/q_+(0.1767766952966368*BstarXdBmag[0]*hamil[4]*rdvpar2)/m_); 
  alphax[1] = 1.732050807568877*rdx2*((jacobTotInv[0]*(0.125*b_y[0]*hamil[7]*rdz2-0.125*b_z[0]*hamil[6]*rdy2))/q_+(0.1767766952966368*BstarXdBmag[1]*hamil[4]*rdvpar2)/m_); 
  alphax[2] = 1.732050807568877*rdx2*((0.125*b_y[0]*jacobTotInv[0]*hamil[8]*rdz2)/q_+(0.1767766952966368*BstarXdBmag[2]*hamil[4]*rdvpar2)/m_); 
  alphax[3] = 1.732050807568877*rdx2*((0.1767766952966368*BstarXdBmag[3]*hamil[4]*rdvpar2)/m_-(0.125*b_z[0]*jacobTotInv[0]*hamil[8]*rdy2)/q_); 
  alphax[6] = 1.732050807568877*rdx2*((0.125*b_y[0]*jacobTotInv[0]*hamil[16]*rdz2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[6]*rdvpar2)/m_); 
  alphax[7] = 1.732050807568877*rdx2*((0.1767766952966368*hamil[4]*BstarXdBmag[7]*rdvpar2)/m_-(0.125*b_z[0]*jacobTotInv[0]*hamil[16]*rdy2)/q_); 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphax[1]-1.414213562373095*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[32]; 
  alphay[0] = 1.732050807568877*rdy2*((jacobTotInv[0]*(0.125*b_z[0]*hamil[1]*rdx2-0.125*b_x[0]*hamil[3]*rdz2))/q_+(0.1767766952966368*BstarYdBmag[0]*hamil[4]*rdvpar2)/m_); 
  alphay[1] = 1.732050807568877*rdy2*((0.1767766952966368*BstarYdBmag[1]*hamil[4]*rdvpar2)/m_-(0.125*b_x[0]*jacobTotInv[0]*hamil[7]*rdz2)/q_); 
  alphay[2] = 1.732050807568877*rdy2*((jacobTotInv[0]*(0.125*b_z[0]*hamil[6]*rdx2-0.125*b_x[0]*hamil[8]*rdz2))/q_+(0.1767766952966368*BstarYdBmag[2]*hamil[4]*rdvpar2)/m_); 
  alphay[3] = 1.732050807568877*((0.125*b_z[0]*jacobTotInv[0]*hamil[7]*rdx2)/q_+(0.1767766952966368*BstarYdBmag[3]*hamil[4]*rdvpar2)/m_)*rdy2; 
  alphay[6] = 1.732050807568877*rdy2*((0.1767766952966368*hamil[4]*BstarYdBmag[6]*rdvpar2)/m_-(0.125*b_x[0]*jacobTotInv[0]*hamil[16]*rdz2)/q_); 
  alphay[8] = 1.732050807568877*((0.125*b_z[0]*jacobTotInv[0]*hamil[16]*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarYdBmag[8]*rdvpar2)/m_)*rdy2; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphay[2]-1.414213562373095*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphay[2]+1.414213562373095*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*(0.3061862178478971*(alphay[8]+alphay[6])-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[8]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[8])+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphay[8]+alphay[6]))+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphay[8]+alphay[6])-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[8]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[8])+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphay[8]+alphay[6]))+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphay[8]+alphay[6])-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[8]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[8])+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphay[8]+alphay[6]))+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphay[8]+alphay[6])-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[8]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[8])+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphay[8]+alphay[6]))+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*((-0.3061862178478971*(alphay[8]+alphay[6]))-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[8])+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[8]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphay[8]+alphay[6])+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphay[8]+alphay[6]))-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[8])+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[8]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphay[8]+alphay[6])+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphay[8]+alphay[6]))-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[8])+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[8]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphay[8]+alphay[6])+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphay[8]+alphay[6]))-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[8])+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[8]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphay[8]+alphay[6])+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphaz[32]; 
  alphaz[0] = 1.732050807568877*((jacobTotInv[0]*(0.125*b_x[0]*hamil[2]*rdy2-0.125*b_y[0]*hamil[1]*rdx2))/q_+(0.1767766952966368*BstarZdBmag[0]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[1] = 1.732050807568877*((0.125*b_x[0]*jacobTotInv[0]*hamil[6]*rdy2)/q_+(0.1767766952966368*BstarZdBmag[1]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[2] = 1.732050807568877*((0.1767766952966368*BstarZdBmag[2]*hamil[4]*rdvpar2)/m_-(0.125*b_y[0]*jacobTotInv[0]*hamil[6]*rdx2)/q_)*rdz2; 
  alphaz[3] = 1.732050807568877*((jacobTotInv[0]*(0.125*b_x[0]*hamil[8]*rdy2-0.125*b_y[0]*hamil[7]*rdx2))/q_+(0.1767766952966368*BstarZdBmag[3]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[7] = 1.732050807568877*((0.125*b_x[0]*jacobTotInv[0]*hamil[16]*rdy2)/q_+(0.1767766952966368*hamil[4]*BstarZdBmag[7]*rdvpar2)/m_)*rdz2; 
  alphaz[8] = 1.732050807568877*((0.1767766952966368*hamil[4]*BstarZdBmag[8]*rdvpar2)/m_-(0.125*b_y[0]*jacobTotInv[0]*hamil[16]*rdx2)/q_)*rdz2; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphaz[3]-1.414213562373095*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphaz[3]+1.414213562373095*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*alphaz[7]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])-0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*alphaz[7]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])-0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*alphaz[7]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])-0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*alphaz[7]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])-0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[8]+alphaz[7]))+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[8]+alphaz[7]))+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[8]+alphaz[7]))+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[8]+alphaz[7]))+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[32]; 
  alphavpar[0] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[2]*hamil[8]+BstarZdBmag[1]*hamil[7]+BstarZdBmag[0]*hamil[3])*rdz2+(BstarYdBmag[3]*hamil[8]+BstarYdBmag[1]*hamil[6]+BstarYdBmag[0]*hamil[2])*rdy2+(BstarXdBmag[3]*hamil[7]+BstarXdBmag[2]*hamil[6]+BstarXdBmag[0]*hamil[1])*rdx2))/m_; 
  alphavpar[1] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[2]*hamil[16]+BstarZdBmag[0]*hamil[7]+BstarZdBmag[1]*hamil[3])*rdz2+(BstarYdBmag[3]*hamil[16]+BstarYdBmag[0]*hamil[6]+BstarYdBmag[1]*hamil[2])*rdy2+(BstarXdBmag[7]*hamil[7]+BstarXdBmag[6]*hamil[6]+BstarXdBmag[1]*hamil[1])*rdx2))/m_; 
  alphavpar[2] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[1]*hamil[16]+BstarZdBmag[0]*hamil[8]+BstarZdBmag[2]*hamil[3])*rdz2+(BstarYdBmag[8]*hamil[8]+BstarYdBmag[6]*hamil[6]+BstarYdBmag[2]*hamil[2])*rdy2+(BstarXdBmag[3]*hamil[16]+BstarXdBmag[0]*hamil[6]+hamil[1]*BstarXdBmag[2])*rdx2))/m_; 
  alphavpar[3] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[8]*hamil[8]+BstarZdBmag[7]*hamil[7]+BstarZdBmag[3]*hamil[3])*rdz2+(BstarYdBmag[1]*hamil[16]+BstarYdBmag[0]*hamil[8]+hamil[2]*BstarYdBmag[3])*rdy2+(BstarXdBmag[2]*hamil[16]+BstarXdBmag[0]*hamil[7]+hamil[1]*BstarXdBmag[3])*rdx2))/m_; 
  alphavpar[6] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[0]*hamil[16]+BstarZdBmag[1]*hamil[8]+BstarZdBmag[2]*hamil[7])*rdz2+(BstarYdBmag[8]*hamil[16]+BstarYdBmag[2]*hamil[6]+hamil[2]*BstarYdBmag[6])*rdy2+(BstarXdBmag[7]*hamil[16]+BstarXdBmag[1]*hamil[6]+hamil[1]*BstarXdBmag[6])*rdx2))/m_; 
  alphavpar[7] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[8]*hamil[16]+BstarZdBmag[3]*hamil[7]+hamil[3]*BstarZdBmag[7])*rdz2+(BstarYdBmag[0]*hamil[16]+BstarYdBmag[1]*hamil[8]+BstarYdBmag[3]*hamil[6])*rdy2+(BstarXdBmag[6]*hamil[16]+BstarXdBmag[1]*hamil[7]+hamil[1]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[8] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[7]*hamil[16]+BstarZdBmag[3]*hamil[8]+hamil[3]*BstarZdBmag[8])*rdz2+(BstarYdBmag[6]*hamil[16]+BstarYdBmag[2]*hamil[8]+hamil[2]*BstarYdBmag[8])*rdy2+(BstarXdBmag[0]*hamil[16]+BstarXdBmag[2]*hamil[7]+BstarXdBmag[3]*hamil[6])*rdx2))/m_; 
  alphavpar[16] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[3]*hamil[16]+BstarZdBmag[7]*hamil[8]+hamil[7]*BstarZdBmag[8])*rdz2+(BstarYdBmag[2]*hamil[16]+BstarYdBmag[6]*hamil[8]+hamil[6]*BstarYdBmag[8])*rdy2+(BstarXdBmag[1]*hamil[16]+BstarXdBmag[6]*hamil[7]+hamil[6]*BstarXdBmag[7])*rdx2))/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.0883883476483184*alphavpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0883883476483184*alphavpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[16])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphavpar[16]+alphavpar[8])-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[3])+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[16]+alphavpar[8]+alphavpar[7]))+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[3])-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[16]+alphavpar[8]))+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[16])+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0055242717280199*(alphavpar[16]+alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[16])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphavpar[16]+alphavpar[8])-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[3])+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[16]+alphavpar[8]+alphavpar[7]))+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[3])-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[16]+alphavpar[8]))+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[16])+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0055242717280199*(alphavpar[16]+alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[16])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphavpar[16]+alphavpar[8])-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[3])+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[16]+alphavpar[8]+alphavpar[7]))+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[3])-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[16]+alphavpar[8]))+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[16])+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0055242717280199*(alphavpar[16]+alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[16])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphavpar[16]+alphavpar[8])-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[3])+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[16]+alphavpar[8]+alphavpar[7]))+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[3])-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[16]+alphavpar[8]))+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[16])+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0055242717280199*(alphavpar[16]+alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.3061862178478971*(alphax[7]*f[7]+alphax[6]*f[6]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.3061862178478971*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.3061862178478971*(alphaz[8]*f[8]+alphaz[7]*f[7]+alphaz[3]*f[3]+alphaz[2]*f[2]+alphaz[1]*f[1]+alphaz[0]*f[0]); 
  out[4] += 0.3061862178478971*(alphavpar[16]*f[16]+alphavpar[8]*f[8]+alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[6] += 0.3061862178478971*((alphay[8]+alphax[7])*f[16]+alphax[3]*f[8]+alphay[3]*f[7]+(alphay[2]+alphax[1])*f[6]+f[2]*alphay[6]+f[1]*alphax[6]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[7] += 0.3061862178478971*((alphaz[8]+alphax[6])*f[16]+alphax[2]*f[8]+(alphaz[3]+alphax[1])*f[7]+f[3]*alphaz[7]+f[1]*alphax[7]+alphaz[2]*f[6]+alphax[0]*f[3]+f[0]*alphax[3]+alphaz[0]*f[1]+f[0]*alphaz[1]); 
  out[8] += 0.3061862178478971*((alphaz[7]+alphay[6])*f[16]+(alphaz[3]+alphay[2])*f[8]+f[3]*alphaz[8]+f[2]*alphay[8]+alphay[1]*f[7]+alphaz[1]*f[6]+alphay[0]*f[3]+f[0]*alphay[3]+alphaz[0]*f[2]+f[0]*alphaz[2]); 
  out[9] += 0.3061862178478971*(alphax[7]*f[18]+alphax[6]*f[17]+alphavpar[8]*f[16]+f[8]*alphavpar[16]+alphax[3]*f[11]+alphax[2]*f[10]+alphax[1]*f[9]+alphavpar[3]*f[7]+f[3]*alphavpar[7]+alphavpar[2]*f[6]+f[2]*alphavpar[6]+alphax[0]*f[4]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[10] += 0.3061862178478971*(alphay[8]*f[19]+alphay[6]*f[17]+alphavpar[7]*f[16]+f[7]*alphavpar[16]+alphay[3]*f[11]+alphay[2]*f[10]+alphay[1]*f[9]+alphavpar[3]*f[8]+f[3]*alphavpar[8]+alphavpar[1]*f[6]+f[1]*alphavpar[6]+alphay[0]*f[4]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[11] += 0.3061862178478971*(alphaz[8]*f[19]+alphaz[7]*f[18]+alphavpar[6]*f[16]+f[6]*alphavpar[16]+alphaz[3]*f[11]+alphaz[2]*f[10]+alphaz[1]*f[9]+alphavpar[2]*f[8]+f[2]*alphavpar[8]+alphavpar[1]*f[7]+f[1]*alphavpar[7]+alphaz[0]*f[4]+alphavpar[0]*f[3]+f[0]*alphavpar[3]); 
  out[12] += 0.3061862178478971*(alphax[7]*f[21]+alphax[6]*f[20]+alphax[3]*f[14]+alphax[2]*f[13]+alphax[1]*f[12]+alphax[0]*f[5]); 
  out[13] += 0.3061862178478971*(alphay[8]*f[22]+alphay[6]*f[20]+alphay[3]*f[14]+alphay[2]*f[13]+alphay[1]*f[12]+alphay[0]*f[5]); 
  out[14] += 0.3061862178478971*(alphaz[8]*f[22]+alphaz[7]*f[21]+alphaz[3]*f[14]+alphaz[2]*f[13]+alphaz[1]*f[12]+alphaz[0]*f[5]); 
  out[15] += 0.3061862178478971*(alphavpar[16]*f[27]+alphavpar[8]*f[22]+alphavpar[7]*f[21]+alphavpar[6]*f[20]+alphavpar[3]*f[14]+alphavpar[2]*f[13]+alphavpar[1]*f[12]+alphavpar[0]*f[5]); 
  out[16] += 0.3061862178478971*((alphaz[3]+alphay[2]+alphax[1])*f[16]+(alphaz[7]+alphay[6]+alphax[0])*f[8]+f[7]*alphaz[8]+f[6]*alphay[8]+(alphax[6]+alphay[0])*f[7]+f[6]*(alphax[7]+alphaz[0])+(alphax[2]+alphay[1])*f[3]+f[1]*alphay[3]+f[2]*(alphax[3]+alphaz[1])+f[1]*alphaz[2]); 
  out[17] += 0.3061862178478971*((alphay[8]+alphax[7])*f[26]+alphax[3]*f[19]+alphay[3]*f[18]+(alphay[2]+alphax[1])*f[17]+alphavpar[3]*f[16]+f[3]*alphavpar[16]+(alphay[6]+alphax[0])*f[10]+(alphax[6]+alphay[0])*f[9]+alphavpar[7]*f[8]+f[7]*alphavpar[8]+alphavpar[0]*f[6]+f[0]*alphavpar[6]+(alphax[2]+alphay[1])*f[4]+alphavpar[1]*f[2]+f[1]*alphavpar[2]); 
  out[18] += 0.3061862178478971*((alphaz[8]+alphax[6])*f[26]+alphax[2]*f[19]+(alphaz[3]+alphax[1])*f[18]+alphaz[2]*f[17]+alphavpar[2]*f[16]+f[2]*alphavpar[16]+(alphaz[7]+alphax[0])*f[11]+(alphax[7]+alphaz[0])*f[9]+alphavpar[6]*f[8]+f[6]*alphavpar[8]+alphavpar[0]*f[7]+f[0]*alphavpar[7]+(alphax[3]+alphaz[1])*f[4]+alphavpar[1]*f[3]+f[1]*alphavpar[3]); 
  out[19] += 0.3061862178478971*((alphaz[7]+alphay[6])*f[26]+(alphaz[3]+alphay[2])*f[19]+alphay[1]*f[18]+alphaz[1]*f[17]+alphavpar[1]*f[16]+f[1]*alphavpar[16]+(alphaz[8]+alphay[0])*f[11]+(alphay[8]+alphaz[0])*f[10]+alphavpar[0]*f[8]+f[0]*alphavpar[8]+alphavpar[6]*f[7]+f[6]*alphavpar[7]+(alphay[3]+alphaz[2])*f[4]+alphavpar[2]*f[3]+f[2]*alphavpar[3]); 
  out[20] += 0.3061862178478971*((alphay[8]+alphax[7])*f[27]+alphax[3]*f[22]+alphay[3]*f[21]+(alphay[2]+alphax[1])*f[20]+(alphay[6]+alphax[0])*f[13]+(alphax[6]+alphay[0])*f[12]+(alphax[2]+alphay[1])*f[5]); 
  out[21] += 0.3061862178478971*((alphaz[8]+alphax[6])*f[27]+alphax[2]*f[22]+(alphaz[3]+alphax[1])*f[21]+alphaz[2]*f[20]+(alphaz[7]+alphax[0])*f[14]+(alphax[7]+alphaz[0])*f[12]+(alphax[3]+alphaz[1])*f[5]); 
  out[22] += 0.3061862178478971*((alphaz[7]+alphay[6])*f[27]+(alphaz[3]+alphay[2])*f[22]+alphay[1]*f[21]+alphaz[1]*f[20]+(alphaz[8]+alphay[0])*f[14]+(alphay[8]+alphaz[0])*f[13]+(alphay[3]+alphaz[2])*f[5]); 
  out[23] += 0.3061862178478971*(alphax[7]*f[29]+alphax[6]*f[28]+alphavpar[8]*f[27]+alphax[3]*f[25]+alphax[2]*f[24]+alphax[1]*f[23]+alphavpar[16]*f[22]+alphavpar[3]*f[21]+alphavpar[2]*f[20]+alphax[0]*f[15]+alphavpar[7]*f[14]+alphavpar[6]*f[13]+alphavpar[0]*f[12]+alphavpar[1]*f[5]); 
  out[24] += 0.3061862178478971*(alphay[8]*f[30]+alphay[6]*f[28]+alphavpar[7]*f[27]+alphay[3]*f[25]+alphay[2]*f[24]+alphay[1]*f[23]+alphavpar[3]*f[22]+alphavpar[16]*f[21]+alphavpar[1]*f[20]+alphay[0]*f[15]+alphavpar[8]*f[14]+alphavpar[0]*f[13]+alphavpar[6]*f[12]+alphavpar[2]*f[5]); 
  out[25] += 0.3061862178478971*(alphaz[8]*f[30]+alphaz[7]*f[29]+alphavpar[6]*f[27]+alphaz[3]*f[25]+alphaz[2]*f[24]+alphaz[1]*f[23]+alphavpar[2]*f[22]+alphavpar[1]*f[21]+alphavpar[16]*f[20]+alphaz[0]*f[15]+alphavpar[0]*f[14]+alphavpar[8]*f[13]+alphavpar[7]*f[12]+alphavpar[3]*f[5]); 
  out[26] += 0.3061862178478971*((alphaz[3]+alphay[2]+alphax[1])*f[26]+(alphaz[7]+alphay[6]+alphax[0])*f[19]+(alphaz[8]+alphax[6]+alphay[0])*f[18]+(alphay[8]+alphax[7]+alphaz[0])*f[17]+alphavpar[0]*f[16]+f[0]*alphavpar[16]+(alphax[2]+alphay[1])*f[11]+(alphax[3]+alphaz[1])*f[10]+(alphay[3]+alphaz[2])*f[9]+alphavpar[1]*f[8]+f[1]*alphavpar[8]+alphavpar[2]*f[7]+f[2]*alphavpar[7]+alphavpar[3]*f[6]+f[3]*alphavpar[6]); 
  out[27] += 0.3061862178478971*((alphaz[3]+alphay[2]+alphax[1])*f[27]+(alphaz[7]+alphay[6]+alphax[0])*f[22]+(alphaz[8]+alphax[6]+alphay[0])*f[21]+(alphay[8]+alphax[7]+alphaz[0])*f[20]+(alphax[2]+alphay[1])*f[14]+(alphax[3]+alphaz[1])*f[13]+(alphay[3]+alphaz[2])*f[12]); 
  out[28] += 0.3061862178478971*((alphay[8]+alphax[7])*f[31]+alphax[3]*f[30]+alphay[3]*f[29]+(alphay[2]+alphax[1])*f[28]+alphavpar[3]*f[27]+(alphay[6]+alphax[0])*f[24]+(alphax[6]+alphay[0])*f[23]+alphavpar[7]*f[22]+alphavpar[8]*f[21]+alphavpar[0]*f[20]+f[14]*alphavpar[16]+(alphax[2]+alphay[1])*f[15]+alphavpar[1]*f[13]+alphavpar[2]*f[12]+f[5]*alphavpar[6]); 
  out[29] += 0.3061862178478971*((alphaz[8]+alphax[6])*f[31]+alphax[2]*f[30]+(alphaz[3]+alphax[1])*f[29]+alphaz[2]*f[28]+alphavpar[2]*f[27]+(alphaz[7]+alphax[0])*f[25]+(alphax[7]+alphaz[0])*f[23]+alphavpar[6]*f[22]+alphavpar[0]*f[21]+alphavpar[8]*f[20]+f[13]*alphavpar[16]+(alphax[3]+alphaz[1])*f[15]+alphavpar[1]*f[14]+alphavpar[3]*f[12]+f[5]*alphavpar[7]); 
  out[30] += 0.3061862178478971*((alphaz[7]+alphay[6])*f[31]+(alphaz[3]+alphay[2])*f[30]+alphay[1]*f[29]+alphaz[1]*f[28]+alphavpar[1]*f[27]+(alphaz[8]+alphay[0])*f[25]+(alphay[8]+alphaz[0])*f[24]+alphavpar[0]*f[22]+alphavpar[6]*f[21]+alphavpar[7]*f[20]+f[12]*alphavpar[16]+(alphay[3]+alphaz[2])*f[15]+alphavpar[2]*f[14]+alphavpar[3]*f[13]+f[5]*alphavpar[8]); 
  out[31] += 0.3061862178478971*((alphaz[3]+alphay[2]+alphax[1])*f[31]+(alphaz[7]+alphay[6]+alphax[0])*f[30]+(alphaz[8]+alphax[6]+alphay[0])*f[29]+(alphay[8]+alphax[7]+alphaz[0])*f[28]+alphavpar[0]*f[27]+(alphax[2]+alphay[1])*f[25]+(alphax[3]+alphaz[1])*f[24]+(alphay[3]+alphaz[2])*f[23]+alphavpar[1]*f[22]+alphavpar[2]*f[21]+alphavpar[3]*f[20]+f[5]*alphavpar[16]+alphavpar[6]*f[14]+alphavpar[7]*f[13]+alphavpar[8]*f[12]); 
  return cflFreq; 
} 
double EmGyrokineticGenGeoVol3x2vSerP1_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f, double *out) 
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
  double wz = w[2];
  double rdz2 = 2.0/dxv[2];
  double wvpar = w[3];
  double rdvpar2 = 2.0/dxv[3];
  double wmu = w[4];
  double rdmu2 = 2.0/dxv[4];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wzSq = w[2]*w[2];
  double rdz2Sq = rdz2*rdz2;
  double wvparSq = w[3]*w[3];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[4]*w[4];
  double rdmu2Sq = rdmu2*rdmu2;

  double hamil[32]; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(4.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+4.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = 2.0*phi[3]*q_; 
  hamil[4] = (3.265986323710906*m_*wvpar)/rdvpar2; 
  hamil[5] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[6] = 2.0*phi[4]*q_; 
  hamil[7] = 2.0*phi[5]*q_; 
  hamil[8] = 2.0*phi[6]*q_; 
  hamil[12] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[16] = 2.0*phi[7]*q_; 

  double BstarXdBmag[32]; 
  BstarXdBmag[0] = -0.4330127018922193*(((b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*Apar[5]+(b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*Apar[3])*rdz2-1.0*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[4]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])*rdy2); 
  BstarXdBmag[1] = -0.08660254037844387*(((9.0*b_y[1]*jacobTotInv[1]+5.0*b_y[0]*jacobTotInv[0])*Apar[5]+5.0*(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*Apar[3])*rdz2-1.0*((9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[4]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])*rdy2); 
  BstarXdBmag[2] = -0.4330127018922193*((b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*Apar[7]+(b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*Apar[6])*rdz2; 
  BstarXdBmag[3] = 0.4330127018922193*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[7]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[6])*rdy2; 
  BstarXdBmag[6] = -0.08660254037844387*((9.0*b_y[1]*jacobTotInv[1]+5.0*b_y[0]*jacobTotInv[0])*Apar[7]+5.0*(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*Apar[6])*rdz2; 
  BstarXdBmag[7] = 0.08660254037844387*((9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[7]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[6])*rdy2; 

  double BstarYdBmag[32]; 
  BstarYdBmag[0] = -(0.3061862178478971*(4.0*jacobTotInv[0]*b_z[1]*m_*rdx2*wvpar+1.414213562373095*q_*((2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*rdx2-1.0*((b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[5]+(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*Apar[3])*rdz2)))/q_; 
  BstarYdBmag[1] = -(0.06123724356957942*(20.0*b_z[1]*jacobTotInv[1]*m_*rdx2*wvpar+1.414213562373095*q_*(5.0*((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])*rdx2-1.0*((9.0*b_x[1]*jacobTotInv[1]+5.0*b_x[0]*jacobTotInv[0])*Apar[5]+5.0*(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[3])*rdz2)))/q_; 
  BstarYdBmag[2] = 0.4330127018922193*(((b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[7]+(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*Apar[6])*rdz2-1.0*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[4]+jacobTotInv[0]*b_z[1]*Apar[2])*rdx2); 
  BstarYdBmag[3] = -0.4330127018922193*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[5]+jacobTotInv[0]*b_z[1]*Apar[3])*rdx2; 
  BstarYdBmag[4] = -(0.7071067811865475*jacobTotInv[0]*b_z[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[6] = 0.08660254037844387*(((9.0*b_x[1]*jacobTotInv[1]+5.0*b_x[0]*jacobTotInv[0])*Apar[7]+5.0*(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[6])*rdz2-5.0*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[4]+b_z[1]*jacobTotInv[1]*Apar[2])*rdx2); 
  BstarYdBmag[7] = -0.4330127018922193*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[5]+b_z[1]*jacobTotInv[1]*Apar[3])*rdx2; 
  BstarYdBmag[8] = -0.4330127018922193*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[7]+jacobTotInv[0]*b_z[1]*Apar[6])*rdx2; 
  BstarYdBmag[9] = -(0.7071067811865475*b_z[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[16] = -0.4330127018922193*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[7]+b_z[1]*jacobTotInv[1]*Apar[6])*rdx2; 

  double BstarZdBmag[32]; 
  BstarZdBmag[0] = (0.1767766952966368*(6.928203230275509*jacobTotInv[0]*b_y[1]*m_*rdx2*wvpar+q_*(2.449489742783178*((2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2-1.0*((b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[4]+(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*Apar[2])*rdy2)+4.0*(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0]))))/q_; 
  BstarZdBmag[1] = (0.03535533905932736*(34.64101615137754*b_y[1]*jacobTotInv[1]*m_*rdx2*wvpar+q_*(5.0*(2.449489742783178*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2+4.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))-2.449489742783178*((9.0*b_x[1]*jacobTotInv[1]+5.0*b_x[0]*jacobTotInv[0])*Apar[4]+5.0*(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[2])*rdy2)))/q_; 
  BstarZdBmag[2] = 0.4330127018922193*((2.0*b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*Apar[4]+jacobTotInv[0]*b_y[1]*Apar[2])*rdx2; 
  BstarZdBmag[3] = -0.4330127018922193*(((b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[7]+(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*Apar[6])*rdy2-1.0*((2.0*b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*Apar[5]+jacobTotInv[0]*b_y[1]*Apar[3])*rdx2); 
  BstarZdBmag[4] = (0.7071067811865475*jacobTotInv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[6] = 0.4330127018922193*((b_y[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_y[1])*Apar[4]+b_y[1]*jacobTotInv[1]*Apar[2])*rdx2; 
  BstarZdBmag[7] = -0.08660254037844387*(((9.0*b_x[1]*jacobTotInv[1]+5.0*b_x[0]*jacobTotInv[0])*Apar[7]+5.0*(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[6])*rdy2-5.0*((b_y[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_y[1])*Apar[5]+b_y[1]*jacobTotInv[1]*Apar[3])*rdx2); 
  BstarZdBmag[8] = 0.4330127018922193*((2.0*b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*Apar[7]+jacobTotInv[0]*b_y[1]*Apar[6])*rdx2; 
  BstarZdBmag[9] = (0.7071067811865475*b_y[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[16] = 0.4330127018922193*((b_y[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_y[1])*Apar[7]+b_y[1]*jacobTotInv[1]*Apar[6])*rdx2; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[32]; 
  alphax[0] = 1.732050807568877*rdx2*((0.125*((b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[7]+(b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[3])*rdz2-0.125*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[6]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[2])*rdy2)/q_+(0.1767766952966368*BstarXdBmag[0]*hamil[4]*rdvpar2)/m_); 
  alphax[1] = 1.732050807568877*rdx2*(((0.125*(b_y[0]*jacobTotInv[0]*hamil[7]+(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[3])+0.225*b_y[1]*jacobTotInv[1]*hamil[7])*rdz2+((-0.125*(b_z[0]*jacobTotInv[0]*hamil[6]+(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[2]))-0.225*b_z[1]*jacobTotInv[1]*hamil[6])*rdy2)/q_+(0.1767766952966368*BstarXdBmag[1]*hamil[4]*rdvpar2)/m_); 
  alphax[2] = 1.732050807568877*rdx2*((0.125*((b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[16]+(b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[8])*rdz2)/q_+(0.1767766952966368*BstarXdBmag[2]*hamil[4]*rdvpar2)/m_); 
  alphax[3] = 1.732050807568877*rdx2*((0.1767766952966368*BstarXdBmag[3]*hamil[4]*rdvpar2)/m_-(0.125*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[16]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[8])*rdy2)/q_); 
  alphax[6] = 1.732050807568877*rdx2*(((0.125*(b_y[0]*jacobTotInv[0]*hamil[16]+(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[8])+0.225*b_y[1]*jacobTotInv[1]*hamil[16])*rdz2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[6]*rdvpar2)/m_); 
  alphax[7] = 1.732050807568877*rdx2*((((-0.125*(b_z[0]*jacobTotInv[0]*hamil[16]+(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[8]))-0.225*b_z[1]*jacobTotInv[1]*hamil[16])*rdy2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[7]*rdvpar2)/m_); 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphax[1]-1.414213562373095*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[7]+alphax[6]))-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[7])+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[32]; 
  alphay[0] = 1.732050807568877*rdy2*((0.125*hamil[1]*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*rdx2-0.125*((b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[7]+(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*hamil[3])*rdz2)/q_+(0.1767766952966368*BstarYdBmag[0]*hamil[4]*rdvpar2)/m_); 
  alphay[1] = 1.732050807568877*rdy2*((((-0.125*(b_x[0]*jacobTotInv[0]*hamil[7]+(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[3]))-0.225*b_x[1]*jacobTotInv[1]*hamil[7])*rdz2+0.125*hamil[1]*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*rdx2)/q_+(0.1767766952966368*BstarYdBmag[1]*hamil[4]*rdvpar2)/m_); 
  alphay[2] = 1.732050807568877*rdy2*((0.125*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[6]*rdx2-0.125*((b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[16]+(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*hamil[8])*rdz2)/q_+(0.1767766952966368*BstarYdBmag[2]*hamil[4]*rdvpar2)/m_); 
  alphay[3] = 1.732050807568877*((0.125*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[7]*rdx2)/q_+(0.1767766952966368*BstarYdBmag[3]*hamil[4]*rdvpar2)/m_)*rdy2; 
  alphay[4] = (0.3061862178478971*BstarYdBmag[4]*hamil[4]*rdvpar2*rdy2)/m_; 
  alphay[5] = (0.2165063509461096*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[12]*rdx2*rdy2)/q_; 
  alphay[6] = 1.732050807568877*rdy2*((((-0.125*(b_x[0]*jacobTotInv[0]*hamil[16]+(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[8]))-0.225*b_x[1]*jacobTotInv[1]*hamil[16])*rdz2+0.125*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[6]*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarYdBmag[6]*rdvpar2)/m_); 
  alphay[7] = 1.732050807568877*((0.125*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[7]*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarYdBmag[7]*rdvpar2)/m_)*rdy2; 
  alphay[8] = 1.732050807568877*((0.125*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[16]*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarYdBmag[8]*rdvpar2)/m_)*rdy2; 
  alphay[9] = (0.3061862178478971*hamil[4]*BstarYdBmag[9]*rdvpar2*rdy2)/m_; 
  alphay[12] = (0.2165063509461096*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[12]*rdx2*rdy2)/q_; 
  alphay[16] = 1.732050807568877*((0.125*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[16]*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarYdBmag[16]*rdvpar2)/m_)*rdy2; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphay[2]-1.414213562373095*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphay[2]+1.414213562373095*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*(alphay[12]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[12]+alphay[9])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[12]+alphay[9])-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*(alphay[12]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*(alphay[12]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[12]+alphay[9])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[12]+alphay[9])-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*(alphay[12]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[12]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*(alphay[12]+alphay[9])-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*(alphay[12]+alphay[9])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[12]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[12]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*(alphay[12]+alphay[9])-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*(alphay[12]+alphay[9])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[12]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphaz[32]; 
  alphaz[0] = 1.732050807568877*((0.125*((b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[6]+(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*hamil[2])*rdy2-0.125*hamil[1]*(b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*rdx2)/q_+(0.1767766952966368*BstarZdBmag[0]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[1] = 1.732050807568877*(((0.125*(b_x[0]*jacobTotInv[0]*hamil[6]+(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[2])+0.225*b_x[1]*jacobTotInv[1]*hamil[6])*rdy2-0.125*hamil[1]*(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*rdx2)/q_+(0.1767766952966368*BstarZdBmag[1]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[2] = 1.732050807568877*((0.1767766952966368*BstarZdBmag[2]*hamil[4]*rdvpar2)/m_-(0.125*(b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[6]*rdx2)/q_)*rdz2; 
  alphaz[3] = 1.732050807568877*((0.125*((b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[16]+(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*hamil[8])*rdy2-0.125*(b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[7]*rdx2)/q_+(0.1767766952966368*BstarZdBmag[3]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[4] = (0.3061862178478971*BstarZdBmag[4]*hamil[4]*rdvpar2*rdz2)/m_; 
  alphaz[5] = -(0.2165063509461096*(b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[12]*rdx2*rdz2)/q_; 
  alphaz[6] = 1.732050807568877*((0.1767766952966368*hamil[4]*BstarZdBmag[6]*rdvpar2)/m_-(0.125*(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[6]*rdx2)/q_)*rdz2; 
  alphaz[7] = 1.732050807568877*(((0.125*(b_x[0]*jacobTotInv[0]*hamil[16]+(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[8])+0.225*b_x[1]*jacobTotInv[1]*hamil[16])*rdy2-0.125*(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[7]*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarZdBmag[7]*rdvpar2)/m_)*rdz2; 
  alphaz[8] = 1.732050807568877*((0.1767766952966368*hamil[4]*BstarZdBmag[8]*rdvpar2)/m_-(0.125*(b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[16]*rdx2)/q_)*rdz2; 
  alphaz[9] = (0.3061862178478971*hamil[4]*BstarZdBmag[9]*rdvpar2*rdz2)/m_; 
  alphaz[12] = -(0.2165063509461096*(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[12]*rdx2*rdz2)/q_; 
  alphaz[16] = 1.732050807568877*((0.1767766952966368*hamil[4]*BstarZdBmag[16]*rdvpar2)/m_-(0.125*(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[16]*rdx2)/q_)*rdz2; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphaz[3]-1.414213562373095*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphaz[3]+1.414213562373095*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.3061862178478971*alphaz[16])+0.1767766952966368*(alphaz[12]+alphaz[9])+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*(alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[16]-0.1767766952966368*(alphaz[12]+alphaz[9])+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[16]+0.1767766952966368*(alphaz[12]+alphaz[9])-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[16])-0.1767766952966368*(alphaz[12]+alphaz[9])-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*(alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[16])+0.1767766952966368*alphaz[12]-0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*alphaz[5]+0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[16]-0.1767766952966368*alphaz[12]+0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5])+0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[16]+0.1767766952966368*alphaz[12]-0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5])+0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[16])-0.1767766952966368*alphaz[12]+0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*alphaz[5]+0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[16])-0.1767766952966368*alphaz[12]+0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5])-0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[16]+0.1767766952966368*alphaz[12]-0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*alphaz[5]-0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[16]-0.1767766952966368*alphaz[12]+0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*alphaz[5]-0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[16])+0.1767766952966368*alphaz[12]-0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5])-0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[16])-0.1767766952966368*(alphaz[12]+alphaz[9])+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[16]+0.1767766952966368*(alphaz[12]+alphaz[9])+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*(alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[16]-0.1767766952966368*(alphaz[12]+alphaz[9])-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*(alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[16])+0.1767766952966368*(alphaz[12]+alphaz[9])-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.3061862178478971*alphaz[16]+0.1767766952966368*(alphaz[12]+alphaz[9])-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*(alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[16])-0.1767766952966368*(alphaz[12]+alphaz[9])-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[16])+0.1767766952966368*(alphaz[12]+alphaz[9])+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[16]-0.1767766952966368*(alphaz[12]+alphaz[9])+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*(alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[16]+0.1767766952966368*alphaz[12]-0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*alphaz[5]+0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[16])-0.1767766952966368*alphaz[12]+0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5])+0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[16])+0.1767766952966368*alphaz[12]-0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5])+0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[16]-0.1767766952966368*alphaz[12]+0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*alphaz[5]+0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[16]-0.1767766952966368*alphaz[12]+0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5])-0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[16])+0.1767766952966368*alphaz[12]-0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*alphaz[5]-0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[16])-0.1767766952966368*alphaz[12]+0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*alphaz[5]-0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[16]+0.1767766952966368*alphaz[12]-0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5])-0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[16]-0.1767766952966368*(alphaz[12]+alphaz[9])-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[16])+0.1767766952966368*(alphaz[12]+alphaz[9])-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*(alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[16])-0.1767766952966368*(alphaz[12]+alphaz[9])+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*(alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[16]+0.1767766952966368*(alphaz[12]+alphaz[9])+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[32]; 
  alphavpar[0] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[6]*hamil[16]+BstarZdBmag[2]*hamil[8]+BstarZdBmag[1]*hamil[7]+BstarZdBmag[0]*hamil[3])*rdz2+(BstarYdBmag[7]*hamil[16]+BstarYdBmag[3]*hamil[8]+BstarYdBmag[1]*hamil[6]+BstarYdBmag[0]*hamil[2])*rdy2+(BstarXdBmag[3]*hamil[7]+BstarXdBmag[2]*hamil[6]+BstarXdBmag[0]*hamil[1])*rdx2))/m_; 
  alphavpar[1] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[2]*hamil[16]+BstarZdBmag[6]*hamil[8]+BstarZdBmag[0]*hamil[7]+BstarZdBmag[1]*hamil[3])*rdz2+(BstarYdBmag[3]*hamil[16]+BstarYdBmag[7]*hamil[8]+BstarYdBmag[0]*hamil[6]+BstarYdBmag[1]*hamil[2])*rdy2+(BstarXdBmag[7]*hamil[7]+BstarXdBmag[6]*hamil[6]+BstarXdBmag[1]*hamil[1])*rdx2))/m_; 
  alphavpar[2] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[1]*hamil[16]+BstarZdBmag[0]*hamil[8]+BstarZdBmag[6]*hamil[7]+BstarZdBmag[2]*hamil[3])*rdz2+(BstarYdBmag[16]*hamil[16]+BstarYdBmag[8]*hamil[8]+BstarYdBmag[6]*hamil[6]+BstarYdBmag[2]*hamil[2])*rdy2+(BstarXdBmag[3]*hamil[16]+BstarXdBmag[0]*hamil[6]+hamil[1]*BstarXdBmag[2])*rdx2))/m_; 
  alphavpar[3] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[16]*hamil[16]+BstarZdBmag[8]*hamil[8]+BstarZdBmag[7]*hamil[7]+BstarZdBmag[3]*hamil[3])*rdz2+(BstarYdBmag[1]*hamil[16]+BstarYdBmag[0]*hamil[8]+hamil[6]*BstarYdBmag[7]+hamil[2]*BstarYdBmag[3])*rdy2+(BstarXdBmag[2]*hamil[16]+BstarXdBmag[0]*hamil[7]+hamil[1]*BstarXdBmag[3])*rdx2))/m_; 
  alphavpar[4] = -(0.3061862178478971*rdvpar2*((hamil[7]*BstarZdBmag[9]+hamil[3]*BstarZdBmag[4])*rdz2+(hamil[6]*BstarYdBmag[9]+hamil[2]*BstarYdBmag[4])*rdy2))/m_; 
  alphavpar[5] = -(0.3061862178478971*BstarXdBmag[0]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[6] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[0]*hamil[16]+BstarZdBmag[1]*hamil[8]+BstarZdBmag[2]*hamil[7]+hamil[3]*BstarZdBmag[6])*rdz2+(BstarYdBmag[8]*hamil[16]+hamil[8]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[6]+hamil[2]*BstarYdBmag[6])*rdy2+(BstarXdBmag[7]*hamil[16]+BstarXdBmag[1]*hamil[6]+hamil[1]*BstarXdBmag[6])*rdx2))/m_; 
  alphavpar[7] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[8]*hamil[16]+hamil[8]*BstarZdBmag[16]+BstarZdBmag[3]*hamil[7]+hamil[3]*BstarZdBmag[7])*rdz2+(BstarYdBmag[0]*hamil[16]+BstarYdBmag[1]*hamil[8]+hamil[2]*BstarYdBmag[7]+BstarYdBmag[3]*hamil[6])*rdy2+(BstarXdBmag[6]*hamil[16]+BstarXdBmag[1]*hamil[7]+hamil[1]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[8] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[7]*hamil[16]+hamil[7]*BstarZdBmag[16]+BstarZdBmag[3]*hamil[8]+hamil[3]*BstarZdBmag[8])*rdz2+(BstarYdBmag[6]*hamil[16]+hamil[6]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[8]+hamil[2]*BstarYdBmag[8])*rdy2+(BstarXdBmag[0]*hamil[16]+BstarXdBmag[2]*hamil[7]+BstarXdBmag[3]*hamil[6])*rdx2))/m_; 
  alphavpar[9] = -(0.3061862178478971*rdvpar2*((hamil[3]*BstarZdBmag[9]+BstarZdBmag[4]*hamil[7])*rdz2+(hamil[2]*BstarYdBmag[9]+BstarYdBmag[4]*hamil[6])*rdy2))/m_; 
  alphavpar[10] = -(0.3061862178478971*(BstarZdBmag[9]*hamil[16]+BstarZdBmag[4]*hamil[8])*rdvpar2*rdz2)/m_; 
  alphavpar[11] = -(0.3061862178478971*(BstarYdBmag[9]*hamil[16]+BstarYdBmag[4]*hamil[8])*rdvpar2*rdy2)/m_; 
  alphavpar[12] = -(0.3061862178478971*BstarXdBmag[1]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[13] = -(0.3061862178478971*BstarXdBmag[2]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[14] = -(0.3061862178478971*BstarXdBmag[3]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[16] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[3]*hamil[16]+hamil[3]*BstarZdBmag[16]+BstarZdBmag[7]*hamil[8]+hamil[7]*BstarZdBmag[8])*rdz2+(BstarYdBmag[2]*hamil[16]+hamil[2]*BstarYdBmag[16]+BstarYdBmag[6]*hamil[8]+hamil[6]*BstarYdBmag[8])*rdy2+(BstarXdBmag[1]*hamil[16]+BstarXdBmag[6]*hamil[7]+hamil[6]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[17] = -(0.3061862178478971*(BstarZdBmag[4]*hamil[16]+hamil[8]*BstarZdBmag[9])*rdvpar2*rdz2)/m_; 
  alphavpar[18] = -(0.3061862178478971*(BstarYdBmag[4]*hamil[16]+hamil[8]*BstarYdBmag[9])*rdvpar2*rdy2)/m_; 
  alphavpar[20] = -(0.3061862178478971*BstarXdBmag[6]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[21] = -(0.3061862178478971*BstarXdBmag[7]*hamil[12]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphavpar[4]-1.414213562373095*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphavpar[4]+1.414213562373095*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[21]+alphavpar[20]))-0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[21])+0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]+0.1767766952966368*(alphavpar[16]+alphavpar[14])-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[21])+0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]-0.1767766952966368*(alphavpar[16]+alphavpar[14])+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[21]+alphavpar[20]))-0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[21]+alphavpar[20]))+0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[21])+0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]-0.1767766952966368*(alphavpar[16]+alphavpar[14])+0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[21])+0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]+0.1767766952966368*(alphavpar[16]+alphavpar[14])-0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[21]+alphavpar[20]))+0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[21]+alphavpar[20]))+0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[21])+0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]+0.1767766952966368*(alphavpar[16]+alphavpar[14])-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[21])+0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]-0.1767766952966368*(alphavpar[16]+alphavpar[14])+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[21]+alphavpar[20]))+0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[21]+alphavpar[20]))-0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[21])+0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]-0.1767766952966368*(alphavpar[16]+alphavpar[14])+0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[21])+0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]+0.1767766952966368*(alphavpar[16]+alphavpar[14])-0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[21]+alphavpar[20]))-0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.3061862178478971*(alphax[7]*f[7]+alphax[6]*f[6]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.3061862178478971*(alphay[16]*f[16]+alphay[12]*f[12]+alphay[9]*f[9]+alphay[8]*f[8]+alphay[7]*f[7]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.3061862178478971*(alphaz[16]*f[16]+alphaz[12]*f[12]+alphaz[9]*f[9]+alphaz[8]*f[8]+alphaz[7]*f[7]+alphaz[6]*f[6]+alphaz[5]*f[5]+alphaz[4]*f[4]+alphaz[3]*f[3]+alphaz[2]*f[2]+alphaz[1]*f[1]+alphaz[0]*f[0]); 
  out[4] += 0.3061862178478971*(alphavpar[21]*f[21]+alphavpar[20]*f[20]+alphavpar[18]*f[18]+alphavpar[17]*f[17]+alphavpar[16]*f[16]+alphavpar[14]*f[14]+alphavpar[13]*f[13]+alphavpar[12]*f[12]+alphavpar[11]*f[11]+alphavpar[10]*f[10]+alphavpar[9]*f[9]+alphavpar[8]*f[8]+alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[6] += 0.3061862178478971*((alphay[8]+alphax[7])*f[16]+f[8]*alphay[16]+alphay[5]*f[12]+f[5]*alphay[12]+alphay[4]*f[9]+f[4]*alphay[9]+alphax[3]*f[8]+alphay[3]*f[7]+f[3]*alphay[7]+(alphay[2]+alphax[1])*f[6]+f[2]*alphay[6]+f[1]*alphax[6]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[7] += 0.3061862178478971*((alphaz[8]+alphax[6])*f[16]+f[8]*alphaz[16]+alphaz[5]*f[12]+f[5]*alphaz[12]+alphaz[4]*f[9]+f[4]*alphaz[9]+alphax[2]*f[8]+(alphaz[3]+alphax[1])*f[7]+f[3]*alphaz[7]+f[1]*alphax[7]+alphaz[2]*f[6]+f[2]*alphaz[6]+alphax[0]*f[3]+f[0]*alphax[3]+alphaz[0]*f[1]+f[0]*alphaz[1]); 
  out[8] += 0.3061862178478971*(alphay[12]*f[21]+alphaz[12]*f[20]+alphay[9]*f[18]+alphaz[9]*f[17]+(alphaz[7]+alphay[6])*f[16]+f[7]*alphaz[16]+f[6]*alphay[16]+alphay[5]*f[14]+alphaz[5]*f[13]+alphay[4]*f[11]+alphaz[4]*f[10]+(alphaz[3]+alphay[2])*f[8]+f[3]*alphaz[8]+f[2]*alphay[8]+alphay[1]*f[7]+f[1]*alphay[7]+alphaz[1]*f[6]+f[1]*alphaz[6]+alphay[0]*f[3]+f[0]*alphay[3]+alphaz[0]*f[2]+f[0]*alphaz[2]); 
  out[9] += 0.3061862178478971*(alphavpar[14]*f[21]+f[14]*alphavpar[21]+alphavpar[13]*f[20]+f[13]*alphavpar[20]+(alphavpar[11]+alphax[7])*f[18]+f[11]*alphavpar[18]+(alphavpar[10]+alphax[6])*f[17]+f[10]*alphavpar[17]+alphavpar[8]*f[16]+f[8]*alphavpar[16]+alphavpar[5]*f[12]+f[5]*alphavpar[12]+alphax[3]*f[11]+alphax[2]*f[10]+(alphavpar[4]+alphax[1])*f[9]+f[4]*alphavpar[9]+alphavpar[3]*f[7]+f[3]*alphavpar[7]+alphavpar[2]*f[6]+f[2]*alphavpar[6]+alphax[0]*f[4]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[10] += 0.3061862178478971*(alphavpar[21]*f[27]+(alphavpar[18]+alphay[16])*f[26]+alphay[12]*f[23]+alphavpar[14]*f[22]+alphavpar[12]*f[20]+f[12]*alphavpar[20]+(alphavpar[11]+alphay[8])*f[19]+alphay[7]*f[18]+(alphavpar[9]+alphay[6])*f[17]+f[9]*alphavpar[17]+alphavpar[7]*f[16]+f[7]*alphavpar[16]+alphay[5]*f[15]+alphavpar[5]*f[13]+f[5]*alphavpar[13]+alphay[3]*f[11]+(alphavpar[4]+alphay[2])*f[10]+f[4]*alphavpar[10]+alphay[1]*f[9]+f[1]*alphay[9]+alphavpar[3]*f[8]+f[3]*alphavpar[8]+alphavpar[1]*f[6]+f[1]*alphavpar[6]+alphay[0]*f[4]+f[0]*alphay[4]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[11] += 0.3061862178478971*(alphavpar[20]*f[27]+(alphavpar[17]+alphaz[16])*f[26]+alphaz[12]*f[23]+alphavpar[13]*f[22]+alphavpar[12]*f[21]+f[12]*alphavpar[21]+(alphavpar[10]+alphaz[8])*f[19]+(alphavpar[9]+alphaz[7])*f[18]+f[9]*alphavpar[18]+alphaz[6]*f[17]+alphavpar[6]*f[16]+f[6]*alphavpar[16]+alphaz[5]*f[15]+alphavpar[5]*f[14]+f[5]*alphavpar[14]+(alphavpar[4]+alphaz[3])*f[11]+f[4]*alphavpar[11]+alphaz[2]*f[10]+alphaz[1]*f[9]+f[1]*alphaz[9]+alphavpar[2]*f[8]+f[2]*alphavpar[8]+alphavpar[1]*f[7]+f[1]*alphavpar[7]+alphaz[0]*f[4]+f[0]*alphaz[4]+alphavpar[0]*f[3]+f[0]*alphavpar[3]); 
  out[12] += 0.3061862178478971*(alphax[7]*f[21]+alphax[6]*f[20]+alphax[3]*f[14]+alphax[2]*f[13]+alphax[1]*f[12]+alphax[0]*f[5]); 
  out[13] += 0.3061862178478971*(alphay[16]*f[27]+alphay[9]*f[23]+alphay[8]*f[22]+alphay[7]*f[21]+alphay[6]*f[20]+alphay[4]*f[15]+alphay[3]*f[14]+alphay[2]*f[13]+alphay[1]*f[12]+f[1]*alphay[12]+alphay[0]*f[5]+f[0]*alphay[5]); 
  out[14] += 0.3061862178478971*(alphaz[16]*f[27]+alphaz[9]*f[23]+alphaz[8]*f[22]+alphaz[7]*f[21]+alphaz[6]*f[20]+alphaz[4]*f[15]+alphaz[3]*f[14]+alphaz[2]*f[13]+alphaz[1]*f[12]+f[1]*alphaz[12]+alphaz[0]*f[5]+f[0]*alphaz[5]); 
  out[15] += 0.3061862178478971*(alphavpar[18]*f[29]+alphavpar[17]*f[28]+alphavpar[16]*f[27]+alphavpar[11]*f[25]+alphavpar[10]*f[24]+alphavpar[9]*f[23]+alphavpar[8]*f[22]+alphavpar[7]*f[21]+f[7]*alphavpar[21]+alphavpar[6]*f[20]+f[6]*alphavpar[20]+alphavpar[4]*f[15]+alphavpar[3]*f[14]+f[3]*alphavpar[14]+alphavpar[2]*f[13]+f[2]*alphavpar[13]+alphavpar[1]*f[12]+f[1]*alphavpar[12]+alphavpar[0]*f[5]+f[0]*alphavpar[5]); 
  out[16] += 0.3061862178478971*(alphay[5]*f[21]+alphaz[5]*f[20]+alphay[4]*f[18]+alphaz[4]*f[17]+(alphaz[3]+alphay[2]+alphax[1])*f[16]+f[3]*alphaz[16]+f[2]*alphay[16]+alphay[12]*f[14]+alphaz[12]*f[13]+alphay[9]*f[11]+alphaz[9]*f[10]+(alphaz[7]+alphay[6]+alphax[0])*f[8]+f[7]*alphaz[8]+f[6]*alphay[8]+(alphax[6]+alphay[0])*f[7]+f[0]*alphay[7]+f[6]*(alphax[7]+alphaz[0])+f[0]*alphaz[6]+(alphax[2]+alphay[1])*f[3]+f[1]*alphay[3]+f[2]*(alphax[3]+alphaz[1])+f[1]*alphaz[2]); 
  out[17] += 0.3061862178478971*(alphavpar[14]*f[27]+(alphavpar[11]+alphay[8]+alphax[7])*f[26]+alphay[5]*f[23]+alphavpar[21]*f[22]+alphavpar[5]*f[20]+f[5]*alphavpar[20]+(alphavpar[18]+alphay[16]+alphax[3])*f[19]+alphay[3]*f[18]+(alphavpar[4]+alphay[2]+alphax[1])*f[17]+f[4]*alphavpar[17]+alphavpar[3]*f[16]+f[3]*alphavpar[16]+alphay[12]*f[15]+alphavpar[12]*f[13]+f[12]*alphavpar[13]+alphay[7]*f[11]+(alphavpar[9]+alphay[6]+alphax[0])*f[10]+f[9]*(alphavpar[10]+alphax[6]+alphay[0])+f[0]*alphay[9]+alphavpar[7]*f[8]+f[7]*alphavpar[8]+alphavpar[0]*f[6]+f[0]*alphavpar[6]+(alphax[2]+alphay[1])*f[4]+f[1]*alphay[4]+alphavpar[1]*f[2]+f[1]*alphavpar[2]); 
  out[18] += 0.3061862178478971*(alphavpar[13]*f[27]+(alphavpar[10]+alphaz[8]+alphax[6])*f[26]+alphaz[5]*f[23]+alphavpar[20]*f[22]+alphavpar[5]*f[21]+f[5]*alphavpar[21]+(alphavpar[17]+alphaz[16]+alphax[2])*f[19]+(alphavpar[4]+alphaz[3]+alphax[1])*f[18]+f[4]*alphavpar[18]+alphaz[2]*f[17]+alphavpar[2]*f[16]+f[2]*alphavpar[16]+alphaz[12]*f[15]+alphavpar[12]*f[14]+f[12]*alphavpar[14]+(alphavpar[9]+alphaz[7]+alphax[0])*f[11]+f[9]*alphavpar[11]+alphaz[6]*f[10]+(alphax[7]+alphaz[0])*f[9]+f[0]*alphaz[9]+alphavpar[6]*f[8]+f[6]*alphavpar[8]+alphavpar[0]*f[7]+f[0]*alphavpar[7]+(alphax[3]+alphaz[1])*f[4]+f[1]*alphaz[4]+alphavpar[1]*f[3]+f[1]*alphavpar[3]); 
  out[19] += 0.3061862178478971*(alphay[12]*f[29]+alphaz[12]*f[28]+alphavpar[12]*f[27]+(alphavpar[9]+alphaz[7]+alphay[6])*f[26]+alphay[5]*f[25]+alphaz[5]*f[24]+alphavpar[5]*f[22]+alphavpar[20]*f[21]+f[20]*alphavpar[21]+(alphavpar[4]+alphaz[3]+alphay[2])*f[19]+(alphavpar[17]+alphaz[16]+alphay[1])*f[18]+f[17]*(alphavpar[18]+alphay[16]+alphaz[1])+alphavpar[1]*f[16]+f[1]*alphavpar[16]+alphavpar[13]*f[14]+f[13]*alphavpar[14]+(alphavpar[10]+alphaz[8]+alphay[0])*f[11]+f[10]*(alphavpar[11]+alphay[8]+alphaz[0])+(alphay[7]+alphaz[6])*f[9]+f[6]*alphaz[9]+f[7]*alphay[9]+alphavpar[0]*f[8]+f[0]*alphavpar[8]+alphavpar[6]*f[7]+f[6]*alphavpar[7]+(alphay[3]+alphaz[2])*f[4]+f[2]*alphaz[4]+f[3]*(alphay[4]+alphavpar[2])+f[2]*alphavpar[3]); 
  out[20] += 0.3061862178478971*((alphay[8]+alphax[7])*f[27]+alphay[4]*f[23]+(alphay[16]+alphax[3])*f[22]+alphay[3]*f[21]+(alphay[2]+alphax[1])*f[20]+alphay[9]*f[15]+alphay[7]*f[14]+(alphay[6]+alphax[0])*f[13]+(alphax[6]+alphay[0])*f[12]+f[0]*alphay[12]+(alphax[2]+alphay[1])*f[5]+f[1]*alphay[5]); 
  out[21] += 0.3061862178478971*((alphaz[8]+alphax[6])*f[27]+alphaz[4]*f[23]+(alphaz[16]+alphax[2])*f[22]+(alphaz[3]+alphax[1])*f[21]+alphaz[2]*f[20]+alphaz[9]*f[15]+(alphaz[7]+alphax[0])*f[14]+alphaz[6]*f[13]+(alphax[7]+alphaz[0])*f[12]+f[0]*alphaz[12]+(alphax[3]+alphaz[1])*f[5]+f[1]*alphaz[5]); 
  out[22] += 0.3061862178478971*(alphay[9]*f[29]+alphaz[9]*f[28]+(alphaz[7]+alphay[6])*f[27]+alphay[4]*f[25]+alphaz[4]*f[24]+(alphaz[3]+alphay[2])*f[22]+(alphaz[16]+alphay[1])*f[21]+(alphay[16]+alphaz[1])*f[20]+(alphaz[8]+alphay[0])*f[14]+(alphay[8]+alphaz[0])*f[13]+(alphay[7]+alphaz[6])*f[12]+f[6]*alphaz[12]+f[7]*alphay[12]+(alphay[3]+alphaz[2])*f[5]+f[2]*alphaz[5]+f[3]*alphay[5]); 
  out[23] += 0.3061862178478971*((alphavpar[11]+alphax[7])*f[29]+(alphavpar[10]+alphax[6])*f[28]+alphavpar[8]*f[27]+(alphavpar[18]+alphax[3])*f[25]+(alphavpar[17]+alphax[2])*f[24]+(alphavpar[4]+alphax[1])*f[23]+alphavpar[16]*f[22]+alphavpar[3]*f[21]+f[3]*alphavpar[21]+alphavpar[2]*f[20]+f[2]*alphavpar[20]+(alphavpar[9]+alphax[0])*f[15]+alphavpar[7]*f[14]+f[7]*alphavpar[14]+alphavpar[6]*f[13]+f[6]*alphavpar[13]+alphavpar[0]*f[12]+f[0]*alphavpar[12]+alphavpar[1]*f[5]+f[1]*alphavpar[5]); 
  out[24] += 0.3061862178478971*((alphavpar[18]+alphay[16])*f[31]+(alphavpar[11]+alphay[8])*f[30]+alphay[7]*f[29]+(alphavpar[9]+alphay[6])*f[28]+alphavpar[7]*f[27]+alphay[3]*f[25]+(alphavpar[4]+alphay[2])*f[24]+(alphavpar[17]+alphay[1])*f[23]+alphavpar[3]*f[22]+alphavpar[16]*f[21]+f[16]*alphavpar[21]+alphavpar[1]*f[20]+f[1]*alphavpar[20]+(alphavpar[10]+alphay[0])*f[15]+alphavpar[8]*f[14]+f[8]*alphavpar[14]+alphavpar[0]*f[13]+f[0]*alphavpar[13]+(alphay[9]+alphavpar[6])*f[12]+f[9]*alphay[12]+f[6]*alphavpar[12]+(alphay[4]+alphavpar[2])*f[5]+f[4]*alphay[5]+f[2]*alphavpar[5]); 
  out[25] += 0.3061862178478971*((alphavpar[17]+alphaz[16])*f[31]+(alphavpar[10]+alphaz[8])*f[30]+(alphavpar[9]+alphaz[7])*f[29]+alphaz[6]*f[28]+alphavpar[6]*f[27]+(alphavpar[4]+alphaz[3])*f[25]+alphaz[2]*f[24]+(alphavpar[18]+alphaz[1])*f[23]+alphavpar[2]*f[22]+alphavpar[1]*f[21]+f[1]*alphavpar[21]+alphavpar[16]*f[20]+f[16]*alphavpar[20]+(alphavpar[11]+alphaz[0])*f[15]+alphavpar[0]*f[14]+f[0]*alphavpar[14]+alphavpar[8]*f[13]+f[8]*alphavpar[13]+(alphaz[9]+alphavpar[7])*f[12]+f[9]*alphaz[12]+f[7]*alphavpar[12]+(alphaz[4]+alphavpar[3])*f[5]+f[4]*alphaz[5]+f[3]*alphavpar[5]); 
  out[26] += 0.3061862178478971*(alphay[5]*f[29]+alphaz[5]*f[28]+alphavpar[5]*f[27]+(alphavpar[4]+alphaz[3]+alphay[2]+alphax[1])*f[26]+alphay[12]*f[25]+alphaz[12]*f[24]+alphavpar[12]*f[22]+alphavpar[13]*f[21]+f[13]*alphavpar[21]+alphavpar[14]*f[20]+f[14]*alphavpar[20]+(alphavpar[9]+alphaz[7]+alphay[6]+alphax[0])*f[19]+(alphavpar[10]+alphaz[8]+alphax[6]+alphay[0])*f[18]+f[10]*alphavpar[18]+(alphavpar[11]+alphay[8]+alphax[7]+alphaz[0])*f[17]+f[11]*alphavpar[17]+alphavpar[0]*f[16]+f[11]*alphaz[16]+f[10]*alphay[16]+f[0]*alphavpar[16]+(alphax[2]+alphay[1])*f[11]+(alphax[3]+alphaz[1])*f[10]+(alphay[3]+alphaz[2])*f[9]+f[2]*alphaz[9]+f[3]*alphay[9]+alphavpar[1]*f[8]+f[1]*alphavpar[8]+(alphay[4]+alphavpar[2])*f[7]+f[4]*alphay[7]+f[2]*alphavpar[7]+(alphaz[4]+alphavpar[3])*f[6]+f[4]*alphaz[6]+f[3]*alphavpar[6]); 
  out[27] += 0.3061862178478971*(alphay[4]*f[29]+alphaz[4]*f[28]+(alphaz[3]+alphay[2]+alphax[1])*f[27]+alphay[9]*f[25]+alphaz[9]*f[24]+(alphaz[7]+alphay[6]+alphax[0])*f[22]+(alphaz[8]+alphax[6]+alphay[0])*f[21]+(alphay[8]+alphax[7]+alphaz[0])*f[20]+f[14]*alphaz[16]+f[13]*alphay[16]+(alphax[2]+alphay[1])*f[14]+(alphax[3]+alphaz[1])*f[13]+(alphay[3]+alphaz[2])*f[12]+f[2]*alphaz[12]+f[3]*alphay[12]+alphay[5]*f[7]+f[5]*alphay[7]+alphaz[5]*f[6]+f[5]*alphaz[6]); 
  out[28] += 0.3061862178478971*((alphavpar[11]+alphay[8]+alphax[7])*f[31]+(alphavpar[18]+alphay[16]+alphax[3])*f[30]+alphay[3]*f[29]+(alphavpar[4]+alphay[2]+alphax[1])*f[28]+alphavpar[3]*f[27]+alphay[7]*f[25]+(alphavpar[9]+alphay[6]+alphax[0])*f[24]+(alphavpar[10]+alphax[6]+alphay[0])*f[23]+alphavpar[7]*f[22]+alphavpar[8]*f[21]+f[8]*alphavpar[21]+alphavpar[0]*f[20]+f[0]*alphavpar[20]+f[15]*alphavpar[17]+alphavpar[14]*f[16]+f[14]*alphavpar[16]+(alphax[2]+alphay[1])*f[15]+alphavpar[1]*f[13]+f[1]*alphavpar[13]+(alphay[4]+alphavpar[2])*f[12]+f[4]*alphay[12]+f[2]*alphavpar[12]+alphay[5]*f[9]+f[5]*alphay[9]+alphavpar[5]*f[6]+f[5]*alphavpar[6]); 
  out[29] += 0.3061862178478971*((alphavpar[10]+alphaz[8]+alphax[6])*f[31]+(alphavpar[17]+alphaz[16]+alphax[2])*f[30]+(alphavpar[4]+alphaz[3]+alphax[1])*f[29]+alphaz[2]*f[28]+alphavpar[2]*f[27]+(alphavpar[9]+alphaz[7]+alphax[0])*f[25]+alphaz[6]*f[24]+(alphavpar[11]+alphax[7]+alphaz[0])*f[23]+alphavpar[6]*f[22]+alphavpar[0]*f[21]+f[0]*alphavpar[21]+alphavpar[8]*f[20]+f[8]*alphavpar[20]+f[15]*alphavpar[18]+alphavpar[13]*f[16]+f[13]*alphavpar[16]+(alphax[3]+alphaz[1])*f[15]+alphavpar[1]*f[14]+f[1]*alphavpar[14]+(alphaz[4]+alphavpar[3])*f[12]+f[4]*alphaz[12]+f[3]*alphavpar[12]+alphaz[5]*f[9]+f[5]*alphaz[9]+alphavpar[5]*f[7]+f[5]*alphavpar[7]); 
  out[30] += 0.3061862178478971*((alphavpar[9]+alphaz[7]+alphay[6])*f[31]+(alphavpar[4]+alphaz[3]+alphay[2])*f[30]+(alphavpar[17]+alphaz[16]+alphay[1])*f[29]+(alphavpar[18]+alphay[16]+alphaz[1])*f[28]+alphavpar[1]*f[27]+(alphavpar[10]+alphaz[8]+alphay[0])*f[25]+(alphavpar[11]+alphay[8]+alphaz[0])*f[24]+(alphay[7]+alphaz[6])*f[23]+alphavpar[0]*f[22]+(alphay[9]+alphavpar[6])*f[21]+f[6]*alphavpar[21]+(alphaz[9]+alphavpar[7])*f[20]+f[7]*alphavpar[20]+alphay[12]*f[18]+alphaz[12]*f[17]+alphavpar[12]*f[16]+f[12]*alphavpar[16]+(alphay[3]+alphaz[2])*f[15]+(alphay[4]+alphavpar[2])*f[14]+f[2]*alphavpar[14]+(alphaz[4]+alphavpar[3])*f[13]+f[3]*alphavpar[13]+alphay[5]*f[11]+alphaz[5]*f[10]+alphavpar[5]*f[8]+f[5]*alphavpar[8]); 
  out[31] += 0.3061862178478971*((alphavpar[4]+alphaz[3]+alphay[2]+alphax[1])*f[31]+(alphavpar[9]+alphaz[7]+alphay[6]+alphax[0])*f[30]+(alphavpar[10]+alphaz[8]+alphax[6]+alphay[0])*f[29]+(alphavpar[11]+alphay[8]+alphax[7]+alphaz[0])*f[28]+alphavpar[0]*f[27]+(alphavpar[17]+alphaz[16]+alphax[2]+alphay[1])*f[25]+(alphavpar[18]+alphay[16]+alphax[3]+alphaz[1])*f[24]+(alphay[3]+alphaz[2])*f[23]+alphavpar[1]*f[22]+(alphay[4]+alphavpar[2])*f[21]+f[2]*alphavpar[21]+(alphaz[4]+alphavpar[3])*f[20]+f[3]*alphavpar[20]+alphay[5]*f[18]+alphaz[5]*f[17]+alphavpar[5]*f[16]+f[5]*alphavpar[16]+(alphay[7]+alphaz[6])*f[15]+(alphay[9]+alphavpar[6])*f[14]+f[6]*alphavpar[14]+(alphaz[9]+alphavpar[7])*f[13]+f[7]*alphavpar[13]+alphavpar[8]*f[12]+f[10]*alphaz[12]+f[11]*alphay[12]+f[8]*alphavpar[12]); 
  return cflFreq; 
} 
double EmGyrokineticGenGeoVol3x2vSerP1_Bvarsz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f, double *out) 
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
  double wz = w[2];
  double rdz2 = 2.0/dxv[2];
  double wvpar = w[3];
  double rdvpar2 = 2.0/dxv[3];
  double wmu = w[4];
  double rdmu2 = 2.0/dxv[4];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wzSq = w[2]*w[2];
  double rdz2Sq = rdz2*rdz2;
  double wvparSq = w[3]*w[3];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[4]*w[4];
  double rdmu2Sq = rdmu2*rdmu2;

  double hamil[32]; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(4.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+4.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = 2.0*(bmag[3]*wmu+phi[3]*q_); 
  hamil[4] = (3.265986323710906*m_*wvpar)/rdvpar2; 
  hamil[5] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[6] = 2.0*phi[4]*q_; 
  hamil[7] = 2.0*phi[5]*q_; 
  hamil[8] = 2.0*phi[6]*q_; 
  hamil[14] = (1.154700538379252*bmag[3])/rdmu2; 
  hamil[16] = 2.0*phi[7]*q_; 

  double BstarXdBmag[32]; 
  BstarXdBmag[0] = -(0.3061862178478971*(4.0*jacobTotInv[0]*b_y[3]*m_*rdz2*wvpar+1.414213562373095*q_*((2.0*Apar[3]*b_y[3]*jacobTotInv[3]+jacobTotInv[0]*(Apar[0]*b_y[3]+b_y[0]*Apar[3]))*rdz2-1.0*((b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[6]+Apar[2]*(b_z[3]*jacobTotInv[3]+b_z[0]*jacobTotInv[0]))*rdy2)))/q_; 
  BstarXdBmag[1] = -0.4330127018922193*(((2.0*b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0])*Apar[5]+jacobTotInv[0]*Apar[1]*b_y[3])*rdz2-1.0*((b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[7]+(b_z[3]*jacobTotInv[3]+b_z[0]*jacobTotInv[0])*Apar[4])*rdy2); 
  BstarXdBmag[2] = -0.4330127018922193*((2.0*b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0])*Apar[6]+jacobTotInv[0]*Apar[2]*b_y[3])*rdz2; 
  BstarXdBmag[3] = -(0.06123724356957942*(20.0*b_y[3]*jacobTotInv[3]*m_*rdz2*wvpar+1.414213562373095*q_*(5.0*((Apar[0]*b_y[3]+b_y[0]*Apar[3])*jacobTotInv[3]+2.0*jacobTotInv[0]*Apar[3]*b_y[3])*rdz2-1.0*((9.0*b_z[3]*jacobTotInv[3]+5.0*b_z[0]*jacobTotInv[0])*Apar[6]+5.0*Apar[2]*(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3]))*rdy2)))/q_; 
  BstarXdBmag[4] = -(0.7071067811865475*jacobTotInv[0]*b_y[3]*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[6] = -0.4330127018922193*((2.0*b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0])*Apar[7]+jacobTotInv[0]*b_y[3]*Apar[4])*rdz2; 
  BstarXdBmag[7] = -0.08660254037844387*(5.0*((b_y[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_y[3])*Apar[5]+Apar[1]*b_y[3]*jacobTotInv[3])*rdz2-1.0*((9.0*b_z[3]*jacobTotInv[3]+5.0*b_z[0]*jacobTotInv[0])*Apar[7]+5.0*(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[4])*rdy2); 
  BstarXdBmag[8] = -0.4330127018922193*((b_y[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_y[3])*Apar[6]+Apar[2]*b_y[3]*jacobTotInv[3])*rdz2; 
  BstarXdBmag[11] = -(0.7071067811865475*b_y[3]*jacobTotInv[3]*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[16] = -0.4330127018922193*((b_y[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_y[3])*Apar[7]+b_y[3]*jacobTotInv[3]*Apar[4])*rdz2; 

  double BstarYdBmag[32]; 
  BstarYdBmag[0] = (0.3061862178478971*(4.0*jacobTotInv[0]*b_x[3]*m_*rdz2*wvpar+1.414213562373095*q_*((2.0*Apar[3]*b_x[3]*jacobTotInv[3]+jacobTotInv[0]*(Apar[0]*b_x[3]+b_x[0]*Apar[3]))*rdz2-1.0*((b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[5]+Apar[1]*(b_z[3]*jacobTotInv[3]+b_z[0]*jacobTotInv[0]))*rdx2)))/q_; 
  BstarYdBmag[1] = 0.4330127018922193*((2.0*b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0])*Apar[5]+jacobTotInv[0]*Apar[1]*b_x[3])*rdz2; 
  BstarYdBmag[2] = 0.4330127018922193*(((2.0*b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0])*Apar[6]+jacobTotInv[0]*Apar[2]*b_x[3])*rdz2-1.0*((b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[7]+(b_z[3]*jacobTotInv[3]+b_z[0]*jacobTotInv[0])*Apar[4])*rdx2); 
  BstarYdBmag[3] = (0.06123724356957942*(20.0*b_x[3]*jacobTotInv[3]*m_*rdz2*wvpar+1.414213562373095*q_*(5.0*((Apar[0]*b_x[3]+b_x[0]*Apar[3])*jacobTotInv[3]+2.0*jacobTotInv[0]*Apar[3]*b_x[3])*rdz2-1.0*((9.0*b_z[3]*jacobTotInv[3]+5.0*b_z[0]*jacobTotInv[0])*Apar[5]+5.0*Apar[1]*(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3]))*rdx2)))/q_; 
  BstarYdBmag[4] = (0.7071067811865475*jacobTotInv[0]*b_x[3]*m_*rdz2)/(q_*rdvpar2); 
  BstarYdBmag[6] = 0.4330127018922193*((2.0*b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0])*Apar[7]+jacobTotInv[0]*b_x[3]*Apar[4])*rdz2; 
  BstarYdBmag[7] = 0.4330127018922193*((b_x[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_x[3])*Apar[5]+Apar[1]*b_x[3]*jacobTotInv[3])*rdz2; 
  BstarYdBmag[8] = 0.08660254037844387*(5.0*((b_x[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_x[3])*Apar[6]+Apar[2]*b_x[3]*jacobTotInv[3])*rdz2-1.0*((9.0*b_z[3]*jacobTotInv[3]+5.0*b_z[0]*jacobTotInv[0])*Apar[7]+5.0*(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[4])*rdx2); 
  BstarYdBmag[11] = (0.7071067811865475*b_x[3]*jacobTotInv[3]*m_*rdz2)/(q_*rdvpar2); 
  BstarYdBmag[16] = 0.4330127018922193*((b_x[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_x[3])*Apar[7]+b_x[3]*jacobTotInv[3]*Apar[4])*rdz2; 

  double BstarZdBmag[32]; 
  BstarZdBmag[0] = -0.1767766952966368*(2.449489742783178*((b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*Apar[6]+Apar[2]*(b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0]))*rdy2-1.0*(2.449489742783178*((b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*Apar[5]+Apar[1]*(b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0]))*rdx2+4.0*(cmag[3]*jacobTotInv[3]+cmag[0]*jacobTotInv[0]))); 
  BstarZdBmag[1] = -0.4330127018922193*((b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*Apar[7]+(b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0])*Apar[4])*rdy2; 
  BstarZdBmag[2] = 0.4330127018922193*((b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*Apar[7]+(b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0])*Apar[4])*rdx2; 
  BstarZdBmag[3] = -0.03535533905932736*(2.449489742783178*((9.0*b_x[3]*jacobTotInv[3]+5.0*b_x[0]*jacobTotInv[0])*Apar[6]+5.0*Apar[2]*(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3]))*rdy2-1.0*(2.449489742783178*((9.0*b_y[3]*jacobTotInv[3]+5.0*b_y[0]*jacobTotInv[0])*Apar[5]+5.0*Apar[1]*(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3]))*rdx2+20.0*(cmag[0]*jacobTotInv[3]+jacobTotInv[0]*cmag[3]))); 
  BstarZdBmag[7] = -0.08660254037844387*((9.0*b_x[3]*jacobTotInv[3]+5.0*b_x[0]*jacobTotInv[0])*Apar[7]+5.0*(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*Apar[4])*rdy2; 
  BstarZdBmag[8] = 0.08660254037844387*((9.0*b_y[3]*jacobTotInv[3]+5.0*b_y[0]*jacobTotInv[0])*Apar[7]+5.0*(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*Apar[4])*rdx2; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[32]; 
  alphax[0] = 1.732050807568877*rdx2*((0.125*hamil[3]*(b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0])*rdz2-0.125*((b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[8]+hamil[2]*(b_z[3]*jacobTotInv[3]+b_z[0]*jacobTotInv[0]))*rdy2)/q_+(0.1767766952966368*BstarXdBmag[0]*hamil[4]*rdvpar2)/m_); 
  alphax[1] = 1.732050807568877*rdx2*((0.125*(b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0])*hamil[7]*rdz2-0.125*((b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[16]+(b_z[3]*jacobTotInv[3]+b_z[0]*jacobTotInv[0])*hamil[6])*rdy2)/q_+(0.1767766952966368*BstarXdBmag[1]*hamil[4]*rdvpar2)/m_); 
  alphax[2] = 1.732050807568877*rdx2*((0.125*(b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0])*hamil[8]*rdz2)/q_+(0.1767766952966368*BstarXdBmag[2]*hamil[4]*rdvpar2)/m_); 
  alphax[3] = 1.732050807568877*rdx2*((0.125*hamil[3]*(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*rdz2+((-0.125*(b_z[0]*jacobTotInv[0]*hamil[8]+hamil[2]*(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])))-0.225*b_z[3]*jacobTotInv[3]*hamil[8])*rdy2)/q_+(0.1767766952966368*BstarXdBmag[3]*hamil[4]*rdvpar2)/m_); 
  alphax[4] = (0.3061862178478971*BstarXdBmag[4]*hamil[4]*rdvpar2*rdx2)/m_; 
  alphax[5] = (0.2165063509461096*(b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0])*hamil[14]*rdx2*rdz2)/q_; 
  alphax[6] = 1.732050807568877*rdx2*((0.125*(b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0])*hamil[16]*rdz2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[6]*rdvpar2)/m_); 
  alphax[7] = 1.732050807568877*rdx2*((0.125*(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[7]*rdz2+((-0.125*(b_z[0]*jacobTotInv[0]*hamil[16]+(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[6]))-0.225*b_z[3]*jacobTotInv[3]*hamil[16])*rdy2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[7]*rdvpar2)/m_); 
  alphax[8] = 1.732050807568877*rdx2*((0.125*(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[8]*rdz2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[8]*rdvpar2)/m_); 
  alphax[11] = (0.3061862178478971*hamil[4]*BstarXdBmag[11]*rdvpar2*rdx2)/m_; 
  alphax[14] = (0.2165063509461096*(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[14]*rdx2*rdz2)/q_; 
  alphax[16] = 1.732050807568877*rdx2*((0.125*(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[16]*rdz2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[16]*rdvpar2)/m_); 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphax[1]-1.414213562373095*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*(alphax[14]+alphax[11]+alphax[8])+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[5]+alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*(alphax[14]+alphax[11])-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[5]+alphax[4]+alphax[3])+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*(alphax[14]+alphax[11]+alphax[8])-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[5]+alphax[4])+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*(alphax[14]+alphax[11])+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[5]+alphax[4])+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[14]-0.1767766952966368*alphax[11]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[5]+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]-0.1767766952966368*(alphax[11]+alphax[8])+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[5]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]+0.1767766952966368*alphax[11]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[5]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[14]+0.1767766952966368*(alphax[11]+alphax[8])-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[5]+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[14]+0.1767766952966368*(alphax[11]+alphax[8])+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[5]-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]+0.1767766952966368*alphax[11]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[5]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]-0.1767766952966368*(alphax[11]+alphax[8])-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[5]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[14]-0.1767766952966368*alphax[11]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[5]-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*(alphax[14]+alphax[11])+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[5]+alphax[4])-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*(alphax[14]+alphax[11]+alphax[8])+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[5]+alphax[4])-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*(alphax[14]+alphax[11])-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[5]+alphax[4]+alphax[3])-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*(alphax[14]+alphax[11]+alphax[8])-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[5]+alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*(alphax[14]+alphax[11]+alphax[8])-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[5]+alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*(alphax[14]+alphax[11])-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[5]+alphax[4]+alphax[3])+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*(alphax[14]+alphax[11]+alphax[8])+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[5]+alphax[4])+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*(alphax[14]+alphax[11])+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[5]+alphax[4])+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]-0.1767766952966368*alphax[11]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[5]+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[14]-0.1767766952966368*(alphax[11]+alphax[8])-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[5]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[14]+0.1767766952966368*alphax[11]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[5]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]+0.1767766952966368*(alphax[11]+alphax[8])+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[5]+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]+0.1767766952966368*(alphax[11]+alphax[8])-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[5]-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[14]+0.1767766952966368*alphax[11]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[5]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[14]-0.1767766952966368*(alphax[11]+alphax[8])+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[5]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]-0.1767766952966368*alphax[11]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[5]-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*(alphax[14]+alphax[11])+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[5]+alphax[4])-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*(alphax[14]+alphax[11]+alphax[8])-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[5]+alphax[4])-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*(alphax[14]+alphax[11])-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[5]+alphax[4]+alphax[3])-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*(alphax[14]+alphax[11]+alphax[8])+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[5]+alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[32]; 
  alphay[0] = 1.732050807568877*rdy2*((0.125*((b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[7]+hamil[1]*(b_z[3]*jacobTotInv[3]+b_z[0]*jacobTotInv[0]))*rdx2-0.125*hamil[3]*(b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0])*rdz2)/q_+(0.1767766952966368*BstarYdBmag[0]*hamil[4]*rdvpar2)/m_); 
  alphay[1] = 1.732050807568877*rdy2*((0.1767766952966368*BstarYdBmag[1]*hamil[4]*rdvpar2)/m_-(0.125*(b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0])*hamil[7]*rdz2)/q_); 
  alphay[2] = 1.732050807568877*rdy2*((0.125*((b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[16]+(b_z[3]*jacobTotInv[3]+b_z[0]*jacobTotInv[0])*hamil[6])*rdx2-0.125*(b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0])*hamil[8]*rdz2)/q_+(0.1767766952966368*BstarYdBmag[2]*hamil[4]*rdvpar2)/m_); 
  alphay[3] = 1.732050807568877*rdy2*(((0.125*(b_z[0]*jacobTotInv[0]*hamil[7]+hamil[1]*(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3]))+0.225*b_z[3]*jacobTotInv[3]*hamil[7])*rdx2-0.125*hamil[3]*(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*rdz2)/q_+(0.1767766952966368*BstarYdBmag[3]*hamil[4]*rdvpar2)/m_); 
  alphay[4] = (0.3061862178478971*BstarYdBmag[4]*hamil[4]*rdvpar2*rdy2)/m_; 
  alphay[5] = -(0.2165063509461096*(b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0])*hamil[14]*rdy2*rdz2)/q_; 
  alphay[6] = 1.732050807568877*rdy2*((0.1767766952966368*hamil[4]*BstarYdBmag[6]*rdvpar2)/m_-(0.125*(b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0])*hamil[16]*rdz2)/q_); 
  alphay[7] = 1.732050807568877*rdy2*((0.1767766952966368*hamil[4]*BstarYdBmag[7]*rdvpar2)/m_-(0.125*(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[7]*rdz2)/q_); 
  alphay[8] = 1.732050807568877*rdy2*(((0.125*(b_z[0]*jacobTotInv[0]*hamil[16]+(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[6])+0.225*b_z[3]*jacobTotInv[3]*hamil[16])*rdx2-0.125*(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[8]*rdz2)/q_+(0.1767766952966368*hamil[4]*BstarYdBmag[8]*rdvpar2)/m_); 
  alphay[11] = (0.3061862178478971*hamil[4]*BstarYdBmag[11]*rdvpar2*rdy2)/m_; 
  alphay[14] = -(0.2165063509461096*(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[14]*rdy2*rdz2)/q_; 
  alphay[16] = 1.732050807568877*rdy2*((0.1767766952966368*hamil[4]*BstarYdBmag[16]*rdvpar2)/m_-(0.125*(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[16]*rdz2)/q_); 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphay[2]-1.414213562373095*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphay[2]+1.414213562373095*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*(alphay[14]+alphay[11])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[11])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[11])-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*(alphay[14]+alphay[11])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[11]+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[11]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[11]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[11]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[11]+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[11]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[11]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[11]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*(alphay[14]+alphay[11])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[11])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[11])-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*(alphay[14]+alphay[11])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[11])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*(alphay[14]+alphay[11])-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*(alphay[14]+alphay[11])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[11])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[11]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[11]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[11]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[11]+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[11]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[11]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[11]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[11]+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[11])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])-0.1767766952966368*(alphay[14]+alphay[11])-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.1767766952966368*(alphay[14]+alphay[11])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[11])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphaz[32]; 
  alphaz[0] = 1.732050807568877*((0.125*((b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[8]+hamil[2]*(b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0]))*rdy2-0.125*((b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[7]+hamil[1]*(b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0]))*rdx2)/q_+(0.1767766952966368*BstarZdBmag[0]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[1] = 1.732050807568877*((0.125*((b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[16]+(b_x[3]*jacobTotInv[3]+b_x[0]*jacobTotInv[0])*hamil[6])*rdy2)/q_+(0.1767766952966368*BstarZdBmag[1]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[2] = 1.732050807568877*((0.1767766952966368*BstarZdBmag[2]*hamil[4]*rdvpar2)/m_-(0.125*((b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[16]+(b_y[3]*jacobTotInv[3]+b_y[0]*jacobTotInv[0])*hamil[6])*rdx2)/q_)*rdz2; 
  alphaz[3] = 1.732050807568877*(((0.125*(b_x[0]*jacobTotInv[0]*hamil[8]+hamil[2]*(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3]))+0.225*b_x[3]*jacobTotInv[3]*hamil[8])*rdy2+((-0.125*(b_y[0]*jacobTotInv[0]*hamil[7]+hamil[1]*(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])))-0.225*b_y[3]*jacobTotInv[3]*hamil[7])*rdx2)/q_+(0.1767766952966368*BstarZdBmag[3]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[7] = 1.732050807568877*(((0.125*(b_x[0]*jacobTotInv[0]*hamil[16]+(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[6])+0.225*b_x[3]*jacobTotInv[3]*hamil[16])*rdy2)/q_+(0.1767766952966368*hamil[4]*BstarZdBmag[7]*rdvpar2)/m_)*rdz2; 
  alphaz[8] = 1.732050807568877*((((-0.125*(b_y[0]*jacobTotInv[0]*hamil[16]+(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[6]))-0.225*b_y[3]*jacobTotInv[3]*hamil[16])*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarZdBmag[8]*rdvpar2)/m_)*rdz2; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphaz[3]-1.414213562373095*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphaz[3]+1.414213562373095*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*alphaz[7]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])-0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*alphaz[7]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])-0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*alphaz[7]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])-0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*alphaz[7]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])-0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[8]+alphaz[7]))+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[8]+alphaz[7]))+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[8]+alphaz[7]))+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[8]+alphaz[7]))+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[8])+0.3061862178478971*(alphaz[7]+alphaz[3])-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[8]+alphaz[7]+alphaz[3])+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[32]; 
  alphavpar[0] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[2]*hamil[8]+BstarZdBmag[1]*hamil[7]+BstarZdBmag[0]*hamil[3])*rdz2+(BstarYdBmag[7]*hamil[16]+BstarYdBmag[3]*hamil[8]+BstarYdBmag[1]*hamil[6]+BstarYdBmag[0]*hamil[2])*rdy2+(BstarXdBmag[8]*hamil[16]+BstarXdBmag[3]*hamil[7]+BstarXdBmag[2]*hamil[6]+BstarXdBmag[0]*hamil[1])*rdx2))/m_; 
  alphavpar[1] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[2]*hamil[16]+BstarZdBmag[0]*hamil[7]+BstarZdBmag[1]*hamil[3])*rdz2+(BstarYdBmag[3]*hamil[16]+BstarYdBmag[7]*hamil[8]+BstarYdBmag[0]*hamil[6]+BstarYdBmag[1]*hamil[2])*rdy2+(BstarXdBmag[16]*hamil[16]+BstarXdBmag[7]*hamil[7]+BstarXdBmag[6]*hamil[6]+BstarXdBmag[1]*hamil[1])*rdx2))/m_; 
  alphavpar[2] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[1]*hamil[16]+BstarZdBmag[0]*hamil[8]+BstarZdBmag[2]*hamil[3])*rdz2+(BstarYdBmag[16]*hamil[16]+BstarYdBmag[8]*hamil[8]+BstarYdBmag[6]*hamil[6]+BstarYdBmag[2]*hamil[2])*rdy2+(BstarXdBmag[3]*hamil[16]+hamil[7]*BstarXdBmag[8]+BstarXdBmag[0]*hamil[6]+hamil[1]*BstarXdBmag[2])*rdx2))/m_; 
  alphavpar[3] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[8]*hamil[8]+BstarZdBmag[7]*hamil[7]+BstarZdBmag[3]*hamil[3])*rdz2+(BstarYdBmag[1]*hamil[16]+BstarYdBmag[0]*hamil[8]+hamil[6]*BstarYdBmag[7]+hamil[2]*BstarYdBmag[3])*rdy2+(BstarXdBmag[2]*hamil[16]+hamil[6]*BstarXdBmag[8]+BstarXdBmag[0]*hamil[7]+hamil[1]*BstarXdBmag[3])*rdx2))/m_; 
  alphavpar[4] = -(0.3061862178478971*rdvpar2*((hamil[8]*BstarYdBmag[11]+hamil[2]*BstarYdBmag[4])*rdy2+(hamil[7]*BstarXdBmag[11]+hamil[1]*BstarXdBmag[4])*rdx2))/m_; 
  alphavpar[5] = -(0.3061862178478971*BstarZdBmag[0]*hamil[14]*rdvpar2*rdz2)/m_; 
  alphavpar[6] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[0]*hamil[16]+BstarZdBmag[1]*hamil[8]+BstarZdBmag[2]*hamil[7])*rdz2+(BstarYdBmag[8]*hamil[16]+hamil[8]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[6]+hamil[2]*BstarYdBmag[6])*rdy2+(BstarXdBmag[7]*hamil[16]+hamil[7]*BstarXdBmag[16]+BstarXdBmag[1]*hamil[6]+hamil[1]*BstarXdBmag[6])*rdx2))/m_; 
  alphavpar[7] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[8]*hamil[16]+BstarZdBmag[3]*hamil[7]+hamil[3]*BstarZdBmag[7])*rdz2+(BstarYdBmag[0]*hamil[16]+BstarYdBmag[1]*hamil[8]+hamil[2]*BstarYdBmag[7]+BstarYdBmag[3]*hamil[6])*rdy2+(BstarXdBmag[6]*hamil[16]+hamil[6]*BstarXdBmag[16]+BstarXdBmag[1]*hamil[7]+hamil[1]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[8] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[7]*hamil[16]+BstarZdBmag[3]*hamil[8]+hamil[3]*BstarZdBmag[8])*rdz2+(BstarYdBmag[6]*hamil[16]+hamil[6]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[8]+hamil[2]*BstarYdBmag[8])*rdy2+(BstarXdBmag[0]*hamil[16]+hamil[1]*BstarXdBmag[8]+BstarXdBmag[2]*hamil[7]+BstarXdBmag[3]*hamil[6])*rdx2))/m_; 
  alphavpar[9] = -(0.3061862178478971*(BstarYdBmag[11]*hamil[16]+BstarYdBmag[4]*hamil[6])*rdvpar2*rdy2)/m_; 
  alphavpar[10] = -(0.3061862178478971*(BstarXdBmag[11]*hamil[16]+BstarXdBmag[4]*hamil[6])*rdvpar2*rdx2)/m_; 
  alphavpar[11] = -(0.3061862178478971*rdvpar2*((hamil[2]*BstarYdBmag[11]+BstarYdBmag[4]*hamil[8])*rdy2+(hamil[1]*BstarXdBmag[11]+BstarXdBmag[4]*hamil[7])*rdx2))/m_; 
  alphavpar[12] = -(0.3061862178478971*BstarZdBmag[1]*hamil[14]*rdvpar2*rdz2)/m_; 
  alphavpar[13] = -(0.3061862178478971*BstarZdBmag[2]*hamil[14]*rdvpar2*rdz2)/m_; 
  alphavpar[14] = -(0.3061862178478971*BstarZdBmag[3]*hamil[14]*rdvpar2*rdz2)/m_; 
  alphavpar[16] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[3]*hamil[16]+BstarZdBmag[7]*hamil[8]+hamil[7]*BstarZdBmag[8])*rdz2+(BstarYdBmag[2]*hamil[16]+hamil[2]*BstarYdBmag[16]+BstarYdBmag[6]*hamil[8]+hamil[6]*BstarYdBmag[8])*rdy2+(BstarXdBmag[1]*hamil[16]+hamil[1]*BstarXdBmag[16]+BstarXdBmag[6]*hamil[7]+hamil[6]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[18] = -(0.3061862178478971*(BstarYdBmag[4]*hamil[16]+hamil[6]*BstarYdBmag[11])*rdvpar2*rdy2)/m_; 
  alphavpar[19] = -(0.3061862178478971*(BstarXdBmag[4]*hamil[16]+hamil[6]*BstarXdBmag[11])*rdvpar2*rdx2)/m_; 
  alphavpar[21] = -(0.3061862178478971*BstarZdBmag[7]*hamil[14]*rdvpar2*rdz2)/m_; 
  alphavpar[22] = -(0.3061862178478971*BstarZdBmag[8]*hamil[14]*rdvpar2*rdz2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphavpar[4]-1.414213562373095*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphavpar[4]+1.414213562373095*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[22]+alphavpar[21]))-0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*alphavpar[16]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[22])+0.1767766952966368*alphavpar[21]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]+0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.1767766952966368*(alphavpar[16]+alphavpar[14])-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*alphavpar[16]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]-0.1767766952966368*(alphavpar[16]+alphavpar[14])+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[22])+0.1767766952966368*alphavpar[21]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[22]+alphavpar[21]))-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[22])+0.1767766952966368*alphavpar[21]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[22]+alphavpar[21]))+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*(alphavpar[16]+alphavpar[14])+0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphavpar[22]+alphavpar[21]))+0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*(alphavpar[16]+alphavpar[14])-0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[22])+0.1767766952966368*alphavpar[21]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[22]+alphavpar[21]))+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*alphavpar[16]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[22])+0.1767766952966368*alphavpar[21]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]+0.1767766952966368*(alphavpar[16]+alphavpar[14])-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*alphavpar[16]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.1767766952966368*(alphavpar[16]+alphavpar[14])+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[22])+0.1767766952966368*alphavpar[21]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]-0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[22]+alphavpar[21]))+0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[22])+0.1767766952966368*alphavpar[21]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[22]+alphavpar[21]))-0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*(alphavpar[16]+alphavpar[14])+0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphavpar[22]+alphavpar[21]))-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*(alphavpar[16]+alphavpar[14])-0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[22])+0.1767766952966368*alphavpar[21]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*(alphavpar[16]+alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.3061862178478971*(alphax[16]*f[16]+alphax[14]*f[14]+alphax[11]*f[11]+alphax[8]*f[8]+alphax[7]*f[7]+alphax[6]*f[6]+alphax[5]*f[5]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.3061862178478971*(alphay[16]*f[16]+alphay[14]*f[14]+alphay[11]*f[11]+alphay[8]*f[8]+alphay[7]*f[7]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.3061862178478971*(alphaz[8]*f[8]+alphaz[7]*f[7]+alphaz[3]*f[3]+alphaz[2]*f[2]+alphaz[1]*f[1]+alphaz[0]*f[0]); 
  out[4] += 0.3061862178478971*(alphavpar[22]*f[22]+alphavpar[21]*f[21]+alphavpar[19]*f[19]+alphavpar[18]*f[18]+alphavpar[16]*f[16]+alphavpar[14]*f[14]+alphavpar[13]*f[13]+alphavpar[12]*f[12]+alphavpar[11]*f[11]+alphavpar[10]*f[10]+alphavpar[9]*f[9]+alphavpar[8]*f[8]+alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[6] += 0.3061862178478971*(alphax[14]*f[22]+alphay[14]*f[21]+alphax[11]*f[19]+alphay[11]*f[18]+(alphay[8]+alphax[7])*f[16]+f[8]*alphay[16]+f[7]*alphax[16]+alphax[5]*f[13]+alphay[5]*f[12]+alphax[4]*f[10]+alphay[4]*f[9]+alphax[3]*f[8]+f[3]*alphax[8]+alphay[3]*f[7]+f[3]*alphay[7]+(alphay[2]+alphax[1])*f[6]+f[2]*alphay[6]+f[1]*alphax[6]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[7] += 0.3061862178478971*((alphaz[8]+alphax[6])*f[16]+f[6]*alphax[16]+alphax[5]*f[14]+f[5]*alphax[14]+alphax[4]*f[11]+f[4]*alphax[11]+alphax[2]*f[8]+f[2]*alphax[8]+(alphaz[3]+alphax[1])*f[7]+f[3]*alphaz[7]+f[1]*alphax[7]+alphaz[2]*f[6]+alphax[0]*f[3]+f[0]*alphax[3]+alphaz[0]*f[1]+f[0]*alphaz[1]); 
  out[8] += 0.3061862178478971*((alphaz[7]+alphay[6])*f[16]+f[6]*alphay[16]+alphay[5]*f[14]+f[5]*alphay[14]+alphay[4]*f[11]+f[4]*alphay[11]+(alphaz[3]+alphay[2])*f[8]+f[3]*alphaz[8]+f[2]*alphay[8]+alphay[1]*f[7]+f[1]*alphay[7]+alphaz[1]*f[6]+alphay[0]*f[3]+f[0]*alphay[3]+alphaz[0]*f[2]+f[0]*alphaz[2]); 
  out[9] += 0.3061862178478971*(alphavpar[22]*f[27]+(alphavpar[19]+alphax[16])*f[26]+alphax[14]*f[25]+alphavpar[14]*f[21]+f[14]*alphavpar[21]+alphavpar[13]*f[20]+alphax[8]*f[19]+(alphavpar[11]+alphax[7])*f[18]+f[11]*alphavpar[18]+(alphavpar[10]+alphax[6])*f[17]+alphavpar[8]*f[16]+f[8]*alphavpar[16]+alphax[5]*f[15]+alphavpar[5]*f[12]+f[5]*alphavpar[12]+alphax[3]*f[11]+f[3]*alphax[11]+alphax[2]*f[10]+(alphavpar[4]+alphax[1])*f[9]+f[4]*alphavpar[9]+alphavpar[3]*f[7]+f[3]*alphavpar[7]+alphavpar[2]*f[6]+f[2]*alphavpar[6]+alphax[0]*f[4]+f[0]*alphax[4]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[10] += 0.3061862178478971*(alphavpar[21]*f[27]+(alphavpar[18]+alphay[16])*f[26]+alphay[14]*f[25]+alphavpar[14]*f[22]+f[14]*alphavpar[22]+alphavpar[12]*f[20]+(alphavpar[11]+alphay[8])*f[19]+f[11]*alphavpar[19]+alphay[7]*f[18]+(alphavpar[9]+alphay[6])*f[17]+alphavpar[7]*f[16]+f[7]*alphavpar[16]+alphay[5]*f[15]+alphavpar[5]*f[13]+f[5]*alphavpar[13]+alphay[3]*f[11]+f[3]*alphay[11]+(alphavpar[4]+alphay[2])*f[10]+f[4]*alphavpar[10]+alphay[1]*f[9]+alphavpar[3]*f[8]+f[3]*alphavpar[8]+alphavpar[1]*f[6]+f[1]*alphavpar[6]+alphay[0]*f[4]+f[0]*alphay[4]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[11] += 0.3061862178478971*(alphavpar[13]*f[22]+f[13]*alphavpar[22]+alphavpar[12]*f[21]+f[12]*alphavpar[21]+(alphavpar[10]+alphaz[8])*f[19]+f[10]*alphavpar[19]+(alphavpar[9]+alphaz[7])*f[18]+f[9]*alphavpar[18]+alphavpar[6]*f[16]+f[6]*alphavpar[16]+alphavpar[5]*f[14]+f[5]*alphavpar[14]+(alphavpar[4]+alphaz[3])*f[11]+f[4]*alphavpar[11]+alphaz[2]*f[10]+alphaz[1]*f[9]+alphavpar[2]*f[8]+f[2]*alphavpar[8]+alphavpar[1]*f[7]+f[1]*alphavpar[7]+alphaz[0]*f[4]+alphavpar[0]*f[3]+f[0]*alphavpar[3]); 
  out[12] += 0.3061862178478971*(alphax[16]*f[27]+alphax[11]*f[25]+alphax[8]*f[22]+alphax[7]*f[21]+alphax[6]*f[20]+alphax[4]*f[15]+alphax[3]*f[14]+f[3]*alphax[14]+alphax[2]*f[13]+alphax[1]*f[12]+alphax[0]*f[5]+f[0]*alphax[5]); 
  out[13] += 0.3061862178478971*(alphay[16]*f[27]+alphay[11]*f[25]+alphay[8]*f[22]+alphay[7]*f[21]+alphay[6]*f[20]+alphay[4]*f[15]+alphay[3]*f[14]+f[3]*alphay[14]+alphay[2]*f[13]+alphay[1]*f[12]+alphay[0]*f[5]+f[0]*alphay[5]); 
  out[14] += 0.3061862178478971*(alphaz[8]*f[22]+alphaz[7]*f[21]+alphaz[3]*f[14]+alphaz[2]*f[13]+alphaz[1]*f[12]+alphaz[0]*f[5]); 
  out[15] += 0.3061862178478971*(alphavpar[19]*f[30]+alphavpar[18]*f[29]+alphavpar[16]*f[27]+alphavpar[11]*f[25]+alphavpar[10]*f[24]+alphavpar[9]*f[23]+alphavpar[8]*f[22]+f[8]*alphavpar[22]+alphavpar[7]*f[21]+f[7]*alphavpar[21]+alphavpar[6]*f[20]+alphavpar[4]*f[15]+alphavpar[3]*f[14]+f[3]*alphavpar[14]+alphavpar[2]*f[13]+f[2]*alphavpar[13]+alphavpar[1]*f[12]+f[1]*alphavpar[12]+alphavpar[0]*f[5]+f[0]*alphavpar[5]); 
  out[16] += 0.3061862178478971*(alphax[5]*f[22]+alphay[5]*f[21]+alphax[4]*f[19]+alphay[4]*f[18]+(alphaz[3]+alphay[2]+alphax[1])*f[16]+f[2]*alphay[16]+f[1]*alphax[16]+f[12]*alphay[14]+f[13]*alphax[14]+f[9]*alphay[11]+f[10]*alphax[11]+(alphaz[7]+alphay[6]+alphax[0])*f[8]+f[7]*alphaz[8]+f[6]*alphay[8]+f[0]*alphax[8]+(alphax[6]+alphay[0])*f[7]+f[0]*alphay[7]+f[6]*(alphax[7]+alphaz[0])+(alphax[2]+alphay[1])*f[3]+f[1]*alphay[3]+f[2]*(alphax[3]+alphaz[1])+f[1]*alphaz[2]); 
  out[17] += 0.3061862178478971*(alphax[14]*f[30]+alphay[14]*f[29]+alphavpar[14]*f[27]+(alphavpar[11]+alphay[8]+alphax[7])*f[26]+alphax[5]*f[24]+alphay[5]*f[23]+alphavpar[21]*f[22]+f[21]*alphavpar[22]+alphavpar[5]*f[20]+(alphavpar[18]+alphay[16]+alphax[3])*f[19]+f[18]*(alphavpar[19]+alphax[16]+alphay[3])+(alphavpar[4]+alphay[2]+alphax[1])*f[17]+alphavpar[3]*f[16]+f[3]*alphavpar[16]+alphavpar[12]*f[13]+f[12]*alphavpar[13]+(alphax[8]+alphay[7])*f[11]+f[7]*alphay[11]+f[8]*alphax[11]+(alphavpar[9]+alphay[6]+alphax[0])*f[10]+f[9]*(alphavpar[10]+alphax[6]+alphay[0])+alphavpar[7]*f[8]+f[7]*alphavpar[8]+alphavpar[0]*f[6]+f[0]*alphavpar[6]+(alphax[2]+alphay[1])*f[4]+f[1]*alphay[4]+f[2]*(alphax[4]+alphavpar[1])+f[1]*alphavpar[2]); 
  out[18] += 0.3061862178478971*(alphavpar[13]*f[27]+(alphavpar[10]+alphaz[8]+alphax[6])*f[26]+alphax[5]*f[25]+f[20]*alphavpar[22]+alphavpar[5]*f[21]+f[5]*alphavpar[21]+alphax[2]*f[19]+f[17]*alphavpar[19]+(alphavpar[4]+alphaz[3]+alphax[1])*f[18]+f[4]*alphavpar[18]+(alphax[16]+alphaz[2])*f[17]+alphavpar[2]*f[16]+f[2]*alphavpar[16]+alphax[14]*f[15]+alphavpar[12]*f[14]+f[12]*alphavpar[14]+(alphavpar[9]+alphaz[7]+alphax[0])*f[11]+f[0]*alphax[11]+f[9]*alphavpar[11]+alphax[8]*f[10]+(alphax[7]+alphaz[0])*f[9]+alphavpar[6]*f[8]+f[6]*alphavpar[8]+alphavpar[0]*f[7]+f[0]*alphavpar[7]+(alphax[3]+alphaz[1])*f[4]+f[3]*(alphax[4]+alphavpar[1])+f[1]*alphavpar[3]); 
  out[19] += 0.3061862178478971*(alphavpar[12]*f[27]+(alphavpar[9]+alphaz[7]+alphay[6])*f[26]+alphay[5]*f[25]+alphavpar[5]*f[22]+f[5]*alphavpar[22]+f[20]*alphavpar[21]+(alphavpar[4]+alphaz[3]+alphay[2])*f[19]+f[4]*alphavpar[19]+alphay[1]*f[18]+f[17]*(alphavpar[18]+alphay[16]+alphaz[1])+alphavpar[1]*f[16]+f[1]*alphavpar[16]+alphay[14]*f[15]+alphavpar[13]*f[14]+f[13]*alphavpar[14]+(alphavpar[10]+alphaz[8]+alphay[0])*f[11]+f[0]*alphay[11]+f[10]*(alphavpar[11]+alphay[8]+alphaz[0])+alphay[7]*f[9]+alphavpar[0]*f[8]+f[0]*alphavpar[8]+alphavpar[6]*f[7]+f[6]*alphavpar[7]+(alphay[3]+alphaz[2])*f[4]+f[3]*(alphay[4]+alphavpar[2])+f[2]*alphavpar[3]); 
  out[20] += 0.3061862178478971*(alphax[11]*f[30]+alphay[11]*f[29]+(alphay[8]+alphax[7])*f[27]+alphax[4]*f[24]+alphay[4]*f[23]+(alphay[16]+alphax[3])*f[22]+(alphax[16]+alphay[3])*f[21]+(alphay[2]+alphax[1])*f[20]+(alphax[8]+alphay[7])*f[14]+f[7]*alphay[14]+f[8]*alphax[14]+(alphay[6]+alphax[0])*f[13]+(alphax[6]+alphay[0])*f[12]+(alphax[2]+alphay[1])*f[5]+f[1]*alphay[5]+f[2]*alphax[5]); 
  out[21] += 0.3061862178478971*((alphaz[8]+alphax[6])*f[27]+alphax[4]*f[25]+alphax[2]*f[22]+(alphaz[3]+alphax[1])*f[21]+(alphax[16]+alphaz[2])*f[20]+alphax[11]*f[15]+(alphaz[7]+alphax[0])*f[14]+f[0]*alphax[14]+alphax[8]*f[13]+(alphax[7]+alphaz[0])*f[12]+(alphax[3]+alphaz[1])*f[5]+f[3]*alphax[5]); 
  out[22] += 0.3061862178478971*((alphaz[7]+alphay[6])*f[27]+alphay[4]*f[25]+(alphaz[3]+alphay[2])*f[22]+alphay[1]*f[21]+(alphay[16]+alphaz[1])*f[20]+alphay[11]*f[15]+(alphaz[8]+alphay[0])*f[14]+f[0]*alphay[14]+(alphay[8]+alphaz[0])*f[13]+alphay[7]*f[12]+(alphay[3]+alphaz[2])*f[5]+f[3]*alphay[5]); 
  out[23] += 0.3061862178478971*((alphavpar[19]+alphax[16])*f[31]+alphax[8]*f[30]+(alphavpar[11]+alphax[7])*f[29]+(alphavpar[10]+alphax[6])*f[28]+alphavpar[8]*f[27]+(alphavpar[18]+alphax[3])*f[25]+alphax[2]*f[24]+(alphavpar[4]+alphax[1])*f[23]+alphavpar[16]*f[22]+f[16]*alphavpar[22]+alphavpar[3]*f[21]+f[3]*alphavpar[21]+alphavpar[2]*f[20]+(alphavpar[9]+alphax[0])*f[15]+(alphax[11]+alphavpar[7])*f[14]+f[11]*alphax[14]+f[7]*alphavpar[14]+alphavpar[6]*f[13]+f[6]*alphavpar[13]+alphavpar[0]*f[12]+f[0]*alphavpar[12]+(alphax[4]+alphavpar[1])*f[5]+f[4]*alphax[5]+f[1]*alphavpar[5]); 
  out[24] += 0.3061862178478971*((alphavpar[18]+alphay[16])*f[31]+(alphavpar[11]+alphay[8])*f[30]+alphay[7]*f[29]+(alphavpar[9]+alphay[6])*f[28]+alphavpar[7]*f[27]+(alphavpar[19]+alphay[3])*f[25]+(alphavpar[4]+alphay[2])*f[24]+alphay[1]*f[23]+alphavpar[3]*f[22]+f[3]*alphavpar[22]+alphavpar[16]*f[21]+f[16]*alphavpar[21]+alphavpar[1]*f[20]+(alphavpar[10]+alphay[0])*f[15]+(alphay[11]+alphavpar[8])*f[14]+f[11]*alphay[14]+f[8]*alphavpar[14]+alphavpar[0]*f[13]+f[0]*alphavpar[13]+alphavpar[6]*f[12]+f[6]*alphavpar[12]+(alphay[4]+alphavpar[2])*f[5]+f[4]*alphay[5]+f[2]*alphavpar[5]); 
  out[25] += 0.3061862178478971*((alphavpar[10]+alphaz[8])*f[30]+(alphavpar[9]+alphaz[7])*f[29]+alphavpar[6]*f[27]+(alphavpar[4]+alphaz[3])*f[25]+(alphavpar[19]+alphaz[2])*f[24]+(alphavpar[18]+alphaz[1])*f[23]+alphavpar[2]*f[22]+f[2]*alphavpar[22]+alphavpar[1]*f[21]+f[1]*alphavpar[21]+alphavpar[16]*f[20]+(alphavpar[11]+alphaz[0])*f[15]+alphavpar[0]*f[14]+f[0]*alphavpar[14]+alphavpar[8]*f[13]+f[8]*alphavpar[13]+alphavpar[7]*f[12]+f[7]*alphavpar[12]+alphavpar[3]*f[5]+f[3]*alphavpar[5]); 
  out[26] += 0.3061862178478971*(alphax[5]*f[30]+alphay[5]*f[29]+alphavpar[5]*f[27]+(alphavpar[4]+alphaz[3]+alphay[2]+alphax[1])*f[26]+alphax[14]*f[24]+alphay[14]*f[23]+alphavpar[12]*f[22]+f[12]*alphavpar[22]+alphavpar[13]*f[21]+f[13]*alphavpar[21]+alphavpar[14]*f[20]+(alphavpar[9]+alphaz[7]+alphay[6]+alphax[0])*f[19]+f[9]*alphavpar[19]+(alphavpar[10]+alphaz[8]+alphax[6]+alphay[0])*f[18]+f[10]*alphavpar[18]+(alphavpar[11]+alphay[8]+alphax[7]+alphaz[0])*f[17]+alphavpar[0]*f[16]+f[10]*alphay[16]+f[9]*alphax[16]+f[0]*alphavpar[16]+(alphax[2]+alphay[1])*f[11]+f[1]*alphay[11]+f[2]*alphax[11]+(alphax[3]+alphaz[1])*f[10]+(alphay[3]+alphaz[2])*f[9]+(alphax[4]+alphavpar[1])*f[8]+f[4]*alphax[8]+f[1]*alphavpar[8]+(alphay[4]+alphavpar[2])*f[7]+f[4]*alphay[7]+f[2]*alphavpar[7]+alphavpar[3]*f[6]+f[3]*alphavpar[6]); 
  out[27] += 0.3061862178478971*(alphax[4]*f[30]+alphay[4]*f[29]+(alphaz[3]+alphay[2]+alphax[1])*f[27]+alphax[11]*f[24]+alphay[11]*f[23]+(alphaz[7]+alphay[6]+alphax[0])*f[22]+(alphaz[8]+alphax[6]+alphay[0])*f[21]+(alphay[8]+alphax[7]+alphaz[0])*f[20]+f[13]*alphay[16]+f[12]*alphax[16]+(alphax[2]+alphay[1])*f[14]+f[1]*alphay[14]+f[2]*alphax[14]+(alphax[3]+alphaz[1])*f[13]+(alphay[3]+alphaz[2])*f[12]+alphax[5]*f[8]+f[5]*alphax[8]+alphay[5]*f[7]+f[5]*alphay[7]); 
  out[28] += 0.3061862178478971*((alphavpar[11]+alphay[8]+alphax[7])*f[31]+(alphavpar[18]+alphay[16]+alphax[3])*f[30]+(alphavpar[19]+alphax[16]+alphay[3])*f[29]+(alphavpar[4]+alphay[2]+alphax[1])*f[28]+alphavpar[3]*f[27]+(alphax[8]+alphay[7])*f[25]+(alphavpar[9]+alphay[6]+alphax[0])*f[24]+(alphavpar[10]+alphax[6]+alphay[0])*f[23]+(alphax[11]+alphavpar[7])*f[22]+f[7]*alphavpar[22]+(alphay[11]+alphavpar[8])*f[21]+f[8]*alphavpar[21]+alphavpar[0]*f[20]+alphax[14]*f[19]+alphay[14]*f[18]+alphavpar[14]*f[16]+f[14]*alphavpar[16]+(alphax[2]+alphay[1])*f[15]+(alphax[4]+alphavpar[1])*f[13]+f[1]*alphavpar[13]+(alphay[4]+alphavpar[2])*f[12]+f[2]*alphavpar[12]+alphax[5]*f[10]+alphay[5]*f[9]+alphavpar[5]*f[6]+f[5]*alphavpar[6]); 
  out[29] += 0.3061862178478971*((alphavpar[10]+alphaz[8]+alphax[6])*f[31]+alphax[2]*f[30]+(alphavpar[4]+alphaz[3]+alphax[1])*f[29]+(alphavpar[19]+alphax[16]+alphaz[2])*f[28]+alphavpar[2]*f[27]+(alphavpar[9]+alphaz[7]+alphax[0])*f[25]+alphax[8]*f[24]+(alphavpar[11]+alphax[7]+alphaz[0])*f[23]+alphavpar[6]*f[22]+f[6]*alphavpar[22]+alphavpar[0]*f[21]+f[0]*alphavpar[21]+alphavpar[8]*f[20]+f[15]*alphavpar[18]+alphavpar[13]*f[16]+f[13]*alphavpar[16]+(alphax[3]+alphaz[1])*f[15]+(alphax[4]+alphavpar[1])*f[14]+f[4]*alphax[14]+f[1]*alphavpar[14]+alphavpar[3]*f[12]+f[3]*alphavpar[12]+alphax[5]*f[11]+f[5]*alphax[11]+alphavpar[5]*f[7]+f[5]*alphavpar[7]); 
  out[30] += 0.3061862178478971*((alphavpar[9]+alphaz[7]+alphay[6])*f[31]+(alphavpar[4]+alphaz[3]+alphay[2])*f[30]+alphay[1]*f[29]+(alphavpar[18]+alphay[16]+alphaz[1])*f[28]+alphavpar[1]*f[27]+(alphavpar[10]+alphaz[8]+alphay[0])*f[25]+(alphavpar[11]+alphay[8]+alphaz[0])*f[24]+alphay[7]*f[23]+alphavpar[0]*f[22]+f[0]*alphavpar[22]+alphavpar[6]*f[21]+f[6]*alphavpar[21]+alphavpar[7]*f[20]+f[15]*alphavpar[19]+alphavpar[12]*f[16]+f[12]*alphavpar[16]+(alphay[3]+alphaz[2])*f[15]+(alphay[4]+alphavpar[2])*f[14]+f[4]*alphay[14]+f[2]*alphavpar[14]+alphavpar[3]*f[13]+f[3]*alphavpar[13]+alphay[5]*f[11]+f[5]*alphay[11]+alphavpar[5]*f[8]+f[5]*alphavpar[8]); 
  out[31] += 0.3061862178478971*((alphavpar[4]+alphaz[3]+alphay[2]+alphax[1])*f[31]+(alphavpar[9]+alphaz[7]+alphay[6]+alphax[0])*f[30]+(alphavpar[10]+alphaz[8]+alphax[6]+alphay[0])*f[29]+(alphavpar[11]+alphay[8]+alphax[7]+alphaz[0])*f[28]+alphavpar[0]*f[27]+(alphax[2]+alphay[1])*f[25]+(alphavpar[18]+alphay[16]+alphax[3]+alphaz[1])*f[24]+(alphavpar[19]+alphax[16]+alphay[3]+alphaz[2])*f[23]+(alphax[4]+alphavpar[1])*f[22]+f[1]*alphavpar[22]+(alphay[4]+alphavpar[2])*f[21]+f[2]*alphavpar[21]+alphavpar[3]*f[20]+alphax[5]*f[19]+alphay[5]*f[18]+alphavpar[5]*f[16]+f[5]*alphavpar[16]+(alphax[8]+alphay[7])*f[15]+alphavpar[6]*f[14]+f[9]*alphay[14]+f[10]*alphax[14]+f[6]*alphavpar[14]+(alphax[11]+alphavpar[7])*f[13]+f[7]*alphavpar[13]+(alphay[11]+alphavpar[8])*f[12]+f[8]*alphavpar[12]); 
  return cflFreq; 
} 
double EmGyrokineticGenGeoVol3x2vSerP1_Bvarsxz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f, double *out) 
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
  double wz = w[2];
  double rdz2 = 2.0/dxv[2];
  double wvpar = w[3];
  double rdvpar2 = 2.0/dxv[3];
  double wmu = w[4];
  double rdmu2 = 2.0/dxv[4];

  double wxSq = w[0]*w[0];
  double rdx2Sq = rdx2*rdx2;
  double wySq = w[1]*w[1];
  double rdy2Sq = rdy2*rdy2;
  double wzSq = w[2]*w[2];
  double rdz2Sq = rdz2*rdz2;
  double wvparSq = w[3]*w[3];
  double rdvpar2Sq = rdvpar2*rdvpar2;
  double wmuSq = w[4]*w[4];
  double rdmu2Sq = rdmu2*rdmu2;

  double hamil[32]; 
  hamil[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(4.0*m_*wvparSq+2.828427124746191*(bmag[0]*wmu+phi[0]*q_))+4.0*m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = 2.0*(bmag[3]*wmu+phi[3]*q_); 
  hamil[4] = (3.265986323710906*m_*wvpar)/rdvpar2; 
  hamil[5] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[6] = 2.0*phi[4]*q_; 
  hamil[7] = 2.0*(bmag[5]*wmu+phi[5]*q_); 
  hamil[8] = 2.0*phi[6]*q_; 
  hamil[12] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[14] = (1.154700538379252*bmag[3])/rdmu2; 
  hamil[16] = 2.0*phi[7]*q_; 
  hamil[21] = (1.154700538379252*bmag[5])/rdmu2; 

  double BstarXdBmag[32]; 
  BstarXdBmag[0] = -(0.3061862178478971*(4.0*(jacobTotInv[1]*b_y[5]+jacobTotInv[0]*b_y[3])*m_*rdz2*wvpar+1.414213562373095*q_*((2.0*(Apar[3]*b_y[5]+b_y[3]*Apar[5])*jacobTotInv[5]+(2.0*jacobTotInv[3]*Apar[5]+Apar[0]*jacobTotInv[1]+jacobTotInv[0]*Apar[1])*b_y[5]+(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*Apar[5]+b_y[3]*(2.0*Apar[3]*jacobTotInv[3]+Apar[1]*jacobTotInv[1]+Apar[0]*jacobTotInv[0])+(b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*Apar[3])*rdz2-1.0*((b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*Apar[7]+(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5]+b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[6]+(Apar[2]*b_z[5]+b_z[3]*Apar[4])*jacobTotInv[5]+Apar[4]*(jacobTotInv[3]*b_z[5]+b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])+Apar[2]*(b_z[3]*jacobTotInv[3]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*rdy2)))/q_; 
  BstarXdBmag[1] = -(0.06123724356957942*(20.0*(jacobTotInv[0]*b_y[5]+jacobTotInv[1]*b_y[3])*m_*rdz2*wvpar+1.414213562373095*q_*((2.0*(9.0*Apar[5]*b_y[5]+5.0*Apar[3]*b_y[3])*jacobTotInv[5]+(10.0*Apar[3]*jacobTotInv[3]+9.0*Apar[1]*jacobTotInv[1]+5.0*Apar[0]*jacobTotInv[0])*b_y[5]+(10.0*b_y[3]*jacobTotInv[3]+9.0*b_y[1]*jacobTotInv[1]+5.0*b_y[0]*jacobTotInv[0])*Apar[5]+5.0*((Apar[0]*jacobTotInv[1]+jacobTotInv[0]*Apar[1])*b_y[3]+(b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*Apar[3]))*rdz2-1.0*((9.0*(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5])+5.0*(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3]))*Apar[7]+5.0*(b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*Apar[6]+(9.0*Apar[4]*b_z[5]+5.0*Apar[2]*b_z[3])*jacobTotInv[5]+5.0*Apar[2]*jacobTotInv[3]*b_z[5]+(5.0*b_z[3]*jacobTotInv[3]+9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[4]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])*rdy2)))/q_; 
  BstarXdBmag[2] = -0.4330127018922193*((2.0*(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5])+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*Apar[7]+(2.0*(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3])+b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*Apar[6]+(jacobTotInv[0]*Apar[4]+jacobTotInv[1]*Apar[2])*b_y[5]+b_y[3]*(jacobTotInv[1]*Apar[4]+jacobTotInv[0]*Apar[2]))*rdz2; 
  BstarXdBmag[3] = -(0.06123724356957942*(20.0*(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3])*m_*rdz2*wvpar+1.414213562373095*q_*(5.0*((Apar[0]*b_y[5]+b_y[0]*Apar[5]+Apar[1]*b_y[3]+b_y[1]*Apar[3])*jacobTotInv[5]+(2.0*jacobTotInv[0]*Apar[5]+Apar[1]*jacobTotInv[3]+2.0*jacobTotInv[1]*Apar[3])*b_y[5]+(b_y[1]*jacobTotInv[3]+2.0*jacobTotInv[1]*b_y[3])*Apar[5]+(Apar[0]*b_y[3]+b_y[0]*Apar[3])*jacobTotInv[3]+2.0*jacobTotInv[0]*Apar[3]*b_y[3])*rdz2-1.0*((9.0*(b_z[3]*jacobTotInv[5]+jacobTotInv[3]*b_z[5])+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[7]+(9.0*(b_z[5]*jacobTotInv[5]+b_z[3]*jacobTotInv[3])+5.0*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*Apar[6]+5.0*((b_z[0]*Apar[4]+b_z[1]*Apar[2])*jacobTotInv[5]+(jacobTotInv[0]*Apar[4]+jacobTotInv[1]*Apar[2])*b_z[5]+(b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*Apar[4]+Apar[2]*(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])))*rdy2)))/q_; 
  BstarXdBmag[4] = -(0.7071067811865475*(jacobTotInv[1]*b_y[5]+jacobTotInv[0]*b_y[3])*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[6] = -0.08660254037844387*((2.0*(9.0*b_y[5]*jacobTotInv[5]+5.0*b_y[3]*jacobTotInv[3])+9.0*b_y[1]*jacobTotInv[1]+5.0*b_y[0]*jacobTotInv[0])*Apar[7]+5.0*(2.0*(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5])+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*Apar[6]+(9.0*jacobTotInv[1]*Apar[4]+5.0*jacobTotInv[0]*Apar[2])*b_y[5]+5.0*b_y[3]*(jacobTotInv[0]*Apar[4]+jacobTotInv[1]*Apar[2]))*rdz2; 
  BstarXdBmag[7] = -(0.01224744871391589*(100.0*(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5])*m_*rdz2*wvpar+1.414213562373095*q_*(5.0*((9.0*(Apar[1]*b_y[5]+b_y[1]*Apar[5])+5.0*(Apar[0]*b_y[3]+b_y[0]*Apar[3]))*jacobTotInv[5]+(18.0*jacobTotInv[1]*Apar[5]+5.0*(Apar[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*Apar[3]))*b_y[5]+5.0*((b_y[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_y[3])*Apar[5]+(Apar[1]*b_y[3]+b_y[1]*Apar[3])*jacobTotInv[3]+2.0*jacobTotInv[1]*Apar[3]*b_y[3]))*rdz2-1.0*((81.0*b_z[5]*jacobTotInv[5]+5.0*(9.0*(b_z[3]*jacobTotInv[3]+b_z[1]*jacobTotInv[1])+5.0*b_z[0]*jacobTotInv[0]))*Apar[7]+5.0*((9.0*(b_z[3]*jacobTotInv[5]+jacobTotInv[3]*b_z[5])+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[6]+(9.0*b_z[1]*Apar[4]+5.0*b_z[0]*Apar[2])*jacobTotInv[5]+(9.0*jacobTotInv[1]*Apar[4]+5.0*jacobTotInv[0]*Apar[2])*b_z[5]+5.0*((b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[4]+Apar[2]*(b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3]))))*rdy2)))/q_; 
  BstarXdBmag[8] = -0.4330127018922193*((b_y[0]*jacobTotInv[5]+2.0*jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+2.0*jacobTotInv[1]*b_y[3])*Apar[7]+(b_y[1]*jacobTotInv[5]+2.0*jacobTotInv[1]*b_y[5]+b_y[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_y[3])*Apar[6]+(Apar[2]*b_y[5]+b_y[3]*Apar[4])*jacobTotInv[5]+jacobTotInv[3]*(Apar[4]*b_y[5]+Apar[2]*b_y[3]))*rdz2; 
  BstarXdBmag[9] = -(0.7071067811865475*(jacobTotInv[0]*b_y[5]+jacobTotInv[1]*b_y[3])*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[11] = -(0.7071067811865475*(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3])*m_*rdz2)/(q_*rdvpar2); 
  BstarXdBmag[16] = -0.08660254037844387*((9.0*(b_y[1]*jacobTotInv[5]+2.0*jacobTotInv[1]*b_y[5])+5.0*(b_y[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_y[3]))*Apar[7]+5.0*(b_y[0]*jacobTotInv[5]+2.0*jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+2.0*jacobTotInv[1]*b_y[3])*Apar[6]+(9.0*Apar[4]*b_y[5]+5.0*Apar[2]*b_y[3])*jacobTotInv[5]+5.0*jacobTotInv[3]*(Apar[2]*b_y[5]+b_y[3]*Apar[4]))*rdz2; 
  BstarXdBmag[18] = -(0.7071067811865475*(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5])*m_*rdz2)/(q_*rdvpar2); 

  double BstarYdBmag[32]; 
  BstarYdBmag[0] = (0.3061862178478971*(4.0*m_*((jacobTotInv[1]*b_x[5]+jacobTotInv[0]*b_x[3])*rdz2-1.0*(jacobTotInv[3]*b_z[5]+jacobTotInv[0]*b_z[1])*rdx2)*wvpar+1.414213562373095*q_*((2.0*(Apar[3]*b_x[5]+b_x[3]*Apar[5])*jacobTotInv[5]+(2.0*jacobTotInv[3]*Apar[5]+Apar[0]*jacobTotInv[1]+jacobTotInv[0]*Apar[1])*b_x[5]+(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[5]+b_x[3]*(2.0*Apar[3]*jacobTotInv[3]+Apar[1]*jacobTotInv[1]+Apar[0]*jacobTotInv[0])+(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*Apar[3])*rdz2-1.0*(2.0*(Apar[1]*b_z[5]+b_z[1]*Apar[5])*jacobTotInv[5]+(2.0*jacobTotInv[1]*Apar[5]+Apar[0]*jacobTotInv[3]+jacobTotInv[0]*Apar[3])*b_z[5]+(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[5]+(Apar[1]*b_z[3]+b_z[1]*Apar[3])*jacobTotInv[3]+2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*rdx2)))/q_; 
  BstarYdBmag[1] = (0.06123724356957942*(20.0*m_*((jacobTotInv[0]*b_x[5]+jacobTotInv[1]*b_x[3])*rdz2-1.0*(b_z[5]*jacobTotInv[5]+b_z[1]*jacobTotInv[1])*rdx2)*wvpar+1.414213562373095*q_*((2.0*(9.0*Apar[5]*b_x[5]+5.0*Apar[3]*b_x[3])*jacobTotInv[5]+(10.0*Apar[3]*jacobTotInv[3]+9.0*Apar[1]*jacobTotInv[1]+5.0*Apar[0]*jacobTotInv[0])*b_x[5]+(10.0*b_x[3]*jacobTotInv[3]+9.0*b_x[1]*jacobTotInv[1]+5.0*b_x[0]*jacobTotInv[0])*Apar[5]+5.0*((Apar[0]*jacobTotInv[1]+jacobTotInv[0]*Apar[1])*b_x[3]+(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[3]))*rdz2-5.0*((Apar[0]*b_z[5]+b_z[0]*Apar[5]+Apar[1]*b_z[3]+b_z[1]*Apar[3])*jacobTotInv[5]+(2.0*(jacobTotInv[0]*Apar[5]+Apar[1]*jacobTotInv[3])+jacobTotInv[1]*Apar[3])*b_z[5]+(2.0*b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*Apar[5]+(Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])*rdx2)))/q_; 
  BstarYdBmag[2] = 0.4330127018922193*(((2.0*(b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5])+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[7]+(2.0*(b_x[5]*jacobTotInv[5]+b_x[3]*jacobTotInv[3])+b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*Apar[6]+(jacobTotInv[0]*Apar[4]+jacobTotInv[1]*Apar[2])*b_x[5]+b_x[3]*(jacobTotInv[1]*Apar[4]+jacobTotInv[0]*Apar[2]))*rdz2-1.0*((2.0*(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5])+b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[7]+(jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3])*Apar[6]+b_z[5]*(2.0*Apar[4]*jacobTotInv[5]+Apar[2]*jacobTotInv[3])+(b_z[3]*jacobTotInv[3]+2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[4]+jacobTotInv[0]*b_z[1]*Apar[2])*rdx2); 
  BstarYdBmag[3] = (0.06123724356957942*(20.0*m_*((b_x[5]*jacobTotInv[5]+b_x[3]*jacobTotInv[3])*rdz2-1.0*(jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3])*rdx2)*wvpar+1.414213562373095*q_*(5.0*((Apar[0]*b_x[5]+b_x[0]*Apar[5]+Apar[1]*b_x[3]+b_x[1]*Apar[3])*jacobTotInv[5]+(2.0*jacobTotInv[0]*Apar[5]+Apar[1]*jacobTotInv[3]+2.0*jacobTotInv[1]*Apar[3])*b_x[5]+(b_x[1]*jacobTotInv[3]+2.0*jacobTotInv[1]*b_x[3])*Apar[5]+(Apar[0]*b_x[3]+b_x[0]*Apar[3])*jacobTotInv[3]+2.0*jacobTotInv[0]*Apar[3]*b_x[3])*rdz2-1.0*(2.0*(9.0*Apar[5]*b_z[5]+5.0*Apar[1]*b_z[1])*jacobTotInv[5]+(9.0*Apar[3]*jacobTotInv[3]+5.0*(2.0*Apar[1]*jacobTotInv[1]+Apar[0]*jacobTotInv[0]))*b_z[5]+(9.0*b_z[3]*jacobTotInv[3]+5.0*(2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*Apar[5]+5.0*((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[3]+jacobTotInv[0]*(Apar[1]*b_z[3]+b_z[1]*Apar[3])))*rdx2)))/q_; 
  BstarYdBmag[4] = (0.7071067811865475*m_*((jacobTotInv[1]*b_x[5]+jacobTotInv[0]*b_x[3])*rdz2-1.0*(jacobTotInv[3]*b_z[5]+jacobTotInv[0]*b_z[1])*rdx2))/(q_*rdvpar2); 
  BstarYdBmag[6] = 0.08660254037844387*(((2.0*(9.0*b_x[5]*jacobTotInv[5]+5.0*b_x[3]*jacobTotInv[3])+9.0*b_x[1]*jacobTotInv[1]+5.0*b_x[0]*jacobTotInv[0])*Apar[7]+5.0*(2.0*(b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5])+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[6]+(9.0*jacobTotInv[1]*Apar[4]+5.0*jacobTotInv[0]*Apar[2])*b_x[5]+5.0*b_x[3]*(jacobTotInv[0]*Apar[4]+jacobTotInv[1]*Apar[2]))*rdz2-5.0*((b_z[0]*jacobTotInv[5]+2.0*(jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3])+jacobTotInv[1]*b_z[3])*Apar[7]+(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5])*Apar[6]+(Apar[2]*b_z[5]+b_z[3]*Apar[4])*jacobTotInv[5]+Apar[4]*(2.0*jacobTotInv[3]*b_z[5]+b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])+b_z[1]*jacobTotInv[1]*Apar[2])*rdx2); 
  BstarYdBmag[7] = (0.06123724356957942*(20.0*m_*((b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5])*rdz2-1.0*(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5])*rdx2)*wvpar+1.414213562373095*q_*(((9.0*(Apar[1]*b_x[5]+b_x[1]*Apar[5])+5.0*(Apar[0]*b_x[3]+b_x[0]*Apar[3]))*jacobTotInv[5]+(18.0*jacobTotInv[1]*Apar[5]+5.0*(Apar[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*Apar[3]))*b_x[5]+5.0*((b_x[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_x[3])*Apar[5]+(Apar[1]*b_x[3]+b_x[1]*Apar[3])*jacobTotInv[3]+2.0*jacobTotInv[1]*Apar[3]*b_x[3]))*rdz2-1.0*((9.0*(Apar[3]*b_z[5]+b_z[3]*Apar[5])+5.0*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*jacobTotInv[5]+(18.0*jacobTotInv[3]*Apar[5]+5.0*(Apar[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]))*b_z[5]+5.0*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[5]+2.0*Apar[1]*b_z[1]*jacobTotInv[3]+jacobTotInv[1]*(Apar[1]*b_z[3]+b_z[1]*Apar[3])))*rdx2)))/q_; 
  BstarYdBmag[8] = 0.08660254037844387*(5.0*((b_x[0]*jacobTotInv[5]+2.0*jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+2.0*jacobTotInv[1]*b_x[3])*Apar[7]+(b_x[1]*jacobTotInv[5]+2.0*jacobTotInv[1]*b_x[5]+b_x[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_x[3])*Apar[6]+(Apar[2]*b_x[5]+b_x[3]*Apar[4])*jacobTotInv[5]+jacobTotInv[3]*(Apar[4]*b_x[5]+Apar[2]*b_x[3]))*rdz2-1.0*((9.0*(2.0*b_z[5]*jacobTotInv[5]+b_z[3]*jacobTotInv[3])+5.0*(2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*Apar[7]+(9.0*jacobTotInv[3]*b_z[5]+5.0*jacobTotInv[0]*b_z[1])*Apar[6]+5.0*(2.0*b_z[1]*Apar[4]*jacobTotInv[5]+(2.0*jacobTotInv[1]*Apar[4]+jacobTotInv[0]*Apar[2])*b_z[5]+(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*Apar[4]+b_z[1]*Apar[2]*jacobTotInv[3]))*rdx2); 
  BstarYdBmag[9] = (0.7071067811865475*m_*((jacobTotInv[0]*b_x[5]+jacobTotInv[1]*b_x[3])*rdz2-1.0*(b_z[5]*jacobTotInv[5]+b_z[1]*jacobTotInv[1])*rdx2))/(q_*rdvpar2); 
  BstarYdBmag[11] = (0.7071067811865475*m_*((b_x[5]*jacobTotInv[5]+b_x[3]*jacobTotInv[3])*rdz2-1.0*(jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3])*rdx2))/(q_*rdvpar2); 
  BstarYdBmag[16] = 0.08660254037844387*(((9.0*(b_x[1]*jacobTotInv[5]+2.0*jacobTotInv[1]*b_x[5])+5.0*(b_x[0]*jacobTotInv[3]+2.0*jacobTotInv[0]*b_x[3]))*Apar[7]+5.0*(b_x[0]*jacobTotInv[5]+2.0*jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+2.0*jacobTotInv[1]*b_x[3])*Apar[6]+(9.0*Apar[4]*b_x[5]+5.0*Apar[2]*b_x[3])*jacobTotInv[5]+5.0*jacobTotInv[3]*(Apar[2]*b_x[5]+b_x[3]*Apar[4]))*rdz2-1.0*((9.0*(b_z[3]*jacobTotInv[5]+2.0*jacobTotInv[3]*b_z[5])+5.0*(b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1]))*Apar[7]+(9.0*b_z[5]*jacobTotInv[5]+5.0*b_z[1]*jacobTotInv[1])*Apar[6]+5.0*((b_z[0]*Apar[4]+b_z[1]*Apar[2])*jacobTotInv[5]+(2.0*jacobTotInv[0]*Apar[4]+jacobTotInv[1]*Apar[2])*b_z[5]+(2.0*b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*Apar[4]))*rdx2); 
  BstarYdBmag[18] = (0.7071067811865475*m_*((b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5])*rdz2-1.0*(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5])*rdx2))/(q_*rdvpar2); 

  double BstarZdBmag[32]; 
  BstarZdBmag[0] = (0.1767766952966368*(6.928203230275509*(jacobTotInv[3]*b_y[5]+jacobTotInv[0]*b_y[1])*m_*rdx2*wvpar+q_*(2.449489742783178*((2.0*(Apar[1]*b_y[5]+b_y[1]*Apar[5])*jacobTotInv[5]+(2.0*jacobTotInv[1]*Apar[5]+Apar[0]*jacobTotInv[3]+jacobTotInv[0]*Apar[3])*b_y[5]+(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*Apar[5]+(Apar[1]*b_y[3]+b_y[1]*Apar[3])*jacobTotInv[3]+2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2-1.0*((b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*Apar[7]+(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5]+b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*Apar[6]+(Apar[2]*b_x[5]+b_x[3]*Apar[4])*jacobTotInv[5]+Apar[4]*(jacobTotInv[3]*b_x[5]+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])+Apar[2]*(b_x[3]*jacobTotInv[3]+b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0]))*rdy2)+4.0*(cmag[5]*jacobTotInv[5]+cmag[3]*jacobTotInv[3]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0]))))/q_; 
  BstarZdBmag[1] = (0.03535533905932736*(34.64101615137754*(b_y[5]*jacobTotInv[5]+b_y[1]*jacobTotInv[1])*m_*rdx2*wvpar+q_*(5.0*(2.449489742783178*((Apar[0]*b_y[5]+b_y[0]*Apar[5]+Apar[1]*b_y[3]+b_y[1]*Apar[3])*jacobTotInv[5]+(2.0*(jacobTotInv[0]*Apar[5]+Apar[1]*jacobTotInv[3])+jacobTotInv[1]*Apar[3])*b_y[5]+(2.0*b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*Apar[5]+(Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2+4.0*(cmag[3]*jacobTotInv[5]+jacobTotInv[3]*cmag[5]+cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))-2.449489742783178*((9.0*(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5])+5.0*(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3]))*Apar[7]+5.0*(b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*Apar[6]+(9.0*Apar[4]*b_x[5]+5.0*Apar[2]*b_x[3])*jacobTotInv[5]+5.0*Apar[2]*jacobTotInv[3]*b_x[5]+(5.0*b_x[3]*jacobTotInv[3]+9.0*b_x[1]*jacobTotInv[1]+5.0*b_x[0]*jacobTotInv[0])*Apar[4]+5.0*(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*Apar[2])*rdy2)))/q_; 
  BstarZdBmag[2] = 0.4330127018922193*((2.0*(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5])+b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*Apar[7]+(jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3])*Apar[6]+b_y[5]*(2.0*Apar[4]*jacobTotInv[5]+Apar[2]*jacobTotInv[3])+(b_y[3]*jacobTotInv[3]+2.0*b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*Apar[4]+jacobTotInv[0]*b_y[1]*Apar[2])*rdx2; 
  BstarZdBmag[3] = (0.03535533905932736*(34.64101615137754*(jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3])*m_*rdx2*wvpar+q_*(2.449489742783178*((2.0*(9.0*Apar[5]*b_y[5]+5.0*Apar[1]*b_y[1])*jacobTotInv[5]+(9.0*Apar[3]*jacobTotInv[3]+5.0*(2.0*Apar[1]*jacobTotInv[1]+Apar[0]*jacobTotInv[0]))*b_y[5]+(9.0*b_y[3]*jacobTotInv[3]+5.0*(2.0*b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0]))*Apar[5]+5.0*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[3]+jacobTotInv[0]*(Apar[1]*b_y[3]+b_y[1]*Apar[3])))*rdx2-1.0*((9.0*(b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5])+5.0*(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1]))*Apar[7]+(9.0*(b_x[5]*jacobTotInv[5]+b_x[3]*jacobTotInv[3])+5.0*(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0]))*Apar[6]+5.0*((b_x[0]*Apar[4]+b_x[1]*Apar[2])*jacobTotInv[5]+(jacobTotInv[0]*Apar[4]+jacobTotInv[1]*Apar[2])*b_x[5]+(b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*Apar[4]+Apar[2]*(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])))*rdy2)+20.0*(cmag[1]*jacobTotInv[5]+jacobTotInv[1]*cmag[5]+cmag[0]*jacobTotInv[3]+jacobTotInv[0]*cmag[3]))))/q_; 
  BstarZdBmag[4] = (0.7071067811865475*(jacobTotInv[3]*b_y[5]+jacobTotInv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[6] = 0.4330127018922193*((b_y[0]*jacobTotInv[5]+2.0*(jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3])+jacobTotInv[1]*b_y[3])*Apar[7]+(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5])*Apar[6]+(Apar[2]*b_y[5]+b_y[3]*Apar[4])*jacobTotInv[5]+Apar[4]*(2.0*jacobTotInv[3]*b_y[5]+b_y[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_y[1])+b_y[1]*jacobTotInv[1]*Apar[2])*rdx2; 
  BstarZdBmag[7] = (0.007071067811865473*(173.2050807568877*(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5])*m_*rdx2*wvpar+q_*(5.0*(2.449489742783178*((9.0*(Apar[3]*b_y[5]+b_y[3]*Apar[5])+5.0*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*jacobTotInv[5]+(18.0*jacobTotInv[3]*Apar[5]+5.0*(Apar[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]))*b_y[5]+5.0*((b_y[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_y[1])*Apar[5]+2.0*Apar[1]*b_y[1]*jacobTotInv[3]+jacobTotInv[1]*(Apar[1]*b_y[3]+b_y[1]*Apar[3])))*rdx2+20.0*(cmag[0]*jacobTotInv[5]+jacobTotInv[0]*cmag[5]+cmag[1]*jacobTotInv[3]+jacobTotInv[1]*cmag[3]))-2.449489742783178*((81.0*b_x[5]*jacobTotInv[5]+5.0*(9.0*(b_x[3]*jacobTotInv[3]+b_x[1]*jacobTotInv[1])+5.0*b_x[0]*jacobTotInv[0]))*Apar[7]+5.0*((9.0*(b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5])+5.0*(b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1]))*Apar[6]+(9.0*b_x[1]*Apar[4]+5.0*b_x[0]*Apar[2])*jacobTotInv[5]+(9.0*jacobTotInv[1]*Apar[4]+5.0*jacobTotInv[0]*Apar[2])*b_x[5]+5.0*((b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*Apar[4]+Apar[2]*(b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3]))))*rdy2)))/q_; 
  BstarZdBmag[8] = 0.08660254037844387*((9.0*(2.0*b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3])+5.0*(2.0*b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0]))*Apar[7]+(9.0*jacobTotInv[3]*b_y[5]+5.0*jacobTotInv[0]*b_y[1])*Apar[6]+5.0*(2.0*b_y[1]*Apar[4]*jacobTotInv[5]+(2.0*jacobTotInv[1]*Apar[4]+jacobTotInv[0]*Apar[2])*b_y[5]+(b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*Apar[4]+b_y[1]*Apar[2]*jacobTotInv[3]))*rdx2; 
  BstarZdBmag[9] = (0.7071067811865475*(b_y[5]*jacobTotInv[5]+b_y[1]*jacobTotInv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[11] = (0.7071067811865475*(jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[16] = 0.08660254037844387*((9.0*(b_y[3]*jacobTotInv[5]+2.0*jacobTotInv[3]*b_y[5])+5.0*(b_y[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_y[1]))*Apar[7]+(9.0*b_y[5]*jacobTotInv[5]+5.0*b_y[1]*jacobTotInv[1])*Apar[6]+5.0*((b_y[0]*Apar[4]+b_y[1]*Apar[2])*jacobTotInv[5]+(2.0*jacobTotInv[0]*Apar[4]+jacobTotInv[1]*Apar[2])*b_y[5]+(2.0*b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*Apar[4]))*rdx2; 
  BstarZdBmag[18] = (0.7071067811865475*(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[32]; 
  alphax[0] = 1.732050807568877*rdx2*((0.125*((b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5]+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[7]+hamil[3]*(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3]+b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0]))*rdz2-0.125*((b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*hamil[16]+(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5]+b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[8]+(b_z[3]*jacobTotInv[5]+jacobTotInv[3]*b_z[5]+b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[6]+hamil[2]*(b_z[5]*jacobTotInv[5]+b_z[3]*jacobTotInv[3]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*rdy2)/q_+(0.1767766952966368*BstarXdBmag[0]*hamil[4]*rdvpar2)/m_); 
  alphax[1] = 1.732050807568877*rdx2*(((0.125*(b_y[0]*jacobTotInv[0]*hamil[7]+hamil[3]*(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5]+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1]))+(0.225*b_y[5]*jacobTotInv[5]+0.125*b_y[3]*jacobTotInv[3]+0.225*b_y[1]*jacobTotInv[1])*hamil[7])*rdz2+((-0.125*((b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[16]+(b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*hamil[8]+b_z[0]*jacobTotInv[0]*hamil[6]+hamil[2]*(b_z[3]*jacobTotInv[5]+jacobTotInv[3]*b_z[5]+b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])))-0.225*(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5])*hamil[16]+((-0.225*b_z[5]*jacobTotInv[5])-0.125*b_z[3]*jacobTotInv[3]-0.225*b_z[1]*jacobTotInv[1])*hamil[6])*rdy2)/q_+(0.1767766952966368*BstarXdBmag[1]*hamil[4]*rdvpar2)/m_); 
  alphax[2] = 1.732050807568877*rdx2*((0.125*((b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5]+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[16]+(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3]+b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[8])*rdz2)/q_+(0.1767766952966368*BstarXdBmag[2]*hamil[4]*rdvpar2)/m_); 
  alphax[3] = 1.732050807568877*rdx2*((0.125*((b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*hamil[7]+hamil[3]*(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5]+b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3]))*rdz2+((-0.225*(jacobTotInv[3]*b_z[5]*hamil[16]+(b_z[5]*jacobTotInv[5]+b_z[3]*jacobTotInv[3])*hamil[8]))-0.125*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[16]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[8]+(b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*hamil[6]+hamil[2]*(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5]+b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3]))-0.225*b_z[3]*jacobTotInv[5]*hamil[16])*rdy2)/q_+(0.1767766952966368*BstarXdBmag[3]*hamil[4]*rdvpar2)/m_); 
  alphax[4] = (0.3061862178478971*BstarXdBmag[4]*hamil[4]*rdvpar2*rdx2)/m_; 
  alphax[5] = (0.2165063509461096*((b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5]+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[21]+(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3]+b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[14])*rdx2*rdz2)/q_; 
  alphax[6] = 1.732050807568877*rdx2*(((0.125*(b_y[0]*jacobTotInv[0]*hamil[16]+(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5]+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[8])+(0.225*b_y[5]*jacobTotInv[5]+0.125*b_y[3]*jacobTotInv[3]+0.225*b_y[1]*jacobTotInv[1])*hamil[16])*rdz2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[6]*rdvpar2)/m_); 
  alphax[7] = 1.732050807568877*rdx2*(((0.125*((b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[7]+hamil[3]*(b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3]))+0.225*(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5])*hamil[7])*rdz2+((9.0*((-0.045*b_z[5]*jacobTotInv[5])-0.025*(b_z[3]*jacobTotInv[3]+b_z[1]*jacobTotInv[1]))-0.125*b_z[0]*jacobTotInv[0])*hamil[16]-0.225*(jacobTotInv[3]*b_z[5]*hamil[8]+(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5])*hamil[6])-0.125*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[8]+(b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[6]+hamil[2]*(b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3]))-0.225*b_z[3]*jacobTotInv[5]*hamil[8])*rdy2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[7]*rdvpar2)/m_); 
  alphax[8] = 1.732050807568877*rdx2*((0.125*((b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*hamil[16]+(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5]+b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[8])*rdz2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[8]*rdvpar2)/m_); 
  alphax[9] = (0.3061862178478971*hamil[4]*BstarXdBmag[9]*rdvpar2*rdx2)/m_; 
  alphax[11] = (0.3061862178478971*hamil[4]*BstarXdBmag[11]*rdvpar2*rdx2)/m_; 
  alphax[12] = (1.732050807568877*(0.125*(b_y[0]*jacobTotInv[0]*hamil[21]+(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5]+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[14])+(0.225*b_y[5]*jacobTotInv[5]+0.125*b_y[3]*jacobTotInv[3]+0.225*b_y[1]*jacobTotInv[1])*hamil[21])*rdx2*rdz2)/q_; 
  alphax[14] = (0.2165063509461096*((b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*hamil[21]+(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5]+b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[14])*rdx2*rdz2)/q_; 
  alphax[16] = 1.732050807568877*rdx2*(((0.125*((b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[16]+(b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*hamil[8])+0.225*(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5])*hamil[16])*rdz2)/q_+(0.1767766952966368*hamil[4]*BstarXdBmag[16]*rdvpar2)/m_); 
  alphax[18] = (0.3061862178478971*hamil[4]*BstarXdBmag[18]*rdvpar2*rdx2)/m_; 
  alphax[21] = (1.732050807568877*(0.125*((b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[21]+(b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*hamil[14])+0.225*(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5])*hamil[21])*rdx2*rdz2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphax[1]-1.414213562373095*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.3061862178478971*(alphax[21]+alphax[18]+alphax[16]))+0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[5]+alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[21]+alphax[18]))+0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[5]+alphax[4]+alphax[3])+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[21]+alphax[18]+alphax[16])-0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[5]+alphax[4])+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[21]+alphax[18])-0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[5]+alphax[4])+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[21])+0.3061862178478971*alphax[18]-0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[5]+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[21])+0.3061862178478971*(alphax[18]+alphax[16])+0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[5]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[21]-0.3061862178478971*alphax[18]+0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[5]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[21]-0.3061862178478971*(alphax[18]+alphax[16])-0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[5]+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[21]-0.3061862178478971*(alphax[18]+alphax[16])-0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[5]-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[21]-0.3061862178478971*alphax[18]+0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[5]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[21])+0.3061862178478971*(alphax[18]+alphax[16])+0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[5]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[21])+0.3061862178478971*alphax[18]-0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[5]-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[21]+alphax[18])-0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[5]+alphax[4])-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[21]+alphax[18]+alphax[16])-0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[5]+alphax[4])-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[21]+alphax[18]))+0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[5]+alphax[4]+alphax[3])-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[21]+alphax[18]+alphax[16]))+0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[5]+alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.3061862178478971*(alphax[21]+alphax[18]+alphax[16])+0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[5]+alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[21]+alphax[18])-0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[5]+alphax[4]+alphax[3])+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[21]+alphax[18]+alphax[16]))-0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[5]+alphax[4])+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[21]+alphax[18]))+0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[5]+alphax[4])+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[21]-0.3061862178478971*alphax[18]+0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[5]+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[21]-0.3061862178478971*(alphax[18]+alphax[16])+0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[5]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[21])+0.3061862178478971*alphax[18]-0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[5]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[21])+0.3061862178478971*(alphax[18]+alphax[16])-0.1767766952966368*alphax[14]-0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[5]+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[21])+0.3061862178478971*(alphax[18]+alphax[16])-0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[5]-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[21])+0.3061862178478971*alphax[18]-0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[5]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[21]-0.3061862178478971*(alphax[18]+alphax[16])+0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[5]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[21]-0.3061862178478971*alphax[18]+0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[5]-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[21]+alphax[18]))+0.3061862178478971*alphax[16]-0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[5]+alphax[4])-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[21]+alphax[18]+alphax[16]))-0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]-0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[5]+alphax[4])-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[21]+alphax[18])-0.3061862178478971*alphax[16]+0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[5]+alphax[4]+alphax[3])-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[21]+alphax[18]+alphax[16])+0.1767766952966368*alphax[14]+0.3061862178478971*alphax[12]+0.1767766952966368*alphax[11]+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[5]+alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[32]; 
  alphay[0] = 1.732050807568877*rdy2*((0.125*((b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5]+b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[7]+hamil[1]*(b_z[5]*jacobTotInv[5]+b_z[3]*jacobTotInv[3]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*rdx2-0.125*((b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5]+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[7]+hamil[3]*(b_x[5]*jacobTotInv[5]+b_x[3]*jacobTotInv[3]+b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0]))*rdz2)/q_+(0.1767766952966368*BstarYdBmag[0]*hamil[4]*rdvpar2)/m_); 
  alphay[1] = 1.732050807568877*rdy2*(((((-0.225*b_x[5]*jacobTotInv[5])-0.125*b_x[3]*jacobTotInv[3]-0.225*b_x[1]*jacobTotInv[1])*hamil[7]-0.125*(b_x[0]*jacobTotInv[0]*hamil[7]+hamil[3]*(b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5]+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])))*rdz2+0.125*((b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*hamil[7]+hamil[1]*(b_z[3]*jacobTotInv[5]+jacobTotInv[3]*b_z[5]+b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*rdx2)/q_+(0.1767766952966368*BstarYdBmag[1]*hamil[4]*rdvpar2)/m_); 
  alphay[2] = 1.732050807568877*rdy2*((0.125*((b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5]+b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[16]+(b_z[5]*jacobTotInv[5]+b_z[3]*jacobTotInv[3]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[6])*rdx2-0.125*((b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5]+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[16]+(b_x[5]*jacobTotInv[5]+b_x[3]*jacobTotInv[3]+b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*hamil[8])*rdz2)/q_+(0.1767766952966368*BstarYdBmag[2]*hamil[4]*rdvpar2)/m_); 
  alphay[3] = 1.732050807568877*rdy2*(((0.125*((b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[7]+hamil[1]*(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5]+b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3]))+0.225*(b_z[5]*jacobTotInv[5]+b_z[3]*jacobTotInv[3])*hamil[7])*rdx2-0.125*((b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*hamil[7]+hamil[3]*(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5]+b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3]))*rdz2)/q_+(0.1767766952966368*BstarYdBmag[3]*hamil[4]*rdvpar2)/m_); 
  alphay[4] = (0.3061862178478971*BstarYdBmag[4]*hamil[4]*rdvpar2*rdy2)/m_; 
  alphay[5] = (1.732050807568877*rdy2*(0.125*((b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5]+b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[21]+(b_z[5]*jacobTotInv[5]+b_z[3]*jacobTotInv[3]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[12])*rdx2-0.125*((b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5]+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[21]+(b_x[5]*jacobTotInv[5]+b_x[3]*jacobTotInv[3]+b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*hamil[14])*rdz2))/q_; 
  alphay[6] = 1.732050807568877*rdy2*(((((-0.225*b_x[5]*jacobTotInv[5])-0.125*b_x[3]*jacobTotInv[3]-0.225*b_x[1]*jacobTotInv[1])*hamil[16]-0.125*(b_x[0]*jacobTotInv[0]*hamil[16]+(b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5]+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[8]))*rdz2+0.125*((b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*hamil[16]+(b_z[3]*jacobTotInv[5]+jacobTotInv[3]*b_z[5]+b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[6])*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarYdBmag[6]*rdvpar2)/m_); 
  alphay[7] = 1.732050807568877*rdy2*((((-0.125*((b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[7]+hamil[3]*(b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])))-0.225*(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5])*hamil[7])*rdz2+(0.125*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[7]+hamil[1]*(b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3]))+0.225*(b_z[3]*jacobTotInv[5]+jacobTotInv[3]*b_z[5])*hamil[7])*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarYdBmag[7]*rdvpar2)/m_); 
  alphay[8] = 1.732050807568877*rdy2*(((0.125*((b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[16]+(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5]+b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[6])+0.225*(b_z[5]*jacobTotInv[5]+b_z[3]*jacobTotInv[3])*hamil[16])*rdx2-0.125*((b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*hamil[16]+(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5]+b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[8])*rdz2)/q_+(0.1767766952966368*hamil[4]*BstarYdBmag[8]*rdvpar2)/m_); 
  alphay[9] = (0.3061862178478971*hamil[4]*BstarYdBmag[9]*rdvpar2*rdy2)/m_; 
  alphay[11] = (0.3061862178478971*hamil[4]*BstarYdBmag[11]*rdvpar2*rdy2)/m_; 
  alphay[12] = (1.732050807568877*rdy2*((((-0.225*b_x[5]*jacobTotInv[5])-0.125*b_x[3]*jacobTotInv[3]-0.225*b_x[1]*jacobTotInv[1])*hamil[21]-0.125*(b_x[0]*jacobTotInv[0]*hamil[21]+(b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5]+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[14]))*rdz2+0.125*((b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*hamil[21]+(b_z[3]*jacobTotInv[5]+jacobTotInv[3]*b_z[5]+b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[12])*rdx2))/q_; 
  alphay[14] = (1.732050807568877*rdy2*((0.125*((b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[21]+(b_z[1]*jacobTotInv[5]+jacobTotInv[1]*b_z[5]+b_z[0]*jacobTotInv[3]+jacobTotInv[0]*b_z[3])*hamil[12])+0.225*(b_z[5]*jacobTotInv[5]+b_z[3]*jacobTotInv[3])*hamil[21])*rdx2-0.125*((b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*hamil[21]+(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5]+b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[14])*rdz2))/q_; 
  alphay[16] = 1.732050807568877*rdy2*((((-0.125*((b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[16]+(b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*hamil[8]))-0.225*(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5])*hamil[16])*rdz2+(0.125*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[16]+(b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*hamil[6])+0.225*(b_z[3]*jacobTotInv[5]+jacobTotInv[3]*b_z[5])*hamil[16])*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarYdBmag[16]*rdvpar2)/m_); 
  alphay[18] = (0.3061862178478971*hamil[4]*BstarYdBmag[18]*rdvpar2*rdy2)/m_; 
  alphay[21] = (1.732050807568877*rdy2*(((-0.125*((b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[21]+(b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*hamil[14]))-0.225*(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5])*hamil[21])*rdz2+(0.125*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[21]+(b_z[0]*jacobTotInv[5]+jacobTotInv[0]*b_z[5]+b_z[1]*jacobTotInv[3]+jacobTotInv[1]*b_z[3])*hamil[12])+0.225*(b_z[3]*jacobTotInv[5]+jacobTotInv[3]*b_z[5])*hamil[21])*rdx2))/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphay[2]-1.414213562373095*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphay[2]+1.414213562373095*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.1767766952966368*(alphay[21]+alphay[18]))-0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[12]+alphay[11]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphay[21]+alphay[18])+0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[11]-0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphay[21]+alphay[18])+0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[11]+0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphay[21]+alphay[18]))-0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[12]+alphay[11]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphay[21])+0.1767766952966368*alphay[18]-0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[12])-0.1767766952966368*(alphay[11]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphay[21]-0.1767766952966368*alphay[18]+0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*(alphay[12]+alphay[11])+0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphay[21]-0.1767766952966368*alphay[18]+0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*(alphay[12]+alphay[11])-0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphay[21])+0.1767766952966368*alphay[18]-0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[12])+0.1767766952966368*(alphay[11]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphay[21]-0.1767766952966368*alphay[18]-0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[12])+0.1767766952966368*(alphay[11]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphay[21])+0.1767766952966368*alphay[18]+0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*(alphay[12]+alphay[11])-0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphay[21])+0.1767766952966368*alphay[18]+0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*(alphay[12]+alphay[11])+0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphay[21]-0.1767766952966368*alphay[18]-0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[12])-0.1767766952966368*(alphay[11]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphay[21]+alphay[18])-0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[12]+alphay[11]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphay[21]+alphay[18]))+0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[11]+0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*(alphay[21]+alphay[18]))+0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[11]-0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*(alphay[21]+alphay[18])-0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[12]+alphay[11]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*((-0.1767766952966368*(alphay[21]+alphay[18]))+0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[12]+alphay[11]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphay[21]+alphay[18])-0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[11]-0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphay[21]+alphay[18])-0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[11]+0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphay[21]+alphay[18]))+0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[12]+alphay[11]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[5]+alphay[4])+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphay[21])+0.1767766952966368*alphay[18]+0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[12])-0.1767766952966368*(alphay[11]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphay[21]-0.1767766952966368*alphay[18]-0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*(alphay[12]+alphay[11])+0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphay[21]-0.1767766952966368*alphay[18]-0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*(alphay[12]+alphay[11])-0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphay[21])+0.1767766952966368*alphay[18]+0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[12])+0.1767766952966368*(alphay[11]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[5]+0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphay[21]-0.1767766952966368*alphay[18]+0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[12])+0.1767766952966368*(alphay[11]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphay[21])+0.1767766952966368*alphay[18]-0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*(alphay[12]+alphay[11])-0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphay[21])+0.1767766952966368*alphay[18]-0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*(alphay[12]+alphay[11])+0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphay[21]-0.1767766952966368*alphay[18]+0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[12])-0.1767766952966368*(alphay[11]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[5]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphay[21]+alphay[18])+0.3061862178478971*alphay[16]-0.1767766952966368*(alphay[14]+alphay[12]+alphay[11]+alphay[9])-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphay[21]+alphay[18]))-0.3061862178478971*alphay[16]-0.1767766952966368*alphay[14]+0.1767766952966368*alphay[12]-0.1767766952966368*alphay[11]+0.1767766952966368*alphay[9]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4])-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*(alphay[21]+alphay[18]))-0.3061862178478971*alphay[16]+0.1767766952966368*alphay[14]-0.1767766952966368*alphay[12]+0.1767766952966368*alphay[11]-0.1767766952966368*alphay[9]+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*(alphay[21]+alphay[18])+0.3061862178478971*alphay[16]+0.1767766952966368*(alphay[14]+alphay[12]+alphay[11]+alphay[9])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[5]+alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphaz[32]; 
  alphaz[0] = 1.732050807568877*((0.125*((b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*hamil[16]+(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5]+b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[8]+(b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5]+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[6]+hamil[2]*(b_x[5]*jacobTotInv[5]+b_x[3]*jacobTotInv[3]+b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0]))*rdy2-0.125*((b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5]+b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[7]+hamil[1]*(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3]+b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0]))*rdx2)/q_+(0.1767766952966368*BstarZdBmag[0]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[1] = 1.732050807568877*(((0.125*((b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[16]+(b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*hamil[8]+b_x[0]*jacobTotInv[0]*hamil[6]+hamil[2]*(b_x[3]*jacobTotInv[5]+jacobTotInv[3]*b_x[5]+b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1]))+0.225*(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5])*hamil[16]+(0.225*b_x[5]*jacobTotInv[5]+0.125*b_x[3]*jacobTotInv[3]+0.225*b_x[1]*jacobTotInv[1])*hamil[6])*rdy2-0.125*((b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*hamil[7]+hamil[1]*(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5]+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1]))*rdx2)/q_+(0.1767766952966368*BstarZdBmag[1]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[2] = 1.732050807568877*((0.1767766952966368*BstarZdBmag[2]*hamil[4]*rdvpar2)/m_-(0.125*((b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5]+b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[16]+(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3]+b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[6])*rdx2)/q_)*rdz2; 
  alphaz[3] = 1.732050807568877*(((0.225*(jacobTotInv[3]*b_x[5]*hamil[16]+(b_x[5]*jacobTotInv[5]+b_x[3]*jacobTotInv[3])*hamil[8])+0.125*((b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[16]+(b_x[1]*jacobTotInv[1]+b_x[0]*jacobTotInv[0])*hamil[8]+(b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3])*hamil[6]+hamil[2]*(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5]+b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3]))+0.225*b_x[3]*jacobTotInv[5]*hamil[16])*rdy2+((-0.125*((b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[7]+hamil[1]*(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5]+b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])))-0.225*(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3])*hamil[7])*rdx2)/q_+(0.1767766952966368*BstarZdBmag[3]*hamil[4]*rdvpar2)/m_)*rdz2; 
  alphaz[4] = (0.3061862178478971*BstarZdBmag[4]*hamil[4]*rdvpar2*rdz2)/m_; 
  alphaz[5] = -(0.2165063509461096*((b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5]+b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[21]+(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3]+b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[12])*rdx2*rdz2)/q_; 
  alphaz[6] = 1.732050807568877*((0.1767766952966368*hamil[4]*BstarZdBmag[6]*rdvpar2)/m_-(0.125*((b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*hamil[16]+(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5]+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[6])*rdx2)/q_)*rdz2; 
  alphaz[7] = 1.732050807568877*((((9.0*(0.045*b_x[5]*jacobTotInv[5]+0.025*(b_x[3]*jacobTotInv[3]+b_x[1]*jacobTotInv[1]))+0.125*b_x[0]*jacobTotInv[0])*hamil[16]+0.225*(jacobTotInv[3]*b_x[5]*hamil[8]+(b_x[1]*jacobTotInv[5]+jacobTotInv[1]*b_x[5])*hamil[6])+0.125*((b_x[0]*jacobTotInv[1]+jacobTotInv[0]*b_x[1])*hamil[8]+(b_x[0]*jacobTotInv[3]+jacobTotInv[0]*b_x[3])*hamil[6]+hamil[2]*(b_x[0]*jacobTotInv[5]+jacobTotInv[0]*b_x[5]+b_x[1]*jacobTotInv[3]+jacobTotInv[1]*b_x[3]))+0.225*b_x[3]*jacobTotInv[5]*hamil[8])*rdy2+((-0.125*((b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[7]+hamil[1]*(b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])))-0.225*(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5])*hamil[7])*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarZdBmag[7]*rdvpar2)/m_)*rdz2; 
  alphaz[8] = 1.732050807568877*((((-0.125*((b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[16]+(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5]+b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[6]))-0.225*(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3])*hamil[16])*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarZdBmag[8]*rdvpar2)/m_)*rdz2; 
  alphaz[9] = (0.3061862178478971*hamil[4]*BstarZdBmag[9]*rdvpar2*rdz2)/m_; 
  alphaz[11] = (0.3061862178478971*hamil[4]*BstarZdBmag[11]*rdvpar2*rdz2)/m_; 
  alphaz[12] = -(0.2165063509461096*((b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*hamil[21]+(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5]+b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[12])*rdx2*rdz2)/q_; 
  alphaz[14] = (1.732050807568877*((-0.125*((b_y[1]*jacobTotInv[1]+b_y[0]*jacobTotInv[0])*hamil[21]+(b_y[1]*jacobTotInv[5]+jacobTotInv[1]*b_y[5]+b_y[0]*jacobTotInv[3]+jacobTotInv[0]*b_y[3])*hamil[12]))-0.225*(b_y[5]*jacobTotInv[5]+b_y[3]*jacobTotInv[3])*hamil[21])*rdx2*rdz2)/q_; 
  alphaz[16] = 1.732050807568877*((((-0.125*((b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[16]+(b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*hamil[6]))-0.225*(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5])*hamil[16])*rdx2)/q_+(0.1767766952966368*hamil[4]*BstarZdBmag[16]*rdvpar2)/m_)*rdz2; 
  alphaz[18] = (0.3061862178478971*hamil[4]*BstarZdBmag[18]*rdvpar2*rdz2)/m_; 
  alphaz[21] = (1.732050807568877*((-0.125*((b_y[0]*jacobTotInv[1]+jacobTotInv[0]*b_y[1])*hamil[21]+(b_y[0]*jacobTotInv[5]+jacobTotInv[0]*b_y[5]+b_y[1]*jacobTotInv[3]+jacobTotInv[1]*b_y[3])*hamil[12]))-0.225*(b_y[3]*jacobTotInv[5]+jacobTotInv[3]*b_y[5])*hamil[21])*rdx2*rdz2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphaz[3]-1.414213562373095*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphaz[3]+1.414213562373095*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.3061862178478971*(alphaz[21]+alphaz[18]+alphaz[16]))+0.3061862178478971*alphaz[14]+0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*(alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphaz[21]+alphaz[18]+alphaz[16]+alphaz[14])-0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphaz[21]+alphaz[18]))+0.3061862178478971*(alphaz[16]+alphaz[14])+0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphaz[21]+alphaz[18])-0.3061862178478971*alphaz[16]+0.3061862178478971*alphaz[14]-0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*(alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[21])+0.3061862178478971*alphaz[18]-0.3061862178478971*alphaz[16]+0.3061862178478971*alphaz[14]+0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*alphaz[5]+0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[21]-0.3061862178478971*alphaz[18]+0.3061862178478971*(alphaz[16]+alphaz[14])-0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5])+0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[21])+0.3061862178478971*(alphaz[18]+alphaz[16]+alphaz[14])+0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5])+0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[21]-0.3061862178478971*(alphaz[18]+alphaz[16])+0.3061862178478971*alphaz[14]-0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*alphaz[5]+0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[21]-0.3061862178478971*(alphaz[18]+alphaz[16]+alphaz[14])-0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5])-0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[21])+0.3061862178478971*(alphaz[18]+alphaz[16])-0.3061862178478971*alphaz[14]+0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*alphaz[5]-0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphaz[21]-0.3061862178478971*alphaz[18]+0.3061862178478971*alphaz[16]-0.3061862178478971*alphaz[14]-0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*alphaz[5]-0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphaz[21])+0.3061862178478971*alphaz[18]-0.3061862178478971*(alphaz[16]+alphaz[14])+0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5])-0.1767766952966368*alphaz[4]-0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphaz[21]+alphaz[18])-0.3061862178478971*(alphaz[16]+alphaz[14])-0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphaz[21]+alphaz[18]))+0.3061862178478971*alphaz[16]-0.3061862178478971*alphaz[14]+0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*(alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphaz[21]+alphaz[18]+alphaz[16])-0.3061862178478971*alphaz[14]-0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*(alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphaz[21]+alphaz[18]+alphaz[16]+alphaz[14]))+0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])-0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.3061862178478971*(alphaz[21]+alphaz[18]+alphaz[16])-0.3061862178478971*alphaz[14]+0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*(alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[21]+alphaz[18]+alphaz[16]+alphaz[14]))-0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[21]+alphaz[18])-0.3061862178478971*(alphaz[16]+alphaz[14])+0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[21]+alphaz[18]))+0.3061862178478971*alphaz[16]-0.3061862178478971*alphaz[14]-0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*(alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[21]-0.3061862178478971*alphaz[18]+0.3061862178478971*alphaz[16]-0.3061862178478971*alphaz[14]+0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*alphaz[5]+0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[21])+0.3061862178478971*alphaz[18]-0.3061862178478971*(alphaz[16]+alphaz[14])-0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5])+0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[21]-0.3061862178478971*(alphaz[18]+alphaz[16]+alphaz[14])+0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*(alphaz[6]+alphaz[5])+0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[21])+0.3061862178478971*(alphaz[18]+alphaz[16])-0.3061862178478971*alphaz[14]-0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*alphaz[6]-0.1767766952966368*alphaz[5]+0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[21])+0.3061862178478971*(alphaz[18]+alphaz[16]+alphaz[14])-0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5])-0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[21]-0.3061862178478971*(alphaz[18]+alphaz[16])+0.3061862178478971*alphaz[14]+0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*alphaz[5]-0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphaz[21])+0.3061862178478971*alphaz[18]-0.3061862178478971*alphaz[16]+0.3061862178478971*alphaz[14]-0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*alphaz[5]-0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphaz[21]-0.3061862178478971*alphaz[18]+0.3061862178478971*(alphaz[16]+alphaz[14])+0.1767766952966368*alphaz[12]-0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5])-0.1767766952966368*alphaz[4]+0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[21]+alphaz[18]))+0.3061862178478971*(alphaz[16]+alphaz[14])-0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]-0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]-0.1767766952966368*(alphaz[2]+alphaz[1])+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[21]+alphaz[18])-0.3061862178478971*alphaz[16]+0.3061862178478971*alphaz[14]+0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]-0.3061862178478971*alphaz[8]+0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*(alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]-0.1767766952966368*alphaz[2]+0.1767766952966368*(alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphaz[21]+alphaz[18]+alphaz[16]))+0.3061862178478971*alphaz[14]-0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]-0.1767766952966368*alphaz[9]+0.3061862178478971*alphaz[8]-0.3061862178478971*alphaz[7]-0.1767766952966368*alphaz[6]+0.1767766952966368*(alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]+0.1767766952966368*alphaz[2]-0.1767766952966368*alphaz[1]+0.1767766952966368*alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphaz[21]+alphaz[18]+alphaz[16]+alphaz[14])+0.1767766952966368*alphaz[12]+0.3061862178478971*alphaz[11]+0.1767766952966368*alphaz[9]+0.3061862178478971*(alphaz[8]+alphaz[7])+0.1767766952966368*(alphaz[6]+alphaz[5]+alphaz[4])+0.3061862178478971*alphaz[3]+0.1767766952966368*(alphaz[2]+alphaz[1]+alphaz[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[32]; 
  alphavpar[0] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[6]*hamil[16]+BstarZdBmag[2]*hamil[8]+BstarZdBmag[1]*hamil[7]+BstarZdBmag[0]*hamil[3])*rdz2+(BstarYdBmag[7]*hamil[16]+BstarYdBmag[3]*hamil[8]+BstarYdBmag[1]*hamil[6]+BstarYdBmag[0]*hamil[2])*rdy2+(BstarXdBmag[8]*hamil[16]+BstarXdBmag[3]*hamil[7]+BstarXdBmag[2]*hamil[6]+BstarXdBmag[0]*hamil[1])*rdx2))/m_; 
  alphavpar[1] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[2]*hamil[16]+BstarZdBmag[6]*hamil[8]+BstarZdBmag[0]*hamil[7]+BstarZdBmag[1]*hamil[3])*rdz2+(BstarYdBmag[3]*hamil[16]+BstarYdBmag[7]*hamil[8]+BstarYdBmag[0]*hamil[6]+BstarYdBmag[1]*hamil[2])*rdy2+(BstarXdBmag[16]*hamil[16]+BstarXdBmag[7]*hamil[7]+BstarXdBmag[6]*hamil[6]+BstarXdBmag[1]*hamil[1])*rdx2))/m_; 
  alphavpar[2] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[1]*hamil[16]+BstarZdBmag[0]*hamil[8]+BstarZdBmag[6]*hamil[7]+BstarZdBmag[2]*hamil[3])*rdz2+(BstarYdBmag[16]*hamil[16]+BstarYdBmag[8]*hamil[8]+BstarYdBmag[6]*hamil[6]+BstarYdBmag[2]*hamil[2])*rdy2+(BstarXdBmag[3]*hamil[16]+hamil[7]*BstarXdBmag[8]+BstarXdBmag[0]*hamil[6]+hamil[1]*BstarXdBmag[2])*rdx2))/m_; 
  alphavpar[3] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[16]*hamil[16]+BstarZdBmag[8]*hamil[8]+BstarZdBmag[7]*hamil[7]+BstarZdBmag[3]*hamil[3])*rdz2+(BstarYdBmag[1]*hamil[16]+BstarYdBmag[0]*hamil[8]+hamil[6]*BstarYdBmag[7]+hamil[2]*BstarYdBmag[3])*rdy2+(BstarXdBmag[2]*hamil[16]+hamil[6]*BstarXdBmag[8]+BstarXdBmag[0]*hamil[7]+hamil[1]*BstarXdBmag[3])*rdx2))/m_; 
  alphavpar[4] = -(0.3061862178478971*rdvpar2*((hamil[7]*BstarZdBmag[9]+hamil[3]*BstarZdBmag[4])*rdz2+(hamil[16]*BstarYdBmag[18]+hamil[8]*BstarYdBmag[11]+hamil[6]*BstarYdBmag[9]+hamil[2]*BstarYdBmag[4])*rdy2+(hamil[7]*BstarXdBmag[11]+hamil[1]*BstarXdBmag[4])*rdx2))/m_; 
  alphavpar[5] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[1]*hamil[21]+BstarZdBmag[0]*hamil[14])*rdz2+(BstarXdBmag[3]*hamil[21]+BstarXdBmag[0]*hamil[12])*rdx2))/m_; 
  alphavpar[6] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[0]*hamil[16]+BstarZdBmag[1]*hamil[8]+BstarZdBmag[2]*hamil[7]+hamil[3]*BstarZdBmag[6])*rdz2+(BstarYdBmag[8]*hamil[16]+hamil[8]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[6]+hamil[2]*BstarYdBmag[6])*rdy2+(BstarXdBmag[7]*hamil[16]+hamil[7]*BstarXdBmag[16]+BstarXdBmag[1]*hamil[6]+hamil[1]*BstarXdBmag[6])*rdx2))/m_; 
  alphavpar[7] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[8]*hamil[16]+hamil[8]*BstarZdBmag[16]+BstarZdBmag[3]*hamil[7]+hamil[3]*BstarZdBmag[7])*rdz2+(BstarYdBmag[0]*hamil[16]+BstarYdBmag[1]*hamil[8]+hamil[2]*BstarYdBmag[7]+BstarYdBmag[3]*hamil[6])*rdy2+(BstarXdBmag[6]*hamil[16]+hamil[6]*BstarXdBmag[16]+BstarXdBmag[1]*hamil[7]+hamil[1]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[8] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[7]*hamil[16]+hamil[7]*BstarZdBmag[16]+BstarZdBmag[3]*hamil[8]+hamil[3]*BstarZdBmag[8])*rdz2+(BstarYdBmag[6]*hamil[16]+hamil[6]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[8]+hamil[2]*BstarYdBmag[8])*rdy2+(BstarXdBmag[0]*hamil[16]+hamil[1]*BstarXdBmag[8]+BstarXdBmag[2]*hamil[7]+BstarXdBmag[3]*hamil[6])*rdx2))/m_; 
  alphavpar[9] = -(0.3061862178478971*rdvpar2*((hamil[3]*BstarZdBmag[9]+BstarZdBmag[4]*hamil[7])*rdz2+(hamil[8]*BstarYdBmag[18]+BstarYdBmag[11]*hamil[16]+hamil[2]*BstarYdBmag[9]+BstarYdBmag[4]*hamil[6])*rdy2+(hamil[7]*BstarXdBmag[18]+hamil[1]*BstarXdBmag[9])*rdx2))/m_; 
  alphavpar[10] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[9]*hamil[16]+BstarZdBmag[4]*hamil[8])*rdz2+(BstarXdBmag[11]*hamil[16]+BstarXdBmag[4]*hamil[6])*rdx2))/m_; 
  alphavpar[11] = -(0.3061862178478971*rdvpar2*((hamil[7]*BstarZdBmag[18]+hamil[3]*BstarZdBmag[11])*rdz2+(hamil[6]*BstarYdBmag[18]+BstarYdBmag[9]*hamil[16]+hamil[2]*BstarYdBmag[11]+BstarYdBmag[4]*hamil[8])*rdy2+(hamil[1]*BstarXdBmag[11]+BstarXdBmag[4]*hamil[7])*rdx2))/m_; 
  alphavpar[12] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[0]*hamil[21]+BstarZdBmag[1]*hamil[14])*rdz2+(BstarXdBmag[7]*hamil[21]+BstarXdBmag[1]*hamil[12])*rdx2))/m_; 
  alphavpar[13] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[6]*hamil[21]+BstarZdBmag[2]*hamil[14])*rdz2+(BstarXdBmag[8]*hamil[21]+BstarXdBmag[2]*hamil[12])*rdx2))/m_; 
  alphavpar[14] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[7]*hamil[21]+BstarZdBmag[3]*hamil[14])*rdz2+(BstarXdBmag[0]*hamil[21]+BstarXdBmag[3]*hamil[12])*rdx2))/m_; 
  alphavpar[15] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[9]*hamil[21]+BstarZdBmag[4]*hamil[14])*rdz2+(BstarXdBmag[11]*hamil[21]+BstarXdBmag[4]*hamil[12])*rdx2))/m_; 
  alphavpar[16] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[3]*hamil[16]+hamil[3]*BstarZdBmag[16]+BstarZdBmag[7]*hamil[8]+hamil[7]*BstarZdBmag[8])*rdz2+(BstarYdBmag[2]*hamil[16]+hamil[2]*BstarYdBmag[16]+BstarYdBmag[6]*hamil[8]+hamil[6]*BstarYdBmag[8])*rdy2+(BstarXdBmag[1]*hamil[16]+hamil[1]*BstarXdBmag[16]+BstarXdBmag[6]*hamil[7]+hamil[6]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[17] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[4]*hamil[16]+hamil[8]*BstarZdBmag[9])*rdz2+(hamil[16]*BstarXdBmag[18]+hamil[6]*BstarXdBmag[9])*rdx2))/m_; 
  alphavpar[18] = -(0.3061862178478971*rdvpar2*((hamil[3]*BstarZdBmag[18]+hamil[7]*BstarZdBmag[11])*rdz2+(hamil[2]*BstarYdBmag[18]+BstarYdBmag[4]*hamil[16]+hamil[6]*BstarYdBmag[11]+hamil[8]*BstarYdBmag[9])*rdy2+(hamil[1]*BstarXdBmag[18]+hamil[7]*BstarXdBmag[9])*rdx2))/m_; 
  alphavpar[19] = -(0.3061862178478971*rdvpar2*((hamil[16]*BstarZdBmag[18]+hamil[8]*BstarZdBmag[11])*rdz2+(BstarXdBmag[4]*hamil[16]+hamil[6]*BstarXdBmag[11])*rdx2))/m_; 
  alphavpar[20] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[2]*hamil[21]+BstarZdBmag[6]*hamil[14])*rdz2+(BstarXdBmag[16]*hamil[21]+BstarXdBmag[6]*hamil[12])*rdx2))/m_; 
  alphavpar[21] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[3]*hamil[21]+BstarZdBmag[7]*hamil[14])*rdz2+(BstarXdBmag[1]*hamil[21]+BstarXdBmag[7]*hamil[12])*rdx2))/m_; 
  alphavpar[22] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[16]*hamil[21]+BstarZdBmag[8]*hamil[14])*rdz2+(BstarXdBmag[2]*hamil[21]+BstarXdBmag[8]*hamil[12])*rdx2))/m_; 
  alphavpar[23] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[4]*hamil[21]+BstarZdBmag[9]*hamil[14])*rdz2+(BstarXdBmag[18]*hamil[21]+BstarXdBmag[9]*hamil[12])*rdx2))/m_; 
  alphavpar[25] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[18]*hamil[21]+BstarZdBmag[11]*hamil[14])*rdz2+(BstarXdBmag[4]*hamil[21]+BstarXdBmag[11]*hamil[12])*rdx2))/m_; 
  alphavpar[26] = -(0.3061862178478971*rdvpar2*((hamil[8]*BstarZdBmag[18]+BstarZdBmag[11]*hamil[16])*rdz2+(hamil[6]*BstarXdBmag[18]+BstarXdBmag[9]*hamil[16])*rdx2))/m_; 
  alphavpar[27] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[8]*hamil[21]+hamil[14]*BstarZdBmag[16])*rdz2+(BstarXdBmag[6]*hamil[21]+hamil[12]*BstarXdBmag[16])*rdx2))/m_; 
  alphavpar[29] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[11]*hamil[21]+hamil[14]*BstarZdBmag[18])*rdz2+(BstarXdBmag[9]*hamil[21]+hamil[12]*BstarXdBmag[18])*rdx2))/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphavpar[4]-1.414213562373095*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphavpar[4]+1.414213562373095*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*(0.3061862178478971*alphavpar[29]+0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*(alphavpar[25]+alphavpar[23])-0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[29])-0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[25])+0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*alphavpar[19]+0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[29]-0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[25]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[29])+0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[25]+0.3061862178478971*alphavpar[23]+0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.1767766952966368*alphavpar[20]+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[29])-0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[25]-0.3061862178478971*alphavpar[23]+0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.1767766952966368*alphavpar[20]+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[29]+0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[25]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[29])+0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[25])-0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*alphavpar[19]+0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[29]-0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*(alphavpar[25]+alphavpar[23])-0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[29])-0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[25]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[29]+0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[25]-0.3061862178478971*alphavpar[23]+0.1767766952966368*alphavpar[22]-0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*alphavpar[19]+0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[29])+0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*(alphavpar[25]+alphavpar[23])-0.1767766952966368*alphavpar[22]+0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[29]-0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[25])-0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.1767766952966368*alphavpar[20]+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[29]+0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[25])+0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.1767766952966368*alphavpar[20]+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[29])-0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*(alphavpar[25]+alphavpar[23])-0.1767766952966368*alphavpar[22]+0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[29]-0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[25]+0.3061862178478971*alphavpar[23]+0.1767766952966368*alphavpar[22]-0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*alphavpar[19]+0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[29])+0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[25]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[29])+0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*(alphavpar[25]+alphavpar[23])-0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[29]-0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[25])-0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*alphavpar[19]-0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[29])-0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[25]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[29]+0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[25]-0.3061862178478971*alphavpar[23]+0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.1767766952966368*alphavpar[20]-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[29]-0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[25]+0.3061862178478971*alphavpar[23]+0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.1767766952966368*alphavpar[20]-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[29])+0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[25]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[29]+0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[25])+0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*alphavpar[19]-0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[29])-0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*(alphavpar[25]+alphavpar[23])-0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[29]-0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[25]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[29])+0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[25]+0.3061862178478971*alphavpar[23]+0.1767766952966368*alphavpar[22]-0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*alphavpar[19]-0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[29]+0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*(alphavpar[25]+alphavpar[23])-0.1767766952966368*alphavpar[22]+0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[29])-0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[25])+0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.1767766952966368*alphavpar[20]-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[29])+0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[25])-0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.1767766952966368*alphavpar[20]-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[29]-0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*(alphavpar[25]+alphavpar[23])-0.1767766952966368*alphavpar[22]+0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[29])-0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[25]-0.3061862178478971*alphavpar[23]+0.1767766952966368*alphavpar[22]-0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*alphavpar[19]-0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[29]+0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[25]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.3061862178478971*(alphax[21]*f[21]+alphax[18]*f[18]+alphax[16]*f[16]+alphax[14]*f[14]+alphax[12]*f[12]+alphax[11]*f[11]+alphax[9]*f[9]+alphax[8]*f[8]+alphax[7]*f[7]+alphax[6]*f[6]+alphax[5]*f[5]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.3061862178478971*(alphay[21]*f[21]+alphay[18]*f[18]+alphay[16]*f[16]+alphay[14]*f[14]+alphay[12]*f[12]+alphay[11]*f[11]+alphay[9]*f[9]+alphay[8]*f[8]+alphay[7]*f[7]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.3061862178478971*(alphaz[21]*f[21]+alphaz[18]*f[18]+alphaz[16]*f[16]+alphaz[14]*f[14]+alphaz[12]*f[12]+alphaz[11]*f[11]+alphaz[9]*f[9]+alphaz[8]*f[8]+alphaz[7]*f[7]+alphaz[6]*f[6]+alphaz[5]*f[5]+alphaz[4]*f[4]+alphaz[3]*f[3]+alphaz[2]*f[2]+alphaz[1]*f[1]+alphaz[0]*f[0]); 
  out[4] += 0.3061862178478971*(alphavpar[29]*f[29]+alphavpar[27]*f[27]+alphavpar[26]*f[26]+alphavpar[25]*f[25]+alphavpar[23]*f[23]+alphavpar[22]*f[22]+alphavpar[21]*f[21]+alphavpar[20]*f[20]+alphavpar[19]*f[19]+alphavpar[18]*f[18]+alphavpar[17]*f[17]+alphavpar[16]*f[16]+alphavpar[15]*f[15]+alphavpar[14]*f[14]+alphavpar[13]*f[13]+alphavpar[12]*f[12]+alphavpar[11]*f[11]+alphavpar[10]*f[10]+alphavpar[9]*f[9]+alphavpar[8]*f[8]+alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[6] += 0.3061862178478971*(alphax[21]*f[27]+alphax[18]*f[26]+alphax[14]*f[22]+alphay[14]*f[21]+f[14]*alphay[21]+alphax[12]*f[20]+alphax[11]*f[19]+alphay[11]*f[18]+f[11]*alphay[18]+alphax[9]*f[17]+(alphay[8]+alphax[7])*f[16]+f[8]*alphay[16]+f[7]*alphax[16]+alphax[5]*f[13]+alphay[5]*f[12]+f[5]*alphay[12]+alphax[4]*f[10]+alphay[4]*f[9]+f[4]*alphay[9]+alphax[3]*f[8]+f[3]*alphax[8]+alphay[3]*f[7]+f[3]*alphay[7]+(alphay[2]+alphax[1])*f[6]+f[2]*alphay[6]+f[1]*alphax[6]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[7] += 0.3061862178478971*((alphaz[14]+alphax[12])*f[21]+f[14]*alphaz[21]+f[12]*alphax[21]+(alphaz[11]+alphax[9])*f[18]+f[11]*alphaz[18]+f[9]*alphax[18]+(alphaz[8]+alphax[6])*f[16]+f[8]*alphaz[16]+f[6]*alphax[16]+alphax[5]*f[14]+f[5]*alphax[14]+alphaz[5]*f[12]+f[5]*alphaz[12]+alphax[4]*f[11]+f[4]*alphax[11]+alphaz[4]*f[9]+f[4]*alphaz[9]+alphax[2]*f[8]+f[2]*alphax[8]+(alphaz[3]+alphax[1])*f[7]+f[3]*alphaz[7]+f[1]*alphax[7]+alphaz[2]*f[6]+f[2]*alphaz[6]+alphax[0]*f[3]+f[0]*alphax[3]+alphaz[0]*f[1]+f[0]*alphaz[1]); 
  out[8] += 0.3061862178478971*(alphaz[21]*f[27]+alphaz[18]*f[26]+alphaz[14]*f[22]+alphay[12]*f[21]+f[12]*alphay[21]+alphaz[12]*f[20]+alphaz[11]*f[19]+alphay[9]*f[18]+f[9]*alphay[18]+alphaz[9]*f[17]+(alphaz[7]+alphay[6])*f[16]+f[7]*alphaz[16]+f[6]*alphay[16]+alphay[5]*f[14]+f[5]*alphay[14]+alphaz[5]*f[13]+alphay[4]*f[11]+f[4]*alphay[11]+alphaz[4]*f[10]+(alphaz[3]+alphay[2])*f[8]+f[3]*alphaz[8]+f[2]*alphay[8]+alphay[1]*f[7]+f[1]*alphay[7]+alphaz[1]*f[6]+f[1]*alphaz[6]+alphay[0]*f[3]+f[0]*alphay[3]+alphaz[0]*f[2]+f[0]*alphaz[2]); 
  out[9] += 0.3061862178478971*((alphavpar[25]+alphax[21])*f[29]+f[25]*alphavpar[29]+alphavpar[22]*f[27]+f[22]*alphavpar[27]+(alphavpar[19]+alphax[16])*f[26]+f[19]*alphavpar[26]+alphax[14]*f[25]+(alphavpar[15]+alphax[12])*f[23]+f[15]*alphavpar[23]+alphavpar[14]*f[21]+f[14]*alphavpar[21]+alphavpar[13]*f[20]+f[13]*alphavpar[20]+alphax[8]*f[19]+(alphavpar[11]+alphax[7])*f[18]+f[7]*alphax[18]+f[11]*alphavpar[18]+(alphavpar[10]+alphax[6])*f[17]+f[10]*alphavpar[17]+alphavpar[8]*f[16]+f[8]*alphavpar[16]+alphax[5]*f[15]+alphavpar[5]*f[12]+f[5]*alphavpar[12]+alphax[3]*f[11]+f[3]*alphax[11]+alphax[2]*f[10]+(alphavpar[4]+alphax[1])*f[9]+f[1]*alphax[9]+f[4]*alphavpar[9]+alphavpar[3]*f[7]+f[3]*alphavpar[7]+alphavpar[2]*f[6]+f[2]*alphavpar[6]+alphax[0]*f[4]+f[0]*alphax[4]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[10] += 0.3061862178478971*(alphavpar[29]*f[31]+alphavpar[25]*f[30]+alphay[21]*f[29]+alphavpar[23]*f[28]+alphavpar[21]*f[27]+f[21]*alphavpar[27]+(alphavpar[18]+alphay[16])*f[26]+f[18]*alphavpar[26]+alphay[14]*f[25]+alphavpar[15]*f[24]+alphay[12]*f[23]+alphavpar[14]*f[22]+f[14]*alphavpar[22]+alphavpar[12]*f[20]+f[12]*alphavpar[20]+(alphavpar[11]+alphay[8])*f[19]+f[11]*alphavpar[19]+alphay[7]*f[18]+f[7]*alphay[18]+(alphavpar[9]+alphay[6])*f[17]+f[9]*alphavpar[17]+alphavpar[7]*f[16]+f[7]*alphavpar[16]+alphay[5]*f[15]+alphavpar[5]*f[13]+f[5]*alphavpar[13]+alphay[3]*f[11]+f[3]*alphay[11]+(alphavpar[4]+alphay[2])*f[10]+f[4]*alphavpar[10]+alphay[1]*f[9]+f[1]*alphay[9]+alphavpar[3]*f[8]+f[3]*alphavpar[8]+alphavpar[1]*f[6]+f[1]*alphavpar[6]+alphay[0]*f[4]+f[0]*alphay[4]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[11] += 0.3061862178478971*((alphavpar[23]+alphaz[21])*f[29]+f[23]*alphavpar[29]+alphavpar[20]*f[27]+f[20]*alphavpar[27]+(alphavpar[17]+alphaz[16])*f[26]+f[17]*alphavpar[26]+(alphavpar[15]+alphaz[14])*f[25]+f[15]*alphavpar[25]+alphaz[12]*f[23]+alphavpar[13]*f[22]+f[13]*alphavpar[22]+alphavpar[12]*f[21]+f[12]*alphavpar[21]+(alphavpar[10]+alphaz[8])*f[19]+f[10]*alphavpar[19]+(alphavpar[9]+alphaz[7])*f[18]+f[7]*alphaz[18]+f[9]*alphavpar[18]+alphaz[6]*f[17]+alphavpar[6]*f[16]+f[6]*alphavpar[16]+alphaz[5]*f[15]+alphavpar[5]*f[14]+f[5]*alphavpar[14]+(alphavpar[4]+alphaz[3])*f[11]+f[3]*alphaz[11]+f[4]*alphavpar[11]+alphaz[2]*f[10]+alphaz[1]*f[9]+f[1]*alphaz[9]+alphavpar[2]*f[8]+f[2]*alphavpar[8]+alphavpar[1]*f[7]+f[1]*alphavpar[7]+alphaz[0]*f[4]+f[0]*alphaz[4]+alphavpar[0]*f[3]+f[0]*alphavpar[3]); 
  out[12] += 0.3061862178478971*(alphax[18]*f[29]+alphax[16]*f[27]+alphax[11]*f[25]+alphax[9]*f[23]+alphax[8]*f[22]+alphax[7]*f[21]+f[7]*alphax[21]+alphax[6]*f[20]+alphax[4]*f[15]+alphax[3]*f[14]+f[3]*alphax[14]+alphax[2]*f[13]+alphax[1]*f[12]+f[1]*alphax[12]+alphax[0]*f[5]+f[0]*alphax[5]); 
  out[13] += 0.3061862178478971*(alphay[18]*f[29]+alphay[16]*f[27]+alphay[11]*f[25]+alphay[9]*f[23]+alphay[8]*f[22]+alphay[7]*f[21]+f[7]*alphay[21]+alphay[6]*f[20]+alphay[4]*f[15]+alphay[3]*f[14]+f[3]*alphay[14]+alphay[2]*f[13]+alphay[1]*f[12]+f[1]*alphay[12]+alphay[0]*f[5]+f[0]*alphay[5]); 
  out[14] += 0.3061862178478971*(alphaz[18]*f[29]+alphaz[16]*f[27]+alphaz[11]*f[25]+alphaz[9]*f[23]+alphaz[8]*f[22]+alphaz[7]*f[21]+f[7]*alphaz[21]+alphaz[6]*f[20]+alphaz[4]*f[15]+alphaz[3]*f[14]+f[3]*alphaz[14]+alphaz[2]*f[13]+alphaz[1]*f[12]+f[1]*alphaz[12]+alphaz[0]*f[5]+f[0]*alphaz[5]); 
  out[15] += 0.3061862178478971*(alphavpar[26]*f[31]+alphavpar[19]*f[30]+alphavpar[18]*f[29]+f[18]*alphavpar[29]+alphavpar[17]*f[28]+alphavpar[16]*f[27]+f[16]*alphavpar[27]+alphavpar[11]*f[25]+f[11]*alphavpar[25]+alphavpar[10]*f[24]+alphavpar[9]*f[23]+f[9]*alphavpar[23]+alphavpar[8]*f[22]+f[8]*alphavpar[22]+alphavpar[7]*f[21]+f[7]*alphavpar[21]+alphavpar[6]*f[20]+f[6]*alphavpar[20]+alphavpar[4]*f[15]+f[4]*alphavpar[15]+alphavpar[3]*f[14]+f[3]*alphavpar[14]+alphavpar[2]*f[13]+f[2]*alphavpar[13]+alphavpar[1]*f[12]+f[1]*alphavpar[12]+alphavpar[0]*f[5]+f[0]*alphavpar[5]); 
  out[16] += 0.3061862178478971*((alphaz[14]+alphax[12])*f[27]+(alphaz[11]+alphax[9])*f[26]+(alphaz[21]+alphax[5])*f[22]+alphay[5]*f[21]+f[5]*alphay[21]+f[20]*(alphax[21]+alphaz[5])+(alphaz[18]+alphax[4])*f[19]+alphay[4]*f[18]+f[4]*alphay[18]+f[17]*(alphax[18]+alphaz[4])+(alphaz[3]+alphay[2]+alphax[1])*f[16]+f[3]*alphaz[16]+f[2]*alphay[16]+f[1]*alphax[16]+alphay[12]*f[14]+f[12]*alphay[14]+f[13]*(alphax[14]+alphaz[12])+alphay[9]*f[11]+f[9]*alphay[11]+f[10]*(alphax[11]+alphaz[9])+(alphaz[7]+alphay[6]+alphax[0])*f[8]+f[7]*alphaz[8]+f[6]*alphay[8]+f[0]*alphax[8]+(alphax[6]+alphay[0])*f[7]+f[0]*alphay[7]+f[6]*(alphax[7]+alphaz[0])+f[0]*alphaz[6]+(alphax[2]+alphay[1])*f[3]+f[1]*alphay[3]+f[2]*(alphax[3]+alphaz[1])+f[1]*alphaz[2]); 
  out[17] += 0.3061862178478971*((alphavpar[25]+alphax[21])*f[31]+(alphavpar[29]+alphax[14])*f[30]+alphay[14]*f[29]+(alphavpar[15]+alphax[12])*f[28]+alphavpar[14]*f[27]+f[14]*alphavpar[27]+(alphavpar[11]+alphay[8]+alphax[7])*f[26]+f[11]*alphavpar[26]+alphay[21]*f[25]+(alphavpar[23]+alphax[5])*f[24]+alphay[5]*f[23]+alphavpar[21]*f[22]+f[21]*alphavpar[22]+alphavpar[5]*f[20]+f[5]*alphavpar[20]+(alphavpar[18]+alphay[16]+alphax[3])*f[19]+f[18]*(alphavpar[19]+alphax[16]+alphay[3])+f[3]*alphay[18]+f[16]*alphax[18]+(alphavpar[4]+alphay[2]+alphax[1])*f[17]+f[4]*alphavpar[17]+alphavpar[3]*f[16]+f[3]*alphavpar[16]+alphay[12]*f[15]+alphavpar[12]*f[13]+f[12]*alphavpar[13]+(alphax[8]+alphay[7])*f[11]+f[7]*alphay[11]+f[8]*alphax[11]+(alphavpar[9]+alphay[6]+alphax[0])*f[10]+f[9]*(alphavpar[10]+alphax[6]+alphay[0])+f[0]*alphay[9]+f[6]*alphax[9]+alphavpar[7]*f[8]+f[7]*alphavpar[8]+alphavpar[0]*f[6]+f[0]*alphavpar[6]+(alphax[2]+alphay[1])*f[4]+f[1]*alphay[4]+f[2]*(alphax[4]+alphavpar[1])+f[1]*alphavpar[2]); 
  out[18] += 0.3061862178478971*((alphavpar[15]+alphaz[14]+alphax[12])*f[29]+f[15]*alphavpar[29]+alphavpar[13]*f[27]+f[13]*alphavpar[27]+(alphavpar[10]+alphaz[8]+alphax[6])*f[26]+f[10]*alphavpar[26]+(alphavpar[23]+alphaz[21]+alphax[5])*f[25]+f[23]*(alphavpar[25]+alphax[21]+alphaz[5])+alphavpar[20]*f[22]+f[20]*alphavpar[22]+alphavpar[5]*f[21]+f[5]*alphavpar[21]+(alphavpar[17]+alphaz[16]+alphax[2])*f[19]+f[17]*alphavpar[19]+(alphavpar[4]+alphaz[3]+alphax[1])*f[18]+f[3]*alphaz[18]+f[1]*alphax[18]+f[4]*alphavpar[18]+(alphax[16]+alphaz[2])*f[17]+alphavpar[2]*f[16]+f[2]*alphavpar[16]+(alphax[14]+alphaz[12])*f[15]+alphavpar[12]*f[14]+f[12]*alphavpar[14]+(alphavpar[9]+alphaz[7]+alphax[0])*f[11]+f[7]*alphaz[11]+f[0]*alphax[11]+f[9]*alphavpar[11]+(alphax[8]+alphaz[6])*f[10]+(alphax[7]+alphaz[0])*f[9]+f[0]*alphaz[9]+f[7]*alphax[9]+alphavpar[6]*f[8]+f[6]*alphavpar[8]+alphavpar[0]*f[7]+f[0]*alphavpar[7]+(alphax[3]+alphaz[1])*f[4]+f[1]*alphaz[4]+f[3]*(alphax[4]+alphavpar[1])+f[1]*alphavpar[3]); 
  out[19] += 0.3061862178478971*((alphavpar[23]+alphaz[21])*f[31]+(alphavpar[15]+alphaz[14])*f[30]+alphay[12]*f[29]+f[28]*(alphavpar[29]+alphaz[12])+alphavpar[12]*f[27]+f[12]*alphavpar[27]+(alphavpar[9]+alphaz[7]+alphay[6])*f[26]+f[9]*alphavpar[26]+alphay[5]*f[25]+f[24]*(alphavpar[25]+alphaz[5])+alphay[21]*f[23]+alphavpar[5]*f[22]+f[5]*alphavpar[22]+alphavpar[20]*f[21]+f[20]*alphavpar[21]+(alphavpar[4]+alphaz[3]+alphay[2])*f[19]+f[4]*alphavpar[19]+(alphavpar[17]+alphaz[16]+alphay[1])*f[18]+f[16]*alphaz[18]+f[1]*alphay[18]+f[17]*(alphavpar[18]+alphay[16]+alphaz[1])+alphavpar[1]*f[16]+f[1]*alphavpar[16]+alphay[14]*f[15]+alphavpar[13]*f[14]+f[13]*alphavpar[14]+(alphavpar[10]+alphaz[8]+alphay[0])*f[11]+f[8]*alphaz[11]+f[0]*alphay[11]+f[10]*(alphavpar[11]+alphay[8]+alphaz[0])+(alphay[7]+alphaz[6])*f[9]+f[6]*alphaz[9]+f[7]*alphay[9]+alphavpar[0]*f[8]+f[0]*alphavpar[8]+alphavpar[6]*f[7]+f[6]*alphavpar[7]+(alphay[3]+alphaz[2])*f[4]+f[2]*alphaz[4]+f[3]*(alphay[4]+alphavpar[2])+f[2]*alphavpar[3]); 
  out[20] += 0.3061862178478971*(alphax[18]*f[31]+alphax[11]*f[30]+alphay[11]*f[29]+alphax[9]*f[28]+(alphay[8]+alphax[7])*f[27]+alphay[18]*f[25]+alphax[4]*f[24]+alphay[4]*f[23]+(alphay[16]+alphax[3])*f[22]+(alphax[16]+alphay[3])*f[21]+f[3]*alphay[21]+f[16]*alphax[21]+(alphay[2]+alphax[1])*f[20]+alphay[9]*f[15]+(alphax[8]+alphay[7])*f[14]+f[7]*alphay[14]+f[8]*alphax[14]+(alphay[6]+alphax[0])*f[13]+(alphax[6]+alphay[0])*f[12]+f[0]*alphay[12]+f[6]*alphax[12]+(alphax[2]+alphay[1])*f[5]+f[1]*alphay[5]+f[2]*alphax[5]); 
  out[21] += 0.3061862178478971*((alphaz[11]+alphax[9])*f[29]+(alphaz[8]+alphax[6])*f[27]+(alphaz[18]+alphax[4])*f[25]+(alphax[18]+alphaz[4])*f[23]+(alphaz[16]+alphax[2])*f[22]+(alphaz[3]+alphax[1])*f[21]+f[3]*alphaz[21]+f[1]*alphax[21]+(alphax[16]+alphaz[2])*f[20]+(alphax[11]+alphaz[9])*f[15]+(alphaz[7]+alphax[0])*f[14]+f[7]*alphaz[14]+f[0]*alphax[14]+(alphax[8]+alphaz[6])*f[13]+(alphax[7]+alphaz[0])*f[12]+f[0]*alphaz[12]+f[7]*alphax[12]+(alphax[3]+alphaz[1])*f[5]+f[1]*alphaz[5]+f[3]*alphax[5]); 
  out[22] += 0.3061862178478971*(alphaz[18]*f[31]+alphaz[11]*f[30]+alphay[9]*f[29]+alphaz[9]*f[28]+(alphaz[7]+alphay[6])*f[27]+alphay[4]*f[25]+alphaz[4]*f[24]+alphay[18]*f[23]+(alphaz[3]+alphay[2])*f[22]+(alphaz[16]+alphay[1])*f[21]+f[16]*alphaz[21]+f[1]*alphay[21]+(alphay[16]+alphaz[1])*f[20]+alphay[11]*f[15]+(alphaz[8]+alphay[0])*f[14]+f[8]*alphaz[14]+f[0]*alphay[14]+(alphay[8]+alphaz[0])*f[13]+(alphay[7]+alphaz[6])*f[12]+f[6]*alphaz[12]+f[7]*alphay[12]+(alphay[3]+alphaz[2])*f[5]+f[2]*alphaz[5]+f[3]*alphay[5]); 
  out[23] += 0.3061862178478971*((alphavpar[19]+alphax[16])*f[31]+(alphavpar[26]+alphax[8])*f[30]+(alphavpar[11]+alphax[7])*f[29]+f[11]*alphavpar[29]+(alphavpar[10]+alphax[6])*f[28]+alphavpar[8]*f[27]+f[8]*alphavpar[27]+(alphavpar[18]+alphax[3])*f[25]+f[18]*alphavpar[25]+(alphavpar[17]+alphax[2])*f[24]+(alphavpar[4]+alphax[1])*f[23]+f[4]*alphavpar[23]+alphavpar[16]*f[22]+f[16]*alphavpar[22]+(alphax[18]+alphavpar[3])*f[21]+f[18]*alphax[21]+f[3]*alphavpar[21]+alphavpar[2]*f[20]+f[2]*alphavpar[20]+(alphavpar[9]+alphax[0])*f[15]+f[9]*alphavpar[15]+(alphax[11]+alphavpar[7])*f[14]+f[11]*alphax[14]+f[7]*alphavpar[14]+alphavpar[6]*f[13]+f[6]*alphavpar[13]+(alphax[9]+alphavpar[0])*f[12]+f[9]*alphax[12]+f[0]*alphavpar[12]+(alphax[4]+alphavpar[1])*f[5]+f[4]*alphax[5]+f[1]*alphavpar[5]); 
  out[24] += 0.3061862178478971*((alphavpar[18]+alphay[16])*f[31]+(alphavpar[11]+alphay[8])*f[30]+(alphavpar[26]+alphay[7])*f[29]+f[26]*alphavpar[29]+(alphavpar[9]+alphay[6])*f[28]+alphavpar[7]*f[27]+f[7]*alphavpar[27]+(alphavpar[19]+alphay[3])*f[25]+f[19]*alphavpar[25]+(alphavpar[4]+alphay[2])*f[24]+(alphavpar[17]+alphay[1])*f[23]+f[17]*alphavpar[23]+alphavpar[3]*f[22]+f[3]*alphavpar[22]+(alphay[18]+alphavpar[16])*f[21]+f[18]*alphay[21]+f[16]*alphavpar[21]+alphavpar[1]*f[20]+f[1]*alphavpar[20]+(alphavpar[10]+alphay[0])*f[15]+f[10]*alphavpar[15]+(alphay[11]+alphavpar[8])*f[14]+f[11]*alphay[14]+f[8]*alphavpar[14]+alphavpar[0]*f[13]+f[0]*alphavpar[13]+(alphay[9]+alphavpar[6])*f[12]+f[9]*alphay[12]+f[6]*alphavpar[12]+(alphay[4]+alphavpar[2])*f[5]+f[4]*alphay[5]+f[2]*alphavpar[5]); 
  out[25] += 0.3061862178478971*((alphavpar[17]+alphaz[16])*f[31]+(alphavpar[10]+alphaz[8])*f[30]+(alphavpar[9]+alphaz[7])*f[29]+f[9]*alphavpar[29]+(alphavpar[26]+alphaz[6])*f[28]+alphavpar[6]*f[27]+f[6]*alphavpar[27]+(alphavpar[4]+alphaz[3])*f[25]+f[4]*alphavpar[25]+(alphavpar[19]+alphaz[2])*f[24]+(alphavpar[18]+alphaz[1])*f[23]+f[18]*alphavpar[23]+alphavpar[2]*f[22]+f[2]*alphavpar[22]+(alphaz[18]+alphavpar[1])*f[21]+f[18]*alphaz[21]+f[1]*alphavpar[21]+alphavpar[16]*f[20]+f[16]*alphavpar[20]+(alphavpar[11]+alphaz[0])*f[15]+f[11]*alphavpar[15]+(alphaz[11]+alphavpar[0])*f[14]+f[11]*alphaz[14]+f[0]*alphavpar[14]+alphavpar[8]*f[13]+f[8]*alphavpar[13]+(alphaz[9]+alphavpar[7])*f[12]+f[9]*alphaz[12]+f[7]*alphavpar[12]+(alphaz[4]+alphavpar[3])*f[5]+f[4]*alphaz[5]+f[3]*alphavpar[5]); 
  out[26] += 0.3061862178478971*((alphavpar[15]+alphaz[14]+alphax[12])*f[31]+(alphavpar[23]+alphaz[21]+alphax[5])*f[30]+alphay[5]*f[29]+f[24]*alphavpar[29]+(alphavpar[25]+alphax[21]+alphaz[5])*f[28]+alphavpar[5]*f[27]+f[5]*alphavpar[27]+(alphavpar[4]+alphaz[3]+alphay[2]+alphax[1])*f[26]+f[4]*alphavpar[26]+alphay[12]*f[25]+(alphax[14]+alphaz[12])*f[24]+alphay[14]*f[23]+alphavpar[12]*f[22]+f[12]*alphavpar[22]+alphavpar[13]*f[21]+f[15]*alphay[21]+f[13]*alphavpar[21]+alphavpar[14]*f[20]+f[14]*alphavpar[20]+(alphavpar[9]+alphaz[7]+alphay[6]+alphax[0])*f[19]+f[9]*alphavpar[19]+(alphavpar[10]+alphaz[8]+alphax[6]+alphay[0])*f[18]+f[8]*alphaz[18]+f[0]*alphay[18]+f[6]*alphax[18]+f[10]*alphavpar[18]+(alphavpar[11]+alphay[8]+alphax[7]+alphaz[0])*f[17]+f[11]*alphavpar[17]+(alphaz[11]+alphax[9]+alphavpar[0])*f[16]+f[11]*alphaz[16]+f[10]*alphay[16]+f[9]*alphax[16]+f[0]*alphavpar[16]+(alphax[2]+alphay[1])*f[11]+f[1]*alphay[11]+f[2]*alphax[11]+(alphax[3]+alphaz[1])*f[10]+(alphay[3]+alphaz[2])*f[9]+f[2]*alphaz[9]+f[3]*alphay[9]+(alphax[4]+alphavpar[1])*f[8]+f[4]*alphax[8]+f[1]*alphavpar[8]+(alphay[4]+alphavpar[2])*f[7]+f[4]*alphay[7]+f[2]*alphavpar[7]+(alphaz[4]+alphavpar[3])*f[6]+f[4]*alphaz[6]+f[3]*alphavpar[6]); 
  out[27] += 0.3061862178478971*((alphaz[11]+alphax[9])*f[31]+(alphaz[18]+alphax[4])*f[30]+alphay[4]*f[29]+(alphax[18]+alphaz[4])*f[28]+(alphaz[3]+alphay[2]+alphax[1])*f[27]+alphay[9]*f[25]+(alphax[11]+alphaz[9])*f[24]+alphay[11]*f[23]+(alphaz[7]+alphay[6]+alphax[0])*f[22]+(alphaz[8]+alphax[6]+alphay[0])*f[21]+f[8]*alphaz[21]+f[0]*alphay[21]+f[6]*alphax[21]+(alphay[8]+alphax[7]+alphaz[0])*f[20]+f[15]*alphay[18]+(alphaz[14]+alphax[12])*f[16]+f[14]*alphaz[16]+f[13]*alphay[16]+f[12]*alphax[16]+(alphax[2]+alphay[1])*f[14]+f[1]*alphay[14]+f[2]*alphax[14]+(alphax[3]+alphaz[1])*f[13]+(alphay[3]+alphaz[2])*f[12]+f[2]*alphaz[12]+f[3]*alphay[12]+alphax[5]*f[8]+f[5]*alphax[8]+alphay[5]*f[7]+f[5]*alphay[7]+alphaz[5]*f[6]+f[5]*alphaz[6]); 
  out[28] += 0.3061862178478971*((alphavpar[11]+alphay[8]+alphax[7])*f[31]+(alphavpar[18]+alphay[16]+alphax[3])*f[30]+(alphavpar[19]+alphax[16]+alphay[3])*f[29]+f[19]*alphavpar[29]+(alphavpar[4]+alphay[2]+alphax[1])*f[28]+(alphax[18]+alphavpar[3])*f[27]+f[3]*alphavpar[27]+(alphavpar[25]+alphax[21])*f[26]+f[25]*(alphavpar[26]+alphax[8]+alphay[7])+(alphavpar[9]+alphay[6]+alphax[0])*f[24]+(alphavpar[10]+alphax[6]+alphay[0])*f[23]+f[10]*alphavpar[23]+(alphax[11]+alphavpar[7])*f[22]+f[7]*alphavpar[22]+(alphay[11]+alphavpar[8])*f[21]+f[11]*alphay[21]+f[8]*alphavpar[21]+(alphax[9]+alphavpar[0])*f[20]+f[0]*alphavpar[20]+alphax[14]*f[19]+alphay[14]*f[18]+f[14]*alphay[18]+(alphavpar[15]+alphax[12])*f[17]+f[15]*alphavpar[17]+alphavpar[14]*f[16]+f[14]*alphavpar[16]+(alphax[2]+alphay[1])*f[15]+(alphax[4]+alphavpar[1])*f[13]+f[1]*alphavpar[13]+(alphay[4]+alphavpar[2])*f[12]+f[4]*alphay[12]+f[2]*alphavpar[12]+alphax[5]*f[10]+alphay[5]*f[9]+f[5]*alphay[9]+alphavpar[5]*f[6]+f[5]*alphavpar[6]); 
  out[29] += 0.3061862178478971*((alphavpar[10]+alphaz[8]+alphax[6])*f[31]+(alphavpar[17]+alphaz[16]+alphax[2])*f[30]+(alphavpar[4]+alphaz[3]+alphax[1])*f[29]+f[4]*alphavpar[29]+(alphavpar[19]+alphax[16]+alphaz[2])*f[28]+alphavpar[2]*f[27]+f[2]*alphavpar[27]+f[24]*alphavpar[26]+(alphavpar[9]+alphaz[7]+alphax[0])*f[25]+f[9]*alphavpar[25]+(alphax[8]+alphaz[6])*f[24]+(alphavpar[11]+alphax[7]+alphaz[0])*f[23]+f[11]*alphavpar[23]+alphavpar[6]*f[22]+f[6]*alphavpar[22]+(alphaz[11]+alphax[9]+alphavpar[0])*f[21]+f[11]*alphaz[21]+f[9]*alphax[21]+f[0]*alphavpar[21]+alphavpar[8]*f[20]+f[8]*alphavpar[20]+(alphavpar[15]+alphaz[14]+alphax[12])*f[18]+f[14]*alphaz[18]+f[12]*alphax[18]+f[15]*alphavpar[18]+alphavpar[13]*f[16]+f[13]*alphavpar[16]+(alphax[3]+alphaz[1])*f[15]+(alphax[4]+alphavpar[1])*f[14]+f[4]*alphax[14]+f[1]*alphavpar[14]+(alphaz[4]+alphavpar[3])*f[12]+f[4]*alphaz[12]+f[3]*alphavpar[12]+alphax[5]*f[11]+f[5]*alphax[11]+alphaz[5]*f[9]+f[5]*alphaz[9]+alphavpar[5]*f[7]+f[5]*alphavpar[7]); 
  out[30] += 0.3061862178478971*((alphavpar[9]+alphaz[7]+alphay[6])*f[31]+(alphavpar[4]+alphaz[3]+alphay[2])*f[30]+(alphavpar[17]+alphaz[16]+alphay[1])*f[29]+f[17]*alphavpar[29]+(alphavpar[18]+alphay[16]+alphaz[1])*f[28]+(alphaz[18]+alphavpar[1])*f[27]+f[1]*alphavpar[27]+(alphavpar[23]+alphaz[21])*f[26]+f[23]*alphavpar[26]+(alphavpar[10]+alphaz[8]+alphay[0])*f[25]+f[10]*alphavpar[25]+(alphavpar[11]+alphay[8]+alphaz[0])*f[24]+(alphay[7]+alphaz[6])*f[23]+(alphaz[11]+alphavpar[0])*f[22]+f[0]*alphavpar[22]+(alphay[9]+alphavpar[6])*f[21]+f[9]*alphay[21]+f[6]*alphavpar[21]+(alphaz[9]+alphavpar[7])*f[20]+f[7]*alphavpar[20]+(alphavpar[15]+alphaz[14])*f[19]+f[15]*alphavpar[19]+alphay[12]*f[18]+f[12]*alphay[18]+alphaz[12]*f[17]+alphavpar[12]*f[16]+f[12]*alphavpar[16]+(alphay[3]+alphaz[2])*f[15]+(alphay[4]+alphavpar[2])*f[14]+f[4]*alphay[14]+f[2]*alphavpar[14]+(alphaz[4]+alphavpar[3])*f[13]+f[3]*alphavpar[13]+alphay[5]*f[11]+f[5]*alphay[11]+alphaz[5]*f[10]+alphavpar[5]*f[8]+f[5]*alphavpar[8]); 
  out[31] += 0.3061862178478971*((alphavpar[4]+alphaz[3]+alphay[2]+alphax[1])*f[31]+(alphavpar[9]+alphaz[7]+alphay[6]+alphax[0])*f[30]+(alphavpar[10]+alphaz[8]+alphax[6]+alphay[0])*f[29]+f[10]*alphavpar[29]+(alphavpar[11]+alphay[8]+alphax[7]+alphaz[0])*f[28]+(alphaz[11]+alphax[9]+alphavpar[0])*f[27]+f[0]*alphavpar[27]+(alphavpar[15]+alphaz[14]+alphax[12])*f[26]+f[15]*alphavpar[26]+(alphavpar[17]+alphaz[16]+alphax[2]+alphay[1])*f[25]+f[17]*alphavpar[25]+(alphavpar[18]+alphay[16]+alphax[3]+alphaz[1])*f[24]+(alphavpar[19]+alphax[16]+alphay[3]+alphaz[2])*f[23]+f[19]*alphavpar[23]+(alphaz[18]+alphax[4]+alphavpar[1])*f[22]+f[1]*alphavpar[22]+(alphay[4]+alphavpar[2])*f[21]+f[19]*alphaz[21]+f[4]*alphay[21]+f[17]*alphax[21]+f[2]*alphavpar[21]+(alphax[18]+alphaz[4]+alphavpar[3])*f[20]+f[3]*alphavpar[20]+alphax[5]*f[19]+alphay[5]*f[18]+f[5]*alphay[18]+alphaz[5]*f[17]+alphavpar[5]*f[16]+f[5]*alphavpar[16]+(alphax[8]+alphay[7]+alphaz[6])*f[15]+(alphay[9]+alphavpar[6])*f[14]+f[9]*alphay[14]+f[10]*alphax[14]+f[6]*alphavpar[14]+(alphax[11]+alphaz[9]+alphavpar[7])*f[13]+f[7]*alphavpar[13]+(alphay[11]+alphavpar[8])*f[12]+f[10]*alphaz[12]+f[11]*alphay[12]+f[8]*alphavpar[12]); 
  return cflFreq; 
} 
double EmGyrokineticGenGeoStep2Vol3x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
{ 
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // dApardt: time derivative of Apar.
  // f: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[3]; 
  out[4] += -(0.6123724356957944*(dApardt[7]*f[16]+dApardt[6]*f[8]+dApardt[5]*f[7]+dApardt[4]*f[6]+dApardt[3]*f[3]+dApardt[2]*f[2]+dApardt[1]*f[1]+dApardt[0]*f[0])*q_*rdvpar2)/m_; 
  out[9] += -(0.6123724356957944*(dApardt[6]*f[16]+dApardt[7]*f[8]+dApardt[3]*f[7]+dApardt[2]*f[6]+f[3]*dApardt[5]+f[2]*dApardt[4]+dApardt[0]*f[1]+f[0]*dApardt[1])*q_*rdvpar2)/m_; 
  out[10] += -(0.6123724356957944*(dApardt[5]*f[16]+dApardt[3]*f[8]+dApardt[7]*f[7]+dApardt[1]*f[6]+f[3]*dApardt[6]+f[1]*dApardt[4]+dApardt[0]*f[2]+f[0]*dApardt[2])*q_*rdvpar2)/m_; 
  out[11] += -(0.6123724356957944*(dApardt[4]*f[16]+dApardt[2]*f[8]+dApardt[1]*f[7]+f[6]*dApardt[7]+f[2]*dApardt[6]+f[1]*dApardt[5]+dApardt[0]*f[3]+f[0]*dApardt[3])*q_*rdvpar2)/m_; 
  out[15] += -(0.6123724356957944*(dApardt[7]*f[27]+dApardt[6]*f[22]+dApardt[5]*f[21]+dApardt[4]*f[20]+dApardt[3]*f[14]+dApardt[2]*f[13]+dApardt[1]*f[12]+dApardt[0]*f[5])*q_*rdvpar2)/m_; 
  out[17] += -(0.6123724356957944*(dApardt[3]*f[16]+dApardt[5]*f[8]+dApardt[6]*f[7]+f[3]*dApardt[7]+dApardt[0]*f[6]+f[0]*dApardt[4]+dApardt[1]*f[2]+f[1]*dApardt[2])*q_*rdvpar2)/m_; 
  out[18] += -(0.6123724356957944*(dApardt[2]*f[16]+dApardt[4]*f[8]+dApardt[0]*f[7]+f[2]*dApardt[7]+dApardt[6]*f[6]+f[0]*dApardt[5]+dApardt[1]*f[3]+f[1]*dApardt[3])*q_*rdvpar2)/m_; 
  out[19] += -(0.6123724356957944*(dApardt[1]*f[16]+dApardt[0]*f[8]+dApardt[4]*f[7]+f[1]*dApardt[7]+dApardt[5]*f[6]+f[0]*dApardt[6]+dApardt[2]*f[3]+f[2]*dApardt[3])*q_*rdvpar2)/m_; 
  out[23] += -(0.6123724356957944*(dApardt[6]*f[27]+dApardt[7]*f[22]+dApardt[3]*f[21]+dApardt[2]*f[20]+dApardt[5]*f[14]+dApardt[4]*f[13]+dApardt[0]*f[12]+dApardt[1]*f[5])*q_*rdvpar2)/m_; 
  out[24] += -(0.6123724356957944*(dApardt[5]*f[27]+dApardt[3]*f[22]+dApardt[7]*f[21]+dApardt[1]*f[20]+dApardt[6]*f[14]+dApardt[0]*f[13]+dApardt[4]*f[12]+dApardt[2]*f[5])*q_*rdvpar2)/m_; 
  out[25] += -(0.6123724356957944*(dApardt[4]*f[27]+dApardt[2]*f[22]+dApardt[1]*f[21]+dApardt[7]*f[20]+dApardt[0]*f[14]+dApardt[6]*f[13]+dApardt[5]*f[12]+dApardt[3]*f[5])*q_*rdvpar2)/m_; 
  out[26] += -(0.6123724356957944*(dApardt[0]*f[16]+dApardt[1]*f[8]+dApardt[2]*f[7]+f[0]*dApardt[7]+dApardt[3]*f[6]+f[1]*dApardt[6]+f[2]*dApardt[5]+f[3]*dApardt[4])*q_*rdvpar2)/m_; 
  out[28] += -(0.6123724356957944*(dApardt[3]*f[27]+dApardt[5]*f[22]+dApardt[6]*f[21]+dApardt[0]*f[20]+dApardt[7]*f[14]+dApardt[1]*f[13]+dApardt[2]*f[12]+dApardt[4]*f[5])*q_*rdvpar2)/m_; 
  out[29] += -(0.6123724356957944*(dApardt[2]*f[27]+dApardt[4]*f[22]+dApardt[0]*f[21]+dApardt[6]*f[20]+dApardt[1]*f[14]+dApardt[7]*f[13]+dApardt[3]*f[12]+dApardt[5]*f[5])*q_*rdvpar2)/m_; 
  out[30] += -(0.6123724356957944*(dApardt[1]*f[27]+dApardt[0]*f[22]+dApardt[4]*f[21]+dApardt[5]*f[20]+dApardt[2]*f[14]+dApardt[3]*f[13]+dApardt[7]*f[12]+f[5]*dApardt[6])*q_*rdvpar2)/m_; 
  out[31] += -(0.6123724356957944*(dApardt[0]*f[27]+dApardt[1]*f[22]+dApardt[2]*f[21]+dApardt[3]*f[20]+dApardt[4]*f[14]+dApardt[5]*f[13]+dApardt[6]*f[12]+f[5]*dApardt[7])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -(0.1767766952966368*dApardt[0]*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = -(0.1767766952966368*dApardt[0]*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = -(0.03125*((-0.3535533905932737*dApardt[7])+0.3535533905932737*(dApardt[6]+dApardt[5]+dApardt[4])-0.3535533905932737*(dApardt[3]+dApardt[2]+dApardt[1])+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(0.3535533905932737*(dApardt[7]+dApardt[6])-0.3535533905932737*(dApardt[5]+dApardt[4]+dApardt[3]+dApardt[2])+0.3535533905932737*(dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(0.3535533905932737*dApardt[7]-0.3535533905932737*dApardt[6]+0.3535533905932737*dApardt[5]-0.3535533905932737*(dApardt[4]+dApardt[3])+0.3535533905932737*dApardt[2]-0.3535533905932737*dApardt[1]+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*((-0.3535533905932737*(dApardt[7]+dApardt[6]+dApardt[5]))+0.3535533905932737*dApardt[4]-0.3535533905932737*dApardt[3]+0.3535533905932737*(dApardt[2]+dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(0.3535533905932737*dApardt[7]-0.3535533905932737*(dApardt[6]+dApardt[5])+0.3535533905932737*(dApardt[4]+dApardt[3])-0.3535533905932737*(dApardt[2]+dApardt[1])+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*((-0.3535533905932737*(dApardt[7]+dApardt[6]))+0.3535533905932737*dApardt[5]-0.3535533905932737*dApardt[4]+0.3535533905932737*dApardt[3]-0.3535533905932737*dApardt[2]+0.3535533905932737*(dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*((-0.3535533905932737*dApardt[7])+0.3535533905932737*dApardt[6]-0.3535533905932737*(dApardt[5]+dApardt[4])+0.3535533905932737*(dApardt[3]+dApardt[2])-0.3535533905932737*dApardt[1]+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0110485434560398*(dApardt[7]+dApardt[6]+dApardt[5]+dApardt[4]+dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*((-0.3535533905932737*dApardt[7])+0.3535533905932737*(dApardt[6]+dApardt[5]+dApardt[4])-0.3535533905932737*(dApardt[3]+dApardt[2]+dApardt[1])+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(0.3535533905932737*(dApardt[7]+dApardt[6])-0.3535533905932737*(dApardt[5]+dApardt[4]+dApardt[3]+dApardt[2])+0.3535533905932737*(dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(0.3535533905932737*dApardt[7]-0.3535533905932737*dApardt[6]+0.3535533905932737*dApardt[5]-0.3535533905932737*(dApardt[4]+dApardt[3])+0.3535533905932737*dApardt[2]-0.3535533905932737*dApardt[1]+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*((-0.3535533905932737*(dApardt[7]+dApardt[6]+dApardt[5]))+0.3535533905932737*dApardt[4]-0.3535533905932737*dApardt[3]+0.3535533905932737*(dApardt[2]+dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*(0.3535533905932737*dApardt[7]-0.3535533905932737*(dApardt[6]+dApardt[5])+0.3535533905932737*(dApardt[4]+dApardt[3])-0.3535533905932737*(dApardt[2]+dApardt[1])+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*((-0.3535533905932737*(dApardt[7]+dApardt[6]))+0.3535533905932737*dApardt[5]-0.3535533905932737*dApardt[4]+0.3535533905932737*dApardt[3]-0.3535533905932737*dApardt[2]+0.3535533905932737*(dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.03125*((-0.3535533905932737*dApardt[7])+0.3535533905932737*dApardt[6]-0.3535533905932737*(dApardt[5]+dApardt[4])+0.3535533905932737*(dApardt[3]+dApardt[2])-0.3535533905932737*dApardt[1]+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = -(0.0110485434560398*(dApardt[7]+dApardt[6]+dApardt[5]+dApardt[4]+dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = -(0.03125*((-0.3535533905932737*dApardt[7])+0.3535533905932737*(dApardt[6]+dApardt[5]+dApardt[4])-0.3535533905932737*(dApardt[3]+dApardt[2]+dApardt[1])+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(0.3535533905932737*(dApardt[7]+dApardt[6])-0.3535533905932737*(dApardt[5]+dApardt[4]+dApardt[3]+dApardt[2])+0.3535533905932737*(dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(0.3535533905932737*dApardt[7]-0.3535533905932737*dApardt[6]+0.3535533905932737*dApardt[5]-0.3535533905932737*(dApardt[4]+dApardt[3])+0.3535533905932737*dApardt[2]-0.3535533905932737*dApardt[1]+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*((-0.3535533905932737*(dApardt[7]+dApardt[6]+dApardt[5]))+0.3535533905932737*dApardt[4]-0.3535533905932737*dApardt[3]+0.3535533905932737*(dApardt[2]+dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(0.3535533905932737*dApardt[7]-0.3535533905932737*(dApardt[6]+dApardt[5])+0.3535533905932737*(dApardt[4]+dApardt[3])-0.3535533905932737*(dApardt[2]+dApardt[1])+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*((-0.3535533905932737*(dApardt[7]+dApardt[6]))+0.3535533905932737*dApardt[5]-0.3535533905932737*dApardt[4]+0.3535533905932737*dApardt[3]-0.3535533905932737*dApardt[2]+0.3535533905932737*(dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*((-0.3535533905932737*dApardt[7])+0.3535533905932737*dApardt[6]-0.3535533905932737*(dApardt[5]+dApardt[4])+0.3535533905932737*(dApardt[3]+dApardt[2])-0.3535533905932737*dApardt[1]+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0110485434560398*(dApardt[7]+dApardt[6]+dApardt[5]+dApardt[4]+dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*((-0.3535533905932737*dApardt[7])+0.3535533905932737*(dApardt[6]+dApardt[5]+dApardt[4])-0.3535533905932737*(dApardt[3]+dApardt[2]+dApardt[1])+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(0.3535533905932737*(dApardt[7]+dApardt[6])-0.3535533905932737*(dApardt[5]+dApardt[4]+dApardt[3]+dApardt[2])+0.3535533905932737*(dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(0.3535533905932737*dApardt[7]-0.3535533905932737*dApardt[6]+0.3535533905932737*dApardt[5]-0.3535533905932737*(dApardt[4]+dApardt[3])+0.3535533905932737*dApardt[2]-0.3535533905932737*dApardt[1]+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*((-0.3535533905932737*(dApardt[7]+dApardt[6]+dApardt[5]))+0.3535533905932737*dApardt[4]-0.3535533905932737*dApardt[3]+0.3535533905932737*(dApardt[2]+dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*(0.3535533905932737*dApardt[7]-0.3535533905932737*(dApardt[6]+dApardt[5])+0.3535533905932737*(dApardt[4]+dApardt[3])-0.3535533905932737*(dApardt[2]+dApardt[1])+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*((-0.3535533905932737*(dApardt[7]+dApardt[6]))+0.3535533905932737*dApardt[5]-0.3535533905932737*dApardt[4]+0.3535533905932737*dApardt[3]-0.3535533905932737*dApardt[2]+0.3535533905932737*(dApardt[1]+dApardt[0]))*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.03125*((-0.3535533905932737*dApardt[7])+0.3535533905932737*dApardt[6]-0.3535533905932737*(dApardt[5]+dApardt[4])+0.3535533905932737*(dApardt[3]+dApardt[2])-0.3535533905932737*dApardt[1]+0.3535533905932737*dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = -(0.0110485434560398*(dApardt[7]+dApardt[6]+dApardt[5]+dApardt[4]+dApardt[3]+dApardt[2]+dApardt[1]+dApardt[0])*q_*rdvpar2)/m_; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

return cflFreq; 
} 
