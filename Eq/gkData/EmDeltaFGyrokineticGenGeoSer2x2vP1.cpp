#include <DeltaFGyrokineticModDecl.h> 
double EmDeltaFGyrokineticGenGeoVol2x2vSerP1_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f0, const double *f1, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
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

  double hamil0[16]; 
  hamil0[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu)+m_))/rdvpar2Sq; 
  hamil0[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil0[4] = (1.154700538379252*bmag[0])/rdmu2; 

  double hamil1[16]; 
  hamil1[0] = 2.0*phi[0]*q_; 
  hamil1[1] = 2.0*phi[1]*q_; 
  hamil1[2] = 2.0*phi[2]*q_; 
  hamil1[5] = 2.0*phi[3]*q_; 

  double BstarXdBmag[16]; 

  double BstarYdBmag[16]; 

  double BstarXdBmagEM[16]; 
  BstarXdBmagEM[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2; 
  BstarXdBmagEM[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2; 

  double BstarYdBmagEM[16]; 
  BstarYdBmagEM[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2; 
  BstarYdBmagEM[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[16]; 
  double alpha1x[16]; 
  alpha1x[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmagEM[0]*hamil0[3]*rdvpar2)/m_-(0.25*b_z[0]*jacobTotInv[0]*hamil1[2]*rdy2)/q_); 
  alpha1x[1] = 1.732050807568877*rdx2*((0.25*BstarXdBmagEM[1]*hamil0[3]*rdvpar2)/m_-(0.25*b_z[0]*jacobTotInv[0]*hamil1[5]*rdy2)/q_); 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0y[16]; 
  double alpha1y[16]; 
  alpha1y[0] = 0.4330127018922193*((b_z[0]*jacobTotInv[0]*hamil1[1]*rdx2)/q_+(BstarYdBmagEM[0]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha1y[2] = 0.4330127018922193*((b_z[0]*jacobTotInv[0]*hamil1[5]*rdx2)/q_+(BstarYdBmagEM[2]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[16]; 
  double alpha1vpar[16]; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.4330127018922193*(alpha1x[1]*f0[1]+alpha1x[0]*f0[0]); 
  out[2] += 0.4330127018922193*(alpha1y[2]*f0[2]+alpha1y[0]*f0[0]); 
  out[5] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f0[5]+alpha1x[0]*f0[2]+alpha1y[0]*f0[1]); 
  out[6] += 0.4330127018922193*(alpha1x[1]*f0[6]+alpha1x[0]*f0[3]); 
  out[7] += 0.4330127018922193*(alpha1y[2]*f0[7]+alpha1y[0]*f0[3]); 
  out[8] += 0.4330127018922193*(alpha1x[1]*f0[8]+alpha1x[0]*f0[4]); 
  out[9] += 0.4330127018922193*(alpha1y[2]*f0[9]+alpha1y[0]*f0[4]); 
  out[11] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f0[11]+alpha1x[0]*f0[7]+alpha1y[0]*f0[6]); 
  out[12] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f0[12]+alpha1x[0]*f0[9]+alpha1y[0]*f0[8]); 
  out[13] += 0.4330127018922193*(alpha1x[1]*f0[13]+alpha1x[0]*f0[10]); 
  out[14] += 0.4330127018922193*(alpha1y[2]*f0[14]+alpha1y[0]*f0[10]); 
  out[15] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f0[15]+alpha1x[0]*f0[14]+alpha1y[0]*f0[13]); 
  out[1] += 0.4330127018922193*(alpha1x[1]*f1[1]+alpha1x[0]*f1[0]); 
  out[2] += 0.4330127018922193*(alpha1y[2]*f1[2]+alpha1y[0]*f1[0]); 
  out[5] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f1[5]+alpha1x[0]*f1[2]+alpha1y[0]*f1[1]); 
  out[6] += 0.4330127018922193*(alpha1x[1]*f1[6]+alpha1x[0]*f1[3]); 
  out[7] += 0.4330127018922193*(alpha1y[2]*f1[7]+alpha1y[0]*f1[3]); 
  out[8] += 0.4330127018922193*(alpha1x[1]*f1[8]+alpha1x[0]*f1[4]); 
  out[9] += 0.4330127018922193*(alpha1y[2]*f1[9]+alpha1y[0]*f1[4]); 
  out[11] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f1[11]+alpha1x[0]*f1[7]+alpha1y[0]*f1[6]); 
  out[12] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f1[12]+alpha1x[0]*f1[9]+alpha1y[0]*f1[8]); 
  out[13] += 0.4330127018922193*(alpha1x[1]*f1[13]+alpha1x[0]*f1[10]); 
  out[14] += 0.4330127018922193*(alpha1y[2]*f1[14]+alpha1y[0]*f1[10]); 
  out[15] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f1[15]+alpha1x[0]*f1[14]+alpha1y[0]*f1[13]); 
  return cflFreq; 
} 
double EmDeltaFGyrokineticGenGeoVol2x2vSerP1_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f0, const double *f1, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
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

  double hamil0[16]; 
  hamil0[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu)+m_))/rdvpar2Sq; 
  hamil0[1] = 2.0*bmag[1]*wmu; 
  hamil0[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil0[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil0[8] = (1.154700538379252*bmag[1])/rdmu2; 

  double hamil1[16]; 
  hamil1[0] = 2.0*phi[0]*q_; 
  hamil1[1] = 2.0*phi[1]*q_; 
  hamil1[2] = 2.0*phi[2]*q_; 
  hamil1[5] = 2.0*phi[3]*q_; 

  double BstarXdBmag[16]; 

  double BstarYdBmag[16]; 
  BstarYdBmag[0] = -(1.732050807568877*jacobTotInv[0]*b_z[1]*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[1] = -(1.732050807568877*b_z[1]*jacobTotInv[1]*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 

  double BstarXdBmagEM[16]; 
  BstarXdBmagEM[0] = 0.8660254037844386*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])*rdy2; 
  BstarXdBmagEM[1] = 0.1732050807568877*((9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])*rdy2; 

  double BstarYdBmagEM[16]; 
  BstarYdBmagEM[0] = -0.8660254037844386*(2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*rdx2; 
  BstarYdBmagEM[1] = -0.8660254037844386*((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])*rdx2; 
  BstarYdBmagEM[2] = -0.8660254037844386*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])*rdx2; 
  BstarYdBmagEM[5] = -0.8660254037844386*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])*rdx2; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[16]; 
  double alpha1x[16]; 
  alpha1x[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmagEM[0]*hamil0[3]*rdvpar2)/m_-(0.25*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil1[5]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil1[2])*rdy2)/q_); 
  alpha1x[1] = 1.732050807568877*rdx2*((((-0.25*(b_z[0]*jacobTotInv[0]*hamil1[5]+(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil1[2]))-0.45*b_z[1]*jacobTotInv[1]*hamil1[5])*rdy2)/q_+(0.25*BstarXdBmagEM[1]*hamil0[3]*rdvpar2)/m_); 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0y[16]; 
  alpha0y[0] = 0.4330127018922193*((hamil0[1]*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*rdx2)/q_+(BstarYdBmag[0]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha0y[1] = 0.4330127018922193*((hamil0[1]*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*rdx2)/q_+(BstarYdBmag[1]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha0y[3] = (0.4330127018922193*BstarYdBmag[3]*hamil0[3]*rdvpar2*rdy2)/m_; 
  alpha0y[4] = (0.4330127018922193*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil0[8]*rdx2*rdy2)/q_; 
  alpha0y[6] = (0.4330127018922193*hamil0[3]*BstarYdBmag[6]*rdvpar2*rdy2)/m_; 
  alpha0y[8] = (0.4330127018922193*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil0[8]*rdx2*rdy2)/q_; 
  double alpha1y[16]; 
  alpha1y[0] = 0.4330127018922193*(((hamil1[1]+hamil0[1])*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*rdx2)/q_+(BstarYdBmagEM[0]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha1y[1] = 0.4330127018922193*(((hamil1[1]+hamil0[1])*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*rdx2)/q_+(BstarYdBmagEM[1]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha1y[2] = 0.4330127018922193*(((b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil1[5]*rdx2)/q_+(BstarYdBmagEM[2]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha1y[4] = (0.4330127018922193*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil0[8]*rdx2*rdy2)/q_; 
  alpha1y[5] = 0.4330127018922193*(((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil1[5]*rdx2)/q_+(hamil0[3]*BstarYdBmagEM[5]*rdvpar2)/m_)*rdy2; 
  alpha1y[8] = (0.4330127018922193*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil0[8]*rdx2*rdy2)/q_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6])+0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])-0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6]))-0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])-0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alpha1y[8]+alpha0y[8])-0.25*alpha0y[6]+0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4])+0.25*alpha0y[3]-0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]))+0.25*alpha0y[6]-0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4])+0.25*alpha0y[3]-0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]))+0.25*alpha0y[6]+0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4])-0.25*alpha0y[3]-0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alpha1y[8]+alpha0y[8])-0.25*alpha0y[6]-0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4])-0.25*alpha0y[3]-0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6]))+0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])-0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6])-0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])-0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6])-0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])+0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6]))+0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])+0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alpha1y[8]+alpha0y[8])-0.25*alpha0y[6]-0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4])+0.25*alpha0y[3]+0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]))+0.25*alpha0y[6]+0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4])+0.25*alpha0y[3]+0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]))+0.25*alpha0y[6]-0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4])-0.25*alpha0y[3]+0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alpha1y[8]+alpha0y[8])-0.25*alpha0y[6]+0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4])-0.25*alpha0y[3]+0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6]))-0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])+0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6])+0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])+0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[16]; 
  double alpha1vpar[16]; 
  alpha1vpar[0] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[1]*hamil1[5]+BstarYdBmag[0]*hamil1[2])*rdy2+BstarXdBmagEM[0]*hamil0[1]*rdx2))/m_; 
  alpha1vpar[1] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[0]*hamil1[5]+BstarYdBmag[1]*hamil1[2])*rdy2+BstarXdBmagEM[1]*hamil0[1]*rdx2))/m_; 
  alpha1vpar[3] = -(0.4330127018922193*(hamil1[5]*BstarYdBmag[6]+hamil1[2]*BstarYdBmag[3])*rdvpar2*rdy2)/m_; 
  alpha1vpar[4] = -(0.4330127018922193*BstarXdBmagEM[0]*hamil0[8]*rdvpar2*rdx2)/m_; 
  alpha1vpar[6] = -(0.4330127018922193*(hamil1[2]*BstarYdBmag[6]+BstarYdBmag[3]*hamil1[5])*rdvpar2*rdy2)/m_; 
  alpha1vpar[8] = -(0.4330127018922193*BstarXdBmagEM[1]*hamil0[8]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alpha1vpar[8]+0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alpha1vpar[8])-0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1vpar[8]+0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alpha1vpar[8])-0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alpha1vpar[8])+0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1vpar[8]-0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alpha1vpar[8])+0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1vpar[8]-0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*alpha1vpar[8]-0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alpha1vpar[8])+0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alpha1vpar[8]-0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alpha1vpar[8])+0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alpha1vpar[8])-0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alpha1vpar[8]+0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alpha1vpar[8])-0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alpha1vpar[8]+0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.4330127018922193*(alpha1x[1]*f0[1]+alpha1x[0]*f0[0]); 
  out[2] += 0.4330127018922193*(alpha0y[8]*f1[8]+alpha1y[8]*f0[8]+alpha0y[6]*f1[6]+alpha1y[5]*f0[5]+alpha0y[4]*f1[4]+alpha1y[4]*f0[4]+alpha0y[3]*f1[3]+alpha1y[2]*f0[2]+alpha0y[1]*f1[1]+alpha1y[1]*f0[1]+alpha0y[0]*f1[0]+alpha1y[0]*f0[0]); 
  out[3] += 0.4330127018922193*(alpha1vpar[8]*f0[8]+alpha1vpar[6]*f0[6]+alpha1vpar[4]*f0[4]+alpha1vpar[3]*f0[3]+alpha1vpar[1]*f0[1]+alpha1vpar[0]*f0[0]); 
  out[5] += 0.4330127018922193*(alpha0y[4]*f1[8]+alpha1y[4]*f0[8]+f0[4]*alpha1y[8]+f1[4]*alpha0y[8]+alpha0y[3]*f1[6]+f1[3]*alpha0y[6]+(alpha1y[2]+alpha1x[1])*f0[5]+f0[2]*(alpha1y[5]+alpha1x[0])+alpha0y[0]*f1[1]+alpha1y[0]*f0[1]+f0[0]*alpha1y[1]+f1[0]*alpha0y[1]); 
  out[6] += 0.4330127018922193*(alpha1vpar[4]*f0[8]+f0[4]*alpha1vpar[8]+(alpha1vpar[3]+alpha1x[1])*f0[6]+f0[3]*(alpha1vpar[6]+alpha1x[0])+alpha1vpar[0]*f0[1]+f0[0]*alpha1vpar[1]); 
  out[7] += 0.4330127018922193*(alpha0y[8]*f1[13]+alpha1y[8]*f0[13]+alpha1vpar[8]*f0[12]+(alpha1vpar[6]+alpha1y[5])*f0[11]+alpha0y[4]*f1[10]+alpha1y[4]*f0[10]+alpha1vpar[4]*f0[9]+(alpha1vpar[3]+alpha1y[2])*f0[7]+alpha0y[1]*f1[6]+alpha1y[1]*f0[6]+f1[1]*alpha0y[6]+alpha1vpar[1]*f0[5]+alpha0y[0]*f1[3]+alpha1y[0]*f0[3]+f1[0]*alpha0y[3]+alpha1vpar[0]*f0[2]); 
  out[8] += 0.4330127018922193*(alpha1x[1]*f0[8]+alpha1x[0]*f0[4]); 
  out[9] += 0.4330127018922193*(alpha0y[6]*f1[13]+alpha1y[5]*f0[12]+alpha0y[3]*f1[10]+alpha1y[2]*f0[9]+alpha0y[1]*f1[8]+alpha1y[1]*f0[8]+f0[1]*alpha1y[8]+f1[1]*alpha0y[8]+alpha0y[0]*f1[4]+alpha1y[0]*f0[4]+f0[0]*alpha1y[4]+f1[0]*alpha0y[4]); 
  out[10] += 0.4330127018922193*(alpha1vpar[6]*f0[13]+alpha1vpar[3]*f0[10]+alpha1vpar[1]*f0[8]+f0[1]*alpha1vpar[8]+alpha1vpar[0]*f0[4]+f0[0]*alpha1vpar[4]); 
  out[11] += 0.4330127018922193*(alpha0y[4]*f1[13]+alpha1y[4]*f0[13]+alpha1vpar[4]*f0[12]+(alpha1vpar[3]+alpha1y[2]+alpha1x[1])*f0[11]+alpha0y[8]*f1[10]+alpha1y[8]*f0[10]+alpha1vpar[8]*f0[9]+(alpha1vpar[6]+alpha1y[5]+alpha1x[0])*f0[7]+alpha0y[0]*f1[6]+alpha1y[0]*f0[6]+f1[0]*alpha0y[6]+alpha1vpar[0]*f0[5]+alpha0y[1]*f1[3]+alpha1y[1]*f0[3]+f1[1]*alpha0y[3]+alpha1vpar[1]*f0[2]); 
  out[12] += 0.4330127018922193*(alpha0y[3]*f1[13]+(alpha1y[2]+alpha1x[1])*f0[12]+alpha0y[6]*f1[10]+(alpha1y[5]+alpha1x[0])*f0[9]+alpha0y[0]*f1[8]+alpha1y[0]*f0[8]+f0[0]*alpha1y[8]+f1[0]*alpha0y[8]+alpha0y[1]*f1[4]+alpha1y[1]*f0[4]+f0[1]*alpha1y[4]+f1[1]*alpha0y[4]); 
  out[13] += 0.4330127018922193*((alpha1vpar[3]+alpha1x[1])*f0[13]+(alpha1vpar[6]+alpha1x[0])*f0[10]+alpha1vpar[0]*f0[8]+f0[0]*alpha1vpar[8]+alpha1vpar[1]*f0[4]+f0[1]*alpha1vpar[4]); 
  out[14] += 0.4330127018922193*((alpha1vpar[6]+alpha1y[5])*f0[15]+(alpha1vpar[3]+alpha1y[2])*f0[14]+alpha0y[1]*f1[13]+alpha1y[1]*f0[13]+alpha1vpar[1]*f0[12]+alpha0y[0]*f1[10]+alpha1y[0]*f0[10]+alpha1vpar[0]*f0[9]+alpha0y[6]*f1[8]+f0[6]*alpha1y[8]+f0[5]*alpha1vpar[8]+f1[6]*alpha0y[8]+alpha0y[3]*f1[4]+f0[3]*alpha1y[4]+f0[2]*alpha1vpar[4]+f1[3]*alpha0y[4]); 
  out[15] += 0.4330127018922193*((alpha1vpar[3]+alpha1y[2]+alpha1x[1])*f0[15]+(alpha1vpar[6]+alpha1y[5]+alpha1x[0])*f0[14]+alpha0y[0]*f1[13]+alpha1y[0]*f0[13]+alpha1vpar[0]*f0[12]+alpha0y[1]*f1[10]+alpha1y[1]*f0[10]+alpha1vpar[1]*f0[9]+alpha0y[3]*f1[8]+f0[3]*alpha1y[8]+f0[2]*alpha1vpar[8]+f1[3]*alpha0y[8]+alpha0y[4]*f1[6]+alpha1y[4]*f0[6]+f1[4]*alpha0y[6]+alpha1vpar[4]*f0[5]); 
  out[1] += 0.4330127018922193*(alpha1x[1]*f1[1]+alpha1x[0]*f1[0]); 
  out[2] += 0.4330127018922193*(alpha1y[8]*f1[8]+alpha1y[5]*f1[5]+alpha1y[4]*f1[4]+alpha1y[2]*f1[2]+alpha1y[1]*f1[1]+alpha1y[0]*f1[0]); 
  out[3] += 0.4330127018922193*(alpha1vpar[8]*f1[8]+alpha1vpar[6]*f1[6]+alpha1vpar[4]*f1[4]+alpha1vpar[3]*f1[3]+alpha1vpar[1]*f1[1]+alpha1vpar[0]*f1[0]); 
  out[5] += 0.4330127018922193*(alpha1y[4]*f1[8]+f1[4]*alpha1y[8]+(alpha1y[2]+alpha1x[1])*f1[5]+f1[2]*(alpha1y[5]+alpha1x[0])+alpha1y[0]*f1[1]+f1[0]*alpha1y[1]); 
  out[6] += 0.4330127018922193*(alpha1vpar[4]*f1[8]+f1[4]*alpha1vpar[8]+(alpha1vpar[3]+alpha1x[1])*f1[6]+f1[3]*(alpha1vpar[6]+alpha1x[0])+alpha1vpar[0]*f1[1]+f1[0]*alpha1vpar[1]); 
  out[7] += 0.4330127018922193*(alpha1y[8]*f1[13]+alpha1vpar[8]*f1[12]+(alpha1vpar[6]+alpha1y[5])*f1[11]+alpha1y[4]*f1[10]+alpha1vpar[4]*f1[9]+(alpha1vpar[3]+alpha1y[2])*f1[7]+alpha1y[1]*f1[6]+alpha1vpar[1]*f1[5]+alpha1y[0]*f1[3]+alpha1vpar[0]*f1[2]); 
  out[8] += 0.4330127018922193*(alpha1x[1]*f1[8]+alpha1x[0]*f1[4]); 
  out[9] += 0.4330127018922193*(alpha1y[5]*f1[12]+alpha1y[2]*f1[9]+alpha1y[1]*f1[8]+f1[1]*alpha1y[8]+alpha1y[0]*f1[4]+f1[0]*alpha1y[4]); 
  out[10] += 0.4330127018922193*(alpha1vpar[6]*f1[13]+alpha1vpar[3]*f1[10]+alpha1vpar[1]*f1[8]+f1[1]*alpha1vpar[8]+alpha1vpar[0]*f1[4]+f1[0]*alpha1vpar[4]); 
  out[11] += 0.4330127018922193*(alpha1y[4]*f1[13]+alpha1vpar[4]*f1[12]+(alpha1vpar[3]+alpha1y[2]+alpha1x[1])*f1[11]+alpha1y[8]*f1[10]+alpha1vpar[8]*f1[9]+(alpha1vpar[6]+alpha1y[5]+alpha1x[0])*f1[7]+alpha1y[0]*f1[6]+alpha1vpar[0]*f1[5]+alpha1y[1]*f1[3]+alpha1vpar[1]*f1[2]); 
  out[12] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f1[12]+(alpha1y[5]+alpha1x[0])*f1[9]+alpha1y[0]*f1[8]+f1[0]*alpha1y[8]+alpha1y[1]*f1[4]+f1[1]*alpha1y[4]); 
  out[13] += 0.4330127018922193*((alpha1vpar[3]+alpha1x[1])*f1[13]+(alpha1vpar[6]+alpha1x[0])*f1[10]+alpha1vpar[0]*f1[8]+f1[0]*alpha1vpar[8]+alpha1vpar[1]*f1[4]+f1[1]*alpha1vpar[4]); 
  out[14] += 0.4330127018922193*((alpha1vpar[6]+alpha1y[5])*f1[15]+(alpha1vpar[3]+alpha1y[2])*f1[14]+alpha1y[1]*f1[13]+alpha1vpar[1]*f1[12]+alpha1y[0]*f1[10]+alpha1vpar[0]*f1[9]+f1[6]*alpha1y[8]+f1[5]*alpha1vpar[8]+f1[3]*alpha1y[4]+f1[2]*alpha1vpar[4]); 
  out[15] += 0.4330127018922193*((alpha1vpar[3]+alpha1y[2]+alpha1x[1])*f1[15]+(alpha1vpar[6]+alpha1y[5]+alpha1x[0])*f1[14]+alpha1y[0]*f1[13]+alpha1vpar[0]*f1[12]+alpha1y[1]*f1[10]+alpha1vpar[1]*f1[9]+f1[3]*alpha1y[8]+f1[2]*alpha1vpar[8]+alpha1y[4]*f1[6]+alpha1vpar[4]*f1[5]); 
  return cflFreq; 
} 
double EmDeltaFGyrokineticGenGeoVol2x2vSerP1_Bvarsz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f0, const double *f1, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
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

  double hamil0[16]; 
  hamil0[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu)+m_))/rdvpar2Sq; 
  hamil0[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil0[4] = (1.154700538379252*bmag[0])/rdmu2; 

  double hamil1[16]; 
  hamil1[0] = 2.0*phi[0]*q_; 
  hamil1[1] = 2.0*phi[1]*q_; 
  hamil1[2] = 2.0*phi[2]*q_; 
  hamil1[5] = 2.0*phi[3]*q_; 

  double BstarXdBmag[16]; 

  double BstarYdBmag[16]; 

  double BstarXdBmagEM[16]; 
  BstarXdBmagEM[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2; 
  BstarXdBmagEM[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2; 

  double BstarYdBmagEM[16]; 
  BstarYdBmagEM[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2; 
  BstarYdBmagEM[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[16]; 
  double alpha1x[16]; 
  alpha1x[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmagEM[0]*hamil0[3]*rdvpar2)/m_-(0.25*b_z[0]*jacobTotInv[0]*hamil1[2]*rdy2)/q_); 
  alpha1x[1] = 1.732050807568877*rdx2*((0.25*BstarXdBmagEM[1]*hamil0[3]*rdvpar2)/m_-(0.25*b_z[0]*jacobTotInv[0]*hamil1[5]*rdy2)/q_); 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0y[16]; 
  double alpha1y[16]; 
  alpha1y[0] = 0.4330127018922193*((b_z[0]*jacobTotInv[0]*hamil1[1]*rdx2)/q_+(BstarYdBmagEM[0]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha1y[2] = 0.4330127018922193*((b_z[0]*jacobTotInv[0]*hamil1[5]*rdx2)/q_+(BstarYdBmagEM[2]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1y[0]-0.4330127018922193*alpha1y[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1y[2]+0.25*alpha1y[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[16]; 
  double alpha1vpar[16]; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.4330127018922193*(alpha1x[1]*f0[1]+alpha1x[0]*f0[0]); 
  out[2] += 0.4330127018922193*(alpha1y[2]*f0[2]+alpha1y[0]*f0[0]); 
  out[5] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f0[5]+alpha1x[0]*f0[2]+alpha1y[0]*f0[1]); 
  out[6] += 0.4330127018922193*(alpha1x[1]*f0[6]+alpha1x[0]*f0[3]); 
  out[7] += 0.4330127018922193*(alpha1y[2]*f0[7]+alpha1y[0]*f0[3]); 
  out[8] += 0.4330127018922193*(alpha1x[1]*f0[8]+alpha1x[0]*f0[4]); 
  out[9] += 0.4330127018922193*(alpha1y[2]*f0[9]+alpha1y[0]*f0[4]); 
  out[11] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f0[11]+alpha1x[0]*f0[7]+alpha1y[0]*f0[6]); 
  out[12] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f0[12]+alpha1x[0]*f0[9]+alpha1y[0]*f0[8]); 
  out[13] += 0.4330127018922193*(alpha1x[1]*f0[13]+alpha1x[0]*f0[10]); 
  out[14] += 0.4330127018922193*(alpha1y[2]*f0[14]+alpha1y[0]*f0[10]); 
  out[15] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f0[15]+alpha1x[0]*f0[14]+alpha1y[0]*f0[13]); 
  out[1] += 0.4330127018922193*(alpha1x[1]*f1[1]+alpha1x[0]*f1[0]); 
  out[2] += 0.4330127018922193*(alpha1y[2]*f1[2]+alpha1y[0]*f1[0]); 
  out[5] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f1[5]+alpha1x[0]*f1[2]+alpha1y[0]*f1[1]); 
  out[6] += 0.4330127018922193*(alpha1x[1]*f1[6]+alpha1x[0]*f1[3]); 
  out[7] += 0.4330127018922193*(alpha1y[2]*f1[7]+alpha1y[0]*f1[3]); 
  out[8] += 0.4330127018922193*(alpha1x[1]*f1[8]+alpha1x[0]*f1[4]); 
  out[9] += 0.4330127018922193*(alpha1y[2]*f1[9]+alpha1y[0]*f1[4]); 
  out[11] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f1[11]+alpha1x[0]*f1[7]+alpha1y[0]*f1[6]); 
  out[12] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f1[12]+alpha1x[0]*f1[9]+alpha1y[0]*f1[8]); 
  out[13] += 0.4330127018922193*(alpha1x[1]*f1[13]+alpha1x[0]*f1[10]); 
  out[14] += 0.4330127018922193*(alpha1y[2]*f1[14]+alpha1y[0]*f1[10]); 
  out[15] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f1[15]+alpha1x[0]*f1[14]+alpha1y[0]*f1[13]); 
  return cflFreq; 
} 
double EmDeltaFGyrokineticGenGeoVol2x2vSerP1_Bvarsxz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double* dApardt, const double *f0, const double *f1, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
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

  double hamil0[16]; 
  hamil0[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu)+m_))/rdvpar2Sq; 
  hamil0[1] = 2.0*bmag[1]*wmu; 
  hamil0[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil0[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil0[8] = (1.154700538379252*bmag[1])/rdmu2; 

  double hamil1[16]; 
  hamil1[0] = 2.0*phi[0]*q_; 
  hamil1[1] = 2.0*phi[1]*q_; 
  hamil1[2] = 2.0*phi[2]*q_; 
  hamil1[5] = 2.0*phi[3]*q_; 

  double BstarXdBmag[16]; 

  double BstarYdBmag[16]; 
  BstarYdBmag[0] = -(1.732050807568877*jacobTotInv[0]*b_z[1]*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[1] = -(1.732050807568877*b_z[1]*jacobTotInv[1]*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 

  double BstarXdBmagEM[16]; 
  BstarXdBmagEM[0] = 0.8660254037844386*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])*rdy2; 
  BstarXdBmagEM[1] = 0.1732050807568877*((9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])*rdy2; 

  double BstarYdBmagEM[16]; 
  BstarYdBmagEM[0] = -0.8660254037844386*(2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*rdx2; 
  BstarYdBmagEM[1] = -0.8660254037844386*((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])*rdx2; 
  BstarYdBmagEM[2] = -0.8660254037844386*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])*rdx2; 
  BstarYdBmagEM[5] = -0.8660254037844386*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])*rdx2; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[16]; 
  double alpha1x[16]; 
  alpha1x[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmagEM[0]*hamil0[3]*rdvpar2)/m_-(0.25*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil1[5]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil1[2])*rdy2)/q_); 
  alpha1x[1] = 1.732050807568877*rdx2*((((-0.25*(b_z[0]*jacobTotInv[0]*hamil1[5]+(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil1[2]))-0.45*b_z[1]*jacobTotInv[1]*hamil1[5])*rdy2)/q_+(0.25*BstarXdBmagEM[1]*hamil0[3]*rdvpar2)/m_); 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1x[0]-0.4330127018922193*alpha1x[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alpha1x[1]+0.25*alpha1x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0y[16]; 
  alpha0y[0] = 0.4330127018922193*((hamil0[1]*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*rdx2)/q_+(BstarYdBmag[0]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha0y[1] = 0.4330127018922193*((hamil0[1]*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*rdx2)/q_+(BstarYdBmag[1]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha0y[3] = (0.4330127018922193*BstarYdBmag[3]*hamil0[3]*rdvpar2*rdy2)/m_; 
  alpha0y[4] = (0.4330127018922193*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil0[8]*rdx2*rdy2)/q_; 
  alpha0y[6] = (0.4330127018922193*hamil0[3]*BstarYdBmag[6]*rdvpar2*rdy2)/m_; 
  alpha0y[8] = (0.4330127018922193*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil0[8]*rdx2*rdy2)/q_; 
  double alpha1y[16]; 
  alpha1y[0] = 0.4330127018922193*(((hamil1[1]+hamil0[1])*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*rdx2)/q_+(BstarYdBmagEM[0]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha1y[1] = 0.4330127018922193*(((hamil1[1]+hamil0[1])*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*rdx2)/q_+(BstarYdBmagEM[1]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha1y[2] = 0.4330127018922193*(((b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil1[5]*rdx2)/q_+(BstarYdBmagEM[2]*hamil0[3]*rdvpar2)/m_)*rdy2; 
  alpha1y[4] = (0.4330127018922193*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil0[8]*rdx2*rdy2)/q_; 
  alpha1y[5] = 0.4330127018922193*(((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil1[5]*rdx2)/q_+(hamil0[3]*BstarYdBmagEM[5]*rdvpar2)/m_)*rdy2; 
  alpha1y[8] = (0.4330127018922193*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil0[8]*rdx2*rdy2)/q_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6])+0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])-0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6]))-0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])-0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alpha1y[8]+alpha0y[8])-0.25*alpha0y[6]+0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4])+0.25*alpha0y[3]-0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]))+0.25*alpha0y[6]-0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4])+0.25*alpha0y[3]-0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]))+0.25*alpha0y[6]+0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4])-0.25*alpha0y[3]-0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alpha1y[8]+alpha0y[8])-0.25*alpha0y[6]-0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4])-0.25*alpha0y[3]-0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6]))+0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])-0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6])-0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])-0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6])-0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])+0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6]))+0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])+0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alpha1y[8]+alpha0y[8])-0.25*alpha0y[6]-0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4])+0.25*alpha0y[3]+0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]))+0.25*alpha0y[6]+0.4330127018922193*alpha1y[5]-0.25*(alpha1y[4]+alpha0y[4])+0.25*alpha0y[3]+0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]))+0.25*alpha0y[6]-0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4])-0.25*alpha0y[3]+0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alpha1y[8]+alpha0y[8])-0.25*alpha0y[6]+0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4])-0.25*alpha0y[3]+0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6]))-0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])+0.4330127018922193*alpha1y[2]-0.25*(alpha1y[1]+alpha0y[1])+0.25*(alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*(alpha1y[8]+alpha0y[8]+alpha0y[6])+0.4330127018922193*alpha1y[5]+0.25*(alpha1y[4]+alpha0y[4]+alpha0y[3])+0.4330127018922193*alpha1y[2]+0.25*(alpha1y[1]+alpha0y[1]+alpha1y[0]+alpha0y[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[16]; 
  double alpha1vpar[16]; 
  alpha1vpar[0] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[1]*hamil1[5]+BstarYdBmag[0]*hamil1[2])*rdy2+BstarXdBmagEM[0]*hamil0[1]*rdx2))/m_; 
  alpha1vpar[1] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[0]*hamil1[5]+BstarYdBmag[1]*hamil1[2])*rdy2+BstarXdBmagEM[1]*hamil0[1]*rdx2))/m_; 
  alpha1vpar[3] = -(0.4330127018922193*(hamil1[5]*BstarYdBmag[6]+hamil1[2]*BstarYdBmag[3])*rdvpar2*rdy2)/m_; 
  alpha1vpar[4] = -(0.4330127018922193*BstarXdBmagEM[0]*hamil0[8]*rdvpar2*rdx2)/m_; 
  alpha1vpar[6] = -(0.4330127018922193*(hamil1[2]*BstarYdBmag[6]+BstarYdBmag[3]*hamil1[5])*rdvpar2*rdy2)/m_; 
  alpha1vpar[8] = -(0.4330127018922193*BstarXdBmagEM[1]*hamil0[8]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.25*alpha1vpar[8]+0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alpha1vpar[8])-0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1vpar[8]+0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alpha1vpar[8])-0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alpha1vpar[8])+0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1vpar[8]-0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.25*alpha1vpar[8])+0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.25*alpha1vpar[8]-0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]-0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.25*alpha1vpar[8]-0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alpha1vpar[8])+0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alpha1vpar[8]-0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alpha1vpar[8])+0.4330127018922193*alpha1vpar[6]-0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alpha1vpar[8])-0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alpha1vpar[8]+0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.25*alpha1vpar[8])-0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]-0.25*alpha1vpar[1]+0.25*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.25*alpha1vpar[8]+0.4330127018922193*alpha1vpar[6]+0.25*alpha1vpar[4]+0.4330127018922193*alpha1vpar[3]+0.25*(alpha1vpar[1]+alpha1vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.4330127018922193*(alpha1x[1]*f0[1]+alpha1x[0]*f0[0]); 
  out[2] += 0.4330127018922193*(alpha0y[8]*f1[8]+alpha1y[8]*f0[8]+alpha0y[6]*f1[6]+alpha1y[5]*f0[5]+alpha0y[4]*f1[4]+alpha1y[4]*f0[4]+alpha0y[3]*f1[3]+alpha1y[2]*f0[2]+alpha0y[1]*f1[1]+alpha1y[1]*f0[1]+alpha0y[0]*f1[0]+alpha1y[0]*f0[0]); 
  out[3] += 0.4330127018922193*(alpha1vpar[8]*f0[8]+alpha1vpar[6]*f0[6]+alpha1vpar[4]*f0[4]+alpha1vpar[3]*f0[3]+alpha1vpar[1]*f0[1]+alpha1vpar[0]*f0[0]); 
  out[5] += 0.4330127018922193*(alpha0y[4]*f1[8]+alpha1y[4]*f0[8]+f0[4]*alpha1y[8]+f1[4]*alpha0y[8]+alpha0y[3]*f1[6]+f1[3]*alpha0y[6]+(alpha1y[2]+alpha1x[1])*f0[5]+f0[2]*(alpha1y[5]+alpha1x[0])+alpha0y[0]*f1[1]+alpha1y[0]*f0[1]+f0[0]*alpha1y[1]+f1[0]*alpha0y[1]); 
  out[6] += 0.4330127018922193*(alpha1vpar[4]*f0[8]+f0[4]*alpha1vpar[8]+(alpha1vpar[3]+alpha1x[1])*f0[6]+f0[3]*(alpha1vpar[6]+alpha1x[0])+alpha1vpar[0]*f0[1]+f0[0]*alpha1vpar[1]); 
  out[7] += 0.4330127018922193*(alpha0y[8]*f1[13]+alpha1y[8]*f0[13]+alpha1vpar[8]*f0[12]+(alpha1vpar[6]+alpha1y[5])*f0[11]+alpha0y[4]*f1[10]+alpha1y[4]*f0[10]+alpha1vpar[4]*f0[9]+(alpha1vpar[3]+alpha1y[2])*f0[7]+alpha0y[1]*f1[6]+alpha1y[1]*f0[6]+f1[1]*alpha0y[6]+alpha1vpar[1]*f0[5]+alpha0y[0]*f1[3]+alpha1y[0]*f0[3]+f1[0]*alpha0y[3]+alpha1vpar[0]*f0[2]); 
  out[8] += 0.4330127018922193*(alpha1x[1]*f0[8]+alpha1x[0]*f0[4]); 
  out[9] += 0.4330127018922193*(alpha0y[6]*f1[13]+alpha1y[5]*f0[12]+alpha0y[3]*f1[10]+alpha1y[2]*f0[9]+alpha0y[1]*f1[8]+alpha1y[1]*f0[8]+f0[1]*alpha1y[8]+f1[1]*alpha0y[8]+alpha0y[0]*f1[4]+alpha1y[0]*f0[4]+f0[0]*alpha1y[4]+f1[0]*alpha0y[4]); 
  out[10] += 0.4330127018922193*(alpha1vpar[6]*f0[13]+alpha1vpar[3]*f0[10]+alpha1vpar[1]*f0[8]+f0[1]*alpha1vpar[8]+alpha1vpar[0]*f0[4]+f0[0]*alpha1vpar[4]); 
  out[11] += 0.4330127018922193*(alpha0y[4]*f1[13]+alpha1y[4]*f0[13]+alpha1vpar[4]*f0[12]+(alpha1vpar[3]+alpha1y[2]+alpha1x[1])*f0[11]+alpha0y[8]*f1[10]+alpha1y[8]*f0[10]+alpha1vpar[8]*f0[9]+(alpha1vpar[6]+alpha1y[5]+alpha1x[0])*f0[7]+alpha0y[0]*f1[6]+alpha1y[0]*f0[6]+f1[0]*alpha0y[6]+alpha1vpar[0]*f0[5]+alpha0y[1]*f1[3]+alpha1y[1]*f0[3]+f1[1]*alpha0y[3]+alpha1vpar[1]*f0[2]); 
  out[12] += 0.4330127018922193*(alpha0y[3]*f1[13]+(alpha1y[2]+alpha1x[1])*f0[12]+alpha0y[6]*f1[10]+(alpha1y[5]+alpha1x[0])*f0[9]+alpha0y[0]*f1[8]+alpha1y[0]*f0[8]+f0[0]*alpha1y[8]+f1[0]*alpha0y[8]+alpha0y[1]*f1[4]+alpha1y[1]*f0[4]+f0[1]*alpha1y[4]+f1[1]*alpha0y[4]); 
  out[13] += 0.4330127018922193*((alpha1vpar[3]+alpha1x[1])*f0[13]+(alpha1vpar[6]+alpha1x[0])*f0[10]+alpha1vpar[0]*f0[8]+f0[0]*alpha1vpar[8]+alpha1vpar[1]*f0[4]+f0[1]*alpha1vpar[4]); 
  out[14] += 0.4330127018922193*((alpha1vpar[6]+alpha1y[5])*f0[15]+(alpha1vpar[3]+alpha1y[2])*f0[14]+alpha0y[1]*f1[13]+alpha1y[1]*f0[13]+alpha1vpar[1]*f0[12]+alpha0y[0]*f1[10]+alpha1y[0]*f0[10]+alpha1vpar[0]*f0[9]+alpha0y[6]*f1[8]+f0[6]*alpha1y[8]+f0[5]*alpha1vpar[8]+f1[6]*alpha0y[8]+alpha0y[3]*f1[4]+f0[3]*alpha1y[4]+f0[2]*alpha1vpar[4]+f1[3]*alpha0y[4]); 
  out[15] += 0.4330127018922193*((alpha1vpar[3]+alpha1y[2]+alpha1x[1])*f0[15]+(alpha1vpar[6]+alpha1y[5]+alpha1x[0])*f0[14]+alpha0y[0]*f1[13]+alpha1y[0]*f0[13]+alpha1vpar[0]*f0[12]+alpha0y[1]*f1[10]+alpha1y[1]*f0[10]+alpha1vpar[1]*f0[9]+alpha0y[3]*f1[8]+f0[3]*alpha1y[8]+f0[2]*alpha1vpar[8]+f1[3]*alpha0y[8]+alpha0y[4]*f1[6]+alpha1y[4]*f0[6]+f1[4]*alpha0y[6]+alpha1vpar[4]*f0[5]); 
  out[1] += 0.4330127018922193*(alpha1x[1]*f1[1]+alpha1x[0]*f1[0]); 
  out[2] += 0.4330127018922193*(alpha1y[8]*f1[8]+alpha1y[5]*f1[5]+alpha1y[4]*f1[4]+alpha1y[2]*f1[2]+alpha1y[1]*f1[1]+alpha1y[0]*f1[0]); 
  out[3] += 0.4330127018922193*(alpha1vpar[8]*f1[8]+alpha1vpar[6]*f1[6]+alpha1vpar[4]*f1[4]+alpha1vpar[3]*f1[3]+alpha1vpar[1]*f1[1]+alpha1vpar[0]*f1[0]); 
  out[5] += 0.4330127018922193*(alpha1y[4]*f1[8]+f1[4]*alpha1y[8]+(alpha1y[2]+alpha1x[1])*f1[5]+f1[2]*(alpha1y[5]+alpha1x[0])+alpha1y[0]*f1[1]+f1[0]*alpha1y[1]); 
  out[6] += 0.4330127018922193*(alpha1vpar[4]*f1[8]+f1[4]*alpha1vpar[8]+(alpha1vpar[3]+alpha1x[1])*f1[6]+f1[3]*(alpha1vpar[6]+alpha1x[0])+alpha1vpar[0]*f1[1]+f1[0]*alpha1vpar[1]); 
  out[7] += 0.4330127018922193*(alpha1y[8]*f1[13]+alpha1vpar[8]*f1[12]+(alpha1vpar[6]+alpha1y[5])*f1[11]+alpha1y[4]*f1[10]+alpha1vpar[4]*f1[9]+(alpha1vpar[3]+alpha1y[2])*f1[7]+alpha1y[1]*f1[6]+alpha1vpar[1]*f1[5]+alpha1y[0]*f1[3]+alpha1vpar[0]*f1[2]); 
  out[8] += 0.4330127018922193*(alpha1x[1]*f1[8]+alpha1x[0]*f1[4]); 
  out[9] += 0.4330127018922193*(alpha1y[5]*f1[12]+alpha1y[2]*f1[9]+alpha1y[1]*f1[8]+f1[1]*alpha1y[8]+alpha1y[0]*f1[4]+f1[0]*alpha1y[4]); 
  out[10] += 0.4330127018922193*(alpha1vpar[6]*f1[13]+alpha1vpar[3]*f1[10]+alpha1vpar[1]*f1[8]+f1[1]*alpha1vpar[8]+alpha1vpar[0]*f1[4]+f1[0]*alpha1vpar[4]); 
  out[11] += 0.4330127018922193*(alpha1y[4]*f1[13]+alpha1vpar[4]*f1[12]+(alpha1vpar[3]+alpha1y[2]+alpha1x[1])*f1[11]+alpha1y[8]*f1[10]+alpha1vpar[8]*f1[9]+(alpha1vpar[6]+alpha1y[5]+alpha1x[0])*f1[7]+alpha1y[0]*f1[6]+alpha1vpar[0]*f1[5]+alpha1y[1]*f1[3]+alpha1vpar[1]*f1[2]); 
  out[12] += 0.4330127018922193*((alpha1y[2]+alpha1x[1])*f1[12]+(alpha1y[5]+alpha1x[0])*f1[9]+alpha1y[0]*f1[8]+f1[0]*alpha1y[8]+alpha1y[1]*f1[4]+f1[1]*alpha1y[4]); 
  out[13] += 0.4330127018922193*((alpha1vpar[3]+alpha1x[1])*f1[13]+(alpha1vpar[6]+alpha1x[0])*f1[10]+alpha1vpar[0]*f1[8]+f1[0]*alpha1vpar[8]+alpha1vpar[1]*f1[4]+f1[1]*alpha1vpar[4]); 
  out[14] += 0.4330127018922193*((alpha1vpar[6]+alpha1y[5])*f1[15]+(alpha1vpar[3]+alpha1y[2])*f1[14]+alpha1y[1]*f1[13]+alpha1vpar[1]*f1[12]+alpha1y[0]*f1[10]+alpha1vpar[0]*f1[9]+f1[6]*alpha1y[8]+f1[5]*alpha1vpar[8]+f1[3]*alpha1y[4]+f1[2]*alpha1vpar[4]); 
  out[15] += 0.4330127018922193*((alpha1vpar[3]+alpha1y[2]+alpha1x[1])*f1[15]+(alpha1vpar[6]+alpha1y[5]+alpha1x[0])*f1[14]+alpha1y[0]*f1[13]+alpha1vpar[0]*f1[12]+alpha1y[1]*f1[10]+alpha1vpar[1]*f1[9]+f1[3]*alpha1y[8]+f1[2]*alpha1vpar[8]+alpha1y[4]*f1[6]+alpha1vpar[4]*f1[5]); 
  return cflFreq; 
} 
double EmDeltaFGyrokineticGenGeoStep2Vol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f0, const double *f1, double *out) 
{ 
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // dApardt: time derivative of Apar.
  // f: Distribution function.
  // out: output increment.

  double rdvpar2 = 2.0/dxv[2]; 
  out[3] += -(0.8660254037844386*(dApardt[3]*f0[5]+dApardt[2]*f0[2]+dApardt[1]*f0[1]+dApardt[0]*f0[0])*q_*rdvpar2)/m_; 
  out[6] += -(0.8660254037844386*(dApardt[2]*f0[5]+f0[2]*dApardt[3]+dApardt[0]*f0[1]+f0[0]*dApardt[1])*q_*rdvpar2)/m_; 
  out[7] += -(0.8660254037844386*(dApardt[1]*f0[5]+f0[1]*dApardt[3]+dApardt[0]*f0[2]+f0[0]*dApardt[2])*q_*rdvpar2)/m_; 
  out[10] += -(0.8660254037844386*(dApardt[3]*f0[12]+dApardt[2]*f0[9]+dApardt[1]*f0[8]+dApardt[0]*f0[4])*q_*rdvpar2)/m_; 
  out[11] += -(0.8660254037844386*(dApardt[0]*f0[5]+f0[0]*dApardt[3]+dApardt[1]*f0[2]+f0[1]*dApardt[2])*q_*rdvpar2)/m_; 
  out[13] += -(0.8660254037844386*(dApardt[2]*f0[12]+dApardt[3]*f0[9]+dApardt[0]*f0[8]+dApardt[1]*f0[4])*q_*rdvpar2)/m_; 
  out[14] += -(0.8660254037844386*(dApardt[1]*f0[12]+dApardt[0]*f0[9]+dApardt[3]*f0[8]+dApardt[2]*f0[4])*q_*rdvpar2)/m_; 
  out[15] += -(0.8660254037844386*(dApardt[0]*f0[12]+dApardt[1]*f0[9]+dApardt[2]*f0[8]+dApardt[3]*f0[4])*q_*rdvpar2)/m_; 
  out[3] += -(0.8660254037844386*(dApardt[3]*f1[5]+dApardt[2]*f1[2]+dApardt[1]*f1[1]+dApardt[0]*f1[0])*q_*rdvpar2)/m_; 
  out[6] += -(0.8660254037844386*(dApardt[2]*f1[5]+f1[2]*dApardt[3]+dApardt[0]*f1[1]+f1[0]*dApardt[1])*q_*rdvpar2)/m_; 
  out[7] += -(0.8660254037844386*(dApardt[1]*f1[5]+f1[1]*dApardt[3]+dApardt[0]*f1[2]+f1[0]*dApardt[2])*q_*rdvpar2)/m_; 
  out[10] += -(0.8660254037844386*(dApardt[3]*f1[12]+dApardt[2]*f1[9]+dApardt[1]*f1[8]+dApardt[0]*f1[4])*q_*rdvpar2)/m_; 
  out[11] += -(0.8660254037844386*(dApardt[0]*f1[5]+f1[0]*dApardt[3]+dApardt[1]*f1[2]+f1[1]*dApardt[2])*q_*rdvpar2)/m_; 
  out[13] += -(0.8660254037844386*(dApardt[2]*f1[12]+dApardt[3]*f1[9]+dApardt[0]*f1[8]+dApardt[1]*f1[4])*q_*rdvpar2)/m_; 
  out[14] += -(0.8660254037844386*(dApardt[1]*f1[12]+dApardt[0]*f1[9]+dApardt[3]*f1[8]+dApardt[2]*f1[4])*q_*rdvpar2)/m_; 
  out[15] += -(0.8660254037844386*(dApardt[0]*f1[12]+dApardt[1]*f1[9]+dApardt[2]*f1[8]+dApardt[3]*f1[4])*q_*rdvpar2)/m_; 
  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
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

return cflFreq; 
} 
