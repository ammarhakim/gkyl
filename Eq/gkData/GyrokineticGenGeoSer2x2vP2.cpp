#include <GyrokineticModDecl.h> 
double GyrokineticGenGeoVol2x2vSerP2_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f, double *out) 
{ 
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

  double hamil[48]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[11] = 2.0*phi[4]*q_; 
  hamil[12] = 2.0*phi[5]*q_; 
  hamil[13] = (0.5962847939999438*m_)/rdvpar2Sq; 
  hamil[19] = 2.0*phi[6]*q_; 
  hamil[20] = 2.0*phi[7]*q_; 

  double BstarXdBmag[48]; 

  double BstarYdBmag[48]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[48]; 
  alphax[0] = -(0.4330127018922193*b_z[0]*jacobTotInv[0]*hamil[2]*rdx2*rdy2)/q_; 
  alphax[1] = -(0.4330127018922193*b_z[0]*jacobTotInv[0]*hamil[5]*rdx2*rdy2)/q_; 
  alphax[2] = -(0.9682458365518543*b_z[0]*jacobTotInv[0]*hamil[12]*rdx2*rdy2)/q_; 
  alphax[5] = -(0.9682458365518543*b_z[0]*jacobTotInv[0]*hamil[20]*rdx2*rdy2)/q_; 
  alphax[11] = -(0.4330127018922193*b_z[0]*jacobTotInv[0]*hamil[19]*rdx2*rdy2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*(2.23606797749979*alphax[11]-1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.23606797749979*alphax[11]+1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[48]; 
  alphay[0] = (0.4330127018922193*b_z[0]*jacobTotInv[0]*hamil[1]*rdx2*rdy2)/q_; 
  alphay[1] = (0.9682458365518543*b_z[0]*jacobTotInv[0]*hamil[11]*rdx2*rdy2)/q_; 
  alphay[2] = (0.4330127018922193*b_z[0]*jacobTotInv[0]*hamil[5]*rdx2*rdy2)/q_; 
  alphay[5] = (0.9682458365518543*b_z[0]*jacobTotInv[0]*hamil[19]*rdx2*rdy2)/q_; 
  alphay[12] = (0.4330127018922194*b_z[0]*jacobTotInv[0]*hamil[20]*rdx2*rdy2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*(2.23606797749979*alphay[12]-1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.23606797749979*alphay[12]+1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[48]; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
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
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[11]*f[11]+alphax[5]*f[5]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[12]*f[12]+alphay[5]*f[5]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[5] += 0.05*((8.660254037844387*alphay[12]+7.745966692414834*alphax[5])*f[20]+8.660254037844387*alphax[11]*f[19]+7.745966692414834*(alphay[5]*f[19]+alphax[2]*f[12]+alphay[1]*f[11])+8.660254037844386*((alphay[2]+alphax[1])*f[5]+f[2]*alphay[5]+f[1]*alphax[5]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1])); 
  out[6] += 0.05*(8.660254037844387*alphax[11]*f[21]+8.660254037844386*(alphax[5]*f[15]+alphax[2]*f[7]+alphax[1]*f[6]+alphax[0]*f[3])); 
  out[7] += 0.05*(8.660254037844387*alphay[12]*f[22]+8.660254037844386*(alphay[5]*f[15]+alphay[2]*f[7]+alphay[1]*f[6]+alphay[0]*f[3])); 
  out[8] += 0.05*(8.660254037844387*alphax[11]*f[25]+8.660254037844386*(alphax[5]*f[16]+alphax[2]*f[9]+alphax[1]*f[8]+alphax[0]*f[4])); 
  out[9] += 0.05*(8.660254037844387*alphay[12]*f[26]+8.660254037844386*(alphay[5]*f[16]+alphay[2]*f[9]+alphay[1]*f[8]+alphay[0]*f[4])); 
  out[11] += 0.05*(17.32050807568877*alphax[5]*f[19]+17.32050807568877*(alphax[1]*f[11]+f[1]*alphax[11])+19.36491673103709*(alphax[2]*f[5]+f[2]*alphax[5]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[12] += 0.05*(17.32050807568877*alphay[5]*f[20]+17.32050807568877*(alphay[2]*f[12]+f[2]*alphay[12])+19.36491673103709*(alphay[1]*f[5]+f[1]*alphay[5]+alphay[0]*f[2]+f[0]*alphay[2])); 
  out[15] += 0.05*((8.660254037844386*alphay[12]+7.745966692414834*alphax[5])*f[33]+8.660254037844386*alphax[11]*f[32]+7.745966692414834*(alphay[5]*f[32]+alphax[2]*f[22]+alphay[1]*f[21])+8.660254037844386*((alphay[2]+alphax[1])*f[15]+(alphay[5]+alphax[0])*f[7]+(alphax[5]+alphay[0])*f[6]+(alphax[2]+alphay[1])*f[3])); 
  out[16] += 0.05*((8.660254037844386*alphay[12]+7.745966692414834*alphax[5])*f[36]+8.660254037844386*alphax[11]*f[35]+7.745966692414834*(alphay[5]*f[35]+alphax[2]*f[26]+alphay[1]*f[25])+8.660254037844386*((alphay[2]+alphax[1])*f[16]+(alphay[5]+alphax[0])*f[9]+(alphax[5]+alphay[0])*f[8]+(alphax[2]+alphay[1])*f[4])); 
  out[17] += 0.4330127018922193*(alphax[11]*f[37]+alphax[5]*f[31]+alphax[2]*f[18]+alphax[1]*f[17]+alphax[0]*f[10]); 
  out[18] += 0.4330127018922193*(alphay[12]*f[38]+alphay[5]*f[31]+alphay[2]*f[18]+alphay[1]*f[17]+alphay[0]*f[10]); 
  out[19] += 0.05*(17.32050807568877*alphax[2]*f[20]+(8.660254037844386*alphay[2]+17.32050807568877*alphax[1])*f[19]+17.32050807568877*alphax[5]*f[12]+(17.32050807568877*alphax[5]+8.660254037844387*alphay[0])*f[11]+f[5]*(17.32050807568877*alphax[11]+7.745966692414834*alphay[5])+19.36491673103708*(alphax[0]*f[5]+f[0]*alphax[5]+alphax[1]*f[2])+f[1]*(19.36491673103708*alphax[2]+7.745966692414834*alphay[1])); 
  out[20] += 0.05*((17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[20]+17.32050807568877*alphay[1]*f[19]+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[12]+17.32050807568877*(f[5]*alphay[12]+alphay[5]*f[11])+7.745966692414834*alphax[5]*f[5]+19.36491673103708*(alphay[0]*f[5]+f[0]*alphay[5])+7.745966692414834*alphax[2]*f[2]+19.36491673103708*(alphay[1]*f[2]+f[1]*alphay[2])); 
  out[21] += 0.05*(17.32050807568877*alphax[5]*f[32]+17.32050807568877*alphax[1]*f[21]+19.36491673103708*alphax[2]*f[15]+17.32050807568877*f[6]*alphax[11]+19.36491673103708*(alphax[5]*f[7]+alphax[0]*f[6]+alphax[1]*f[3])); 
  out[22] += 0.05*(17.32050807568877*alphay[5]*f[33]+17.32050807568877*alphay[2]*f[22]+19.36491673103708*alphay[1]*f[15]+17.32050807568877*f[7]*alphay[12]+19.36491673103708*(alphay[0]*f[7]+alphay[5]*f[6]+alphay[2]*f[3])); 
  out[23] += 0.05*(8.660254037844387*alphax[5]*f[34]+8.660254037844386*(alphax[2]*f[24]+alphax[1]*f[23])+8.660254037844387*alphax[0]*f[13]); 
  out[24] += 0.05*(8.660254037844387*alphay[5]*f[34]+8.660254037844386*(alphay[2]*f[24]+alphay[1]*f[23])+8.660254037844387*alphay[0]*f[13]); 
  out[25] += 0.05*(17.32050807568877*alphax[5]*f[35]+17.32050807568877*alphax[1]*f[25]+19.36491673103708*alphax[2]*f[16]+17.32050807568877*f[8]*alphax[11]+19.36491673103708*(alphax[5]*f[9]+alphax[0]*f[8]+alphax[1]*f[4])); 
  out[26] += 0.05*(17.32050807568877*alphay[5]*f[36]+17.32050807568877*alphay[2]*f[26]+19.36491673103708*alphay[1]*f[16]+17.32050807568877*f[9]*alphay[12]+19.36491673103708*(alphay[0]*f[9]+alphay[5]*f[8]+alphay[2]*f[4])); 
  out[28] += 0.05*(8.660254037844387*alphax[5]*f[41]+8.660254037844386*(alphax[2]*f[29]+alphax[1]*f[28])+8.660254037844387*alphax[0]*f[14]); 
  out[29] += 0.05*(8.660254037844387*alphay[5]*f[41]+8.660254037844386*(alphay[2]*f[29]+alphay[1]*f[28])+8.660254037844387*alphay[0]*f[14]); 
  out[31] += 0.05*((8.660254037844387*alphay[12]+7.745966692414834*alphax[5])*f[45]+8.660254037844387*alphax[11]*f[44]+7.745966692414834*(alphay[5]*f[44]+alphax[2]*f[38]+alphay[1]*f[37])+8.660254037844386*((alphay[2]+alphax[1])*f[31]+(alphay[5]+alphax[0])*f[18]+(alphax[5]+alphay[0])*f[17]+(alphax[2]+alphay[1])*f[10])); 
  out[32] += 0.05*(17.32050807568877*alphax[2]*f[33]+(8.660254037844386*alphay[2]+17.32050807568877*alphax[1])*f[32]+17.32050807568877*alphax[5]*f[22]+(17.32050807568877*alphax[5]+8.660254037844387*alphay[0])*f[21]+(17.32050807568877*alphax[11]+7.745966692414834*alphay[5])*f[15]+19.36491673103709*(alphax[0]*f[15]+alphax[1]*f[7])+(19.36491673103709*alphax[2]+7.745966692414834*alphay[1])*f[6]+19.36491673103709*f[3]*alphax[5]); 
  out[33] += 0.05*((17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[33]+17.32050807568877*alphay[1]*f[32]+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[22]+17.32050807568877*alphay[5]*f[21]+(17.32050807568877*alphay[12]+7.745966692414834*alphax[5]+19.36491673103709*alphay[0])*f[15]+7.745966692414834*alphax[2]*f[7]+19.36491673103709*(alphay[1]*f[7]+alphay[2]*f[6]+f[3]*alphay[5])); 
  out[34] += 0.05*(8.660254037844386*(alphay[2]+alphax[1])*f[34]+8.660254037844387*((alphay[5]+alphax[0])*f[24]+(alphax[5]+alphay[0])*f[23])+8.660254037844386*(alphax[2]+alphay[1])*f[13]); 
  out[35] += 0.05*(17.32050807568877*alphax[2]*f[36]+(8.660254037844386*alphay[2]+17.32050807568877*alphax[1])*f[35]+17.32050807568877*alphax[5]*f[26]+(17.32050807568877*alphax[5]+8.660254037844387*alphay[0])*f[25]+(17.32050807568877*alphax[11]+7.745966692414834*alphay[5])*f[16]+19.36491673103709*(alphax[0]*f[16]+alphax[1]*f[9])+(19.36491673103709*alphax[2]+7.745966692414834*alphay[1])*f[8]+19.36491673103709*f[4]*alphax[5]); 
  out[36] += 0.05*((17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[36]+17.32050807568877*alphay[1]*f[35]+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[26]+17.32050807568877*alphay[5]*f[25]+(17.32050807568877*alphay[12]+7.745966692414834*alphax[5]+19.36491673103709*alphay[0])*f[16]+7.745966692414834*alphax[2]*f[9]+19.36491673103709*(alphay[1]*f[9]+alphay[2]*f[8]+f[4]*alphay[5])); 
  out[37] += 0.05*(17.32050807568877*alphax[5]*f[44]+17.32050807568877*alphax[1]*f[37]+19.36491673103709*(alphax[2]*f[31]+alphax[5]*f[18])+17.32050807568877*alphax[11]*f[17]+19.36491673103709*(alphax[0]*f[17]+alphax[1]*f[10])); 
  out[38] += 0.05*(17.32050807568877*alphay[5]*f[45]+17.32050807568877*alphay[2]*f[38]+19.36491673103709*alphay[1]*f[31]+17.32050807568877*alphay[12]*f[18]+19.36491673103709*(alphay[0]*f[18]+alphay[5]*f[17]+alphay[2]*f[10])); 
  out[39] += 0.05*(8.660254037844387*alphax[5]*f[46]+8.660254037844386*(alphax[2]*f[40]+alphax[1]*f[39])+8.660254037844387*alphax[0]*f[27]); 
  out[40] += 0.05*(8.660254037844387*alphay[5]*f[46]+8.660254037844386*(alphay[2]*f[40]+alphay[1]*f[39])+8.660254037844387*alphay[0]*f[27]); 
  out[41] += 0.05*(8.660254037844386*(alphay[2]+alphax[1])*f[41]+8.660254037844387*((alphay[5]+alphax[0])*f[29]+(alphax[5]+alphay[0])*f[28])+8.660254037844386*(alphax[2]+alphay[1])*f[14]); 
  out[42] += 0.05*(8.660254037844387*alphax[5]*f[47]+8.660254037844386*(alphax[2]*f[43]+alphax[1]*f[42])+8.660254037844387*alphax[0]*f[30]); 
  out[43] += 0.05*(8.660254037844387*alphay[5]*f[47]+8.660254037844386*(alphay[2]*f[43]+alphay[1]*f[42])+8.660254037844387*alphay[0]*f[30]); 
  out[44] += 0.05*(17.32050807568877*alphax[2]*f[45]+(8.660254037844386*alphay[2]+17.32050807568877*alphax[1])*f[44]+17.32050807568877*alphax[5]*f[38]+(17.32050807568877*alphax[5]+8.660254037844387*alphay[0])*f[37]+(17.32050807568877*alphax[11]+7.745966692414834*alphay[5])*f[31]+19.36491673103708*(alphax[0]*f[31]+alphax[1]*f[18])+(19.36491673103708*alphax[2]+7.745966692414834*alphay[1])*f[17]+19.36491673103708*alphax[5]*f[10]); 
  out[45] += 0.05*((17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[45]+17.32050807568877*alphay[1]*f[44]+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[38]+17.32050807568877*alphay[5]*f[37]+(17.32050807568877*alphay[12]+7.745966692414834*alphax[5]+19.36491673103708*alphay[0])*f[31]+7.745966692414834*alphax[2]*f[18]+19.36491673103708*(alphay[1]*f[18]+alphay[2]*f[17]+alphay[5]*f[10])); 
  out[46] += 0.05*(8.660254037844386*(alphay[2]+alphax[1])*f[46]+8.660254037844387*((alphay[5]+alphax[0])*f[40]+(alphax[5]+alphay[0])*f[39])+8.660254037844386*(alphax[2]+alphay[1])*f[27]); 
  out[47] += 0.05*(8.660254037844386*(alphay[2]+alphax[1])*f[47]+8.660254037844387*((alphay[5]+alphax[0])*f[43]+(alphax[5]+alphay[0])*f[42])+8.660254037844386*(alphax[2]+alphay[1])*f[30]); 
  return cflFreq; 
} 
double GyrokineticGenGeoVol2x2vSerP2_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f, double *out) 
{ 
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

  double hamil[48]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[11] = 2.0*(bmag[4]*wmu+phi[4]*q_); 
  hamil[12] = 2.0*phi[5]*q_; 
  hamil[13] = (0.5962847939999438*m_)/rdvpar2Sq; 
  hamil[19] = 2.0*phi[6]*q_; 
  hamil[20] = 2.0*phi[7]*q_; 
  hamil[25] = (1.154700538379251*bmag[4])/rdmu2; 

  double BstarXdBmag[48]; 

  double BstarYdBmag[48]; 
  BstarYdBmag[0] = -(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[1] = -(1.732050807568877*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[3] = -(1.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[6] = -(1.0*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[11] = -(1.732050807568877*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[21] = -(1.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[48]; 
  alphax[0] = ((3.872983346207417*((-0.1118033988749895*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.07142857142857142*b_z[4]*jacobTotInv[4]-0.1*b_z[1]*jacobTotInv[1])*hamil[19]+1.732050807568877*((-0.25*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[5]+hamil[2]*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])))-0.223606797749979*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])*hamil[5]))*rdx2*rdy2)/q_; 
  alphax[1] = ((3.872983346207417*((-0.1756910553749835*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-0.1*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[19]+1.732050807568877*(((-0.223606797749979*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.3928571428571428*b_z[4]*jacobTotInv[4]-0.45*b_z[1]*jacobTotInv[1]-0.25*b_z[0]*jacobTotInv[0])*hamil[5]+hamil[2]*((-0.223606797749979*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-0.25*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))))*rdx2*rdy2)/q_; 
  alphax[2] = (((-0.9682458365518543*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[20]+(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[12]))-0.8660254037844387*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])*hamil[20])*rdx2*rdy2)/q_; 
  alphax[5] = ((3.872983346207417*((-0.223606797749979*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.3928571428571428*b_z[4]*jacobTotInv[4]-0.45*b_z[1]*jacobTotInv[1]-0.25*b_z[0]*jacobTotInv[0])*hamil[20]+1.732050807568877*((-0.5*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-0.5590169943749475*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[12])*rdx2*rdy2)/q_; 
  alphax[11] = ((3.872983346207417*((-0.07142857142857142*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.2395787118749775*b_z[4]*jacobTotInv[4]+0.4472135954999579*((-0.3928571428571428*b_z[1]*jacobTotInv[1])-0.25*b_z[0]*jacobTotInv[0]))*hamil[19]+1.732050807568877*(((-0.3928571428571428*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-0.223606797749979*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[5]+hamil[2]*((-0.25*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.159719141249985*b_z[4]*jacobTotInv[4]-0.223606797749979*b_z[1]*jacobTotInv[1])))*rdx2*rdy2)/q_; 
  alphax[19] = ((1.732050807568877*((-0.8784552768749174*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-0.5*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[20]+3.872983346207417*((-0.25*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.159719141249985*b_z[4]*jacobTotInv[4]-0.223606797749979*b_z[1]*jacobTotInv[1])*hamil[12])*rdx2*rdy2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*(2.23606797749979*alphax[11]-1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.23606797749979*alphax[11]+1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[48]; 
  alphay[0] = 1.732050807568877*(0.25*(((2.23606797749979*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[11]+hamil[1]*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*rdx2)/q_+((2.23606797749979*BstarYdBmag[3]*hamil[13]+BstarYdBmag[0]*hamil[3])*rdvpar2)/m_)+(0.5*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])*hamil[11]*rdx2)/q_)*rdy2; 
  alphay[1] = 1.732050807568877*(0.25*((hamil[1]*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*rdx2)/q_+((2.23606797749979*BstarYdBmag[6]*hamil[13]+BstarYdBmag[1]*hamil[3])*rdvpar2)/m_)+(((0.5*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+0.8784552768749174*b_z[4]*jacobTotInv[4]+0.25*(4.024922359499621*b_z[1]*jacobTotInv[1]+2.23606797749979*b_z[0]*jacobTotInv[0]))*hamil[11]+0.223606797749979*hamil[1]*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))*rdx2)/q_)*rdy2; 
  alphay[2] = ((0.25*(3.872983346207417*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[19]+1.732050807568877*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[5])+0.8660254037844387*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])*hamil[19])*rdx2*rdy2)/q_; 
  alphay[3] = (0.4330127018922193*(2.23606797749979*BstarYdBmag[0]*hamil[13]+BstarYdBmag[3]*hamil[3])*rdvpar2*rdy2)/m_; 
  alphay[4] = ((0.25*(3.872983346207417*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[25]+1.732050807568877*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[8])+0.8660254037844387*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])*hamil[25])*rdx2*rdy2)/q_; 
  alphay[5] = ((3.872983346207417*(0.223606797749979*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+0.3928571428571428*b_z[4]*jacobTotInv[4]+0.45*b_z[1]*jacobTotInv[1]+0.25*b_z[0]*jacobTotInv[0])*hamil[19]+1.732050807568877*(0.223606797749979*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.25*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[5])*rdx2*rdy2)/q_; 
  alphay[6] = (0.4330127018922193*(2.23606797749979*BstarYdBmag[1]*hamil[13]+hamil[3]*BstarYdBmag[6])*rdvpar2*rdy2)/m_; 
  alphay[8] = ((3.872983346207417*(0.223606797749979*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+0.3928571428571428*b_z[4]*jacobTotInv[4]+0.45*b_z[1]*jacobTotInv[1]+0.25*b_z[0]*jacobTotInv[0])*hamil[25]+1.732050807568877*(0.223606797749979*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.25*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[8])*rdx2*rdy2)/q_; 
  alphay[11] = ((1.732050807568877*((0.8784552768749174*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.5*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[11]+hamil[1]*(0.25*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+0.159719141249985*b_z[4]*jacobTotInv[4]+0.223606797749979*b_z[1]*jacobTotInv[1]))*rdx2)/q_+(0.25*(3.872983346207417*hamil[13]*BstarYdBmag[21]+1.732050807568877*hamil[3]*BstarYdBmag[11])*rdvpar2)/m_)*rdy2; 
  alphay[12] = (0.4330127018922194*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[20]*rdx2*rdy2)/q_; 
  alphay[13] = (0.8660254037844386*BstarYdBmag[3]*hamil[13]*rdvpar2*rdy2)/m_; 
  alphay[19] = ((1.732050807568877*(0.8784552768749174*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.5*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[19]+(0.2581988897471611*(1.677050983124842*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+1.5*b_z[1]*jacobTotInv[1])+0.276641667586244*b_z[4]*jacobTotInv[4])*hamil[5])*rdx2*rdy2)/q_; 
  alphay[20] = (1.732050807568877*(0.223606797749979*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.25*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[20]*rdx2*rdy2)/q_; 
  alphay[21] = (0.25*(1.732050807568877*hamil[3]*BstarYdBmag[21]+3.872983346207417*BstarYdBmag[11]*hamil[13])*rdvpar2*rdy2)/m_; 
  alphay[23] = (0.8660254037844387*BstarYdBmag[6]*hamil[13]*rdvpar2*rdy2)/m_; 
  alphay[25] = ((1.732050807568877*(0.8784552768749174*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.5*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[25]+3.872983346207417*(0.1118033988749895*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+0.07142857142857142*b_z[4]*jacobTotInv[4]+0.1*b_z[1]*jacobTotInv[1])*hamil[8])*rdx2*rdy2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*(2.23606797749979*alphay[12]-1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.23606797749979*alphay[12]+1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])-0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*(alphay[8]+alphay[6])+0.5809475019311124*alphay[5]-0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.375*(alphay[25]+alphay[21])+0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]-0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])+0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*(alphay[8]+alphay[6])-0.5809475019311124*alphay[5]-0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])+0.375*alphay[23]-0.75*alphay[20]-0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])-0.375*alphay[23]+0.75*alphay[20]-0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])-0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.45*alphay[6]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.375*alphay[25]-0.375*alphay[21]+0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])+0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.45*alphay[6]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999998*alphay[23])-0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[6]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[6]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999998*alphay[23])+0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[6]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[6]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]-0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.45*alphay[6]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.375*alphay[25])+0.375*alphay[21]+0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]+0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.45*alphay[6]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]+0.375*alphay[23]-0.75*alphay[20]-0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]-0.375*alphay[23]+0.75*alphay[20]-0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]-0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*(alphay[8]+alphay[6])+0.5809475019311124*alphay[5]+0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.375*(alphay[25]+alphay[21]))+0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]+0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]+0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*(alphay[8]+alphay[6])-0.5809475019311124*alphay[5]+0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])-0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*(alphay[8]+alphay[6])-0.5809475019311124*alphay[5]-0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.375*(alphay[25]+alphay[21])-0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]-0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])+0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*(alphay[8]+alphay[6])+0.5809475019311124*alphay[5]-0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])+0.375*alphay[23]-0.75*alphay[20]+0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])-0.375*alphay[23]+0.75*alphay[20]+0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])-0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.45*alphay[6]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.375*alphay[25]-0.375*alphay[21]-0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])+0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.45*alphay[6]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999998*alphay[23])-0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[6]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[6]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999998*alphay[23])+0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[6]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[6]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]-0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.45*alphay[6]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.375*alphay[25])+0.375*alphay[21]-0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]+0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.45*alphay[6]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]+0.375*alphay[23]-0.75*alphay[20]+0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]-0.375*alphay[23]+0.75*alphay[20]+0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]-0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*(alphay[8]+alphay[6])-0.5809475019311124*alphay[5]+0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.375*(alphay[25]+alphay[21]))-0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]+0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]+0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*(alphay[8]+alphay[6])+0.5809475019311124*alphay[5]+0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[48]; 
  alphavpar[0] = -(0.25*(1.732050807568877*BstarYdBmag[11]*hamil[19]+1.732050807568877*(BstarYdBmag[1]*hamil[5]+BstarYdBmag[0]*hamil[2]))*rdvpar2*rdy2)/m_; 
  alphavpar[1] = ((1.732050807568877*((-0.223606797749979*hamil[5]*BstarYdBmag[11])-0.25*(BstarYdBmag[0]*hamil[5]+BstarYdBmag[1]*hamil[2]))-0.3872983346207417*BstarYdBmag[1]*hamil[19])*rdvpar2*rdy2)/m_; 
  alphavpar[2] = -(0.9682458365518543*(BstarYdBmag[1]*hamil[20]+BstarYdBmag[0]*hamil[12])*rdvpar2*rdy2)/m_; 
  alphavpar[3] = -(0.4330127018922193*(hamil[19]*BstarYdBmag[21]+hamil[5]*BstarYdBmag[6]+hamil[2]*BstarYdBmag[3])*rdvpar2*rdy2)/m_; 
  alphavpar[5] = (((-0.9682458365518543*(BstarYdBmag[0]*hamil[20]+BstarYdBmag[1]*hamil[12]))-0.8660254037844387*BstarYdBmag[11]*hamil[20])*rdvpar2*rdy2)/m_; 
  alphavpar[6] = (((-0.3872983346207417*(hamil[5]*BstarYdBmag[21]+BstarYdBmag[6]*hamil[19]))-0.4330127018922193*(hamil[2]*BstarYdBmag[6]+BstarYdBmag[3]*hamil[5]))*rdvpar2*rdy2)/m_; 
  alphavpar[7] = -(0.9682458365518543*(BstarYdBmag[6]*hamil[20]+BstarYdBmag[3]*hamil[12])*rdvpar2*rdy2)/m_; 
  alphavpar[11] = ((3.872983346207417*((-0.07142857142857142*BstarYdBmag[11])-0.1118033988749895*BstarYdBmag[0])*hamil[19]+1.732050807568877*((-0.25*hamil[2]*BstarYdBmag[11])-0.223606797749979*BstarYdBmag[1]*hamil[5]))*rdvpar2*rdy2)/m_; 
  alphavpar[15] = (((-0.8660254037844386*hamil[20]*BstarYdBmag[21])-0.9682458365518543*(BstarYdBmag[3]*hamil[20]+BstarYdBmag[6]*hamil[12]))*rdvpar2*rdy2)/m_; 
  alphavpar[19] = (((-0.8660254037844386*BstarYdBmag[1]*hamil[20])-0.9682458365518543*BstarYdBmag[11]*hamil[12])*rdvpar2*rdy2)/m_; 
  alphavpar[21] = (((-0.4330127018922193*(hamil[2]*BstarYdBmag[21]+BstarYdBmag[3]*hamil[19]))-0.276641667586244*hamil[19]*BstarYdBmag[21]-0.3872983346207417*hamil[5]*BstarYdBmag[6])*rdvpar2*rdy2)/m_; 
  alphavpar[32] = (3.872983346207417*((-0.25*hamil[12]*BstarYdBmag[21])-0.223606797749979*BstarYdBmag[6]*hamil[20])*rdvpar2*rdy2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphavpar[3]-1.0*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphavpar[3]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.6495190528383289*alphavpar[32])+0.4841229182759271*alphavpar[21]+0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.3872983346207416*alphavpar[21])+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[6]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.3872983346207416*alphavpar[21])+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[6]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.6495190528383289*alphavpar[32]+0.4841229182759271*alphavpar[21]-0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]-0.5809475019311124*alphavpar[7]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.6495190528383289*alphavpar[32])+0.4841229182759271*alphavpar[21]+0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.3872983346207416*alphavpar[21])+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[6]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.3872983346207416*alphavpar[21])+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[6]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.6495190528383289*alphavpar[32]+0.4841229182759271*alphavpar[21]-0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]-0.5809475019311124*alphavpar[7]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.6495190528383289*alphavpar[32]-0.4841229182759271*alphavpar[21]+0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.3872983346207416*alphavpar[21]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[6]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.3872983346207416*alphavpar[21]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[6]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.6495190528383289*alphavpar[32])-0.4841229182759271*alphavpar[21]-0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]+0.5809475019311124*alphavpar[7]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.6495190528383289*alphavpar[32]-0.4841229182759271*alphavpar[21]+0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.3872983346207416*alphavpar[21]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[6]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.3872983346207416*alphavpar[21]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[6]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.6495190528383289*alphavpar[32])-0.4841229182759271*alphavpar[21]-0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]+0.5809475019311124*alphavpar[7]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[19]*f[19]+alphax[11]*f[11]+alphax[5]*f[5]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[25]*f[25]+alphay[23]*f[23]+alphay[21]*f[21]+alphay[20]*f[20]+alphay[19]*f[19]+alphay[13]*f[13]+alphay[12]*f[12]+alphay[11]*f[11]+alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphavpar[32]*f[32]+alphavpar[21]*f[21]+alphavpar[19]*f[19]+alphavpar[15]*f[15]+alphavpar[11]*f[11]+alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[5] += 0.05*(7.745966692414834*(alphay[8]*f[25]+f[8]*alphay[25])+8.660254037844387*(alphay[13]*f[23]+f[13]*alphay[23])+7.745966692414834*(alphay[6]*f[21]+f[6]*alphay[21])+(8.660254037844387*alphay[12]+7.745966692414834*alphax[5])*f[20]+8.660254037844387*(f[12]*alphay[20]+alphax[11]*f[19])+7.745966692414834*(alphay[5]*f[19]+f[5]*alphay[19])+8.660254037844387*f[11]*alphax[19]+7.745966692414834*(alphax[2]*f[12]+alphay[1]*f[11]+f[1]*alphay[11])+8.660254037844386*(alphay[4]*f[8]+f[4]*alphay[8]+alphay[3]*f[6]+f[3]*alphay[6]+(alphay[2]+alphax[1])*f[5]+f[2]*alphay[5]+f[1]*alphax[5]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1])); 
  out[6] += 0.05*(8.660254037844387*alphax[19]*f[32]+7.745966692414834*(alphavpar[15]*f[32]+f[15]*alphavpar[32])+8.660254037844387*alphax[11]*f[21]+7.745966692414834*(alphavpar[6]*f[21]+f[6]*alphavpar[21]+alphavpar[5]*f[19]+f[5]*alphavpar[19])+8.660254037844386*((alphavpar[7]+alphax[5])*f[15]+f[7]*alphavpar[15])+7.745966692414834*(alphavpar[1]*f[11]+f[1]*alphavpar[11])+8.660254037844386*(alphax[2]*f[7]+(alphavpar[3]+alphax[1])*f[6]+f[3]*alphavpar[6]+alphavpar[2]*f[5]+f[2]*alphavpar[5]+alphax[0]*f[3]+alphavpar[0]*f[1]+f[0]*alphavpar[1])); 
  out[7] += 0.05*(8.660254037844387*alphay[25]*f[37]+(8.660254037844387*alphay[20]+7.745966692414834*alphavpar[15])*f[33]+8.660254037844387*((alphavpar[21]+alphay[19])*f[32]+f[21]*alphavpar[32])+7.745966692414834*(alphay[6]*f[23]+f[6]*alphay[23])+(8.660254037844387*alphay[12]+7.745966692414834*alphavpar[7])*f[22]+8.660254037844387*(alphay[11]*f[21]+f[11]*alphay[21])+7.745966692414834*alphavpar[5]*f[20]+8.660254037844387*(alphavpar[11]*f[19]+f[11]*alphavpar[19])+8.660254037844386*(alphay[8]*f[17]+(alphavpar[6]+alphay[5])*f[15]+f[6]*alphavpar[15])+7.745966692414834*(alphay[3]*f[13]+f[3]*alphay[13]+alphavpar[2]*f[12])+8.660254037844386*(alphay[4]*f[10]+(alphavpar[3]+alphay[2])*f[7]+f[3]*alphavpar[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphay[0]*f[3]+f[0]*alphay[3]+alphavpar[0]*f[2]+f[0]*alphavpar[2])); 
  out[8] += 0.05*(8.660254037844387*(alphax[19]*f[35]+alphax[11]*f[25])+8.660254037844386*(alphax[5]*f[16]+alphax[2]*f[9]+alphax[1]*f[8]+alphax[0]*f[4])); 
  out[9] += 0.05*(8.660254037844387*(alphay[23]*f[39]+alphay[21]*f[37]+alphay[20]*f[36]+alphay[19]*f[35])+7.745966692414834*alphay[8]*f[28]+8.660254037844387*(alphay[13]*f[27]+alphay[12]*f[26]+alphay[11]*f[25]+f[11]*alphay[25])+8.660254037844386*(alphay[6]*f[17]+alphay[5]*f[16])+7.745966692414834*alphay[4]*f[14]+8.660254037844386*(alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4])); 
  out[10] += 0.05*(8.660254037844387*(alphavpar[32]*f[44]+alphavpar[21]*f[37]+alphavpar[19]*f[35])+8.660254037844386*alphavpar[15]*f[31]+8.660254037844387*alphavpar[11]*f[25]+8.660254037844386*(alphavpar[7]*f[18]+alphavpar[6]*f[17]+alphavpar[5]*f[16]+alphavpar[3]*f[10]+alphavpar[2]*f[9]+alphavpar[1]*f[8]+alphavpar[0]*f[4])); 
  out[11] += 0.05*(17.32050807568877*(alphax[5]*f[19]+f[5]*alphax[19])+17.32050807568877*(alphax[1]*f[11]+f[1]*alphax[11])+19.36491673103709*(alphax[2]*f[5]+f[2]*alphax[5]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[12] += 0.05*(19.36491673103708*(alphay[25]*f[35]+alphay[23]*f[34]+alphay[21]*f[32]+alphay[13]*f[24])+17.32050807568877*(alphay[5]*f[20]+f[5]*alphay[20])+19.36491673103708*(alphay[11]*f[19]+f[11]*alphay[19])+19.36491673103709*(alphay[8]*f[16]+alphay[6]*f[15])+17.32050807568877*(alphay[2]*f[12]+f[2]*alphay[12])+19.36491673103709*(alphay[4]*f[9]+alphay[3]*f[7]+alphay[1]*f[5]+f[1]*alphay[5]+alphay[0]*f[2]+f[0]*alphay[2])); 
  out[13] += 0.05*(17.32050807568877*alphavpar[15]*f[34]+19.36491673103708*(alphavpar[19]*f[32]+f[19]*alphavpar[32])+17.32050807568877*(alphavpar[7]*f[24]+alphavpar[6]*f[23])+19.36491673103708*(alphavpar[11]*f[21]+f[11]*alphavpar[21])+19.36491673103709*(alphavpar[5]*f[15]+f[5]*alphavpar[15])+17.32050807568877*alphavpar[3]*f[13]+19.36491673103709*(alphavpar[2]*f[7]+f[2]*alphavpar[7]+alphavpar[1]*f[6]+f[1]*alphavpar[6]+alphavpar[0]*f[3]+f[0]*alphavpar[3])); 
  out[15] += 0.05*(7.745966692414834*alphay[8]*f[37]+(6.928203230275509*alphavpar[32]+8.660254037844386*alphay[12]+7.745966692414834*(alphavpar[7]+alphax[5]))*f[33]+8.660254037844386*alphax[11]*f[32]+7.745966692414834*((alphavpar[6]+alphay[5])*f[32]+f[6]*alphavpar[32]+f[17]*alphay[25])+(6.928203230275509*alphay[21]+7.745966692414834*alphay[3])*f[23]+(6.928203230275509*f[21]+7.745966692414834*f[3])*alphay[23]+(8.660254037844386*alphay[20]+7.745966692414834*(alphavpar[15]+alphax[2]))*f[22]+8.660254037844386*alphax[19]*f[21]+7.745966692414834*((alphavpar[15]+alphay[1])*f[21]+f[1]*alphay[21]+f[15]*alphavpar[21])+6.928203230275509*alphavpar[19]*f[20]+7.745966692414834*(alphavpar[2]*f[20]+alphavpar[1]*f[19]+f[15]*alphay[19]+f[1]*alphavpar[19])+8.660254037844386*(alphay[4]*f[17]+(alphavpar[3]+alphay[2]+alphax[1])*f[15]+f[3]*alphavpar[15])+7.745966692414834*(alphay[6]*f[13]+f[6]*alphay[13]+alphavpar[5]*f[12]+(alphay[6]+alphavpar[5])*f[11]+f[6]*alphay[11]+f[5]*alphavpar[11])+8.660254037844386*(alphay[8]*f[10]+(alphavpar[6]+alphay[5]+alphax[0])*f[7]+f[6]*(alphavpar[7]+alphax[5]+alphay[0])+f[0]*alphay[6]+alphavpar[0]*f[5]+f[0]*alphavpar[5]+(alphax[2]+alphay[1])*f[3]+f[1]*alphay[3]+alphavpar[1]*f[2]+f[1]*alphavpar[2])); 
  out[16] += 0.05*(8.660254037844386*alphay[13]*f[39]+7.745966692414834*alphay[6]*f[37]+(8.660254037844386*alphay[12]+7.745966692414834*alphax[5])*f[36]+(8.660254037844386*alphax[11]+7.745966692414834*alphay[5])*f[35]+(6.928203230275509*alphay[25]+7.745966692414834*alphay[4])*f[28]+8.660254037844386*alphay[23]*f[27]+(8.660254037844386*alphay[20]+7.745966692414834*alphax[2])*f[26]+8.660254037844386*alphax[19]*f[25]+7.745966692414834*(alphay[1]*f[25]+f[1]*alphay[25]+f[17]*alphay[21]+f[16]*alphay[19])+8.660254037844386*(alphay[3]*f[17]+(alphay[2]+alphax[1])*f[16])+7.745966692414834*(alphay[8]*(f[14]+f[11])+f[8]*alphay[11])+8.660254037844386*(alphay[6]*f[10]+(alphay[5]+alphax[0])*f[9]+(alphax[5]+alphay[0])*f[8]+f[0]*alphay[8]+(alphax[2]+alphay[1])*f[4]+f[1]*alphay[4])); 
  out[17] += 0.05*((8.660254037844386*alphax[19]+7.745966692414834*alphavpar[15])*f[44]+8.660254037844386*alphax[11]*f[37]+7.745966692414834*(alphavpar[6]*f[37]+alphavpar[5]*f[35])+f[31]*(7.745966692414834*alphavpar[32]+8.660254037844386*(alphavpar[7]+alphax[5]))+7.745966692414834*(alphavpar[1]*f[25]+f[17]*alphavpar[21]+f[16]*alphavpar[19])+8.660254037844386*((alphavpar[15]+alphax[2])*f[18]+(alphavpar[3]+alphax[1])*f[17]+alphavpar[2]*f[16])+7.745966692414834*f[8]*alphavpar[11]+8.660254037844386*((alphavpar[6]+alphax[0])*f[10]+alphavpar[5]*f[9]+alphavpar[0]*f[8]+alphavpar[1]*f[4])); 
  out[18] += 0.05*((8.660254037844386*alphay[20]+7.745966692414834*alphavpar[15])*f[45]+8.660254037844386*(alphavpar[21]+alphay[19])*f[44]+7.745966692414834*(alphay[8]*f[42]+alphay[6]*f[39])+(8.660254037844386*alphay[12]+7.745966692414834*alphavpar[7])*f[38]+8.660254037844386*(alphavpar[32]+alphay[11])*f[37]+7.745966692414834*alphavpar[5]*f[36]+8.660254037844386*(alphavpar[11]*f[35]+(alphavpar[6]+alphay[5])*f[31])+7.745966692414834*(alphay[4]*f[30]+alphay[3]*f[27]+alphavpar[2]*f[26])+8.660254037844386*((alphay[21]+alphavpar[19])*f[25]+f[21]*alphay[25])+7.745966692414834*f[17]*alphay[23]+8.660254037844386*((alphavpar[3]+alphay[2])*f[18]+(alphavpar[15]+alphay[1])*f[17]+alphavpar[1]*f[16])+7.745966692414834*f[10]*alphay[13]+8.660254037844386*((alphavpar[7]+alphay[0])*f[10]+alphavpar[0]*f[9]+(alphay[6]+alphavpar[5])*f[8]+f[6]*alphay[8]+(alphay[3]+alphavpar[2])*f[4]+f[3]*alphay[4])); 
  out[19] += 0.007142857142857143*(38.72983346207417*alphay[25]*f[25]+60.6217782649107*(alphay[4]*f[25]+f[4]*alphay[25])+54.22176684690384*alphay[23]*f[23]+38.72983346207417*alphay[21]*f[21]+60.6217782649107*(alphay[3]*f[21]+f[3]*alphay[21])+(54.22176684690384*alphay[20]+108.4435336938077*alphax[19]+121.2435565298214*alphax[2])*f[20]+(38.72983346207417*alphay[19]+60.6217782649107*alphay[2]+121.2435565298214*alphax[1])*f[19]+60.6217782649107*f[2]*alphay[19]+121.2435565298214*f[1]*alphax[19]+121.2435565298214*alphax[5]*f[12]+(38.72983346207417*alphay[11]+121.2435565298214*alphax[5])*f[11]+60.62177826491071*(alphay[0]*f[11]+f[0]*alphay[11])+121.2435565298214*f[5]*alphax[11]+54.22176684690384*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5])+135.5544171172596*(alphax[0]*f[5]+f[0]*alphax[5]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphay[1])); 
  out[20] += 0.05*(17.32050807568877*alphay[8]*f[35]+19.36491673103708*alphay[13]*f[34]+17.32050807568877*alphay[6]*f[32]+17.32050807568877*f[16]*alphay[25]+19.36491673103708*alphay[23]*f[24]+17.32050807568877*f[15]*alphay[21]+(15.49193338482967*alphay[19]+17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[20]+(15.49193338482967*f[19]+17.32050807568877*f[2])*alphay[20]+7.745966692414834*alphax[19]*f[19]+17.32050807568877*(alphay[1]*f[19]+f[1]*alphay[19])+19.36491673103708*(alphay[4]*f[16]+alphay[3]*f[15])+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[12]+17.32050807568877*(f[5]*alphay[12]+alphay[5]*f[11]+f[5]*alphay[11])+19.36491673103708*(alphay[8]*f[9]+alphay[6]*f[7])+7.745966692414834*alphax[5]*f[5]+19.36491673103708*(alphay[0]*f[5]+f[0]*alphay[5])+7.745966692414834*alphax[2]*f[2]+19.36491673103708*(alphay[1]*f[2]+f[1]*alphay[2])); 
  out[21] += 0.007142857142857143*((38.72983346207417*alphavpar[32]+60.62177826491071*alphavpar[7]+121.2435565298214*alphax[5])*f[32]+60.62177826491071*f[7]*alphavpar[32]+(38.72983346207417*alphavpar[21]+60.6217782649107*alphavpar[3]+121.2435565298214*alphax[1])*f[21]+60.6217782649107*f[3]*alphavpar[21]+(38.72983346207417*alphavpar[19]+60.6217782649107*alphavpar[2])*f[19]+121.2435565298214*f[15]*alphax[19]+60.6217782649107*f[2]*alphavpar[19]+(54.22176684690384*alphavpar[15]+135.5544171172596*alphax[2])*f[15]+(38.72983346207417*alphavpar[11]+60.62177826491071*alphavpar[0])*f[11]+121.2435565298214*f[6]*alphax[11]+60.62177826491071*f[0]*alphavpar[11]+135.5544171172596*alphax[5]*f[7]+(54.22176684690384*alphavpar[6]+135.5544171172596*alphax[0])*f[6]+54.22176684690384*alphavpar[5]*f[5]+135.5544171172596*alphax[1]*f[3]+54.22176684690384*alphavpar[1]*f[1]); 
  out[22] += 0.05*(19.36491673103708*alphay[25]*f[44]+17.32050807568877*alphay[6]*f[34]+(8.660254037844387*alphavpar[6]+17.32050807568877*alphay[5])*f[33]+7.745966692414834*alphavpar[32]*f[32]+19.36491673103708*(alphay[11]*f[32]+alphay[8]*f[31])+17.32050807568877*(alphay[3]*f[24]+f[15]*alphay[23])+(8.660254037844386*alphavpar[3]+17.32050807568877*alphay[2])*f[22]+19.36491673103708*(alphay[19]*f[21]+f[19]*alphay[21])+8.660254037844386*alphavpar[1]*f[20]+17.32050807568877*f[15]*alphay[20]+7.745966692414834*alphavpar[19]*f[19]+19.36491673103708*alphay[4]*f[18]+(7.745966692414834*alphavpar[15]+19.36491673103708*alphay[1])*f[15]+17.32050807568877*f[7]*alphay[13]+8.660254037844387*alphavpar[0]*f[12]+f[7]*(17.32050807568877*alphay[12]+7.745966692414834*alphavpar[7])+19.36491673103708*(alphay[0]*f[7]+alphay[5]*f[6])+f[5]*(19.36491673103708*alphay[6]+7.745966692414834*alphavpar[5])+19.36491673103708*alphay[2]*f[3]+f[2]*(19.36491673103708*alphay[3]+7.745966692414834*alphavpar[2])); 
  out[23] += 0.05*((15.49193338482967*alphavpar[32]+17.32050807568877*alphavpar[7]+8.660254037844387*alphax[5])*f[34]+17.32050807568877*(alphavpar[5]*f[32]+f[5]*alphavpar[32])+(17.32050807568877*alphavpar[15]+8.660254037844386*alphax[2])*f[24]+(15.49193338482967*alphavpar[21]+17.32050807568877*alphavpar[3]+8.660254037844386*alphax[1])*f[23]+17.32050807568877*(alphavpar[1]*f[21]+f[1]*alphavpar[21]+alphavpar[15]*f[19]+f[15]*alphavpar[19])+19.36491673103708*(alphavpar[2]*f[15]+f[2]*alphavpar[15])+(17.32050807568877*alphavpar[6]+8.660254037844387*alphax[0])*f[13]+17.32050807568877*(alphavpar[6]*f[11]+f[6]*alphavpar[11])+19.36491673103708*(alphavpar[5]*f[7]+f[5]*alphavpar[7]+alphavpar[0]*f[6]+f[0]*alphavpar[6]+alphavpar[1]*f[3]+f[1]*alphavpar[3])); 
  out[24] += 0.007142857142857143*(60.62177826491071*alphay[8]*f[39]+(121.2435565298214*alphavpar[6]+60.62177826491071*alphay[5])*f[34]+121.2435565298214*alphavpar[5]*f[33]+135.5544171172596*(alphavpar[11]*f[32]+f[11]*alphavpar[32])+60.6217782649107*alphay[4]*f[27]+(121.2435565298214*alphavpar[3]+60.6217782649107*alphay[2])*f[24]+(38.72983346207417*alphay[23]+121.2435565298214*alphavpar[15])*f[23]+60.6217782649107*(alphay[1]*f[23]+f[1]*alphay[23])+121.2435565298214*alphavpar[2]*f[22]+54.22176684690384*alphay[21]*f[21]+135.5544171172596*(alphavpar[19]*f[21]+f[19]*alphavpar[21])+121.2435565298214*alphavpar[15]*f[20]+135.5544171172596*(alphavpar[1]*f[15]+f[1]*alphavpar[15])+(38.72983346207417*alphay[13]+121.2435565298214*alphavpar[7])*f[13]+60.62177826491071*(alphay[0]*f[13]+f[0]*alphay[13])+121.2435565298214*alphavpar[7]*f[12]+135.5544171172596*(alphavpar[0]*f[7]+f[0]*alphavpar[7])+54.22176684690384*alphay[6]*f[6]+135.5544171172596*(alphavpar[5]*f[6]+f[5]*alphavpar[6])+54.22176684690384*alphay[3]*f[3]+135.5544171172596*(alphavpar[2]*f[3]+f[2]*alphavpar[3])); 
  out[25] += 0.05*(17.32050807568877*alphax[5]*f[35]+17.32050807568877*alphax[1]*f[25]+f[16]*(17.32050807568877*alphax[19]+19.36491673103708*alphax[2])+17.32050807568877*f[8]*alphax[11]+19.36491673103708*(alphax[5]*f[9]+alphax[0]*f[8]+alphax[1]*f[4])); 
  out[26] += 0.05*(19.36491673103708*(alphay[23]*f[46]+alphay[21]*f[44])+17.32050807568877*alphay[8]*f[41]+19.36491673103708*alphay[13]*f[40]+17.32050807568877*alphay[5]*f[36]+19.36491673103708*(alphay[11]*f[35]+alphay[6]*f[31])+17.32050807568877*(alphay[4]*f[29]+alphay[2]*f[26])+19.36491673103708*(alphay[19]*f[25]+f[19]*alphay[25])+17.32050807568877*f[16]*alphay[20]+19.36491673103708*(alphay[3]*f[18]+alphay[1]*f[16])+17.32050807568877*f[9]*alphay[12]+19.36491673103708*(alphay[0]*f[9]+alphay[5]*f[8]+f[5]*alphay[8]+alphay[2]*f[4]+f[2]*alphay[4])); 
  out[27] += 0.05*(17.32050807568877*alphavpar[15]*f[46]+19.36491673103708*alphavpar[19]*f[44]+17.32050807568877*(alphavpar[7]*f[40]+alphavpar[6]*f[39])+19.36491673103708*(alphavpar[11]*f[37]+alphavpar[32]*f[35]+alphavpar[5]*f[31])+17.32050807568877*alphavpar[3]*f[27]+19.36491673103708*(alphavpar[21]*f[25]+alphavpar[2]*f[18]+alphavpar[1]*f[17]+alphavpar[15]*f[16]+alphavpar[0]*f[10]+alphavpar[7]*f[9]+alphavpar[6]*f[8]+alphavpar[3]*f[4])); 
  out[28] += 0.05*(8.660254037844387*alphax[5]*f[41]+8.660254037844386*(alphax[2]*f[29]+alphax[1]*f[28])+8.660254037844387*alphax[0]*f[14]); 
  out[29] += 0.05*(8.660254037844387*(alphay[6]*f[42]+alphay[5]*f[41])+8.660254037844386*(alphay[3]*f[30]+alphay[2]*f[29]+alphay[1]*f[28])+7.745966692414834*alphay[25]*f[25]+8.660254037844387*alphay[0]*f[14]+7.745966692414834*(alphay[8]*f[8]+alphay[4]*f[4])); 
  out[30] += 0.05*(8.660254037844386*alphavpar[15]*f[47]+8.660254037844387*(alphavpar[7]*f[43]+alphavpar[6]*f[42]+alphavpar[5]*f[41])+8.660254037844386*(alphavpar[3]*f[30]+alphavpar[2]*f[29]+alphavpar[1]*f[28])+8.660254037844387*alphavpar[0]*f[14]); 
  out[31] += 0.01*((34.64101615137755*alphavpar[32]+43.30127018922195*alphay[12]+38.72983346207417*(alphavpar[7]+alphax[5]))*f[45]+(43.30127018922195*alphax[11]+38.72983346207417*(alphavpar[6]+alphay[5]))*f[44]+(34.64101615137755*alphay[25]+38.72983346207418*alphay[4])*f[42]+(34.64101615137755*alphay[21]+38.72983346207418*alphay[3])*f[39]+(43.30127018922195*alphay[20]+38.72983346207418*(alphavpar[15]+alphax[2]))*f[38]+(34.64101615137755*alphay[23]+43.30127018922195*alphax[19]+38.72983346207418*(alphavpar[15]+alphay[1]))*f[37]+34.64101615137755*alphavpar[19]*f[36]+38.72983346207418*(alphavpar[2]*f[36]+alphavpar[1]*f[35]+f[17]*alphavpar[32])+(38.72983346207417*(alphavpar[21]+alphay[19])+43.30127018922193*(alphavpar[3]+alphay[2]+alphax[1]))*f[31]+38.72983346207417*(alphay[8]*f[30]+alphay[6]*f[27]+alphavpar[5]*f[26]+(alphay[6]+alphavpar[5])*f[25]+f[6]*alphay[25]+f[10]*alphay[23]+alphay[8]*f[21]+f[8]*(alphay[21]+alphavpar[19]))+43.30127018922193*(alphavpar[6]+alphay[5]+alphax[0])*f[18]+(38.72983346207418*(alphay[13]+alphay[11])+43.30127018922193*(alphavpar[7]+alphax[5]+alphay[0]))*f[17]+38.72983346207418*alphavpar[11]*f[16]+43.30127018922193*(alphavpar[0]*f[16]+f[10]*(alphavpar[15]+alphax[2]+alphay[1])+alphavpar[1]*f[9]+(alphay[3]+alphavpar[2])*f[8]+f[3]*alphay[8]+alphay[4]*f[6]+f[4]*(alphay[6]+alphavpar[5]))); 
  out[32] += 0.001428571428571429*((193.6491673103708*alphay[25]+303.1088913245535*alphay[4])*f[37]+(271.1088342345192*alphay[20]+542.2176684690384*alphax[19]+242.4871130596428*alphavpar[15]+606.217782649107*alphax[2])*f[33]+(193.6491673103708*(alphavpar[21]+alphay[19])+303.1088913245535*(alphavpar[3]+alphay[2])+606.217782649107*alphax[1])*f[32]+(271.1088342345192*f[22]+193.6491673103708*f[21]+303.1088913245535*f[3])*alphavpar[32]+303.1088913245536*f[10]*alphay[25]+242.4871130596428*(alphay[6]*f[23]+f[6]*alphay[23])+606.2177826491072*alphax[5]*f[22]+(271.1088342345192*alphay[13]+193.6491673103708*alphay[11]+303.1088913245536*alphavpar[7]+606.2177826491072*alphax[5]+303.1088913245536*alphay[0])*f[21]+(271.1088342345192*f[13]+193.6491673103708*f[11])*alphay[21]+303.1088913245536*(f[0]*alphay[21]+f[7]*alphavpar[21])+242.4871130596428*alphavpar[5]*f[20]+193.6491673103708*alphavpar[11]*f[19]+303.1088913245536*(alphavpar[0]*f[19]+f[7]*alphay[19])+606.2177826491072*f[6]*alphax[19]+(271.1088342345192*f[12]+193.6491673103708*f[11]+303.1088913245536*f[0])*alphavpar[19]+271.1088342345192*alphay[8]*f[17]+(606.217782649107*alphax[11]+271.1088342345192*(alphavpar[6]+alphay[5])+677.7720855862981*alphax[0])*f[15]+271.1088342345192*f[6]*alphavpar[15]+303.1088913245535*((alphay[3]+alphavpar[2])*f[11]+f[3]*alphay[11]+f[2]*alphavpar[11])+677.7720855862981*(alphax[1]*f[7]+alphax[2]*f[6])+271.1088342345192*(alphay[1]*f[6]+f[1]*alphay[6]+alphavpar[1]*f[5])+677.7720855862981*f[3]*alphax[5]+271.1088342345192*f[1]*alphavpar[5]); 
  out[33] += 0.01*(86.60254037844389*alphay[8]*f[44]+(77.45966692414834*alphay[21]+86.60254037844386*alphay[3])*f[34]+(38.72983346207417*alphavpar[21]+77.45966692414834*alphay[19]+43.30127018922193*alphavpar[3]+86.60254037844386*alphay[2]+43.30127018922193*alphax[1])*f[33]+(77.45966692414834*(alphay[23]+alphay[20])+38.72983346207417*alphax[19]+34.64101615137754*alphavpar[15]+86.60254037844386*alphay[1])*f[32]+34.64101615137754*f[15]*alphavpar[32]+(86.60254037844389*alphay[25]+96.82458365518544*alphay[4])*f[31]+86.60254037844389*(alphay[6]*f[24]+f[7]*alphay[23])+(43.30127018922195*alphavpar[6]+86.60254037844389*alphay[5]+43.30127018922195*alphax[0])*f[22]+86.60254037844389*(alphay[5]*f[21]+f[5]*alphay[21])+(38.72983346207417*alphavpar[11]+43.30127018922195*alphavpar[0])*f[20]+86.60254037844389*f[7]*alphay[20]+(86.60254037844389*alphay[6]+34.64101615137755*alphavpar[5])*f[19]+86.60254037844389*f[6]*alphay[19]+34.64101615137755*f[5]*alphavpar[19]+96.82458365518544*alphay[8]*f[18]+(86.60254037844386*(alphay[13]+alphay[12]+alphay[11])+38.72983346207418*(alphavpar[7]+alphax[5])+96.82458365518544*alphay[0])*f[15]+38.72983346207418*f[7]*alphavpar[15]+43.30127018922193*alphavpar[1]*f[12]+38.72983346207418*alphax[2]*f[7]+96.82458365518544*(alphay[1]*f[7]+alphay[2]*f[6]+f[2]*alphay[6])+(96.82458365518544*alphay[3]+38.72983346207418*alphavpar[2])*f[5]+96.82458365518544*f[3]*alphay[5]+38.72983346207418*f[2]*alphavpar[5]); 
  out[34] += 0.001428571428571429*((271.1088342345192*alphay[25]+303.1088913245535*alphay[4])*f[39]+(542.2176684690384*alphavpar[21]+271.1088342345192*alphay[19]+606.217782649107*alphavpar[3]+303.1088913245535*(alphay[2]+alphax[1]))*f[34]+542.2176684690384*alphavpar[19]*f[33]+606.217782649107*(alphavpar[2]*f[33]+alphavpar[1]*f[32])+(542.2176684690384*(f[23]+f[20])+606.217782649107*f[1])*alphavpar[32]+303.1088913245536*alphay[8]*f[27]+(606.2177826491072*alphavpar[6]+303.1088913245536*(alphay[5]+alphax[0]))*f[24]+(193.6491673103708*alphay[13]+271.1088342345192*alphay[11]+606.2177826491072*alphavpar[7]+303.1088913245536*(alphax[5]+alphay[0]))*f[23]+(193.6491673103708*f[13]+271.1088342345192*f[11]+303.1088913245536*f[0])*alphay[23]+606.2177826491072*alphavpar[5]*f[22]+(242.4871130596428*alphay[6]+606.2177826491072*alphavpar[5])*f[21]+242.4871130596428*f[6]*alphay[21]+606.2177826491072*(f[5]*alphavpar[21]+alphavpar[7]*f[20]+alphavpar[6]*f[19]+f[6]*alphavpar[19])+(606.217782649107*alphavpar[11]+677.7720855862981*alphavpar[0])*f[15]+(606.217782649107*(f[13]+f[12]+f[11])+677.7720855862981*f[0])*alphavpar[15]+303.1088913245535*((alphax[2]+alphay[1])*f[13]+f[1]*alphay[13])+677.7720855862981*(alphavpar[1]*f[7]+f[1]*alphavpar[7])+(271.1088342345192*alphay[3]+677.7720855862981*alphavpar[2])*f[6]+271.1088342345192*f[3]*alphay[6]+677.7720855862981*(f[2]*alphavpar[6]+alphavpar[3]*f[5]+f[3]*alphavpar[5])); 
  out[35] += 0.001428571428571429*(271.1088342345192*alphay[23]*f[39]+(193.6491673103708*alphay[21]+303.1088913245535*alphay[3])*f[37]+(271.1088342345192*alphay[20]+542.2176684690384*alphax[19]+606.217782649107*alphax[2])*f[36]+(193.6491673103708*alphay[19]+303.1088913245535*alphay[2]+606.217782649107*alphax[1])*f[35]+242.4871130596428*alphay[8]*f[28]+606.2177826491072*alphax[5]*f[26]+(193.6491673103708*alphay[11]+606.2177826491072*alphax[5]+303.1088913245536*alphay[0])*f[25]+(271.1088342345192*f[14]+193.6491673103708*f[11])*alphay[25]+303.1088913245536*(f[0]*alphay[25]+f[10]*alphay[21]+f[9]*alphay[19])+606.2177826491072*f[8]*alphax[19]+271.1088342345192*alphay[6]*f[17]+(606.217782649107*alphax[11]+271.1088342345192*alphay[5]+677.7720855862981*alphax[0])*f[16]+303.1088913245535*(alphay[4]*f[11]+f[4]*alphay[11])+677.7720855862981*(alphax[1]*f[9]+alphax[2]*f[8])+271.1088342345192*(alphay[1]*f[8]+f[1]*alphay[8])+677.7720855862981*f[4]*alphax[5]); 
  out[36] += 0.05*(19.36491673103708*alphay[13]*f[46]+17.32050807568877*alphay[6]*f[44]+(15.49193338482967*alphay[25]+17.32050807568877*alphay[4])*f[41]+19.36491673103708*alphay[23]*f[40]+(15.49193338482967*alphay[19]+17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[36]+(15.49193338482967*alphay[20]+7.745966692414834*alphax[19]+17.32050807568877*alphay[1])*f[35]+(17.32050807568877*alphay[21]+19.36491673103709*alphay[3])*f[31]+17.32050807568877*alphay[8]*f[29]+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[26]+17.32050807568877*(alphay[5]*f[25]+f[5]*alphay[25]+f[9]*alphay[20]+alphay[8]*f[19]+f[8]*alphay[19])+19.36491673103709*alphay[6]*f[18]+(17.32050807568877*(alphay[12]+alphay[11])+7.745966692414834*alphax[5]+19.36491673103709*alphay[0])*f[16]+7.745966692414834*alphax[2]*f[9]+19.36491673103709*(alphay[1]*f[9]+alphay[2]*f[8]+f[2]*alphay[8]+alphay[4]*f[5]+f[4]*alphay[5])); 
  out[37] += 0.007142857142857143*((38.72983346207417*alphavpar[32]+60.62177826491071*alphavpar[7]+121.2435565298214*alphax[5])*f[44]+(38.72983346207417*alphavpar[21]+60.6217782649107*alphavpar[3]+121.2435565298214*alphax[1])*f[37]+38.72983346207417*alphavpar[19]*f[35]+60.6217782649107*(alphavpar[2]*f[35]+f[18]*alphavpar[32])+(121.2435565298214*alphax[19]+54.22176684690384*alphavpar[15]+135.5544171172596*alphax[2])*f[31]+38.72983346207417*alphavpar[11]*f[25]+60.62177826491071*(alphavpar[0]*f[25]+f[10]*alphavpar[21]+f[9]*alphavpar[19])+135.5544171172596*alphax[5]*f[18]+(121.2435565298214*alphax[11]+54.22176684690384*alphavpar[6]+135.5544171172596*alphax[0])*f[17]+54.22176684690384*alphavpar[5]*f[16]+60.6217782649107*f[4]*alphavpar[11]+135.5544171172596*alphax[1]*f[10]+54.22176684690384*alphavpar[1]*f[8]); 
  out[38] += 0.05*(17.32050807568877*(alphay[8]*f[47]+alphay[6]*f[46])+(8.660254037844387*alphavpar[6]+17.32050807568877*alphay[5])*f[45]+(7.745966692414834*alphavpar[32]+19.36491673103708*alphay[11])*f[44]+17.32050807568877*(alphay[4]*f[43]+alphay[3]*f[40])+(8.660254037844386*alphavpar[3]+17.32050807568877*alphay[2])*f[38]+19.36491673103708*alphay[19]*f[37]+8.660254037844386*alphavpar[1]*f[36]+(19.36491673103708*alphay[21]+7.745966692414834*alphavpar[19])*f[35]+19.36491673103708*alphay[25]*f[32]+(17.32050807568877*(alphay[23]+alphay[20])+7.745966692414834*alphavpar[15]+19.36491673103709*alphay[1])*f[31]+8.660254037844387*alphavpar[0]*f[26]+(17.32050807568877*(alphay[13]+alphay[12])+7.745966692414834*alphavpar[7])*f[18]+19.36491673103709*(alphay[0]*f[18]+alphay[5]*f[17])+(19.36491673103709*alphay[6]+7.745966692414834*alphavpar[5])*f[16]+19.36491673103709*(alphay[8]*f[15]+alphay[2]*f[10])+(19.36491673103709*alphay[3]+7.745966692414834*alphavpar[2])*f[9]+19.36491673103709*alphay[4]*f[7]); 
  out[39] += 0.05*((15.49193338482967*alphavpar[32]+17.32050807568877*alphavpar[7]+8.660254037844387*alphax[5])*f[46]+17.32050807568877*alphavpar[5]*f[44]+(17.32050807568877*alphavpar[15]+8.660254037844386*alphax[2])*f[40]+(15.49193338482967*alphavpar[21]+17.32050807568877*alphavpar[3]+8.660254037844386*alphax[1])*f[39]+17.32050807568877*(alphavpar[1]*f[37]+alphavpar[15]*f[35]+f[16]*alphavpar[32])+(17.32050807568877*alphavpar[19]+19.36491673103709*alphavpar[2])*f[31]+(17.32050807568877*alphavpar[6]+8.660254037844387*alphax[0])*f[27]+17.32050807568877*(alphavpar[6]*f[25]+f[8]*alphavpar[21])+19.36491673103709*alphavpar[5]*f[18]+17.32050807568877*alphavpar[11]*f[17]+19.36491673103709*(alphavpar[0]*f[17]+alphavpar[7]*f[16]+f[9]*alphavpar[15]+alphavpar[1]*f[10]+alphavpar[3]*f[8]+f[4]*alphavpar[6])); 
  out[40] += 0.007142857142857143*((121.2435565298214*alphavpar[6]+60.62177826491071*alphay[5])*f[46]+121.2435565298214*alphavpar[5]*f[45]+135.5544171172596*alphavpar[11]*f[44]+(121.2435565298214*alphavpar[3]+60.6217782649107*alphay[2])*f[40]+(38.72983346207417*alphay[23]+121.2435565298214*alphavpar[15]+60.6217782649107*alphay[1])*f[39]+121.2435565298214*alphavpar[2]*f[38]+(54.22176684690384*alphay[21]+135.5544171172596*alphavpar[19])*f[37]+121.2435565298214*alphavpar[15]*f[36]+135.5544171172596*(alphavpar[21]*f[35]+f[25]*alphavpar[32])+135.5544171172596*alphavpar[1]*f[31]+(38.72983346207417*alphay[13]+121.2435565298214*alphavpar[7]+60.62177826491071*alphay[0])*f[27]+121.2435565298214*alphavpar[7]*f[26]+60.62177826491071*(alphay[8]*f[23]+f[8]*alphay[23])+135.5544171172596*alphavpar[0]*f[18]+54.22176684690384*alphay[6]*f[17]+135.5544171172596*(alphavpar[5]*f[17]+alphavpar[6]*f[16]+f[8]*alphavpar[15])+60.6217782649107*(alphay[4]*f[13]+f[4]*alphay[13])+54.22176684690384*alphay[3]*f[10]+135.5544171172596*(alphavpar[2]*f[10]+alphavpar[3]*f[9]+f[4]*alphavpar[7])); 
  out[41] += 0.01*((38.72983346207417*alphay[21]+43.30127018922193*alphay[3])*f[42]+(38.72983346207417*alphay[19]+43.30127018922193*(alphay[2]+alphax[1]))*f[41]+43.30127018922195*(alphay[6]*f[30]+(alphay[5]+alphax[0])*f[29])+(38.72983346207417*alphay[11]+43.30127018922195*(alphax[5]+alphay[0]))*f[28]+34.64101615137755*(alphay[8]*f[25]+f[8]*alphay[25])+43.30127018922193*(alphax[2]+alphay[1])*f[14]+38.72983346207418*(alphay[4]*f[8]+f[4]*alphay[8])); 
  out[42] += 0.05*((7.745966692414834*alphavpar[32]+8.660254037844387*(alphavpar[7]+alphax[5]))*f[47]+8.660254037844386*(alphavpar[15]+alphax[2])*f[43]+(7.745966692414834*alphavpar[21]+8.660254037844386*(alphavpar[3]+alphax[1]))*f[42]+(7.745966692414834*alphavpar[19]+8.660254037844386*alphavpar[2])*f[41]+8.660254037844387*((alphavpar[6]+alphax[0])*f[30]+alphavpar[5]*f[29])+(7.745966692414834*alphavpar[11]+8.660254037844387*alphavpar[0])*f[28]+8.660254037844386*alphavpar[1]*f[14]); 
  out[43] += 0.05*(8.660254037844387*(alphavpar[6]+alphay[5])*f[47]+8.660254037844386*(alphavpar[3]+alphay[2])*f[43]+7.745966692414834*alphay[23]*f[42]+8.660254037844386*((alphavpar[15]+alphay[1])*f[42]+alphavpar[1]*f[41])+7.745966692414834*(alphay[25]*f[37]+alphay[13]*f[30])+8.660254037844387*((alphavpar[7]+alphay[0])*f[30]+alphavpar[0]*f[29]+(alphay[6]+alphavpar[5])*f[28])+7.745966692414834*alphay[8]*f[17]+8.660254037844386*(alphay[3]+alphavpar[2])*f[14]+7.745966692414834*alphay[4]*f[10]); 
  out[44] += 0.001428571428571429*((271.1088342345192*alphay[20]+542.2176684690384*alphax[19]+242.4871130596428*alphavpar[15]+606.217782649107*alphax[2])*f[45]+(193.6491673103708*(alphavpar[21]+alphay[19])+303.1088913245535*(alphavpar[3]+alphay[2])+606.217782649107*alphax[1])*f[44]+242.4871130596428*(alphay[8]*f[42]+alphay[6]*f[39])+(271.1088342345192*alphavpar[32]+606.2177826491072*alphax[5])*f[38]+(193.6491673103708*alphavpar[32]+271.1088342345192*alphay[13]+193.6491673103708*alphay[11]+303.1088913245536*alphavpar[7]+606.2177826491072*alphax[5]+303.1088913245536*alphay[0])*f[37]+242.4871130596428*alphavpar[5]*f[36]+193.6491673103708*alphavpar[11]*f[35]+303.1088913245536*(alphavpar[0]*f[35]+f[10]*alphavpar[32])+(606.2177826491072*alphax[11]+271.1088342345192*(alphavpar[6]+alphay[5])+677.772085586298*alphax[0])*f[31]+271.1088342345192*(alphay[25]*f[30]+alphay[21]*f[27]+alphavpar[19]*f[26])+(193.6491673103708*(alphay[21]+alphavpar[19])+303.1088913245535*(alphay[3]+alphavpar[2]))*f[25]+(193.6491673103708*f[21]+303.1088913245535*f[3])*alphay[25]+242.4871130596428*f[17]*alphay[23]+303.1088913245535*(alphay[4]*f[21]+f[4]*alphay[21]+f[18]*(alphavpar[21]+alphay[19]))+606.217782649107*f[17]*alphax[19]+303.1088913245535*f[4]*alphavpar[19]+677.772085586298*alphax[1]*f[18]+(271.1088342345192*alphavpar[15]+677.772085586298*alphax[2])*f[17]+271.1088342345192*(alphay[1]*f[17]+alphavpar[1]*f[16])+303.1088913245536*(f[10]*alphay[11]+f[9]*alphavpar[11])+677.772085586298*alphax[5]*f[10]+271.1088342345192*((alphay[6]+alphavpar[5])*f[8]+f[6]*alphay[8])); 
  out[45] += 0.01*((77.45966692414834*alphay[25]+86.60254037844386*alphay[4])*f[47]+(77.45966692414834*alphay[21]+86.60254037844386*alphay[3])*f[46]+(38.72983346207417*alphavpar[21]+77.45966692414834*alphay[19]+43.30127018922193*alphavpar[3]+86.60254037844386*alphay[2]+43.30127018922193*alphax[1])*f[45]+(77.45966692414834*(alphay[23]+alphay[20])+38.72983346207417*alphax[19]+34.64101615137754*alphavpar[15]+86.60254037844386*alphay[1])*f[44]+86.60254037844389*(alphay[8]*f[43]+alphay[6]*f[40])+(43.30127018922195*alphavpar[6]+86.60254037844389*alphay[5]+43.30127018922195*alphax[0])*f[38]+86.60254037844389*alphay[5]*f[37]+(38.72983346207417*alphavpar[11]+43.30127018922195*alphavpar[0])*f[36]+(86.60254037844389*alphay[6]+34.64101615137755*alphavpar[5])*f[35]+86.60254037844389*alphay[8]*f[32]+f[31]*(34.64101615137755*alphavpar[32]+86.60254037844389*(alphay[13]+alphay[12]+alphay[11])+38.72983346207417*(alphavpar[7]+alphax[5])+96.82458365518542*alphay[0])+43.30127018922193*alphavpar[1]*f[26]+86.60254037844386*(f[15]*alphay[25]+f[18]*alphay[23]+f[16]*alphay[21]+f[18]*alphay[20]+f[17]*alphay[19])+34.64101615137754*f[16]*alphavpar[19]+38.72983346207417*(alphavpar[15]+alphax[2])*f[18]+96.82458365518542*(alphay[1]*f[18]+alphay[2]*f[17])+(96.82458365518542*alphay[3]+38.72983346207417*alphavpar[2])*f[16]+96.82458365518542*(alphay[4]*f[15]+alphay[5]*f[10])+(96.82458365518542*alphay[6]+38.72983346207417*alphavpar[5])*f[9]+96.82458365518542*f[7]*alphay[8]); 
  out[46] += 0.001428571428571429*((542.2176684690384*alphavpar[21]+271.1088342345192*alphay[19]+606.217782649107*alphavpar[3]+303.1088913245535*(alphay[2]+alphax[1]))*f[46]+542.2176684690384*alphavpar[19]*f[45]+606.217782649107*(alphavpar[2]*f[45]+alphavpar[1]*f[44])+(606.2177826491072*alphavpar[6]+303.1088913245536*(alphay[5]+alphax[0]))*f[40]+(542.2176684690384*alphavpar[32]+193.6491673103708*alphay[13]+271.1088342345192*alphay[11]+606.2177826491072*alphavpar[7]+303.1088913245536*(alphax[5]+alphay[0]))*f[39]+606.2177826491072*alphavpar[5]*f[38]+(242.4871130596428*alphay[6]+606.2177826491072*alphavpar[5])*f[37]+542.2176684690384*alphavpar[32]*f[36]+606.2177826491072*(alphavpar[7]*f[36]+alphavpar[6]*f[35]+f[8]*alphavpar[32])+(606.2177826491072*alphavpar[11]+677.772085586298*alphavpar[0])*f[31]+(193.6491673103708*alphay[23]+606.217782649107*alphavpar[15]+303.1088913245535*(alphax[2]+alphay[1]))*f[27]+606.217782649107*alphavpar[15]*f[26]+(271.1088342345192*alphay[23]+606.217782649107*alphavpar[15])*f[25]+271.1088342345192*f[23]*alphay[25]+303.1088913245535*(alphay[4]*f[23]+f[4]*alphay[23])+242.4871130596428*f[17]*alphay[21]+606.217782649107*(f[16]*alphavpar[21]+f[17]*alphavpar[19])+677.772085586298*alphavpar[1]*f[18]+271.1088342345192*alphay[3]*f[17]+677.772085586298*(alphavpar[2]*f[17]+alphavpar[3]*f[16]+f[4]*alphavpar[15])+303.1088913245536*(alphay[8]*f[13]+f[8]*alphay[13])+271.1088342345192*alphay[6]*f[10]+677.772085586298*(alphavpar[5]*f[10]+alphavpar[6]*f[9]+alphavpar[7]*f[8])); 
  out[47] += 0.01*((38.72983346207417*(alphavpar[21]+alphay[19])+43.30127018922193*(alphavpar[3]+alphay[2]+alphax[1]))*f[47]+43.30127018922195*(alphavpar[6]+alphay[5]+alphax[0])*f[43]+(38.72983346207417*(alphavpar[32]+alphay[13]+alphay[11])+43.30127018922195*(alphavpar[7]+alphax[5]+alphay[0]))*f[42]+(38.72983346207417*alphavpar[11]+43.30127018922195*alphavpar[0])*f[41]+34.64101615137755*alphay[8]*f[37]+38.72983346207417*alphay[23]*f[30]+43.30127018922193*((alphavpar[15]+alphax[2]+alphay[1])*f[30]+alphavpar[1]*f[29])+(38.72983346207417*(alphay[21]+alphavpar[19])+43.30127018922193*(alphay[3]+alphavpar[2]))*f[28]+f[17]*(34.64101615137754*alphay[25]+38.72983346207417*alphay[4])+43.30127018922195*(alphay[6]+alphavpar[5])*f[14]+38.72983346207417*alphay[8]*f[10]); 
  return cflFreq; 
} 
double GyrokineticGenGeoVol2x2vSerP2_Bvarsz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f, double *out) 
{ 
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

  double hamil[48]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[11] = 2.0*phi[4]*q_; 
  hamil[12] = 2.0*phi[5]*q_; 
  hamil[13] = (0.5962847939999438*m_)/rdvpar2Sq; 
  hamil[19] = 2.0*phi[6]*q_; 
  hamil[20] = 2.0*phi[7]*q_; 

  double BstarXdBmag[48]; 

  double BstarYdBmag[48]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[48]; 
  alphax[0] = -(0.4330127018922193*b_z[0]*jacobTotInv[0]*hamil[2]*rdx2*rdy2)/q_; 
  alphax[1] = -(0.4330127018922193*b_z[0]*jacobTotInv[0]*hamil[5]*rdx2*rdy2)/q_; 
  alphax[2] = -(0.9682458365518543*b_z[0]*jacobTotInv[0]*hamil[12]*rdx2*rdy2)/q_; 
  alphax[5] = -(0.9682458365518543*b_z[0]*jacobTotInv[0]*hamil[20]*rdx2*rdy2)/q_; 
  alphax[11] = -(0.4330127018922193*b_z[0]*jacobTotInv[0]*hamil[19]*rdx2*rdy2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*(2.23606797749979*alphax[11]-1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.23606797749979*alphax[11]+1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[48]; 
  alphay[0] = (0.4330127018922193*b_z[0]*jacobTotInv[0]*hamil[1]*rdx2*rdy2)/q_; 
  alphay[1] = (0.9682458365518543*b_z[0]*jacobTotInv[0]*hamil[11]*rdx2*rdy2)/q_; 
  alphay[2] = (0.4330127018922193*b_z[0]*jacobTotInv[0]*hamil[5]*rdx2*rdy2)/q_; 
  alphay[5] = (0.9682458365518543*b_z[0]*jacobTotInv[0]*hamil[19]*rdx2*rdy2)/q_; 
  alphay[12] = (0.4330127018922194*b_z[0]*jacobTotInv[0]*hamil[20]*rdx2*rdy2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*(2.23606797749979*alphay[12]-1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.23606797749979*alphay[12]+1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]-0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphay[12]+0.5809475019311124*alphay[5]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[48]; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.0; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
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
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[11]*f[11]+alphax[5]*f[5]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[12]*f[12]+alphay[5]*f[5]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[5] += 0.05*((8.660254037844387*alphay[12]+7.745966692414834*alphax[5])*f[20]+8.660254037844387*alphax[11]*f[19]+7.745966692414834*(alphay[5]*f[19]+alphax[2]*f[12]+alphay[1]*f[11])+8.660254037844386*((alphay[2]+alphax[1])*f[5]+f[2]*alphay[5]+f[1]*alphax[5]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1])); 
  out[6] += 0.05*(8.660254037844387*alphax[11]*f[21]+8.660254037844386*(alphax[5]*f[15]+alphax[2]*f[7]+alphax[1]*f[6]+alphax[0]*f[3])); 
  out[7] += 0.05*(8.660254037844387*alphay[12]*f[22]+8.660254037844386*(alphay[5]*f[15]+alphay[2]*f[7]+alphay[1]*f[6]+alphay[0]*f[3])); 
  out[8] += 0.05*(8.660254037844387*alphax[11]*f[25]+8.660254037844386*(alphax[5]*f[16]+alphax[2]*f[9]+alphax[1]*f[8]+alphax[0]*f[4])); 
  out[9] += 0.05*(8.660254037844387*alphay[12]*f[26]+8.660254037844386*(alphay[5]*f[16]+alphay[2]*f[9]+alphay[1]*f[8]+alphay[0]*f[4])); 
  out[11] += 0.05*(17.32050807568877*alphax[5]*f[19]+17.32050807568877*(alphax[1]*f[11]+f[1]*alphax[11])+19.36491673103709*(alphax[2]*f[5]+f[2]*alphax[5]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[12] += 0.05*(17.32050807568877*alphay[5]*f[20]+17.32050807568877*(alphay[2]*f[12]+f[2]*alphay[12])+19.36491673103709*(alphay[1]*f[5]+f[1]*alphay[5]+alphay[0]*f[2]+f[0]*alphay[2])); 
  out[15] += 0.05*((8.660254037844386*alphay[12]+7.745966692414834*alphax[5])*f[33]+8.660254037844386*alphax[11]*f[32]+7.745966692414834*(alphay[5]*f[32]+alphax[2]*f[22]+alphay[1]*f[21])+8.660254037844386*((alphay[2]+alphax[1])*f[15]+(alphay[5]+alphax[0])*f[7]+(alphax[5]+alphay[0])*f[6]+(alphax[2]+alphay[1])*f[3])); 
  out[16] += 0.05*((8.660254037844386*alphay[12]+7.745966692414834*alphax[5])*f[36]+8.660254037844386*alphax[11]*f[35]+7.745966692414834*(alphay[5]*f[35]+alphax[2]*f[26]+alphay[1]*f[25])+8.660254037844386*((alphay[2]+alphax[1])*f[16]+(alphay[5]+alphax[0])*f[9]+(alphax[5]+alphay[0])*f[8]+(alphax[2]+alphay[1])*f[4])); 
  out[17] += 0.4330127018922193*(alphax[11]*f[37]+alphax[5]*f[31]+alphax[2]*f[18]+alphax[1]*f[17]+alphax[0]*f[10]); 
  out[18] += 0.4330127018922193*(alphay[12]*f[38]+alphay[5]*f[31]+alphay[2]*f[18]+alphay[1]*f[17]+alphay[0]*f[10]); 
  out[19] += 0.05*(17.32050807568877*alphax[2]*f[20]+(8.660254037844386*alphay[2]+17.32050807568877*alphax[1])*f[19]+17.32050807568877*alphax[5]*f[12]+(17.32050807568877*alphax[5]+8.660254037844387*alphay[0])*f[11]+f[5]*(17.32050807568877*alphax[11]+7.745966692414834*alphay[5])+19.36491673103708*(alphax[0]*f[5]+f[0]*alphax[5]+alphax[1]*f[2])+f[1]*(19.36491673103708*alphax[2]+7.745966692414834*alphay[1])); 
  out[20] += 0.05*((17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[20]+17.32050807568877*alphay[1]*f[19]+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[12]+17.32050807568877*(f[5]*alphay[12]+alphay[5]*f[11])+7.745966692414834*alphax[5]*f[5]+19.36491673103708*(alphay[0]*f[5]+f[0]*alphay[5])+7.745966692414834*alphax[2]*f[2]+19.36491673103708*(alphay[1]*f[2]+f[1]*alphay[2])); 
  out[21] += 0.05*(17.32050807568877*alphax[5]*f[32]+17.32050807568877*alphax[1]*f[21]+19.36491673103708*alphax[2]*f[15]+17.32050807568877*f[6]*alphax[11]+19.36491673103708*(alphax[5]*f[7]+alphax[0]*f[6]+alphax[1]*f[3])); 
  out[22] += 0.05*(17.32050807568877*alphay[5]*f[33]+17.32050807568877*alphay[2]*f[22]+19.36491673103708*alphay[1]*f[15]+17.32050807568877*f[7]*alphay[12]+19.36491673103708*(alphay[0]*f[7]+alphay[5]*f[6]+alphay[2]*f[3])); 
  out[23] += 0.05*(8.660254037844387*alphax[5]*f[34]+8.660254037844386*(alphax[2]*f[24]+alphax[1]*f[23])+8.660254037844387*alphax[0]*f[13]); 
  out[24] += 0.05*(8.660254037844387*alphay[5]*f[34]+8.660254037844386*(alphay[2]*f[24]+alphay[1]*f[23])+8.660254037844387*alphay[0]*f[13]); 
  out[25] += 0.05*(17.32050807568877*alphax[5]*f[35]+17.32050807568877*alphax[1]*f[25]+19.36491673103708*alphax[2]*f[16]+17.32050807568877*f[8]*alphax[11]+19.36491673103708*(alphax[5]*f[9]+alphax[0]*f[8]+alphax[1]*f[4])); 
  out[26] += 0.05*(17.32050807568877*alphay[5]*f[36]+17.32050807568877*alphay[2]*f[26]+19.36491673103708*alphay[1]*f[16]+17.32050807568877*f[9]*alphay[12]+19.36491673103708*(alphay[0]*f[9]+alphay[5]*f[8]+alphay[2]*f[4])); 
  out[28] += 0.05*(8.660254037844387*alphax[5]*f[41]+8.660254037844386*(alphax[2]*f[29]+alphax[1]*f[28])+8.660254037844387*alphax[0]*f[14]); 
  out[29] += 0.05*(8.660254037844387*alphay[5]*f[41]+8.660254037844386*(alphay[2]*f[29]+alphay[1]*f[28])+8.660254037844387*alphay[0]*f[14]); 
  out[31] += 0.05*((8.660254037844387*alphay[12]+7.745966692414834*alphax[5])*f[45]+8.660254037844387*alphax[11]*f[44]+7.745966692414834*(alphay[5]*f[44]+alphax[2]*f[38]+alphay[1]*f[37])+8.660254037844386*((alphay[2]+alphax[1])*f[31]+(alphay[5]+alphax[0])*f[18]+(alphax[5]+alphay[0])*f[17]+(alphax[2]+alphay[1])*f[10])); 
  out[32] += 0.05*(17.32050807568877*alphax[2]*f[33]+(8.660254037844386*alphay[2]+17.32050807568877*alphax[1])*f[32]+17.32050807568877*alphax[5]*f[22]+(17.32050807568877*alphax[5]+8.660254037844387*alphay[0])*f[21]+(17.32050807568877*alphax[11]+7.745966692414834*alphay[5])*f[15]+19.36491673103709*(alphax[0]*f[15]+alphax[1]*f[7])+(19.36491673103709*alphax[2]+7.745966692414834*alphay[1])*f[6]+19.36491673103709*f[3]*alphax[5]); 
  out[33] += 0.05*((17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[33]+17.32050807568877*alphay[1]*f[32]+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[22]+17.32050807568877*alphay[5]*f[21]+(17.32050807568877*alphay[12]+7.745966692414834*alphax[5]+19.36491673103709*alphay[0])*f[15]+7.745966692414834*alphax[2]*f[7]+19.36491673103709*(alphay[1]*f[7]+alphay[2]*f[6]+f[3]*alphay[5])); 
  out[34] += 0.05*(8.660254037844386*(alphay[2]+alphax[1])*f[34]+8.660254037844387*((alphay[5]+alphax[0])*f[24]+(alphax[5]+alphay[0])*f[23])+8.660254037844386*(alphax[2]+alphay[1])*f[13]); 
  out[35] += 0.05*(17.32050807568877*alphax[2]*f[36]+(8.660254037844386*alphay[2]+17.32050807568877*alphax[1])*f[35]+17.32050807568877*alphax[5]*f[26]+(17.32050807568877*alphax[5]+8.660254037844387*alphay[0])*f[25]+(17.32050807568877*alphax[11]+7.745966692414834*alphay[5])*f[16]+19.36491673103709*(alphax[0]*f[16]+alphax[1]*f[9])+(19.36491673103709*alphax[2]+7.745966692414834*alphay[1])*f[8]+19.36491673103709*f[4]*alphax[5]); 
  out[36] += 0.05*((17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[36]+17.32050807568877*alphay[1]*f[35]+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[26]+17.32050807568877*alphay[5]*f[25]+(17.32050807568877*alphay[12]+7.745966692414834*alphax[5]+19.36491673103709*alphay[0])*f[16]+7.745966692414834*alphax[2]*f[9]+19.36491673103709*(alphay[1]*f[9]+alphay[2]*f[8]+f[4]*alphay[5])); 
  out[37] += 0.05*(17.32050807568877*alphax[5]*f[44]+17.32050807568877*alphax[1]*f[37]+19.36491673103709*(alphax[2]*f[31]+alphax[5]*f[18])+17.32050807568877*alphax[11]*f[17]+19.36491673103709*(alphax[0]*f[17]+alphax[1]*f[10])); 
  out[38] += 0.05*(17.32050807568877*alphay[5]*f[45]+17.32050807568877*alphay[2]*f[38]+19.36491673103709*alphay[1]*f[31]+17.32050807568877*alphay[12]*f[18]+19.36491673103709*(alphay[0]*f[18]+alphay[5]*f[17]+alphay[2]*f[10])); 
  out[39] += 0.05*(8.660254037844387*alphax[5]*f[46]+8.660254037844386*(alphax[2]*f[40]+alphax[1]*f[39])+8.660254037844387*alphax[0]*f[27]); 
  out[40] += 0.05*(8.660254037844387*alphay[5]*f[46]+8.660254037844386*(alphay[2]*f[40]+alphay[1]*f[39])+8.660254037844387*alphay[0]*f[27]); 
  out[41] += 0.05*(8.660254037844386*(alphay[2]+alphax[1])*f[41]+8.660254037844387*((alphay[5]+alphax[0])*f[29]+(alphax[5]+alphay[0])*f[28])+8.660254037844386*(alphax[2]+alphay[1])*f[14]); 
  out[42] += 0.05*(8.660254037844387*alphax[5]*f[47]+8.660254037844386*(alphax[2]*f[43]+alphax[1]*f[42])+8.660254037844387*alphax[0]*f[30]); 
  out[43] += 0.05*(8.660254037844387*alphay[5]*f[47]+8.660254037844386*(alphay[2]*f[43]+alphay[1]*f[42])+8.660254037844387*alphay[0]*f[30]); 
  out[44] += 0.05*(17.32050807568877*alphax[2]*f[45]+(8.660254037844386*alphay[2]+17.32050807568877*alphax[1])*f[44]+17.32050807568877*alphax[5]*f[38]+(17.32050807568877*alphax[5]+8.660254037844387*alphay[0])*f[37]+(17.32050807568877*alphax[11]+7.745966692414834*alphay[5])*f[31]+19.36491673103708*(alphax[0]*f[31]+alphax[1]*f[18])+(19.36491673103708*alphax[2]+7.745966692414834*alphay[1])*f[17]+19.36491673103708*alphax[5]*f[10]); 
  out[45] += 0.05*((17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[45]+17.32050807568877*alphay[1]*f[44]+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[38]+17.32050807568877*alphay[5]*f[37]+(17.32050807568877*alphay[12]+7.745966692414834*alphax[5]+19.36491673103708*alphay[0])*f[31]+7.745966692414834*alphax[2]*f[18]+19.36491673103708*(alphay[1]*f[18]+alphay[2]*f[17]+alphay[5]*f[10])); 
  out[46] += 0.05*(8.660254037844386*(alphay[2]+alphax[1])*f[46]+8.660254037844387*((alphay[5]+alphax[0])*f[40]+(alphax[5]+alphay[0])*f[39])+8.660254037844386*(alphax[2]+alphay[1])*f[27]); 
  out[47] += 0.05*(8.660254037844386*(alphay[2]+alphax[1])*f[47]+8.660254037844387*((alphay[5]+alphax[0])*f[43]+(alphax[5]+alphay[0])*f[42])+8.660254037844386*(alphax[2]+alphay[1])*f[30]); 
  return cflFreq; 
} 
double GyrokineticGenGeoVol2x2vSerP2_Bvarsxz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f, double *out) 
{ 
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

  double hamil[48]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil[11] = 2.0*(bmag[4]*wmu+phi[4]*q_); 
  hamil[12] = 2.0*phi[5]*q_; 
  hamil[13] = (0.5962847939999438*m_)/rdvpar2Sq; 
  hamil[19] = 2.0*phi[6]*q_; 
  hamil[20] = 2.0*phi[7]*q_; 
  hamil[25] = (1.154700538379251*bmag[4])/rdmu2; 

  double BstarXdBmag[48]; 

  double BstarYdBmag[48]; 
  BstarYdBmag[0] = -(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[1] = -(1.732050807568877*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[3] = -(1.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[6] = -(1.0*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarYdBmag[11] = -(1.732050807568877*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2*wvpar)/q_; 
  BstarYdBmag[21] = -(1.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[48]; 
  alphax[0] = ((3.872983346207417*((-0.1118033988749895*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.07142857142857142*b_z[4]*jacobTotInv[4]-0.1*b_z[1]*jacobTotInv[1])*hamil[19]+1.732050807568877*((-0.25*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[5]+hamil[2]*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])))-0.223606797749979*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])*hamil[5]))*rdx2*rdy2)/q_; 
  alphax[1] = ((3.872983346207417*((-0.1756910553749835*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-0.1*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[19]+1.732050807568877*(((-0.223606797749979*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.3928571428571428*b_z[4]*jacobTotInv[4]-0.45*b_z[1]*jacobTotInv[1]-0.25*b_z[0]*jacobTotInv[0])*hamil[5]+hamil[2]*((-0.223606797749979*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-0.25*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))))*rdx2*rdy2)/q_; 
  alphax[2] = (((-0.9682458365518543*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[20]+(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[12]))-0.8660254037844387*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])*hamil[20])*rdx2*rdy2)/q_; 
  alphax[5] = ((3.872983346207417*((-0.223606797749979*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.3928571428571428*b_z[4]*jacobTotInv[4]-0.45*b_z[1]*jacobTotInv[1]-0.25*b_z[0]*jacobTotInv[0])*hamil[20]+1.732050807568877*((-0.5*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-0.5590169943749475*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[12])*rdx2*rdy2)/q_; 
  alphax[11] = ((3.872983346207417*((-0.07142857142857142*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.2395787118749775*b_z[4]*jacobTotInv[4]+0.4472135954999579*((-0.3928571428571428*b_z[1]*jacobTotInv[1])-0.25*b_z[0]*jacobTotInv[0]))*hamil[19]+1.732050807568877*(((-0.3928571428571428*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-0.223606797749979*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[5]+hamil[2]*((-0.25*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.159719141249985*b_z[4]*jacobTotInv[4]-0.223606797749979*b_z[1]*jacobTotInv[1])))*rdx2*rdy2)/q_; 
  alphax[19] = ((1.732050807568877*((-0.8784552768749174*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-0.5*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[20]+3.872983346207417*((-0.25*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4]))-0.159719141249985*b_z[4]*jacobTotInv[4]-0.223606797749979*b_z[1]*jacobTotInv[1])*hamil[12])*rdx2*rdy2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*(2.23606797749979*alphax[11]-1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.23606797749979*alphax[11]+1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5590169943749475*alphax[11]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.75*alphax[19])+0.5590169943749475*alphax[11]-0.5809475019311124*alphax[5]-0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5590169943749475*alphax[11]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.75*alphax[19]+0.5590169943749475*alphax[11]+0.5809475019311124*alphax[5]+0.3354101966249685*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[48]; 
  alphay[0] = 1.732050807568877*(0.25*(((2.23606797749979*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[11]+hamil[1]*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*rdx2)/q_+((2.23606797749979*BstarYdBmag[3]*hamil[13]+BstarYdBmag[0]*hamil[3])*rdvpar2)/m_)+(0.5*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])*hamil[11]*rdx2)/q_)*rdy2; 
  alphay[1] = 1.732050807568877*(0.25*((hamil[1]*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*rdx2)/q_+((2.23606797749979*BstarYdBmag[6]*hamil[13]+BstarYdBmag[1]*hamil[3])*rdvpar2)/m_)+(((0.5*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+0.8784552768749174*b_z[4]*jacobTotInv[4]+0.25*(4.024922359499621*b_z[1]*jacobTotInv[1]+2.23606797749979*b_z[0]*jacobTotInv[0]))*hamil[11]+0.223606797749979*hamil[1]*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))*rdx2)/q_)*rdy2; 
  alphay[2] = ((0.25*(3.872983346207417*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[19]+1.732050807568877*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[5])+0.8660254037844387*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])*hamil[19])*rdx2*rdy2)/q_; 
  alphay[3] = (0.4330127018922193*(2.23606797749979*BstarYdBmag[0]*hamil[13]+BstarYdBmag[3]*hamil[3])*rdvpar2*rdy2)/m_; 
  alphay[4] = ((0.25*(3.872983346207417*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil[25]+1.732050807568877*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[8])+0.8660254037844387*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])*hamil[25])*rdx2*rdy2)/q_; 
  alphay[5] = ((3.872983346207417*(0.223606797749979*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+0.3928571428571428*b_z[4]*jacobTotInv[4]+0.45*b_z[1]*jacobTotInv[1]+0.25*b_z[0]*jacobTotInv[0])*hamil[19]+1.732050807568877*(0.223606797749979*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.25*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[5])*rdx2*rdy2)/q_; 
  alphay[6] = (0.4330127018922193*(2.23606797749979*BstarYdBmag[1]*hamil[13]+hamil[3]*BstarYdBmag[6])*rdvpar2*rdy2)/m_; 
  alphay[8] = ((3.872983346207417*(0.223606797749979*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+0.3928571428571428*b_z[4]*jacobTotInv[4]+0.45*b_z[1]*jacobTotInv[1]+0.25*b_z[0]*jacobTotInv[0])*hamil[25]+1.732050807568877*(0.223606797749979*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.25*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[8])*rdx2*rdy2)/q_; 
  alphay[11] = ((1.732050807568877*((0.8784552768749174*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.5*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[11]+hamil[1]*(0.25*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+0.159719141249985*b_z[4]*jacobTotInv[4]+0.223606797749979*b_z[1]*jacobTotInv[1]))*rdx2)/q_+(0.25*(3.872983346207417*hamil[13]*BstarYdBmag[21]+1.732050807568877*hamil[3]*BstarYdBmag[11])*rdvpar2)/m_)*rdy2; 
  alphay[12] = (0.4330127018922194*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil[20]*rdx2*rdy2)/q_; 
  alphay[13] = (0.8660254037844386*BstarYdBmag[3]*hamil[13]*rdvpar2*rdy2)/m_; 
  alphay[19] = ((1.732050807568877*(0.8784552768749174*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.5*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[19]+(0.2581988897471611*(1.677050983124842*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+1.5*b_z[1]*jacobTotInv[1])+0.276641667586244*b_z[4]*jacobTotInv[4])*hamil[5])*rdx2*rdy2)/q_; 
  alphay[20] = (1.732050807568877*(0.223606797749979*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.25*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[20]*rdx2*rdy2)/q_; 
  alphay[21] = (0.25*(1.732050807568877*hamil[3]*BstarYdBmag[21]+3.872983346207417*BstarYdBmag[11]*hamil[13])*rdvpar2*rdy2)/m_; 
  alphay[23] = (0.8660254037844387*BstarYdBmag[6]*hamil[13]*rdvpar2*rdy2)/m_; 
  alphay[25] = ((1.732050807568877*(0.8784552768749174*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+0.5*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamil[25]+3.872983346207417*(0.1118033988749895*(b_z[0]*jacobTotInv[4]+jacobTotInv[0]*b_z[4])+0.07142857142857142*b_z[4]*jacobTotInv[4]+0.1*b_z[1]*jacobTotInv[1])*hamil[8])*rdx2*rdy2)/q_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.125*(2.23606797749979*alphay[12]-1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(2.23606797749979*alphay[12]+1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])-0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*(alphay[8]+alphay[6])+0.5809475019311124*alphay[5]-0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.375*(alphay[25]+alphay[21])+0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]-0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])+0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*(alphay[8]+alphay[6])-0.5809475019311124*alphay[5]-0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])+0.375*alphay[23]-0.75*alphay[20]-0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])-0.375*alphay[23]+0.75*alphay[20]-0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])-0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.45*alphay[6]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.375*alphay[25]-0.375*alphay[21]+0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999999*alphay[25])+0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.45*alphay[6]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999998*alphay[23])-0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[6]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[6]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.2999999999999998*alphay[23])+0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[6]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[6]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]-0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.45*alphay[6]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.375*alphay[25])+0.375*alphay[21]+0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]+0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.45*alphay[6]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]+0.375*alphay[23]-0.75*alphay[20]-0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]-0.375*alphay[23]+0.75*alphay[20]-0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]-0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]-0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*(alphay[8]+alphay[6])+0.5809475019311124*alphay[5]+0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.375*(alphay[25]+alphay[21]))+0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]+0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.2999999999999999*alphay[25]+0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]-0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*(alphay[8]+alphay[6])-0.5809475019311124*alphay[5]+0.3354101966249685*(alphay[4]+alphay[3])-0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])-0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*(alphay[8]+alphay[6])-0.5809475019311124*alphay[5]-0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.375*(alphay[25]+alphay[21])-0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]-0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])+0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*(alphay[8]+alphay[6])+0.5809475019311124*alphay[5]-0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])+0.375*alphay[23]-0.75*alphay[20]+0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])-0.375*alphay[23]+0.75*alphay[20]+0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])-0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.45*alphay[6]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.375*alphay[25]-0.375*alphay[21]-0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999999*alphay[25])+0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.45*alphay[6]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[4]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999998*alphay[23])-0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[6]-0.5809475019311124*alphay[5]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[6]+0.5809475019311124*alphay[5]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.2999999999999998*alphay[23])+0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[6]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[6]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]-0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]+0.45*alphay[6]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.375*alphay[25])+0.375*alphay[21]-0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]+0.2999999999999998*alphay[23]-0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]-0.45*alphay[6]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]-0.3354101966249685*alphay[3]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]+0.375*alphay[23]-0.75*alphay[20]+0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*alphay[8]-0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]-0.375*alphay[23]+0.75*alphay[20]+0.3872983346207416*alphay[19]-0.2795084971874737*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*alphay[8]+0.5809475019311124*alphay[5]+0.3354101966249685*alphay[4]+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]-0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]-0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]-0.45*(alphay[8]+alphay[6])-0.5809475019311124*alphay[5]+0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]-0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.375*(alphay[25]+alphay[21]))-0.4841229182759271*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]-0.2795084971874737*alphay[11]+0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.2999999999999999*alphay[25]+0.2999999999999998*alphay[23]+0.2999999999999999*alphay[21]+0.75*alphay[20]+0.3872983346207416*alphay[19]+0.223606797749979*alphay[13]+0.5590169943749475*alphay[12]+0.223606797749979*alphay[11]+0.45*(alphay[8]+alphay[6])+0.5809475019311124*alphay[5]+0.3354101966249685*(alphay[4]+alphay[3])+0.4330127018922193*alphay[2]+0.3354101966249685*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[48]; 
  alphavpar[0] = -(0.25*(1.732050807568877*BstarYdBmag[11]*hamil[19]+1.732050807568877*(BstarYdBmag[1]*hamil[5]+BstarYdBmag[0]*hamil[2]))*rdvpar2*rdy2)/m_; 
  alphavpar[1] = ((1.732050807568877*((-0.223606797749979*hamil[5]*BstarYdBmag[11])-0.25*(BstarYdBmag[0]*hamil[5]+BstarYdBmag[1]*hamil[2]))-0.3872983346207417*BstarYdBmag[1]*hamil[19])*rdvpar2*rdy2)/m_; 
  alphavpar[2] = -(0.9682458365518543*(BstarYdBmag[1]*hamil[20]+BstarYdBmag[0]*hamil[12])*rdvpar2*rdy2)/m_; 
  alphavpar[3] = -(0.4330127018922193*(hamil[19]*BstarYdBmag[21]+hamil[5]*BstarYdBmag[6]+hamil[2]*BstarYdBmag[3])*rdvpar2*rdy2)/m_; 
  alphavpar[5] = (((-0.9682458365518543*(BstarYdBmag[0]*hamil[20]+BstarYdBmag[1]*hamil[12]))-0.8660254037844387*BstarYdBmag[11]*hamil[20])*rdvpar2*rdy2)/m_; 
  alphavpar[6] = (((-0.3872983346207417*(hamil[5]*BstarYdBmag[21]+BstarYdBmag[6]*hamil[19]))-0.4330127018922193*(hamil[2]*BstarYdBmag[6]+BstarYdBmag[3]*hamil[5]))*rdvpar2*rdy2)/m_; 
  alphavpar[7] = -(0.9682458365518543*(BstarYdBmag[6]*hamil[20]+BstarYdBmag[3]*hamil[12])*rdvpar2*rdy2)/m_; 
  alphavpar[11] = ((3.872983346207417*((-0.07142857142857142*BstarYdBmag[11])-0.1118033988749895*BstarYdBmag[0])*hamil[19]+1.732050807568877*((-0.25*hamil[2]*BstarYdBmag[11])-0.223606797749979*BstarYdBmag[1]*hamil[5]))*rdvpar2*rdy2)/m_; 
  alphavpar[15] = (((-0.8660254037844386*hamil[20]*BstarYdBmag[21])-0.9682458365518543*(BstarYdBmag[3]*hamil[20]+BstarYdBmag[6]*hamil[12]))*rdvpar2*rdy2)/m_; 
  alphavpar[19] = (((-0.8660254037844386*BstarYdBmag[1]*hamil[20])-0.9682458365518543*BstarYdBmag[11]*hamil[12])*rdvpar2*rdy2)/m_; 
  alphavpar[21] = (((-0.4330127018922193*(hamil[2]*BstarYdBmag[21]+BstarYdBmag[3]*hamil[19]))-0.276641667586244*hamil[19]*BstarYdBmag[21]-0.3872983346207417*hamil[5]*BstarYdBmag[6])*rdvpar2*rdy2)/m_; 
  alphavpar[32] = (3.872983346207417*((-0.25*hamil[12]*BstarYdBmag[21])-0.223606797749979*BstarYdBmag[6]*hamil[20])*rdvpar2*rdy2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphavpar[3]-1.0*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphavpar[3]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.6495190528383289*alphavpar[32])+0.4841229182759271*alphavpar[21]+0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.3872983346207416*alphavpar[21])+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[6]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.3872983346207416*alphavpar[21])+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[6]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.6495190528383289*alphavpar[32]+0.4841229182759271*alphavpar[21]-0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]-0.5809475019311124*alphavpar[7]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.6495190528383289*alphavpar[32])+0.4841229182759271*alphavpar[21]+0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.5196152422706631*alphavpar[32]-0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.3872983346207416*alphavpar[21])+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[6]-0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.3872983346207416*alphavpar[21])+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[6]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.6495190528383289*alphavpar[32]+0.4841229182759271*alphavpar[21]-0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]-0.5809475019311124*alphavpar[7]-0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.5196152422706631*alphavpar[32])-0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.6495190528383289*alphavpar[32]-0.4841229182759271*alphavpar[21]+0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.3872983346207416*alphavpar[21]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[6]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.3872983346207416*alphavpar[21]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[6]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.6495190528383289*alphavpar[32])-0.4841229182759271*alphavpar[21]-0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]+0.5809475019311124*alphavpar[7]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.6495190528383289*alphavpar[32]-0.4841229182759271*alphavpar[21]+0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.5196152422706631*alphavpar[32])+0.3872983346207416*alphavpar[21]-0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[7]+0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[2]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.3872983346207416*alphavpar[21]+0.223606797749979*alphavpar[11]-0.5809475019311124*alphavpar[6]+0.4330127018922193*alphavpar[3]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.3872983346207416*alphavpar[21]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[6]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]-0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*alphavpar[7]-0.5809475019311124*alphavpar[6]-0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]-0.3354101966249685*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.6495190528383289*alphavpar[32])-0.4841229182759271*alphavpar[21]-0.375*alphavpar[19]-0.2795084971874737*alphavpar[11]+0.5809475019311124*alphavpar[7]+0.4330127018922193*alphavpar[3]+0.3354101966249685*alphavpar[2]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.5196152422706631*alphavpar[32]+0.3872983346207416*alphavpar[21]+0.2999999999999999*alphavpar[19]+0.7794228634059945*alphavpar[15]+0.223606797749979*alphavpar[11]+0.5809475019311124*(alphavpar[7]+alphavpar[6])+0.45*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.3354101966249685*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[19]*f[19]+alphax[11]*f[11]+alphax[5]*f[5]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[25]*f[25]+alphay[23]*f[23]+alphay[21]*f[21]+alphay[20]*f[20]+alphay[19]*f[19]+alphay[13]*f[13]+alphay[12]*f[12]+alphay[11]*f[11]+alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphavpar[32]*f[32]+alphavpar[21]*f[21]+alphavpar[19]*f[19]+alphavpar[15]*f[15]+alphavpar[11]*f[11]+alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[5] += 0.05*(7.745966692414834*(alphay[8]*f[25]+f[8]*alphay[25])+8.660254037844387*(alphay[13]*f[23]+f[13]*alphay[23])+7.745966692414834*(alphay[6]*f[21]+f[6]*alphay[21])+(8.660254037844387*alphay[12]+7.745966692414834*alphax[5])*f[20]+8.660254037844387*(f[12]*alphay[20]+alphax[11]*f[19])+7.745966692414834*(alphay[5]*f[19]+f[5]*alphay[19])+8.660254037844387*f[11]*alphax[19]+7.745966692414834*(alphax[2]*f[12]+alphay[1]*f[11]+f[1]*alphay[11])+8.660254037844386*(alphay[4]*f[8]+f[4]*alphay[8]+alphay[3]*f[6]+f[3]*alphay[6]+(alphay[2]+alphax[1])*f[5]+f[2]*alphay[5]+f[1]*alphax[5]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1])); 
  out[6] += 0.05*(8.660254037844387*alphax[19]*f[32]+7.745966692414834*(alphavpar[15]*f[32]+f[15]*alphavpar[32])+8.660254037844387*alphax[11]*f[21]+7.745966692414834*(alphavpar[6]*f[21]+f[6]*alphavpar[21]+alphavpar[5]*f[19]+f[5]*alphavpar[19])+8.660254037844386*((alphavpar[7]+alphax[5])*f[15]+f[7]*alphavpar[15])+7.745966692414834*(alphavpar[1]*f[11]+f[1]*alphavpar[11])+8.660254037844386*(alphax[2]*f[7]+(alphavpar[3]+alphax[1])*f[6]+f[3]*alphavpar[6]+alphavpar[2]*f[5]+f[2]*alphavpar[5]+alphax[0]*f[3]+alphavpar[0]*f[1]+f[0]*alphavpar[1])); 
  out[7] += 0.05*(8.660254037844387*alphay[25]*f[37]+(8.660254037844387*alphay[20]+7.745966692414834*alphavpar[15])*f[33]+8.660254037844387*((alphavpar[21]+alphay[19])*f[32]+f[21]*alphavpar[32])+7.745966692414834*(alphay[6]*f[23]+f[6]*alphay[23])+(8.660254037844387*alphay[12]+7.745966692414834*alphavpar[7])*f[22]+8.660254037844387*(alphay[11]*f[21]+f[11]*alphay[21])+7.745966692414834*alphavpar[5]*f[20]+8.660254037844387*(alphavpar[11]*f[19]+f[11]*alphavpar[19])+8.660254037844386*(alphay[8]*f[17]+(alphavpar[6]+alphay[5])*f[15]+f[6]*alphavpar[15])+7.745966692414834*(alphay[3]*f[13]+f[3]*alphay[13]+alphavpar[2]*f[12])+8.660254037844386*(alphay[4]*f[10]+(alphavpar[3]+alphay[2])*f[7]+f[3]*alphavpar[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphay[0]*f[3]+f[0]*alphay[3]+alphavpar[0]*f[2]+f[0]*alphavpar[2])); 
  out[8] += 0.05*(8.660254037844387*(alphax[19]*f[35]+alphax[11]*f[25])+8.660254037844386*(alphax[5]*f[16]+alphax[2]*f[9]+alphax[1]*f[8]+alphax[0]*f[4])); 
  out[9] += 0.05*(8.660254037844387*(alphay[23]*f[39]+alphay[21]*f[37]+alphay[20]*f[36]+alphay[19]*f[35])+7.745966692414834*alphay[8]*f[28]+8.660254037844387*(alphay[13]*f[27]+alphay[12]*f[26]+alphay[11]*f[25]+f[11]*alphay[25])+8.660254037844386*(alphay[6]*f[17]+alphay[5]*f[16])+7.745966692414834*alphay[4]*f[14]+8.660254037844386*(alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4])); 
  out[10] += 0.05*(8.660254037844387*(alphavpar[32]*f[44]+alphavpar[21]*f[37]+alphavpar[19]*f[35])+8.660254037844386*alphavpar[15]*f[31]+8.660254037844387*alphavpar[11]*f[25]+8.660254037844386*(alphavpar[7]*f[18]+alphavpar[6]*f[17]+alphavpar[5]*f[16]+alphavpar[3]*f[10]+alphavpar[2]*f[9]+alphavpar[1]*f[8]+alphavpar[0]*f[4])); 
  out[11] += 0.05*(17.32050807568877*(alphax[5]*f[19]+f[5]*alphax[19])+17.32050807568877*(alphax[1]*f[11]+f[1]*alphax[11])+19.36491673103709*(alphax[2]*f[5]+f[2]*alphax[5]+alphax[0]*f[1]+f[0]*alphax[1])); 
  out[12] += 0.05*(19.36491673103708*(alphay[25]*f[35]+alphay[23]*f[34]+alphay[21]*f[32]+alphay[13]*f[24])+17.32050807568877*(alphay[5]*f[20]+f[5]*alphay[20])+19.36491673103708*(alphay[11]*f[19]+f[11]*alphay[19])+19.36491673103709*(alphay[8]*f[16]+alphay[6]*f[15])+17.32050807568877*(alphay[2]*f[12]+f[2]*alphay[12])+19.36491673103709*(alphay[4]*f[9]+alphay[3]*f[7]+alphay[1]*f[5]+f[1]*alphay[5]+alphay[0]*f[2]+f[0]*alphay[2])); 
  out[13] += 0.05*(17.32050807568877*alphavpar[15]*f[34]+19.36491673103708*(alphavpar[19]*f[32]+f[19]*alphavpar[32])+17.32050807568877*(alphavpar[7]*f[24]+alphavpar[6]*f[23])+19.36491673103708*(alphavpar[11]*f[21]+f[11]*alphavpar[21])+19.36491673103709*(alphavpar[5]*f[15]+f[5]*alphavpar[15])+17.32050807568877*alphavpar[3]*f[13]+19.36491673103709*(alphavpar[2]*f[7]+f[2]*alphavpar[7]+alphavpar[1]*f[6]+f[1]*alphavpar[6]+alphavpar[0]*f[3]+f[0]*alphavpar[3])); 
  out[15] += 0.05*(7.745966692414834*alphay[8]*f[37]+(6.928203230275509*alphavpar[32]+8.660254037844386*alphay[12]+7.745966692414834*(alphavpar[7]+alphax[5]))*f[33]+8.660254037844386*alphax[11]*f[32]+7.745966692414834*((alphavpar[6]+alphay[5])*f[32]+f[6]*alphavpar[32]+f[17]*alphay[25])+(6.928203230275509*alphay[21]+7.745966692414834*alphay[3])*f[23]+(6.928203230275509*f[21]+7.745966692414834*f[3])*alphay[23]+(8.660254037844386*alphay[20]+7.745966692414834*(alphavpar[15]+alphax[2]))*f[22]+8.660254037844386*alphax[19]*f[21]+7.745966692414834*((alphavpar[15]+alphay[1])*f[21]+f[1]*alphay[21]+f[15]*alphavpar[21])+6.928203230275509*alphavpar[19]*f[20]+7.745966692414834*(alphavpar[2]*f[20]+alphavpar[1]*f[19]+f[15]*alphay[19]+f[1]*alphavpar[19])+8.660254037844386*(alphay[4]*f[17]+(alphavpar[3]+alphay[2]+alphax[1])*f[15]+f[3]*alphavpar[15])+7.745966692414834*(alphay[6]*f[13]+f[6]*alphay[13]+alphavpar[5]*f[12]+(alphay[6]+alphavpar[5])*f[11]+f[6]*alphay[11]+f[5]*alphavpar[11])+8.660254037844386*(alphay[8]*f[10]+(alphavpar[6]+alphay[5]+alphax[0])*f[7]+f[6]*(alphavpar[7]+alphax[5]+alphay[0])+f[0]*alphay[6]+alphavpar[0]*f[5]+f[0]*alphavpar[5]+(alphax[2]+alphay[1])*f[3]+f[1]*alphay[3]+alphavpar[1]*f[2]+f[1]*alphavpar[2])); 
  out[16] += 0.05*(8.660254037844386*alphay[13]*f[39]+7.745966692414834*alphay[6]*f[37]+(8.660254037844386*alphay[12]+7.745966692414834*alphax[5])*f[36]+(8.660254037844386*alphax[11]+7.745966692414834*alphay[5])*f[35]+(6.928203230275509*alphay[25]+7.745966692414834*alphay[4])*f[28]+8.660254037844386*alphay[23]*f[27]+(8.660254037844386*alphay[20]+7.745966692414834*alphax[2])*f[26]+8.660254037844386*alphax[19]*f[25]+7.745966692414834*(alphay[1]*f[25]+f[1]*alphay[25]+f[17]*alphay[21]+f[16]*alphay[19])+8.660254037844386*(alphay[3]*f[17]+(alphay[2]+alphax[1])*f[16])+7.745966692414834*(alphay[8]*(f[14]+f[11])+f[8]*alphay[11])+8.660254037844386*(alphay[6]*f[10]+(alphay[5]+alphax[0])*f[9]+(alphax[5]+alphay[0])*f[8]+f[0]*alphay[8]+(alphax[2]+alphay[1])*f[4]+f[1]*alphay[4])); 
  out[17] += 0.05*((8.660254037844386*alphax[19]+7.745966692414834*alphavpar[15])*f[44]+8.660254037844386*alphax[11]*f[37]+7.745966692414834*(alphavpar[6]*f[37]+alphavpar[5]*f[35])+f[31]*(7.745966692414834*alphavpar[32]+8.660254037844386*(alphavpar[7]+alphax[5]))+7.745966692414834*(alphavpar[1]*f[25]+f[17]*alphavpar[21]+f[16]*alphavpar[19])+8.660254037844386*((alphavpar[15]+alphax[2])*f[18]+(alphavpar[3]+alphax[1])*f[17]+alphavpar[2]*f[16])+7.745966692414834*f[8]*alphavpar[11]+8.660254037844386*((alphavpar[6]+alphax[0])*f[10]+alphavpar[5]*f[9]+alphavpar[0]*f[8]+alphavpar[1]*f[4])); 
  out[18] += 0.05*((8.660254037844386*alphay[20]+7.745966692414834*alphavpar[15])*f[45]+8.660254037844386*(alphavpar[21]+alphay[19])*f[44]+7.745966692414834*(alphay[8]*f[42]+alphay[6]*f[39])+(8.660254037844386*alphay[12]+7.745966692414834*alphavpar[7])*f[38]+8.660254037844386*(alphavpar[32]+alphay[11])*f[37]+7.745966692414834*alphavpar[5]*f[36]+8.660254037844386*(alphavpar[11]*f[35]+(alphavpar[6]+alphay[5])*f[31])+7.745966692414834*(alphay[4]*f[30]+alphay[3]*f[27]+alphavpar[2]*f[26])+8.660254037844386*((alphay[21]+alphavpar[19])*f[25]+f[21]*alphay[25])+7.745966692414834*f[17]*alphay[23]+8.660254037844386*((alphavpar[3]+alphay[2])*f[18]+(alphavpar[15]+alphay[1])*f[17]+alphavpar[1]*f[16])+7.745966692414834*f[10]*alphay[13]+8.660254037844386*((alphavpar[7]+alphay[0])*f[10]+alphavpar[0]*f[9]+(alphay[6]+alphavpar[5])*f[8]+f[6]*alphay[8]+(alphay[3]+alphavpar[2])*f[4]+f[3]*alphay[4])); 
  out[19] += 0.007142857142857143*(38.72983346207417*alphay[25]*f[25]+60.6217782649107*(alphay[4]*f[25]+f[4]*alphay[25])+54.22176684690384*alphay[23]*f[23]+38.72983346207417*alphay[21]*f[21]+60.6217782649107*(alphay[3]*f[21]+f[3]*alphay[21])+(54.22176684690384*alphay[20]+108.4435336938077*alphax[19]+121.2435565298214*alphax[2])*f[20]+(38.72983346207417*alphay[19]+60.6217782649107*alphay[2]+121.2435565298214*alphax[1])*f[19]+60.6217782649107*f[2]*alphay[19]+121.2435565298214*f[1]*alphax[19]+121.2435565298214*alphax[5]*f[12]+(38.72983346207417*alphay[11]+121.2435565298214*alphax[5])*f[11]+60.62177826491071*(alphay[0]*f[11]+f[0]*alphay[11])+121.2435565298214*f[5]*alphax[11]+54.22176684690384*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5])+135.5544171172596*(alphax[0]*f[5]+f[0]*alphax[5]+alphax[1]*f[2])+f[1]*(135.5544171172596*alphax[2]+54.22176684690384*alphay[1])); 
  out[20] += 0.05*(17.32050807568877*alphay[8]*f[35]+19.36491673103708*alphay[13]*f[34]+17.32050807568877*alphay[6]*f[32]+17.32050807568877*f[16]*alphay[25]+19.36491673103708*alphay[23]*f[24]+17.32050807568877*f[15]*alphay[21]+(15.49193338482967*alphay[19]+17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[20]+(15.49193338482967*f[19]+17.32050807568877*f[2])*alphay[20]+7.745966692414834*alphax[19]*f[19]+17.32050807568877*(alphay[1]*f[19]+f[1]*alphay[19])+19.36491673103708*(alphay[4]*f[16]+alphay[3]*f[15])+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[12]+17.32050807568877*(f[5]*alphay[12]+alphay[5]*f[11]+f[5]*alphay[11])+19.36491673103708*(alphay[8]*f[9]+alphay[6]*f[7])+7.745966692414834*alphax[5]*f[5]+19.36491673103708*(alphay[0]*f[5]+f[0]*alphay[5])+7.745966692414834*alphax[2]*f[2]+19.36491673103708*(alphay[1]*f[2]+f[1]*alphay[2])); 
  out[21] += 0.007142857142857143*((38.72983346207417*alphavpar[32]+60.62177826491071*alphavpar[7]+121.2435565298214*alphax[5])*f[32]+60.62177826491071*f[7]*alphavpar[32]+(38.72983346207417*alphavpar[21]+60.6217782649107*alphavpar[3]+121.2435565298214*alphax[1])*f[21]+60.6217782649107*f[3]*alphavpar[21]+(38.72983346207417*alphavpar[19]+60.6217782649107*alphavpar[2])*f[19]+121.2435565298214*f[15]*alphax[19]+60.6217782649107*f[2]*alphavpar[19]+(54.22176684690384*alphavpar[15]+135.5544171172596*alphax[2])*f[15]+(38.72983346207417*alphavpar[11]+60.62177826491071*alphavpar[0])*f[11]+121.2435565298214*f[6]*alphax[11]+60.62177826491071*f[0]*alphavpar[11]+135.5544171172596*alphax[5]*f[7]+(54.22176684690384*alphavpar[6]+135.5544171172596*alphax[0])*f[6]+54.22176684690384*alphavpar[5]*f[5]+135.5544171172596*alphax[1]*f[3]+54.22176684690384*alphavpar[1]*f[1]); 
  out[22] += 0.05*(19.36491673103708*alphay[25]*f[44]+17.32050807568877*alphay[6]*f[34]+(8.660254037844387*alphavpar[6]+17.32050807568877*alphay[5])*f[33]+7.745966692414834*alphavpar[32]*f[32]+19.36491673103708*(alphay[11]*f[32]+alphay[8]*f[31])+17.32050807568877*(alphay[3]*f[24]+f[15]*alphay[23])+(8.660254037844386*alphavpar[3]+17.32050807568877*alphay[2])*f[22]+19.36491673103708*(alphay[19]*f[21]+f[19]*alphay[21])+8.660254037844386*alphavpar[1]*f[20]+17.32050807568877*f[15]*alphay[20]+7.745966692414834*alphavpar[19]*f[19]+19.36491673103708*alphay[4]*f[18]+(7.745966692414834*alphavpar[15]+19.36491673103708*alphay[1])*f[15]+17.32050807568877*f[7]*alphay[13]+8.660254037844387*alphavpar[0]*f[12]+f[7]*(17.32050807568877*alphay[12]+7.745966692414834*alphavpar[7])+19.36491673103708*(alphay[0]*f[7]+alphay[5]*f[6])+f[5]*(19.36491673103708*alphay[6]+7.745966692414834*alphavpar[5])+19.36491673103708*alphay[2]*f[3]+f[2]*(19.36491673103708*alphay[3]+7.745966692414834*alphavpar[2])); 
  out[23] += 0.05*((15.49193338482967*alphavpar[32]+17.32050807568877*alphavpar[7]+8.660254037844387*alphax[5])*f[34]+17.32050807568877*(alphavpar[5]*f[32]+f[5]*alphavpar[32])+(17.32050807568877*alphavpar[15]+8.660254037844386*alphax[2])*f[24]+(15.49193338482967*alphavpar[21]+17.32050807568877*alphavpar[3]+8.660254037844386*alphax[1])*f[23]+17.32050807568877*(alphavpar[1]*f[21]+f[1]*alphavpar[21]+alphavpar[15]*f[19]+f[15]*alphavpar[19])+19.36491673103708*(alphavpar[2]*f[15]+f[2]*alphavpar[15])+(17.32050807568877*alphavpar[6]+8.660254037844387*alphax[0])*f[13]+17.32050807568877*(alphavpar[6]*f[11]+f[6]*alphavpar[11])+19.36491673103708*(alphavpar[5]*f[7]+f[5]*alphavpar[7]+alphavpar[0]*f[6]+f[0]*alphavpar[6]+alphavpar[1]*f[3]+f[1]*alphavpar[3])); 
  out[24] += 0.007142857142857143*(60.62177826491071*alphay[8]*f[39]+(121.2435565298214*alphavpar[6]+60.62177826491071*alphay[5])*f[34]+121.2435565298214*alphavpar[5]*f[33]+135.5544171172596*(alphavpar[11]*f[32]+f[11]*alphavpar[32])+60.6217782649107*alphay[4]*f[27]+(121.2435565298214*alphavpar[3]+60.6217782649107*alphay[2])*f[24]+(38.72983346207417*alphay[23]+121.2435565298214*alphavpar[15])*f[23]+60.6217782649107*(alphay[1]*f[23]+f[1]*alphay[23])+121.2435565298214*alphavpar[2]*f[22]+54.22176684690384*alphay[21]*f[21]+135.5544171172596*(alphavpar[19]*f[21]+f[19]*alphavpar[21])+121.2435565298214*alphavpar[15]*f[20]+135.5544171172596*(alphavpar[1]*f[15]+f[1]*alphavpar[15])+(38.72983346207417*alphay[13]+121.2435565298214*alphavpar[7])*f[13]+60.62177826491071*(alphay[0]*f[13]+f[0]*alphay[13])+121.2435565298214*alphavpar[7]*f[12]+135.5544171172596*(alphavpar[0]*f[7]+f[0]*alphavpar[7])+54.22176684690384*alphay[6]*f[6]+135.5544171172596*(alphavpar[5]*f[6]+f[5]*alphavpar[6])+54.22176684690384*alphay[3]*f[3]+135.5544171172596*(alphavpar[2]*f[3]+f[2]*alphavpar[3])); 
  out[25] += 0.05*(17.32050807568877*alphax[5]*f[35]+17.32050807568877*alphax[1]*f[25]+f[16]*(17.32050807568877*alphax[19]+19.36491673103708*alphax[2])+17.32050807568877*f[8]*alphax[11]+19.36491673103708*(alphax[5]*f[9]+alphax[0]*f[8]+alphax[1]*f[4])); 
  out[26] += 0.05*(19.36491673103708*(alphay[23]*f[46]+alphay[21]*f[44])+17.32050807568877*alphay[8]*f[41]+19.36491673103708*alphay[13]*f[40]+17.32050807568877*alphay[5]*f[36]+19.36491673103708*(alphay[11]*f[35]+alphay[6]*f[31])+17.32050807568877*(alphay[4]*f[29]+alphay[2]*f[26])+19.36491673103708*(alphay[19]*f[25]+f[19]*alphay[25])+17.32050807568877*f[16]*alphay[20]+19.36491673103708*(alphay[3]*f[18]+alphay[1]*f[16])+17.32050807568877*f[9]*alphay[12]+19.36491673103708*(alphay[0]*f[9]+alphay[5]*f[8]+f[5]*alphay[8]+alphay[2]*f[4]+f[2]*alphay[4])); 
  out[27] += 0.05*(17.32050807568877*alphavpar[15]*f[46]+19.36491673103708*alphavpar[19]*f[44]+17.32050807568877*(alphavpar[7]*f[40]+alphavpar[6]*f[39])+19.36491673103708*(alphavpar[11]*f[37]+alphavpar[32]*f[35]+alphavpar[5]*f[31])+17.32050807568877*alphavpar[3]*f[27]+19.36491673103708*(alphavpar[21]*f[25]+alphavpar[2]*f[18]+alphavpar[1]*f[17]+alphavpar[15]*f[16]+alphavpar[0]*f[10]+alphavpar[7]*f[9]+alphavpar[6]*f[8]+alphavpar[3]*f[4])); 
  out[28] += 0.05*(8.660254037844387*alphax[5]*f[41]+8.660254037844386*(alphax[2]*f[29]+alphax[1]*f[28])+8.660254037844387*alphax[0]*f[14]); 
  out[29] += 0.05*(8.660254037844387*(alphay[6]*f[42]+alphay[5]*f[41])+8.660254037844386*(alphay[3]*f[30]+alphay[2]*f[29]+alphay[1]*f[28])+7.745966692414834*alphay[25]*f[25]+8.660254037844387*alphay[0]*f[14]+7.745966692414834*(alphay[8]*f[8]+alphay[4]*f[4])); 
  out[30] += 0.05*(8.660254037844386*alphavpar[15]*f[47]+8.660254037844387*(alphavpar[7]*f[43]+alphavpar[6]*f[42]+alphavpar[5]*f[41])+8.660254037844386*(alphavpar[3]*f[30]+alphavpar[2]*f[29]+alphavpar[1]*f[28])+8.660254037844387*alphavpar[0]*f[14]); 
  out[31] += 0.01*((34.64101615137755*alphavpar[32]+43.30127018922195*alphay[12]+38.72983346207417*(alphavpar[7]+alphax[5]))*f[45]+(43.30127018922195*alphax[11]+38.72983346207417*(alphavpar[6]+alphay[5]))*f[44]+(34.64101615137755*alphay[25]+38.72983346207418*alphay[4])*f[42]+(34.64101615137755*alphay[21]+38.72983346207418*alphay[3])*f[39]+(43.30127018922195*alphay[20]+38.72983346207418*(alphavpar[15]+alphax[2]))*f[38]+(34.64101615137755*alphay[23]+43.30127018922195*alphax[19]+38.72983346207418*(alphavpar[15]+alphay[1]))*f[37]+34.64101615137755*alphavpar[19]*f[36]+38.72983346207418*(alphavpar[2]*f[36]+alphavpar[1]*f[35]+f[17]*alphavpar[32])+(38.72983346207417*(alphavpar[21]+alphay[19])+43.30127018922193*(alphavpar[3]+alphay[2]+alphax[1]))*f[31]+38.72983346207417*(alphay[8]*f[30]+alphay[6]*f[27]+alphavpar[5]*f[26]+(alphay[6]+alphavpar[5])*f[25]+f[6]*alphay[25]+f[10]*alphay[23]+alphay[8]*f[21]+f[8]*(alphay[21]+alphavpar[19]))+43.30127018922193*(alphavpar[6]+alphay[5]+alphax[0])*f[18]+(38.72983346207418*(alphay[13]+alphay[11])+43.30127018922193*(alphavpar[7]+alphax[5]+alphay[0]))*f[17]+38.72983346207418*alphavpar[11]*f[16]+43.30127018922193*(alphavpar[0]*f[16]+f[10]*(alphavpar[15]+alphax[2]+alphay[1])+alphavpar[1]*f[9]+(alphay[3]+alphavpar[2])*f[8]+f[3]*alphay[8]+alphay[4]*f[6]+f[4]*(alphay[6]+alphavpar[5]))); 
  out[32] += 0.001428571428571429*((193.6491673103708*alphay[25]+303.1088913245535*alphay[4])*f[37]+(271.1088342345192*alphay[20]+542.2176684690384*alphax[19]+242.4871130596428*alphavpar[15]+606.217782649107*alphax[2])*f[33]+(193.6491673103708*(alphavpar[21]+alphay[19])+303.1088913245535*(alphavpar[3]+alphay[2])+606.217782649107*alphax[1])*f[32]+(271.1088342345192*f[22]+193.6491673103708*f[21]+303.1088913245535*f[3])*alphavpar[32]+303.1088913245536*f[10]*alphay[25]+242.4871130596428*(alphay[6]*f[23]+f[6]*alphay[23])+606.2177826491072*alphax[5]*f[22]+(271.1088342345192*alphay[13]+193.6491673103708*alphay[11]+303.1088913245536*alphavpar[7]+606.2177826491072*alphax[5]+303.1088913245536*alphay[0])*f[21]+(271.1088342345192*f[13]+193.6491673103708*f[11])*alphay[21]+303.1088913245536*(f[0]*alphay[21]+f[7]*alphavpar[21])+242.4871130596428*alphavpar[5]*f[20]+193.6491673103708*alphavpar[11]*f[19]+303.1088913245536*(alphavpar[0]*f[19]+f[7]*alphay[19])+606.2177826491072*f[6]*alphax[19]+(271.1088342345192*f[12]+193.6491673103708*f[11]+303.1088913245536*f[0])*alphavpar[19]+271.1088342345192*alphay[8]*f[17]+(606.217782649107*alphax[11]+271.1088342345192*(alphavpar[6]+alphay[5])+677.7720855862981*alphax[0])*f[15]+271.1088342345192*f[6]*alphavpar[15]+303.1088913245535*((alphay[3]+alphavpar[2])*f[11]+f[3]*alphay[11]+f[2]*alphavpar[11])+677.7720855862981*(alphax[1]*f[7]+alphax[2]*f[6])+271.1088342345192*(alphay[1]*f[6]+f[1]*alphay[6]+alphavpar[1]*f[5])+677.7720855862981*f[3]*alphax[5]+271.1088342345192*f[1]*alphavpar[5]); 
  out[33] += 0.01*(86.60254037844389*alphay[8]*f[44]+(77.45966692414834*alphay[21]+86.60254037844386*alphay[3])*f[34]+(38.72983346207417*alphavpar[21]+77.45966692414834*alphay[19]+43.30127018922193*alphavpar[3]+86.60254037844386*alphay[2]+43.30127018922193*alphax[1])*f[33]+(77.45966692414834*(alphay[23]+alphay[20])+38.72983346207417*alphax[19]+34.64101615137754*alphavpar[15]+86.60254037844386*alphay[1])*f[32]+34.64101615137754*f[15]*alphavpar[32]+(86.60254037844389*alphay[25]+96.82458365518544*alphay[4])*f[31]+86.60254037844389*(alphay[6]*f[24]+f[7]*alphay[23])+(43.30127018922195*alphavpar[6]+86.60254037844389*alphay[5]+43.30127018922195*alphax[0])*f[22]+86.60254037844389*(alphay[5]*f[21]+f[5]*alphay[21])+(38.72983346207417*alphavpar[11]+43.30127018922195*alphavpar[0])*f[20]+86.60254037844389*f[7]*alphay[20]+(86.60254037844389*alphay[6]+34.64101615137755*alphavpar[5])*f[19]+86.60254037844389*f[6]*alphay[19]+34.64101615137755*f[5]*alphavpar[19]+96.82458365518544*alphay[8]*f[18]+(86.60254037844386*(alphay[13]+alphay[12]+alphay[11])+38.72983346207418*(alphavpar[7]+alphax[5])+96.82458365518544*alphay[0])*f[15]+38.72983346207418*f[7]*alphavpar[15]+43.30127018922193*alphavpar[1]*f[12]+38.72983346207418*alphax[2]*f[7]+96.82458365518544*(alphay[1]*f[7]+alphay[2]*f[6]+f[2]*alphay[6])+(96.82458365518544*alphay[3]+38.72983346207418*alphavpar[2])*f[5]+96.82458365518544*f[3]*alphay[5]+38.72983346207418*f[2]*alphavpar[5]); 
  out[34] += 0.001428571428571429*((271.1088342345192*alphay[25]+303.1088913245535*alphay[4])*f[39]+(542.2176684690384*alphavpar[21]+271.1088342345192*alphay[19]+606.217782649107*alphavpar[3]+303.1088913245535*(alphay[2]+alphax[1]))*f[34]+542.2176684690384*alphavpar[19]*f[33]+606.217782649107*(alphavpar[2]*f[33]+alphavpar[1]*f[32])+(542.2176684690384*(f[23]+f[20])+606.217782649107*f[1])*alphavpar[32]+303.1088913245536*alphay[8]*f[27]+(606.2177826491072*alphavpar[6]+303.1088913245536*(alphay[5]+alphax[0]))*f[24]+(193.6491673103708*alphay[13]+271.1088342345192*alphay[11]+606.2177826491072*alphavpar[7]+303.1088913245536*(alphax[5]+alphay[0]))*f[23]+(193.6491673103708*f[13]+271.1088342345192*f[11]+303.1088913245536*f[0])*alphay[23]+606.2177826491072*alphavpar[5]*f[22]+(242.4871130596428*alphay[6]+606.2177826491072*alphavpar[5])*f[21]+242.4871130596428*f[6]*alphay[21]+606.2177826491072*(f[5]*alphavpar[21]+alphavpar[7]*f[20]+alphavpar[6]*f[19]+f[6]*alphavpar[19])+(606.217782649107*alphavpar[11]+677.7720855862981*alphavpar[0])*f[15]+(606.217782649107*(f[13]+f[12]+f[11])+677.7720855862981*f[0])*alphavpar[15]+303.1088913245535*((alphax[2]+alphay[1])*f[13]+f[1]*alphay[13])+677.7720855862981*(alphavpar[1]*f[7]+f[1]*alphavpar[7])+(271.1088342345192*alphay[3]+677.7720855862981*alphavpar[2])*f[6]+271.1088342345192*f[3]*alphay[6]+677.7720855862981*(f[2]*alphavpar[6]+alphavpar[3]*f[5]+f[3]*alphavpar[5])); 
  out[35] += 0.001428571428571429*(271.1088342345192*alphay[23]*f[39]+(193.6491673103708*alphay[21]+303.1088913245535*alphay[3])*f[37]+(271.1088342345192*alphay[20]+542.2176684690384*alphax[19]+606.217782649107*alphax[2])*f[36]+(193.6491673103708*alphay[19]+303.1088913245535*alphay[2]+606.217782649107*alphax[1])*f[35]+242.4871130596428*alphay[8]*f[28]+606.2177826491072*alphax[5]*f[26]+(193.6491673103708*alphay[11]+606.2177826491072*alphax[5]+303.1088913245536*alphay[0])*f[25]+(271.1088342345192*f[14]+193.6491673103708*f[11])*alphay[25]+303.1088913245536*(f[0]*alphay[25]+f[10]*alphay[21]+f[9]*alphay[19])+606.2177826491072*f[8]*alphax[19]+271.1088342345192*alphay[6]*f[17]+(606.217782649107*alphax[11]+271.1088342345192*alphay[5]+677.7720855862981*alphax[0])*f[16]+303.1088913245535*(alphay[4]*f[11]+f[4]*alphay[11])+677.7720855862981*(alphax[1]*f[9]+alphax[2]*f[8])+271.1088342345192*(alphay[1]*f[8]+f[1]*alphay[8])+677.7720855862981*f[4]*alphax[5]); 
  out[36] += 0.05*(19.36491673103708*alphay[13]*f[46]+17.32050807568877*alphay[6]*f[44]+(15.49193338482967*alphay[25]+17.32050807568877*alphay[4])*f[41]+19.36491673103708*alphay[23]*f[40]+(15.49193338482967*alphay[19]+17.32050807568877*alphay[2]+8.660254037844386*alphax[1])*f[36]+(15.49193338482967*alphay[20]+7.745966692414834*alphax[19]+17.32050807568877*alphay[1])*f[35]+(17.32050807568877*alphay[21]+19.36491673103709*alphay[3])*f[31]+17.32050807568877*alphay[8]*f[29]+(17.32050807568877*alphay[5]+8.660254037844387*alphax[0])*f[26]+17.32050807568877*(alphay[5]*f[25]+f[5]*alphay[25]+f[9]*alphay[20]+alphay[8]*f[19]+f[8]*alphay[19])+19.36491673103709*alphay[6]*f[18]+(17.32050807568877*(alphay[12]+alphay[11])+7.745966692414834*alphax[5]+19.36491673103709*alphay[0])*f[16]+7.745966692414834*alphax[2]*f[9]+19.36491673103709*(alphay[1]*f[9]+alphay[2]*f[8]+f[2]*alphay[8]+alphay[4]*f[5]+f[4]*alphay[5])); 
  out[37] += 0.007142857142857143*((38.72983346207417*alphavpar[32]+60.62177826491071*alphavpar[7]+121.2435565298214*alphax[5])*f[44]+(38.72983346207417*alphavpar[21]+60.6217782649107*alphavpar[3]+121.2435565298214*alphax[1])*f[37]+38.72983346207417*alphavpar[19]*f[35]+60.6217782649107*(alphavpar[2]*f[35]+f[18]*alphavpar[32])+(121.2435565298214*alphax[19]+54.22176684690384*alphavpar[15]+135.5544171172596*alphax[2])*f[31]+38.72983346207417*alphavpar[11]*f[25]+60.62177826491071*(alphavpar[0]*f[25]+f[10]*alphavpar[21]+f[9]*alphavpar[19])+135.5544171172596*alphax[5]*f[18]+(121.2435565298214*alphax[11]+54.22176684690384*alphavpar[6]+135.5544171172596*alphax[0])*f[17]+54.22176684690384*alphavpar[5]*f[16]+60.6217782649107*f[4]*alphavpar[11]+135.5544171172596*alphax[1]*f[10]+54.22176684690384*alphavpar[1]*f[8]); 
  out[38] += 0.05*(17.32050807568877*(alphay[8]*f[47]+alphay[6]*f[46])+(8.660254037844387*alphavpar[6]+17.32050807568877*alphay[5])*f[45]+(7.745966692414834*alphavpar[32]+19.36491673103708*alphay[11])*f[44]+17.32050807568877*(alphay[4]*f[43]+alphay[3]*f[40])+(8.660254037844386*alphavpar[3]+17.32050807568877*alphay[2])*f[38]+19.36491673103708*alphay[19]*f[37]+8.660254037844386*alphavpar[1]*f[36]+(19.36491673103708*alphay[21]+7.745966692414834*alphavpar[19])*f[35]+19.36491673103708*alphay[25]*f[32]+(17.32050807568877*(alphay[23]+alphay[20])+7.745966692414834*alphavpar[15]+19.36491673103709*alphay[1])*f[31]+8.660254037844387*alphavpar[0]*f[26]+(17.32050807568877*(alphay[13]+alphay[12])+7.745966692414834*alphavpar[7])*f[18]+19.36491673103709*(alphay[0]*f[18]+alphay[5]*f[17])+(19.36491673103709*alphay[6]+7.745966692414834*alphavpar[5])*f[16]+19.36491673103709*(alphay[8]*f[15]+alphay[2]*f[10])+(19.36491673103709*alphay[3]+7.745966692414834*alphavpar[2])*f[9]+19.36491673103709*alphay[4]*f[7]); 
  out[39] += 0.05*((15.49193338482967*alphavpar[32]+17.32050807568877*alphavpar[7]+8.660254037844387*alphax[5])*f[46]+17.32050807568877*alphavpar[5]*f[44]+(17.32050807568877*alphavpar[15]+8.660254037844386*alphax[2])*f[40]+(15.49193338482967*alphavpar[21]+17.32050807568877*alphavpar[3]+8.660254037844386*alphax[1])*f[39]+17.32050807568877*(alphavpar[1]*f[37]+alphavpar[15]*f[35]+f[16]*alphavpar[32])+(17.32050807568877*alphavpar[19]+19.36491673103709*alphavpar[2])*f[31]+(17.32050807568877*alphavpar[6]+8.660254037844387*alphax[0])*f[27]+17.32050807568877*(alphavpar[6]*f[25]+f[8]*alphavpar[21])+19.36491673103709*alphavpar[5]*f[18]+17.32050807568877*alphavpar[11]*f[17]+19.36491673103709*(alphavpar[0]*f[17]+alphavpar[7]*f[16]+f[9]*alphavpar[15]+alphavpar[1]*f[10]+alphavpar[3]*f[8]+f[4]*alphavpar[6])); 
  out[40] += 0.007142857142857143*((121.2435565298214*alphavpar[6]+60.62177826491071*alphay[5])*f[46]+121.2435565298214*alphavpar[5]*f[45]+135.5544171172596*alphavpar[11]*f[44]+(121.2435565298214*alphavpar[3]+60.6217782649107*alphay[2])*f[40]+(38.72983346207417*alphay[23]+121.2435565298214*alphavpar[15]+60.6217782649107*alphay[1])*f[39]+121.2435565298214*alphavpar[2]*f[38]+(54.22176684690384*alphay[21]+135.5544171172596*alphavpar[19])*f[37]+121.2435565298214*alphavpar[15]*f[36]+135.5544171172596*(alphavpar[21]*f[35]+f[25]*alphavpar[32])+135.5544171172596*alphavpar[1]*f[31]+(38.72983346207417*alphay[13]+121.2435565298214*alphavpar[7]+60.62177826491071*alphay[0])*f[27]+121.2435565298214*alphavpar[7]*f[26]+60.62177826491071*(alphay[8]*f[23]+f[8]*alphay[23])+135.5544171172596*alphavpar[0]*f[18]+54.22176684690384*alphay[6]*f[17]+135.5544171172596*(alphavpar[5]*f[17]+alphavpar[6]*f[16]+f[8]*alphavpar[15])+60.6217782649107*(alphay[4]*f[13]+f[4]*alphay[13])+54.22176684690384*alphay[3]*f[10]+135.5544171172596*(alphavpar[2]*f[10]+alphavpar[3]*f[9]+f[4]*alphavpar[7])); 
  out[41] += 0.01*((38.72983346207417*alphay[21]+43.30127018922193*alphay[3])*f[42]+(38.72983346207417*alphay[19]+43.30127018922193*(alphay[2]+alphax[1]))*f[41]+43.30127018922195*(alphay[6]*f[30]+(alphay[5]+alphax[0])*f[29])+(38.72983346207417*alphay[11]+43.30127018922195*(alphax[5]+alphay[0]))*f[28]+34.64101615137755*(alphay[8]*f[25]+f[8]*alphay[25])+43.30127018922193*(alphax[2]+alphay[1])*f[14]+38.72983346207418*(alphay[4]*f[8]+f[4]*alphay[8])); 
  out[42] += 0.05*((7.745966692414834*alphavpar[32]+8.660254037844387*(alphavpar[7]+alphax[5]))*f[47]+8.660254037844386*(alphavpar[15]+alphax[2])*f[43]+(7.745966692414834*alphavpar[21]+8.660254037844386*(alphavpar[3]+alphax[1]))*f[42]+(7.745966692414834*alphavpar[19]+8.660254037844386*alphavpar[2])*f[41]+8.660254037844387*((alphavpar[6]+alphax[0])*f[30]+alphavpar[5]*f[29])+(7.745966692414834*alphavpar[11]+8.660254037844387*alphavpar[0])*f[28]+8.660254037844386*alphavpar[1]*f[14]); 
  out[43] += 0.05*(8.660254037844387*(alphavpar[6]+alphay[5])*f[47]+8.660254037844386*(alphavpar[3]+alphay[2])*f[43]+7.745966692414834*alphay[23]*f[42]+8.660254037844386*((alphavpar[15]+alphay[1])*f[42]+alphavpar[1]*f[41])+7.745966692414834*(alphay[25]*f[37]+alphay[13]*f[30])+8.660254037844387*((alphavpar[7]+alphay[0])*f[30]+alphavpar[0]*f[29]+(alphay[6]+alphavpar[5])*f[28])+7.745966692414834*alphay[8]*f[17]+8.660254037844386*(alphay[3]+alphavpar[2])*f[14]+7.745966692414834*alphay[4]*f[10]); 
  out[44] += 0.001428571428571429*((271.1088342345192*alphay[20]+542.2176684690384*alphax[19]+242.4871130596428*alphavpar[15]+606.217782649107*alphax[2])*f[45]+(193.6491673103708*(alphavpar[21]+alphay[19])+303.1088913245535*(alphavpar[3]+alphay[2])+606.217782649107*alphax[1])*f[44]+242.4871130596428*(alphay[8]*f[42]+alphay[6]*f[39])+(271.1088342345192*alphavpar[32]+606.2177826491072*alphax[5])*f[38]+(193.6491673103708*alphavpar[32]+271.1088342345192*alphay[13]+193.6491673103708*alphay[11]+303.1088913245536*alphavpar[7]+606.2177826491072*alphax[5]+303.1088913245536*alphay[0])*f[37]+242.4871130596428*alphavpar[5]*f[36]+193.6491673103708*alphavpar[11]*f[35]+303.1088913245536*(alphavpar[0]*f[35]+f[10]*alphavpar[32])+(606.2177826491072*alphax[11]+271.1088342345192*(alphavpar[6]+alphay[5])+677.772085586298*alphax[0])*f[31]+271.1088342345192*(alphay[25]*f[30]+alphay[21]*f[27]+alphavpar[19]*f[26])+(193.6491673103708*(alphay[21]+alphavpar[19])+303.1088913245535*(alphay[3]+alphavpar[2]))*f[25]+(193.6491673103708*f[21]+303.1088913245535*f[3])*alphay[25]+242.4871130596428*f[17]*alphay[23]+303.1088913245535*(alphay[4]*f[21]+f[4]*alphay[21]+f[18]*(alphavpar[21]+alphay[19]))+606.217782649107*f[17]*alphax[19]+303.1088913245535*f[4]*alphavpar[19]+677.772085586298*alphax[1]*f[18]+(271.1088342345192*alphavpar[15]+677.772085586298*alphax[2])*f[17]+271.1088342345192*(alphay[1]*f[17]+alphavpar[1]*f[16])+303.1088913245536*(f[10]*alphay[11]+f[9]*alphavpar[11])+677.772085586298*alphax[5]*f[10]+271.1088342345192*((alphay[6]+alphavpar[5])*f[8]+f[6]*alphay[8])); 
  out[45] += 0.01*((77.45966692414834*alphay[25]+86.60254037844386*alphay[4])*f[47]+(77.45966692414834*alphay[21]+86.60254037844386*alphay[3])*f[46]+(38.72983346207417*alphavpar[21]+77.45966692414834*alphay[19]+43.30127018922193*alphavpar[3]+86.60254037844386*alphay[2]+43.30127018922193*alphax[1])*f[45]+(77.45966692414834*(alphay[23]+alphay[20])+38.72983346207417*alphax[19]+34.64101615137754*alphavpar[15]+86.60254037844386*alphay[1])*f[44]+86.60254037844389*(alphay[8]*f[43]+alphay[6]*f[40])+(43.30127018922195*alphavpar[6]+86.60254037844389*alphay[5]+43.30127018922195*alphax[0])*f[38]+86.60254037844389*alphay[5]*f[37]+(38.72983346207417*alphavpar[11]+43.30127018922195*alphavpar[0])*f[36]+(86.60254037844389*alphay[6]+34.64101615137755*alphavpar[5])*f[35]+86.60254037844389*alphay[8]*f[32]+f[31]*(34.64101615137755*alphavpar[32]+86.60254037844389*(alphay[13]+alphay[12]+alphay[11])+38.72983346207417*(alphavpar[7]+alphax[5])+96.82458365518542*alphay[0])+43.30127018922193*alphavpar[1]*f[26]+86.60254037844386*(f[15]*alphay[25]+f[18]*alphay[23]+f[16]*alphay[21]+f[18]*alphay[20]+f[17]*alphay[19])+34.64101615137754*f[16]*alphavpar[19]+38.72983346207417*(alphavpar[15]+alphax[2])*f[18]+96.82458365518542*(alphay[1]*f[18]+alphay[2]*f[17])+(96.82458365518542*alphay[3]+38.72983346207417*alphavpar[2])*f[16]+96.82458365518542*(alphay[4]*f[15]+alphay[5]*f[10])+(96.82458365518542*alphay[6]+38.72983346207417*alphavpar[5])*f[9]+96.82458365518542*f[7]*alphay[8]); 
  out[46] += 0.001428571428571429*((542.2176684690384*alphavpar[21]+271.1088342345192*alphay[19]+606.217782649107*alphavpar[3]+303.1088913245535*(alphay[2]+alphax[1]))*f[46]+542.2176684690384*alphavpar[19]*f[45]+606.217782649107*(alphavpar[2]*f[45]+alphavpar[1]*f[44])+(606.2177826491072*alphavpar[6]+303.1088913245536*(alphay[5]+alphax[0]))*f[40]+(542.2176684690384*alphavpar[32]+193.6491673103708*alphay[13]+271.1088342345192*alphay[11]+606.2177826491072*alphavpar[7]+303.1088913245536*(alphax[5]+alphay[0]))*f[39]+606.2177826491072*alphavpar[5]*f[38]+(242.4871130596428*alphay[6]+606.2177826491072*alphavpar[5])*f[37]+542.2176684690384*alphavpar[32]*f[36]+606.2177826491072*(alphavpar[7]*f[36]+alphavpar[6]*f[35]+f[8]*alphavpar[32])+(606.2177826491072*alphavpar[11]+677.772085586298*alphavpar[0])*f[31]+(193.6491673103708*alphay[23]+606.217782649107*alphavpar[15]+303.1088913245535*(alphax[2]+alphay[1]))*f[27]+606.217782649107*alphavpar[15]*f[26]+(271.1088342345192*alphay[23]+606.217782649107*alphavpar[15])*f[25]+271.1088342345192*f[23]*alphay[25]+303.1088913245535*(alphay[4]*f[23]+f[4]*alphay[23])+242.4871130596428*f[17]*alphay[21]+606.217782649107*(f[16]*alphavpar[21]+f[17]*alphavpar[19])+677.772085586298*alphavpar[1]*f[18]+271.1088342345192*alphay[3]*f[17]+677.772085586298*(alphavpar[2]*f[17]+alphavpar[3]*f[16]+f[4]*alphavpar[15])+303.1088913245536*(alphay[8]*f[13]+f[8]*alphay[13])+271.1088342345192*alphay[6]*f[10]+677.772085586298*(alphavpar[5]*f[10]+alphavpar[6]*f[9]+alphavpar[7]*f[8])); 
  out[47] += 0.01*((38.72983346207417*(alphavpar[21]+alphay[19])+43.30127018922193*(alphavpar[3]+alphay[2]+alphax[1]))*f[47]+43.30127018922195*(alphavpar[6]+alphay[5]+alphax[0])*f[43]+(38.72983346207417*(alphavpar[32]+alphay[13]+alphay[11])+43.30127018922195*(alphavpar[7]+alphax[5]+alphay[0]))*f[42]+(38.72983346207417*alphavpar[11]+43.30127018922195*alphavpar[0])*f[41]+34.64101615137755*alphay[8]*f[37]+38.72983346207417*alphay[23]*f[30]+43.30127018922193*((alphavpar[15]+alphax[2]+alphay[1])*f[30]+alphavpar[1]*f[29])+(38.72983346207417*(alphay[21]+alphavpar[19])+43.30127018922193*(alphay[3]+alphavpar[2]))*f[28]+f[17]*(34.64101615137754*alphay[25]+38.72983346207417*alphay[4])+43.30127018922195*(alphay[6]+alphavpar[5])*f[14]+38.72983346207417*alphay[8]*f[10]); 
  return cflFreq; 
} 
