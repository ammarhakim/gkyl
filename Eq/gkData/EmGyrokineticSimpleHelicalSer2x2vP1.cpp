#include <GyrokineticModDecl.h> 
double EmGyrokineticSimpleHelicalVol2x2vSerP1_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
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

  double hamil[16]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*phi[1]*q_; 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 

  double BstarXdBmag[16]; 
  BstarXdBmag[0] = (2.0*BdriftX[0]*m_*wvpar+q_*(1.732050807568877*bmagInv[0]*Apar[2]*rdy2+Apar[0]*BdriftX[0]))/q_; 
  BstarXdBmag[1] = 1.732050807568877*bmagInv[0]*Apar[3]*rdy2+BdriftX[0]*Apar[1]; 
  BstarXdBmag[2] = BdriftX[0]*Apar[2]; 
  BstarXdBmag[3] = (1.154700538379252*BdriftX[0]*m_)/(q_*rdvpar2); 
  BstarXdBmag[5] = BdriftX[0]*Apar[3]; 

  double BstarYdBmag[16]; 
  BstarYdBmag[0] = (2.0*BdriftY[0]*m_*wvpar+q_*(Apar[0]*BdriftY[0]-1.732050807568877*bmagInv[0]*Apar[1]*rdx2))/q_; 
  BstarYdBmag[1] = BdriftY[0]*Apar[1]; 
  BstarYdBmag[2] = -0.5773502691896258*(3.0*bmagInv[0]*Apar[3]*rdx2-1.732050807568877*BdriftY[0]*Apar[2]); 
  BstarYdBmag[3] = (1.154700538379252*BdriftY[0]*m_)/(q_*rdvpar2); 
  BstarYdBmag[5] = BdriftY[0]*Apar[3]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[0]*hamil[3]*rdvpar2)/m_-(0.5*bmagInv[0]*hamil[2]*rdy2)/q_); 
  alphax[1] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[1]*hamil[3]*rdvpar2)/m_-(0.5*bmagInv[0]*hamil[5]*rdy2)/q_); 
  alphax[2] = (0.4330127018922193*BstarXdBmag[2]*hamil[3]*rdvpar2*rdx2)/m_; 
  alphax[3] = (0.4330127018922193*BstarXdBmag[3]*hamil[3]*rdvpar2*rdx2)/m_; 
  alphax[5] = (0.4330127018922193*hamil[3]*BstarXdBmag[5]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.4330127018922193*alphax[5]-0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[5])-0.25*alphax[3]+0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[5])+0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[5]-0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[5])-0.25*alphax[3]+0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[5])+0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.4330127018922193*alphax[5])-0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[5])+0.25*alphax[3]-0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[5]+0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[5])-0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[5])+0.25*alphax[3]-0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[5]+0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[16]; 
  alphay[0] = 1.732050807568877*((0.5*bmagInv[0]*hamil[1]*rdx2)/q_+(0.25*BstarYdBmag[0]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[1] = (0.4330127018922193*BstarYdBmag[1]*hamil[3]*rdvpar2*rdy2)/m_; 
  alphay[2] = 1.732050807568877*((0.5*bmagInv[0]*hamil[5]*rdx2)/q_+(0.25*BstarYdBmag[2]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[3] = (0.4330127018922193*BstarYdBmag[3]*hamil[3]*rdvpar2*rdy2)/m_; 
  alphay[5] = (0.4330127018922193*hamil[3]*BstarYdBmag[5]*rdvpar2*rdy2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphay[2]-1.0*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphay[2]+alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.4330127018922193*alphay[5]-0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphay[5])-0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphay[5]+0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphay[5])+0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphay[5]-0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphay[5])-0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphay[5]+0.25*alphay[3]-0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphay[5])+0.25*alphay[3]-0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.4330127018922193*alphay[5])-0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[5]-0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphay[5])+0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[5]+0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphay[5])-0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[5]-0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphay[5])+0.25*alphay[3]+0.4330127018922193*alphay[2]-0.25*alphay[1]+0.25*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphay[5]+0.25*alphay[3]+0.4330127018922193*alphay[2]+0.25*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[16]; 
  alphavpar[0] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[1]*hamil[5]+BstarYdBmag[0]*hamil[2])*rdy2+(BstarXdBmag[2]*hamil[5]+BstarXdBmag[0]*hamil[1])*rdx2))/m_; 
  alphavpar[1] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[0]*hamil[5]+BstarYdBmag[1]*hamil[2])*rdy2+(BstarXdBmag[5]*hamil[5]+BstarXdBmag[1]*hamil[1])*rdx2))/m_; 
  alphavpar[2] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[5]*hamil[5]+BstarYdBmag[2]*hamil[2])*rdy2+(BstarXdBmag[0]*hamil[5]+hamil[1]*BstarXdBmag[2])*rdx2))/m_; 
  alphavpar[3] = -(0.4330127018922193*rdvpar2*(hamil[2]*BstarYdBmag[3]*rdy2+hamil[1]*BstarXdBmag[3]*rdx2))/m_; 
  alphavpar[5] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[2]*hamil[5]+hamil[2]*BstarYdBmag[5])*rdy2+(BstarXdBmag[1]*hamil[5]+hamil[1]*BstarXdBmag[5])*rdx2))/m_; 
  alphavpar[6] = -(0.4330127018922193*BstarYdBmag[3]*hamil[5]*rdvpar2*rdy2)/m_; 
  alphavpar[7] = -(0.4330127018922193*BstarXdBmag[3]*hamil[5]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphavpar[3]-1.0*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphavpar[3]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphavpar[7]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphavpar[7])+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(alphavpar[7]+alphavpar[6]))+0.25*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphavpar[7]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphavpar[7])+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(alphavpar[7]+alphavpar[6]))+0.25*alphavpar[5]-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.4330127018922193*(alphavpar[7]+alphavpar[6]))+0.25*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphavpar[7])+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphavpar[7]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(alphavpar[7]+alphavpar[6]))+0.25*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphavpar[7])+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphavpar[7]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*alphavpar[5]+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[5]*f[5]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[5]*f[5]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[5] += 0.4330127018922193*(alphax[3]*f[7]+alphay[3]*f[6]+(alphay[2]+alphax[1])*f[5]+f[2]*alphay[5]+f[1]*alphax[5]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[6] += 0.4330127018922193*((alphavpar[7]+alphax[5])*f[11]+alphax[2]*f[7]+(alphavpar[3]+alphax[1])*f[6]+f[3]*alphavpar[6]+alphavpar[2]*f[5]+f[2]*alphavpar[5]+alphax[0]*f[3]+f[0]*alphax[3]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[7] += 0.4330127018922193*((alphavpar[6]+alphay[5])*f[11]+(alphavpar[3]+alphay[2])*f[7]+f[3]*alphavpar[7]+alphay[1]*f[6]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphay[0]*f[3]+f[0]*alphay[3]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[8] += 0.4330127018922193*(alphax[5]*f[12]+alphax[3]*f[10]+alphax[2]*f[9]+alphax[1]*f[8]+alphax[0]*f[4]); 
  out[9] += 0.4330127018922193*(alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+alphay[0]*f[4]); 
  out[10] += 0.4330127018922193*(alphavpar[7]*f[14]+alphavpar[6]*f[13]+alphavpar[5]*f[12]+alphavpar[3]*f[10]+alphavpar[2]*f[9]+alphavpar[1]*f[8]+alphavpar[0]*f[4]); 
  out[11] += 0.4330127018922193*((alphavpar[3]+alphay[2]+alphax[1])*f[11]+(alphavpar[6]+alphay[5]+alphax[0])*f[7]+f[6]*(alphavpar[7]+alphax[5]+alphay[0])+alphavpar[0]*f[5]+f[0]*alphavpar[5]+(alphax[2]+alphay[1])*f[3]+f[1]*alphay[3]+f[2]*(alphax[3]+alphavpar[1])+f[1]*alphavpar[2]); 
  out[12] += 0.4330127018922193*(alphax[3]*f[14]+alphay[3]*f[13]+(alphay[2]+alphax[1])*f[12]+(alphay[5]+alphax[0])*f[9]+(alphax[5]+alphay[0])*f[8]+(alphax[2]+alphay[1])*f[4]); 
  out[13] += 0.4330127018922193*((alphavpar[7]+alphax[5])*f[15]+alphax[2]*f[14]+(alphavpar[3]+alphax[1])*f[13]+alphavpar[2]*f[12]+(alphavpar[6]+alphax[0])*f[10]+alphavpar[5]*f[9]+alphavpar[0]*f[8]+(alphax[3]+alphavpar[1])*f[4]); 
  out[14] += 0.4330127018922193*((alphavpar[6]+alphay[5])*f[15]+(alphavpar[3]+alphay[2])*f[14]+alphay[1]*f[13]+alphavpar[1]*f[12]+(alphavpar[7]+alphay[0])*f[10]+alphavpar[0]*f[9]+alphavpar[5]*f[8]+(alphay[3]+alphavpar[2])*f[4]); 
  out[15] += 0.4330127018922193*((alphavpar[3]+alphay[2]+alphax[1])*f[15]+(alphavpar[6]+alphay[5]+alphax[0])*f[14]+(alphavpar[7]+alphax[5]+alphay[0])*f[13]+alphavpar[0]*f[12]+(alphax[2]+alphay[1])*f[10]+(alphax[3]+alphavpar[1])*f[9]+(alphay[3]+alphavpar[2])*f[8]+f[4]*alphavpar[5]); 
  return cflFreq; 
} 
double EmGyrokineticSimpleHelicalVol2x2vSerP1_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
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

  double hamil[16]; 
  hamil[0] = (0.6666666666666666*(3.0*rdvpar2Sq*(m_*wvparSq+bmag[0]*wmu+phi[0]*q_)+m_))/rdvpar2Sq; 
  hamil[1] = 2.0*(bmag[1]*wmu+phi[1]*q_); 
  hamil[2] = 2.0*phi[2]*q_; 
  hamil[3] = (2.309401076758503*m_*wvpar)/rdvpar2; 
  hamil[4] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil[5] = 2.0*phi[3]*q_; 
  hamil[8] = (1.154700538379252*bmag[1])/rdmu2; 

  double BstarXdBmag[16]; 
  BstarXdBmag[0] = (2.0*BdriftX[0]*m_*wvpar+q_*(1.732050807568877*(bmagInv[1]*Apar[3]+bmagInv[0]*Apar[2])*rdy2+Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0]))/q_; 
  BstarXdBmag[1] = (2.0*BdriftX[1]*m_*wvpar+q_*(1.732050807568877*(bmagInv[0]*Apar[3]+bmagInv[1]*Apar[2])*rdy2+Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1]))/q_; 
  BstarXdBmag[2] = BdriftX[1]*Apar[3]+BdriftX[0]*Apar[2]; 
  BstarXdBmag[3] = (1.154700538379252*BdriftX[0]*m_)/(q_*rdvpar2); 
  BstarXdBmag[5] = BdriftX[0]*Apar[3]+BdriftX[1]*Apar[2]; 
  BstarXdBmag[6] = (1.154700538379252*BdriftX[1]*m_)/(q_*rdvpar2); 

  double BstarYdBmag[16]; 
  BstarYdBmag[0] = (2.0*BdriftY[0]*m_*wvpar+q_*(Apar[1]*(BdriftY[1]-1.732050807568877*bmagInv[0]*rdx2)+Apar[0]*BdriftY[0]))/q_; 
  BstarYdBmag[1] = (2.0*BdriftY[1]*m_*wvpar+q_*((-1.732050807568877*Apar[1]*bmagInv[1]*rdx2)+Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1]))/q_; 
  BstarYdBmag[2] = -0.5773502691896258*(3.0*bmagInv[0]*Apar[3]*rdx2-1.732050807568877*(BdriftY[1]*Apar[3]+BdriftY[0]*Apar[2])); 
  BstarYdBmag[3] = (1.154700538379252*BdriftY[0]*m_)/(q_*rdvpar2); 
  BstarYdBmag[5] = -1.0*(1.732050807568877*bmagInv[1]*Apar[3]*rdx2-1.0*(BdriftY[0]*Apar[3]+BdriftY[1]*Apar[2])); 
  BstarYdBmag[6] = (1.154700538379252*BdriftY[1]*m_)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[16]; 
  alphax[0] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[0]*hamil[3]*rdvpar2)/m_-(0.5*(bmagInv[1]*hamil[5]+bmagInv[0]*hamil[2])*rdy2)/q_); 
  alphax[1] = 1.732050807568877*rdx2*((0.25*BstarXdBmag[1]*hamil[3]*rdvpar2)/m_-(0.5*(bmagInv[0]*hamil[5]+bmagInv[1]*hamil[2])*rdy2)/q_); 
  alphax[2] = (0.4330127018922193*BstarXdBmag[2]*hamil[3]*rdvpar2*rdx2)/m_; 
  alphax[3] = (0.4330127018922193*BstarXdBmag[3]*hamil[3]*rdvpar2*rdx2)/m_; 
  alphax[5] = (0.4330127018922193*hamil[3]*BstarXdBmag[5]*rdvpar2*rdx2)/m_; 
  alphax[6] = (0.4330127018922193*hamil[3]*BstarXdBmag[6]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphax[1]-1.0*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphax[1]+alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*(0.4330127018922193*(alphax[6]+alphax[5])-0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[6])+0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(alphax[6]+alphax[5]))+0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*(alphax[6]+alphax[5])-0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphax[6])+0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*(alphax[6]+alphax[5]))+0.25*(alphax[3]+alphax[2])-0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*((-0.4330127018922193*(alphax[6]+alphax[5]))-0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[6])+0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphax[6]+alphax[5])+0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*(alphax[6]+alphax[5]))-0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphax[6])+0.4330127018922193*alphax[5]-0.25*alphax[3]+0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphax[6]-0.4330127018922193*alphax[5]+0.25*alphax[3]-0.25*alphax[2]+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*(alphax[6]+alphax[5])+0.25*(alphax[3]+alphax[2])+0.4330127018922193*alphax[1]+0.25*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[16]; 
  alphay[0] = 1.732050807568877*((0.5*bmagInv[0]*hamil[1]*rdx2)/q_+(0.25*BstarYdBmag[0]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[1] = 1.732050807568877*((0.5*bmagInv[1]*hamil[1]*rdx2)/q_+(0.25*BstarYdBmag[1]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[2] = 1.732050807568877*((0.5*bmagInv[0]*hamil[5]*rdx2)/q_+(0.25*BstarYdBmag[2]*hamil[3]*rdvpar2)/m_)*rdy2; 
  alphay[3] = (0.4330127018922193*BstarYdBmag[3]*hamil[3]*rdvpar2*rdy2)/m_; 
  alphay[4] = (0.8660254037844386*bmagInv[0]*hamil[8]*rdx2*rdy2)/q_; 
  alphay[5] = 1.732050807568877*((0.5*bmagInv[1]*hamil[5]*rdx2)/q_+(0.25*hamil[3]*BstarYdBmag[5]*rdvpar2)/m_)*rdy2; 
  alphay[6] = (0.4330127018922193*hamil[3]*BstarYdBmag[6]*rdvpar2*rdy2)/m_; 
  alphay[8] = (0.8660254037844386*bmagInv[1]*hamil[8]*rdx2*rdy2)/q_; 
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
  alphavpar[0] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[1]*hamil[5]+BstarYdBmag[0]*hamil[2])*rdy2+(BstarXdBmag[2]*hamil[5]+BstarXdBmag[0]*hamil[1])*rdx2))/m_; 
  alphavpar[1] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[0]*hamil[5]+BstarYdBmag[1]*hamil[2])*rdy2+(BstarXdBmag[5]*hamil[5]+BstarXdBmag[1]*hamil[1])*rdx2))/m_; 
  alphavpar[2] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[5]*hamil[5]+BstarYdBmag[2]*hamil[2])*rdy2+(BstarXdBmag[0]*hamil[5]+hamil[1]*BstarXdBmag[2])*rdx2))/m_; 
  alphavpar[3] = -(0.4330127018922193*rdvpar2*((hamil[5]*BstarYdBmag[6]+hamil[2]*BstarYdBmag[3])*rdy2+hamil[1]*BstarXdBmag[3]*rdx2))/m_; 
  alphavpar[4] = -(0.4330127018922193*BstarXdBmag[0]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphavpar[5] = -(0.4330127018922193*rdvpar2*((BstarYdBmag[2]*hamil[5]+hamil[2]*BstarYdBmag[5])*rdy2+(BstarXdBmag[1]*hamil[5]+hamil[1]*BstarXdBmag[5])*rdx2))/m_; 
  alphavpar[6] = -(0.4330127018922193*rdvpar2*((hamil[2]*BstarYdBmag[6]+BstarYdBmag[3]*hamil[5])*rdy2+hamil[1]*BstarXdBmag[6]*rdx2))/m_; 
  alphavpar[7] = -(0.4330127018922193*BstarXdBmag[3]*hamil[5]*rdvpar2*rdx2)/m_; 
  alphavpar[8] = -(0.4330127018922193*BstarXdBmag[1]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphavpar[9] = -(0.4330127018922193*BstarXdBmag[2]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphavpar[10] = -(0.4330127018922193*BstarXdBmag[3]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphavpar[11] = -(0.4330127018922193*hamil[5]*BstarXdBmag[6]*rdvpar2*rdx2)/m_; 
  alphavpar[12] = -(0.4330127018922193*BstarXdBmag[5]*hamil[8]*rdvpar2*rdx2)/m_; 
  alphavpar[13] = -(0.4330127018922193*BstarXdBmag[6]*hamil[8]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.125*(1.732050807568877*alphavpar[3]-1.0*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.125*(1.732050807568877*alphavpar[3]+alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0625*((-0.4330127018922193*alphavpar[13])-0.25*alphavpar[12]-0.4330127018922193*alphavpar[11]+0.4330127018922193*alphavpar[10]+0.25*(alphavpar[9]+alphavpar[8])+0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphavpar[13]+0.25*alphavpar[12]+0.4330127018922193*(alphavpar[11]+alphavpar[10])+0.25*alphavpar[9]-0.25*alphavpar[8]+0.4330127018922193*alphavpar[7]-0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphavpar[13])+0.25*alphavpar[12]+0.4330127018922193*(alphavpar[11]+alphavpar[10])-0.25*alphavpar[9]+0.25*alphavpar[8]-0.4330127018922193*alphavpar[7]+0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphavpar[13]-0.25*alphavpar[12]-0.4330127018922193*alphavpar[11]+0.4330127018922193*alphavpar[10]-0.25*(alphavpar[9]+alphavpar[8])-0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*alphavpar[5]-0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphavpar[13]+0.25*alphavpar[12]-0.4330127018922193*(alphavpar[11]+alphavpar[10])-0.25*(alphavpar[9]+alphavpar[8])+0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphavpar[13])-0.25*alphavpar[12]+0.4330127018922193*alphavpar[11]-0.4330127018922193*alphavpar[10]-0.25*alphavpar[9]+0.25*alphavpar[8]+0.4330127018922193*alphavpar[7]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*(0.4330127018922193*alphavpar[13]-0.25*alphavpar[12]+0.4330127018922193*alphavpar[11]-0.4330127018922193*alphavpar[10]+0.25*alphavpar[9]-0.25*alphavpar[8]-0.4330127018922193*alphavpar[7]+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]-0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0625*((-0.4330127018922193*alphavpar[13])+0.25*alphavpar[12]-0.4330127018922193*(alphavpar[11]+alphavpar[10])+0.25*(alphavpar[9]+alphavpar[8])-0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*(alphavpar[5]+alphavpar[4])-0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0625*(0.4330127018922193*alphavpar[13]-0.25*alphavpar[12]+0.4330127018922193*alphavpar[11]-0.4330127018922193*alphavpar[10]+0.25*(alphavpar[9]+alphavpar[8])-0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphavpar[13])+0.25*alphavpar[12]-0.4330127018922193*(alphavpar[11]+alphavpar[10])+0.25*alphavpar[9]-0.25*alphavpar[8]-0.4330127018922193*alphavpar[7]+0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphavpar[13]+0.25*alphavpar[12]-0.4330127018922193*(alphavpar[11]+alphavpar[10])-0.25*alphavpar[9]+0.25*alphavpar[8]+0.4330127018922193*alphavpar[7]-0.4330127018922193*alphavpar[6]-0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphavpar[13])-0.25*alphavpar[12]+0.4330127018922193*alphavpar[11]-0.4330127018922193*alphavpar[10]-0.25*(alphavpar[9]+alphavpar[8])+0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*alphavpar[5]-0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphavpar[13])+0.25*alphavpar[12]+0.4330127018922193*(alphavpar[11]+alphavpar[10])-0.25*(alphavpar[9]+alphavpar[8])-0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]-0.25*(alphavpar[2]+alphavpar[1])+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphavpar[13]-0.25*alphavpar[12]-0.4330127018922193*alphavpar[11]+0.4330127018922193*alphavpar[10]-0.25*alphavpar[9]+0.25*alphavpar[8]-0.4330127018922193*alphavpar[7]+0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]-0.25*alphavpar[2]+0.25*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*((-0.4330127018922193*alphavpar[13])-0.25*alphavpar[12]-0.4330127018922193*alphavpar[11]+0.4330127018922193*alphavpar[10]+0.25*alphavpar[9]-0.25*alphavpar[8]+0.4330127018922193*alphavpar[7]-0.4330127018922193*alphavpar[6]-0.25*alphavpar[5]+0.25*alphavpar[4]+0.4330127018922193*alphavpar[3]+0.25*alphavpar[2]-0.25*alphavpar[1]+0.25*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0625*(0.4330127018922193*alphavpar[13]+0.25*alphavpar[12]+0.4330127018922193*(alphavpar[11]+alphavpar[10])+0.25*(alphavpar[9]+alphavpar[8])+0.4330127018922193*(alphavpar[7]+alphavpar[6])+0.25*(alphavpar[5]+alphavpar[4])+0.4330127018922193*alphavpar[3]+0.25*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.4330127018922193*(alphax[6]*f[6]+alphax[5]*f[5]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.4330127018922193*(alphay[8]*f[8]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.4330127018922193*(alphavpar[13]*f[13]+alphavpar[12]*f[12]+alphavpar[11]*f[11]+alphavpar[10]*f[10]+alphavpar[9]*f[9]+alphavpar[8]*f[8]+alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[5] += 0.4330127018922193*(alphax[6]*f[11]+alphay[4]*f[8]+f[4]*alphay[8]+alphax[3]*f[7]+alphay[3]*f[6]+f[3]*alphay[6]+(alphay[2]+alphax[1])*f[5]+f[2]*alphay[5]+f[1]*alphax[5]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[6] += 0.4330127018922193*(alphavpar[10]*f[13]+f[10]*alphavpar[13]+alphavpar[9]*f[12]+f[9]*alphavpar[12]+(alphavpar[7]+alphax[5])*f[11]+f[7]*alphavpar[11]+alphavpar[4]*f[8]+f[4]*alphavpar[8]+alphax[2]*f[7]+(alphavpar[3]+alphax[1])*f[6]+f[1]*alphax[6]+f[3]*alphavpar[6]+alphavpar[2]*f[5]+f[2]*alphavpar[5]+alphax[0]*f[3]+f[0]*alphax[3]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[7] += 0.4330127018922193*(alphavpar[13]*f[15]+alphavpar[10]*f[14]+alphay[8]*f[13]+alphavpar[8]*f[12]+f[8]*alphavpar[12]+(alphavpar[6]+alphay[5])*f[11]+f[6]*alphavpar[11]+alphay[4]*f[10]+alphavpar[4]*f[9]+f[4]*alphavpar[9]+(alphavpar[3]+alphay[2])*f[7]+f[3]*alphavpar[7]+alphay[1]*f[6]+f[1]*alphay[6]+alphavpar[1]*f[5]+f[1]*alphavpar[5]+alphay[0]*f[3]+f[0]*alphay[3]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[8] += 0.4330127018922193*(alphax[6]*f[13]+alphax[5]*f[12]+alphax[3]*f[10]+alphax[2]*f[9]+alphax[1]*f[8]+alphax[0]*f[4]); 
  out[9] += 0.4330127018922193*(alphay[6]*f[13]+alphay[5]*f[12]+alphay[3]*f[10]+alphay[2]*f[9]+alphay[1]*f[8]+f[1]*alphay[8]+alphay[0]*f[4]+f[0]*alphay[4]); 
  out[10] += 0.4330127018922193*(alphavpar[11]*f[15]+alphavpar[7]*f[14]+alphavpar[6]*f[13]+f[6]*alphavpar[13]+alphavpar[5]*f[12]+f[5]*alphavpar[12]+alphavpar[3]*f[10]+f[3]*alphavpar[10]+alphavpar[2]*f[9]+f[2]*alphavpar[9]+alphavpar[1]*f[8]+f[1]*alphavpar[8]+alphavpar[0]*f[4]+f[0]*alphavpar[4]); 
  out[11] += 0.4330127018922193*(alphavpar[10]*f[15]+alphavpar[13]*f[14]+alphay[4]*f[13]+alphavpar[4]*f[12]+f[4]*alphavpar[12]+(alphavpar[3]+alphay[2]+alphax[1])*f[11]+f[3]*alphavpar[11]+alphay[8]*f[10]+alphavpar[8]*f[9]+f[8]*alphavpar[9]+(alphavpar[6]+alphay[5]+alphax[0])*f[7]+f[6]*(alphavpar[7]+alphax[5]+alphay[0])+f[0]*alphay[6]+f[5]*(alphax[6]+alphavpar[0])+f[0]*alphavpar[5]+(alphax[2]+alphay[1])*f[3]+f[1]*alphay[3]+f[2]*(alphax[3]+alphavpar[1])+f[1]*alphavpar[2]); 
  out[12] += 0.4330127018922193*(alphax[6]*f[15]+alphax[3]*f[14]+alphay[3]*f[13]+(alphay[2]+alphax[1])*f[12]+alphay[6]*f[10]+(alphay[5]+alphax[0])*f[9]+(alphax[5]+alphay[0])*f[8]+f[0]*alphay[8]+(alphax[2]+alphay[1])*f[4]+f[1]*alphay[4]); 
  out[13] += 0.4330127018922193*((alphavpar[7]+alphax[5])*f[15]+(alphavpar[11]+alphax[2])*f[14]+(alphavpar[3]+alphax[1])*f[13]+f[3]*alphavpar[13]+alphavpar[2]*f[12]+f[2]*alphavpar[12]+(alphavpar[6]+alphax[0])*f[10]+f[6]*alphavpar[10]+alphavpar[5]*f[9]+f[5]*alphavpar[9]+(alphax[6]+alphavpar[0])*f[8]+f[0]*alphavpar[8]+(alphax[3]+alphavpar[1])*f[4]+f[1]*alphavpar[4]); 
  out[14] += 0.4330127018922193*((alphavpar[6]+alphay[5])*f[15]+(alphavpar[3]+alphay[2])*f[14]+(alphavpar[11]+alphay[1])*f[13]+f[11]*alphavpar[13]+alphavpar[1]*f[12]+f[1]*alphavpar[12]+(alphavpar[7]+alphay[0])*f[10]+f[7]*alphavpar[10]+alphavpar[0]*f[9]+f[0]*alphavpar[9]+(alphay[6]+alphavpar[5])*f[8]+f[6]*alphay[8]+f[5]*alphavpar[8]+(alphay[3]+alphavpar[2])*f[4]+f[3]*alphay[4]+f[2]*alphavpar[4]); 
  out[15] += 0.4330127018922193*((alphavpar[3]+alphay[2]+alphax[1])*f[15]+(alphavpar[6]+alphay[5]+alphax[0])*f[14]+(alphavpar[7]+alphax[5]+alphay[0])*f[13]+f[7]*alphavpar[13]+(alphax[6]+alphavpar[0])*f[12]+f[0]*alphavpar[12]+alphavpar[10]*f[11]+f[10]*(alphavpar[11]+alphax[2]+alphay[1])+(alphax[3]+alphavpar[1])*f[9]+f[1]*alphavpar[9]+(alphay[3]+alphavpar[2])*f[8]+f[3]*alphay[8]+f[2]*alphavpar[8]+alphay[4]*f[6]+f[4]*alphay[6]+alphavpar[4]*f[5]+f[4]*alphavpar[5]); 
  return cflFreq; 
} 
double EmGyrokineticSimpleHelicalStep2Vol2x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
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
