#include <GyrokineticSimpleHelicalModDecl.h> 
double EmGyrokineticSimpleHelicalVol3x2vSerP1_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *gradPar, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // bmagInv: 1/bmag.
  // gradPar: coefficient multiplying parallel gradient.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
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
  BstarXdBmag[0] = (0.7071067811865475*(2.828427124746191*BdriftX[0]*m_*wvpar+q_*(1.732050807568877*bmagInv[0]*Apar[2]*rdy2+Apar[0]*BdriftX[0])))/q_; 
  BstarXdBmag[1] = 0.7071067811865475*(1.732050807568877*bmagInv[0]*Apar[4]*rdy2+BdriftX[0]*Apar[1]); 
  BstarXdBmag[2] = 0.7071067811865475*BdriftX[0]*Apar[2]; 
  BstarXdBmag[3] = 0.408248290463863*(3.0*bmagInv[0]*Apar[6]*rdy2+1.732050807568877*BdriftX[0]*Apar[3]); 
  BstarXdBmag[4] = (1.154700538379252*BdriftX[0]*m_)/(q_*rdvpar2); 
  BstarXdBmag[6] = 0.7071067811865475*BdriftX[0]*Apar[4]; 
  BstarXdBmag[7] = 0.7071067811865475*(1.732050807568877*bmagInv[0]*Apar[7]*rdy2+BdriftX[0]*Apar[5]); 
  BstarXdBmag[8] = 0.7071067811865475*BdriftX[0]*Apar[6]; 
  BstarXdBmag[16] = 0.7071067811865475*BdriftX[0]*Apar[7]; 

  double BstarYdBmag[32]; 
  BstarYdBmag[0] = (0.7071067811865475*(2.828427124746191*BdriftY[0]*m_*wvpar+q_*(Apar[0]*BdriftY[0]-1.732050807568877*bmagInv[0]*Apar[1]*rdx2)))/q_; 
  BstarYdBmag[1] = 0.7071067811865475*BdriftY[0]*Apar[1]; 
  BstarYdBmag[2] = -0.7071067811865475*(1.732050807568877*bmagInv[0]*Apar[4]*rdx2-1.0*BdriftY[0]*Apar[2]); 
  BstarYdBmag[3] = -0.408248290463863*(3.0*bmagInv[0]*Apar[5]*rdx2-1.732050807568877*BdriftY[0]*Apar[3]); 
  BstarYdBmag[4] = (1.154700538379252*BdriftY[0]*m_)/(q_*rdvpar2); 
  BstarYdBmag[6] = 0.7071067811865475*BdriftY[0]*Apar[4]; 
  BstarYdBmag[7] = 0.7071067811865475*BdriftY[0]*Apar[5]; 
  BstarYdBmag[8] = -0.7071067811865475*(1.732050807568877*bmagInv[0]*Apar[7]*rdx2-1.0*BdriftY[0]*Apar[6]); 
  BstarYdBmag[16] = 0.7071067811865475*BdriftY[0]*Apar[7]; 

  double BstarZdBmag[32]; 
  BstarZdBmag[0] = 2.0*gradPar[0]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[32]; 
  alphax[0] = 0.3061862178478971*rdx2*((BstarXdBmag[0]*hamil[4]*rdvpar2)/m_-(2.0*bmagInv[0]*hamil[2]*rdy2)/q_); 
  alphax[1] = 0.3061862178478971*rdx2*((BstarXdBmag[1]*hamil[4]*rdvpar2)/m_-(2.0*bmagInv[0]*hamil[6]*rdy2)/q_); 
  alphax[2] = (0.3061862178478971*BstarXdBmag[2]*hamil[4]*rdvpar2*rdx2)/m_; 
  alphax[3] = 0.3061862178478971*rdx2*((BstarXdBmag[3]*hamil[4]*rdvpar2)/m_-(2.0*bmagInv[0]*hamil[8]*rdy2)/q_); 
  alphax[4] = (0.3061862178478971*BstarXdBmag[4]*hamil[4]*rdvpar2*rdx2)/m_; 
  alphax[6] = (0.3061862178478971*hamil[4]*BstarXdBmag[6]*rdvpar2*rdx2)/m_; 
  alphax[7] = 0.3061862178478971*rdx2*((hamil[4]*BstarXdBmag[7]*rdvpar2)/m_-(2.0*bmagInv[0]*hamil[16]*rdy2)/q_); 
  alphax[8] = (0.3061862178478971*hamil[4]*BstarXdBmag[8]*rdvpar2*rdx2)/m_; 
  alphax[16] = (0.3061862178478971*hamil[4]*BstarXdBmag[16]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphax[1]-1.414213562373095*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[32]; 
  alphay[0] = 0.3061862178478971*((2.0*bmagInv[0]*hamil[1]*rdx2)/q_+(BstarYdBmag[0]*hamil[4]*rdvpar2)/m_)*rdy2; 
  alphay[1] = (0.3061862178478971*BstarYdBmag[1]*hamil[4]*rdvpar2*rdy2)/m_; 
  alphay[2] = 0.3061862178478971*((2.0*bmagInv[0]*hamil[6]*rdx2)/q_+(BstarYdBmag[2]*hamil[4]*rdvpar2)/m_)*rdy2; 
  alphay[3] = 0.3061862178478971*((2.0*bmagInv[0]*hamil[7]*rdx2)/q_+(BstarYdBmag[3]*hamil[4]*rdvpar2)/m_)*rdy2; 
  alphay[4] = (0.3061862178478971*BstarYdBmag[4]*hamil[4]*rdvpar2*rdy2)/m_; 
  alphay[6] = (0.3061862178478971*hamil[4]*BstarYdBmag[6]*rdvpar2*rdy2)/m_; 
  alphay[7] = (0.3061862178478971*hamil[4]*BstarYdBmag[7]*rdvpar2*rdy2)/m_; 
  alphay[8] = 0.3061862178478971*((2.0*bmagInv[0]*hamil[16]*rdx2)/q_+(hamil[4]*BstarYdBmag[8]*rdvpar2)/m_)*rdy2; 
  alphay[16] = (0.3061862178478971*hamil[4]*BstarYdBmag[16]*rdvpar2*rdy2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphay[2]-1.414213562373095*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphay[2]+1.414213562373095*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphay[16]+alphay[8])-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphay[16]+alphay[8]))+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphay[16]+alphay[8])-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphay[16]+alphay[8]))+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphay[16]+alphay[8])-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphay[16]+alphay[8]))+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphay[16])+0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphay[16]+alphay[8])-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphay[16]-0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphay[16]+alphay[8]))+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[4]+alphay[3])-0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphay[16]+alphay[8]))-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphay[16]+alphay[8])+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphay[16]+alphay[8]))-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphay[16]+alphay[8])+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphay[16]+alphay[8]))-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphay[16]+alphay[8])+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]-0.1767766952966368*alphay[4]+0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphay[16]-0.3061862178478971*alphay[8]+0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphay[16]+alphay[8]))-0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*alphay[4]-0.1767766952966368*alphay[3]+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphay[16])+0.3061862178478971*alphay[8]-0.1767766952966368*alphay[7]-0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]-0.1767766952966368*alphay[1]+0.1767766952966368*alphay[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphay[16]+alphay[8])+0.1767766952966368*alphay[7]+0.3061862178478971*alphay[6]+0.1767766952966368*(alphay[4]+alphay[3])+0.3061862178478971*alphay[2]+0.1767766952966368*(alphay[1]+alphay[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphaz[32]; 
  alphaz[0] = (0.3061862178478971*BstarZdBmag[0]*hamil[4]*rdvpar2*rdz2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.0883883476483184*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0883883476483184*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.005524271728019897*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.005524271728019897*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[32]; 
  alphavpar[0] = -(0.3061862178478971*rdvpar2*(BstarZdBmag[0]*hamil[3]*rdz2+(BstarYdBmag[7]*hamil[16]+BstarYdBmag[3]*hamil[8]+BstarYdBmag[1]*hamil[6]+BstarYdBmag[0]*hamil[2])*rdy2+(BstarXdBmag[8]*hamil[16]+BstarXdBmag[3]*hamil[7]+BstarXdBmag[2]*hamil[6]+BstarXdBmag[0]*hamil[1])*rdx2))/m_; 
  alphavpar[1] = -(0.3061862178478971*rdvpar2*(BstarZdBmag[0]*hamil[7]*rdz2+(BstarYdBmag[3]*hamil[16]+BstarYdBmag[7]*hamil[8]+BstarYdBmag[0]*hamil[6]+BstarYdBmag[1]*hamil[2])*rdy2+(BstarXdBmag[16]*hamil[16]+BstarXdBmag[7]*hamil[7]+BstarXdBmag[6]*hamil[6]+BstarXdBmag[1]*hamil[1])*rdx2))/m_; 
  alphavpar[2] = -(0.3061862178478971*rdvpar2*(BstarZdBmag[0]*hamil[8]*rdz2+(BstarYdBmag[16]*hamil[16]+BstarYdBmag[8]*hamil[8]+BstarYdBmag[6]*hamil[6]+BstarYdBmag[2]*hamil[2])*rdy2+(BstarXdBmag[3]*hamil[16]+hamil[7]*BstarXdBmag[8]+BstarXdBmag[0]*hamil[6]+hamil[1]*BstarXdBmag[2])*rdx2))/m_; 
  alphavpar[3] = -(0.3061862178478971*rdvpar2*((BstarYdBmag[1]*hamil[16]+BstarYdBmag[0]*hamil[8]+hamil[6]*BstarYdBmag[7]+hamil[2]*BstarYdBmag[3])*rdy2+(BstarXdBmag[2]*hamil[16]+hamil[6]*BstarXdBmag[8]+BstarXdBmag[0]*hamil[7]+hamil[1]*BstarXdBmag[3])*rdx2))/m_; 
  alphavpar[4] = -(0.3061862178478971*rdvpar2*(hamil[2]*BstarYdBmag[4]*rdy2+hamil[1]*BstarXdBmag[4]*rdx2))/m_; 
  alphavpar[6] = -(0.3061862178478971*rdvpar2*(BstarZdBmag[0]*hamil[16]*rdz2+(BstarYdBmag[8]*hamil[16]+hamil[8]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[6]+hamil[2]*BstarYdBmag[6])*rdy2+(BstarXdBmag[7]*hamil[16]+hamil[7]*BstarXdBmag[16]+BstarXdBmag[1]*hamil[6]+hamil[1]*BstarXdBmag[6])*rdx2))/m_; 
  alphavpar[7] = -(0.3061862178478971*rdvpar2*((BstarYdBmag[0]*hamil[16]+BstarYdBmag[1]*hamil[8]+hamil[2]*BstarYdBmag[7]+BstarYdBmag[3]*hamil[6])*rdy2+(BstarXdBmag[6]*hamil[16]+hamil[6]*BstarXdBmag[16]+BstarXdBmag[1]*hamil[7]+hamil[1]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[8] = -(0.3061862178478971*rdvpar2*((BstarYdBmag[6]*hamil[16]+hamil[6]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[8]+hamil[2]*BstarYdBmag[8])*rdy2+(BstarXdBmag[0]*hamil[16]+hamil[1]*BstarXdBmag[8]+BstarXdBmag[2]*hamil[7]+BstarXdBmag[3]*hamil[6])*rdx2))/m_; 
  alphavpar[9] = -(0.3061862178478971*BstarYdBmag[4]*hamil[6]*rdvpar2*rdy2)/m_; 
  alphavpar[10] = -(0.3061862178478971*BstarXdBmag[4]*hamil[6]*rdvpar2*rdx2)/m_; 
  alphavpar[11] = -(0.3061862178478971*rdvpar2*(BstarYdBmag[4]*hamil[8]*rdy2+BstarXdBmag[4]*hamil[7]*rdx2))/m_; 
  alphavpar[16] = -(0.3061862178478971*rdvpar2*((BstarYdBmag[2]*hamil[16]+hamil[2]*BstarYdBmag[16]+BstarYdBmag[6]*hamil[8]+hamil[6]*BstarYdBmag[8])*rdy2+(BstarXdBmag[1]*hamil[16]+hamil[1]*BstarXdBmag[16]+BstarXdBmag[6]*hamil[7]+hamil[6]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[18] = -(0.3061862178478971*BstarYdBmag[4]*hamil[16]*rdvpar2*rdy2)/m_; 
  alphavpar[19] = -(0.3061862178478971*BstarXdBmag[4]*hamil[16]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphavpar[4]-1.414213562373095*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphavpar[4]+1.414213562373095*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.3061862178478971*(alphavpar[19]+alphavpar[18]))-0.1767766952966368*alphavpar[16]+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[19])+0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[19])+0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphavpar[19]+alphavpar[18]))+0.1767766952966368*alphavpar[16]-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphavpar[19]+alphavpar[18]))-0.1767766952966368*alphavpar[16]+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[19])+0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphavpar[19])+0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphavpar[19]+alphavpar[18]))+0.1767766952966368*alphavpar[16]-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*alphavpar[16]-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[19])+0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphavpar[19]+alphavpar[18]))-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphavpar[19]+alphavpar[18]))+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[19])+0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*alphavpar[16]+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.1767766952966368*alphavpar[16]-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[19])+0.3061862178478971*alphavpar[18]+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphavpar[19]+alphavpar[18]))-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphavpar[19]+alphavpar[18]))+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphavpar[19])+0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]-0.1767766952966368*alphavpar[16]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.1767766952966368*alphavpar[16]+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.3061862178478971*(alphax[16]*f[16]+alphax[8]*f[8]+alphax[7]*f[7]+alphax[6]*f[6]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.3061862178478971*(alphay[16]*f[16]+alphay[8]*f[8]+alphay[7]*f[7]+alphay[6]*f[6]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.3061862178478971*alphaz[0]*f[0]; 
  out[4] += 0.3061862178478971*(alphavpar[19]*f[19]+alphavpar[18]*f[18]+alphavpar[16]*f[16]+alphavpar[11]*f[11]+alphavpar[10]*f[10]+alphavpar[9]*f[9]+alphavpar[8]*f[8]+alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[6] += 0.3061862178478971*((alphay[8]+alphax[7])*f[16]+f[8]*alphay[16]+f[7]*alphax[16]+alphax[4]*f[10]+alphay[4]*f[9]+alphax[3]*f[8]+f[3]*alphax[8]+alphay[3]*f[7]+f[3]*alphay[7]+(alphay[2]+alphax[1])*f[6]+f[2]*alphay[6]+f[1]*alphax[6]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[7] += 0.3061862178478971*(alphax[6]*f[16]+f[6]*alphax[16]+alphax[4]*f[11]+alphax[2]*f[8]+f[2]*alphax[8]+alphax[1]*f[7]+f[1]*alphax[7]+alphax[0]*f[3]+f[0]*alphax[3]+alphaz[0]*f[1]); 
  out[8] += 0.3061862178478971*(alphay[6]*f[16]+f[6]*alphay[16]+alphay[4]*f[11]+alphay[2]*f[8]+f[2]*alphay[8]+alphay[1]*f[7]+f[1]*alphay[7]+alphay[0]*f[3]+f[0]*alphay[3]+alphaz[0]*f[2]); 
  out[9] += 0.3061862178478971*((alphavpar[19]+alphax[16])*f[26]+alphax[8]*f[19]+(alphavpar[11]+alphax[7])*f[18]+f[11]*alphavpar[18]+(alphavpar[10]+alphax[6])*f[17]+alphavpar[8]*f[16]+f[8]*alphavpar[16]+alphax[3]*f[11]+alphax[2]*f[10]+(alphavpar[4]+alphax[1])*f[9]+f[4]*alphavpar[9]+alphavpar[3]*f[7]+f[3]*alphavpar[7]+alphavpar[2]*f[6]+f[2]*alphavpar[6]+alphax[0]*f[4]+f[0]*alphax[4]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[10] += 0.3061862178478971*((alphavpar[18]+alphay[16])*f[26]+(alphavpar[11]+alphay[8])*f[19]+f[11]*alphavpar[19]+alphay[7]*f[18]+(alphavpar[9]+alphay[6])*f[17]+alphavpar[7]*f[16]+f[7]*alphavpar[16]+alphay[3]*f[11]+(alphavpar[4]+alphay[2])*f[10]+f[4]*alphavpar[10]+alphay[1]*f[9]+alphavpar[3]*f[8]+f[3]*alphavpar[8]+alphavpar[1]*f[6]+f[1]*alphavpar[6]+alphay[0]*f[4]+f[0]*alphay[4]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[11] += 0.3061862178478971*(alphavpar[10]*f[19]+f[10]*alphavpar[19]+alphavpar[9]*f[18]+f[9]*alphavpar[18]+alphavpar[6]*f[16]+f[6]*alphavpar[16]+alphavpar[4]*f[11]+f[4]*alphavpar[11]+alphavpar[2]*f[8]+f[2]*alphavpar[8]+alphavpar[1]*f[7]+f[1]*alphavpar[7]+alphaz[0]*f[4]+alphavpar[0]*f[3]+f[0]*alphavpar[3]); 
  out[12] += 0.3061862178478971*(alphax[16]*f[27]+alphax[8]*f[22]+alphax[7]*f[21]+alphax[6]*f[20]+alphax[4]*f[15]+alphax[3]*f[14]+alphax[2]*f[13]+alphax[1]*f[12]+alphax[0]*f[5]); 
  out[13] += 0.3061862178478971*(alphay[16]*f[27]+alphay[8]*f[22]+alphay[7]*f[21]+alphay[6]*f[20]+alphay[4]*f[15]+alphay[3]*f[14]+alphay[2]*f[13]+alphay[1]*f[12]+alphay[0]*f[5]); 
  out[14] += 0.3061862178478971*alphaz[0]*f[5]; 
  out[15] += 0.3061862178478971*(alphavpar[19]*f[30]+alphavpar[18]*f[29]+alphavpar[16]*f[27]+alphavpar[11]*f[25]+alphavpar[10]*f[24]+alphavpar[9]*f[23]+alphavpar[8]*f[22]+alphavpar[7]*f[21]+alphavpar[6]*f[20]+alphavpar[4]*f[15]+alphavpar[3]*f[14]+alphavpar[2]*f[13]+alphavpar[1]*f[12]+alphavpar[0]*f[5]); 
  out[16] += 0.3061862178478971*(alphax[4]*f[19]+alphay[4]*f[18]+(alphay[2]+alphax[1])*f[16]+f[2]*alphay[16]+f[1]*alphax[16]+(alphay[6]+alphax[0])*f[8]+f[6]*alphay[8]+f[0]*alphax[8]+(alphax[6]+alphay[0])*f[7]+f[0]*alphay[7]+f[6]*(alphax[7]+alphaz[0])+(alphax[2]+alphay[1])*f[3]+f[1]*alphay[3]+f[2]*alphax[3]); 
  out[17] += 0.3061862178478971*((alphavpar[11]+alphay[8]+alphax[7])*f[26]+(alphavpar[18]+alphay[16]+alphax[3])*f[19]+f[18]*(alphavpar[19]+alphax[16]+alphay[3])+(alphavpar[4]+alphay[2]+alphax[1])*f[17]+alphavpar[3]*f[16]+f[3]*alphavpar[16]+(alphax[8]+alphay[7])*f[11]+(alphavpar[9]+alphay[6]+alphax[0])*f[10]+f[9]*(alphavpar[10]+alphax[6]+alphay[0])+alphavpar[7]*f[8]+f[7]*alphavpar[8]+alphavpar[0]*f[6]+f[0]*alphavpar[6]+(alphax[2]+alphay[1])*f[4]+f[1]*alphay[4]+f[2]*(alphax[4]+alphavpar[1])+f[1]*alphavpar[2]); 
  out[18] += 0.3061862178478971*((alphavpar[10]+alphax[6])*f[26]+alphax[2]*f[19]+f[17]*alphavpar[19]+(alphavpar[4]+alphax[1])*f[18]+f[4]*alphavpar[18]+alphax[16]*f[17]+alphavpar[2]*f[16]+f[2]*alphavpar[16]+(alphavpar[9]+alphax[0])*f[11]+f[9]*alphavpar[11]+alphax[8]*f[10]+(alphax[7]+alphaz[0])*f[9]+alphavpar[6]*f[8]+f[6]*alphavpar[8]+alphavpar[0]*f[7]+f[0]*alphavpar[7]+alphax[3]*f[4]+f[3]*(alphax[4]+alphavpar[1])+f[1]*alphavpar[3]); 
  out[19] += 0.3061862178478971*((alphavpar[9]+alphay[6])*f[26]+(alphavpar[4]+alphay[2])*f[19]+f[4]*alphavpar[19]+alphay[1]*f[18]+f[17]*(alphavpar[18]+alphay[16])+alphavpar[1]*f[16]+f[1]*alphavpar[16]+(alphavpar[10]+alphay[0])*f[11]+f[10]*(alphavpar[11]+alphay[8]+alphaz[0])+alphay[7]*f[9]+alphavpar[0]*f[8]+f[0]*alphavpar[8]+alphavpar[6]*f[7]+f[6]*alphavpar[7]+alphay[3]*f[4]+f[3]*(alphay[4]+alphavpar[2])+f[2]*alphavpar[3]); 
  out[20] += 0.3061862178478971*((alphay[8]+alphax[7])*f[27]+alphax[4]*f[24]+alphay[4]*f[23]+(alphay[16]+alphax[3])*f[22]+(alphax[16]+alphay[3])*f[21]+(alphay[2]+alphax[1])*f[20]+(alphax[8]+alphay[7])*f[14]+(alphay[6]+alphax[0])*f[13]+(alphax[6]+alphay[0])*f[12]+(alphax[2]+alphay[1])*f[5]); 
  out[21] += 0.3061862178478971*(alphax[6]*f[27]+alphax[4]*f[25]+alphax[2]*f[22]+alphax[1]*f[21]+alphax[16]*f[20]+alphax[0]*f[14]+alphax[8]*f[13]+(alphax[7]+alphaz[0])*f[12]+alphax[3]*f[5]); 
  out[22] += 0.3061862178478971*(alphay[6]*f[27]+alphay[4]*f[25]+alphay[2]*f[22]+alphay[1]*f[21]+alphay[16]*f[20]+alphay[0]*f[14]+(alphay[8]+alphaz[0])*f[13]+alphay[7]*f[12]+alphay[3]*f[5]); 
  out[23] += 0.3061862178478971*((alphavpar[19]+alphax[16])*f[31]+alphax[8]*f[30]+(alphavpar[11]+alphax[7])*f[29]+(alphavpar[10]+alphax[6])*f[28]+alphavpar[8]*f[27]+(alphavpar[18]+alphax[3])*f[25]+alphax[2]*f[24]+(alphavpar[4]+alphax[1])*f[23]+alphavpar[16]*f[22]+alphavpar[3]*f[21]+alphavpar[2]*f[20]+(alphavpar[9]+alphax[0])*f[15]+alphavpar[7]*f[14]+alphavpar[6]*f[13]+alphavpar[0]*f[12]+(alphax[4]+alphavpar[1])*f[5]); 
  out[24] += 0.3061862178478971*((alphavpar[18]+alphay[16])*f[31]+(alphavpar[11]+alphay[8])*f[30]+alphay[7]*f[29]+(alphavpar[9]+alphay[6])*f[28]+alphavpar[7]*f[27]+(alphavpar[19]+alphay[3])*f[25]+(alphavpar[4]+alphay[2])*f[24]+alphay[1]*f[23]+alphavpar[3]*f[22]+alphavpar[16]*f[21]+alphavpar[1]*f[20]+(alphavpar[10]+alphay[0])*f[15]+alphavpar[8]*f[14]+alphavpar[0]*f[13]+alphavpar[6]*f[12]+(alphay[4]+alphavpar[2])*f[5]); 
  out[25] += 0.3061862178478971*(alphavpar[10]*f[30]+alphavpar[9]*f[29]+alphavpar[6]*f[27]+alphavpar[4]*f[25]+alphavpar[19]*f[24]+alphavpar[18]*f[23]+alphavpar[2]*f[22]+alphavpar[1]*f[21]+alphavpar[16]*f[20]+(alphavpar[11]+alphaz[0])*f[15]+alphavpar[0]*f[14]+alphavpar[8]*f[13]+alphavpar[7]*f[12]+alphavpar[3]*f[5]); 
  out[26] += 0.3061862178478971*((alphavpar[4]+alphay[2]+alphax[1])*f[26]+(alphavpar[9]+alphay[6]+alphax[0])*f[19]+f[9]*alphavpar[19]+(alphavpar[10]+alphax[6]+alphay[0])*f[18]+f[10]*alphavpar[18]+(alphavpar[11]+alphay[8]+alphax[7]+alphaz[0])*f[17]+alphavpar[0]*f[16]+f[10]*alphay[16]+f[9]*alphax[16]+f[0]*alphavpar[16]+(alphax[2]+alphay[1])*f[11]+alphax[3]*f[10]+alphay[3]*f[9]+(alphax[4]+alphavpar[1])*f[8]+f[4]*alphax[8]+f[1]*alphavpar[8]+(alphay[4]+alphavpar[2])*f[7]+f[4]*alphay[7]+f[2]*alphavpar[7]+alphavpar[3]*f[6]+f[3]*alphavpar[6]); 
  out[27] += 0.3061862178478971*(alphax[4]*f[30]+alphay[4]*f[29]+(alphay[2]+alphax[1])*f[27]+(alphay[6]+alphax[0])*f[22]+(alphax[6]+alphay[0])*f[21]+(alphay[8]+alphax[7]+alphaz[0])*f[20]+f[13]*alphay[16]+f[12]*alphax[16]+(alphax[2]+alphay[1])*f[14]+alphax[3]*f[13]+alphay[3]*f[12]+f[5]*(alphax[8]+alphay[7])); 
  out[28] += 0.3061862178478971*((alphavpar[11]+alphay[8]+alphax[7])*f[31]+(alphavpar[18]+alphay[16]+alphax[3])*f[30]+(alphavpar[19]+alphax[16]+alphay[3])*f[29]+(alphavpar[4]+alphay[2]+alphax[1])*f[28]+alphavpar[3]*f[27]+(alphax[8]+alphay[7])*f[25]+(alphavpar[9]+alphay[6]+alphax[0])*f[24]+(alphavpar[10]+alphax[6]+alphay[0])*f[23]+alphavpar[7]*f[22]+alphavpar[8]*f[21]+alphavpar[0]*f[20]+f[14]*alphavpar[16]+(alphax[2]+alphay[1])*f[15]+(alphax[4]+alphavpar[1])*f[13]+(alphay[4]+alphavpar[2])*f[12]+f[5]*alphavpar[6]); 
  out[29] += 0.3061862178478971*((alphavpar[10]+alphax[6])*f[31]+alphax[2]*f[30]+(alphavpar[4]+alphax[1])*f[29]+(alphavpar[19]+alphax[16])*f[28]+alphavpar[2]*f[27]+(alphavpar[9]+alphax[0])*f[25]+alphax[8]*f[24]+(alphavpar[11]+alphax[7]+alphaz[0])*f[23]+alphavpar[6]*f[22]+alphavpar[0]*f[21]+alphavpar[8]*f[20]+f[15]*alphavpar[18]+f[13]*alphavpar[16]+alphax[3]*f[15]+(alphax[4]+alphavpar[1])*f[14]+alphavpar[3]*f[12]+f[5]*alphavpar[7]); 
  out[30] += 0.3061862178478971*((alphavpar[9]+alphay[6])*f[31]+(alphavpar[4]+alphay[2])*f[30]+alphay[1]*f[29]+(alphavpar[18]+alphay[16])*f[28]+alphavpar[1]*f[27]+(alphavpar[10]+alphay[0])*f[25]+(alphavpar[11]+alphay[8]+alphaz[0])*f[24]+alphay[7]*f[23]+alphavpar[0]*f[22]+alphavpar[6]*f[21]+alphavpar[7]*f[20]+f[15]*alphavpar[19]+f[12]*alphavpar[16]+alphay[3]*f[15]+(alphay[4]+alphavpar[2])*f[14]+alphavpar[3]*f[13]+f[5]*alphavpar[8]); 
  out[31] += 0.3061862178478971*((alphavpar[4]+alphay[2]+alphax[1])*f[31]+(alphavpar[9]+alphay[6]+alphax[0])*f[30]+(alphavpar[10]+alphax[6]+alphay[0])*f[29]+(alphavpar[11]+alphay[8]+alphax[7]+alphaz[0])*f[28]+alphavpar[0]*f[27]+(alphax[2]+alphay[1])*f[25]+(alphavpar[18]+alphay[16]+alphax[3])*f[24]+(alphavpar[19]+alphax[16]+alphay[3])*f[23]+(alphax[4]+alphavpar[1])*f[22]+(alphay[4]+alphavpar[2])*f[21]+alphavpar[3]*f[20]+f[5]*alphavpar[16]+(alphax[8]+alphay[7])*f[15]+alphavpar[6]*f[14]+alphavpar[7]*f[13]+alphavpar[8]*f[12]); 
  return cflFreq; 
} 
double EmGyrokineticSimpleHelicalVol3x2vSerP1_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *bmagInv, const double *gradPar, const double *BdriftX, const double *BdriftY, const double *phi, const double *Apar, const double *dApardt, const double *f, double *out) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // q_,m_: species charge and mass.
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // bmag: magnetic field amplitude.
  // bmagInv: 1/bmag.
  // gradPar: coefficient multiplying parallel gradient.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
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
  BstarXdBmag[0] = (0.7071067811865475*(2.828427124746191*BdriftX[0]*m_*wvpar+q_*(1.732050807568877*(bmagInv[1]*Apar[4]+bmagInv[0]*Apar[2])*rdy2+Apar[1]*BdriftX[1]+Apar[0]*BdriftX[0])))/q_; 
  BstarXdBmag[1] = (0.7071067811865475*(2.828427124746191*BdriftX[1]*m_*wvpar+q_*(1.732050807568877*(bmagInv[0]*Apar[4]+bmagInv[1]*Apar[2])*rdy2+Apar[0]*BdriftX[1]+BdriftX[0]*Apar[1])))/q_; 
  BstarXdBmag[2] = 0.7071067811865475*(BdriftX[1]*Apar[4]+BdriftX[0]*Apar[2]); 
  BstarXdBmag[3] = 0.408248290463863*(3.0*(bmagInv[1]*Apar[7]+bmagInv[0]*Apar[6])*rdy2+1.732050807568877*(BdriftX[1]*Apar[5]+BdriftX[0]*Apar[3])); 
  BstarXdBmag[4] = (1.154700538379252*BdriftX[0]*m_)/(q_*rdvpar2); 
  BstarXdBmag[6] = 0.7071067811865475*(BdriftX[0]*Apar[4]+BdriftX[1]*Apar[2]); 
  BstarXdBmag[7] = 0.7071067811865475*(1.732050807568877*(bmagInv[0]*Apar[7]+bmagInv[1]*Apar[6])*rdy2+BdriftX[0]*Apar[5]+BdriftX[1]*Apar[3]); 
  BstarXdBmag[8] = 0.7071067811865475*(BdriftX[1]*Apar[7]+BdriftX[0]*Apar[6]); 
  BstarXdBmag[9] = (1.154700538379252*BdriftX[1]*m_)/(q_*rdvpar2); 
  BstarXdBmag[16] = 0.7071067811865475*(BdriftX[0]*Apar[7]+BdriftX[1]*Apar[6]); 

  double BstarYdBmag[32]; 
  BstarYdBmag[0] = (0.7071067811865475*(2.828427124746191*BdriftY[0]*m_*wvpar+q_*(Apar[1]*(BdriftY[1]-1.732050807568877*bmagInv[0]*rdx2)+Apar[0]*BdriftY[0])))/q_; 
  BstarYdBmag[1] = (0.7071067811865475*(2.828427124746191*BdriftY[1]*m_*wvpar+q_*((-1.732050807568877*Apar[1]*bmagInv[1]*rdx2)+Apar[0]*BdriftY[1]+BdriftY[0]*Apar[1])))/q_; 
  BstarYdBmag[2] = -0.7071067811865475*(1.732050807568877*bmagInv[0]*Apar[4]*rdx2-1.0*(BdriftY[1]*Apar[4]+BdriftY[0]*Apar[2])); 
  BstarYdBmag[3] = -0.408248290463863*(3.0*bmagInv[0]*Apar[5]*rdx2-1.732050807568877*(BdriftY[1]*Apar[5]+BdriftY[0]*Apar[3])); 
  BstarYdBmag[4] = (1.154700538379252*BdriftY[0]*m_)/(q_*rdvpar2); 
  BstarYdBmag[6] = -0.7071067811865475*(1.732050807568877*bmagInv[1]*Apar[4]*rdx2-1.0*(BdriftY[0]*Apar[4]+BdriftY[1]*Apar[2])); 
  BstarYdBmag[7] = -0.7071067811865475*(1.732050807568877*bmagInv[1]*Apar[5]*rdx2-1.0*(BdriftY[0]*Apar[5]+BdriftY[1]*Apar[3])); 
  BstarYdBmag[8] = -0.7071067811865475*(1.732050807568877*bmagInv[0]*Apar[7]*rdx2-1.0*(BdriftY[1]*Apar[7]+BdriftY[0]*Apar[6])); 
  BstarYdBmag[9] = (1.154700538379252*BdriftY[1]*m_)/(q_*rdvpar2); 
  BstarYdBmag[16] = -0.7071067811865475*(1.732050807568877*bmagInv[1]*Apar[7]*rdx2-1.0*(BdriftY[0]*Apar[7]+BdriftY[1]*Apar[6])); 

  double BstarZdBmag[32]; 
  BstarZdBmag[0] = 2.0*gradPar[0]; 
  BstarZdBmag[1] = 2.0*gradPar[1]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alphax[32]; 
  alphax[0] = -0.3061862178478971*rdx2*((2.0*(bmagInv[1]*hamil[6]+bmagInv[0]*hamil[2])*rdy2)/q_-(1.0*BstarXdBmag[0]*hamil[4]*rdvpar2)/m_); 
  alphax[1] = -0.3061862178478971*rdx2*((2.0*(bmagInv[0]*hamil[6]+bmagInv[1]*hamil[2])*rdy2)/q_-(1.0*BstarXdBmag[1]*hamil[4]*rdvpar2)/m_); 
  alphax[2] = (0.3061862178478971*BstarXdBmag[2]*hamil[4]*rdvpar2*rdx2)/m_; 
  alphax[3] = -0.3061862178478971*rdx2*((2.0*(bmagInv[1]*hamil[16]+bmagInv[0]*hamil[8])*rdy2)/q_-(1.0*BstarXdBmag[3]*hamil[4]*rdvpar2)/m_); 
  alphax[4] = (0.3061862178478971*BstarXdBmag[4]*hamil[4]*rdvpar2*rdx2)/m_; 
  alphax[6] = (0.3061862178478971*hamil[4]*BstarXdBmag[6]*rdvpar2*rdx2)/m_; 
  alphax[7] = -0.3061862178478971*rdx2*((2.0*(bmagInv[0]*hamil[16]+bmagInv[1]*hamil[8])*rdy2)/q_-(1.0*hamil[4]*BstarXdBmag[7]*rdvpar2)/m_); 
  alphax[8] = (0.3061862178478971*hamil[4]*BstarXdBmag[8]*rdvpar2*rdx2)/m_; 
  alphax[9] = (0.3061862178478971*hamil[4]*BstarXdBmag[9]*rdvpar2*rdx2)/m_; 
  alphax[16] = (0.3061862178478971*hamil[4]*BstarXdBmag[16]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphax[1]-1.414213562373095*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphax[1]+1.414213562373095*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[16]+alphax[9])-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[16]+alphax[9])-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[16]+alphax[9]))+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[16]+alphax[9]))+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[16]+alphax[9])-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*(alphax[16]+alphax[9])-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*alphax[16])+0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[16]+alphax[9]))+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.3061862178478971*alphax[16]-0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.3061862178478971*(alphax[16]+alphax[9]))+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])-0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.3061862178478971*alphax[16]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[16]+alphax[9]))-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[16]+alphax[9]))-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[16]+alphax[9])+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[16]+alphax[9])+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[16]+alphax[9]))-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]-0.1767766952966368*(alphax[4]+alphax[3])+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*(alphax[16]+alphax[9]))-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]-0.1767766952966368*alphax[4]+0.1767766952966368*alphax[3]-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*alphax[16]-0.3061862178478971*alphax[9]+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])-0.1767766952966368*alphax[4]+0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[16]+alphax[9])+0.1767766952966368*alphax[8]-0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*alphax[4]-0.1767766952966368*(alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]-0.3061862178478971*alphax[7]+0.3061862178478971*alphax[6]+0.1767766952966368*alphax[4]-0.1767766952966368*alphax[3]+0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.3061862178478971*alphax[16])+0.3061862178478971*alphax[9]-0.1767766952966368*alphax[8]+0.3061862178478971*alphax[7]-0.3061862178478971*alphax[6]+0.1767766952966368*(alphax[4]+alphax[3])-0.1767766952966368*alphax[2]+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.3061862178478971*(alphax[16]+alphax[9])+0.1767766952966368*alphax[8]+0.3061862178478971*(alphax[7]+alphax[6])+0.1767766952966368*(alphax[4]+alphax[3]+alphax[2])+0.3061862178478971*alphax[1]+0.1767766952966368*alphax[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphay[32]; 
  alphay[0] = 0.3061862178478971*((2.0*bmagInv[0]*hamil[1]*rdx2)/q_+(BstarYdBmag[0]*hamil[4]*rdvpar2)/m_)*rdy2; 
  alphay[1] = 0.3061862178478971*((2.0*bmagInv[1]*hamil[1]*rdx2)/q_+(BstarYdBmag[1]*hamil[4]*rdvpar2)/m_)*rdy2; 
  alphay[2] = 0.3061862178478971*((2.0*bmagInv[0]*hamil[6]*rdx2)/q_+(BstarYdBmag[2]*hamil[4]*rdvpar2)/m_)*rdy2; 
  alphay[3] = 0.3061862178478971*((2.0*bmagInv[0]*hamil[7]*rdx2)/q_+(BstarYdBmag[3]*hamil[4]*rdvpar2)/m_)*rdy2; 
  alphay[4] = (0.3061862178478971*BstarYdBmag[4]*hamil[4]*rdvpar2*rdy2)/m_; 
  alphay[5] = (0.6123724356957944*bmagInv[0]*hamil[12]*rdx2*rdy2)/q_; 
  alphay[6] = 0.3061862178478971*((2.0*bmagInv[1]*hamil[6]*rdx2)/q_+(hamil[4]*BstarYdBmag[6]*rdvpar2)/m_)*rdy2; 
  alphay[7] = 0.3061862178478971*((2.0*bmagInv[1]*hamil[7]*rdx2)/q_+(hamil[4]*BstarYdBmag[7]*rdvpar2)/m_)*rdy2; 
  alphay[8] = 0.3061862178478971*((2.0*bmagInv[0]*hamil[16]*rdx2)/q_+(hamil[4]*BstarYdBmag[8]*rdvpar2)/m_)*rdy2; 
  alphay[9] = (0.3061862178478971*hamil[4]*BstarYdBmag[9]*rdvpar2*rdy2)/m_; 
  alphay[12] = (0.6123724356957944*bmagInv[1]*hamil[12]*rdx2*rdy2)/q_; 
  alphay[16] = 0.3061862178478971*((2.0*bmagInv[1]*hamil[16]*rdx2)/q_+(hamil[4]*BstarYdBmag[16]*rdvpar2)/m_)*rdy2; 
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
  alphaz[0] = (0.3061862178478971*BstarZdBmag[0]*hamil[4]*rdvpar2*rdz2)/m_; 
  alphaz[1] = (0.3061862178478971*BstarZdBmag[1]*hamil[4]*rdvpar2*rdz2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = 0.0883883476483184*alphaz[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0883883476483184*alphaz[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphaz[0]-0.1767766952966368*alphaz[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0055242717280199*(alphaz[1]+alphaz[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  double alphavpar[32]; 
  alphavpar[0] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[1]*hamil[7]+BstarZdBmag[0]*hamil[3])*rdz2+(BstarYdBmag[7]*hamil[16]+BstarYdBmag[3]*hamil[8]+BstarYdBmag[1]*hamil[6]+BstarYdBmag[0]*hamil[2])*rdy2+(BstarXdBmag[8]*hamil[16]+BstarXdBmag[3]*hamil[7]+BstarXdBmag[2]*hamil[6]+BstarXdBmag[0]*hamil[1])*rdx2))/m_; 
  alphavpar[1] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[0]*hamil[7]+BstarZdBmag[1]*hamil[3])*rdz2+(BstarYdBmag[3]*hamil[16]+BstarYdBmag[7]*hamil[8]+BstarYdBmag[0]*hamil[6]+BstarYdBmag[1]*hamil[2])*rdy2+(BstarXdBmag[16]*hamil[16]+BstarXdBmag[7]*hamil[7]+BstarXdBmag[6]*hamil[6]+BstarXdBmag[1]*hamil[1])*rdx2))/m_; 
  alphavpar[2] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[1]*hamil[16]+BstarZdBmag[0]*hamil[8])*rdz2+(BstarYdBmag[16]*hamil[16]+BstarYdBmag[8]*hamil[8]+BstarYdBmag[6]*hamil[6]+BstarYdBmag[2]*hamil[2])*rdy2+(BstarXdBmag[3]*hamil[16]+hamil[7]*BstarXdBmag[8]+BstarXdBmag[0]*hamil[6]+hamil[1]*BstarXdBmag[2])*rdx2))/m_; 
  alphavpar[3] = -(0.3061862178478971*rdvpar2*((BstarYdBmag[1]*hamil[16]+BstarYdBmag[0]*hamil[8]+hamil[6]*BstarYdBmag[7]+hamil[2]*BstarYdBmag[3])*rdy2+(BstarXdBmag[2]*hamil[16]+hamil[6]*BstarXdBmag[8]+BstarXdBmag[0]*hamil[7]+hamil[1]*BstarXdBmag[3])*rdx2))/m_; 
  alphavpar[4] = -(0.3061862178478971*rdvpar2*((hamil[6]*BstarYdBmag[9]+hamil[2]*BstarYdBmag[4])*rdy2+hamil[1]*BstarXdBmag[4]*rdx2))/m_; 
  alphavpar[5] = -(0.3061862178478971*BstarXdBmag[0]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[6] = -(0.3061862178478971*rdvpar2*((BstarZdBmag[0]*hamil[16]+BstarZdBmag[1]*hamil[8])*rdz2+(BstarYdBmag[8]*hamil[16]+hamil[8]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[6]+hamil[2]*BstarYdBmag[6])*rdy2+(BstarXdBmag[7]*hamil[16]+hamil[7]*BstarXdBmag[16]+BstarXdBmag[1]*hamil[6]+hamil[1]*BstarXdBmag[6])*rdx2))/m_; 
  alphavpar[7] = -(0.3061862178478971*rdvpar2*((BstarYdBmag[0]*hamil[16]+BstarYdBmag[1]*hamil[8]+hamil[2]*BstarYdBmag[7]+BstarYdBmag[3]*hamil[6])*rdy2+(BstarXdBmag[6]*hamil[16]+hamil[6]*BstarXdBmag[16]+BstarXdBmag[1]*hamil[7]+hamil[1]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[8] = -(0.3061862178478971*rdvpar2*((BstarYdBmag[6]*hamil[16]+hamil[6]*BstarYdBmag[16]+BstarYdBmag[2]*hamil[8]+hamil[2]*BstarYdBmag[8])*rdy2+(BstarXdBmag[0]*hamil[16]+hamil[1]*BstarXdBmag[8]+BstarXdBmag[2]*hamil[7]+BstarXdBmag[3]*hamil[6])*rdx2))/m_; 
  alphavpar[9] = -(0.3061862178478971*rdvpar2*((hamil[2]*BstarYdBmag[9]+BstarYdBmag[4]*hamil[6])*rdy2+hamil[1]*BstarXdBmag[9]*rdx2))/m_; 
  alphavpar[10] = -(0.3061862178478971*BstarXdBmag[4]*hamil[6]*rdvpar2*rdx2)/m_; 
  alphavpar[11] = -(0.3061862178478971*rdvpar2*((BstarYdBmag[9]*hamil[16]+BstarYdBmag[4]*hamil[8])*rdy2+BstarXdBmag[4]*hamil[7]*rdx2))/m_; 
  alphavpar[12] = -(0.3061862178478971*BstarXdBmag[1]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[13] = -(0.3061862178478971*BstarXdBmag[2]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[14] = -(0.3061862178478971*BstarXdBmag[3]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[15] = -(0.3061862178478971*BstarXdBmag[4]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[16] = -(0.3061862178478971*rdvpar2*((BstarYdBmag[2]*hamil[16]+hamil[2]*BstarYdBmag[16]+BstarYdBmag[6]*hamil[8]+hamil[6]*BstarYdBmag[8])*rdy2+(BstarXdBmag[1]*hamil[16]+hamil[1]*BstarXdBmag[16]+BstarXdBmag[6]*hamil[7]+hamil[6]*BstarXdBmag[7])*rdx2))/m_; 
  alphavpar[17] = -(0.3061862178478971*hamil[6]*BstarXdBmag[9]*rdvpar2*rdx2)/m_; 
  alphavpar[18] = -(0.3061862178478971*rdvpar2*((BstarYdBmag[4]*hamil[16]+hamil[8]*BstarYdBmag[9])*rdy2+hamil[7]*BstarXdBmag[9]*rdx2))/m_; 
  alphavpar[19] = -(0.3061862178478971*BstarXdBmag[4]*hamil[16]*rdvpar2*rdx2)/m_; 
  alphavpar[20] = -(0.3061862178478971*BstarXdBmag[6]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[21] = -(0.3061862178478971*BstarXdBmag[7]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[22] = -(0.3061862178478971*BstarXdBmag[8]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[23] = -(0.3061862178478971*BstarXdBmag[9]*hamil[12]*rdvpar2*rdx2)/m_; 
  alphavpar[26] = -(0.3061862178478971*BstarXdBmag[9]*hamil[16]*rdvpar2*rdx2)/m_; 
  alphavpar[27] = -(0.3061862178478971*hamil[12]*BstarXdBmag[16]*rdvpar2*rdx2)/m_; 
#if cflType == SURFAVG 
  // Evaluate surface-averaged alpha on left.
  alphaL = -0.0625*(2.449489742783178*alphavpar[4]-1.414213562373095*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate surface-averaged alpha on right.
  alphaR = 0.0625*(2.449489742783178*alphavpar[4]+1.414213562373095*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#elif cflType == QUAD 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.03125*(0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[27])-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*alphavpar[19]+0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[27])-0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.1767766952966368*alphavpar[20]+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[27])-0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.1767766952966368*alphavpar[20]+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*alphavpar[19]+0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[27])-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[27])+0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*alphavpar[19]+0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[27])+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.1767766952966368*alphavpar[20]+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.1767766952966368*alphavpar[20]+0.3061862178478971*(alphavpar[19]+alphavpar[18])-0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[27])+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]+0.3061862178478971*alphavpar[19]-0.3061862178478971*alphavpar[18]+0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*((-0.1767766952966368*alphavpar[27])+0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*(alphavpar[21]+alphavpar[20])-0.3061862178478971*alphavpar[19]+0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.03125*(0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])-0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])-0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.03125*(0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[27])+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*alphavpar[19]-0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[27])+0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.1767766952966368*alphavpar[20]-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[27])+0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21])-0.1767766952966368*alphavpar[20]-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*alphavpar[6]-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[27]-0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*alphavpar[21]+0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[27]-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*alphavpar[19]-0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[27])+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]-0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6])-0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[27])-0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])-0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*alphavpar[19]-0.3061862178478971*(alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*(alphavpar[14]+alphavpar[13])+0.1767766952966368*alphavpar[12]-0.3061862178478971*(alphavpar[11]+alphavpar[10])+0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*(alphavpar[3]+alphavpar[2])+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*alphavpar[13]-0.1767766952966368*alphavpar[12]-0.3061862178478971*alphavpar[11]+0.3061862178478971*alphavpar[10]-0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*alphavpar[2]-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[27])-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.1767766952966368*alphavpar[20]-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]-0.1767766952966368*alphavpar[14]+0.1767766952966368*(alphavpar[13]+alphavpar[12])-0.3061862178478971*alphavpar[11]+0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]-0.1767766952966368*alphavpar[3]+0.1767766952966368*(alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[27]+0.3061862178478971*alphavpar[26]-0.3061862178478971*alphavpar[23]-0.1767766952966368*(alphavpar[22]+alphavpar[21])+0.1767766952966368*alphavpar[20]-0.3061862178478971*(alphavpar[19]+alphavpar[18])+0.3061862178478971*alphavpar[17]+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*(alphavpar[13]+alphavpar[12])+0.3061862178478971*alphavpar[11]-0.3061862178478971*(alphavpar[10]+alphavpar[9])-0.1767766952966368*(alphavpar[8]+alphavpar[7])+0.1767766952966368*(alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*(alphavpar[2]+alphavpar[1])+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[27])-0.3061862178478971*alphavpar[26]+0.3061862178478971*alphavpar[23]-0.1767766952966368*alphavpar[22]+0.1767766952966368*alphavpar[21]-0.1767766952966368*alphavpar[20]-0.3061862178478971*alphavpar[19]+0.3061862178478971*alphavpar[18]-0.3061862178478971*alphavpar[17]-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*alphavpar[14]-0.1767766952966368*alphavpar[13]+0.1767766952966368*alphavpar[12]+0.3061862178478971*alphavpar[11]-0.3061862178478971*alphavpar[10]+0.3061862178478971*alphavpar[9]-0.1767766952966368*alphavpar[8]+0.1767766952966368*alphavpar[7]-0.1767766952966368*alphavpar[6]+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*alphavpar[3]-0.1767766952966368*alphavpar[2]+0.1767766952966368*(alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*((-0.1767766952966368*alphavpar[27])-0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*alphavpar[22]-0.1767766952966368*(alphavpar[21]+alphavpar[20])+0.3061862178478971*alphavpar[19]-0.3061862178478971*(alphavpar[18]+alphavpar[17])-0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13])-0.1767766952966368*alphavpar[12]+0.3061862178478971*(alphavpar[11]+alphavpar[10])-0.3061862178478971*alphavpar[9]+0.1767766952966368*alphavpar[8]-0.1767766952966368*(alphavpar[7]+alphavpar[6])+0.1767766952966368*alphavpar[5]+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2])-0.1767766952966368*alphavpar[1]+0.1767766952966368*alphavpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.03125*(0.1767766952966368*alphavpar[27]+0.3061862178478971*(alphavpar[26]+alphavpar[23])+0.1767766952966368*(alphavpar[22]+alphavpar[21]+alphavpar[20])+0.3061862178478971*(alphavpar[19]+alphavpar[18]+alphavpar[17])+0.1767766952966368*alphavpar[16]+0.3061862178478971*alphavpar[15]+0.1767766952966368*(alphavpar[14]+alphavpar[13]+alphavpar[12])+0.3061862178478971*(alphavpar[11]+alphavpar[10]+alphavpar[9])+0.1767766952966368*(alphavpar[8]+alphavpar[7]+alphavpar[6]+alphavpar[5])+0.3061862178478971*alphavpar[4]+0.1767766952966368*(alphavpar[3]+alphavpar[2]+alphavpar[1]+alphavpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
#endif 

  out[1] += 0.3061862178478971*(alphax[16]*f[16]+alphax[9]*f[9]+alphax[8]*f[8]+alphax[7]*f[7]+alphax[6]*f[6]+alphax[4]*f[4]+alphax[3]*f[3]+alphax[2]*f[2]+alphax[1]*f[1]+alphax[0]*f[0]); 
  out[2] += 0.3061862178478971*(alphay[16]*f[16]+alphay[12]*f[12]+alphay[9]*f[9]+alphay[8]*f[8]+alphay[7]*f[7]+alphay[6]*f[6]+alphay[5]*f[5]+alphay[4]*f[4]+alphay[3]*f[3]+alphay[2]*f[2]+alphay[1]*f[1]+alphay[0]*f[0]); 
  out[3] += 0.3061862178478971*(alphaz[1]*f[1]+alphaz[0]*f[0]); 
  out[4] += 0.3061862178478971*(alphavpar[27]*f[27]+alphavpar[26]*f[26]+alphavpar[23]*f[23]+alphavpar[22]*f[22]+alphavpar[21]*f[21]+alphavpar[20]*f[20]+alphavpar[19]*f[19]+alphavpar[18]*f[18]+alphavpar[17]*f[17]+alphavpar[16]*f[16]+alphavpar[15]*f[15]+alphavpar[14]*f[14]+alphavpar[13]*f[13]+alphavpar[12]*f[12]+alphavpar[11]*f[11]+alphavpar[10]*f[10]+alphavpar[9]*f[9]+alphavpar[8]*f[8]+alphavpar[7]*f[7]+alphavpar[6]*f[6]+alphavpar[5]*f[5]+alphavpar[4]*f[4]+alphavpar[3]*f[3]+alphavpar[2]*f[2]+alphavpar[1]*f[1]+alphavpar[0]*f[0]); 
  out[6] += 0.3061862178478971*(alphax[9]*f[17]+(alphay[8]+alphax[7])*f[16]+f[8]*alphay[16]+f[7]*alphax[16]+alphay[5]*f[12]+f[5]*alphay[12]+alphax[4]*f[10]+alphay[4]*f[9]+f[4]*alphay[9]+alphax[3]*f[8]+f[3]*alphax[8]+alphay[3]*f[7]+f[3]*alphay[7]+(alphay[2]+alphax[1])*f[6]+f[2]*alphay[6]+f[1]*alphax[6]+alphax[0]*f[2]+f[0]*alphax[2]+alphay[0]*f[1]+f[0]*alphay[1]); 
  out[7] += 0.3061862178478971*(alphax[9]*f[18]+alphax[6]*f[16]+f[6]*alphax[16]+alphax[4]*f[11]+alphax[2]*f[8]+f[2]*alphax[8]+alphax[1]*f[7]+f[1]*alphax[7]+alphax[0]*f[3]+f[0]*alphax[3]+alphaz[0]*f[1]+f[0]*alphaz[1]); 
  out[8] += 0.3061862178478971*(alphay[12]*f[21]+alphay[9]*f[18]+alphay[6]*f[16]+f[6]*alphay[16]+alphay[5]*f[14]+alphay[4]*f[11]+alphay[2]*f[8]+f[2]*alphay[8]+alphay[1]*f[7]+f[1]*alphay[7]+alphaz[1]*f[6]+alphay[0]*f[3]+f[0]*alphay[3]+alphaz[0]*f[2]); 
  out[9] += 0.3061862178478971*(alphavpar[22]*f[27]+f[22]*alphavpar[27]+(alphavpar[19]+alphax[16])*f[26]+f[19]*alphavpar[26]+alphavpar[15]*f[23]+f[15]*alphavpar[23]+alphavpar[14]*f[21]+f[14]*alphavpar[21]+alphavpar[13]*f[20]+f[13]*alphavpar[20]+alphax[8]*f[19]+(alphavpar[11]+alphax[7])*f[18]+f[11]*alphavpar[18]+(alphavpar[10]+alphax[6])*f[17]+f[10]*alphavpar[17]+alphavpar[8]*f[16]+f[8]*alphavpar[16]+alphavpar[5]*f[12]+f[5]*alphavpar[12]+alphax[3]*f[11]+alphax[2]*f[10]+(alphavpar[4]+alphax[1])*f[9]+f[1]*alphax[9]+f[4]*alphavpar[9]+alphavpar[3]*f[7]+f[3]*alphavpar[7]+alphavpar[2]*f[6]+f[2]*alphavpar[6]+alphax[0]*f[4]+f[0]*alphax[4]+alphavpar[0]*f[1]+f[0]*alphavpar[1]); 
  out[10] += 0.3061862178478971*(alphavpar[23]*f[28]+alphavpar[21]*f[27]+f[21]*alphavpar[27]+(alphavpar[18]+alphay[16])*f[26]+f[18]*alphavpar[26]+alphavpar[15]*f[24]+alphay[12]*f[23]+alphavpar[14]*f[22]+f[14]*alphavpar[22]+alphavpar[12]*f[20]+f[12]*alphavpar[20]+(alphavpar[11]+alphay[8])*f[19]+f[11]*alphavpar[19]+alphay[7]*f[18]+(alphavpar[9]+alphay[6])*f[17]+f[9]*alphavpar[17]+alphavpar[7]*f[16]+f[7]*alphavpar[16]+alphay[5]*f[15]+alphavpar[5]*f[13]+f[5]*alphavpar[13]+alphay[3]*f[11]+(alphavpar[4]+alphay[2])*f[10]+f[4]*alphavpar[10]+alphay[1]*f[9]+f[1]*alphay[9]+alphavpar[3]*f[8]+f[3]*alphavpar[8]+alphavpar[1]*f[6]+f[1]*alphavpar[6]+alphay[0]*f[4]+f[0]*alphay[4]+alphavpar[0]*f[2]+f[0]*alphavpar[2]); 
  out[11] += 0.3061862178478971*(alphavpar[23]*f[29]+alphavpar[20]*f[27]+f[20]*alphavpar[27]+alphavpar[17]*f[26]+f[17]*alphavpar[26]+alphavpar[15]*f[25]+alphavpar[13]*f[22]+f[13]*alphavpar[22]+alphavpar[12]*f[21]+f[12]*alphavpar[21]+alphavpar[10]*f[19]+f[10]*alphavpar[19]+alphavpar[9]*f[18]+f[9]*alphavpar[18]+alphavpar[6]*f[16]+f[6]*alphavpar[16]+alphavpar[5]*f[14]+f[5]*alphavpar[14]+alphavpar[4]*f[11]+f[4]*alphavpar[11]+alphaz[1]*f[9]+alphavpar[2]*f[8]+f[2]*alphavpar[8]+alphavpar[1]*f[7]+f[1]*alphavpar[7]+alphaz[0]*f[4]+alphavpar[0]*f[3]+f[0]*alphavpar[3]); 
  out[12] += 0.3061862178478971*(alphax[16]*f[27]+alphax[9]*f[23]+alphax[8]*f[22]+alphax[7]*f[21]+alphax[6]*f[20]+alphax[4]*f[15]+alphax[3]*f[14]+alphax[2]*f[13]+alphax[1]*f[12]+alphax[0]*f[5]); 
  out[13] += 0.3061862178478971*(alphay[16]*f[27]+alphay[9]*f[23]+alphay[8]*f[22]+alphay[7]*f[21]+alphay[6]*f[20]+alphay[4]*f[15]+alphay[3]*f[14]+alphay[2]*f[13]+alphay[1]*f[12]+f[1]*alphay[12]+alphay[0]*f[5]+f[0]*alphay[5]); 
  out[14] += 0.3061862178478971*(alphaz[1]*f[12]+alphaz[0]*f[5]); 
  out[15] += 0.3061862178478971*(alphavpar[26]*f[31]+alphavpar[19]*f[30]+alphavpar[18]*f[29]+alphavpar[17]*f[28]+alphavpar[16]*f[27]+f[16]*alphavpar[27]+alphavpar[11]*f[25]+alphavpar[10]*f[24]+alphavpar[9]*f[23]+f[9]*alphavpar[23]+alphavpar[8]*f[22]+f[8]*alphavpar[22]+alphavpar[7]*f[21]+f[7]*alphavpar[21]+alphavpar[6]*f[20]+f[6]*alphavpar[20]+alphavpar[4]*f[15]+f[4]*alphavpar[15]+alphavpar[3]*f[14]+f[3]*alphavpar[14]+alphavpar[2]*f[13]+f[2]*alphavpar[13]+alphavpar[1]*f[12]+f[1]*alphavpar[12]+alphavpar[0]*f[5]+f[0]*alphavpar[5]); 
  out[16] += 0.3061862178478971*(alphax[9]*f[26]+alphay[5]*f[21]+alphax[4]*f[19]+alphay[4]*f[18]+(alphay[2]+alphax[1])*f[16]+f[2]*alphay[16]+f[1]*alphax[16]+alphay[12]*f[14]+alphay[9]*f[11]+(alphay[6]+alphax[0])*f[8]+f[6]*alphay[8]+f[0]*alphax[8]+(alphax[6]+alphay[0])*f[7]+f[0]*alphay[7]+f[6]*(alphax[7]+alphaz[0])+(alphax[2]+alphay[1])*f[3]+f[1]*alphay[3]+f[2]*(alphax[3]+alphaz[1])); 
  out[17] += 0.3061862178478971*(alphavpar[15]*f[28]+alphavpar[14]*f[27]+f[14]*alphavpar[27]+(alphavpar[11]+alphay[8]+alphax[7])*f[26]+f[11]*alphavpar[26]+alphavpar[23]*f[24]+alphay[5]*f[23]+alphavpar[21]*f[22]+f[21]*alphavpar[22]+alphavpar[5]*f[20]+f[5]*alphavpar[20]+(alphavpar[18]+alphay[16]+alphax[3])*f[19]+f[18]*(alphavpar[19]+alphax[16]+alphay[3])+(alphavpar[4]+alphay[2]+alphax[1])*f[17]+f[4]*alphavpar[17]+alphavpar[3]*f[16]+f[3]*alphavpar[16]+alphay[12]*f[15]+alphavpar[12]*f[13]+f[12]*alphavpar[13]+(alphax[8]+alphay[7])*f[11]+(alphavpar[9]+alphay[6]+alphax[0])*f[10]+f[9]*(alphavpar[10]+alphax[6]+alphay[0])+f[0]*alphay[9]+f[6]*alphax[9]+alphavpar[7]*f[8]+f[7]*alphavpar[8]+alphavpar[0]*f[6]+f[0]*alphavpar[6]+(alphax[2]+alphay[1])*f[4]+f[1]*alphay[4]+f[2]*(alphax[4]+alphavpar[1])+f[1]*alphavpar[2]); 
  out[18] += 0.3061862178478971*(alphavpar[15]*f[29]+alphavpar[13]*f[27]+f[13]*alphavpar[27]+(alphavpar[10]+alphax[6])*f[26]+f[10]*alphavpar[26]+alphavpar[23]*f[25]+alphavpar[20]*f[22]+f[20]*alphavpar[22]+alphavpar[5]*f[21]+f[5]*alphavpar[21]+(alphavpar[17]+alphax[2])*f[19]+f[17]*alphavpar[19]+(alphavpar[4]+alphax[1])*f[18]+f[4]*alphavpar[18]+alphax[16]*f[17]+alphavpar[2]*f[16]+f[2]*alphavpar[16]+alphavpar[12]*f[14]+f[12]*alphavpar[14]+(alphavpar[9]+alphax[0])*f[11]+f[9]*alphavpar[11]+alphax[8]*f[10]+(alphax[7]+alphaz[0])*f[9]+f[7]*alphax[9]+alphavpar[6]*f[8]+f[6]*alphavpar[8]+alphavpar[0]*f[7]+f[0]*alphavpar[7]+(alphax[3]+alphaz[1])*f[4]+f[3]*(alphax[4]+alphavpar[1])+f[1]*alphavpar[3]); 
  out[19] += 0.3061862178478971*(alphavpar[23]*f[31]+alphavpar[15]*f[30]+alphay[12]*f[29]+alphavpar[12]*f[27]+f[12]*alphavpar[27]+(alphavpar[9]+alphay[6])*f[26]+f[9]*alphavpar[26]+alphay[5]*f[25]+alphavpar[5]*f[22]+f[5]*alphavpar[22]+alphavpar[20]*f[21]+f[20]*alphavpar[21]+(alphavpar[4]+alphay[2])*f[19]+f[4]*alphavpar[19]+(alphavpar[17]+alphay[1])*f[18]+f[17]*(alphavpar[18]+alphay[16]+alphaz[1])+alphavpar[1]*f[16]+f[1]*alphavpar[16]+alphavpar[13]*f[14]+f[13]*alphavpar[14]+(alphavpar[10]+alphay[0])*f[11]+f[10]*(alphavpar[11]+alphay[8]+alphaz[0])+alphay[7]*f[9]+f[7]*alphay[9]+alphavpar[0]*f[8]+f[0]*alphavpar[8]+alphavpar[6]*f[7]+f[6]*alphavpar[7]+alphay[3]*f[4]+f[3]*(alphay[4]+alphavpar[2])+f[2]*alphavpar[3]); 
  out[20] += 0.3061862178478971*(alphax[9]*f[28]+(alphay[8]+alphax[7])*f[27]+alphax[4]*f[24]+alphay[4]*f[23]+(alphay[16]+alphax[3])*f[22]+(alphax[16]+alphay[3])*f[21]+(alphay[2]+alphax[1])*f[20]+alphay[9]*f[15]+(alphax[8]+alphay[7])*f[14]+(alphay[6]+alphax[0])*f[13]+(alphax[6]+alphay[0])*f[12]+f[0]*alphay[12]+(alphax[2]+alphay[1])*f[5]+f[1]*alphay[5]); 
  out[21] += 0.3061862178478971*(alphax[9]*f[29]+alphax[6]*f[27]+alphax[4]*f[25]+alphax[2]*f[22]+alphax[1]*f[21]+alphax[16]*f[20]+alphax[0]*f[14]+alphax[8]*f[13]+(alphax[7]+alphaz[0])*f[12]+(alphax[3]+alphaz[1])*f[5]); 
  out[22] += 0.3061862178478971*(alphay[9]*f[29]+alphay[6]*f[27]+alphay[4]*f[25]+alphay[2]*f[22]+alphay[1]*f[21]+(alphay[16]+alphaz[1])*f[20]+alphay[0]*f[14]+(alphay[8]+alphaz[0])*f[13]+alphay[7]*f[12]+f[7]*alphay[12]+alphay[3]*f[5]+f[3]*alphay[5]); 
  out[23] += 0.3061862178478971*((alphavpar[19]+alphax[16])*f[31]+(alphavpar[26]+alphax[8])*f[30]+(alphavpar[11]+alphax[7])*f[29]+(alphavpar[10]+alphax[6])*f[28]+alphavpar[8]*f[27]+f[8]*alphavpar[27]+(alphavpar[18]+alphax[3])*f[25]+(alphavpar[17]+alphax[2])*f[24]+(alphavpar[4]+alphax[1])*f[23]+f[4]*alphavpar[23]+alphavpar[16]*f[22]+f[16]*alphavpar[22]+alphavpar[3]*f[21]+f[3]*alphavpar[21]+alphavpar[2]*f[20]+f[2]*alphavpar[20]+(alphavpar[9]+alphax[0])*f[15]+f[9]*alphavpar[15]+alphavpar[7]*f[14]+f[7]*alphavpar[14]+alphavpar[6]*f[13]+f[6]*alphavpar[13]+(alphax[9]+alphavpar[0])*f[12]+f[0]*alphavpar[12]+(alphax[4]+alphavpar[1])*f[5]+f[1]*alphavpar[5]); 
  out[24] += 0.3061862178478971*((alphavpar[18]+alphay[16])*f[31]+(alphavpar[11]+alphay[8])*f[30]+(alphavpar[26]+alphay[7])*f[29]+(alphavpar[9]+alphay[6])*f[28]+alphavpar[7]*f[27]+f[7]*alphavpar[27]+(alphavpar[19]+alphay[3])*f[25]+(alphavpar[4]+alphay[2])*f[24]+(alphavpar[17]+alphay[1])*f[23]+f[17]*alphavpar[23]+alphavpar[3]*f[22]+f[3]*alphavpar[22]+alphavpar[16]*f[21]+f[16]*alphavpar[21]+alphavpar[1]*f[20]+f[1]*alphavpar[20]+(alphavpar[10]+alphay[0])*f[15]+f[10]*alphavpar[15]+alphavpar[8]*f[14]+f[8]*alphavpar[14]+alphavpar[0]*f[13]+f[0]*alphavpar[13]+(alphay[9]+alphavpar[6])*f[12]+f[9]*alphay[12]+f[6]*alphavpar[12]+(alphay[4]+alphavpar[2])*f[5]+f[4]*alphay[5]+f[2]*alphavpar[5]); 
  out[25] += 0.3061862178478971*(alphavpar[17]*f[31]+alphavpar[10]*f[30]+alphavpar[9]*f[29]+alphavpar[26]*f[28]+alphavpar[6]*f[27]+f[6]*alphavpar[27]+alphavpar[4]*f[25]+alphavpar[19]*f[24]+(alphavpar[18]+alphaz[1])*f[23]+f[18]*alphavpar[23]+alphavpar[2]*f[22]+f[2]*alphavpar[22]+alphavpar[1]*f[21]+f[1]*alphavpar[21]+alphavpar[16]*f[20]+f[16]*alphavpar[20]+(alphavpar[11]+alphaz[0])*f[15]+f[11]*alphavpar[15]+alphavpar[0]*f[14]+f[0]*alphavpar[14]+alphavpar[8]*f[13]+f[8]*alphavpar[13]+alphavpar[7]*f[12]+f[7]*alphavpar[12]+alphavpar[3]*f[5]+f[3]*alphavpar[5]); 
  out[26] += 0.3061862178478971*(alphavpar[15]*f[31]+alphavpar[23]*f[30]+alphay[5]*f[29]+alphavpar[5]*f[27]+f[5]*alphavpar[27]+(alphavpar[4]+alphay[2]+alphax[1])*f[26]+f[4]*alphavpar[26]+alphay[12]*f[25]+alphavpar[12]*f[22]+f[12]*alphavpar[22]+alphavpar[13]*f[21]+f[13]*alphavpar[21]+alphavpar[14]*f[20]+f[14]*alphavpar[20]+(alphavpar[9]+alphay[6]+alphax[0])*f[19]+f[9]*alphavpar[19]+(alphavpar[10]+alphax[6]+alphay[0])*f[18]+f[10]*alphavpar[18]+(alphavpar[11]+alphay[8]+alphax[7]+alphaz[0])*f[17]+f[11]*alphavpar[17]+(alphax[9]+alphavpar[0])*f[16]+f[10]*alphay[16]+f[9]*alphax[16]+f[0]*alphavpar[16]+(alphax[2]+alphay[1])*f[11]+(alphax[3]+alphaz[1])*f[10]+alphay[3]*f[9]+f[3]*alphay[9]+(alphax[4]+alphavpar[1])*f[8]+f[4]*alphax[8]+f[1]*alphavpar[8]+(alphay[4]+alphavpar[2])*f[7]+f[4]*alphay[7]+f[2]*alphavpar[7]+alphavpar[3]*f[6]+f[3]*alphavpar[6]); 
  out[27] += 0.3061862178478971*(alphax[9]*f[31]+alphax[4]*f[30]+alphay[4]*f[29]+(alphay[2]+alphax[1])*f[27]+alphay[9]*f[25]+(alphay[6]+alphax[0])*f[22]+(alphax[6]+alphay[0])*f[21]+(alphay[8]+alphax[7]+alphaz[0])*f[20]+f[13]*alphay[16]+f[12]*alphax[16]+(alphax[2]+alphay[1])*f[14]+(alphax[3]+alphaz[1])*f[13]+alphay[3]*f[12]+f[3]*alphay[12]+f[5]*alphax[8]+alphay[5]*f[7]+f[5]*alphay[7]); 
  out[28] += 0.3061862178478971*((alphavpar[11]+alphay[8]+alphax[7])*f[31]+(alphavpar[18]+alphay[16]+alphax[3])*f[30]+(alphavpar[19]+alphax[16]+alphay[3])*f[29]+(alphavpar[4]+alphay[2]+alphax[1])*f[28]+alphavpar[3]*f[27]+f[3]*alphavpar[27]+f[25]*(alphavpar[26]+alphax[8]+alphay[7])+(alphavpar[9]+alphay[6]+alphax[0])*f[24]+(alphavpar[10]+alphax[6]+alphay[0])*f[23]+f[10]*alphavpar[23]+alphavpar[7]*f[22]+f[7]*alphavpar[22]+alphavpar[8]*f[21]+f[8]*alphavpar[21]+(alphax[9]+alphavpar[0])*f[20]+f[0]*alphavpar[20]+alphavpar[15]*f[17]+f[15]*alphavpar[17]+alphavpar[14]*f[16]+f[14]*alphavpar[16]+(alphax[2]+alphay[1])*f[15]+(alphax[4]+alphavpar[1])*f[13]+f[1]*alphavpar[13]+(alphay[4]+alphavpar[2])*f[12]+f[4]*alphay[12]+f[2]*alphavpar[12]+alphay[5]*f[9]+f[5]*alphay[9]+alphavpar[5]*f[6]+f[5]*alphavpar[6]); 
  out[29] += 0.3061862178478971*((alphavpar[10]+alphax[6])*f[31]+(alphavpar[17]+alphax[2])*f[30]+(alphavpar[4]+alphax[1])*f[29]+(alphavpar[19]+alphax[16])*f[28]+alphavpar[2]*f[27]+f[2]*alphavpar[27]+f[24]*alphavpar[26]+(alphavpar[9]+alphax[0])*f[25]+alphax[8]*f[24]+(alphavpar[11]+alphax[7]+alphaz[0])*f[23]+f[11]*alphavpar[23]+alphavpar[6]*f[22]+f[6]*alphavpar[22]+(alphax[9]+alphavpar[0])*f[21]+f[0]*alphavpar[21]+alphavpar[8]*f[20]+f[8]*alphavpar[20]+alphavpar[15]*f[18]+f[15]*alphavpar[18]+alphavpar[13]*f[16]+f[13]*alphavpar[16]+(alphax[3]+alphaz[1])*f[15]+(alphax[4]+alphavpar[1])*f[14]+f[1]*alphavpar[14]+alphavpar[3]*f[12]+f[3]*alphavpar[12]+alphavpar[5]*f[7]+f[5]*alphavpar[7]); 
  out[30] += 0.3061862178478971*((alphavpar[9]+alphay[6])*f[31]+(alphavpar[4]+alphay[2])*f[30]+(alphavpar[17]+alphay[1])*f[29]+(alphavpar[18]+alphay[16]+alphaz[1])*f[28]+alphavpar[1]*f[27]+f[1]*alphavpar[27]+alphavpar[23]*f[26]+f[23]*alphavpar[26]+(alphavpar[10]+alphay[0])*f[25]+(alphavpar[11]+alphay[8]+alphaz[0])*f[24]+alphay[7]*f[23]+alphavpar[0]*f[22]+f[0]*alphavpar[22]+(alphay[9]+alphavpar[6])*f[21]+f[6]*alphavpar[21]+alphavpar[7]*f[20]+f[7]*alphavpar[20]+alphavpar[15]*f[19]+f[15]*alphavpar[19]+alphay[12]*f[18]+alphavpar[12]*f[16]+f[12]*alphavpar[16]+alphay[3]*f[15]+(alphay[4]+alphavpar[2])*f[14]+f[2]*alphavpar[14]+alphavpar[3]*f[13]+f[3]*alphavpar[13]+alphay[5]*f[11]+alphavpar[5]*f[8]+f[5]*alphavpar[8]); 
  out[31] += 0.3061862178478971*((alphavpar[4]+alphay[2]+alphax[1])*f[31]+(alphavpar[9]+alphay[6]+alphax[0])*f[30]+(alphavpar[10]+alphax[6]+alphay[0])*f[29]+(alphavpar[11]+alphay[8]+alphax[7]+alphaz[0])*f[28]+(alphax[9]+alphavpar[0])*f[27]+f[0]*alphavpar[27]+alphavpar[15]*f[26]+f[15]*alphavpar[26]+(alphavpar[17]+alphax[2]+alphay[1])*f[25]+(alphavpar[18]+alphay[16]+alphax[3]+alphaz[1])*f[24]+(alphavpar[19]+alphax[16]+alphay[3])*f[23]+f[19]*alphavpar[23]+(alphax[4]+alphavpar[1])*f[22]+f[1]*alphavpar[22]+(alphay[4]+alphavpar[2])*f[21]+f[2]*alphavpar[21]+alphavpar[3]*f[20]+f[3]*alphavpar[20]+alphay[5]*f[18]+alphavpar[5]*f[16]+f[5]*alphavpar[16]+(alphax[8]+alphay[7])*f[15]+(alphay[9]+alphavpar[6])*f[14]+f[6]*alphavpar[14]+alphavpar[7]*f[13]+f[7]*alphavpar[13]+alphavpar[8]*f[12]+f[11]*alphay[12]+f[8]*alphavpar[12]); 
  return cflFreq; 
} 
double EmGyrokineticSimpleHelicalStep2Vol3x2vSerP1(const double q_, const double m_, const double *w, const double *dxv, const double *dApardt, const double *f, double *out) 
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
