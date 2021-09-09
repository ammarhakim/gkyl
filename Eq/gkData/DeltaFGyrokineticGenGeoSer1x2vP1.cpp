#include <GyrokineticModDecl.h> 
double DeltaFGyrokineticGenGeoVol1x2vSerP1_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0, const double *f1, double *out) 
{ 
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

  double hamil0[8]; 
  hamil0[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*bmag[0]*wmu)+2.0*m_))/rdvpar2Sq; 
  hamil0[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil0[3] = (1.154700538379252*bmag[0])/rdmu2; 

  double hamil1[8]; 
  hamil1[0] = 2.0*phi[0]*q_; 
  hamil1[1] = 2.0*phi[1]*q_; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[8]; 
  alpha0x[0] = (0.6123724356957944*BstarZdBmag[0]*hamil0[2]*rdvpar2*rdx2)/m_; 
  double alpha1x[8]; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0441941738241592*alpha0x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha0x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha0x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha0x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0441941738241592*alpha0x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha0x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha0x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha0x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[8]; 
  double alpha1vpar[8]; 
  alpha1vpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamil1[1]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.6123724356957944*alpha0x[0]*f1[0]; 
  out[2] += 0.6123724356957944*alpha1vpar[0]*f0[0]; 
  out[4] += 0.6123724356957944*(alpha0x[0]*f1[2]+alpha1vpar[0]*f0[1]); 
  out[5] += 0.6123724356957944*alpha0x[0]*f1[3]; 
  out[6] += 0.6123724356957944*alpha1vpar[0]*f0[3]; 
  out[7] += 0.6123724356957944*(alpha0x[0]*f1[6]+alpha1vpar[0]*f0[5]); 
  out[2] += 0.6123724356957944*alpha1vpar[0]*f1[0]; 
  out[4] += 0.6123724356957944*alpha1vpar[0]*f1[1]; 
  out[6] += 0.6123724356957944*alpha1vpar[0]*f1[3]; 
  out[7] += 0.6123724356957944*alpha1vpar[0]*f1[5]; 
  return cflFreq; 
} 
double DeltaFGyrokineticGenGeoVol1x2vSerP1_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0, const double *f1, double *out) 
{ 
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

  double hamil0[8]; 
  hamil0[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*bmag[0]*wmu)+2.0*m_))/rdvpar2Sq; 
  hamil0[1] = 2.0*bmag[1]*wmu; 
  hamil0[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil0[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil0[5] = (1.154700538379252*bmag[1])/rdmu2; 

  double hamil1[8]; 
  hamil1[0] = 2.0*phi[0]*q_; 
  hamil1[1] = 2.0*phi[1]*q_; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmag[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2*wvpar+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[8]; 
  alpha0x[0] = (0.6123724356957944*BstarZdBmag[0]*hamil0[2]*rdvpar2*rdx2)/m_; 
  alpha0x[1] = (0.6123724356957944*BstarZdBmag[1]*hamil0[2]*rdvpar2*rdx2)/m_; 
  alpha0x[2] = (0.6123724356957944*BstarZdBmag[2]*hamil0[2]*rdvpar2*rdx2)/m_; 
  alpha0x[4] = (0.6123724356957944*hamil0[2]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
  double alpha1x[8]; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.6123724356957944*alpha0x[4]-0.3535533905932737*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*alpha0x[4])+0.3535533905932737*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*alpha0x[4]-0.3535533905932737*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*alpha0x[4])+0.3535533905932737*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*((-0.6123724356957944*alpha0x[4])-0.3535533905932737*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alpha0x[4]+0.3535533905932737*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*alpha0x[4])-0.3535533905932737*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alpha0x[4]+0.3535533905932737*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[8]; 
  alpha0vpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamil0[1]*rdvpar2*rdx2)/m_; 
  alpha0vpar[1] = -(0.6123724356957944*BstarZdBmag[1]*hamil0[1]*rdvpar2*rdx2)/m_; 
  alpha0vpar[2] = -(0.6123724356957944*hamil0[1]*BstarZdBmag[2]*rdvpar2*rdx2)/m_; 
  alpha0vpar[3] = -(0.6123724356957944*BstarZdBmag[0]*hamil0[5]*rdvpar2*rdx2)/m_; 
  alpha0vpar[4] = -(0.6123724356957944*hamil0[1]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
  alpha0vpar[5] = -(0.6123724356957944*BstarZdBmag[1]*hamil0[5]*rdvpar2*rdx2)/m_; 
  alpha0vpar[6] = -(0.6123724356957944*BstarZdBmag[2]*hamil0[5]*rdvpar2*rdx2)/m_; 
  alpha0vpar[7] = -(0.6123724356957944*BstarZdBmag[4]*hamil0[5]*rdvpar2*rdx2)/m_; 
  double alpha1vpar[8]; 
  alpha1vpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamil1[1]*rdvpar2*rdx2)/m_; 
  alpha1vpar[1] = -(0.6123724356957944*BstarZdBmag[1]*hamil1[1]*rdvpar2*rdx2)/m_; 
  alpha1vpar[2] = -(0.6123724356957944*hamil1[1]*BstarZdBmag[2]*rdvpar2*rdx2)/m_; 
  alpha1vpar[4] = -(0.6123724356957944*hamil1[1]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*((-0.6123724356957944*alpha0vpar[7])+0.6123724356957944*alpha0vpar[6]+0.3535533905932737*alpha0vpar[5]+0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])-0.3535533905932737*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*(alpha0vpar[7]+alpha0vpar[6])-0.3535533905932737*alpha0vpar[5]-0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])-0.3535533905932737*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1]+alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*alpha0vpar[7]-0.6123724356957944*alpha0vpar[6]-0.3535533905932737*alpha0vpar[5]+0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])+0.3535533905932737*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*(alpha0vpar[7]+alpha0vpar[6]))+0.3535533905932737*alpha0vpar[5]-0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])+0.3535533905932737*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1]+alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.6123724356957944*alpha0vpar[7]-0.6123724356957944*alpha0vpar[6]+0.3535533905932737*alpha0vpar[5]-0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])-0.3535533905932737*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*(alpha0vpar[7]+alpha0vpar[6]))-0.3535533905932737*alpha0vpar[5]+0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])-0.3535533905932737*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1]+alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*alpha0vpar[7])+0.6123724356957944*alpha0vpar[6]-0.3535533905932737*alpha0vpar[5]-0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])+0.3535533905932737*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*(alpha0vpar[7]+alpha0vpar[6])+0.3535533905932737*alpha0vpar[5]+0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])+0.3535533905932737*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1]+alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.6123724356957944*(alpha0x[4]*f1[4]+alpha0x[2]*f1[2]+alpha0x[1]*f1[1]+alpha0x[0]*f1[0]); 
  out[2] += 0.6123724356957944*(alpha0vpar[7]*f1[7]+alpha0vpar[6]*f1[6]+alpha0vpar[5]*f1[5]+alpha0vpar[4]*f1[4]+alpha1vpar[4]*f0[4]+alpha0vpar[3]*f1[3]+alpha0vpar[2]*f1[2]+alpha1vpar[2]*f0[2]+alpha0vpar[1]*f1[1]+alpha1vpar[1]*f0[1]+alpha0vpar[0]*f1[0]+alpha1vpar[0]*f0[0]); 
  out[4] += 0.6123724356957944*(alpha0vpar[6]*f1[7]+f1[6]*alpha0vpar[7]+alpha0vpar[3]*f1[5]+f1[3]*alpha0vpar[5]+(alpha0vpar[2]+alpha0x[1])*f1[4]+alpha1vpar[2]*f0[4]+f0[2]*alpha1vpar[4]+f1[1]*alpha0x[4]+f1[2]*(alpha0vpar[4]+alpha0x[0])+f1[0]*alpha0x[2]+alpha0vpar[0]*f1[1]+alpha1vpar[0]*f0[1]+f0[0]*alpha1vpar[1]+f1[0]*alpha0vpar[1]); 
  out[5] += 0.6123724356957944*(alpha0x[4]*f1[7]+alpha0x[2]*f1[6]+alpha0x[1]*f1[5]+alpha0x[0]*f1[3]); 
  out[6] += 0.6123724356957944*(alpha0vpar[4]*f1[7]+alpha1vpar[4]*f0[7]+f1[4]*alpha0vpar[7]+alpha0vpar[2]*f1[6]+alpha1vpar[2]*f0[6]+f1[2]*alpha0vpar[6]+alpha0vpar[1]*f1[5]+alpha1vpar[1]*f0[5]+f1[1]*alpha0vpar[5]+alpha0vpar[0]*f1[3]+alpha1vpar[0]*f0[3]+f1[0]*alpha0vpar[3]); 
  out[7] += 0.6123724356957944*((alpha0vpar[2]+alpha0x[1])*f1[7]+alpha1vpar[2]*f0[7]+f1[2]*alpha0vpar[7]+(alpha0vpar[4]+alpha0x[0])*f1[6]+alpha1vpar[4]*f0[6]+f1[4]*alpha0vpar[6]+(alpha0x[4]+alpha0vpar[0])*f1[5]+alpha1vpar[0]*f0[5]+f1[0]*alpha0vpar[5]+(alpha0x[2]+alpha0vpar[1])*f1[3]+alpha1vpar[1]*f0[3]+f1[1]*alpha0vpar[3]); 
  out[2] += 0.6123724356957944*(alpha1vpar[4]*f1[4]+alpha1vpar[2]*f1[2]+alpha1vpar[1]*f1[1]+alpha1vpar[0]*f1[0]); 
  out[4] += 0.6123724356957944*(alpha1vpar[2]*f1[4]+f1[2]*alpha1vpar[4]+alpha1vpar[0]*f1[1]+f1[0]*alpha1vpar[1]); 
  out[6] += 0.6123724356957944*(alpha1vpar[4]*f1[7]+alpha1vpar[2]*f1[6]+alpha1vpar[1]*f1[5]+alpha1vpar[0]*f1[3]); 
  out[7] += 0.6123724356957944*(alpha1vpar[2]*f1[7]+alpha1vpar[4]*f1[6]+alpha1vpar[0]*f1[5]+alpha1vpar[1]*f1[3]); 
  return cflFreq; 
} 
double DeltaFGyrokineticGenGeoVol1x2vSerP1_Bvarsz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0, const double *f1, double *out) 
{ 
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

  double hamil0[8]; 
  hamil0[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*bmag[0]*wmu)+2.0*m_))/rdvpar2Sq; 
  hamil0[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil0[3] = (1.154700538379252*bmag[0])/rdmu2; 

  double hamil1[8]; 
  hamil1[0] = 2.0*phi[0]*q_; 
  hamil1[1] = 2.0*phi[1]*q_; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[8]; 
  alpha0x[0] = (0.6123724356957944*BstarZdBmag[0]*hamil0[2]*rdvpar2*rdx2)/m_; 
  double alpha1x[8]; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0441941738241592*alpha0x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha0x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha0x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha0x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0441941738241592*alpha0x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha0x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha0x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha0x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[8]; 
  double alpha1vpar[8]; 
  alpha1vpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamil1[1]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.6123724356957944*alpha0x[0]*f1[0]; 
  out[2] += 0.6123724356957944*alpha1vpar[0]*f0[0]; 
  out[4] += 0.6123724356957944*(alpha0x[0]*f1[2]+alpha1vpar[0]*f0[1]); 
  out[5] += 0.6123724356957944*alpha0x[0]*f1[3]; 
  out[6] += 0.6123724356957944*alpha1vpar[0]*f0[3]; 
  out[7] += 0.6123724356957944*(alpha0x[0]*f1[6]+alpha1vpar[0]*f0[5]); 
  out[2] += 0.6123724356957944*alpha1vpar[0]*f1[0]; 
  out[4] += 0.6123724356957944*alpha1vpar[0]*f1[1]; 
  out[6] += 0.6123724356957944*alpha1vpar[0]*f1[3]; 
  out[7] += 0.6123724356957944*alpha1vpar[0]*f1[5]; 
  return cflFreq; 
} 
double DeltaFGyrokineticGenGeoVol1x2vSerP1_Bvarsxz(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0, const double *f1, double *out) 
{ 
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

  double hamil0[8]; 
  hamil0[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*bmag[0]*wmu)+2.0*m_))/rdvpar2Sq; 
  hamil0[1] = 2.0*bmag[1]*wmu; 
  hamil0[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil0[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil0[5] = (1.154700538379252*bmag[1])/rdmu2; 

  double hamil1[8]; 
  hamil1[0] = 2.0*phi[0]*q_; 
  hamil1[1] = 2.0*phi[1]*q_; 

  double BstarZdBmag[8]; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2*wvpar+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmag[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2*wvpar+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[8]; 
  alpha0x[0] = (0.6123724356957944*BstarZdBmag[0]*hamil0[2]*rdvpar2*rdx2)/m_; 
  alpha0x[1] = (0.6123724356957944*BstarZdBmag[1]*hamil0[2]*rdvpar2*rdx2)/m_; 
  alpha0x[2] = (0.6123724356957944*BstarZdBmag[2]*hamil0[2]*rdvpar2*rdx2)/m_; 
  alpha0x[4] = (0.6123724356957944*hamil0[2]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
  double alpha1x[8]; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.6123724356957944*alpha0x[4]-0.3535533905932737*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*alpha0x[4])+0.3535533905932737*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*alpha0x[4]-0.3535533905932737*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*alpha0x[4])+0.3535533905932737*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*((-0.6123724356957944*alpha0x[4])-0.3535533905932737*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alpha0x[4]+0.3535533905932737*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*alpha0x[4])-0.3535533905932737*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*alpha0x[4]+0.3535533905932737*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[8]; 
  alpha0vpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamil0[1]*rdvpar2*rdx2)/m_; 
  alpha0vpar[1] = -(0.6123724356957944*BstarZdBmag[1]*hamil0[1]*rdvpar2*rdx2)/m_; 
  alpha0vpar[2] = -(0.6123724356957944*hamil0[1]*BstarZdBmag[2]*rdvpar2*rdx2)/m_; 
  alpha0vpar[3] = -(0.6123724356957944*BstarZdBmag[0]*hamil0[5]*rdvpar2*rdx2)/m_; 
  alpha0vpar[4] = -(0.6123724356957944*hamil0[1]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
  alpha0vpar[5] = -(0.6123724356957944*BstarZdBmag[1]*hamil0[5]*rdvpar2*rdx2)/m_; 
  alpha0vpar[6] = -(0.6123724356957944*BstarZdBmag[2]*hamil0[5]*rdvpar2*rdx2)/m_; 
  alpha0vpar[7] = -(0.6123724356957944*BstarZdBmag[4]*hamil0[5]*rdvpar2*rdx2)/m_; 
  double alpha1vpar[8]; 
  alpha1vpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamil1[1]*rdvpar2*rdx2)/m_; 
  alpha1vpar[1] = -(0.6123724356957944*BstarZdBmag[1]*hamil1[1]*rdvpar2*rdx2)/m_; 
  alpha1vpar[2] = -(0.6123724356957944*hamil1[1]*BstarZdBmag[2]*rdvpar2*rdx2)/m_; 
  alpha1vpar[4] = -(0.6123724356957944*hamil1[1]*BstarZdBmag[4]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*((-0.6123724356957944*alpha0vpar[7])+0.6123724356957944*alpha0vpar[6]+0.3535533905932737*alpha0vpar[5]+0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])-0.3535533905932737*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*(alpha0vpar[7]+alpha0vpar[6])-0.3535533905932737*alpha0vpar[5]-0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])-0.3535533905932737*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1]+alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6123724356957944*alpha0vpar[7]-0.6123724356957944*alpha0vpar[6]-0.3535533905932737*alpha0vpar[5]+0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])+0.3535533905932737*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.6123724356957944*(alpha0vpar[7]+alpha0vpar[6]))+0.3535533905932737*alpha0vpar[5]-0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])+0.3535533905932737*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1]+alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.6123724356957944*alpha0vpar[7]-0.6123724356957944*alpha0vpar[6]+0.3535533905932737*alpha0vpar[5]-0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])-0.3535533905932737*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*(alpha0vpar[7]+alpha0vpar[6]))-0.3535533905932737*alpha0vpar[5]+0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])-0.3535533905932737*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1]+alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6123724356957944*alpha0vpar[7])+0.6123724356957944*alpha0vpar[6]-0.3535533905932737*alpha0vpar[5]-0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])+0.3535533905932737*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.6123724356957944*(alpha0vpar[7]+alpha0vpar[6])+0.3535533905932737*alpha0vpar[5]+0.6123724356957944*(alpha1vpar[4]+alpha0vpar[4])+0.3535533905932737*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[1]+alpha0vpar[1]+alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.6123724356957944*(alpha0x[4]*f1[4]+alpha0x[2]*f1[2]+alpha0x[1]*f1[1]+alpha0x[0]*f1[0]); 
  out[2] += 0.6123724356957944*(alpha0vpar[7]*f1[7]+alpha0vpar[6]*f1[6]+alpha0vpar[5]*f1[5]+alpha0vpar[4]*f1[4]+alpha1vpar[4]*f0[4]+alpha0vpar[3]*f1[3]+alpha0vpar[2]*f1[2]+alpha1vpar[2]*f0[2]+alpha0vpar[1]*f1[1]+alpha1vpar[1]*f0[1]+alpha0vpar[0]*f1[0]+alpha1vpar[0]*f0[0]); 
  out[4] += 0.6123724356957944*(alpha0vpar[6]*f1[7]+f1[6]*alpha0vpar[7]+alpha0vpar[3]*f1[5]+f1[3]*alpha0vpar[5]+(alpha0vpar[2]+alpha0x[1])*f1[4]+alpha1vpar[2]*f0[4]+f0[2]*alpha1vpar[4]+f1[1]*alpha0x[4]+f1[2]*(alpha0vpar[4]+alpha0x[0])+f1[0]*alpha0x[2]+alpha0vpar[0]*f1[1]+alpha1vpar[0]*f0[1]+f0[0]*alpha1vpar[1]+f1[0]*alpha0vpar[1]); 
  out[5] += 0.6123724356957944*(alpha0x[4]*f1[7]+alpha0x[2]*f1[6]+alpha0x[1]*f1[5]+alpha0x[0]*f1[3]); 
  out[6] += 0.6123724356957944*(alpha0vpar[4]*f1[7]+alpha1vpar[4]*f0[7]+f1[4]*alpha0vpar[7]+alpha0vpar[2]*f1[6]+alpha1vpar[2]*f0[6]+f1[2]*alpha0vpar[6]+alpha0vpar[1]*f1[5]+alpha1vpar[1]*f0[5]+f1[1]*alpha0vpar[5]+alpha0vpar[0]*f1[3]+alpha1vpar[0]*f0[3]+f1[0]*alpha0vpar[3]); 
  out[7] += 0.6123724356957944*((alpha0vpar[2]+alpha0x[1])*f1[7]+alpha1vpar[2]*f0[7]+f1[2]*alpha0vpar[7]+(alpha0vpar[4]+alpha0x[0])*f1[6]+alpha1vpar[4]*f0[6]+f1[4]*alpha0vpar[6]+(alpha0x[4]+alpha0vpar[0])*f1[5]+alpha1vpar[0]*f0[5]+f1[0]*alpha0vpar[5]+(alpha0x[2]+alpha0vpar[1])*f1[3]+alpha1vpar[1]*f0[3]+f1[1]*alpha0vpar[3]); 
  out[2] += 0.6123724356957944*(alpha1vpar[4]*f1[4]+alpha1vpar[2]*f1[2]+alpha1vpar[1]*f1[1]+alpha1vpar[0]*f1[0]); 
  out[4] += 0.6123724356957944*(alpha1vpar[2]*f1[4]+f1[2]*alpha1vpar[4]+alpha1vpar[0]*f1[1]+f1[0]*alpha1vpar[1]); 
  out[6] += 0.6123724356957944*(alpha1vpar[4]*f1[7]+alpha1vpar[2]*f1[6]+alpha1vpar[1]*f1[5]+alpha1vpar[0]*f1[3]); 
  out[7] += 0.6123724356957944*(alpha1vpar[2]*f1[7]+alpha1vpar[4]*f1[6]+alpha1vpar[0]*f1[5]+alpha1vpar[1]*f1[3]); 
  return cflFreq; 
} 
