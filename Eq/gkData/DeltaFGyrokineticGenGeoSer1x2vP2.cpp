#include <DeltaFGyrokineticModDecl.h> 
double DeltaFGyrokineticGenGeoVol1x2vSerP2_Bvars(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0, const double *f1, double *out) 
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

  double hamil0[20]; 
  hamil0[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*bmag[0]*wmu)+2.0*m_))/rdvpar2Sq; 
  hamil0[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil0[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil0[8] = (0.421637021355784*m_)/rdvpar2Sq; 

  double hamil1[20]; 
  hamil1[0] = 2.0*phi[0]*q_; 
  hamil1[1] = 2.0*phi[1]*q_; 
  hamil1[7] = 2.0*phi[2]*q_; 

  double BstarZdBmag[20]; 
  BstarZdBmag[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[20]; 
  alpha0x[0] = (0.6123724356957944*BstarZdBmag[0]*hamil0[2]*rdvpar2*rdx2)/m_; 
  alpha0x[2] = (1.369306393762915*BstarZdBmag[0]*hamil0[8]*rdvpar2*rdx2)/m_; 
  double alpha1x[20]; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.3535533905932737*alpha0x[0]-0.4743416490252568*alpha0x[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha0x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alpha0x[2]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alpha0x[0]-0.4743416490252568*alpha0x[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alpha0x[2]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alpha0x[0]-0.4743416490252568*alpha0x[2]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha0x[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alpha0x[2]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.3535533905932737*alpha0x[0]-0.4743416490252568*alpha0x[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha0x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alpha0x[2]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alpha0x[0]-0.4743416490252568*alpha0x[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alpha0x[2]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alpha0x[0]-0.4743416490252568*alpha0x[2]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha0x[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alpha0x[2]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[20]; 
  double alpha1vpar[20]; 
  alpha1vpar[0] = -(0.6123724356957944*BstarZdBmag[0]*hamil1[1]*rdvpar2*rdx2)/m_; 
  alpha1vpar[1] = -(1.369306393762915*BstarZdBmag[0]*hamil1[7]*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.3535533905932737*alpha1vpar[0]-0.4743416490252568*alpha1vpar[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alpha1vpar[1]+0.3535533905932737*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alpha1vpar[0]-0.4743416490252568*alpha1vpar[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alpha1vpar[1]+0.3535533905932737*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.3535533905932737*alpha1vpar[0]-0.4743416490252568*alpha1vpar[1]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.4743416490252568*alpha1vpar[1]+0.3535533905932737*alpha1vpar[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.3535533905932737*alpha1vpar[0]-0.4743416490252568*alpha1vpar[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alpha1vpar[1]+0.3535533905932737*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alpha1vpar[0]-0.4743416490252568*alpha1vpar[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alpha1vpar[1]+0.3535533905932737*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.3535533905932737*alpha1vpar[0]-0.4743416490252568*alpha1vpar[1]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.0441941738241592*alpha1vpar[0]; 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.4743416490252568*alpha1vpar[1]+0.3535533905932737*alpha1vpar[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.6123724356957944*alpha0x[2]*f1[2]+0.6123724356957944*alpha0x[0]*f1[0]; 
  out[2] += 0.6123724356957944*alpha1vpar[1]*f0[1]+0.6123724356957944*alpha1vpar[0]*f0[0]; 
  out[4] += 0.5477225575051661*alpha0x[2]*f1[8]+0.5477225575051661*alpha1vpar[1]*f0[7]+0.6123724356957944*alpha0x[0]*f1[2]+0.6123724356957944*f1[0]*alpha0x[2]+0.6123724356957944*alpha1vpar[0]*f0[1]+0.6123724356957944*f0[0]*alpha1vpar[1]; 
  out[5] += 0.6123724356957944*alpha0x[2]*f1[6]+0.6123724356957944*alpha0x[0]*f1[3]; 
  out[6] += 0.6123724356957944*alpha1vpar[1]*f0[5]+0.6123724356957944*alpha1vpar[0]*f0[3]; 
  out[7] += 1.369306393762915*alpha0x[2]*f1[4]+1.369306393762915*alpha0x[0]*f1[1]; 
  out[8] += 1.369306393762915*alpha1vpar[1]*f0[4]+1.369306393762915*alpha1vpar[0]*f0[2]; 
  out[10] += 0.5477225575051661*alpha0x[2]*f1[14]+0.5477225575051661*alpha1vpar[1]*f0[13]+0.6123724356957944*alpha0x[0]*f1[6]+0.6123724356957944*alpha1vpar[0]*f0[5]+0.6123724356957944*alpha0x[2]*f1[3]+0.6123724356957944*alpha1vpar[1]*f0[3]; 
  out[11] += 1.224744871391589*alpha0x[2]*f1[12]+0.6123724356957944*alpha1vpar[0]*f0[7]+1.369306393762915*alpha0x[0]*f1[4]+1.369306393762915*f1[1]*alpha0x[2]+0.5477225575051661*alpha1vpar[1]*f0[1]; 
  out[12] += 1.224744871391589*alpha1vpar[1]*f0[11]+0.6123724356957944*alpha0x[0]*f1[8]+1.369306393762915*alpha1vpar[0]*f0[4]+0.5477225575051661*alpha0x[2]*f1[2]+1.369306393762915*alpha1vpar[1]*f0[2]; 
  out[13] += 1.369306393762915*alpha0x[2]*f1[10]+1.369306393762915*alpha0x[0]*f1[5]; 
  out[14] += 1.369306393762915*alpha1vpar[1]*f0[10]+1.369306393762915*alpha1vpar[0]*f0[6]; 
  out[15] += 0.6123724356957944*alpha0x[2]*f1[16]+0.6123724356957944*alpha0x[0]*f1[9]; 
  out[16] += 0.6123724356957944*alpha1vpar[1]*f0[15]+0.6123724356957944*alpha1vpar[0]*f0[9]; 
  out[17] += 1.224744871391589*alpha0x[2]*f1[18]+0.6123724356957944*alpha1vpar[0]*f0[13]+1.369306393762915*alpha0x[0]*f1[10]+1.369306393762915*alpha0x[2]*f1[5]+0.5477225575051661*alpha1vpar[1]*f0[5]; 
  out[18] += 1.224744871391589*alpha1vpar[1]*f0[17]+0.6123724356957944*alpha0x[0]*f1[14]+1.369306393762915*alpha1vpar[0]*f0[10]+0.5477225575051661*alpha0x[2]*f1[6]+1.369306393762915*alpha1vpar[1]*f0[6]; 
  out[19] += 0.6123724356957944*alpha0x[0]*f1[16]+0.6123724356957944*alpha1vpar[0]*f0[15]+0.6123724356957944*alpha0x[2]*f1[9]+0.6123724356957944*alpha1vpar[1]*f0[9]; 
  out[2] += 0.6123724356957944*alpha1vpar[1]*f1[1]+0.6123724356957944*alpha1vpar[0]*f1[0]; 
  out[4] += 0.5477225575051661*alpha1vpar[1]*f1[7]+0.6123724356957944*alpha1vpar[0]*f1[1]+0.6123724356957944*f1[0]*alpha1vpar[1]; 
  out[6] += 0.6123724356957944*alpha1vpar[1]*f1[5]+0.6123724356957944*alpha1vpar[0]*f1[3]; 
  out[8] += 1.369306393762915*alpha1vpar[1]*f1[4]+1.369306393762915*alpha1vpar[0]*f1[2]; 
  out[10] += 0.5477225575051661*alpha1vpar[1]*f1[13]+0.6123724356957944*alpha1vpar[0]*f1[5]+0.6123724356957944*alpha1vpar[1]*f1[3]; 
  out[11] += 0.6123724356957944*alpha1vpar[0]*f1[7]+0.5477225575051661*alpha1vpar[1]*f1[1]; 
  out[12] += 1.224744871391589*alpha1vpar[1]*f1[11]+1.369306393762915*alpha1vpar[0]*f1[4]+1.369306393762915*alpha1vpar[1]*f1[2]; 
  out[14] += 1.369306393762915*alpha1vpar[1]*f1[10]+1.369306393762915*alpha1vpar[0]*f1[6]; 
  out[16] += 0.6123724356957944*alpha1vpar[1]*f1[15]+0.6123724356957944*alpha1vpar[0]*f1[9]; 
  out[17] += 0.6123724356957944*alpha1vpar[0]*f1[13]+0.5477225575051661*alpha1vpar[1]*f1[5]; 
  out[18] += 1.224744871391589*alpha1vpar[1]*f1[17]+1.369306393762915*alpha1vpar[0]*f1[10]+1.369306393762915*alpha1vpar[1]*f1[6]; 
  out[19] += 0.6123724356957944*alpha1vpar[0]*f1[15]+0.6123724356957944*alpha1vpar[1]*f1[9]; 
  return cflFreq; 
} 
double DeltaFGyrokineticGenGeoVol1x2vSerP2_Bvarsx(const double q_, const double m_, const double *w, const double *dxv, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0, const double *f1, double *out) 
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

  double hamil0[20]; 
  hamil0[0] = (0.2357022603955158*(3.0*rdvpar2Sq*(2.0*m_*wvparSq+2.828427124746191*bmag[0]*wmu)+2.0*m_))/rdvpar2Sq; 
  hamil0[1] = 2.0*bmag[1]*wmu; 
  hamil0[2] = (1.632993161855453*m_*wvpar)/rdvpar2; 
  hamil0[3] = (1.154700538379252*bmag[0])/rdmu2; 
  hamil0[5] = (1.154700538379252*bmag[1])/rdmu2; 
  hamil0[7] = 2.0*bmag[2]*wmu; 
  hamil0[8] = (0.421637021355784*m_)/rdvpar2Sq; 
  hamil0[13] = (1.154700538379251*bmag[2])/rdmu2; 

  double hamil1[20]; 
  hamil1[0] = 2.0*phi[0]*q_; 
  hamil1[1] = 2.0*phi[1]*q_; 
  hamil1[7] = 2.0*phi[2]*q_; 

  double BstarZdBmag[20]; 
  BstarZdBmag[0] = (1.414213562373095*(1.732050807568877*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2*wvpar+(cmag[2]*jacobTotInv[2]+cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmag[1] = (0.2828427124746191*(1.732050807568877*(b_y[2]*(10.0*jacobTotInv[2]+11.18033988749895*jacobTotInv[0])+5.0*b_y[1]*jacobTotInv[1])*m_*rdx2*wvpar+(4.47213595499958*(cmag[1]*jacobTotInv[2]+jacobTotInv[1]*cmag[2])+5.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))*q_))/q_; 
  BstarZdBmag[2] = (1.414213562373095*(2.23606797749979*jacobTotInv[1]*b_y[2]+jacobTotInv[0]*b_y[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[4] = (1.414213562373095*(b_y[2]*(2.0*jacobTotInv[2]+2.23606797749979*jacobTotInv[0])+b_y[1]*jacobTotInv[1])*m_*rdx2)/(q_*rdvpar2); 
  BstarZdBmag[7] = (0.04040610178208843*(60.6217782649107*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2*wvpar+((22.3606797749979*cmag[2]+35.0*cmag[0])*jacobTotInv[2]+7.0*(5.0*jacobTotInv[0]*cmag[2]+4.47213595499958*cmag[1]*jacobTotInv[1]))*q_))/q_; 
  BstarZdBmag[11] = (1.414213562373095*(b_y[1]*jacobTotInv[2]+2.0*jacobTotInv[1]*b_y[2])*m_*rdx2)/(q_*rdvpar2); 

  double cflFreq = 0.0; 
  double alphaL = 0.0; 
  double alphaR = 0.0; 
  double alpha0x[20]; 
  alpha0x[0] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[2]*hamil0[8]+BstarZdBmag[0]*hamil0[2])*rdvpar2*rdx2)/m_; 
  alpha0x[1] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[4]*hamil0[8]+BstarZdBmag[1]*hamil0[2])*rdvpar2*rdx2)/m_; 
  alpha0x[2] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[0]*hamil0[8]+BstarZdBmag[2]*hamil0[2])*rdvpar2*rdx2)/m_; 
  alpha0x[4] = (0.6123724356957944*(2.23606797749979*BstarZdBmag[1]*hamil0[8]+hamil0[2]*BstarZdBmag[4])*rdvpar2*rdx2)/m_; 
  alpha0x[7] = (0.3535533905932737*(3.872983346207417*hamil0[8]*BstarZdBmag[11]+1.732050807568877*hamil0[2]*BstarZdBmag[7])*rdvpar2*rdx2)/m_; 
  alpha0x[8] = (1.224744871391589*BstarZdBmag[2]*hamil0[8]*rdvpar2*rdx2)/m_; 
  alpha0x[11] = (0.3535533905932737*(1.732050807568877*hamil0[2]*BstarZdBmag[11]+3.872983346207417*BstarZdBmag[7]*hamil0[8])*rdvpar2*rdx2)/m_; 
  alpha0x[12] = (1.224744871391589*BstarZdBmag[4]*hamil0[8]*rdvpar2*rdx2)/m_; 
  double alpha1x[20]; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*((-0.5477225575051661*alpha0x[12])-1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]+0.8215838362577489*alpha0x[4]-0.4743416490252568*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6846531968814574*alpha0x[12]-0.3952847075210473*alpha0x[8]+0.7905694150420947*alpha0x[7]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.5477225575051661*alpha0x[12])+1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]-0.8215838362577489*alpha0x[4]+0.4743416490252568*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.5477225575051661*alpha0x[12])-1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]+0.8215838362577489*alpha0x[4]-0.4743416490252568*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.5477225575051661*alpha0x[12])+1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]-0.8215838362577489*alpha0x[4]+0.4743416490252568*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.5477225575051661*alpha0x[12])-1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]+0.8215838362577489*alpha0x[4]-0.4743416490252568*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.6846531968814574*alpha0x[12]-0.3952847075210473*alpha0x[8]+0.7905694150420947*alpha0x[7]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.5477225575051661*alpha0x[12])+1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]-0.8215838362577489*alpha0x[4]+0.4743416490252568*alpha0x[2]-0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*(0.5477225575051661*alpha0x[12]-1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]-0.8215838362577489*alpha0x[4]-0.4743416490252568*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6846531968814574*alpha0x[12])-0.3952847075210473*alpha0x[8]+0.7905694150420947*alpha0x[7]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.5477225575051661*alpha0x[12]+1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]+0.8215838362577489*alpha0x[4]+0.4743416490252568*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.5477225575051661*alpha0x[12]-1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]-0.8215838362577489*alpha0x[4]-0.4743416490252568*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.5477225575051661*alpha0x[12]+1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]+0.8215838362577489*alpha0x[4]+0.4743416490252568*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.5477225575051661*alpha0x[12]-1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]-0.8215838362577489*alpha0x[4]-0.4743416490252568*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.6846531968814574*alpha0x[12])-0.3952847075210473*alpha0x[8]+0.7905694150420947*alpha0x[7]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.5477225575051661*alpha0x[12]+1.060660171779821*alpha0x[11]+0.3162277660168379*alpha0x[8]+0.7905694150420947*alpha0x[7]+0.8215838362577489*alpha0x[4]+0.4743416490252568*alpha0x[2]+0.6123724356957944*alpha0x[1]+0.3535533905932737*alpha0x[0]); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  double alpha0vpar[20]; 
  alpha0vpar[0] = -(0.6123724356957944*(2.23606797749979*BstarZdBmag[1]*hamil0[7]+BstarZdBmag[0]*hamil0[1])*rdvpar2*rdx2)/m_; 
  alpha0vpar[1] = -(0.6123724356957944*((2.0*BstarZdBmag[7]+2.23606797749979*BstarZdBmag[0])*hamil0[7]+BstarZdBmag[1]*hamil0[1])*rdvpar2*rdx2)/m_; 
  alpha0vpar[2] = -(0.6123724356957944*(2.23606797749979*BstarZdBmag[4]*hamil0[7]+hamil0[1]*BstarZdBmag[2])*rdvpar2*rdx2)/m_; 
  alpha0vpar[3] = -(0.3535533905932737*(3.872983346207417*BstarZdBmag[1]*hamil0[13]+1.732050807568877*BstarZdBmag[0]*hamil0[5])*rdvpar2*rdx2)/m_; 
  alpha0vpar[4] = -(0.3535533905932737*(3.464101615137755*hamil0[7]*BstarZdBmag[11]+1.732050807568877*(2.23606797749979*BstarZdBmag[2]*hamil0[7]+hamil0[1]*BstarZdBmag[4]))*rdvpar2*rdx2)/m_; 
  alpha0vpar[5] = -(0.3535533905932737*(3.872983346207417*(0.8944271909999159*BstarZdBmag[7]+BstarZdBmag[0])*hamil0[13]+1.732050807568877*BstarZdBmag[1]*hamil0[5])*rdvpar2*rdx2)/m_; 
  alpha0vpar[6] = -(0.3535533905932737*(3.872983346207417*BstarZdBmag[4]*hamil0[13]+1.732050807568877*BstarZdBmag[2]*hamil0[5])*rdvpar2*rdx2)/m_; 
  alpha0vpar[7] = -(0.6123724356957944*(2.0*BstarZdBmag[1]*hamil0[7]+hamil0[1]*BstarZdBmag[7])*rdvpar2*rdx2)/m_; 
  alpha0vpar[10] = -(0.3535533905932737*((3.464101615137754*BstarZdBmag[11]+3.872983346207417*BstarZdBmag[2])*hamil0[13]+1.732050807568877*BstarZdBmag[4]*hamil0[5])*rdvpar2*rdx2)/m_; 
  alpha0vpar[11] = -(0.3535533905932737*(1.732050807568877*hamil0[1]*BstarZdBmag[11]+3.464101615137755*BstarZdBmag[4]*hamil0[7])*rdvpar2*rdx2)/m_; 
  alpha0vpar[13] = -(0.3535533905932737*(3.464101615137754*BstarZdBmag[1]*hamil0[13]+1.732050807568877*hamil0[5]*BstarZdBmag[7])*rdvpar2*rdx2)/m_; 
  alpha0vpar[17] = -(0.6123724356957944*(2.0*BstarZdBmag[4]*hamil0[13]+hamil0[5]*BstarZdBmag[11])*rdvpar2*rdx2)/m_; 
  double alpha1vpar[20]; 
  alpha1vpar[0] = -(0.6123724356957944*(2.23606797749979*BstarZdBmag[1]*hamil1[7]+BstarZdBmag[0]*hamil1[1])*rdvpar2*rdx2)/m_; 
  alpha1vpar[1] = -(0.6123724356957944*((2.0*BstarZdBmag[7]+2.23606797749979*BstarZdBmag[0])*hamil1[7]+BstarZdBmag[1]*hamil1[1])*rdvpar2*rdx2)/m_; 
  alpha1vpar[2] = -(0.6123724356957944*(2.23606797749979*BstarZdBmag[4]*hamil1[7]+hamil1[1]*BstarZdBmag[2])*rdvpar2*rdx2)/m_; 
  alpha1vpar[4] = -(0.3535533905932737*(3.464101615137755*hamil1[7]*BstarZdBmag[11]+1.732050807568877*(2.23606797749979*BstarZdBmag[2]*hamil1[7]+hamil1[1]*BstarZdBmag[4]))*rdvpar2*rdx2)/m_; 
  alpha1vpar[7] = -(0.6123724356957944*(2.0*BstarZdBmag[1]*hamil1[7]+hamil1[1]*BstarZdBmag[7])*rdvpar2*rdx2)/m_; 
  alpha1vpar[11] = -(0.3535533905932737*(1.732050807568877*hamil1[1]*BstarZdBmag[11]+3.464101615137755*BstarZdBmag[4]*hamil1[7])*rdvpar2*rdx2)/m_; 
  // Evaluate alpha at left surface quadrature points.
  alphaL = 0.125*(0.7348469228349533*alpha0vpar[17]-0.4242640687119285*alpha0vpar[13]-0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11])-1.10227038425243*alpha0vpar[10]+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])+0.8215838362577489*alpha0vpar[6]+0.6363961030678926*alpha0vpar[5]+0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])-0.4743416490252568*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.9185586535436913*alpha0vpar[17])+0.5303300858899104*alpha0vpar[13]+0.6846531968814574*(alpha1vpar[11]+alpha0vpar[11])-0.3952847075210473*(alpha1vpar[7]+alpha0vpar[7])+0.8215838362577489*alpha0vpar[6]-0.4743416490252568*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.7348469228349533*alpha0vpar[17]-0.4242640687119285*alpha0vpar[13]-0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11])+1.10227038425243*alpha0vpar[10]+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])+0.8215838362577489*alpha0vpar[6]-0.6363961030678926*alpha0vpar[5]-0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])-0.4743416490252568*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11]))+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])+0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11]))+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])-0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.7348469228349533*alpha0vpar[17])+0.4242640687119285*alpha0vpar[13]-0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11])+1.10227038425243*alpha0vpar[10]+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])-0.8215838362577489*alpha0vpar[6]-0.6363961030678926*alpha0vpar[5]+0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])+0.4743416490252568*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*(0.9185586535436913*alpha0vpar[17]-0.5303300858899104*alpha0vpar[13]+0.6846531968814574*(alpha1vpar[11]+alpha0vpar[11])-0.3952847075210473*(alpha1vpar[7]+alpha0vpar[7])-0.8215838362577489*alpha0vpar[6]+0.4743416490252568*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  alphaL = 0.125*((-0.7348469228349533*alpha0vpar[17])+0.4242640687119285*alpha0vpar[13]-0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11])-1.10227038425243*alpha0vpar[10]+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])-0.8215838362577489*alpha0vpar[6]+0.6363961030678926*alpha0vpar[5]-0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])+0.4743416490252568*alpha0vpar[3]-0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += -0.5*(alphaL-std::abs(alphaL)); 
  // Evaluate alpha at right surface quadrature points.
  alphaR = 0.125*((-0.7348469228349533*alpha0vpar[17])-0.4242640687119285*alpha0vpar[13]+0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11])+1.10227038425243*alpha0vpar[10]+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])-0.8215838362577489*alpha0vpar[6]+0.6363961030678926*alpha0vpar[5]-0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])-0.4743416490252568*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.9185586535436913*alpha0vpar[17]+0.5303300858899104*alpha0vpar[13]-0.6846531968814574*(alpha1vpar[11]+alpha0vpar[11])-0.3952847075210473*(alpha1vpar[7]+alpha0vpar[7])-0.8215838362577489*alpha0vpar[6]-0.4743416490252568*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.7348469228349533*alpha0vpar[17])-0.4242640687119285*alpha0vpar[13]+0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11])-1.10227038425243*alpha0vpar[10]+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])-0.8215838362577489*alpha0vpar[6]-0.6363961030678926*alpha0vpar[5]+0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])-0.4743416490252568*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11])+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])-0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11])+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])+0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.7348469228349533*alpha0vpar[17]+0.4242640687119285*alpha0vpar[13]+0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11])-1.10227038425243*alpha0vpar[10]+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])+0.8215838362577489*alpha0vpar[6]-0.6363961030678926*alpha0vpar[5]-0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])+0.4743416490252568*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])-0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*((-0.9185586535436913*alpha0vpar[17])-0.5303300858899104*alpha0vpar[13]-0.6846531968814574*(alpha1vpar[11]+alpha0vpar[11])-0.3952847075210473*(alpha1vpar[7]+alpha0vpar[7])+0.8215838362577489*alpha0vpar[6]+0.4743416490252568*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 
  alphaR = 0.125*(0.7348469228349533*alpha0vpar[17]+0.4242640687119285*alpha0vpar[13]+0.5477225575051661*(alpha1vpar[11]+alpha0vpar[11])+1.10227038425243*alpha0vpar[10]+0.3162277660168379*(alpha1vpar[7]+alpha0vpar[7])+0.8215838362577489*alpha0vpar[6]+0.6363961030678926*alpha0vpar[5]+0.8215838362577489*(alpha1vpar[4]+alpha0vpar[4])+0.4743416490252568*alpha0vpar[3]+0.6123724356957944*(alpha1vpar[2]+alpha0vpar[2])+0.4743416490252568*(alpha1vpar[1]+alpha0vpar[1])+0.3535533905932737*(alpha1vpar[0]+alpha0vpar[0])); 
  cflFreq += 0.5*(alphaR+std::abs(alphaR)); 

  out[1] += 0.6123724356957944*alpha0x[12]*f1[12]+0.6123724356957944*alpha0x[11]*f1[11]+0.6123724356957944*alpha0x[8]*f1[8]+0.6123724356957944*alpha0x[7]*f1[7]+0.6123724356957944*alpha0x[4]*f1[4]+0.6123724356957944*alpha0x[2]*f1[2]+0.6123724356957944*alpha0x[1]*f1[1]+0.6123724356957944*alpha0x[0]*f1[0]; 
  out[2] += 0.6123724356957944*alpha0vpar[17]*f1[17]+0.6123724356957944*alpha0vpar[13]*f1[13]+0.6123724356957944*alpha0vpar[11]*f1[11]+0.6123724356957944*alpha1vpar[11]*f0[11]+0.6123724356957944*alpha0vpar[10]*f1[10]+0.6123724356957944*alpha0vpar[7]*f1[7]+0.6123724356957944*alpha1vpar[7]*f0[7]+0.6123724356957944*alpha0vpar[6]*f1[6]+0.6123724356957944*alpha0vpar[5]*f1[5]+0.6123724356957944*alpha0vpar[4]*f1[4]+0.6123724356957944*alpha1vpar[4]*f0[4]+0.6123724356957944*alpha0vpar[3]*f1[3]+0.6123724356957944*alpha0vpar[2]*f1[2]+0.6123724356957944*alpha1vpar[2]*f0[2]+0.6123724356957944*alpha0vpar[1]*f1[1]+0.6123724356957944*alpha1vpar[1]*f0[1]+0.6123724356957944*alpha0vpar[0]*f1[0]+0.6123724356957944*alpha1vpar[0]*f0[0]; 
  out[4] += 0.5477225575051661*alpha0vpar[10]*f1[17]+0.5477225575051661*f1[10]*alpha0vpar[17]+0.5477225575051661*alpha0vpar[5]*f1[13]+0.5477225575051661*f1[5]*alpha0vpar[13]+0.5477225575051661*alpha0x[4]*f1[12]+0.5477225575051661*f1[4]*alpha0x[12]+0.6123724356957944*alpha0x[7]*f1[11]+0.5477225575051661*alpha0vpar[4]*f1[11]+0.5477225575051661*alpha1vpar[4]*f0[11]+0.5477225575051661*f0[4]*alpha1vpar[11]+0.6123724356957944*f1[7]*alpha0x[11]+0.5477225575051661*f1[4]*alpha0vpar[11]+0.6123724356957944*alpha0vpar[6]*f1[10]+0.6123724356957944*f1[6]*alpha0vpar[10]+0.5477225575051661*alpha0x[2]*f1[8]+0.5477225575051661*f1[2]*alpha0x[8]+0.5477225575051661*alpha0vpar[1]*f1[7]+0.5477225575051661*alpha1vpar[1]*f0[7]+0.5477225575051661*f0[1]*alpha1vpar[7]+0.5477225575051661*f1[1]*alpha0vpar[7]+0.6123724356957944*alpha0vpar[3]*f1[5]+0.6123724356957944*f1[3]*alpha0vpar[5]+0.6123724356957944*alpha0vpar[2]*f1[4]+0.6123724356957944*alpha0x[1]*f1[4]+0.6123724356957944*alpha1vpar[2]*f0[4]+0.6123724356957944*f0[2]*alpha1vpar[4]+0.6123724356957944*f1[1]*alpha0x[4]+0.6123724356957944*f1[2]*alpha0vpar[4]+0.6123724356957944*alpha0x[0]*f1[2]+0.6123724356957944*f1[0]*alpha0x[2]+0.6123724356957944*alpha0vpar[0]*f1[1]+0.6123724356957944*alpha1vpar[0]*f0[1]+0.6123724356957944*f0[0]*alpha1vpar[1]+0.6123724356957944*f1[0]*alpha0vpar[1]; 
  out[5] += 0.6123724356957944*alpha0x[12]*f1[18]+0.6123724356957944*alpha0x[11]*f1[17]+0.6123724356957944*alpha0x[8]*f1[14]+0.6123724356957944*alpha0x[7]*f1[13]+0.6123724356957944*alpha0x[4]*f1[10]+0.6123724356957944*alpha0x[2]*f1[6]+0.6123724356957944*alpha0x[1]*f1[5]+0.6123724356957944*alpha0x[0]*f1[3]; 
  out[6] += 0.5477225575051661*alpha0vpar[10]*f1[19]+0.6123724356957944*alpha0vpar[11]*f1[17]+0.6123724356957944*alpha1vpar[11]*f0[17]+0.6123724356957944*f1[11]*alpha0vpar[17]+0.5477225575051661*alpha0vpar[6]*f1[16]+0.5477225575051661*alpha0vpar[5]*f1[15]+0.6123724356957944*alpha0vpar[7]*f1[13]+0.6123724356957944*alpha1vpar[7]*f0[13]+0.6123724356957944*f1[7]*alpha0vpar[13]+0.6123724356957944*alpha0vpar[4]*f1[10]+0.6123724356957944*alpha1vpar[4]*f0[10]+0.6123724356957944*f1[4]*alpha0vpar[10]+0.5477225575051661*alpha0vpar[3]*f1[9]+0.6123724356957944*alpha0vpar[2]*f1[6]+0.6123724356957944*alpha1vpar[2]*f0[6]+0.6123724356957944*f1[2]*alpha0vpar[6]+0.6123724356957944*alpha0vpar[1]*f1[5]+0.6123724356957944*alpha1vpar[1]*f0[5]+0.6123724356957944*f1[1]*alpha0vpar[5]+0.6123724356957944*alpha0vpar[0]*f1[3]+0.6123724356957944*alpha1vpar[0]*f0[3]+0.6123724356957944*f1[0]*alpha0vpar[3]; 
  out[7] += 1.369306393762915*alpha0x[8]*f1[12]+1.369306393762915*f1[8]*alpha0x[12]+1.224744871391589*alpha0x[4]*f1[11]+1.224744871391589*f1[4]*alpha0x[11]+1.224744871391589*alpha0x[1]*f1[7]+1.224744871391589*f1[1]*alpha0x[7]+1.369306393762915*alpha0x[2]*f1[4]+1.369306393762915*f1[2]*alpha0x[4]+1.369306393762915*alpha0x[0]*f1[1]+1.369306393762915*f1[0]*alpha0x[1]; 
  out[8] += 1.224744871391589*alpha0vpar[10]*f1[18]+1.369306393762915*alpha0vpar[13]*f1[17]+1.369306393762915*f1[13]*alpha0vpar[17]+1.224744871391589*alpha0vpar[6]*f1[14]+1.224744871391589*alpha0vpar[4]*f1[12]+1.224744871391589*alpha1vpar[4]*f0[12]+1.369306393762915*alpha0vpar[7]*f1[11]+1.369306393762915*alpha1vpar[7]*f0[11]+1.369306393762915*f0[7]*alpha1vpar[11]+1.369306393762915*f1[7]*alpha0vpar[11]+1.369306393762915*alpha0vpar[5]*f1[10]+1.369306393762915*f1[5]*alpha0vpar[10]+1.224744871391589*alpha0vpar[2]*f1[8]+1.224744871391589*alpha1vpar[2]*f0[8]+1.369306393762915*alpha0vpar[3]*f1[6]+1.369306393762915*f1[3]*alpha0vpar[6]+1.369306393762915*alpha0vpar[1]*f1[4]+1.369306393762915*alpha1vpar[1]*f0[4]+1.369306393762915*f0[1]*alpha1vpar[4]+1.369306393762915*f1[1]*alpha0vpar[4]+1.369306393762915*alpha0vpar[0]*f1[2]+1.369306393762915*alpha1vpar[0]*f0[2]+1.369306393762915*f0[0]*alpha1vpar[2]+1.369306393762915*f1[0]*alpha0vpar[2]; 
  out[10] += 0.4898979485566357*alpha0vpar[17]*f1[19]+0.5477225575051661*alpha0vpar[6]*f1[19]+0.5477225575051661*alpha0x[4]*f1[18]+0.6123724356957944*alpha0x[7]*f1[17]+0.5477225575051661*alpha0vpar[4]*f1[17]+0.5477225575051661*alpha1vpar[4]*f0[17]+0.5477225575051661*f1[4]*alpha0vpar[17]+0.5477225575051661*alpha0vpar[10]*f1[16]+0.4898979485566357*alpha0vpar[13]*f1[15]+0.5477225575051661*alpha0vpar[3]*f1[15]+0.5477225575051661*alpha0x[2]*f1[14]+0.6123724356957944*alpha0x[11]*f1[13]+0.5477225575051661*alpha0vpar[1]*f1[13]+0.5477225575051661*alpha1vpar[1]*f0[13]+0.5477225575051661*f1[1]*alpha0vpar[13]+0.5477225575051661*f1[10]*alpha0x[12]+0.5477225575051661*alpha0vpar[10]*f1[11]+0.5477225575051661*f0[10]*alpha1vpar[11]+0.5477225575051661*f1[10]*alpha0vpar[11]+0.6123724356957944*alpha0vpar[2]*f1[10]+0.6123724356957944*alpha0x[1]*f1[10]+0.6123724356957944*alpha1vpar[2]*f0[10]+0.6123724356957944*f1[2]*alpha0vpar[10]+0.5477225575051661*alpha0vpar[5]*f1[9]+0.5477225575051661*f1[6]*alpha0x[8]+0.5477225575051661*alpha0vpar[5]*f1[7]+0.5477225575051661*f0[5]*alpha1vpar[7]+0.5477225575051661*f1[5]*alpha0vpar[7]+0.6123724356957944*alpha0vpar[4]*f1[6]+0.6123724356957944*alpha0x[0]*f1[6]+0.6123724356957944*alpha1vpar[4]*f0[6]+0.6123724356957944*f1[4]*alpha0vpar[6]+0.6123724356957944*alpha0x[4]*f1[5]+0.6123724356957944*alpha0vpar[0]*f1[5]+0.6123724356957944*alpha1vpar[0]*f0[5]+0.6123724356957944*f1[0]*alpha0vpar[5]+0.6123724356957944*alpha0x[2]*f1[3]+0.6123724356957944*alpha0vpar[1]*f1[3]+0.6123724356957944*alpha1vpar[1]*f0[3]+0.6123724356957944*f1[1]*alpha0vpar[3]; 
  out[11] += 0.3912303982179757*alpha0vpar[17]*f1[17]+0.6123724356957944*alpha0vpar[6]*f1[17]+0.6123724356957944*f1[6]*alpha0vpar[17]+0.3912303982179757*alpha0vpar[13]*f1[13]+0.6123724356957944*alpha0vpar[3]*f1[13]+0.6123724356957944*f1[3]*alpha0vpar[13]+1.095445115010332*alpha0x[11]*f1[12]+1.224744871391589*alpha0x[2]*f1[12]+1.095445115010332*f1[11]*alpha0x[12]+1.224744871391589*f1[2]*alpha0x[12]+0.3912303982179757*alpha0vpar[11]*f1[11]+0.6123724356957944*alpha0vpar[2]*f1[11]+1.224744871391589*alpha0x[1]*f1[11]+0.3912303982179757*alpha1vpar[11]*f0[11]+0.6123724356957944*alpha1vpar[2]*f0[11]+0.6123724356957944*f0[2]*alpha1vpar[11]+1.224744871391589*f1[1]*alpha0x[11]+0.6123724356957944*f1[2]*alpha0vpar[11]+0.5477225575051661*alpha0vpar[10]*f1[10]+1.224744871391589*alpha0x[4]*f1[8]+1.224744871391589*f1[4]*alpha0x[8]+0.3912303982179757*alpha0vpar[7]*f1[7]+1.224744871391589*alpha0x[4]*f1[7]+0.6123724356957944*alpha0vpar[0]*f1[7]+0.3912303982179757*alpha1vpar[7]*f0[7]+0.6123724356957944*alpha1vpar[0]*f0[7]+0.6123724356957944*f0[0]*alpha1vpar[7]+1.224744871391589*f1[4]*alpha0x[7]+0.6123724356957944*f1[0]*alpha0vpar[7]+0.5477225575051661*alpha0vpar[5]*f1[5]+0.5477225575051661*alpha0vpar[4]*f1[4]+1.369306393762915*alpha0x[0]*f1[4]+0.5477225575051661*alpha1vpar[4]*f0[4]+1.369306393762915*f1[0]*alpha0x[4]+1.369306393762915*alpha0x[1]*f1[2]+1.369306393762915*f1[1]*alpha0x[2]+0.5477225575051661*alpha0vpar[1]*f1[1]+0.5477225575051661*alpha1vpar[1]*f0[1]; 
  out[12] += 1.095445115010332*alpha0vpar[17]*f1[18]+1.224744871391589*alpha0vpar[6]*f1[18]+1.224744871391589*alpha0vpar[5]*f1[17]+1.224744871391589*f1[5]*alpha0vpar[17]+1.224744871391589*alpha0vpar[10]*f1[14]+1.224744871391589*alpha0vpar[10]*f1[13]+1.224744871391589*f1[10]*alpha0vpar[13]+0.3912303982179757*alpha0x[12]*f1[12]+1.095445115010332*alpha0vpar[11]*f1[12]+1.224744871391589*alpha0vpar[2]*f1[12]+0.6123724356957944*alpha0x[1]*f1[12]+1.095445115010332*alpha1vpar[11]*f0[12]+1.224744871391589*alpha1vpar[2]*f0[12]+0.6123724356957944*f1[1]*alpha0x[12]+0.5477225575051661*alpha0x[11]*f1[11]+1.224744871391589*alpha0vpar[1]*f1[11]+1.224744871391589*alpha1vpar[1]*f0[11]+1.224744871391589*f0[1]*alpha1vpar[11]+1.224744871391589*f1[1]*alpha0vpar[11]+1.369306393762915*alpha0vpar[3]*f1[10]+1.369306393762915*f1[3]*alpha0vpar[10]+0.3912303982179757*alpha0x[8]*f1[8]+1.224744871391589*alpha0vpar[4]*f1[8]+0.6123724356957944*alpha0x[0]*f1[8]+1.224744871391589*alpha1vpar[4]*f0[8]+0.6123724356957944*f1[0]*alpha0x[8]+1.224744871391589*alpha0vpar[4]*f1[7]+1.224744871391589*alpha1vpar[4]*f0[7]+1.224744871391589*f0[4]*alpha1vpar[7]+1.224744871391589*f1[4]*alpha0vpar[7]+1.369306393762915*alpha0vpar[5]*f1[6]+1.369306393762915*f1[5]*alpha0vpar[6]+0.5477225575051661*alpha0x[4]*f1[4]+1.369306393762915*alpha0vpar[0]*f1[4]+1.369306393762915*alpha1vpar[0]*f0[4]+1.369306393762915*f0[0]*alpha1vpar[4]+1.369306393762915*f1[0]*alpha0vpar[4]+0.5477225575051661*alpha0x[2]*f1[2]+1.369306393762915*alpha0vpar[1]*f1[2]+1.369306393762915*alpha1vpar[1]*f0[2]+1.369306393762915*f0[1]*alpha1vpar[2]+1.369306393762915*f1[1]*alpha0vpar[2]; 
  out[13] += 1.369306393762915*alpha0x[8]*f1[18]+1.224744871391589*alpha0x[4]*f1[17]+1.369306393762915*alpha0x[12]*f1[14]+1.224744871391589*alpha0x[1]*f1[13]+1.224744871391589*f1[10]*alpha0x[11]+1.369306393762915*alpha0x[2]*f1[10]+1.224744871391589*f1[5]*alpha0x[7]+1.369306393762915*alpha0x[4]*f1[6]+1.369306393762915*alpha0x[0]*f1[5]+1.369306393762915*alpha0x[1]*f1[3]; 
  out[14] += 1.224744871391589*alpha0vpar[5]*f1[19]+1.224744871391589*alpha0vpar[4]*f1[18]+1.224744871391589*alpha1vpar[4]*f0[18]+1.369306393762915*alpha0vpar[7]*f1[17]+1.369306393762915*alpha1vpar[7]*f0[17]+1.369306393762915*f1[7]*alpha0vpar[17]+1.224744871391589*alpha0vpar[3]*f1[16]+1.224744871391589*alpha0vpar[10]*f1[15]+1.224744871391589*alpha0vpar[2]*f1[14]+1.224744871391589*alpha1vpar[2]*f0[14]+1.369306393762915*alpha0vpar[11]*f1[13]+1.369306393762915*alpha1vpar[11]*f0[13]+1.369306393762915*f1[11]*alpha0vpar[13]+1.224744871391589*alpha0vpar[10]*f1[12]+1.369306393762915*alpha0vpar[1]*f1[10]+1.369306393762915*alpha1vpar[1]*f0[10]+1.369306393762915*f1[1]*alpha0vpar[10]+1.224744871391589*alpha0vpar[6]*f1[9]+1.224744871391589*alpha0vpar[6]*f1[8]+1.369306393762915*alpha0vpar[0]*f1[6]+1.369306393762915*alpha1vpar[0]*f0[6]+1.369306393762915*f1[0]*alpha0vpar[6]+1.369306393762915*alpha0vpar[4]*f1[5]+1.369306393762915*alpha1vpar[4]*f0[5]+1.369306393762915*f1[4]*alpha0vpar[5]+1.369306393762915*alpha0vpar[2]*f1[3]+1.369306393762915*alpha1vpar[2]*f0[3]+1.369306393762915*f1[2]*alpha0vpar[3]; 
  out[15] += 0.6123724356957944*alpha0x[4]*f1[19]+0.6123724356957944*alpha0x[2]*f1[16]+0.6123724356957944*alpha0x[1]*f1[15]+0.6123724356957944*alpha0x[0]*f1[9]; 
  out[16] += 0.6123724356957944*alpha0vpar[4]*f1[19]+0.6123724356957944*alpha1vpar[4]*f0[19]+0.5477225575051661*alpha0vpar[17]*f1[17]+0.6123724356957944*alpha0vpar[2]*f1[16]+0.6123724356957944*alpha1vpar[2]*f0[16]+0.6123724356957944*alpha0vpar[1]*f1[15]+0.6123724356957944*alpha1vpar[1]*f0[15]+0.5477225575051661*alpha0vpar[13]*f1[13]+0.5477225575051661*alpha0vpar[10]*f1[10]+0.6123724356957944*alpha0vpar[0]*f1[9]+0.6123724356957944*alpha1vpar[0]*f0[9]+0.5477225575051661*alpha0vpar[6]*f1[6]+0.5477225575051661*alpha0vpar[5]*f1[5]+0.5477225575051661*alpha0vpar[3]*f1[3]; 
  out[17] += 0.4898979485566357*alpha0vpar[10]*f1[19]+1.095445115010332*alpha0x[11]*f1[18]+1.224744871391589*alpha0x[2]*f1[18]+1.095445115010332*alpha0x[12]*f1[17]+0.3912303982179757*alpha0vpar[11]*f1[17]+0.6123724356957944*alpha0vpar[2]*f1[17]+1.224744871391589*alpha0x[1]*f1[17]+0.3912303982179757*alpha1vpar[11]*f0[17]+0.6123724356957944*alpha1vpar[2]*f0[17]+0.5477225575051661*f1[16]*alpha0vpar[17]+0.3912303982179757*f1[11]*alpha0vpar[17]+0.6123724356957944*f1[2]*alpha0vpar[17]+0.4898979485566356*alpha0vpar[5]*f1[15]+1.224744871391589*alpha0x[4]*f1[14]+0.3912303982179757*alpha0vpar[7]*f1[13]+1.224744871391589*alpha0x[4]*f1[13]+0.6123724356957944*alpha0vpar[0]*f1[13]+0.3912303982179757*alpha1vpar[7]*f0[13]+0.6123724356957944*alpha1vpar[0]*f0[13]+0.5477225575051661*f1[9]*alpha0vpar[13]+0.3912303982179757*f1[7]*alpha0vpar[13]+0.6123724356957944*f1[0]*alpha0vpar[13]+1.224744871391589*f1[6]*alpha0x[12]+0.6123724356957944*alpha0vpar[6]*f1[11]+0.6123724356957944*f0[6]*alpha1vpar[11]+1.224744871391589*f1[5]*alpha0x[11]+0.6123724356957944*f1[6]*alpha0vpar[11]+1.224744871391589*alpha0x[8]*f1[10]+1.224744871391589*alpha0x[7]*f1[10]+0.5477225575051661*alpha0vpar[4]*f1[10]+1.369306393762915*alpha0x[0]*f1[10]+0.5477225575051661*alpha1vpar[4]*f0[10]+0.5477225575051661*f1[4]*alpha0vpar[10]+0.6123724356957944*alpha0vpar[3]*f1[7]+0.6123724356957944*f0[3]*alpha1vpar[7]+0.6123724356957944*f1[3]*alpha0vpar[7]+1.369306393762915*alpha0x[1]*f1[6]+1.369306393762915*alpha0x[2]*f1[5]+0.5477225575051661*alpha0vpar[1]*f1[5]+0.5477225575051661*alpha1vpar[1]*f0[5]+0.5477225575051661*f1[1]*alpha0vpar[5]+1.369306393762915*f1[3]*alpha0x[4]; 
  out[18] += 1.095445115010332*alpha0vpar[13]*f1[19]+1.224744871391589*alpha0vpar[3]*f1[19]+0.3912303982179757*alpha0x[12]*f1[18]+1.095445115010332*alpha0vpar[11]*f1[18]+1.224744871391589*alpha0vpar[2]*f1[18]+0.6123724356957944*alpha0x[1]*f1[18]+1.095445115010332*alpha1vpar[11]*f0[18]+1.224744871391589*alpha1vpar[2]*f0[18]+0.5477225575051661*alpha0x[11]*f1[17]+1.224744871391589*alpha0vpar[1]*f1[17]+1.224744871391589*alpha1vpar[1]*f0[17]+1.095445115010332*f1[15]*alpha0vpar[17]+1.095445115010332*f1[12]*alpha0vpar[17]+1.224744871391589*f1[1]*alpha0vpar[17]+1.224744871391589*alpha0vpar[5]*f1[16]+1.224744871391589*alpha0vpar[6]*f1[15]+0.3912303982179757*alpha0x[8]*f1[14]+1.224744871391589*alpha0vpar[4]*f1[14]+0.6123724356957944*alpha0x[0]*f1[14]+1.224744871391589*alpha1vpar[4]*f0[14]+1.224744871391589*alpha0vpar[4]*f1[13]+1.224744871391589*alpha1vpar[4]*f0[13]+1.224744871391589*f1[4]*alpha0vpar[13]+1.224744871391589*alpha0vpar[6]*f1[12]+0.6123724356957944*f1[5]*alpha0x[12]+1.224744871391589*alpha0vpar[5]*f1[11]+1.224744871391589*f0[5]*alpha1vpar[11]+1.224744871391589*f1[5]*alpha0vpar[11]+1.224744871391589*alpha0vpar[7]*f1[10]+0.5477225575051661*alpha0x[4]*f1[10]+1.369306393762915*alpha0vpar[0]*f1[10]+1.224744871391589*alpha1vpar[7]*f0[10]+1.369306393762915*alpha1vpar[0]*f0[10]+1.224744871391589*f1[9]*alpha0vpar[10]+1.224744871391589*f1[8]*alpha0vpar[10]+1.224744871391589*f1[7]*alpha0vpar[10]+1.369306393762915*f1[0]*alpha0vpar[10]+0.6123724356957944*f1[3]*alpha0x[8]+0.5477225575051661*alpha0x[2]*f1[6]+1.369306393762915*alpha0vpar[1]*f1[6]+1.369306393762915*alpha1vpar[1]*f0[6]+1.369306393762915*f1[1]*alpha0vpar[6]+1.369306393762915*alpha0vpar[2]*f1[5]+1.369306393762915*alpha1vpar[2]*f0[5]+1.369306393762915*f1[2]*alpha0vpar[5]+1.369306393762915*alpha0vpar[3]*f1[4]+1.369306393762915*f0[3]*alpha1vpar[4]+1.369306393762915*f1[3]*alpha0vpar[4]; 
  out[19] += 0.5477225575051661*alpha0x[12]*f1[19]+0.5477225575051661*alpha0vpar[11]*f1[19]+0.6123724356957944*alpha0vpar[2]*f1[19]+0.6123724356957944*alpha0x[1]*f1[19]+0.5477225575051661*alpha1vpar[11]*f0[19]+0.6123724356957944*alpha1vpar[2]*f0[19]+0.4898979485566357*alpha0vpar[10]*f1[17]+0.4898979485566357*f1[10]*alpha0vpar[17]+0.5477225575051661*alpha0x[8]*f1[16]+0.6123724356957944*alpha0vpar[4]*f1[16]+0.6123724356957944*alpha0x[0]*f1[16]+0.6123724356957944*alpha1vpar[4]*f0[16]+0.5477225575051661*alpha0vpar[7]*f1[15]+0.6123724356957944*alpha0x[4]*f1[15]+0.6123724356957944*alpha0vpar[0]*f1[15]+0.5477225575051661*alpha1vpar[7]*f0[15]+0.6123724356957944*alpha1vpar[0]*f0[15]+0.4898979485566356*alpha0vpar[5]*f1[13]+0.4898979485566356*f1[5]*alpha0vpar[13]+0.5477225575051661*alpha0vpar[6]*f1[10]+0.5477225575051661*f1[6]*alpha0vpar[10]+0.6123724356957944*alpha0x[2]*f1[9]+0.6123724356957944*alpha0vpar[1]*f1[9]+0.6123724356957944*alpha1vpar[1]*f0[9]+0.5477225575051661*alpha0vpar[3]*f1[5]+0.5477225575051661*f1[3]*alpha0vpar[5]; 
  out[2] += 0.6123724356957944*alpha1vpar[11]*f1[11]+0.6123724356957944*alpha1vpar[7]*f1[7]+0.6123724356957944*alpha1vpar[4]*f1[4]+0.6123724356957944*alpha1vpar[2]*f1[2]+0.6123724356957944*alpha1vpar[1]*f1[1]+0.6123724356957944*alpha1vpar[0]*f1[0]; 
  out[4] += 0.5477225575051661*alpha1vpar[4]*f1[11]+0.5477225575051661*f1[4]*alpha1vpar[11]+0.5477225575051661*alpha1vpar[1]*f1[7]+0.5477225575051661*f1[1]*alpha1vpar[7]+0.6123724356957944*alpha1vpar[2]*f1[4]+0.6123724356957944*f1[2]*alpha1vpar[4]+0.6123724356957944*alpha1vpar[0]*f1[1]+0.6123724356957944*f1[0]*alpha1vpar[1]; 
  out[6] += 0.6123724356957944*alpha1vpar[11]*f1[17]+0.6123724356957944*alpha1vpar[7]*f1[13]+0.6123724356957944*alpha1vpar[4]*f1[10]+0.6123724356957944*alpha1vpar[2]*f1[6]+0.6123724356957944*alpha1vpar[1]*f1[5]+0.6123724356957944*alpha1vpar[0]*f1[3]; 
  out[8] += 1.224744871391589*alpha1vpar[4]*f1[12]+1.369306393762915*alpha1vpar[7]*f1[11]+1.369306393762915*f1[7]*alpha1vpar[11]+1.224744871391589*alpha1vpar[2]*f1[8]+1.369306393762915*alpha1vpar[1]*f1[4]+1.369306393762915*f1[1]*alpha1vpar[4]+1.369306393762915*alpha1vpar[0]*f1[2]+1.369306393762915*f1[0]*alpha1vpar[2]; 
  out[10] += 0.5477225575051661*alpha1vpar[4]*f1[17]+0.5477225575051661*alpha1vpar[1]*f1[13]+0.5477225575051661*f1[10]*alpha1vpar[11]+0.6123724356957944*alpha1vpar[2]*f1[10]+0.5477225575051661*f1[5]*alpha1vpar[7]+0.6123724356957944*alpha1vpar[4]*f1[6]+0.6123724356957944*alpha1vpar[0]*f1[5]+0.6123724356957944*alpha1vpar[1]*f1[3]; 
  out[11] += 0.3912303982179757*alpha1vpar[11]*f1[11]+0.6123724356957944*alpha1vpar[2]*f1[11]+0.6123724356957944*f1[2]*alpha1vpar[11]+0.3912303982179757*alpha1vpar[7]*f1[7]+0.6123724356957944*alpha1vpar[0]*f1[7]+0.6123724356957944*f1[0]*alpha1vpar[7]+0.5477225575051661*alpha1vpar[4]*f1[4]+0.5477225575051661*alpha1vpar[1]*f1[1]; 
  out[12] += 1.095445115010332*alpha1vpar[11]*f1[12]+1.224744871391589*alpha1vpar[2]*f1[12]+1.224744871391589*alpha1vpar[1]*f1[11]+1.224744871391589*f1[1]*alpha1vpar[11]+1.224744871391589*alpha1vpar[4]*f1[8]+1.224744871391589*alpha1vpar[4]*f1[7]+1.224744871391589*f1[4]*alpha1vpar[7]+1.369306393762915*alpha1vpar[0]*f1[4]+1.369306393762915*f1[0]*alpha1vpar[4]+1.369306393762915*alpha1vpar[1]*f1[2]+1.369306393762915*f1[1]*alpha1vpar[2]; 
  out[14] += 1.224744871391589*alpha1vpar[4]*f1[18]+1.369306393762915*alpha1vpar[7]*f1[17]+1.224744871391589*alpha1vpar[2]*f1[14]+1.369306393762915*alpha1vpar[11]*f1[13]+1.369306393762915*alpha1vpar[1]*f1[10]+1.369306393762915*alpha1vpar[0]*f1[6]+1.369306393762915*alpha1vpar[4]*f1[5]+1.369306393762915*alpha1vpar[2]*f1[3]; 
  out[16] += 0.6123724356957944*alpha1vpar[4]*f1[19]+0.6123724356957944*alpha1vpar[2]*f1[16]+0.6123724356957944*alpha1vpar[1]*f1[15]+0.6123724356957944*alpha1vpar[0]*f1[9]; 
  out[17] += 0.3912303982179757*alpha1vpar[11]*f1[17]+0.6123724356957944*alpha1vpar[2]*f1[17]+0.3912303982179757*alpha1vpar[7]*f1[13]+0.6123724356957944*alpha1vpar[0]*f1[13]+0.6123724356957944*f1[6]*alpha1vpar[11]+0.5477225575051661*alpha1vpar[4]*f1[10]+0.6123724356957944*f1[3]*alpha1vpar[7]+0.5477225575051661*alpha1vpar[1]*f1[5]; 
  out[18] += 1.095445115010332*alpha1vpar[11]*f1[18]+1.224744871391589*alpha1vpar[2]*f1[18]+1.224744871391589*alpha1vpar[1]*f1[17]+1.224744871391589*alpha1vpar[4]*f1[14]+1.224744871391589*alpha1vpar[4]*f1[13]+1.224744871391589*f1[5]*alpha1vpar[11]+1.224744871391589*alpha1vpar[7]*f1[10]+1.369306393762915*alpha1vpar[0]*f1[10]+1.369306393762915*alpha1vpar[1]*f1[6]+1.369306393762915*alpha1vpar[2]*f1[5]+1.369306393762915*f1[3]*alpha1vpar[4]; 
  out[19] += 0.5477225575051661*alpha1vpar[11]*f1[19]+0.6123724356957944*alpha1vpar[2]*f1[19]+0.6123724356957944*alpha1vpar[4]*f1[16]+0.5477225575051661*alpha1vpar[7]*f1[15]+0.6123724356957944*alpha1vpar[0]*f1[15]+0.6123724356957944*alpha1vpar[1]*f1[9]; 
  return cflFreq; 
} 
