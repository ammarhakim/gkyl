#include <DeltaFGyrokineticModDecl.h>
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_x_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagR[16]; 

  double alpha0R[8]; 

  double alpha1R[8]; 
  alpha1R[0] = (0.1767766952966368*b_z[0]*jacobTotInv[0]*(3.0*hamil1R[5]-1.732050807568877*hamil1R[2])*rdy2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*(3.0*b_z[0]*jacobTotInv[0]*hamil1R[5]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamil1R[2])*rdy2R)/q_; 

  double incr[16]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])+0.3535533905932737*(f0R[3]+f0R[2])-0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*(f0R[5]+f0R[4])-0.6123724356957944*(f0L[5]+f0L[4])-0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[1] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])+0.3535533905932737*(f0R[3]+f0R[2])-0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*(f0R[5]+f0R[4])-0.6123724356957944*(f0L[5]+f0L[4])-0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5])+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*(f0L[3]+f0R[2])+0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*f0R[5]-0.6123724356957944*(f0L[5]+f0R[4])+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])+0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[3] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5])+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*(f0L[3]+f0R[2])+0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*f0R[5]-0.6123724356957944*(f0L[5]+f0R[4])+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])+0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[4] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5])-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*(f0L[3]+f0R[2])-0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*f0R[5]+0.6123724356957944*(f0L[5]+f0R[4])-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[5] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5])-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*(f0L[3]+f0R[2])-0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*f0R[5]+0.6123724356957944*(f0L[5]+f0R[4])-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[6] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.3535533905932737*(f0R[3]+f0R[2])+0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*(f0R[5]+f0R[4])+0.6123724356957944*(f0L[5]+f0L[4])+0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[7] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.3535533905932737*(f0R[3]+f0R[2])+0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*(f0R[5]+f0R[4])+0.6123724356957944*(f0L[5]+f0L[4])+0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*alpha1R[0]*fUp[0]; 
  incr[1] = -0.4330127018922193*alpha1R[0]*fUp[0]; 
  incr[2] = 0.25*alpha1R[0]*fUp[1]; 
  incr[3] = 0.25*alpha1R[0]*fUp[2]; 
  incr[4] = 0.25*alpha1R[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alpha1R[0]*fUp[1]; 
  incr[6] = -0.4330127018922193*alpha1R[0]*fUp[2]; 
  incr[7] = 0.25*alpha1R[0]*fUp[4]; 
  incr[8] = -0.4330127018922193*alpha1R[0]*fUp[3]; 
  incr[9] = 0.25*alpha1R[0]*fUp[5]; 
  incr[10] = 0.25*alpha1R[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alpha1R[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alpha1R[0]*fUp[5]; 
  incr[13] = -0.4330127018922193*alpha1R[0]*fUp[6]; 
  incr[14] = 0.25*alpha1R[0]*fUp[7]; 
  incr[15] = -0.4330127018922193*alpha1R[0]*fUp[7]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 
  outR[8] += incr[8]*rdx2R; 
  outR[9] += incr[9]*rdx2R; 
  outR[10] += incr[10]*rdx2R; 
  outR[11] += incr[11]*rdx2R; 
  outR[12] += incr[12]*rdx2R; 
  outR[13] += incr[13]*rdx2R; 
  outR[14] += incr[14]*rdx2R; 
  outR[15] += incr[15]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += -1.0*incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += incr[6]*rdx2L; 
  outL[7] += -1.0*incr[7]*rdx2L; 
  outL[8] += incr[8]*rdx2L; 
  outL[9] += -1.0*incr[9]*rdx2L; 
  outL[10] += -1.0*incr[10]*rdx2L; 
  outL[11] += incr[11]*rdx2L; 
  outL[12] += incr[12]*rdx2L; 
  outL[13] += incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_y_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarYdBmagR[16]; 

  double alpha0R[8]; 

  double alpha1R[8]; 
  alpha1R[0] = -(0.1767766952966368*b_z[0]*jacobTotInv[0]*(3.0*hamil1R[5]-1.732050807568877*hamil1R[1])*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(3.0*b_z[0]*jacobTotInv[0]*hamil1R[5]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamil1R[1])*rdx2R)/q_; 

  double incr[16]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[0] = 0.5*((0.3535533905932737*f0R[7]-0.3535533905932737*(f0L[7]+f0R[6]+f0R[5]+f0R[4])+0.3535533905932737*(f0L[6]+f0L[5]+f0L[4]+f0R[3]+f0R[2]+f0R[1])-0.3535533905932737*(f0L[3]+f0L[2]+f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.3535533905932737*(f0R[7]+f0L[7])+0.3535533905932737*(f0R[6]+f0L[6]+f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2]+f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[1] = 0.5*(((-0.3535533905932737*(f0R[7]+f0R[6]))+0.3535533905932737*(f0L[7]+f0L[6]+f0R[5]+f0R[4]+f0R[3]+f0R[2])-0.3535533905932737*(f0L[5]+f0L[4]+f0L[3]+f0L[2]+f0R[1]+f0R[0])+0.3535533905932737*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.3535533905932737*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.3535533905932737*(f0R[5]+f0L[5]+f0R[4]+f0L[4]+f0R[3]+f0L[3]+f0R[2]+f0L[2])+0.3535533905932737*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[2] = 0.5*(((-0.3535533905932737*f0R[7])+0.3535533905932737*(f0L[7]+f0R[6])-0.3535533905932737*(f0L[6]+f0R[5])+0.3535533905932737*(f0L[5]+f0R[4]+f0R[3])-0.3535533905932737*(f0L[4]+f0L[3]+f0R[2])+0.3535533905932737*(f0L[2]+f0R[1])-0.3535533905932737*(f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.3535533905932737*(f0R[7]+f0L[7])-0.3535533905932737*(f0R[6]+f0L[6])+0.3535533905932737*(f0R[5]+f0L[5])-0.3535533905932737*(f0R[4]+f0L[4]+f0R[3]+f0L[3])+0.3535533905932737*(f0R[2]+f0L[2])-0.3535533905932737*(f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[3] = 0.5*((0.3535533905932737*(f0R[7]+f0R[6]+f0R[5])-0.3535533905932737*(f0L[7]+f0L[6]+f0L[5]+f0R[4])+0.3535533905932737*(f0L[4]+f0R[3])-0.3535533905932737*(f0L[3]+f0R[2]+f0R[1]+f0R[0])+0.3535533905932737*(f0L[2]+f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.3535533905932737*(f0R[7]+f0L[7]+f0R[6]+f0L[6]+f0R[5]+f0L[5])+0.3535533905932737*(f0R[4]+f0L[4])-0.3535533905932737*(f0R[3]+f0L[3])+0.3535533905932737*(f0R[2]+f0L[2]+f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[4] = 0.5*(((-0.3535533905932737*f0R[7])+0.3535533905932737*(f0L[7]+f0R[6]+f0R[5])-0.3535533905932737*(f0L[6]+f0L[5]+f0R[4]+f0R[3])+0.3535533905932737*(f0L[4]+f0L[3]+f0R[2]+f0R[1])-0.3535533905932737*(f0L[2]+f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.3535533905932737*(f0R[7]+f0L[7])-0.3535533905932737*(f0R[6]+f0L[6]+f0R[5]+f0L[5])+0.3535533905932737*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.3535533905932737*(f0R[2]+f0L[2]+f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[5] = 0.5*((0.3535533905932737*(f0R[7]+f0R[6])-0.3535533905932737*(f0L[7]+f0L[6]+f0R[5])+0.3535533905932737*(f0L[5]+f0R[4])-0.3535533905932737*(f0L[4]+f0R[3])+0.3535533905932737*(f0L[3]+f0R[2])-0.3535533905932737*(f0L[2]+f0R[1]+f0R[0])+0.3535533905932737*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.3535533905932737*(f0R[7]+f0L[7]+f0R[6]+f0L[6])+0.3535533905932737*(f0R[5]+f0L[5])-0.3535533905932737*(f0R[4]+f0L[4])+0.3535533905932737*(f0R[3]+f0L[3])-0.3535533905932737*(f0R[2]+f0L[2])+0.3535533905932737*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[6] = 0.5*((0.3535533905932737*f0R[7]-0.3535533905932737*(f0L[7]+f0R[6])+0.3535533905932737*(f0L[6]+f0R[5]+f0R[4])-0.3535533905932737*(f0L[5]+f0L[4]+f0R[3]+f0R[2])+0.3535533905932737*(f0L[3]+f0L[2]+f0R[1])-0.3535533905932737*(f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.3535533905932737*(f0R[7]+f0L[7])+0.3535533905932737*(f0R[6]+f0L[6])-0.3535533905932737*(f0R[5]+f0L[5]+f0R[4]+f0L[4])+0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.3535533905932737*(f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[7] = 0.5*(0.3535533905932737*(f0L[7]+f0L[6]+f0L[5]+f0L[4]+f0L[3]+f0L[2]+f0L[1]+f0L[0])-0.3535533905932737*(f0R[7]+f0R[6]+f0R[5]+f0R[4]+f0R[3]+f0R[2]+f0R[1]+f0R[0]))*sgn(alphaOrdR)+0.1767766952966368*(f0R[7]+f0L[7]+f0R[6]+f0L[6]+f0R[5]+f0L[5]+f0R[4]+f0L[4]+f0R[3]+f0L[3]+f0R[2]+f0L[2]+f0R[1]+f0L[1]+f0R[0]+f0L[0]); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*alpha1R[0]*fUp[0]; 
  incr[1] = 0.25*alpha1R[0]*fUp[1]; 
  incr[2] = -0.4330127018922193*alpha1R[0]*fUp[0]; 
  incr[3] = 0.25*alpha1R[0]*fUp[2]; 
  incr[4] = 0.25*alpha1R[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alpha1R[0]*fUp[1]; 
  incr[6] = 0.25*alpha1R[0]*fUp[4]; 
  incr[7] = -0.4330127018922193*alpha1R[0]*fUp[2]; 
  incr[8] = 0.25*alpha1R[0]*fUp[5]; 
  incr[9] = -0.4330127018922193*alpha1R[0]*fUp[3]; 
  incr[10] = 0.25*alpha1R[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alpha1R[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alpha1R[0]*fUp[5]; 
  incr[13] = 0.25*alpha1R[0]*fUp[7]; 
  incr[14] = -0.4330127018922193*alpha1R[0]*fUp[6]; 
  incr[15] = -0.4330127018922193*alpha1R[0]*fUp[7]; 

  outR[0] += incr[0]*rdy2R; 
  outR[1] += incr[1]*rdy2R; 
  outR[2] += incr[2]*rdy2R; 
  outR[3] += incr[3]*rdy2R; 
  outR[4] += incr[4]*rdy2R; 
  outR[5] += incr[5]*rdy2R; 
  outR[6] += incr[6]*rdy2R; 
  outR[7] += incr[7]*rdy2R; 
  outR[8] += incr[8]*rdy2R; 
  outR[9] += incr[9]*rdy2R; 
  outR[10] += incr[10]*rdy2R; 
  outR[11] += incr[11]*rdy2R; 
  outR[12] += incr[12]*rdy2R; 
  outR[13] += incr[13]*rdy2R; 
  outR[14] += incr[14]*rdy2R; 
  outR[15] += incr[15]*rdy2R; 

  outL[0] += -1.0*incr[0]*rdy2L; 
  outL[1] += -1.0*incr[1]*rdy2L; 
  outL[2] += incr[2]*rdy2L; 
  outL[3] += -1.0*incr[3]*rdy2L; 
  outL[4] += -1.0*incr[4]*rdy2L; 
  outL[5] += incr[5]*rdy2L; 
  outL[6] += -1.0*incr[6]*rdy2L; 
  outL[7] += incr[7]*rdy2L; 
  outL[8] += -1.0*incr[8]*rdy2L; 
  outL[9] += incr[9]*rdy2L; 
  outL[10] += -1.0*incr[10]*rdy2L; 
  outL[11] += incr[11]*rdy2L; 
  outL[12] += incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 

  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagR[16]; 

  double BstarYdBmagR[16]; 

  double alpha0R[8]; 

  double alpha1R[8]; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = 0.0; 

  double incr[16]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // alpha1 == 0, so nothing to do for alpha1*f0 and alpha1*f1 terms.
  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_x_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[8] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagR[16]; 

  double alpha0R[8]; 

  double alpha1R[8]; 
  alpha1R[0] = (0.125*(((12.72792206135786*b_z[1]-7.348469228349534*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(4.242640687119286*b_z[0]-7.348469228349534*b_z[1]))*hamil1R[5]+((4.242640687119286*b_z[0]-7.348469228349534*b_z[1])*jacobTotInv[1]+jacobTotInv[0]*(4.242640687119286*b_z[1]-2.449489742783178*b_z[0]))*hamil1R[2])*rdy2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*(((9.0*b_z[1]-5.196152422706631*b_z[0])*jacobTotInv[1]-5.196152422706631*jacobTotInv[0]*b_z[1]+3.0*b_z[0]*jacobTotInv[0])*hamil1R[5]+((3.0*b_z[0]-5.196152422706631*b_z[1])*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1]-1.732050807568877*b_z[0]*jacobTotInv[0])*hamil1R[2])*rdy2R)/q_; 

  double incr[16]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])+0.3535533905932737*(f0R[3]+f0R[2])-0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*(f0R[5]+f0R[4])-0.6123724356957944*(f0L[5]+f0L[4])-0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[1] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])+0.3535533905932737*(f0R[3]+f0R[2])-0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*(f0R[5]+f0R[4])-0.6123724356957944*(f0L[5]+f0L[4])-0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5])+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*(f0L[3]+f0R[2])+0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*f0R[5]-0.6123724356957944*(f0L[5]+f0R[4])+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])+0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[3] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5])+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*(f0L[3]+f0R[2])+0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*f0R[5]-0.6123724356957944*(f0L[5]+f0R[4])+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])+0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[4] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5])-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*(f0L[3]+f0R[2])-0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*f0R[5]+0.6123724356957944*(f0L[5]+f0R[4])-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[5] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5])-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*(f0L[3]+f0R[2])-0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*f0R[5]+0.6123724356957944*(f0L[5]+f0R[4])-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[6] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.3535533905932737*(f0R[3]+f0R[2])+0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*(f0R[5]+f0R[4])+0.6123724356957944*(f0L[5]+f0L[4])+0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[7] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.3535533905932737*(f0R[3]+f0R[2])+0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*(f0R[5]+f0R[4])+0.6123724356957944*(f0L[5]+f0L[4])+0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*alpha1R[0]*fUp[0]; 
  incr[1] = -0.4330127018922193*alpha1R[0]*fUp[0]; 
  incr[2] = 0.25*alpha1R[0]*fUp[1]; 
  incr[3] = 0.25*alpha1R[0]*fUp[2]; 
  incr[4] = 0.25*alpha1R[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alpha1R[0]*fUp[1]; 
  incr[6] = -0.4330127018922193*alpha1R[0]*fUp[2]; 
  incr[7] = 0.25*alpha1R[0]*fUp[4]; 
  incr[8] = -0.4330127018922193*alpha1R[0]*fUp[3]; 
  incr[9] = 0.25*alpha1R[0]*fUp[5]; 
  incr[10] = 0.25*alpha1R[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alpha1R[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alpha1R[0]*fUp[5]; 
  incr[13] = -0.4330127018922193*alpha1R[0]*fUp[6]; 
  incr[14] = 0.25*alpha1R[0]*fUp[7]; 
  incr[15] = -0.4330127018922193*alpha1R[0]*fUp[7]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 
  outR[8] += incr[8]*rdx2R; 
  outR[9] += incr[9]*rdx2R; 
  outR[10] += incr[10]*rdx2R; 
  outR[11] += incr[11]*rdx2R; 
  outR[12] += incr[12]*rdx2R; 
  outR[13] += incr[13]*rdx2R; 
  outR[14] += incr[14]*rdx2R; 
  outR[15] += incr[15]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += -1.0*incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += incr[6]*rdx2L; 
  outL[7] += -1.0*incr[7]*rdx2L; 
  outL[8] += incr[8]*rdx2L; 
  outL[9] += -1.0*incr[9]*rdx2L; 
  outL[10] += -1.0*incr[10]*rdx2L; 
  outL[11] += incr[11]*rdx2L; 
  outL[12] += incr[12]*rdx2L; 
  outL[13] += incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_y_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[8] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = -(1.732050807568877*jacobTotInv[0]*b_z[1]*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[1] = -(1.732050807568877*b_z[1]*jacobTotInv[1]*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alpha0R[8]; 
  alpha0R[0] = (0.3061862178478971*(hamil0R[1]*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*m_*rdx2R+BstarYdBmagR[0]*hamil0R[3]*q_*rdvpar2R))/(m_*q_); 
  alpha0R[1] = (0.3061862178478971*(hamil0R[1]*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*m_*rdx2R+BstarYdBmagR[1]*hamil0R[3]*q_*rdvpar2R))/(m_*q_); 
  alpha0R[2] = (0.3061862178478971*BstarYdBmagR[3]*hamil0R[3]*rdvpar2R)/m_; 
  alpha0R[3] = (0.3061862178478971*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil0R[8]*rdx2R)/q_; 
  alpha0R[4] = (0.3061862178478971*hamil0R[3]*BstarYdBmagR[6]*rdvpar2R)/m_; 
  alpha0R[5] = (0.3061862178478971*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil0R[8]*rdx2R)/q_; 

  double alpha1R[8]; 
  alpha1R[0] = -(0.125*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*(4.242640687119286*hamil1R[5]-2.449489742783178*hamil1R[1])*rdx2R)/q_; 
  alpha1R[1] = -(0.125*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*(4.242640687119286*hamil1R[5]-2.449489742783178*hamil1R[1])*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(((3.0*b_z[1]*jacobTotInv[1]+3.0*b_z[0]*jacobTotInv[0])*hamil1R[5]+((-1.732050807568877*b_z[1]*hamil1R[1])-1.732050807568877*b_z[1]*hamil0R[1])*jacobTotInv[1]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamil1R[1]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamil0R[1])*m_*rdx2R-1.732050807568877*BstarYdBmagR[0]*hamil0R[3]*q_*rdvpar2R))/(m_*q_); 

  double incr[16]; 
  // linear term alpha0*f1 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*(alpha0R[5]+alpha0R[4])-0.3535533905932737*(alpha0R[3]+alpha0R[2]+alpha0R[1])+0.3535533905932737*alpha0R[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(f1R[15]+f1L[15]))+0.4330127018922193*(f1R[14]+f1L[14])+0.25*f1R[13]-0.25*f1L[13]+0.4330127018922193*(f1R[12]+f1L[12]+f1R[11]+f1L[11])-0.25*f1R[10]+0.25*f1L[10]-0.4330127018922193*(f1R[9]+f1L[9])-0.25*f1R[8]+0.25*f1L[8]-0.4330127018922193*(f1R[7]+f1L[7])-0.25*f1R[6]+0.25*f1L[6]-0.4330127018922193*(f1R[5]+f1L[5])+0.25*(f1R[4]+f1R[3])-0.25*(f1L[4]+f1L[3])+0.4330127018922193*(f1R[2]+f1L[2])+0.25*f1R[1]-0.25*(f1L[1]+f1R[0])+0.25*f1L[0])*sgn(alphaOrdR)+0.4330127018922193*f1R[15]-0.4330127018922193*(f1L[15]+f1R[14])+0.4330127018922193*f1L[14]-0.25*(f1R[13]+f1L[13])-0.4330127018922193*(f1R[12]+f1R[11])+0.4330127018922193*(f1L[12]+f1L[11])+0.25*(f1R[10]+f1L[10])+0.4330127018922193*f1R[9]-0.4330127018922193*f1L[9]+0.25*(f1R[8]+f1L[8])+0.4330127018922193*f1R[7]-0.4330127018922193*f1L[7]+0.25*(f1R[6]+f1L[6])+0.4330127018922193*f1R[5]-0.4330127018922193*f1L[5]-0.25*(f1R[4]+f1L[4]+f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]-0.25*(f1R[1]+f1L[1])+0.25*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha0R[1]+alpha0R[0])-0.3535533905932737*(alpha0R[5]+alpha0R[4]+alpha0R[3]+alpha0R[2]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(f1R[15]+f1L[15]+f1R[14]+f1L[14])-0.25*f1R[13]+0.25*f1L[13]-0.4330127018922193*(f1R[12]+f1L[12]+f1R[11]+f1L[11])-0.25*f1R[10]+0.25*f1L[10]-0.4330127018922193*(f1R[9]+f1L[9])+0.25*f1R[8]-0.25*f1L[8]-0.4330127018922193*(f1R[7]+f1L[7])+0.25*f1R[6]-0.25*f1L[6]+0.4330127018922193*(f1R[5]+f1L[5])+0.25*(f1R[4]+f1R[3])-0.25*(f1L[4]+f1L[3])+0.4330127018922193*(f1R[2]+f1L[2])-0.25*(f1R[1]+f1R[0])+0.25*(f1L[1]+f1L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f1R[15]+f1R[14])+0.4330127018922193*(f1L[15]+f1L[14])+0.25*(f1R[13]+f1L[13])+0.4330127018922193*(f1R[12]+f1R[11])-0.4330127018922193*(f1L[12]+f1L[11])+0.25*(f1R[10]+f1L[10])+0.4330127018922193*f1R[9]-0.4330127018922193*f1L[9]-0.25*(f1R[8]+f1L[8])+0.4330127018922193*f1R[7]-0.4330127018922193*f1L[7]-0.25*(f1R[6]+f1L[6])-0.4330127018922193*f1R[5]+0.4330127018922193*f1L[5]-0.25*(f1R[4]+f1L[4]+f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]+0.25*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 
  alphaOrdR = 0.3535533905932737*alpha0R[5]-0.3535533905932737*(alpha0R[4]+alpha0R[3])+0.3535533905932737*alpha0R[2]-0.3535533905932737*alpha0R[1]+0.3535533905932737*alpha0R[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(f1R[15]+f1L[15])-0.4330127018922193*(f1R[14]+f1L[14])-0.25*f1R[13]+0.25*f1L[13]+0.4330127018922193*(f1R[12]+f1L[12])-0.4330127018922193*(f1R[11]+f1L[11])+0.25*f1R[10]-0.25*f1L[10]-0.4330127018922193*(f1R[9]+f1L[9])-0.25*f1R[8]+0.25*f1L[8]+0.4330127018922193*(f1R[7]+f1L[7])+0.25*f1R[6]-0.25*f1L[6]-0.4330127018922193*(f1R[5]+f1L[5])+0.25*f1R[4]-0.25*(f1L[4]+f1R[3])+0.25*f1L[3]+0.4330127018922193*(f1R[2]+f1L[2])+0.25*f1R[1]-0.25*(f1L[1]+f1R[0])+0.25*f1L[0])*sgn(alphaOrdR)-0.4330127018922193*f1R[15]+0.4330127018922193*(f1L[15]+f1R[14])-0.4330127018922193*f1L[14]+0.25*(f1R[13]+f1L[13])-0.4330127018922193*f1R[12]+0.4330127018922193*(f1L[12]+f1R[11])-0.4330127018922193*f1L[11]-0.25*(f1R[10]+f1L[10])+0.4330127018922193*f1R[9]-0.4330127018922193*f1L[9]+0.25*(f1R[8]+f1L[8])-0.4330127018922193*f1R[7]+0.4330127018922193*f1L[7]-0.25*(f1R[6]+f1L[6])+0.4330127018922193*f1R[5]-0.4330127018922193*f1L[5]-0.25*(f1R[4]+f1L[4])+0.25*(f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]-0.25*(f1R[1]+f1L[1])+0.25*(f1R[0]+f1L[0])); 
  alphaOrdR = (-0.3535533905932737*alpha0R[5])+0.3535533905932737*alpha0R[4]-0.3535533905932737*alpha0R[3]+0.3535533905932737*(alpha0R[2]+alpha0R[1]+alpha0R[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(f1R[15]+f1L[15]+f1R[14]+f1L[14]))+0.25*f1R[13]-0.25*f1L[13]-0.4330127018922193*(f1R[12]+f1L[12])+0.4330127018922193*(f1R[11]+f1L[11])+0.25*f1R[10]-0.25*f1L[10]-0.4330127018922193*(f1R[9]+f1L[9])+0.25*f1R[8]-0.25*f1L[8]+0.4330127018922193*(f1R[7]+f1L[7])-0.25*f1R[6]+0.25*f1L[6]+0.4330127018922193*(f1R[5]+f1L[5])+0.25*f1R[4]-0.25*(f1L[4]+f1R[3])+0.25*f1L[3]+0.4330127018922193*(f1R[2]+f1L[2])-0.25*(f1R[1]+f1R[0])+0.25*(f1L[1]+f1L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f1R[15]+f1R[14])-0.4330127018922193*(f1L[15]+f1L[14])-0.25*(f1R[13]+f1L[13])+0.4330127018922193*f1R[12]-0.4330127018922193*(f1L[12]+f1R[11])+0.4330127018922193*f1L[11]-0.25*(f1R[10]+f1L[10])+0.4330127018922193*f1R[9]-0.4330127018922193*f1L[9]-0.25*(f1R[8]+f1L[8])-0.4330127018922193*f1R[7]+0.4330127018922193*f1L[7]+0.25*(f1R[6]+f1L[6])-0.4330127018922193*f1R[5]+0.4330127018922193*f1L[5]-0.25*(f1R[4]+f1L[4])+0.25*(f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]+0.25*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 
  alphaOrdR = (-0.3535533905932737*alpha0R[5])+0.3535533905932737*(alpha0R[4]+alpha0R[3])-0.3535533905932737*(alpha0R[2]+alpha0R[1])+0.3535533905932737*alpha0R[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(f1R[15]+f1L[15])-0.4330127018922193*(f1R[14]+f1L[14])-0.25*f1R[13]+0.25*f1L[13]-0.4330127018922193*(f1R[12]+f1L[12])+0.4330127018922193*(f1R[11]+f1L[11])+0.25*f1R[10]-0.25*f1L[10]+0.4330127018922193*(f1R[9]+f1L[9])+0.25*f1R[8]-0.25*f1L[8]-0.4330127018922193*(f1R[7]+f1L[7])-0.25*f1R[6]+0.25*f1L[6]-0.4330127018922193*(f1R[5]+f1L[5])-0.25*f1R[4]+0.25*(f1L[4]+f1R[3])-0.25*f1L[3]+0.4330127018922193*(f1R[2]+f1L[2])+0.25*f1R[1]-0.25*(f1L[1]+f1R[0])+0.25*f1L[0])*sgn(alphaOrdR)-0.4330127018922193*f1R[15]+0.4330127018922193*(f1L[15]+f1R[14])-0.4330127018922193*f1L[14]+0.25*(f1R[13]+f1L[13])+0.4330127018922193*f1R[12]-0.4330127018922193*(f1L[12]+f1R[11])+0.4330127018922193*f1L[11]-0.25*(f1R[10]+f1L[10])-0.4330127018922193*f1R[9]+0.4330127018922193*f1L[9]-0.25*(f1R[8]+f1L[8])+0.4330127018922193*f1R[7]-0.4330127018922193*f1L[7]+0.25*(f1R[6]+f1L[6])+0.4330127018922193*f1R[5]-0.4330127018922193*f1L[5]+0.25*(f1R[4]+f1L[4])-0.25*(f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]-0.25*(f1R[1]+f1L[1])+0.25*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.3535533905932737*alpha0R[5]-0.3535533905932737*alpha0R[4]+0.3535533905932737*alpha0R[3]-0.3535533905932737*alpha0R[2]+0.3535533905932737*(alpha0R[1]+alpha0R[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(f1R[15]+f1L[15]+f1R[14]+f1L[14]))+0.25*f1R[13]-0.25*f1L[13]+0.4330127018922193*(f1R[12]+f1L[12])-0.4330127018922193*(f1R[11]+f1L[11])+0.25*f1R[10]-0.25*f1L[10]+0.4330127018922193*(f1R[9]+f1L[9])-0.25*f1R[8]+0.25*f1L[8]-0.4330127018922193*(f1R[7]+f1L[7])+0.25*f1R[6]-0.25*f1L[6]+0.4330127018922193*(f1R[5]+f1L[5])-0.25*f1R[4]+0.25*(f1L[4]+f1R[3])-0.25*f1L[3]+0.4330127018922193*(f1R[2]+f1L[2])-0.25*(f1R[1]+f1R[0])+0.25*(f1L[1]+f1L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f1R[15]+f1R[14])-0.4330127018922193*(f1L[15]+f1L[14])-0.25*(f1R[13]+f1L[13])-0.4330127018922193*f1R[12]+0.4330127018922193*(f1L[12]+f1R[11])-0.4330127018922193*f1L[11]-0.25*(f1R[10]+f1L[10])-0.4330127018922193*f1R[9]+0.4330127018922193*f1L[9]+0.25*(f1R[8]+f1L[8])+0.4330127018922193*f1R[7]-0.4330127018922193*f1L[7]-0.25*(f1R[6]+f1L[6])-0.4330127018922193*f1R[5]+0.4330127018922193*f1L[5]+0.25*(f1R[4]+f1L[4])-0.25*(f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]+0.25*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 
  alphaOrdR = (-0.3535533905932737*(alpha0R[5]+alpha0R[4]))+0.3535533905932737*(alpha0R[3]+alpha0R[2])-0.3535533905932737*alpha0R[1]+0.3535533905932737*alpha0R[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(f1R[15]+f1L[15]))+0.4330127018922193*(f1R[14]+f1L[14])+0.25*f1R[13]-0.25*f1L[13]-0.4330127018922193*(f1R[12]+f1L[12]+f1R[11]+f1L[11])-0.25*f1R[10]+0.25*f1L[10]+0.4330127018922193*(f1R[9]+f1L[9])+0.25*f1R[8]-0.25*f1L[8]+0.4330127018922193*(f1R[7]+f1L[7])+0.25*f1R[6]-0.25*f1L[6]-0.4330127018922193*(f1R[5]+f1L[5])-0.25*(f1R[4]+f1R[3])+0.25*(f1L[4]+f1L[3])+0.4330127018922193*(f1R[2]+f1L[2])+0.25*f1R[1]-0.25*(f1L[1]+f1R[0])+0.25*f1L[0])*sgn(alphaOrdR)+0.4330127018922193*f1R[15]-0.4330127018922193*(f1L[15]+f1R[14])+0.4330127018922193*f1L[14]-0.25*(f1R[13]+f1L[13])+0.4330127018922193*(f1R[12]+f1R[11])-0.4330127018922193*(f1L[12]+f1L[11])+0.25*(f1R[10]+f1L[10])-0.4330127018922193*f1R[9]+0.4330127018922193*f1L[9]-0.25*(f1R[8]+f1L[8])-0.4330127018922193*f1R[7]+0.4330127018922193*f1L[7]-0.25*(f1R[6]+f1L[6])+0.4330127018922193*f1R[5]-0.4330127018922193*f1L[5]+0.25*(f1R[4]+f1L[4]+f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]-0.25*(f1R[1]+f1L[1])+0.25*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha0R[5]+alpha0R[4]+alpha0R[3]+alpha0R[2]+alpha0R[1]+alpha0R[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(f1R[15]+f1L[15]+f1R[14]+f1L[14])-0.25*f1R[13]+0.25*f1L[13]+0.4330127018922193*(f1R[12]+f1L[12]+f1R[11]+f1L[11])-0.25*f1R[10]+0.25*f1L[10]+0.4330127018922193*(f1R[9]+f1L[9])-0.25*f1R[8]+0.25*f1L[8]+0.4330127018922193*(f1R[7]+f1L[7])-0.25*f1R[6]+0.25*f1L[6]+0.4330127018922193*(f1R[5]+f1L[5])-0.25*(f1R[4]+f1R[3])+0.25*(f1L[4]+f1L[3])+0.4330127018922193*(f1R[2]+f1L[2])-0.25*(f1R[1]+f1R[0])+0.25*(f1L[1]+f1L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f1R[15]+f1R[14])+0.4330127018922193*(f1L[15]+f1L[14])+0.25*(f1R[13]+f1L[13])-0.4330127018922193*(f1R[12]+f1R[11])+0.4330127018922193*(f1L[12]+f1L[11])+0.25*(f1R[10]+f1L[10])-0.4330127018922193*f1R[9]+0.4330127018922193*f1L[9]+0.25*(f1R[8]+f1L[8])-0.4330127018922193*f1R[7]+0.4330127018922193*f1L[7]+0.25*(f1R[6]+f1L[6])-0.4330127018922193*f1R[5]+0.4330127018922193*f1L[5]+0.25*(f1R[4]+f1L[4]+f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]+0.25*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alpha0R[5]*fUp[5]+alpha0R[4]*fUp[4]+alpha0R[3]*fUp[3]+alpha0R[2]*fUp[2]+alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[1] = 0.25*(alpha0R[3]*fUp[5]+fUp[3]*alpha0R[5]+alpha0R[2]*fUp[4]+fUp[2]*alpha0R[4]+alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1]); 
  incr[2] = -0.4330127018922193*(alpha0R[5]*fUp[5]+alpha0R[4]*fUp[4]+alpha0R[3]*fUp[3]+alpha0R[2]*fUp[2]+alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[3] = 0.25*(alpha0R[5]*fUp[7]+alpha0R[3]*fUp[6]+alpha0R[1]*fUp[4]+fUp[1]*alpha0R[4]+alpha0R[0]*fUp[2]+fUp[0]*alpha0R[2]); 
  incr[4] = 0.25*(alpha0R[4]*fUp[7]+alpha0R[2]*fUp[6]+alpha0R[1]*fUp[5]+fUp[1]*alpha0R[5]+alpha0R[0]*fUp[3]+fUp[0]*alpha0R[3]); 
  incr[5] = -0.4330127018922193*(alpha0R[3]*fUp[5]+fUp[3]*alpha0R[5]+alpha0R[2]*fUp[4]+fUp[2]*alpha0R[4]+alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1]); 
  incr[6] = 0.25*(alpha0R[3]*fUp[7]+alpha0R[5]*fUp[6]+alpha0R[0]*fUp[4]+fUp[0]*alpha0R[4]+alpha0R[1]*fUp[2]+fUp[1]*alpha0R[2]); 
  incr[7] = -0.4330127018922193*(alpha0R[5]*fUp[7]+alpha0R[3]*fUp[6]+alpha0R[1]*fUp[4]+fUp[1]*alpha0R[4]+alpha0R[0]*fUp[2]+fUp[0]*alpha0R[2]); 
  incr[8] = 0.25*(alpha0R[2]*fUp[7]+alpha0R[4]*fUp[6]+alpha0R[0]*fUp[5]+fUp[0]*alpha0R[5]+alpha0R[1]*fUp[3]+fUp[1]*alpha0R[3]); 
  incr[9] = -0.4330127018922193*(alpha0R[4]*fUp[7]+alpha0R[2]*fUp[6]+alpha0R[1]*fUp[5]+fUp[1]*alpha0R[5]+alpha0R[0]*fUp[3]+fUp[0]*alpha0R[3]); 
  incr[10] = 0.25*(alpha0R[1]*fUp[7]+alpha0R[0]*fUp[6]+alpha0R[4]*fUp[5]+fUp[4]*alpha0R[5]+alpha0R[2]*fUp[3]+fUp[2]*alpha0R[3]); 
  incr[11] = -0.4330127018922193*(alpha0R[3]*fUp[7]+alpha0R[5]*fUp[6]+alpha0R[0]*fUp[4]+fUp[0]*alpha0R[4]+alpha0R[1]*fUp[2]+fUp[1]*alpha0R[2]); 
  incr[12] = -0.4330127018922193*(alpha0R[2]*fUp[7]+alpha0R[4]*fUp[6]+alpha0R[0]*fUp[5]+fUp[0]*alpha0R[5]+alpha0R[1]*fUp[3]+fUp[1]*alpha0R[3]); 
  incr[13] = 0.25*(alpha0R[0]*fUp[7]+alpha0R[1]*fUp[6]+alpha0R[2]*fUp[5]+fUp[2]*alpha0R[5]+alpha0R[3]*fUp[4]+fUp[3]*alpha0R[4]); 
  incr[14] = -0.4330127018922193*(alpha0R[1]*fUp[7]+alpha0R[0]*fUp[6]+alpha0R[4]*fUp[5]+fUp[4]*alpha0R[5]+alpha0R[2]*fUp[3]+fUp[2]*alpha0R[3]); 
  incr[15] = -0.4330127018922193*(alpha0R[0]*fUp[7]+alpha0R[1]*fUp[6]+alpha0R[2]*fUp[5]+fUp[2]*alpha0R[5]+alpha0R[3]*fUp[4]+fUp[3]*alpha0R[4]); 

  outR[0] += incr[0]*rdy2R; 
  outR[1] += incr[1]*rdy2R; 
  outR[2] += incr[2]*rdy2R; 
  outR[3] += incr[3]*rdy2R; 
  outR[4] += incr[4]*rdy2R; 
  outR[5] += incr[5]*rdy2R; 
  outR[6] += incr[6]*rdy2R; 
  outR[7] += incr[7]*rdy2R; 
  outR[8] += incr[8]*rdy2R; 
  outR[9] += incr[9]*rdy2R; 
  outR[10] += incr[10]*rdy2R; 
  outR[11] += incr[11]*rdy2R; 
  outR[12] += incr[12]*rdy2R; 
  outR[13] += incr[13]*rdy2R; 
  outR[14] += incr[14]*rdy2R; 
  outR[15] += incr[15]*rdy2R; 

  outL[0] += -1.0*incr[0]*rdy2L; 
  outL[1] += -1.0*incr[1]*rdy2L; 
  outL[2] += incr[2]*rdy2L; 
  outL[3] += -1.0*incr[3]*rdy2L; 
  outL[4] += -1.0*incr[4]*rdy2L; 
  outL[5] += incr[5]*rdy2L; 
  outL[6] += -1.0*incr[6]*rdy2L; 
  outL[7] += incr[7]*rdy2L; 
  outL[8] += -1.0*incr[8]*rdy2L; 
  outL[9] += incr[9]*rdy2L; 
  outL[10] += -1.0*incr[10]*rdy2L; 
  outL[11] += incr[11]*rdy2L; 
  outL[12] += incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.4330127018922193*(f0R[14]+f0L[14])+0.25*f0R[13]-0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0R[3])-0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*(f0L[15]+f0R[14])+0.4330127018922193*f0L[14]-0.25*(f0R[13]+f0L[13])-0.4330127018922193*(f0R[12]+f0R[11])+0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0R[3])-0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f0R[15]+f0R[14])+0.4330127018922193*(f0L[15]+f0L[14])+0.25*(f0R[13]+f0L[13])+0.4330127018922193*(f0R[12]+f0R[11])-0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.4330127018922193*(f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])+0.25*f0R[4]-0.25*(f0L[4]+f0R[3])+0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*(f0L[15]+f0R[14])-0.4330127018922193*f0L[14]+0.25*(f0R[13]+f0L[13])-0.4330127018922193*f0R[12]+0.4330127018922193*(f0L[12]+f0R[11])-0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4])+0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]))+0.25*f0R[13]-0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])+0.25*f0R[4]-0.25*(f0L[4]+f0R[3])+0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f0R[15]+f0R[14])-0.4330127018922193*(f0L[15]+f0L[14])-0.25*(f0R[13]+f0L[13])+0.4330127018922193*f0R[12]-0.4330127018922193*(f0L[12]+f0R[11])+0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4])+0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.4330127018922193*(f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])-0.25*f0R[4]+0.25*(f0L[4]+f0R[3])-0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*(f0L[15]+f0R[14])-0.4330127018922193*f0L[14]+0.25*(f0R[13]+f0L[13])+0.4330127018922193*f0R[12]-0.4330127018922193*(f0L[12]+f0R[11])+0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4])-0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]))+0.25*f0R[13]-0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])-0.25*f0R[4]+0.25*(f0L[4]+f0R[3])-0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f0R[15]+f0R[14])-0.4330127018922193*(f0L[15]+f0L[14])-0.25*(f0R[13]+f0L[13])-0.4330127018922193*f0R[12]+0.4330127018922193*(f0L[12]+f0R[11])-0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4])-0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.4330127018922193*(f0R[14]+f0L[14])+0.25*f0R[13]-0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0R[3])+0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*(f0L[15]+f0R[14])+0.4330127018922193*f0L[14]-0.25*(f0R[13]+f0L[13])+0.4330127018922193*(f0R[12]+f0R[11])-0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0R[3])+0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f0R[15]+f0R[14])+0.4330127018922193*(f0L[15]+f0L[14])+0.25*(f0R[13]+f0L[13])-0.4330127018922193*(f0R[12]+f0R[11])+0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.25*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.4330127018922193*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = 0.25*(alpha1R[1]*fUp[4]+alpha1R[0]*fUp[2]); 
  incr[4] = 0.25*(alpha1R[1]*fUp[5]+alpha1R[0]*fUp[3]); 
  incr[5] = -0.4330127018922193*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[6] = 0.25*(alpha1R[0]*fUp[4]+alpha1R[1]*fUp[2]); 
  incr[7] = -0.4330127018922193*(alpha1R[1]*fUp[4]+alpha1R[0]*fUp[2]); 
  incr[8] = 0.25*(alpha1R[0]*fUp[5]+alpha1R[1]*fUp[3]); 
  incr[9] = -0.4330127018922193*(alpha1R[1]*fUp[5]+alpha1R[0]*fUp[3]); 
  incr[10] = 0.25*(alpha1R[1]*fUp[7]+alpha1R[0]*fUp[6]); 
  incr[11] = -0.4330127018922193*(alpha1R[0]*fUp[4]+alpha1R[1]*fUp[2]); 
  incr[12] = -0.4330127018922193*(alpha1R[0]*fUp[5]+alpha1R[1]*fUp[3]); 
  incr[13] = 0.25*(alpha1R[0]*fUp[7]+alpha1R[1]*fUp[6]); 
  incr[14] = -0.4330127018922193*(alpha1R[1]*fUp[7]+alpha1R[0]*fUp[6]); 
  incr[15] = -0.4330127018922193*(alpha1R[0]*fUp[7]+alpha1R[1]*fUp[6]); 

  outR[0] += incr[0]*rdy2R; 
  outR[1] += incr[1]*rdy2R; 
  outR[2] += incr[2]*rdy2R; 
  outR[3] += incr[3]*rdy2R; 
  outR[4] += incr[4]*rdy2R; 
  outR[5] += incr[5]*rdy2R; 
  outR[6] += incr[6]*rdy2R; 
  outR[7] += incr[7]*rdy2R; 
  outR[8] += incr[8]*rdy2R; 
  outR[9] += incr[9]*rdy2R; 
  outR[10] += incr[10]*rdy2R; 
  outR[11] += incr[11]*rdy2R; 
  outR[12] += incr[12]*rdy2R; 
  outR[13] += incr[13]*rdy2R; 
  outR[14] += incr[14]*rdy2R; 
  outR[15] += incr[15]*rdy2R; 

  outL[0] += -1.0*incr[0]*rdy2L; 
  outL[1] += -1.0*incr[1]*rdy2L; 
  outL[2] += incr[2]*rdy2L; 
  outL[3] += -1.0*incr[3]*rdy2L; 
  outL[4] += -1.0*incr[4]*rdy2L; 
  outL[5] += incr[5]*rdy2L; 
  outL[6] += -1.0*incr[6]*rdy2L; 
  outL[7] += incr[7]*rdy2L; 
  outL[8] += -1.0*incr[8]*rdy2L; 
  outL[9] += incr[9]*rdy2L; 
  outL[10] += -1.0*incr[10]*rdy2L; 
  outL[11] += incr[11]*rdy2L; 
  outL[12] += incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 

  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[8] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagR[16]; 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = -(1.732050807568877*jacobTotInv[0]*b_z[1]*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[1] = -(1.732050807568877*b_z[1]*jacobTotInv[1]*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alpha0R[8]; 

  double alpha1R[8]; 
  alpha1R[0] = (0.125*(hamil1R[5]*(4.242640687119286*BstarYdBmagR[6]-2.449489742783178*BstarYdBmagR[1])+hamil1R[2]*(4.242640687119286*BstarYdBmagR[3]-2.449489742783178*BstarYdBmagR[0]))*rdy2R)/m_; 
  alpha1R[1] = (0.125*(4.242640687119286*hamil1R[2]*BstarYdBmagR[6]+(4.242640687119286*BstarYdBmagR[3]-2.449489742783178*BstarYdBmagR[0])*hamil1R[5]-2.449489742783178*BstarYdBmagR[1]*hamil1R[2])*rdy2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*(3.0*hamil1R[5]*BstarYdBmagR[6]-1.732050807568877*BstarYdBmagR[1]*hamil1R[5]+3.0*hamil1R[2]*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0]*hamil1R[2])*rdy2R)/m_; 

  double incr[16]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.4330127018922193*(f0R[14]+f0L[14]+f0R[13]+f0L[13])+0.25*f0R[12]-0.25*f0L[12]+0.4330127018922193*(f0R[11]+f0L[11])-0.4330127018922193*(f0R[10]+f0L[10])-0.25*(f0R[9]+f0R[8])+0.25*(f0L[9]+f0L[8])-0.4330127018922193*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.25*f0R[5]+0.25*(f0L[5]+f0R[4])-0.25*f0L[4]+0.4330127018922193*(f0R[3]+f0L[3])+0.25*(f0R[2]+f0R[1])-0.25*(f0L[2]+f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*(f0L[15]+f0R[14]+f0R[13])+0.4330127018922193*(f0L[14]+f0L[13])-0.25*(f0R[12]+f0L[12])-0.4330127018922193*f0R[11]+0.4330127018922193*(f0L[11]+f0R[10])-0.4330127018922193*f0L[10]+0.25*(f0R[9]+f0L[9]+f0R[8]+f0L[8])+0.4330127018922193*(f0R[7]+f0R[6])-0.4330127018922193*(f0L[7]+f0L[6])+0.25*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]-0.25*(f0R[2]+f0L[2]+f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14])-0.4330127018922193*(f0R[13]+f0L[13])-0.25*f0R[12]+0.25*f0L[12]-0.4330127018922193*(f0R[11]+f0L[11]+f0R[10]+f0L[10])-0.25*f0R[9]+0.25*(f0L[9]+f0R[8])-0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])+0.4330127018922193*(f0R[6]+f0L[6])+0.25*(f0R[5]+f0R[4])-0.25*(f0L[5]+f0L[4])+0.4330127018922193*(f0R[3]+f0L[3])+0.25*f0R[2]-0.25*(f0L[2]+f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f0R[15]+f0R[14])+0.4330127018922193*(f0L[15]+f0L[14]+f0R[13])-0.4330127018922193*f0L[13]+0.25*(f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0R[10])-0.4330127018922193*(f0L[11]+f0L[10])+0.25*(f0R[9]+f0L[9])-0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*(f0L[7]+f0R[6])+0.4330127018922193*f0L[6]-0.25*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]-0.25*(f0R[2]+f0L[2])+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.4330127018922193*(f0R[14]+f0L[14])+0.4330127018922193*(f0R[13]+f0L[13])-0.25*f0R[12]+0.25*f0L[12]-0.4330127018922193*(f0R[11]+f0L[11]+f0R[10]+f0L[10])+0.25*f0R[9]-0.25*(f0L[9]+f0R[8])+0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])-0.4330127018922193*(f0R[6]+f0L[6])+0.25*(f0R[5]+f0R[4])-0.25*(f0L[5]+f0L[4])+0.4330127018922193*(f0R[3]+f0L[3])-0.25*f0R[2]+0.25*(f0L[2]+f0R[1])-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*(f0L[15]+f0R[14])-0.4330127018922193*(f0L[14]+f0R[13])+0.4330127018922193*f0L[13]+0.25*(f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0R[10])-0.4330127018922193*(f0L[11]+f0L[10])-0.25*(f0R[9]+f0L[9])+0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*(f0L[7]+f0R[6])-0.4330127018922193*f0L[6]-0.25*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]+0.25*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]+f0R[13]+f0L[13]))+0.25*f0R[12]-0.25*f0L[12]+0.4330127018922193*(f0R[11]+f0L[11])-0.4330127018922193*(f0R[10]+f0L[10])+0.25*(f0R[9]+f0R[8])-0.25*(f0L[9]+f0L[8])+0.4330127018922193*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.25*f0R[5]+0.25*(f0L[5]+f0R[4])-0.25*f0L[4]+0.4330127018922193*(f0R[3]+f0L[3])-0.25*(f0R[2]+f0R[1]+f0R[0])+0.25*(f0L[2]+f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f0R[15]+f0R[14]+f0R[13])-0.4330127018922193*(f0L[15]+f0L[14]+f0L[13])-0.25*(f0R[12]+f0L[12])-0.4330127018922193*f0R[11]+0.4330127018922193*(f0L[11]+f0R[10])-0.4330127018922193*f0L[10]-0.25*(f0R[9]+f0L[9]+f0R[8]+f0L[8])-0.4330127018922193*(f0R[7]+f0R[6])+0.4330127018922193*(f0L[7]+f0L[6])+0.25*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]+0.25*(f0R[2]+f0L[2]+f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.4330127018922193*(f0R[14]+f0L[14]+f0R[13]+f0L[13])-0.25*f0R[12]+0.25*f0L[12]+0.4330127018922193*(f0R[11]+f0L[11]+f0R[10]+f0L[10])+0.25*(f0R[9]+f0R[8])-0.25*(f0L[9]+f0L[8])-0.4330127018922193*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.25*(f0R[5]+f0R[4])+0.25*(f0L[5]+f0L[4])+0.4330127018922193*(f0R[3]+f0L[3])+0.25*(f0R[2]+f0R[1])-0.25*(f0L[2]+f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*(f0L[15]+f0R[14]+f0R[13])-0.4330127018922193*(f0L[14]+f0L[13])+0.25*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0R[10])+0.4330127018922193*(f0L[11]+f0L[10])-0.25*(f0R[9]+f0L[9]+f0R[8]+f0L[8])+0.4330127018922193*(f0R[7]+f0R[6])-0.4330127018922193*(f0L[7]+f0L[6])+0.25*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]-0.25*(f0R[2]+f0L[2]+f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]))+0.4330127018922193*(f0R[13]+f0L[13])+0.25*f0R[12]-0.25*f0L[12]-0.4330127018922193*(f0R[11]+f0L[11])+0.4330127018922193*(f0R[10]+f0L[10])+0.25*f0R[9]-0.25*(f0L[9]+f0R[8])+0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])+0.4330127018922193*(f0R[6]+f0L[6])+0.25*f0R[5]-0.25*(f0L[5]+f0R[4])+0.25*f0L[4]+0.4330127018922193*(f0R[3]+f0L[3])+0.25*f0R[2]-0.25*(f0L[2]+f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f0R[15]+f0R[14])-0.4330127018922193*(f0L[15]+f0L[14]+f0R[13])+0.4330127018922193*f0L[13]-0.25*(f0R[12]+f0L[12])+0.4330127018922193*f0R[11]-0.4330127018922193*(f0L[11]+f0R[10])+0.4330127018922193*f0L[10]-0.25*(f0R[9]+f0L[9])+0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*(f0L[7]+f0R[6])+0.4330127018922193*f0L[6]-0.25*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]-0.25*(f0R[2]+f0L[2])+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.4330127018922193*(f0R[14]+f0L[14])-0.4330127018922193*(f0R[13]+f0L[13])+0.25*f0R[12]-0.25*f0L[12]-0.4330127018922193*(f0R[11]+f0L[11])+0.4330127018922193*(f0R[10]+f0L[10])-0.25*f0R[9]+0.25*(f0L[9]+f0R[8])-0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])-0.4330127018922193*(f0R[6]+f0L[6])+0.25*f0R[5]-0.25*(f0L[5]+f0R[4])+0.25*f0L[4]+0.4330127018922193*(f0R[3]+f0L[3])-0.25*f0R[2]+0.25*(f0L[2]+f0R[1])-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*(f0L[15]+f0R[14])+0.4330127018922193*(f0L[14]+f0R[13])-0.4330127018922193*f0L[13]-0.25*(f0R[12]+f0L[12])+0.4330127018922193*f0R[11]-0.4330127018922193*(f0L[11]+f0R[10])+0.4330127018922193*f0L[10]+0.25*(f0R[9]+f0L[9])-0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*(f0L[7]+f0R[6])-0.4330127018922193*f0L[6]-0.25*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]+0.25*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]+f0R[13]+f0L[13])-0.25*f0R[12]+0.25*f0L[12]+0.4330127018922193*(f0R[11]+f0L[11]+f0R[10]+f0L[10])-0.25*(f0R[9]+f0R[8])+0.25*(f0L[9]+f0L[8])+0.4330127018922193*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.25*(f0R[5]+f0R[4])+0.25*(f0L[5]+f0L[4])+0.4330127018922193*(f0R[3]+f0L[3])-0.25*(f0R[2]+f0R[1]+f0R[0])+0.25*(f0L[2]+f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f0R[15]+f0R[14]+f0R[13])+0.4330127018922193*(f0L[15]+f0L[14]+f0L[13])+0.25*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0R[10])+0.4330127018922193*(f0L[11]+f0L[10])+0.25*(f0R[9]+f0L[9]+f0R[8]+f0L[8])-0.4330127018922193*(f0R[7]+f0R[6])+0.4330127018922193*(f0L[7]+f0L[6])+0.25*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]+0.25*(f0R[2]+f0L[2]+f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.25*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = 0.25*(alpha1R[1]*fUp[4]+alpha1R[0]*fUp[2]); 
  incr[3] = -0.4330127018922193*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[4] = 0.25*(alpha1R[1]*fUp[5]+alpha1R[0]*fUp[3]); 
  incr[5] = 0.25*(alpha1R[0]*fUp[4]+alpha1R[1]*fUp[2]); 
  incr[6] = -0.4330127018922193*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[7] = -0.4330127018922193*(alpha1R[1]*fUp[4]+alpha1R[0]*fUp[2]); 
  incr[8] = 0.25*(alpha1R[0]*fUp[5]+alpha1R[1]*fUp[3]); 
  incr[9] = 0.25*(alpha1R[1]*fUp[7]+alpha1R[0]*fUp[6]); 
  incr[10] = -0.4330127018922193*(alpha1R[1]*fUp[5]+alpha1R[0]*fUp[3]); 
  incr[11] = -0.4330127018922193*(alpha1R[0]*fUp[4]+alpha1R[1]*fUp[2]); 
  incr[12] = 0.25*(alpha1R[0]*fUp[7]+alpha1R[1]*fUp[6]); 
  incr[13] = -0.4330127018922193*(alpha1R[0]*fUp[5]+alpha1R[1]*fUp[3]); 
  incr[14] = -0.4330127018922193*(alpha1R[1]*fUp[7]+alpha1R[0]*fUp[6]); 
  incr[15] = -0.4330127018922193*(alpha1R[0]*fUp[7]+alpha1R[1]*fUp[6]); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 
  outR[8] += incr[8]*rdvpar2R; 
  outR[9] += incr[9]*rdvpar2R; 
  outR[10] += incr[10]*rdvpar2R; 
  outR[11] += incr[11]*rdvpar2R; 
  outR[12] += incr[12]*rdvpar2R; 
  outR[13] += incr[13]*rdvpar2R; 
  outR[14] += incr[14]*rdvpar2R; 
  outR[15] += incr[15]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += -1.0*incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 
  outL[4] += -1.0*incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 
  outL[8] += -1.0*incr[8]*rdvpar2L; 
  outL[9] += -1.0*incr[9]*rdvpar2L; 
  outL[10] += incr[10]*rdvpar2L; 
  outL[11] += incr[11]*rdvpar2L; 
  outL[12] += -1.0*incr[12]*rdvpar2L; 
  outL[13] += incr[13]*rdvpar2L; 
  outL[14] += incr[14]*rdvpar2L; 
  outL[15] += incr[15]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_x_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagR[16]; 

  double alpha0R[8]; 

  double alpha1R[8]; 
  alpha1R[0] = (0.1767766952966368*b_z[0]*jacobTotInv[0]*(3.0*hamil1R[5]-1.732050807568877*hamil1R[2])*rdy2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*(3.0*b_z[0]*jacobTotInv[0]*hamil1R[5]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamil1R[2])*rdy2R)/q_; 

  double incr[16]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.25*f0R[14]-0.25*f0L[14]+0.4330127018922193*(f0R[13]+f0L[13]+f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*(f0R[10]+f0R[9])+0.25*(f0L[10]+f0L[9])-0.4330127018922193*(f0R[8]+f0L[8])-0.25*f0R[7]+0.25*f0L[7]-0.4330127018922193*(f0R[6]+f0L[6]+f0R[5]+f0L[5])+0.25*(f0R[4]+f0R[3]+f0R[2])-0.25*(f0L[4]+f0L[3]+f0L[2])+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*f0L[15]-0.25*(f0R[14]+f0L[14])-0.4330127018922193*(f0R[13]+f0R[12]+f0R[11])+0.4330127018922193*(f0L[13]+f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10]+f0R[9]+f0L[9])+0.4330127018922193*f0R[8]-0.4330127018922193*f0L[8]+0.25*(f0R[7]+f0L[7])+0.4330127018922193*(f0R[6]+f0R[5])-0.4330127018922193*(f0L[6]+f0L[5])-0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.25*f0R[14]+0.25*f0L[14]+0.4330127018922193*(f0R[13]+f0L[13])-0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*(f0L[10]+f0R[9])-0.25*f0L[9]-0.4330127018922193*(f0R[8]+f0L[8])+0.25*f0R[7]-0.25*f0L[7]-0.4330127018922193*(f0R[6]+f0L[6])+0.4330127018922193*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0R[3])-0.25*(f0L[4]+f0L[3]+f0R[2])+0.25*f0L[2]+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*f0L[15]+0.25*(f0R[14]+f0L[14])-0.4330127018922193*f0R[13]+0.4330127018922193*(f0L[13]+f0R[12]+f0R[11])-0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])-0.25*(f0R[9]+f0L[9])+0.4330127018922193*f0R[8]-0.4330127018922193*f0L[8]-0.25*(f0R[7]+f0L[7])+0.4330127018922193*f0R[6]-0.4330127018922193*(f0L[6]+f0R[5])+0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])+0.25*(f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.25*f0R[14]+0.25*f0L[14]-0.4330127018922193*(f0R[13]+f0L[13])+0.4330127018922193*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*(f0L[10]+f0R[9])+0.25*f0L[9]-0.4330127018922193*(f0R[8]+f0L[8])+0.25*f0R[7]-0.25*f0L[7]+0.4330127018922193*(f0R[6]+f0L[6])-0.4330127018922193*(f0R[5]+f0L[5])+0.25*f0R[4]-0.25*(f0L[4]+f0R[3])+0.25*(f0L[3]+f0R[2])-0.25*f0L[2]+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*f0L[15]+0.25*(f0R[14]+f0L[14])+0.4330127018922193*f0R[13]-0.4330127018922193*(f0L[13]+f0R[12])+0.4330127018922193*(f0L[12]+f0R[11])-0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])+0.25*(f0R[9]+f0L[9])+0.4330127018922193*f0R[8]-0.4330127018922193*f0L[8]-0.25*(f0R[7]+f0L[7])-0.4330127018922193*f0R[6]+0.4330127018922193*(f0L[6]+f0R[5])-0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4])+0.25*(f0R[3]+f0L[3])-0.25*(f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.25*f0R[14]-0.25*f0L[14]-0.4330127018922193*(f0R[13]+f0L[13]+f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0L[11])+0.25*(f0R[10]+f0R[9])-0.25*(f0L[10]+f0L[9])-0.4330127018922193*(f0R[8]+f0L[8])-0.25*f0R[7]+0.25*f0L[7]+0.4330127018922193*(f0R[6]+f0L[6]+f0R[5]+f0L[5])+0.25*f0R[4]-0.25*(f0L[4]+f0R[3]+f0R[2])+0.25*(f0L[3]+f0L[2])+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*f0L[15]-0.25*(f0R[14]+f0L[14])+0.4330127018922193*(f0R[13]+f0R[12])-0.4330127018922193*(f0L[13]+f0L[12]+f0R[11])+0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10]+f0R[9]+f0L[9])+0.4330127018922193*f0R[8]-0.4330127018922193*f0L[8]+0.25*(f0R[7]+f0L[7])-0.4330127018922193*(f0R[6]+f0R[5])+0.4330127018922193*(f0L[6]+f0L[5])-0.25*(f0R[4]+f0L[4])+0.25*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.25*f0R[14]+0.25*f0L[14]-0.4330127018922193*(f0R[13]+f0L[13]+f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0L[11])+0.25*(f0R[10]+f0R[9])-0.25*(f0L[10]+f0L[9])+0.4330127018922193*(f0R[8]+f0L[8])-0.25*f0R[7]+0.25*f0L[7]-0.4330127018922193*(f0R[6]+f0L[6]+f0R[5]+f0L[5])-0.25*f0R[4]+0.25*(f0L[4]+f0R[3]+f0R[2])-0.25*(f0L[3]+f0L[2])+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*f0L[15]+0.25*(f0R[14]+f0L[14])+0.4330127018922193*(f0R[13]+f0R[12])-0.4330127018922193*(f0L[13]+f0L[12]+f0R[11])+0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10]+f0R[9]+f0L[9])-0.4330127018922193*f0R[8]+0.4330127018922193*f0L[8]+0.25*(f0R[7]+f0L[7])+0.4330127018922193*(f0R[6]+f0R[5])-0.4330127018922193*(f0L[6]+f0L[5])+0.25*(f0R[4]+f0L[4])-0.25*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.25*f0R[14]-0.25*f0L[14]-0.4330127018922193*(f0R[13]+f0L[13])+0.4330127018922193*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*(f0L[10]+f0R[9])+0.25*f0L[9]+0.4330127018922193*(f0R[8]+f0L[8])+0.25*f0R[7]-0.25*f0L[7]-0.4330127018922193*(f0R[6]+f0L[6])+0.4330127018922193*(f0R[5]+f0L[5])-0.25*f0R[4]+0.25*(f0L[4]+f0R[3])-0.25*(f0L[3]+f0R[2])+0.25*f0L[2]+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*f0L[15]-0.25*(f0R[14]+f0L[14])+0.4330127018922193*f0R[13]-0.4330127018922193*(f0L[13]+f0R[12])+0.4330127018922193*(f0L[12]+f0R[11])-0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])+0.25*(f0R[9]+f0L[9])-0.4330127018922193*f0R[8]+0.4330127018922193*f0L[8]-0.25*(f0R[7]+f0L[7])+0.4330127018922193*f0R[6]-0.4330127018922193*(f0L[6]+f0R[5])+0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4])-0.25*(f0R[3]+f0L[3])+0.25*(f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.25*f0R[14]-0.25*f0L[14]+0.4330127018922193*(f0R[13]+f0L[13])-0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*(f0L[10]+f0R[9])-0.25*f0L[9]+0.4330127018922193*(f0R[8]+f0L[8])+0.25*f0R[7]-0.25*f0L[7]+0.4330127018922193*(f0R[6]+f0L[6])-0.4330127018922193*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0R[3])+0.25*(f0L[4]+f0L[3]+f0R[2])-0.25*f0L[2]+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*f0L[15]-0.25*(f0R[14]+f0L[14])-0.4330127018922193*f0R[13]+0.4330127018922193*(f0L[13]+f0R[12]+f0R[11])-0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])-0.25*(f0R[9]+f0L[9])-0.4330127018922193*f0R[8]+0.4330127018922193*f0L[8]-0.25*(f0R[7]+f0L[7])-0.4330127018922193*f0R[6]+0.4330127018922193*(f0L[6]+f0R[5])-0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.25*(f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[7] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.25*f0R[14]+0.25*f0L[14]+0.4330127018922193*(f0R[13]+f0L[13]+f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*(f0R[10]+f0R[9])+0.25*(f0L[10]+f0L[9])+0.4330127018922193*(f0R[8]+f0L[8])-0.25*f0R[7]+0.25*f0L[7]+0.4330127018922193*(f0R[6]+f0L[6]+f0R[5]+f0L[5])-0.25*(f0R[4]+f0R[3]+f0R[2])+0.25*(f0L[4]+f0L[3]+f0L[2])+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*f0L[15]+0.25*(f0R[14]+f0L[14])-0.4330127018922193*(f0R[13]+f0R[12]+f0R[11])+0.4330127018922193*(f0L[13]+f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10]+f0R[9]+f0L[9])-0.4330127018922193*f0R[8]+0.4330127018922193*f0L[8]+0.25*(f0R[7]+f0L[7])-0.4330127018922193*(f0R[6]+f0R[5])+0.4330127018922193*(f0L[6]+f0L[5])+0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*alpha1R[0]*fUp[0]; 
  incr[1] = -0.4330127018922193*alpha1R[0]*fUp[0]; 
  incr[2] = 0.25*alpha1R[0]*fUp[1]; 
  incr[3] = 0.25*alpha1R[0]*fUp[2]; 
  incr[4] = 0.25*alpha1R[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alpha1R[0]*fUp[1]; 
  incr[6] = -0.4330127018922193*alpha1R[0]*fUp[2]; 
  incr[7] = 0.25*alpha1R[0]*fUp[4]; 
  incr[8] = -0.4330127018922193*alpha1R[0]*fUp[3]; 
  incr[9] = 0.25*alpha1R[0]*fUp[5]; 
  incr[10] = 0.25*alpha1R[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alpha1R[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alpha1R[0]*fUp[5]; 
  incr[13] = -0.4330127018922193*alpha1R[0]*fUp[6]; 
  incr[14] = 0.25*alpha1R[0]*fUp[7]; 
  incr[15] = -0.4330127018922193*alpha1R[0]*fUp[7]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 
  outR[8] += incr[8]*rdx2R; 
  outR[9] += incr[9]*rdx2R; 
  outR[10] += incr[10]*rdx2R; 
  outR[11] += incr[11]*rdx2R; 
  outR[12] += incr[12]*rdx2R; 
  outR[13] += incr[13]*rdx2R; 
  outR[14] += incr[14]*rdx2R; 
  outR[15] += incr[15]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += -1.0*incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += incr[6]*rdx2L; 
  outL[7] += -1.0*incr[7]*rdx2L; 
  outL[8] += incr[8]*rdx2L; 
  outL[9] += -1.0*incr[9]*rdx2L; 
  outL[10] += -1.0*incr[10]*rdx2L; 
  outL[11] += incr[11]*rdx2L; 
  outL[12] += incr[12]*rdx2L; 
  outL[13] += incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_y_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarYdBmagR[16]; 

  double alpha0R[8]; 

  double alpha1R[8]; 
  alpha1R[0] = -(0.1767766952966368*b_z[0]*jacobTotInv[0]*(3.0*hamil1R[5]-1.732050807568877*hamil1R[1])*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(3.0*b_z[0]*jacobTotInv[0]*hamil1R[5]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamil1R[1])*rdx2R)/q_; 

  double incr[16]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.4330127018922193*(f0R[14]+f0L[14])+0.25*f0R[13]-0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0R[3])-0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*(f0L[15]+f0R[14])+0.4330127018922193*f0L[14]-0.25*(f0R[13]+f0L[13])-0.4330127018922193*(f0R[12]+f0R[11])+0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0R[3])-0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f0R[15]+f0R[14])+0.4330127018922193*(f0L[15]+f0L[14])+0.25*(f0R[13]+f0L[13])+0.4330127018922193*(f0R[12]+f0R[11])-0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.4330127018922193*(f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])+0.25*f0R[4]-0.25*(f0L[4]+f0R[3])+0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*(f0L[15]+f0R[14])-0.4330127018922193*f0L[14]+0.25*(f0R[13]+f0L[13])-0.4330127018922193*f0R[12]+0.4330127018922193*(f0L[12]+f0R[11])-0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4])+0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]))+0.25*f0R[13]-0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])+0.25*f0R[4]-0.25*(f0L[4]+f0R[3])+0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f0R[15]+f0R[14])-0.4330127018922193*(f0L[15]+f0L[14])-0.25*(f0R[13]+f0L[13])+0.4330127018922193*f0R[12]-0.4330127018922193*(f0L[12]+f0R[11])+0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4])+0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.4330127018922193*(f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])-0.25*f0R[4]+0.25*(f0L[4]+f0R[3])-0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*(f0L[15]+f0R[14])-0.4330127018922193*f0L[14]+0.25*(f0R[13]+f0L[13])+0.4330127018922193*f0R[12]-0.4330127018922193*(f0L[12]+f0R[11])+0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4])-0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]))+0.25*f0R[13]-0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])-0.25*f0R[4]+0.25*(f0L[4]+f0R[3])-0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f0R[15]+f0R[14])-0.4330127018922193*(f0L[15]+f0L[14])-0.25*(f0R[13]+f0L[13])-0.4330127018922193*f0R[12]+0.4330127018922193*(f0L[12]+f0R[11])-0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4])-0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.4330127018922193*(f0R[14]+f0L[14])+0.25*f0R[13]-0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0R[3])+0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*(f0L[15]+f0R[14])+0.4330127018922193*f0L[14]-0.25*(f0R[13]+f0L[13])+0.4330127018922193*(f0R[12]+f0R[11])-0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[7] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0R[3])+0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f0R[15]+f0R[14])+0.4330127018922193*(f0L[15]+f0L[14])+0.25*(f0R[13]+f0L[13])-0.4330127018922193*(f0R[12]+f0R[11])+0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*alpha1R[0]*fUp[0]; 
  incr[1] = 0.25*alpha1R[0]*fUp[1]; 
  incr[2] = -0.4330127018922193*alpha1R[0]*fUp[0]; 
  incr[3] = 0.25*alpha1R[0]*fUp[2]; 
  incr[4] = 0.25*alpha1R[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alpha1R[0]*fUp[1]; 
  incr[6] = 0.25*alpha1R[0]*fUp[4]; 
  incr[7] = -0.4330127018922193*alpha1R[0]*fUp[2]; 
  incr[8] = 0.25*alpha1R[0]*fUp[5]; 
  incr[9] = -0.4330127018922193*alpha1R[0]*fUp[3]; 
  incr[10] = 0.25*alpha1R[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alpha1R[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alpha1R[0]*fUp[5]; 
  incr[13] = 0.25*alpha1R[0]*fUp[7]; 
  incr[14] = -0.4330127018922193*alpha1R[0]*fUp[6]; 
  incr[15] = -0.4330127018922193*alpha1R[0]*fUp[7]; 

  outR[0] += incr[0]*rdy2R; 
  outR[1] += incr[1]*rdy2R; 
  outR[2] += incr[2]*rdy2R; 
  outR[3] += incr[3]*rdy2R; 
  outR[4] += incr[4]*rdy2R; 
  outR[5] += incr[5]*rdy2R; 
  outR[6] += incr[6]*rdy2R; 
  outR[7] += incr[7]*rdy2R; 
  outR[8] += incr[8]*rdy2R; 
  outR[9] += incr[9]*rdy2R; 
  outR[10] += incr[10]*rdy2R; 
  outR[11] += incr[11]*rdy2R; 
  outR[12] += incr[12]*rdy2R; 
  outR[13] += incr[13]*rdy2R; 
  outR[14] += incr[14]*rdy2R; 
  outR[15] += incr[15]*rdy2R; 

  outL[0] += -1.0*incr[0]*rdy2L; 
  outL[1] += -1.0*incr[1]*rdy2L; 
  outL[2] += incr[2]*rdy2L; 
  outL[3] += -1.0*incr[3]*rdy2L; 
  outL[4] += -1.0*incr[4]*rdy2L; 
  outL[5] += incr[5]*rdy2L; 
  outL[6] += -1.0*incr[6]*rdy2L; 
  outL[7] += incr[7]*rdy2L; 
  outL[8] += -1.0*incr[8]*rdy2L; 
  outL[9] += incr[9]*rdy2L; 
  outL[10] += -1.0*incr[10]*rdy2L; 
  outL[11] += incr[11]*rdy2L; 
  outL[12] += incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 

  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_vpar_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagR[16]; 

  double BstarYdBmagR[16]; 

  double alpha0R[8]; 

  double alpha1R[8]; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = 0.0; 

  double incr[16]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // alpha1 == 0, so nothing to do for alpha1*f0 and alpha1*f1 terms.
  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_x_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[8] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagR[16]; 

  double alpha0R[8]; 

  double alpha1R[8]; 
  alpha1R[0] = (0.125*(((12.72792206135786*b_z[1]-7.348469228349534*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(4.242640687119286*b_z[0]-7.348469228349534*b_z[1]))*hamil1R[5]+((4.242640687119286*b_z[0]-7.348469228349534*b_z[1])*jacobTotInv[1]+jacobTotInv[0]*(4.242640687119286*b_z[1]-2.449489742783178*b_z[0]))*hamil1R[2])*rdy2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*(((9.0*b_z[1]-5.196152422706631*b_z[0])*jacobTotInv[1]-5.196152422706631*jacobTotInv[0]*b_z[1]+3.0*b_z[0]*jacobTotInv[0])*hamil1R[5]+((3.0*b_z[0]-5.196152422706631*b_z[1])*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1]-1.732050807568877*b_z[0]*jacobTotInv[0])*hamil1R[2])*rdy2R)/q_; 

  double incr[16]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.25*f0R[14]-0.25*f0L[14]+0.4330127018922193*(f0R[13]+f0L[13]+f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*(f0R[10]+f0R[9])+0.25*(f0L[10]+f0L[9])-0.4330127018922193*(f0R[8]+f0L[8])-0.25*f0R[7]+0.25*f0L[7]-0.4330127018922193*(f0R[6]+f0L[6]+f0R[5]+f0L[5])+0.25*(f0R[4]+f0R[3]+f0R[2])-0.25*(f0L[4]+f0L[3]+f0L[2])+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*f0L[15]-0.25*(f0R[14]+f0L[14])-0.4330127018922193*(f0R[13]+f0R[12]+f0R[11])+0.4330127018922193*(f0L[13]+f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10]+f0R[9]+f0L[9])+0.4330127018922193*f0R[8]-0.4330127018922193*f0L[8]+0.25*(f0R[7]+f0L[7])+0.4330127018922193*(f0R[6]+f0R[5])-0.4330127018922193*(f0L[6]+f0L[5])-0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.25*f0R[14]+0.25*f0L[14]+0.4330127018922193*(f0R[13]+f0L[13])-0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*(f0L[10]+f0R[9])-0.25*f0L[9]-0.4330127018922193*(f0R[8]+f0L[8])+0.25*f0R[7]-0.25*f0L[7]-0.4330127018922193*(f0R[6]+f0L[6])+0.4330127018922193*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0R[3])-0.25*(f0L[4]+f0L[3]+f0R[2])+0.25*f0L[2]+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*f0L[15]+0.25*(f0R[14]+f0L[14])-0.4330127018922193*f0R[13]+0.4330127018922193*(f0L[13]+f0R[12]+f0R[11])-0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])-0.25*(f0R[9]+f0L[9])+0.4330127018922193*f0R[8]-0.4330127018922193*f0L[8]-0.25*(f0R[7]+f0L[7])+0.4330127018922193*f0R[6]-0.4330127018922193*(f0L[6]+f0R[5])+0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])+0.25*(f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.25*f0R[14]+0.25*f0L[14]-0.4330127018922193*(f0R[13]+f0L[13])+0.4330127018922193*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*(f0L[10]+f0R[9])+0.25*f0L[9]-0.4330127018922193*(f0R[8]+f0L[8])+0.25*f0R[7]-0.25*f0L[7]+0.4330127018922193*(f0R[6]+f0L[6])-0.4330127018922193*(f0R[5]+f0L[5])+0.25*f0R[4]-0.25*(f0L[4]+f0R[3])+0.25*(f0L[3]+f0R[2])-0.25*f0L[2]+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*f0L[15]+0.25*(f0R[14]+f0L[14])+0.4330127018922193*f0R[13]-0.4330127018922193*(f0L[13]+f0R[12])+0.4330127018922193*(f0L[12]+f0R[11])-0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])+0.25*(f0R[9]+f0L[9])+0.4330127018922193*f0R[8]-0.4330127018922193*f0L[8]-0.25*(f0R[7]+f0L[7])-0.4330127018922193*f0R[6]+0.4330127018922193*(f0L[6]+f0R[5])-0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4])+0.25*(f0R[3]+f0L[3])-0.25*(f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.25*f0R[14]-0.25*f0L[14]-0.4330127018922193*(f0R[13]+f0L[13]+f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0L[11])+0.25*(f0R[10]+f0R[9])-0.25*(f0L[10]+f0L[9])-0.4330127018922193*(f0R[8]+f0L[8])-0.25*f0R[7]+0.25*f0L[7]+0.4330127018922193*(f0R[6]+f0L[6]+f0R[5]+f0L[5])+0.25*f0R[4]-0.25*(f0L[4]+f0R[3]+f0R[2])+0.25*(f0L[3]+f0L[2])+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*f0L[15]-0.25*(f0R[14]+f0L[14])+0.4330127018922193*(f0R[13]+f0R[12])-0.4330127018922193*(f0L[13]+f0L[12]+f0R[11])+0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10]+f0R[9]+f0L[9])+0.4330127018922193*f0R[8]-0.4330127018922193*f0L[8]+0.25*(f0R[7]+f0L[7])-0.4330127018922193*(f0R[6]+f0R[5])+0.4330127018922193*(f0L[6]+f0L[5])-0.25*(f0R[4]+f0L[4])+0.25*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.25*f0R[14]+0.25*f0L[14]-0.4330127018922193*(f0R[13]+f0L[13]+f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0L[11])+0.25*(f0R[10]+f0R[9])-0.25*(f0L[10]+f0L[9])+0.4330127018922193*(f0R[8]+f0L[8])-0.25*f0R[7]+0.25*f0L[7]-0.4330127018922193*(f0R[6]+f0L[6]+f0R[5]+f0L[5])-0.25*f0R[4]+0.25*(f0L[4]+f0R[3]+f0R[2])-0.25*(f0L[3]+f0L[2])+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*f0L[15]+0.25*(f0R[14]+f0L[14])+0.4330127018922193*(f0R[13]+f0R[12])-0.4330127018922193*(f0L[13]+f0L[12]+f0R[11])+0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10]+f0R[9]+f0L[9])-0.4330127018922193*f0R[8]+0.4330127018922193*f0L[8]+0.25*(f0R[7]+f0L[7])+0.4330127018922193*(f0R[6]+f0R[5])-0.4330127018922193*(f0L[6]+f0L[5])+0.25*(f0R[4]+f0L[4])-0.25*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.25*f0R[14]-0.25*f0L[14]-0.4330127018922193*(f0R[13]+f0L[13])+0.4330127018922193*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*(f0L[10]+f0R[9])+0.25*f0L[9]+0.4330127018922193*(f0R[8]+f0L[8])+0.25*f0R[7]-0.25*f0L[7]-0.4330127018922193*(f0R[6]+f0L[6])+0.4330127018922193*(f0R[5]+f0L[5])-0.25*f0R[4]+0.25*(f0L[4]+f0R[3])-0.25*(f0L[3]+f0R[2])+0.25*f0L[2]+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*f0L[15]-0.25*(f0R[14]+f0L[14])+0.4330127018922193*f0R[13]-0.4330127018922193*(f0L[13]+f0R[12])+0.4330127018922193*(f0L[12]+f0R[11])-0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])+0.25*(f0R[9]+f0L[9])-0.4330127018922193*f0R[8]+0.4330127018922193*f0L[8]-0.25*(f0R[7]+f0L[7])+0.4330127018922193*f0R[6]-0.4330127018922193*(f0L[6]+f0R[5])+0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4])-0.25*(f0R[3]+f0L[3])+0.25*(f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.25*f0R[14]-0.25*f0L[14]+0.4330127018922193*(f0R[13]+f0L[13])-0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*(f0L[10]+f0R[9])-0.25*f0L[9]+0.4330127018922193*(f0R[8]+f0L[8])+0.25*f0R[7]-0.25*f0L[7]+0.4330127018922193*(f0R[6]+f0L[6])-0.4330127018922193*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0R[3])+0.25*(f0L[4]+f0L[3]+f0R[2])-0.25*f0L[2]+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*f0L[15]-0.25*(f0R[14]+f0L[14])-0.4330127018922193*f0R[13]+0.4330127018922193*(f0L[13]+f0R[12]+f0R[11])-0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])-0.25*(f0R[9]+f0L[9])-0.4330127018922193*f0R[8]+0.4330127018922193*f0L[8]-0.25*(f0R[7]+f0L[7])-0.4330127018922193*f0R[6]+0.4330127018922193*(f0L[6]+f0R[5])-0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.25*(f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]; 
  fUpOrd[7] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.25*f0R[14]+0.25*f0L[14]+0.4330127018922193*(f0R[13]+f0L[13]+f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*(f0R[10]+f0R[9])+0.25*(f0L[10]+f0L[9])+0.4330127018922193*(f0R[8]+f0L[8])-0.25*f0R[7]+0.25*f0L[7]+0.4330127018922193*(f0R[6]+f0L[6]+f0R[5]+f0L[5])-0.25*(f0R[4]+f0R[3]+f0R[2])+0.25*(f0L[4]+f0L[3]+f0L[2])+0.4330127018922193*(f0R[1]+f0L[1])-0.25*f0R[0]+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*f0L[15]+0.25*(f0R[14]+f0L[14])-0.4330127018922193*(f0R[13]+f0R[12]+f0R[11])+0.4330127018922193*(f0L[13]+f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10]+f0R[9]+f0L[9])-0.4330127018922193*f0R[8]+0.4330127018922193*f0L[8]+0.25*(f0R[7]+f0L[7])-0.4330127018922193*(f0R[6]+f0R[5])+0.4330127018922193*(f0L[6]+f0L[5])+0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.4330127018922193*f0R[1]+0.4330127018922193*f0L[1]+0.25*(f0R[0]+f0L[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*alpha1R[0]*fUp[0]; 
  incr[1] = -0.4330127018922193*alpha1R[0]*fUp[0]; 
  incr[2] = 0.25*alpha1R[0]*fUp[1]; 
  incr[3] = 0.25*alpha1R[0]*fUp[2]; 
  incr[4] = 0.25*alpha1R[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alpha1R[0]*fUp[1]; 
  incr[6] = -0.4330127018922193*alpha1R[0]*fUp[2]; 
  incr[7] = 0.25*alpha1R[0]*fUp[4]; 
  incr[8] = -0.4330127018922193*alpha1R[0]*fUp[3]; 
  incr[9] = 0.25*alpha1R[0]*fUp[5]; 
  incr[10] = 0.25*alpha1R[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alpha1R[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alpha1R[0]*fUp[5]; 
  incr[13] = -0.4330127018922193*alpha1R[0]*fUp[6]; 
  incr[14] = 0.25*alpha1R[0]*fUp[7]; 
  incr[15] = -0.4330127018922193*alpha1R[0]*fUp[7]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 
  outR[8] += incr[8]*rdx2R; 
  outR[9] += incr[9]*rdx2R; 
  outR[10] += incr[10]*rdx2R; 
  outR[11] += incr[11]*rdx2R; 
  outR[12] += incr[12]*rdx2R; 
  outR[13] += incr[13]*rdx2R; 
  outR[14] += incr[14]*rdx2R; 
  outR[15] += incr[15]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += -1.0*incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += incr[6]*rdx2L; 
  outL[7] += -1.0*incr[7]*rdx2L; 
  outL[8] += incr[8]*rdx2L; 
  outL[9] += -1.0*incr[9]*rdx2L; 
  outL[10] += -1.0*incr[10]*rdx2L; 
  outL[11] += incr[11]*rdx2L; 
  outL[12] += incr[12]*rdx2L; 
  outL[13] += incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_y_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[8] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = -(1.732050807568877*jacobTotInv[0]*b_z[1]*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[1] = -(1.732050807568877*b_z[1]*jacobTotInv[1]*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alpha0R[8]; 
  alpha0R[0] = (0.3061862178478971*(hamil0R[1]*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*m_*rdx2R+BstarYdBmagR[0]*hamil0R[3]*q_*rdvpar2R))/(m_*q_); 
  alpha0R[1] = (0.3061862178478971*(hamil0R[1]*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*m_*rdx2R+BstarYdBmagR[1]*hamil0R[3]*q_*rdvpar2R))/(m_*q_); 
  alpha0R[2] = (0.3061862178478971*BstarYdBmagR[3]*hamil0R[3]*rdvpar2R)/m_; 
  alpha0R[3] = (0.3061862178478971*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamil0R[8]*rdx2R)/q_; 
  alpha0R[4] = (0.3061862178478971*hamil0R[3]*BstarYdBmagR[6]*rdvpar2R)/m_; 
  alpha0R[5] = (0.3061862178478971*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamil0R[8]*rdx2R)/q_; 

  double alpha1R[8]; 
  alpha1R[0] = -(0.125*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*(4.242640687119286*hamil1R[5]-2.449489742783178*hamil1R[1])*rdx2R)/q_; 
  alpha1R[1] = -(0.125*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*(4.242640687119286*hamil1R[5]-2.449489742783178*hamil1R[1])*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(((3.0*b_z[1]*jacobTotInv[1]+3.0*b_z[0]*jacobTotInv[0])*hamil1R[5]+((-1.732050807568877*b_z[1]*hamil1R[1])-1.732050807568877*b_z[1]*hamil0R[1])*jacobTotInv[1]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamil1R[1]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamil0R[1])*m_*rdx2R-1.732050807568877*BstarYdBmagR[0]*hamil0R[3]*q_*rdvpar2R))/(m_*q_); 

  double incr[16]; 
  // linear term alpha0*f1 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*(alpha0R[5]+alpha0R[4])-0.3535533905932737*(alpha0R[3]+alpha0R[2]+alpha0R[1])+0.3535533905932737*alpha0R[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(f1R[15]+f1L[15]))+0.4330127018922193*(f1R[14]+f1L[14])+0.25*f1R[13]-0.25*f1L[13]+0.4330127018922193*(f1R[12]+f1L[12]+f1R[11]+f1L[11])-0.25*f1R[10]+0.25*f1L[10]-0.4330127018922193*(f1R[9]+f1L[9])-0.25*f1R[8]+0.25*f1L[8]-0.4330127018922193*(f1R[7]+f1L[7])-0.25*f1R[6]+0.25*f1L[6]-0.4330127018922193*(f1R[5]+f1L[5])+0.25*(f1R[4]+f1R[3])-0.25*(f1L[4]+f1L[3])+0.4330127018922193*(f1R[2]+f1L[2])+0.25*f1R[1]-0.25*(f1L[1]+f1R[0])+0.25*f1L[0])*sgn(alphaOrdR)+0.4330127018922193*f1R[15]-0.4330127018922193*(f1L[15]+f1R[14])+0.4330127018922193*f1L[14]-0.25*(f1R[13]+f1L[13])-0.4330127018922193*(f1R[12]+f1R[11])+0.4330127018922193*(f1L[12]+f1L[11])+0.25*(f1R[10]+f1L[10])+0.4330127018922193*f1R[9]-0.4330127018922193*f1L[9]+0.25*(f1R[8]+f1L[8])+0.4330127018922193*f1R[7]-0.4330127018922193*f1L[7]+0.25*(f1R[6]+f1L[6])+0.4330127018922193*f1R[5]-0.4330127018922193*f1L[5]-0.25*(f1R[4]+f1L[4]+f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]-0.25*(f1R[1]+f1L[1])+0.25*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha0R[1]+alpha0R[0])-0.3535533905932737*(alpha0R[5]+alpha0R[4]+alpha0R[3]+alpha0R[2]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(f1R[15]+f1L[15]+f1R[14]+f1L[14])-0.25*f1R[13]+0.25*f1L[13]-0.4330127018922193*(f1R[12]+f1L[12]+f1R[11]+f1L[11])-0.25*f1R[10]+0.25*f1L[10]-0.4330127018922193*(f1R[9]+f1L[9])+0.25*f1R[8]-0.25*f1L[8]-0.4330127018922193*(f1R[7]+f1L[7])+0.25*f1R[6]-0.25*f1L[6]+0.4330127018922193*(f1R[5]+f1L[5])+0.25*(f1R[4]+f1R[3])-0.25*(f1L[4]+f1L[3])+0.4330127018922193*(f1R[2]+f1L[2])-0.25*(f1R[1]+f1R[0])+0.25*(f1L[1]+f1L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f1R[15]+f1R[14])+0.4330127018922193*(f1L[15]+f1L[14])+0.25*(f1R[13]+f1L[13])+0.4330127018922193*(f1R[12]+f1R[11])-0.4330127018922193*(f1L[12]+f1L[11])+0.25*(f1R[10]+f1L[10])+0.4330127018922193*f1R[9]-0.4330127018922193*f1L[9]-0.25*(f1R[8]+f1L[8])+0.4330127018922193*f1R[7]-0.4330127018922193*f1L[7]-0.25*(f1R[6]+f1L[6])-0.4330127018922193*f1R[5]+0.4330127018922193*f1L[5]-0.25*(f1R[4]+f1L[4]+f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]+0.25*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 
  alphaOrdR = 0.3535533905932737*alpha0R[5]-0.3535533905932737*(alpha0R[4]+alpha0R[3])+0.3535533905932737*alpha0R[2]-0.3535533905932737*alpha0R[1]+0.3535533905932737*alpha0R[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(f1R[15]+f1L[15])-0.4330127018922193*(f1R[14]+f1L[14])-0.25*f1R[13]+0.25*f1L[13]+0.4330127018922193*(f1R[12]+f1L[12])-0.4330127018922193*(f1R[11]+f1L[11])+0.25*f1R[10]-0.25*f1L[10]-0.4330127018922193*(f1R[9]+f1L[9])-0.25*f1R[8]+0.25*f1L[8]+0.4330127018922193*(f1R[7]+f1L[7])+0.25*f1R[6]-0.25*f1L[6]-0.4330127018922193*(f1R[5]+f1L[5])+0.25*f1R[4]-0.25*(f1L[4]+f1R[3])+0.25*f1L[3]+0.4330127018922193*(f1R[2]+f1L[2])+0.25*f1R[1]-0.25*(f1L[1]+f1R[0])+0.25*f1L[0])*sgn(alphaOrdR)-0.4330127018922193*f1R[15]+0.4330127018922193*(f1L[15]+f1R[14])-0.4330127018922193*f1L[14]+0.25*(f1R[13]+f1L[13])-0.4330127018922193*f1R[12]+0.4330127018922193*(f1L[12]+f1R[11])-0.4330127018922193*f1L[11]-0.25*(f1R[10]+f1L[10])+0.4330127018922193*f1R[9]-0.4330127018922193*f1L[9]+0.25*(f1R[8]+f1L[8])-0.4330127018922193*f1R[7]+0.4330127018922193*f1L[7]-0.25*(f1R[6]+f1L[6])+0.4330127018922193*f1R[5]-0.4330127018922193*f1L[5]-0.25*(f1R[4]+f1L[4])+0.25*(f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]-0.25*(f1R[1]+f1L[1])+0.25*(f1R[0]+f1L[0])); 
  alphaOrdR = (-0.3535533905932737*alpha0R[5])+0.3535533905932737*alpha0R[4]-0.3535533905932737*alpha0R[3]+0.3535533905932737*(alpha0R[2]+alpha0R[1]+alpha0R[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(f1R[15]+f1L[15]+f1R[14]+f1L[14]))+0.25*f1R[13]-0.25*f1L[13]-0.4330127018922193*(f1R[12]+f1L[12])+0.4330127018922193*(f1R[11]+f1L[11])+0.25*f1R[10]-0.25*f1L[10]-0.4330127018922193*(f1R[9]+f1L[9])+0.25*f1R[8]-0.25*f1L[8]+0.4330127018922193*(f1R[7]+f1L[7])-0.25*f1R[6]+0.25*f1L[6]+0.4330127018922193*(f1R[5]+f1L[5])+0.25*f1R[4]-0.25*(f1L[4]+f1R[3])+0.25*f1L[3]+0.4330127018922193*(f1R[2]+f1L[2])-0.25*(f1R[1]+f1R[0])+0.25*(f1L[1]+f1L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f1R[15]+f1R[14])-0.4330127018922193*(f1L[15]+f1L[14])-0.25*(f1R[13]+f1L[13])+0.4330127018922193*f1R[12]-0.4330127018922193*(f1L[12]+f1R[11])+0.4330127018922193*f1L[11]-0.25*(f1R[10]+f1L[10])+0.4330127018922193*f1R[9]-0.4330127018922193*f1L[9]-0.25*(f1R[8]+f1L[8])-0.4330127018922193*f1R[7]+0.4330127018922193*f1L[7]+0.25*(f1R[6]+f1L[6])-0.4330127018922193*f1R[5]+0.4330127018922193*f1L[5]-0.25*(f1R[4]+f1L[4])+0.25*(f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]+0.25*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 
  alphaOrdR = (-0.3535533905932737*alpha0R[5])+0.3535533905932737*(alpha0R[4]+alpha0R[3])-0.3535533905932737*(alpha0R[2]+alpha0R[1])+0.3535533905932737*alpha0R[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(f1R[15]+f1L[15])-0.4330127018922193*(f1R[14]+f1L[14])-0.25*f1R[13]+0.25*f1L[13]-0.4330127018922193*(f1R[12]+f1L[12])+0.4330127018922193*(f1R[11]+f1L[11])+0.25*f1R[10]-0.25*f1L[10]+0.4330127018922193*(f1R[9]+f1L[9])+0.25*f1R[8]-0.25*f1L[8]-0.4330127018922193*(f1R[7]+f1L[7])-0.25*f1R[6]+0.25*f1L[6]-0.4330127018922193*(f1R[5]+f1L[5])-0.25*f1R[4]+0.25*(f1L[4]+f1R[3])-0.25*f1L[3]+0.4330127018922193*(f1R[2]+f1L[2])+0.25*f1R[1]-0.25*(f1L[1]+f1R[0])+0.25*f1L[0])*sgn(alphaOrdR)-0.4330127018922193*f1R[15]+0.4330127018922193*(f1L[15]+f1R[14])-0.4330127018922193*f1L[14]+0.25*(f1R[13]+f1L[13])+0.4330127018922193*f1R[12]-0.4330127018922193*(f1L[12]+f1R[11])+0.4330127018922193*f1L[11]-0.25*(f1R[10]+f1L[10])-0.4330127018922193*f1R[9]+0.4330127018922193*f1L[9]-0.25*(f1R[8]+f1L[8])+0.4330127018922193*f1R[7]-0.4330127018922193*f1L[7]+0.25*(f1R[6]+f1L[6])+0.4330127018922193*f1R[5]-0.4330127018922193*f1L[5]+0.25*(f1R[4]+f1L[4])-0.25*(f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]-0.25*(f1R[1]+f1L[1])+0.25*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.3535533905932737*alpha0R[5]-0.3535533905932737*alpha0R[4]+0.3535533905932737*alpha0R[3]-0.3535533905932737*alpha0R[2]+0.3535533905932737*(alpha0R[1]+alpha0R[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(f1R[15]+f1L[15]+f1R[14]+f1L[14]))+0.25*f1R[13]-0.25*f1L[13]+0.4330127018922193*(f1R[12]+f1L[12])-0.4330127018922193*(f1R[11]+f1L[11])+0.25*f1R[10]-0.25*f1L[10]+0.4330127018922193*(f1R[9]+f1L[9])-0.25*f1R[8]+0.25*f1L[8]-0.4330127018922193*(f1R[7]+f1L[7])+0.25*f1R[6]-0.25*f1L[6]+0.4330127018922193*(f1R[5]+f1L[5])-0.25*f1R[4]+0.25*(f1L[4]+f1R[3])-0.25*f1L[3]+0.4330127018922193*(f1R[2]+f1L[2])-0.25*(f1R[1]+f1R[0])+0.25*(f1L[1]+f1L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f1R[15]+f1R[14])-0.4330127018922193*(f1L[15]+f1L[14])-0.25*(f1R[13]+f1L[13])-0.4330127018922193*f1R[12]+0.4330127018922193*(f1L[12]+f1R[11])-0.4330127018922193*f1L[11]-0.25*(f1R[10]+f1L[10])-0.4330127018922193*f1R[9]+0.4330127018922193*f1L[9]+0.25*(f1R[8]+f1L[8])+0.4330127018922193*f1R[7]-0.4330127018922193*f1L[7]-0.25*(f1R[6]+f1L[6])-0.4330127018922193*f1R[5]+0.4330127018922193*f1L[5]+0.25*(f1R[4]+f1L[4])-0.25*(f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]+0.25*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 
  alphaOrdR = (-0.3535533905932737*(alpha0R[5]+alpha0R[4]))+0.3535533905932737*(alpha0R[3]+alpha0R[2])-0.3535533905932737*alpha0R[1]+0.3535533905932737*alpha0R[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(f1R[15]+f1L[15]))+0.4330127018922193*(f1R[14]+f1L[14])+0.25*f1R[13]-0.25*f1L[13]-0.4330127018922193*(f1R[12]+f1L[12]+f1R[11]+f1L[11])-0.25*f1R[10]+0.25*f1L[10]+0.4330127018922193*(f1R[9]+f1L[9])+0.25*f1R[8]-0.25*f1L[8]+0.4330127018922193*(f1R[7]+f1L[7])+0.25*f1R[6]-0.25*f1L[6]-0.4330127018922193*(f1R[5]+f1L[5])-0.25*(f1R[4]+f1R[3])+0.25*(f1L[4]+f1L[3])+0.4330127018922193*(f1R[2]+f1L[2])+0.25*f1R[1]-0.25*(f1L[1]+f1R[0])+0.25*f1L[0])*sgn(alphaOrdR)+0.4330127018922193*f1R[15]-0.4330127018922193*(f1L[15]+f1R[14])+0.4330127018922193*f1L[14]-0.25*(f1R[13]+f1L[13])+0.4330127018922193*(f1R[12]+f1R[11])-0.4330127018922193*(f1L[12]+f1L[11])+0.25*(f1R[10]+f1L[10])-0.4330127018922193*f1R[9]+0.4330127018922193*f1L[9]-0.25*(f1R[8]+f1L[8])-0.4330127018922193*f1R[7]+0.4330127018922193*f1L[7]-0.25*(f1R[6]+f1L[6])+0.4330127018922193*f1R[5]-0.4330127018922193*f1L[5]+0.25*(f1R[4]+f1L[4]+f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]-0.25*(f1R[1]+f1L[1])+0.25*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha0R[5]+alpha0R[4]+alpha0R[3]+alpha0R[2]+alpha0R[1]+alpha0R[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(f1R[15]+f1L[15]+f1R[14]+f1L[14])-0.25*f1R[13]+0.25*f1L[13]+0.4330127018922193*(f1R[12]+f1L[12]+f1R[11]+f1L[11])-0.25*f1R[10]+0.25*f1L[10]+0.4330127018922193*(f1R[9]+f1L[9])-0.25*f1R[8]+0.25*f1L[8]+0.4330127018922193*(f1R[7]+f1L[7])-0.25*f1R[6]+0.25*f1L[6]+0.4330127018922193*(f1R[5]+f1L[5])-0.25*(f1R[4]+f1R[3])+0.25*(f1L[4]+f1L[3])+0.4330127018922193*(f1R[2]+f1L[2])-0.25*(f1R[1]+f1R[0])+0.25*(f1L[1]+f1L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f1R[15]+f1R[14])+0.4330127018922193*(f1L[15]+f1L[14])+0.25*(f1R[13]+f1L[13])-0.4330127018922193*(f1R[12]+f1R[11])+0.4330127018922193*(f1L[12]+f1L[11])+0.25*(f1R[10]+f1L[10])-0.4330127018922193*f1R[9]+0.4330127018922193*f1L[9]+0.25*(f1R[8]+f1L[8])-0.4330127018922193*f1R[7]+0.4330127018922193*f1L[7]+0.25*(f1R[6]+f1L[6])-0.4330127018922193*f1R[5]+0.4330127018922193*f1L[5]+0.25*(f1R[4]+f1L[4]+f1R[3]+f1L[3])-0.4330127018922193*f1R[2]+0.4330127018922193*f1L[2]+0.25*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alpha0R[5]*fUp[5]+alpha0R[4]*fUp[4]+alpha0R[3]*fUp[3]+alpha0R[2]*fUp[2]+alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[1] = 0.25*(alpha0R[3]*fUp[5]+fUp[3]*alpha0R[5]+alpha0R[2]*fUp[4]+fUp[2]*alpha0R[4]+alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1]); 
  incr[2] = -0.4330127018922193*(alpha0R[5]*fUp[5]+alpha0R[4]*fUp[4]+alpha0R[3]*fUp[3]+alpha0R[2]*fUp[2]+alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[3] = 0.25*(alpha0R[5]*fUp[7]+alpha0R[3]*fUp[6]+alpha0R[1]*fUp[4]+fUp[1]*alpha0R[4]+alpha0R[0]*fUp[2]+fUp[0]*alpha0R[2]); 
  incr[4] = 0.25*(alpha0R[4]*fUp[7]+alpha0R[2]*fUp[6]+alpha0R[1]*fUp[5]+fUp[1]*alpha0R[5]+alpha0R[0]*fUp[3]+fUp[0]*alpha0R[3]); 
  incr[5] = -0.4330127018922193*(alpha0R[3]*fUp[5]+fUp[3]*alpha0R[5]+alpha0R[2]*fUp[4]+fUp[2]*alpha0R[4]+alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1]); 
  incr[6] = 0.25*(alpha0R[3]*fUp[7]+alpha0R[5]*fUp[6]+alpha0R[0]*fUp[4]+fUp[0]*alpha0R[4]+alpha0R[1]*fUp[2]+fUp[1]*alpha0R[2]); 
  incr[7] = -0.4330127018922193*(alpha0R[5]*fUp[7]+alpha0R[3]*fUp[6]+alpha0R[1]*fUp[4]+fUp[1]*alpha0R[4]+alpha0R[0]*fUp[2]+fUp[0]*alpha0R[2]); 
  incr[8] = 0.25*(alpha0R[2]*fUp[7]+alpha0R[4]*fUp[6]+alpha0R[0]*fUp[5]+fUp[0]*alpha0R[5]+alpha0R[1]*fUp[3]+fUp[1]*alpha0R[3]); 
  incr[9] = -0.4330127018922193*(alpha0R[4]*fUp[7]+alpha0R[2]*fUp[6]+alpha0R[1]*fUp[5]+fUp[1]*alpha0R[5]+alpha0R[0]*fUp[3]+fUp[0]*alpha0R[3]); 
  incr[10] = 0.25*(alpha0R[1]*fUp[7]+alpha0R[0]*fUp[6]+alpha0R[4]*fUp[5]+fUp[4]*alpha0R[5]+alpha0R[2]*fUp[3]+fUp[2]*alpha0R[3]); 
  incr[11] = -0.4330127018922193*(alpha0R[3]*fUp[7]+alpha0R[5]*fUp[6]+alpha0R[0]*fUp[4]+fUp[0]*alpha0R[4]+alpha0R[1]*fUp[2]+fUp[1]*alpha0R[2]); 
  incr[12] = -0.4330127018922193*(alpha0R[2]*fUp[7]+alpha0R[4]*fUp[6]+alpha0R[0]*fUp[5]+fUp[0]*alpha0R[5]+alpha0R[1]*fUp[3]+fUp[1]*alpha0R[3]); 
  incr[13] = 0.25*(alpha0R[0]*fUp[7]+alpha0R[1]*fUp[6]+alpha0R[2]*fUp[5]+fUp[2]*alpha0R[5]+alpha0R[3]*fUp[4]+fUp[3]*alpha0R[4]); 
  incr[14] = -0.4330127018922193*(alpha0R[1]*fUp[7]+alpha0R[0]*fUp[6]+alpha0R[4]*fUp[5]+fUp[4]*alpha0R[5]+alpha0R[2]*fUp[3]+fUp[2]*alpha0R[3]); 
  incr[15] = -0.4330127018922193*(alpha0R[0]*fUp[7]+alpha0R[1]*fUp[6]+alpha0R[2]*fUp[5]+fUp[2]*alpha0R[5]+alpha0R[3]*fUp[4]+fUp[3]*alpha0R[4]); 

  outR[0] += incr[0]*rdy2R; 
  outR[1] += incr[1]*rdy2R; 
  outR[2] += incr[2]*rdy2R; 
  outR[3] += incr[3]*rdy2R; 
  outR[4] += incr[4]*rdy2R; 
  outR[5] += incr[5]*rdy2R; 
  outR[6] += incr[6]*rdy2R; 
  outR[7] += incr[7]*rdy2R; 
  outR[8] += incr[8]*rdy2R; 
  outR[9] += incr[9]*rdy2R; 
  outR[10] += incr[10]*rdy2R; 
  outR[11] += incr[11]*rdy2R; 
  outR[12] += incr[12]*rdy2R; 
  outR[13] += incr[13]*rdy2R; 
  outR[14] += incr[14]*rdy2R; 
  outR[15] += incr[15]*rdy2R; 

  outL[0] += -1.0*incr[0]*rdy2L; 
  outL[1] += -1.0*incr[1]*rdy2L; 
  outL[2] += incr[2]*rdy2L; 
  outL[3] += -1.0*incr[3]*rdy2L; 
  outL[4] += -1.0*incr[4]*rdy2L; 
  outL[5] += incr[5]*rdy2L; 
  outL[6] += -1.0*incr[6]*rdy2L; 
  outL[7] += incr[7]*rdy2L; 
  outL[8] += -1.0*incr[8]*rdy2L; 
  outL[9] += incr[9]*rdy2L; 
  outL[10] += -1.0*incr[10]*rdy2L; 
  outL[11] += incr[11]*rdy2L; 
  outL[12] += incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.4330127018922193*(f0R[14]+f0L[14])+0.25*f0R[13]-0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0R[3])-0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*(f0L[15]+f0R[14])+0.4330127018922193*f0L[14]-0.25*(f0R[13]+f0L[13])-0.4330127018922193*(f0R[12]+f0R[11])+0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0R[3])-0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f0R[15]+f0R[14])+0.4330127018922193*(f0L[15]+f0L[14])+0.25*(f0R[13]+f0L[13])+0.4330127018922193*(f0R[12]+f0R[11])-0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.4330127018922193*(f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])+0.25*f0R[4]-0.25*(f0L[4]+f0R[3])+0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*(f0L[15]+f0R[14])-0.4330127018922193*f0L[14]+0.25*(f0R[13]+f0L[13])-0.4330127018922193*f0R[12]+0.4330127018922193*(f0L[12]+f0R[11])-0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4])+0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]))+0.25*f0R[13]-0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]-0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])+0.25*f0R[4]-0.25*(f0L[4]+f0R[3])+0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f0R[15]+f0R[14])-0.4330127018922193*(f0L[15]+f0L[14])-0.25*(f0R[13]+f0L[13])+0.4330127018922193*f0R[12]-0.4330127018922193*(f0L[12]+f0R[11])+0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])+0.4330127018922193*f0R[9]-0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]-0.25*(f0R[4]+f0L[4])+0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.4330127018922193*(f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])-0.25*f0R[4]+0.25*(f0L[4]+f0R[3])-0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*(f0L[15]+f0R[14])-0.4330127018922193*f0L[14]+0.25*(f0R[13]+f0L[13])+0.4330127018922193*f0R[12]-0.4330127018922193*(f0L[12]+f0R[11])+0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4])-0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]))+0.25*f0R[13]-0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0L[11])+0.25*f0R[10]-0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])-0.25*f0R[4]+0.25*(f0L[4]+f0R[3])-0.25*f0L[3]+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f0R[15]+f0R[14])-0.4330127018922193*(f0L[15]+f0L[14])-0.25*(f0R[13]+f0L[13])-0.4330127018922193*f0R[12]+0.4330127018922193*(f0L[12]+f0R[11])-0.4330127018922193*f0L[11]-0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4])-0.25*(f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.4330127018922193*(f0R[14]+f0L[14])+0.25*f0R[13]-0.25*f0L[13]-0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])+0.25*f0R[8]-0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])+0.25*f0R[6]-0.25*f0L[6]-0.4330127018922193*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0R[3])+0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])+0.25*f0R[1]-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*(f0L[15]+f0R[14])+0.4330127018922193*f0L[14]-0.25*(f0R[13]+f0L[13])+0.4330127018922193*(f0R[12]+f0R[11])-0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]-0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]-0.25*(f0R[6]+f0L[6])+0.4330127018922193*f0R[5]-0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14])-0.25*f0R[13]+0.25*f0L[13]+0.4330127018922193*(f0R[12]+f0L[12]+f0R[11]+f0L[11])-0.25*f0R[10]+0.25*f0L[10]+0.4330127018922193*(f0R[9]+f0L[9])-0.25*f0R[8]+0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])-0.25*f0R[6]+0.25*f0L[6]+0.4330127018922193*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0R[3])+0.25*(f0L[4]+f0L[3])+0.4330127018922193*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f0R[15]+f0R[14])+0.4330127018922193*(f0L[15]+f0L[14])+0.25*(f0R[13]+f0L[13])-0.4330127018922193*(f0R[12]+f0R[11])+0.4330127018922193*(f0L[12]+f0L[11])+0.25*(f0R[10]+f0L[10])-0.4330127018922193*f0R[9]+0.4330127018922193*f0L[9]+0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*f0L[7]+0.25*(f0R[6]+f0L[6])-0.4330127018922193*f0R[5]+0.4330127018922193*f0L[5]+0.25*(f0R[4]+f0L[4]+f0R[3]+f0L[3])-0.4330127018922193*f0R[2]+0.4330127018922193*f0L[2]+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.25*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.4330127018922193*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = 0.25*(alpha1R[1]*fUp[4]+alpha1R[0]*fUp[2]); 
  incr[4] = 0.25*(alpha1R[1]*fUp[5]+alpha1R[0]*fUp[3]); 
  incr[5] = -0.4330127018922193*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[6] = 0.25*(alpha1R[0]*fUp[4]+alpha1R[1]*fUp[2]); 
  incr[7] = -0.4330127018922193*(alpha1R[1]*fUp[4]+alpha1R[0]*fUp[2]); 
  incr[8] = 0.25*(alpha1R[0]*fUp[5]+alpha1R[1]*fUp[3]); 
  incr[9] = -0.4330127018922193*(alpha1R[1]*fUp[5]+alpha1R[0]*fUp[3]); 
  incr[10] = 0.25*(alpha1R[1]*fUp[7]+alpha1R[0]*fUp[6]); 
  incr[11] = -0.4330127018922193*(alpha1R[0]*fUp[4]+alpha1R[1]*fUp[2]); 
  incr[12] = -0.4330127018922193*(alpha1R[0]*fUp[5]+alpha1R[1]*fUp[3]); 
  incr[13] = 0.25*(alpha1R[0]*fUp[7]+alpha1R[1]*fUp[6]); 
  incr[14] = -0.4330127018922193*(alpha1R[1]*fUp[7]+alpha1R[0]*fUp[6]); 
  incr[15] = -0.4330127018922193*(alpha1R[0]*fUp[7]+alpha1R[1]*fUp[6]); 

  outR[0] += incr[0]*rdy2R; 
  outR[1] += incr[1]*rdy2R; 
  outR[2] += incr[2]*rdy2R; 
  outR[3] += incr[3]*rdy2R; 
  outR[4] += incr[4]*rdy2R; 
  outR[5] += incr[5]*rdy2R; 
  outR[6] += incr[6]*rdy2R; 
  outR[7] += incr[7]*rdy2R; 
  outR[8] += incr[8]*rdy2R; 
  outR[9] += incr[9]*rdy2R; 
  outR[10] += incr[10]*rdy2R; 
  outR[11] += incr[11]*rdy2R; 
  outR[12] += incr[12]*rdy2R; 
  outR[13] += incr[13]*rdy2R; 
  outR[14] += incr[14]*rdy2R; 
  outR[15] += incr[15]*rdy2R; 

  outL[0] += -1.0*incr[0]*rdy2L; 
  outL[1] += -1.0*incr[1]*rdy2L; 
  outL[2] += incr[2]*rdy2L; 
  outL[3] += -1.0*incr[3]*rdy2L; 
  outL[4] += -1.0*incr[4]*rdy2L; 
  outL[5] += incr[5]*rdy2L; 
  outL[6] += -1.0*incr[6]*rdy2L; 
  outL[7] += incr[7]*rdy2L; 
  outL[8] += -1.0*incr[8]*rdy2L; 
  outL[9] += incr[9]*rdy2L; 
  outL[10] += -1.0*incr[10]*rdy2L; 
  outL[11] += incr[11]*rdy2L; 
  outL[12] += incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 

  return std::abs(alphaSurfAvgR); 
} 
double LinearDeltaFGyrokineticGenGeoSurf2x2vSer_vpar_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR) 
{ 
  // jacobTotInv: reciprocal of the conf-space jacobian time the guiding center coordinate Jacobian.
  // b_x,b_y,b_z: covariant components of the field aligned unit vector.
  // q_,m_: species charge and mass.
  // cflL,cflR: CFL rate in left and right cells.
  // wL[NDIM],wR[NDIM]: cell-center in left and right cells.
  // dxvL[NDIM],dxvR[NDIM]: Cell length in left and right cells.
  // amax_in: maximum phase-space speed.
  // bmag: magnetic field amplitude.
  // cmag: coefficient multiplying parallel gradient.
  // phi: electrostatic potential .
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wyL = wL[1];
  double wyR = wR[1];
  double rdy2L = 2.0/dxvL[1];
  double rdy2R = 2.0/dxvR[1];
  double wvparL = wL[2];
  double wvparR = wR[2];
  double rdvpar2L = 2.0/dxvL[2];
  double rdvpar2R = 2.0/dxvR[2];
  double wmuL = wL[3];
  double wmuR = wR[3];
  double rdmu2L = 2.0/dxvL[3];
  double rdmu2R = 2.0/dxvR[3];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wySqL = wL[1]*wL[1];
  double wySqR = wR[1]*wR[1];
  double rdy2SqL = rdy2L*rdy2L;
  double rdy2SqR = rdy2R*rdy2R;
  double wvparSqL = wL[2]*wL[2];
  double wvparSqR = wR[2]*wR[2];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[3]*wL[3];
  double wmuSqR = wR[3]*wR[3];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[16]; 
  hamil0R[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR)+m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamil0R[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[8] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[16]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 
  hamil1R[2] = 2.0*phi[2]*q_; 
  hamil1R[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagR[16]; 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = -(1.732050807568877*jacobTotInv[0]*b_z[1]*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[1] = -(1.732050807568877*b_z[1]*jacobTotInv[1]*m_*rdx2R*wvparR)/q_; 
  BstarYdBmagR[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alpha0R[8]; 

  double alpha1R[8]; 
  alpha1R[0] = (0.125*(hamil1R[5]*(4.242640687119286*BstarYdBmagR[6]-2.449489742783178*BstarYdBmagR[1])+hamil1R[2]*(4.242640687119286*BstarYdBmagR[3]-2.449489742783178*BstarYdBmagR[0]))*rdy2R)/m_; 
  alpha1R[1] = (0.125*(4.242640687119286*hamil1R[2]*BstarYdBmagR[6]+(4.242640687119286*BstarYdBmagR[3]-2.449489742783178*BstarYdBmagR[0])*hamil1R[5]-2.449489742783178*BstarYdBmagR[1]*hamil1R[2])*rdy2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*(3.0*hamil1R[5]*BstarYdBmagR[6]-1.732050807568877*BstarYdBmagR[1]*hamil1R[5]+3.0*hamil1R[2]*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0]*hamil1R[2])*rdy2R)/m_; 

  double incr[16]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.4330127018922193*(f0R[14]+f0L[14]+f0R[13]+f0L[13])+0.25*f0R[12]-0.25*f0L[12]+0.4330127018922193*(f0R[11]+f0L[11])-0.4330127018922193*(f0R[10]+f0L[10])-0.25*(f0R[9]+f0R[8])+0.25*(f0L[9]+f0L[8])-0.4330127018922193*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.25*f0R[5]+0.25*(f0L[5]+f0R[4])-0.25*f0L[4]+0.4330127018922193*(f0R[3]+f0L[3])+0.25*(f0R[2]+f0R[1])-0.25*(f0L[2]+f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*(f0L[15]+f0R[14]+f0R[13])+0.4330127018922193*(f0L[14]+f0L[13])-0.25*(f0R[12]+f0L[12])-0.4330127018922193*f0R[11]+0.4330127018922193*(f0L[11]+f0R[10])-0.4330127018922193*f0L[10]+0.25*(f0R[9]+f0L[9]+f0R[8]+f0L[8])+0.4330127018922193*(f0R[7]+f0R[6])-0.4330127018922193*(f0L[7]+f0L[6])+0.25*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]-0.25*(f0R[2]+f0L[2]+f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14])-0.4330127018922193*(f0R[13]+f0L[13])-0.25*f0R[12]+0.25*f0L[12]-0.4330127018922193*(f0R[11]+f0L[11]+f0R[10]+f0L[10])-0.25*f0R[9]+0.25*(f0L[9]+f0R[8])-0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])+0.4330127018922193*(f0R[6]+f0L[6])+0.25*(f0R[5]+f0R[4])-0.25*(f0L[5]+f0L[4])+0.4330127018922193*(f0R[3]+f0L[3])+0.25*f0R[2]-0.25*(f0L[2]+f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f0R[15]+f0R[14])+0.4330127018922193*(f0L[15]+f0L[14]+f0R[13])-0.4330127018922193*f0L[13]+0.25*(f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0R[10])-0.4330127018922193*(f0L[11]+f0L[10])+0.25*(f0R[9]+f0L[9])-0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*(f0L[7]+f0R[6])+0.4330127018922193*f0L[6]-0.25*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]-0.25*(f0R[2]+f0L[2])+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.4330127018922193*(f0R[14]+f0L[14])+0.4330127018922193*(f0R[13]+f0L[13])-0.25*f0R[12]+0.25*f0L[12]-0.4330127018922193*(f0R[11]+f0L[11]+f0R[10]+f0L[10])+0.25*f0R[9]-0.25*(f0L[9]+f0R[8])+0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])-0.4330127018922193*(f0R[6]+f0L[6])+0.25*(f0R[5]+f0R[4])-0.25*(f0L[5]+f0L[4])+0.4330127018922193*(f0R[3]+f0L[3])-0.25*f0R[2]+0.25*(f0L[2]+f0R[1])-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*(f0L[15]+f0R[14])-0.4330127018922193*(f0L[14]+f0R[13])+0.4330127018922193*f0L[13]+0.25*(f0R[12]+f0L[12])+0.4330127018922193*(f0R[11]+f0R[10])-0.4330127018922193*(f0L[11]+f0L[10])-0.25*(f0R[9]+f0L[9])+0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*(f0L[7]+f0R[6])-0.4330127018922193*f0L[6]-0.25*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]+0.25*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]+f0R[13]+f0L[13]))+0.25*f0R[12]-0.25*f0L[12]+0.4330127018922193*(f0R[11]+f0L[11])-0.4330127018922193*(f0R[10]+f0L[10])+0.25*(f0R[9]+f0R[8])-0.25*(f0L[9]+f0L[8])+0.4330127018922193*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.25*f0R[5]+0.25*(f0L[5]+f0R[4])-0.25*f0L[4]+0.4330127018922193*(f0R[3]+f0L[3])-0.25*(f0R[2]+f0R[1]+f0R[0])+0.25*(f0L[2]+f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f0R[15]+f0R[14]+f0R[13])-0.4330127018922193*(f0L[15]+f0L[14]+f0L[13])-0.25*(f0R[12]+f0L[12])-0.4330127018922193*f0R[11]+0.4330127018922193*(f0L[11]+f0R[10])-0.4330127018922193*f0L[10]-0.25*(f0R[9]+f0L[9]+f0R[8]+f0L[8])-0.4330127018922193*(f0R[7]+f0R[6])+0.4330127018922193*(f0L[7]+f0L[6])+0.25*(f0R[5]+f0L[5])-0.25*(f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]+0.25*(f0R[2]+f0L[2]+f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15])-0.4330127018922193*(f0R[14]+f0L[14]+f0R[13]+f0L[13])-0.25*f0R[12]+0.25*f0L[12]+0.4330127018922193*(f0R[11]+f0L[11]+f0R[10]+f0L[10])+0.25*(f0R[9]+f0R[8])-0.25*(f0L[9]+f0L[8])-0.4330127018922193*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.25*(f0R[5]+f0R[4])+0.25*(f0L[5]+f0L[4])+0.4330127018922193*(f0R[3]+f0L[3])+0.25*(f0R[2]+f0R[1])-0.25*(f0L[2]+f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)-0.4330127018922193*f0R[15]+0.4330127018922193*(f0L[15]+f0R[14]+f0R[13])-0.4330127018922193*(f0L[14]+f0L[13])+0.25*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0R[10])+0.4330127018922193*(f0L[11]+f0L[10])-0.25*(f0R[9]+f0L[9]+f0R[8]+f0L[8])+0.4330127018922193*(f0R[7]+f0R[6])-0.4330127018922193*(f0L[7]+f0L[6])+0.25*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]-0.25*(f0R[2]+f0L[2]+f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]))+0.4330127018922193*(f0R[13]+f0L[13])+0.25*f0R[12]-0.25*f0L[12]-0.4330127018922193*(f0R[11]+f0L[11])+0.4330127018922193*(f0R[10]+f0L[10])+0.25*f0R[9]-0.25*(f0L[9]+f0R[8])+0.25*f0L[8]-0.4330127018922193*(f0R[7]+f0L[7])+0.4330127018922193*(f0R[6]+f0L[6])+0.25*f0R[5]-0.25*(f0L[5]+f0R[4])+0.25*f0L[4]+0.4330127018922193*(f0R[3]+f0L[3])+0.25*f0R[2]-0.25*(f0L[2]+f0R[1]+f0R[0])+0.25*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.4330127018922193*(f0R[15]+f0R[14])-0.4330127018922193*(f0L[15]+f0L[14]+f0R[13])+0.4330127018922193*f0L[13]-0.25*(f0R[12]+f0L[12])+0.4330127018922193*f0R[11]-0.4330127018922193*(f0L[11]+f0R[10])+0.4330127018922193*f0L[10]-0.25*(f0R[9]+f0L[9])+0.25*(f0R[8]+f0L[8])+0.4330127018922193*f0R[7]-0.4330127018922193*(f0L[7]+f0R[6])+0.4330127018922193*f0L[6]-0.25*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]-0.25*(f0R[2]+f0L[2])+0.25*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*alpha1R[0]-0.3535533905932737*alpha1R[1]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(f0R[15]+f0L[15]))+0.4330127018922193*(f0R[14]+f0L[14])-0.4330127018922193*(f0R[13]+f0L[13])+0.25*f0R[12]-0.25*f0L[12]-0.4330127018922193*(f0R[11]+f0L[11])+0.4330127018922193*(f0R[10]+f0L[10])-0.25*f0R[9]+0.25*(f0L[9]+f0R[8])-0.25*f0L[8]+0.4330127018922193*(f0R[7]+f0L[7])-0.4330127018922193*(f0R[6]+f0L[6])+0.25*f0R[5]-0.25*(f0L[5]+f0R[4])+0.25*f0L[4]+0.4330127018922193*(f0R[3]+f0L[3])-0.25*f0R[2]+0.25*(f0L[2]+f0R[1])-0.25*(f0L[1]+f0R[0])+0.25*f0L[0])*sgn(alphaOrdR)+0.4330127018922193*f0R[15]-0.4330127018922193*(f0L[15]+f0R[14])+0.4330127018922193*(f0L[14]+f0R[13])-0.4330127018922193*f0L[13]-0.25*(f0R[12]+f0L[12])+0.4330127018922193*f0R[11]-0.4330127018922193*(f0L[11]+f0R[10])+0.4330127018922193*f0L[10]+0.25*(f0R[9]+f0L[9])-0.25*(f0R[8]+f0L[8])-0.4330127018922193*f0R[7]+0.4330127018922193*(f0L[7]+f0R[6])-0.4330127018922193*f0L[6]-0.25*(f0R[5]+f0L[5])+0.25*(f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]+0.25*(f0R[2]+f0L[2])-0.25*(f0R[1]+f0L[1])+0.25*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.3535533905932737*(alpha1R[1]+alpha1R[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(f0R[15]+f0L[15]+f0R[14]+f0L[14]+f0R[13]+f0L[13])-0.25*f0R[12]+0.25*f0L[12]+0.4330127018922193*(f0R[11]+f0L[11]+f0R[10]+f0L[10])-0.25*(f0R[9]+f0R[8])+0.25*(f0L[9]+f0L[8])+0.4330127018922193*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.25*(f0R[5]+f0R[4])+0.25*(f0L[5]+f0L[4])+0.4330127018922193*(f0R[3]+f0L[3])-0.25*(f0R[2]+f0R[1]+f0R[0])+0.25*(f0L[2]+f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.4330127018922193*(f0R[15]+f0R[14]+f0R[13])+0.4330127018922193*(f0L[15]+f0L[14]+f0L[13])+0.25*(f0R[12]+f0L[12])-0.4330127018922193*(f0R[11]+f0R[10])+0.4330127018922193*(f0L[11]+f0L[10])+0.25*(f0R[9]+f0L[9]+f0R[8]+f0L[8])-0.4330127018922193*(f0R[7]+f0R[6])+0.4330127018922193*(f0L[7]+f0L[6])+0.25*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.4330127018922193*f0R[3]+0.4330127018922193*f0L[3]+0.25*(f0R[2]+f0L[2]+f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.25*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = 0.25*(alpha1R[1]*fUp[4]+alpha1R[0]*fUp[2]); 
  incr[3] = -0.4330127018922193*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[4] = 0.25*(alpha1R[1]*fUp[5]+alpha1R[0]*fUp[3]); 
  incr[5] = 0.25*(alpha1R[0]*fUp[4]+alpha1R[1]*fUp[2]); 
  incr[6] = -0.4330127018922193*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[7] = -0.4330127018922193*(alpha1R[1]*fUp[4]+alpha1R[0]*fUp[2]); 
  incr[8] = 0.25*(alpha1R[0]*fUp[5]+alpha1R[1]*fUp[3]); 
  incr[9] = 0.25*(alpha1R[1]*fUp[7]+alpha1R[0]*fUp[6]); 
  incr[10] = -0.4330127018922193*(alpha1R[1]*fUp[5]+alpha1R[0]*fUp[3]); 
  incr[11] = -0.4330127018922193*(alpha1R[0]*fUp[4]+alpha1R[1]*fUp[2]); 
  incr[12] = 0.25*(alpha1R[0]*fUp[7]+alpha1R[1]*fUp[6]); 
  incr[13] = -0.4330127018922193*(alpha1R[0]*fUp[5]+alpha1R[1]*fUp[3]); 
  incr[14] = -0.4330127018922193*(alpha1R[1]*fUp[7]+alpha1R[0]*fUp[6]); 
  incr[15] = -0.4330127018922193*(alpha1R[0]*fUp[7]+alpha1R[1]*fUp[6]); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 
  outR[8] += incr[8]*rdvpar2R; 
  outR[9] += incr[9]*rdvpar2R; 
  outR[10] += incr[10]*rdvpar2R; 
  outR[11] += incr[11]*rdvpar2R; 
  outR[12] += incr[12]*rdvpar2R; 
  outR[13] += incr[13]*rdvpar2R; 
  outR[14] += incr[14]*rdvpar2R; 
  outR[15] += incr[15]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += -1.0*incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 
  outL[4] += -1.0*incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 
  outL[8] += -1.0*incr[8]*rdvpar2L; 
  outL[9] += -1.0*incr[9]*rdvpar2L; 
  outL[10] += incr[10]*rdvpar2L; 
  outL[11] += incr[11]*rdvpar2L; 
  outL[12] += -1.0*incr[12]*rdvpar2L; 
  outL[13] += incr[13]*rdvpar2L; 
  outL[14] += incr[14]*rdvpar2L; 
  outL[15] += incr[15]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
