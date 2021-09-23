#include <DeltaFGyrokineticModDecl.h>
double EmDeltaFGyrokineticGenGeoSurf1x2vSer_x_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // AparL: Apar in the left cell.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
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
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double BstarZdBmagEML[8]; 
  BstarZdBmagEML[0] = -1.732050807568877*b_y[0]*jacobTotInv[0]*AparL[1]*rdx2R; 

  double alpha0R[4]; 
  alpha0R[0] = (0.4330127018922193*BstarZdBmagR[0]*hamil0R[2]*rdvpar2R)/m_; 

  double alpha1R[4]; 
  alpha1R[0] = (0.2165063509461096*(BstarZdBmagEMR[0]+BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*(3.464101615137754*BstarZdBmagR[0]+1.732050807568877*BstarZdBmagEMR[0]+1.732050807568877*BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  double incr[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0]); 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.3535533905932737*f1R[6]+0.3535533905932737*f1L[6]-0.6123724356957944*(f1R[5]+f1L[5]+f1R[4]+f1L[4])+0.3535533905932737*(f1R[3]+f1R[2])-0.3535533905932737*(f1L[3]+f1L[2])+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*f1L[7]+0.3535533905932737*(f1R[6]+f1L[6])+0.6123724356957944*(f1R[5]+f1R[4])-0.6123724356957944*(f1L[5]+f1L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.3535533905932737*f1R[6]-0.3535533905932737*f1L[6]-0.6123724356957944*(f1R[5]+f1L[5])+0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*(f1L[3]+f1R[2])+0.3535533905932737*f1L[2]+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*f1L[7]-0.3535533905932737*(f1R[6]+f1L[6])+0.6123724356957944*f1R[5]-0.6123724356957944*(f1L[5]+f1R[4])+0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])+0.3535533905932737*(f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0]); 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.3535533905932737*f1R[6]-0.3535533905932737*f1L[6]+0.6123724356957944*(f1R[5]+f1L[5])-0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*(f1L[3]+f1R[2])-0.3535533905932737*f1L[2]+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*f1L[7]-0.3535533905932737*(f1R[6]+f1L[6])-0.6123724356957944*f1R[5]+0.6123724356957944*(f1L[5]+f1R[4])-0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.3535533905932737*(f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.3535533905932737*f1R[6]+0.3535533905932737*f1L[6]+0.6123724356957944*(f1R[5]+f1L[5]+f1R[4]+f1L[4])-0.3535533905932737*(f1R[3]+f1R[2])+0.3535533905932737*(f1L[3]+f1L[2])+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*f1L[7]+0.3535533905932737*(f1R[6]+f1L[6])-0.6123724356957944*(f1R[5]+f1R[4])+0.6123724356957944*(f1L[5]+f1L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alpha1R[0]+alpha0R[0])*fUp[0]; 
  incr[1] = -0.6123724356957944*(alpha1R[0]+alpha0R[0])*fUp[0]; 
  incr[2] = 0.3535533905932737*(alpha1R[0]+alpha0R[0])*fUp[1]; 
  incr[3] = 0.3535533905932737*(alpha1R[0]+alpha0R[0])*fUp[2]; 
  incr[4] = -0.6123724356957944*(alpha1R[0]+alpha0R[0])*fUp[1]; 
  incr[5] = -0.6123724356957944*(alpha1R[0]+alpha0R[0])*fUp[2]; 
  incr[6] = 0.3535533905932737*(alpha1R[0]+alpha0R[0])*fUp[3]; 
  incr[7] = -0.6123724356957944*(alpha1R[0]+alpha0R[0])*fUp[3]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += incr[7]*rdx2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])+0.3535533905932737*(f0R[3]+f0R[2])-0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*(f0R[5]+f0R[4])-0.6123724356957944*(f0L[5]+f0L[4])-0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5])+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*(f0L[3]+f0R[2])+0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*f0R[5]-0.6123724356957944*(f0L[5]+f0R[4])+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])+0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5])-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*(f0L[3]+f0R[2])-0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*f0R[5]+0.6123724356957944*(f0L[5]+f0R[4])-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.3535533905932737*(f0R[3]+f0R[2])+0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*(f0R[5]+f0R[4])+0.6123724356957944*(f0L[5]+f0L[4])+0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 

  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*alpha1R[0]*fUp[0]; 
  incr[1] = -0.6123724356957944*alpha1R[0]*fUp[0]; 
  incr[2] = 0.3535533905932737*alpha1R[0]*fUp[1]; 
  incr[3] = 0.3535533905932737*alpha1R[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alpha1R[0]*fUp[1]; 
  incr[5] = -0.6123724356957944*alpha1R[0]*fUp[2]; 
  incr[6] = 0.3535533905932737*alpha1R[0]*fUp[3]; 
  incr[7] = -0.6123724356957944*alpha1R[0]*fUp[3]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += incr[7]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSer_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // AparL: Apar in the left cell.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
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
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double alpha0R[4]; 

  double alpha0UpR[4]; 

  double alpha1R[4]; 
  alpha1R[0] = -(0.4330127018922193*BstarZdBmagR[0]*hamil1R[1]*rdx2R)/m_; 

  double alpha1UpR[4]; 
  alpha1UpR[0] = -(0.25*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+5.656854249492382*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = -(1.414213562373095*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+5.656854249492382*dApardtPrev[0]*q_))/m_; 

  double incr[8]; 
  double incrEmMod[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.6123724356957944*(f1R[6]+f1L[6])-0.3535533905932737*f1R[5]+0.3535533905932737*f1L[5]-0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])+0.3535533905932737*f1R[1]-0.3535533905932737*(f1L[1]+f1R[0])+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*(f1L[7]+f1R[6])-0.6123724356957944*f1L[6]+0.3535533905932737*(f1R[5]+f1L[5])+0.6123724356957944*f1R[4]-0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]-0.3535533905932737*(f1R[1]+f1L[1])+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f1R[6]+f1L[6]))+0.3535533905932737*f1R[5]-0.3535533905932737*f1L[5]+0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])-0.3535533905932737*(f1R[1]+f1R[0])+0.3535533905932737*(f1L[1]+f1L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f1R[6])-0.6123724356957944*(f1L[7]+f1L[6])-0.3535533905932737*(f1R[5]+f1L[5])-0.6123724356957944*f1R[4]+0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]+0.3535533905932737*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.6123724356957944*(f1R[6]+f1L[6])+0.3535533905932737*f1R[5]-0.3535533905932737*f1L[5]-0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])+0.3535533905932737*f1R[1]-0.3535533905932737*(f1L[1]+f1R[0])+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*(f1L[7]+f1R[6])+0.6123724356957944*f1L[6]-0.3535533905932737*(f1R[5]+f1L[5])+0.6123724356957944*f1R[4]-0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]-0.3535533905932737*(f1R[1]+f1L[1])+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f1R[6]+f1L[6])-0.3535533905932737*f1R[5]+0.3535533905932737*f1L[5]+0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])-0.3535533905932737*(f1R[1]+f1R[0])+0.3535533905932737*(f1L[1]+f1L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f1R[6])+0.6123724356957944*(f1L[7]+f1L[6])+0.3535533905932737*(f1R[5]+f1L[5])-0.6123724356957944*f1R[4]+0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]+0.3535533905932737*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = 0.7071067811865475*fUp[2]; 
  incrEmMod[4] = -1.224744871391589*fUp[1]; 
  incrEmMod[5] = 0.7071067811865475*fUp[3]; 
  incrEmMod[6] = -1.224744871391589*fUp[2]; 
  incrEmMod[7] = -1.224744871391589*fUp[3]; 

  incr[0] = 0.3535533905932737*alpha1R[0]*fUp[0]; 
  incr[1] = 0.3535533905932737*alpha1R[0]*fUp[1]; 
  incr[2] = -0.6123724356957944*alpha1R[0]*fUp[0]; 
  incr[3] = 0.3535533905932737*alpha1R[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alpha1R[0]*fUp[1]; 
  incr[5] = 0.3535533905932737*alpha1R[0]*fUp[3]; 
  incr[6] = -0.6123724356957944*alpha1R[0]*fUp[2]; 
  incr[7] = -0.6123724356957944*alpha1R[0]*fUp[3]; 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModR[4] += incrEmMod[4]*rdvpar2R; 
  emModR[5] += incrEmMod[5]*rdvpar2R; 
  emModR[6] += incrEmMod[6]*rdvpar2R; 
  emModR[7] += incrEmMod[7]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += -1.0*incrEmMod[3]*rdvpar2L; 
  emModL[4] += incrEmMod[4]*rdvpar2L; 
  emModL[5] += -1.0*incrEmMod[5]*rdvpar2L; 
  emModL[6] += incrEmMod[6]*rdvpar2L; 
  emModL[7] += incrEmMod[7]*rdvpar2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.6123724356957944*(f0R[6]+f0L[6])-0.3535533905932737*f0R[5]+0.3535533905932737*f0L[5]-0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])+0.3535533905932737*f0R[1]-0.3535533905932737*(f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*(f0L[7]+f0R[6])-0.6123724356957944*f0L[6]+0.3535533905932737*(f0R[5]+f0L[5])+0.6123724356957944*f0R[4]-0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]-0.3535533905932737*(f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]+f0R[6]+f0L[6]))+0.3535533905932737*f0R[5]-0.3535533905932737*f0L[5]+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])-0.3535533905932737*(f0R[1]+f0R[0])+0.3535533905932737*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f0R[7]+f0R[6])-0.6123724356957944*(f0L[7]+f0L[6])-0.3535533905932737*(f0R[5]+f0L[5])-0.6123724356957944*f0R[4]+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]+0.3535533905932737*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.6123724356957944*(f0R[6]+f0L[6])+0.3535533905932737*f0R[5]-0.3535533905932737*f0L[5]-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])+0.3535533905932737*f0R[1]-0.3535533905932737*(f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*(f0L[7]+f0R[6])+0.6123724356957944*f0L[6]-0.3535533905932737*(f0R[5]+f0L[5])+0.6123724356957944*f0R[4]-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]-0.3535533905932737*(f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.3535533905932737*f0R[5]+0.3535533905932737*f0L[5]+0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])-0.3535533905932737*(f0R[1]+f0R[0])+0.3535533905932737*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f0R[7]+f0R[6])+0.6123724356957944*(f0L[7]+f0L[6])+0.3535533905932737*(f0R[5]+f0L[5])-0.6123724356957944*f0R[4]+0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]+0.3535533905932737*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = 0.7071067811865475*fUp[2]; 
  incrEmMod[4] = -1.224744871391589*fUp[1]; 
  incrEmMod[5] = 0.7071067811865475*fUp[3]; 
  incrEmMod[6] = -1.224744871391589*fUp[2]; 
  incrEmMod[7] = -1.224744871391589*fUp[3]; 

  incr[0] = 0.3535533905932737*alpha1R[0]*fUp[0]; 
  incr[1] = 0.3535533905932737*alpha1R[0]*fUp[1]; 
  incr[2] = -0.6123724356957944*alpha1R[0]*fUp[0]; 
  incr[3] = 0.3535533905932737*alpha1R[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alpha1R[0]*fUp[1]; 
  incr[5] = 0.3535533905932737*alpha1R[0]*fUp[3]; 
  incr[6] = -0.6123724356957944*alpha1R[0]*fUp[2]; 
  incr[7] = -0.6123724356957944*alpha1R[0]*fUp[3]; 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModR[4] += incrEmMod[4]*rdvpar2R; 
  emModR[5] += incrEmMod[5]*rdvpar2R; 
  emModR[6] += incrEmMod[6]*rdvpar2R; 
  emModR[7] += incrEmMod[7]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += -1.0*incrEmMod[3]*rdvpar2L; 
  emModL[4] += incrEmMod[4]*rdvpar2L; 
  emModL[5] += -1.0*incrEmMod[5]*rdvpar2L; 
  emModL[6] += incrEmMod[6]*rdvpar2L; 
  emModL[7] += incrEmMod[7]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSerStep2_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.
  // emModL,emModR: .

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double alpha0R[4]; 

  double alpha0UpR[4]; 

  double alpha1R[4]; 
  alpha1R[0] = -(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha1R[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 

  double alpha1UpR[4]; 
  alpha1UpR[0] = -(0.25*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+5.656854249492382*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = -(1.414213562373095*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+5.656854249492382*dApardtPrev[0]*q_))/m_; 

  double incr[8]; 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7])-0.6123724356957944*(f1R[6]+f1L[6]+f0R[6]+f0L[6])-0.3535533905932737*(f1R[5]+f0R[5])+0.3535533905932737*(f1L[5]+f0L[5])-0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])+0.3535533905932737*(f1R[3]+f0R[3])-0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])+0.3535533905932737*(f1R[1]+f0R[1])-0.3535533905932737*(f1L[1]+f0L[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[0]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f0R[7])+0.6123724356957944*(f1L[7]+f0L[7]+f1R[6]+f0R[6])-0.6123724356957944*(f1L[6]+f0L[6])+0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])+0.6123724356957944*(f1R[4]+f0R[4])-0.6123724356957944*(f1L[4]+f0L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])-0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1])+0.3535533905932737*(f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]+f1R[6]+f1L[6]+f0R[6]+f0L[6]))+0.3535533905932737*(f1R[5]+f0R[5])-0.3535533905932737*(f1L[5]+f0L[5])+0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])+0.3535533905932737*(f1R[3]+f0R[3])-0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])-0.3535533905932737*(f1R[1]+f0R[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[1]+f0L[1]+f1L[0]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f0R[7]+f1R[6]+f0R[6])-0.6123724356957944*(f1L[7]+f0L[7]+f1L[6]+f0L[6])-0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])-0.6123724356957944*(f1R[4]+f0R[4])+0.6123724356957944*(f1L[4]+f0L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])+0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1]+f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]))+0.6123724356957944*(f1R[6]+f1L[6]+f0R[6]+f0L[6])+0.3535533905932737*(f1R[5]+f0R[5])-0.3535533905932737*(f1L[5]+f0L[5])-0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])-0.3535533905932737*(f1R[3]+f0R[3])+0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])+0.3535533905932737*(f1R[1]+f0R[1])-0.3535533905932737*(f1L[1]+f0L[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[0]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f0R[7])-0.6123724356957944*(f1L[7]+f0L[7]+f1R[6]+f0R[6])+0.6123724356957944*(f1L[6]+f0L[6])-0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])+0.6123724356957944*(f1R[4]+f0R[4])-0.6123724356957944*(f1L[4]+f0L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])-0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1])+0.3535533905932737*(f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]+f1R[6]+f1L[6]+f0R[6]+f0L[6])-0.3535533905932737*(f1R[5]+f0R[5])+0.3535533905932737*(f1L[5]+f0L[5])+0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])-0.3535533905932737*(f1R[3]+f0R[3])+0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])-0.3535533905932737*(f1R[1]+f0R[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[1]+f0L[1]+f1L[0]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f0R[7]+f1R[6]+f0R[6])+0.6123724356957944*(f1L[7]+f0L[7]+f1L[6]+f0L[6])+0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])-0.6123724356957944*(f1R[4]+f0R[4])+0.6123724356957944*(f1L[4]+f0L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])+0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1]+f1R[0]+f1L[0]+f0R[0]+f0L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.3535533905932737*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.6123724356957944*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = 0.3535533905932737*(alpha1R[1]*fUp[3]+alpha1R[0]*fUp[2]); 
  incr[4] = -0.6123724356957944*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[5] = 0.3535533905932737*(alpha1R[0]*fUp[3]+alpha1R[1]*fUp[2]); 
  incr[6] = -0.6123724356957944*(alpha1R[1]*fUp[3]+alpha1R[0]*fUp[2]); 
  incr[7] = -0.6123724356957944*(alpha1R[0]*fUp[3]+alpha1R[1]*fUp[2]); 


  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSer_x_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // AparL: Apar in the left cell.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
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
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[5] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.732050807568877*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double BstarZdBmagEML[8]; 
  BstarZdBmagEML[0] = -1.732050807568877*(2.0*AparL[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(b_y[0]*AparL[1]-1.0*AparL[0]*b_y[1]))*rdx2R; 
  BstarZdBmagEML[1] = 1.732050807568877*((AparL[0]*b_y[1]-1.0*b_y[0]*AparL[1])*jacobTotInv[1]-2.0*jacobTotInv[0]*AparL[1]*b_y[1])*rdx2R; 

  double alpha0R[4]; 
  alpha0R[0] = -(0.25*(3.0*BstarZdBmagR[1]-1.732050807568877*BstarZdBmagR[0])*hamil0R[2]*rdvpar2R)/m_; 
  alpha0R[1] = -(0.25*hamil0R[2]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[2])*rdvpar2R)/m_; 

  double alpha1R[4]; 
  alpha1R[0] = -(0.125*(3.0*(BstarZdBmagEMR[1]+BstarZdBmagEML[1])-1.732050807568877*(BstarZdBmagEMR[0]+BstarZdBmagEML[0]))*hamil0R[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(6.0*BstarZdBmagR[1]+3.0*BstarZdBmagEMR[1]+3.0*BstarZdBmagEML[1]-3.464101615137754*BstarZdBmagR[0]-1.732050807568877*BstarZdBmagEMR[0]-1.732050807568877*BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  double incr[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0])-0.5*alpha0R[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.3535533905932737*f1R[6]+0.3535533905932737*f1L[6]-0.6123724356957944*(f1R[5]+f1L[5]+f1R[4]+f1L[4])+0.3535533905932737*(f1R[3]+f1R[2])-0.3535533905932737*(f1L[3]+f1L[2])+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*f1L[7]+0.3535533905932737*(f1R[6]+f1L[6])+0.6123724356957944*(f1R[5]+f1R[4])-0.6123724356957944*(f1L[5]+f1L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha0R[1]+alpha1R[0]+alpha0R[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.3535533905932737*f1R[6]-0.3535533905932737*f1L[6]-0.6123724356957944*(f1R[5]+f1L[5])+0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*(f1L[3]+f1R[2])+0.3535533905932737*f1L[2]+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*f1L[7]-0.3535533905932737*(f1R[6]+f1L[6])+0.6123724356957944*f1R[5]-0.6123724356957944*(f1L[5]+f1R[4])+0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])+0.3535533905932737*(f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0])-0.5*alpha0R[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.3535533905932737*f1R[6]-0.3535533905932737*f1L[6]+0.6123724356957944*(f1R[5]+f1L[5])-0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*(f1L[3]+f1R[2])-0.3535533905932737*f1L[2]+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*f1L[7]-0.3535533905932737*(f1R[6]+f1L[6])-0.6123724356957944*f1R[5]+0.6123724356957944*(f1L[5]+f1R[4])-0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.3535533905932737*(f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha0R[1]+alpha1R[0]+alpha0R[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.3535533905932737*f1R[6]+0.3535533905932737*f1L[6]+0.6123724356957944*(f1R[5]+f1L[5]+f1R[4]+f1L[4])-0.3535533905932737*(f1R[3]+f1R[2])+0.3535533905932737*(f1L[3]+f1L[2])+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*f1L[7]+0.3535533905932737*(f1R[6]+f1L[6])-0.6123724356957944*(f1R[5]+f1R[4])+0.6123724356957944*(f1L[5]+f1L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alpha0R[1]*fUp[1]+(alpha1R[0]+alpha0R[0])*fUp[0]); 
  incr[1] = -0.6123724356957944*(alpha0R[1]*fUp[1]+(alpha1R[0]+alpha0R[0])*fUp[0]); 
  incr[2] = 0.3535533905932737*((alpha1R[0]+alpha0R[0])*fUp[1]+fUp[0]*alpha0R[1]); 
  incr[3] = 0.3535533905932737*(alpha0R[1]*fUp[3]+(alpha1R[0]+alpha0R[0])*fUp[2]); 
  incr[4] = -0.6123724356957944*((alpha1R[0]+alpha0R[0])*fUp[1]+fUp[0]*alpha0R[1]); 
  incr[5] = -0.6123724356957944*(alpha0R[1]*fUp[3]+(alpha1R[0]+alpha0R[0])*fUp[2]); 
  incr[6] = 0.3535533905932737*((alpha1R[0]+alpha0R[0])*fUp[3]+alpha0R[1]*fUp[2]); 
  incr[7] = -0.6123724356957944*((alpha1R[0]+alpha0R[0])*fUp[3]+alpha0R[1]*fUp[2]); 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += incr[7]*rdx2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])+0.3535533905932737*(f0R[3]+f0R[2])-0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*(f0R[5]+f0R[4])-0.6123724356957944*(f0L[5]+f0L[4])-0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5])+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*(f0L[3]+f0R[2])+0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*f0R[5]-0.6123724356957944*(f0L[5]+f0R[4])+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])+0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5])-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*(f0L[3]+f0R[2])-0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*f0R[5]+0.6123724356957944*(f0L[5]+f0R[4])-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.3535533905932737*(f0R[3]+f0R[2])+0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*(f0R[5]+f0R[4])+0.6123724356957944*(f0L[5]+f0L[4])+0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 

  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*alpha1R[0]*fUp[0]; 
  incr[1] = -0.6123724356957944*alpha1R[0]*fUp[0]; 
  incr[2] = 0.3535533905932737*alpha1R[0]*fUp[1]; 
  incr[3] = 0.3535533905932737*alpha1R[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alpha1R[0]*fUp[1]; 
  incr[5] = -0.6123724356957944*alpha1R[0]*fUp[2]; 
  incr[6] = 0.3535533905932737*alpha1R[0]*fUp[3]; 
  incr[7] = -0.6123724356957944*alpha1R[0]*fUp[3]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += incr[7]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSer_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // AparL: Apar in the left cell.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
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
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[5] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.732050807568877*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double alpha0R[4]; 
  alpha0R[0] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alpha0R[1] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 
  alpha0R[2] = (0.25*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*hamil0R[5]*rdx2R)/m_; 
  alpha0R[3] = (0.25*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*hamil0R[5]*rdx2R)/m_; 

  double alpha0UpR[4]; 
  alpha0UpR[0] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alpha0UpR[1] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 
  alpha0UpR[2] = (0.25*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*hamil0R[5]*rdx2R)/m_; 
  alpha0UpR[3] = (0.25*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*hamil0R[5]*rdx2R)/m_; 

  double alpha1R[4]; 
  alpha1R[0] = (0.25*(3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*(BstarZdBmagR[0]*hamil1R[1]+BstarZdBmagEMR[0]*hamil0R[1]))*rdx2R)/m_; 
  alpha1R[1] = (0.25*(3.0*hamil1R[1]*BstarZdBmagR[4]-1.732050807568877*(BstarZdBmagR[1]*hamil1R[1]+BstarZdBmagEMR[1]*hamil0R[1]))*rdx2R)/m_; 
  alpha1R[2] = -(0.4330127018922193*BstarZdBmagEMR[0]*hamil0R[5]*rdx2R)/m_; 
  alpha1R[3] = -(0.4330127018922193*BstarZdBmagEMR[1]*hamil0R[5]*rdx2R)/m_; 

  double alpha1UpR[4]; 
  alpha1UpR[0] = (0.25*((3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*(BstarZdBmagR[0]*hamil1R[1]+BstarZdBmagEMR[0]*hamil0R[1]))*rdx2R-5.656854249492382*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = (0.25*((3.0*hamil1R[1]*BstarZdBmagR[4]-1.732050807568877*(BstarZdBmagR[1]*hamil1R[1]+BstarZdBmagEMR[1]*hamil0R[1]))*rdx2R-5.656854249492382*dApardtPrev[1]*q_))/m_; 
  alpha1UpR[2] = -(0.4330127018922193*BstarZdBmagEMR[0]*hamil0R[5]*rdx2R)/m_; 
  alpha1UpR[3] = -(0.4330127018922193*BstarZdBmagEMR[1]*hamil0R[5]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(((3.0*hamil1R[1]+3.0*hamil0R[1])*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamil1R[1]+((-1.732050807568877*BstarZdBmagR[0])-1.732050807568877*BstarZdBmagEMR[0])*hamil0R[1])*rdx2R-5.656854249492382*dApardtPrev[0]*q_))/m_; 

  double incr[8]; 
  double incrEmMod[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*(alpha1UpR[3]+alpha0UpR[3])-0.5*(alpha1UpR[2]+alpha0UpR[2]+alpha1UpR[1]+alpha0UpR[1])+0.5*(alpha1UpR[0]+alpha0UpR[0]); 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.6123724356957944*(f1R[6]+f1L[6])-0.3535533905932737*f1R[5]+0.3535533905932737*f1L[5]-0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])+0.3535533905932737*f1R[1]-0.3535533905932737*(f1L[1]+f1R[0])+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*(f1L[7]+f1R[6])-0.6123724356957944*f1L[6]+0.3535533905932737*(f1R[5]+f1L[5])+0.6123724356957944*f1R[4]-0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]-0.3535533905932737*(f1R[1]+f1L[1])+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha0UpR[1]+alpha1UpR[0]+alpha0UpR[0])-0.5*(alpha1UpR[3]+alpha0UpR[3]+alpha1UpR[2]+alpha0UpR[2]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f1R[6]+f1L[6]))+0.3535533905932737*f1R[5]-0.3535533905932737*f1L[5]+0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])-0.3535533905932737*(f1R[1]+f1R[0])+0.3535533905932737*(f1L[1]+f1L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f1R[6])-0.6123724356957944*(f1L[7]+f1L[6])-0.3535533905932737*(f1R[5]+f1L[5])-0.6123724356957944*f1R[4]+0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]+0.3535533905932737*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 
  alphaOrdR = (-0.5*(alpha1UpR[3]+alpha0UpR[3]))+0.5*(alpha1UpR[2]+alpha0UpR[2])-0.5*(alpha1UpR[1]+alpha0UpR[1])+0.5*(alpha1UpR[0]+alpha0UpR[0]); 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.6123724356957944*(f1R[6]+f1L[6])+0.3535533905932737*f1R[5]-0.3535533905932737*f1L[5]-0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])+0.3535533905932737*f1R[1]-0.3535533905932737*(f1L[1]+f1R[0])+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*(f1L[7]+f1R[6])+0.6123724356957944*f1L[6]-0.3535533905932737*(f1R[5]+f1L[5])+0.6123724356957944*f1R[4]-0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]-0.3535533905932737*(f1R[1]+f1L[1])+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[3]+alpha0UpR[3]+alpha1UpR[2]+alpha0UpR[2]+alpha1UpR[1]+alpha0UpR[1]+alpha1UpR[0]+alpha0UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f1R[6]+f1L[6])-0.3535533905932737*f1R[5]+0.3535533905932737*f1L[5]+0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])-0.3535533905932737*(f1R[1]+f1R[0])+0.3535533905932737*(f1L[1]+f1L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f1R[6])+0.6123724356957944*(f1L[7]+f1L[6])+0.3535533905932737*(f1R[5]+f1L[5])-0.6123724356957944*f1R[4]+0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]+0.3535533905932737*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = 0.7071067811865475*fUp[2]; 
  incrEmMod[4] = -1.224744871391589*fUp[1]; 
  incrEmMod[5] = 0.7071067811865475*fUp[3]; 
  incrEmMod[6] = -1.224744871391589*fUp[2]; 
  incrEmMod[7] = -1.224744871391589*fUp[3]; 

  incr[0] = 0.3535533905932737*((alpha1R[3]+alpha0R[3])*fUp[3]+(alpha1R[2]+alpha0R[2])*fUp[2]+(alpha1R[1]+alpha0R[1])*fUp[1]+(alpha1R[0]+alpha0R[0])*fUp[0]); 
  incr[1] = 0.3535533905932737*((alpha1R[2]+alpha0R[2])*fUp[3]+fUp[2]*(alpha1R[3]+alpha0R[3])+(alpha1R[0]+alpha0R[0])*fUp[1]+fUp[0]*(alpha1R[1]+alpha0R[1])); 
  incr[2] = -0.6123724356957944*((alpha1R[3]+alpha0R[3])*fUp[3]+(alpha1R[2]+alpha0R[2])*fUp[2]+(alpha1R[1]+alpha0R[1])*fUp[1]+(alpha1R[0]+alpha0R[0])*fUp[0]); 
  incr[3] = 0.3535533905932737*((alpha1R[1]+alpha0R[1])*fUp[3]+fUp[1]*(alpha1R[3]+alpha0R[3])+(alpha1R[0]+alpha0R[0])*fUp[2]+fUp[0]*(alpha1R[2]+alpha0R[2])); 
  incr[4] = -0.6123724356957944*((alpha1R[2]+alpha0R[2])*fUp[3]+fUp[2]*(alpha1R[3]+alpha0R[3])+(alpha1R[0]+alpha0R[0])*fUp[1]+fUp[0]*(alpha1R[1]+alpha0R[1])); 
  incr[5] = 0.3535533905932737*((alpha1R[0]+alpha0R[0])*fUp[3]+fUp[0]*(alpha1R[3]+alpha0R[3])+(alpha1R[1]+alpha0R[1])*fUp[2]+fUp[1]*(alpha1R[2]+alpha0R[2])); 
  incr[6] = -0.6123724356957944*((alpha1R[1]+alpha0R[1])*fUp[3]+fUp[1]*(alpha1R[3]+alpha0R[3])+(alpha1R[0]+alpha0R[0])*fUp[2]+fUp[0]*(alpha1R[2]+alpha0R[2])); 
  incr[7] = -0.6123724356957944*((alpha1R[0]+alpha0R[0])*fUp[3]+fUp[0]*(alpha1R[3]+alpha0R[3])+(alpha1R[1]+alpha0R[1])*fUp[2]+fUp[1]*(alpha1R[2]+alpha0R[2])); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModR[4] += incrEmMod[4]*rdvpar2R; 
  emModR[5] += incrEmMod[5]*rdvpar2R; 
  emModR[6] += incrEmMod[6]*rdvpar2R; 
  emModR[7] += incrEmMod[7]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += -1.0*incrEmMod[3]*rdvpar2L; 
  emModL[4] += incrEmMod[4]*rdvpar2L; 
  emModL[5] += -1.0*incrEmMod[5]*rdvpar2L; 
  emModL[6] += incrEmMod[6]*rdvpar2L; 
  emModL[7] += incrEmMod[7]*rdvpar2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.5*alpha1UpR[3]-0.5*(alpha1UpR[2]+alpha1UpR[1])+0.5*alpha1UpR[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.6123724356957944*(f0R[6]+f0L[6])-0.3535533905932737*f0R[5]+0.3535533905932737*f0L[5]-0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])+0.3535533905932737*f0R[1]-0.3535533905932737*(f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*(f0L[7]+f0R[6])-0.6123724356957944*f0L[6]+0.3535533905932737*(f0R[5]+f0L[5])+0.6123724356957944*f0R[4]-0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]-0.3535533905932737*(f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0])-0.5*(alpha1UpR[3]+alpha1UpR[2]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]+f0R[6]+f0L[6]))+0.3535533905932737*f0R[5]-0.3535533905932737*f0L[5]+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])-0.3535533905932737*(f0R[1]+f0R[0])+0.3535533905932737*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f0R[7]+f0R[6])-0.6123724356957944*(f0L[7]+f0L[6])-0.3535533905932737*(f0R[5]+f0L[5])-0.6123724356957944*f0R[4]+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]+0.3535533905932737*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = (-0.5*alpha1UpR[3])+0.5*alpha1UpR[2]-0.5*alpha1UpR[1]+0.5*alpha1UpR[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.6123724356957944*(f0R[6]+f0L[6])+0.3535533905932737*f0R[5]-0.3535533905932737*f0L[5]-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])+0.3535533905932737*f0R[1]-0.3535533905932737*(f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*(f0L[7]+f0R[6])+0.6123724356957944*f0L[6]-0.3535533905932737*(f0R[5]+f0L[5])+0.6123724356957944*f0R[4]-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]-0.3535533905932737*(f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[3]+alpha1UpR[2]+alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.3535533905932737*f0R[5]+0.3535533905932737*f0L[5]+0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])-0.3535533905932737*(f0R[1]+f0R[0])+0.3535533905932737*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f0R[7]+f0R[6])+0.6123724356957944*(f0L[7]+f0L[6])+0.3535533905932737*(f0R[5]+f0L[5])-0.6123724356957944*f0R[4]+0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]+0.3535533905932737*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = 0.7071067811865475*fUp[2]; 
  incrEmMod[4] = -1.224744871391589*fUp[1]; 
  incrEmMod[5] = 0.7071067811865475*fUp[3]; 
  incrEmMod[6] = -1.224744871391589*fUp[2]; 
  incrEmMod[7] = -1.224744871391589*fUp[3]; 

  incr[0] = 0.3535533905932737*(alpha1R[3]*fUp[3]+alpha1R[2]*fUp[2]+alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.3535533905932737*(alpha1R[2]*fUp[3]+fUp[2]*alpha1R[3]+alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.6123724356957944*(alpha1R[3]*fUp[3]+alpha1R[2]*fUp[2]+alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = 0.3535533905932737*(alpha1R[1]*fUp[3]+fUp[1]*alpha1R[3]+alpha1R[0]*fUp[2]+fUp[0]*alpha1R[2]); 
  incr[4] = -0.6123724356957944*(alpha1R[2]*fUp[3]+fUp[2]*alpha1R[3]+alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[5] = 0.3535533905932737*(alpha1R[0]*fUp[3]+fUp[0]*alpha1R[3]+alpha1R[1]*fUp[2]+fUp[1]*alpha1R[2]); 
  incr[6] = -0.6123724356957944*(alpha1R[1]*fUp[3]+fUp[1]*alpha1R[3]+alpha1R[0]*fUp[2]+fUp[0]*alpha1R[2]); 
  incr[7] = -0.6123724356957944*(alpha1R[0]*fUp[3]+fUp[0]*alpha1R[3]+alpha1R[1]*fUp[2]+fUp[1]*alpha1R[2]); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModR[4] += incrEmMod[4]*rdvpar2R; 
  emModR[5] += incrEmMod[5]*rdvpar2R; 
  emModR[6] += incrEmMod[6]*rdvpar2R; 
  emModR[7] += incrEmMod[7]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += -1.0*incrEmMod[3]*rdvpar2L; 
  emModL[4] += incrEmMod[4]*rdvpar2L; 
  emModL[5] += -1.0*incrEmMod[5]*rdvpar2L; 
  emModL[6] += incrEmMod[6]*rdvpar2L; 
  emModL[7] += incrEmMod[7]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSerStep2_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.
  // emModL,emModR: .

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[5] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.732050807568877*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double alpha0R[4]; 

  double alpha0UpR[4]; 
  alpha0UpR[0] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alpha0UpR[1] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 
  alpha0UpR[2] = (0.25*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*hamil0R[5]*rdx2R)/m_; 
  alpha0UpR[3] = (0.25*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*hamil0R[5]*rdx2R)/m_; 

  double alpha1R[4]; 
  alpha1R[0] = -(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha1R[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 

  double alpha1UpR[4]; 
  alpha1UpR[0] = (0.25*((3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*(BstarZdBmagR[0]*hamil1R[1]+BstarZdBmagEMR[0]*hamil0R[1]))*rdx2R-5.656854249492382*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = (0.25*((3.0*hamil1R[1]*BstarZdBmagR[4]-1.732050807568877*(BstarZdBmagR[1]*hamil1R[1]+BstarZdBmagEMR[1]*hamil0R[1]))*rdx2R-5.656854249492382*dApardtPrev[1]*q_))/m_; 
  alpha1UpR[2] = -(0.4330127018922193*BstarZdBmagEMR[0]*hamil0R[5]*rdx2R)/m_; 
  alpha1UpR[3] = -(0.4330127018922193*BstarZdBmagEMR[1]*hamil0R[5]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(((3.0*hamil1R[1]+3.0*hamil0R[1])*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamil1R[1]+((-1.732050807568877*BstarZdBmagR[0])-1.732050807568877*BstarZdBmagEMR[0])*hamil0R[1])*rdx2R-5.656854249492382*dApardtPrev[0]*q_))/m_; 

  double incr[8]; 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alpha1UpR[3]-0.5*(alpha1UpR[2]+alpha1UpR[1])+0.5*alpha1UpR[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7])-0.6123724356957944*(f1R[6]+f1L[6]+f0R[6]+f0L[6])-0.3535533905932737*(f1R[5]+f0R[5])+0.3535533905932737*(f1L[5]+f0L[5])-0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])+0.3535533905932737*(f1R[3]+f0R[3])-0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])+0.3535533905932737*(f1R[1]+f0R[1])-0.3535533905932737*(f1L[1]+f0L[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[0]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f0R[7])+0.6123724356957944*(f1L[7]+f0L[7]+f1R[6]+f0R[6])-0.6123724356957944*(f1L[6]+f0L[6])+0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])+0.6123724356957944*(f1R[4]+f0R[4])-0.6123724356957944*(f1L[4]+f0L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])-0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1])+0.3535533905932737*(f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0])-0.5*(alpha1UpR[3]+alpha1UpR[2]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]+f1R[6]+f1L[6]+f0R[6]+f0L[6]))+0.3535533905932737*(f1R[5]+f0R[5])-0.3535533905932737*(f1L[5]+f0L[5])+0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])+0.3535533905932737*(f1R[3]+f0R[3])-0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])-0.3535533905932737*(f1R[1]+f0R[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[1]+f0L[1]+f1L[0]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f0R[7]+f1R[6]+f0R[6])-0.6123724356957944*(f1L[7]+f0L[7]+f1L[6]+f0L[6])-0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])-0.6123724356957944*(f1R[4]+f0R[4])+0.6123724356957944*(f1L[4]+f0L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])+0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1]+f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = (-0.5*alpha1UpR[3])+0.5*alpha1UpR[2]-0.5*alpha1UpR[1]+0.5*alpha1UpR[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]))+0.6123724356957944*(f1R[6]+f1L[6]+f0R[6]+f0L[6])+0.3535533905932737*(f1R[5]+f0R[5])-0.3535533905932737*(f1L[5]+f0L[5])-0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])-0.3535533905932737*(f1R[3]+f0R[3])+0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])+0.3535533905932737*(f1R[1]+f0R[1])-0.3535533905932737*(f1L[1]+f0L[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[0]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f0R[7])-0.6123724356957944*(f1L[7]+f0L[7]+f1R[6]+f0R[6])+0.6123724356957944*(f1L[6]+f0L[6])-0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])+0.6123724356957944*(f1R[4]+f0R[4])-0.6123724356957944*(f1L[4]+f0L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])-0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1])+0.3535533905932737*(f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[3]+alpha1UpR[2]+alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]+f1R[6]+f1L[6]+f0R[6]+f0L[6])-0.3535533905932737*(f1R[5]+f0R[5])+0.3535533905932737*(f1L[5]+f0L[5])+0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])-0.3535533905932737*(f1R[3]+f0R[3])+0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])-0.3535533905932737*(f1R[1]+f0R[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[1]+f0L[1]+f1L[0]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f0R[7]+f1R[6]+f0R[6])+0.6123724356957944*(f1L[7]+f0L[7]+f1L[6]+f0L[6])+0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])-0.6123724356957944*(f1R[4]+f0R[4])+0.6123724356957944*(f1L[4]+f0L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])+0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1]+f1R[0]+f1L[0]+f0R[0]+f0L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.3535533905932737*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.6123724356957944*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = 0.3535533905932737*(alpha1R[1]*fUp[3]+alpha1R[0]*fUp[2]); 
  incr[4] = -0.6123724356957944*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[5] = 0.3535533905932737*(alpha1R[0]*fUp[3]+alpha1R[1]*fUp[2]); 
  incr[6] = -0.6123724356957944*(alpha1R[1]*fUp[3]+alpha1R[0]*fUp[2]); 
  incr[7] = -0.6123724356957944*(alpha1R[0]*fUp[3]+alpha1R[1]*fUp[2]); 


  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSer_x_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // AparL: Apar in the left cell.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
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
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double BstarZdBmagEML[8]; 
  BstarZdBmagEML[0] = -1.732050807568877*b_y[0]*jacobTotInv[0]*AparL[1]*rdx2R; 

  double alpha0R[4]; 
  alpha0R[0] = (0.4330127018922193*BstarZdBmagR[0]*hamil0R[2]*rdvpar2R)/m_; 

  double alpha1R[4]; 
  alpha1R[0] = (0.2165063509461096*(BstarZdBmagEMR[0]+BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*(3.464101615137754*BstarZdBmagR[0]+1.732050807568877*BstarZdBmagEMR[0]+1.732050807568877*BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  double incr[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0]); 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.3535533905932737*f1R[6]+0.3535533905932737*f1L[6]-0.6123724356957944*(f1R[5]+f1L[5]+f1R[4]+f1L[4])+0.3535533905932737*(f1R[3]+f1R[2])-0.3535533905932737*(f1L[3]+f1L[2])+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*f1L[7]+0.3535533905932737*(f1R[6]+f1L[6])+0.6123724356957944*(f1R[5]+f1R[4])-0.6123724356957944*(f1L[5]+f1L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.3535533905932737*f1R[6]-0.3535533905932737*f1L[6]-0.6123724356957944*(f1R[5]+f1L[5])+0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*(f1L[3]+f1R[2])+0.3535533905932737*f1L[2]+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*f1L[7]-0.3535533905932737*(f1R[6]+f1L[6])+0.6123724356957944*f1R[5]-0.6123724356957944*(f1L[5]+f1R[4])+0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])+0.3535533905932737*(f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0]); 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.3535533905932737*f1R[6]-0.3535533905932737*f1L[6]+0.6123724356957944*(f1R[5]+f1L[5])-0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*(f1L[3]+f1R[2])-0.3535533905932737*f1L[2]+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*f1L[7]-0.3535533905932737*(f1R[6]+f1L[6])-0.6123724356957944*f1R[5]+0.6123724356957944*(f1L[5]+f1R[4])-0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.3535533905932737*(f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.3535533905932737*f1R[6]+0.3535533905932737*f1L[6]+0.6123724356957944*(f1R[5]+f1L[5]+f1R[4]+f1L[4])-0.3535533905932737*(f1R[3]+f1R[2])+0.3535533905932737*(f1L[3]+f1L[2])+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*f1L[7]+0.3535533905932737*(f1R[6]+f1L[6])-0.6123724356957944*(f1R[5]+f1R[4])+0.6123724356957944*(f1L[5]+f1L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alpha1R[0]+alpha0R[0])*fUp[0]; 
  incr[1] = -0.6123724356957944*(alpha1R[0]+alpha0R[0])*fUp[0]; 
  incr[2] = 0.3535533905932737*(alpha1R[0]+alpha0R[0])*fUp[1]; 
  incr[3] = 0.3535533905932737*(alpha1R[0]+alpha0R[0])*fUp[2]; 
  incr[4] = -0.6123724356957944*(alpha1R[0]+alpha0R[0])*fUp[1]; 
  incr[5] = -0.6123724356957944*(alpha1R[0]+alpha0R[0])*fUp[2]; 
  incr[6] = 0.3535533905932737*(alpha1R[0]+alpha0R[0])*fUp[3]; 
  incr[7] = -0.6123724356957944*(alpha1R[0]+alpha0R[0])*fUp[3]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += incr[7]*rdx2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])+0.3535533905932737*(f0R[3]+f0R[2])-0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*(f0R[5]+f0R[4])-0.6123724356957944*(f0L[5]+f0L[4])-0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5])+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*(f0L[3]+f0R[2])+0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*f0R[5]-0.6123724356957944*(f0L[5]+f0R[4])+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])+0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5])-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*(f0L[3]+f0R[2])-0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*f0R[5]+0.6123724356957944*(f0L[5]+f0R[4])-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.3535533905932737*(f0R[3]+f0R[2])+0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*(f0R[5]+f0R[4])+0.6123724356957944*(f0L[5]+f0L[4])+0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 

  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*alpha1R[0]*fUp[0]; 
  incr[1] = -0.6123724356957944*alpha1R[0]*fUp[0]; 
  incr[2] = 0.3535533905932737*alpha1R[0]*fUp[1]; 
  incr[3] = 0.3535533905932737*alpha1R[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alpha1R[0]*fUp[1]; 
  incr[5] = -0.6123724356957944*alpha1R[0]*fUp[2]; 
  incr[6] = 0.3535533905932737*alpha1R[0]*fUp[3]; 
  incr[7] = -0.6123724356957944*alpha1R[0]*fUp[3]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += incr[7]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSer_vpar_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // AparL: Apar in the left cell.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
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
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double alpha0R[4]; 

  double alpha0UpR[4]; 

  double alpha1R[4]; 
  alpha1R[0] = -(0.4330127018922193*BstarZdBmagR[0]*hamil1R[1]*rdx2R)/m_; 

  double alpha1UpR[4]; 
  alpha1UpR[0] = -(0.25*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+5.656854249492382*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = -(1.414213562373095*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+5.656854249492382*dApardtPrev[0]*q_))/m_; 

  double incr[8]; 
  double incrEmMod[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.6123724356957944*(f1R[6]+f1L[6])-0.3535533905932737*f1R[5]+0.3535533905932737*f1L[5]-0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])+0.3535533905932737*f1R[1]-0.3535533905932737*(f1L[1]+f1R[0])+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*(f1L[7]+f1R[6])-0.6123724356957944*f1L[6]+0.3535533905932737*(f1R[5]+f1L[5])+0.6123724356957944*f1R[4]-0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]-0.3535533905932737*(f1R[1]+f1L[1])+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f1R[6]+f1L[6]))+0.3535533905932737*f1R[5]-0.3535533905932737*f1L[5]+0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])-0.3535533905932737*(f1R[1]+f1R[0])+0.3535533905932737*(f1L[1]+f1L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f1R[6])-0.6123724356957944*(f1L[7]+f1L[6])-0.3535533905932737*(f1R[5]+f1L[5])-0.6123724356957944*f1R[4]+0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]+0.3535533905932737*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.6123724356957944*(f1R[6]+f1L[6])+0.3535533905932737*f1R[5]-0.3535533905932737*f1L[5]-0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])+0.3535533905932737*f1R[1]-0.3535533905932737*(f1L[1]+f1R[0])+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*(f1L[7]+f1R[6])+0.6123724356957944*f1L[6]-0.3535533905932737*(f1R[5]+f1L[5])+0.6123724356957944*f1R[4]-0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]-0.3535533905932737*(f1R[1]+f1L[1])+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f1R[6]+f1L[6])-0.3535533905932737*f1R[5]+0.3535533905932737*f1L[5]+0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])-0.3535533905932737*(f1R[1]+f1R[0])+0.3535533905932737*(f1L[1]+f1L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f1R[6])+0.6123724356957944*(f1L[7]+f1L[6])+0.3535533905932737*(f1R[5]+f1L[5])-0.6123724356957944*f1R[4]+0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]+0.3535533905932737*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = 0.7071067811865475*fUp[2]; 
  incrEmMod[4] = -1.224744871391589*fUp[1]; 
  incrEmMod[5] = 0.7071067811865475*fUp[3]; 
  incrEmMod[6] = -1.224744871391589*fUp[2]; 
  incrEmMod[7] = -1.224744871391589*fUp[3]; 

  incr[0] = 0.3535533905932737*alpha1R[0]*fUp[0]; 
  incr[1] = 0.3535533905932737*alpha1R[0]*fUp[1]; 
  incr[2] = -0.6123724356957944*alpha1R[0]*fUp[0]; 
  incr[3] = 0.3535533905932737*alpha1R[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alpha1R[0]*fUp[1]; 
  incr[5] = 0.3535533905932737*alpha1R[0]*fUp[3]; 
  incr[6] = -0.6123724356957944*alpha1R[0]*fUp[2]; 
  incr[7] = -0.6123724356957944*alpha1R[0]*fUp[3]; 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModR[4] += incrEmMod[4]*rdvpar2R; 
  emModR[5] += incrEmMod[5]*rdvpar2R; 
  emModR[6] += incrEmMod[6]*rdvpar2R; 
  emModR[7] += incrEmMod[7]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += -1.0*incrEmMod[3]*rdvpar2L; 
  emModL[4] += incrEmMod[4]*rdvpar2L; 
  emModL[5] += -1.0*incrEmMod[5]*rdvpar2L; 
  emModL[6] += incrEmMod[6]*rdvpar2L; 
  emModL[7] += incrEmMod[7]*rdvpar2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.6123724356957944*(f0R[6]+f0L[6])-0.3535533905932737*f0R[5]+0.3535533905932737*f0L[5]-0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])+0.3535533905932737*f0R[1]-0.3535533905932737*(f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*(f0L[7]+f0R[6])-0.6123724356957944*f0L[6]+0.3535533905932737*(f0R[5]+f0L[5])+0.6123724356957944*f0R[4]-0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]-0.3535533905932737*(f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]+f0R[6]+f0L[6]))+0.3535533905932737*f0R[5]-0.3535533905932737*f0L[5]+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])-0.3535533905932737*(f0R[1]+f0R[0])+0.3535533905932737*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f0R[7]+f0R[6])-0.6123724356957944*(f0L[7]+f0L[6])-0.3535533905932737*(f0R[5]+f0L[5])-0.6123724356957944*f0R[4]+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]+0.3535533905932737*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.6123724356957944*(f0R[6]+f0L[6])+0.3535533905932737*f0R[5]-0.3535533905932737*f0L[5]-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])+0.3535533905932737*f0R[1]-0.3535533905932737*(f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*(f0L[7]+f0R[6])+0.6123724356957944*f0L[6]-0.3535533905932737*(f0R[5]+f0L[5])+0.6123724356957944*f0R[4]-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]-0.3535533905932737*(f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.3535533905932737*f0R[5]+0.3535533905932737*f0L[5]+0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])-0.3535533905932737*(f0R[1]+f0R[0])+0.3535533905932737*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f0R[7]+f0R[6])+0.6123724356957944*(f0L[7]+f0L[6])+0.3535533905932737*(f0R[5]+f0L[5])-0.6123724356957944*f0R[4]+0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]+0.3535533905932737*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = 0.7071067811865475*fUp[2]; 
  incrEmMod[4] = -1.224744871391589*fUp[1]; 
  incrEmMod[5] = 0.7071067811865475*fUp[3]; 
  incrEmMod[6] = -1.224744871391589*fUp[2]; 
  incrEmMod[7] = -1.224744871391589*fUp[3]; 

  incr[0] = 0.3535533905932737*alpha1R[0]*fUp[0]; 
  incr[1] = 0.3535533905932737*alpha1R[0]*fUp[1]; 
  incr[2] = -0.6123724356957944*alpha1R[0]*fUp[0]; 
  incr[3] = 0.3535533905932737*alpha1R[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alpha1R[0]*fUp[1]; 
  incr[5] = 0.3535533905932737*alpha1R[0]*fUp[3]; 
  incr[6] = -0.6123724356957944*alpha1R[0]*fUp[2]; 
  incr[7] = -0.6123724356957944*alpha1R[0]*fUp[3]; 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModR[4] += incrEmMod[4]*rdvpar2R; 
  emModR[5] += incrEmMod[5]*rdvpar2R; 
  emModR[6] += incrEmMod[6]*rdvpar2R; 
  emModR[7] += incrEmMod[7]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += -1.0*incrEmMod[3]*rdvpar2L; 
  emModL[4] += incrEmMod[4]*rdvpar2L; 
  emModL[5] += -1.0*incrEmMod[5]*rdvpar2L; 
  emModL[6] += incrEmMod[6]*rdvpar2L; 
  emModL[7] += incrEmMod[7]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSerStep2_vpar_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.
  // emModL,emModR: .

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double alpha0R[4]; 

  double alpha0UpR[4]; 

  double alpha1R[4]; 
  alpha1R[0] = -(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha1R[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 

  double alpha1UpR[4]; 
  alpha1UpR[0] = -(0.25*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+5.656854249492382*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = -(1.414213562373095*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+5.656854249492382*dApardtPrev[0]*q_))/m_; 

  double incr[8]; 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7])-0.6123724356957944*(f1R[6]+f1L[6]+f0R[6]+f0L[6])-0.3535533905932737*(f1R[5]+f0R[5])+0.3535533905932737*(f1L[5]+f0L[5])-0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])+0.3535533905932737*(f1R[3]+f0R[3])-0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])+0.3535533905932737*(f1R[1]+f0R[1])-0.3535533905932737*(f1L[1]+f0L[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[0]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f0R[7])+0.6123724356957944*(f1L[7]+f0L[7]+f1R[6]+f0R[6])-0.6123724356957944*(f1L[6]+f0L[6])+0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])+0.6123724356957944*(f1R[4]+f0R[4])-0.6123724356957944*(f1L[4]+f0L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])-0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1])+0.3535533905932737*(f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]+f1R[6]+f1L[6]+f0R[6]+f0L[6]))+0.3535533905932737*(f1R[5]+f0R[5])-0.3535533905932737*(f1L[5]+f0L[5])+0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])+0.3535533905932737*(f1R[3]+f0R[3])-0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])-0.3535533905932737*(f1R[1]+f0R[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[1]+f0L[1]+f1L[0]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f0R[7]+f1R[6]+f0R[6])-0.6123724356957944*(f1L[7]+f0L[7]+f1L[6]+f0L[6])-0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])-0.6123724356957944*(f1R[4]+f0R[4])+0.6123724356957944*(f1L[4]+f0L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])+0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1]+f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1UpR[0]-0.5*alpha1UpR[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]))+0.6123724356957944*(f1R[6]+f1L[6]+f0R[6]+f0L[6])+0.3535533905932737*(f1R[5]+f0R[5])-0.3535533905932737*(f1L[5]+f0L[5])-0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])-0.3535533905932737*(f1R[3]+f0R[3])+0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])+0.3535533905932737*(f1R[1]+f0R[1])-0.3535533905932737*(f1L[1]+f0L[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[0]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f0R[7])-0.6123724356957944*(f1L[7]+f0L[7]+f1R[6]+f0R[6])+0.6123724356957944*(f1L[6]+f0L[6])-0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])+0.6123724356957944*(f1R[4]+f0R[4])-0.6123724356957944*(f1L[4]+f0L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])-0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1])+0.3535533905932737*(f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]+f1R[6]+f1L[6]+f0R[6]+f0L[6])-0.3535533905932737*(f1R[5]+f0R[5])+0.3535533905932737*(f1L[5]+f0L[5])+0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])-0.3535533905932737*(f1R[3]+f0R[3])+0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])-0.3535533905932737*(f1R[1]+f0R[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[1]+f0L[1]+f1L[0]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f0R[7]+f1R[6]+f0R[6])+0.6123724356957944*(f1L[7]+f0L[7]+f1L[6]+f0L[6])+0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])-0.6123724356957944*(f1R[4]+f0R[4])+0.6123724356957944*(f1L[4]+f0L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])+0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1]+f1R[0]+f1L[0]+f0R[0]+f0L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.3535533905932737*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.6123724356957944*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = 0.3535533905932737*(alpha1R[1]*fUp[3]+alpha1R[0]*fUp[2]); 
  incr[4] = -0.6123724356957944*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[5] = 0.3535533905932737*(alpha1R[0]*fUp[3]+alpha1R[1]*fUp[2]); 
  incr[6] = -0.6123724356957944*(alpha1R[1]*fUp[3]+alpha1R[0]*fUp[2]); 
  incr[7] = -0.6123724356957944*(alpha1R[0]*fUp[3]+alpha1R[1]*fUp[2]); 


  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSer_x_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // AparL: Apar in the left cell.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
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
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[5] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.732050807568877*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double BstarZdBmagEML[8]; 
  BstarZdBmagEML[0] = -1.732050807568877*(2.0*AparL[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(b_y[0]*AparL[1]-1.0*AparL[0]*b_y[1]))*rdx2R; 
  BstarZdBmagEML[1] = 1.732050807568877*((AparL[0]*b_y[1]-1.0*b_y[0]*AparL[1])*jacobTotInv[1]-2.0*jacobTotInv[0]*AparL[1]*b_y[1])*rdx2R; 

  double alpha0R[4]; 
  alpha0R[0] = -(0.25*(3.0*BstarZdBmagR[1]-1.732050807568877*BstarZdBmagR[0])*hamil0R[2]*rdvpar2R)/m_; 
  alpha0R[1] = -(0.25*hamil0R[2]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[2])*rdvpar2R)/m_; 

  double alpha1R[4]; 
  alpha1R[0] = -(0.125*(3.0*(BstarZdBmagEMR[1]+BstarZdBmagEML[1])-1.732050807568877*(BstarZdBmagEMR[0]+BstarZdBmagEML[0]))*hamil0R[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(6.0*BstarZdBmagR[1]+3.0*BstarZdBmagEMR[1]+3.0*BstarZdBmagEML[1]-3.464101615137754*BstarZdBmagR[0]-1.732050807568877*BstarZdBmagEMR[0]-1.732050807568877*BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  double incr[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0])-0.5*alpha0R[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.3535533905932737*f1R[6]+0.3535533905932737*f1L[6]-0.6123724356957944*(f1R[5]+f1L[5]+f1R[4]+f1L[4])+0.3535533905932737*(f1R[3]+f1R[2])-0.3535533905932737*(f1L[3]+f1L[2])+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*f1L[7]+0.3535533905932737*(f1R[6]+f1L[6])+0.6123724356957944*(f1R[5]+f1R[4])-0.6123724356957944*(f1L[5]+f1L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha0R[1]+alpha1R[0]+alpha0R[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.3535533905932737*f1R[6]-0.3535533905932737*f1L[6]-0.6123724356957944*(f1R[5]+f1L[5])+0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*(f1L[3]+f1R[2])+0.3535533905932737*f1L[2]+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*f1L[7]-0.3535533905932737*(f1R[6]+f1L[6])+0.6123724356957944*f1R[5]-0.6123724356957944*(f1L[5]+f1R[4])+0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])+0.3535533905932737*(f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1R[0]+alpha0R[0])-0.5*alpha0R[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.3535533905932737*f1R[6]-0.3535533905932737*f1L[6]+0.6123724356957944*(f1R[5]+f1L[5])-0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*(f1L[3]+f1R[2])-0.3535533905932737*f1L[2]+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*f1L[7]-0.3535533905932737*(f1R[6]+f1L[6])-0.6123724356957944*f1R[5]+0.6123724356957944*(f1L[5]+f1R[4])-0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.3535533905932737*(f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha0R[1]+alpha1R[0]+alpha0R[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.3535533905932737*f1R[6]+0.3535533905932737*f1L[6]+0.6123724356957944*(f1R[5]+f1L[5]+f1R[4]+f1L[4])-0.3535533905932737*(f1R[3]+f1R[2])+0.3535533905932737*(f1L[3]+f1L[2])+0.6123724356957944*(f1R[1]+f1L[1])-0.3535533905932737*f1R[0]+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*f1L[7]+0.3535533905932737*(f1R[6]+f1L[6])-0.6123724356957944*(f1R[5]+f1R[4])+0.6123724356957944*(f1L[5]+f1L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f1R[2]+f1L[2])-0.6123724356957944*f1R[1]+0.6123724356957944*f1L[1]+0.3535533905932737*(f1R[0]+f1L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alpha0R[1]*fUp[1]+(alpha1R[0]+alpha0R[0])*fUp[0]); 
  incr[1] = -0.6123724356957944*(alpha0R[1]*fUp[1]+(alpha1R[0]+alpha0R[0])*fUp[0]); 
  incr[2] = 0.3535533905932737*((alpha1R[0]+alpha0R[0])*fUp[1]+fUp[0]*alpha0R[1]); 
  incr[3] = 0.3535533905932737*(alpha0R[1]*fUp[3]+(alpha1R[0]+alpha0R[0])*fUp[2]); 
  incr[4] = -0.6123724356957944*((alpha1R[0]+alpha0R[0])*fUp[1]+fUp[0]*alpha0R[1]); 
  incr[5] = -0.6123724356957944*(alpha0R[1]*fUp[3]+(alpha1R[0]+alpha0R[0])*fUp[2]); 
  incr[6] = 0.3535533905932737*((alpha1R[0]+alpha0R[0])*fUp[3]+alpha0R[1]*fUp[2]); 
  incr[7] = -0.6123724356957944*((alpha1R[0]+alpha0R[0])*fUp[3]+alpha0R[1]*fUp[2]); 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += incr[7]*rdx2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])+0.3535533905932737*(f0R[3]+f0R[2])-0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*(f0R[5]+f0R[4])-0.6123724356957944*(f0L[5]+f0L[4])-0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]-0.6123724356957944*(f0R[5]+f0L[5])+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*(f0L[3]+f0R[2])+0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])+0.6123724356957944*f0R[5]-0.6123724356957944*(f0L[5]+f0R[4])+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])+0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.3535533905932737*f0R[6]-0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5])-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*(f0L[3]+f0R[2])-0.3535533905932737*f0L[2]+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*f0L[7]-0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*f0R[5]+0.6123724356957944*(f0L[5]+f0R[4])-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.3535533905932737*(f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*alpha1R[0]; 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.3535533905932737*f0R[6]+0.3535533905932737*f0L[6]+0.6123724356957944*(f0R[5]+f0L[5]+f0R[4]+f0L[4])-0.3535533905932737*(f0R[3]+f0R[2])+0.3535533905932737*(f0L[3]+f0L[2])+0.6123724356957944*(f0R[1]+f0L[1])-0.3535533905932737*f0R[0]+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*f0L[7]+0.3535533905932737*(f0R[6]+f0L[6])-0.6123724356957944*(f0R[5]+f0R[4])+0.6123724356957944*(f0L[5]+f0L[4])+0.3535533905932737*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.6123724356957944*f0R[1]+0.6123724356957944*f0L[1]+0.3535533905932737*(f0R[0]+f0L[0])); 

  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*alpha1R[0]*fUp[0]; 
  incr[1] = -0.6123724356957944*alpha1R[0]*fUp[0]; 
  incr[2] = 0.3535533905932737*alpha1R[0]*fUp[1]; 
  incr[3] = 0.3535533905932737*alpha1R[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alpha1R[0]*fUp[1]; 
  incr[5] = -0.6123724356957944*alpha1R[0]*fUp[2]; 
  incr[6] = 0.3535533905932737*alpha1R[0]*fUp[3]; 
  incr[7] = -0.6123724356957944*alpha1R[0]*fUp[3]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 
  outR[4] += incr[4]*rdx2R; 
  outR[5] += incr[5]*rdx2R; 
  outR[6] += incr[6]*rdx2R; 
  outR[7] += incr[7]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += -1.0*incr[3]*rdx2L; 
  outL[4] += incr[4]*rdx2L; 
  outL[5] += incr[5]*rdx2L; 
  outL[6] += -1.0*incr[6]*rdx2L; 
  outL[7] += incr[7]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSer_vpar_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
  // AparL: Apar in the left cell.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // emModL,emModR: .
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
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[5] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.732050807568877*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double alpha0R[4]; 
  alpha0R[0] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alpha0R[1] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 
  alpha0R[2] = (0.25*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*hamil0R[5]*rdx2R)/m_; 
  alpha0R[3] = (0.25*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*hamil0R[5]*rdx2R)/m_; 

  double alpha0UpR[4]; 
  alpha0UpR[0] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alpha0UpR[1] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 
  alpha0UpR[2] = (0.25*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*hamil0R[5]*rdx2R)/m_; 
  alpha0UpR[3] = (0.25*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*hamil0R[5]*rdx2R)/m_; 

  double alpha1R[4]; 
  alpha1R[0] = (0.25*(3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*(BstarZdBmagR[0]*hamil1R[1]+BstarZdBmagEMR[0]*hamil0R[1]))*rdx2R)/m_; 
  alpha1R[1] = (0.25*(3.0*hamil1R[1]*BstarZdBmagR[4]-1.732050807568877*(BstarZdBmagR[1]*hamil1R[1]+BstarZdBmagEMR[1]*hamil0R[1]))*rdx2R)/m_; 
  alpha1R[2] = -(0.4330127018922193*BstarZdBmagEMR[0]*hamil0R[5]*rdx2R)/m_; 
  alpha1R[3] = -(0.4330127018922193*BstarZdBmagEMR[1]*hamil0R[5]*rdx2R)/m_; 

  double alpha1UpR[4]; 
  alpha1UpR[0] = (0.25*((3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*(BstarZdBmagR[0]*hamil1R[1]+BstarZdBmagEMR[0]*hamil0R[1]))*rdx2R-5.656854249492382*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = (0.25*((3.0*hamil1R[1]*BstarZdBmagR[4]-1.732050807568877*(BstarZdBmagR[1]*hamil1R[1]+BstarZdBmagEMR[1]*hamil0R[1]))*rdx2R-5.656854249492382*dApardtPrev[1]*q_))/m_; 
  alpha1UpR[2] = -(0.4330127018922193*BstarZdBmagEMR[0]*hamil0R[5]*rdx2R)/m_; 
  alpha1UpR[3] = -(0.4330127018922193*BstarZdBmagEMR[1]*hamil0R[5]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(((3.0*hamil1R[1]+3.0*hamil0R[1])*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamil1R[1]+((-1.732050807568877*BstarZdBmagR[0])-1.732050807568877*BstarZdBmagEMR[0])*hamil0R[1])*rdx2R-5.656854249492382*dApardtPrev[0]*q_))/m_; 

  double incr[8]; 
  double incrEmMod[8]; 
  // linear + nonlinear terms (alpha0+alpha1)*f1 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*(alpha1UpR[3]+alpha0UpR[3])-0.5*(alpha1UpR[2]+alpha0UpR[2]+alpha1UpR[1]+alpha0UpR[1])+0.5*(alpha1UpR[0]+alpha0UpR[0]); 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7])-0.6123724356957944*(f1R[6]+f1L[6])-0.3535533905932737*f1R[5]+0.3535533905932737*f1L[5]-0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])+0.3535533905932737*f1R[1]-0.3535533905932737*(f1L[1]+f1R[0])+0.3535533905932737*f1L[0])*sgn(alphaOrdR)-0.6123724356957944*f1R[7]+0.6123724356957944*(f1L[7]+f1R[6])-0.6123724356957944*f1L[6]+0.3535533905932737*(f1R[5]+f1L[5])+0.6123724356957944*f1R[4]-0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]-0.3535533905932737*(f1R[1]+f1L[1])+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha0UpR[1]+alpha1UpR[0]+alpha0UpR[0])-0.5*(alpha1UpR[3]+alpha0UpR[3]+alpha1UpR[2]+alpha0UpR[2]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f1R[6]+f1L[6]))+0.3535533905932737*f1R[5]-0.3535533905932737*f1L[5]+0.6123724356957944*(f1R[4]+f1L[4])+0.3535533905932737*f1R[3]-0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])-0.3535533905932737*(f1R[1]+f1R[0])+0.3535533905932737*(f1L[1]+f1L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f1R[6])-0.6123724356957944*(f1L[7]+f1L[6])-0.3535533905932737*(f1R[5]+f1L[5])-0.6123724356957944*f1R[4]+0.6123724356957944*f1L[4]-0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]+0.3535533905932737*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 
  alphaOrdR = (-0.5*(alpha1UpR[3]+alpha0UpR[3]))+0.5*(alpha1UpR[2]+alpha0UpR[2])-0.5*(alpha1UpR[1]+alpha0UpR[1])+0.5*(alpha1UpR[0]+alpha0UpR[0]); 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]))+0.6123724356957944*(f1R[6]+f1L[6])+0.3535533905932737*f1R[5]-0.3535533905932737*f1L[5]-0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])+0.3535533905932737*f1R[1]-0.3535533905932737*(f1L[1]+f1R[0])+0.3535533905932737*f1L[0])*sgn(alphaOrdR)+0.6123724356957944*f1R[7]-0.6123724356957944*(f1L[7]+f1R[6])+0.6123724356957944*f1L[6]-0.3535533905932737*(f1R[5]+f1L[5])+0.6123724356957944*f1R[4]-0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]-0.3535533905932737*(f1R[1]+f1L[1])+0.3535533905932737*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[3]+alpha0UpR[3]+alpha1UpR[2]+alpha0UpR[2]+alpha1UpR[1]+alpha0UpR[1]+alpha1UpR[0]+alpha0UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f1R[6]+f1L[6])-0.3535533905932737*f1R[5]+0.3535533905932737*f1L[5]+0.6123724356957944*(f1R[4]+f1L[4])-0.3535533905932737*f1R[3]+0.3535533905932737*f1L[3]+0.6123724356957944*(f1R[2]+f1L[2])-0.3535533905932737*(f1R[1]+f1R[0])+0.3535533905932737*(f1L[1]+f1L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f1R[6])+0.6123724356957944*(f1L[7]+f1L[6])+0.3535533905932737*(f1R[5]+f1L[5])-0.6123724356957944*f1R[4]+0.6123724356957944*f1L[4]+0.3535533905932737*(f1R[3]+f1L[3])-0.6123724356957944*f1R[2]+0.6123724356957944*f1L[2]+0.3535533905932737*(f1R[1]+f1L[1]+f1R[0]+f1L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = 0.7071067811865475*fUp[2]; 
  incrEmMod[4] = -1.224744871391589*fUp[1]; 
  incrEmMod[5] = 0.7071067811865475*fUp[3]; 
  incrEmMod[6] = -1.224744871391589*fUp[2]; 
  incrEmMod[7] = -1.224744871391589*fUp[3]; 

  incr[0] = 0.3535533905932737*((alpha1R[3]+alpha0R[3])*fUp[3]+(alpha1R[2]+alpha0R[2])*fUp[2]+(alpha1R[1]+alpha0R[1])*fUp[1]+(alpha1R[0]+alpha0R[0])*fUp[0]); 
  incr[1] = 0.3535533905932737*((alpha1R[2]+alpha0R[2])*fUp[3]+fUp[2]*(alpha1R[3]+alpha0R[3])+(alpha1R[0]+alpha0R[0])*fUp[1]+fUp[0]*(alpha1R[1]+alpha0R[1])); 
  incr[2] = -0.6123724356957944*((alpha1R[3]+alpha0R[3])*fUp[3]+(alpha1R[2]+alpha0R[2])*fUp[2]+(alpha1R[1]+alpha0R[1])*fUp[1]+(alpha1R[0]+alpha0R[0])*fUp[0]); 
  incr[3] = 0.3535533905932737*((alpha1R[1]+alpha0R[1])*fUp[3]+fUp[1]*(alpha1R[3]+alpha0R[3])+(alpha1R[0]+alpha0R[0])*fUp[2]+fUp[0]*(alpha1R[2]+alpha0R[2])); 
  incr[4] = -0.6123724356957944*((alpha1R[2]+alpha0R[2])*fUp[3]+fUp[2]*(alpha1R[3]+alpha0R[3])+(alpha1R[0]+alpha0R[0])*fUp[1]+fUp[0]*(alpha1R[1]+alpha0R[1])); 
  incr[5] = 0.3535533905932737*((alpha1R[0]+alpha0R[0])*fUp[3]+fUp[0]*(alpha1R[3]+alpha0R[3])+(alpha1R[1]+alpha0R[1])*fUp[2]+fUp[1]*(alpha1R[2]+alpha0R[2])); 
  incr[6] = -0.6123724356957944*((alpha1R[1]+alpha0R[1])*fUp[3]+fUp[1]*(alpha1R[3]+alpha0R[3])+(alpha1R[0]+alpha0R[0])*fUp[2]+fUp[0]*(alpha1R[2]+alpha0R[2])); 
  incr[7] = -0.6123724356957944*((alpha1R[0]+alpha0R[0])*fUp[3]+fUp[0]*(alpha1R[3]+alpha0R[3])+(alpha1R[1]+alpha0R[1])*fUp[2]+fUp[1]*(alpha1R[2]+alpha0R[2])); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModR[4] += incrEmMod[4]*rdvpar2R; 
  emModR[5] += incrEmMod[5]*rdvpar2R; 
  emModR[6] += incrEmMod[6]*rdvpar2R; 
  emModR[7] += incrEmMod[7]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += -1.0*incrEmMod[3]*rdvpar2L; 
  emModL[4] += incrEmMod[4]*rdvpar2L; 
  emModL[5] += -1.0*incrEmMod[5]*rdvpar2L; 
  emModL[6] += incrEmMod[6]*rdvpar2L; 
  emModL[7] += incrEmMod[7]*rdvpar2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.5*alpha1UpR[3]-0.5*(alpha1UpR[2]+alpha1UpR[1])+0.5*alpha1UpR[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7])-0.6123724356957944*(f0R[6]+f0L[6])-0.3535533905932737*f0R[5]+0.3535533905932737*f0L[5]-0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])+0.3535533905932737*f0R[1]-0.3535533905932737*(f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)-0.6123724356957944*f0R[7]+0.6123724356957944*(f0L[7]+f0R[6])-0.6123724356957944*f0L[6]+0.3535533905932737*(f0R[5]+f0L[5])+0.6123724356957944*f0R[4]-0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]-0.3535533905932737*(f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0])-0.5*(alpha1UpR[3]+alpha1UpR[2]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]+f0R[6]+f0L[6]))+0.3535533905932737*f0R[5]-0.3535533905932737*f0L[5]+0.6123724356957944*(f0R[4]+f0L[4])+0.3535533905932737*f0R[3]-0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])-0.3535533905932737*(f0R[1]+f0R[0])+0.3535533905932737*(f0L[1]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f0R[7]+f0R[6])-0.6123724356957944*(f0L[7]+f0L[6])-0.3535533905932737*(f0R[5]+f0L[5])-0.6123724356957944*f0R[4]+0.6123724356957944*f0L[4]-0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]+0.3535533905932737*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 
  alphaOrdR = (-0.5*alpha1UpR[3])+0.5*alpha1UpR[2]-0.5*alpha1UpR[1]+0.5*alpha1UpR[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f0R[7]+f0L[7]))+0.6123724356957944*(f0R[6]+f0L[6])+0.3535533905932737*f0R[5]-0.3535533905932737*f0L[5]-0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])+0.3535533905932737*f0R[1]-0.3535533905932737*(f0L[1]+f0R[0])+0.3535533905932737*f0L[0])*sgn(alphaOrdR)+0.6123724356957944*f0R[7]-0.6123724356957944*(f0L[7]+f0R[6])+0.6123724356957944*f0L[6]-0.3535533905932737*(f0R[5]+f0L[5])+0.6123724356957944*f0R[4]-0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]-0.3535533905932737*(f0R[1]+f0L[1])+0.3535533905932737*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[3]+alpha1UpR[2]+alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f0R[7]+f0L[7]+f0R[6]+f0L[6])-0.3535533905932737*f0R[5]+0.3535533905932737*f0L[5]+0.6123724356957944*(f0R[4]+f0L[4])-0.3535533905932737*f0R[3]+0.3535533905932737*f0L[3]+0.6123724356957944*(f0R[2]+f0L[2])-0.3535533905932737*(f0R[1]+f0R[0])+0.3535533905932737*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f0R[7]+f0R[6])+0.6123724356957944*(f0L[7]+f0L[6])+0.3535533905932737*(f0R[5]+f0L[5])-0.6123724356957944*f0R[4]+0.6123724356957944*f0L[4]+0.3535533905932737*(f0R[3]+f0L[3])-0.6123724356957944*f0R[2]+0.6123724356957944*f0L[2]+0.3535533905932737*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = 0.7071067811865475*fUp[2]; 
  incrEmMod[4] = -1.224744871391589*fUp[1]; 
  incrEmMod[5] = 0.7071067811865475*fUp[3]; 
  incrEmMod[6] = -1.224744871391589*fUp[2]; 
  incrEmMod[7] = -1.224744871391589*fUp[3]; 

  incr[0] = 0.3535533905932737*(alpha1R[3]*fUp[3]+alpha1R[2]*fUp[2]+alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.3535533905932737*(alpha1R[2]*fUp[3]+fUp[2]*alpha1R[3]+alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.6123724356957944*(alpha1R[3]*fUp[3]+alpha1R[2]*fUp[2]+alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = 0.3535533905932737*(alpha1R[1]*fUp[3]+fUp[1]*alpha1R[3]+alpha1R[0]*fUp[2]+fUp[0]*alpha1R[2]); 
  incr[4] = -0.6123724356957944*(alpha1R[2]*fUp[3]+fUp[2]*alpha1R[3]+alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[5] = 0.3535533905932737*(alpha1R[0]*fUp[3]+fUp[0]*alpha1R[3]+alpha1R[1]*fUp[2]+fUp[1]*alpha1R[2]); 
  incr[6] = -0.6123724356957944*(alpha1R[1]*fUp[3]+fUp[1]*alpha1R[3]+alpha1R[0]*fUp[2]+fUp[0]*alpha1R[2]); 
  incr[7] = -0.6123724356957944*(alpha1R[0]*fUp[3]+fUp[0]*alpha1R[3]+alpha1R[1]*fUp[2]+fUp[1]*alpha1R[2]); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModR[4] += incrEmMod[4]*rdvpar2R; 
  emModR[5] += incrEmMod[5]*rdvpar2R; 
  emModR[6] += incrEmMod[6]*rdvpar2R; 
  emModR[7] += incrEmMod[7]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += -1.0*incrEmMod[3]*rdvpar2L; 
  emModL[4] += incrEmMod[4]*rdvpar2L; 
  emModL[5] += -1.0*incrEmMod[5]*rdvpar2L; 
  emModL[6] += incrEmMod[6]*rdvpar2L; 
  emModL[7] += incrEmMod[7]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmDeltaFGyrokineticGenGeoSurf1x2vSerStep2_vpar_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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
  // Apar: parallel component of magnetic vector potential.
  // dApardt: time derivative of Apar.
  // dApardtPrev: previous dApardt.
  // fL,fR: Distribution function in left and right cells.
  // outL/outR: Output increment in left and right cells.
  // emModL,emModR: .

  double wxL = wL[0];
  double wxR = wR[0];
  double rdx2L = 2.0/dxvL[0];
  double rdx2R = 2.0/dxvR[0];
  double wvparL = wL[1];
  double wvparR = wR[1];
  double rdvpar2L = 2.0/dxvL[1];
  double rdvpar2R = 2.0/dxvR[1];
  double wmuL = wL[2];
  double wmuR = wR[2];
  double rdmu2L = 2.0/dxvL[2];
  double rdmu2R = 2.0/dxvR[2];

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;
  double wmuSqL = wL[2]*wL[2];
  double wmuSqR = wR[2]*wR[2];
  double rdmu2SqL = rdmu2L*rdmu2L;
  double rdmu2SqR = rdmu2R*rdmu2R;

  double hamil0R[8]; 
  hamil0R[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*bmag[0]*wmuR)+2.0*m_))/rdvpar2SqR; 
  hamil0R[1] = 2.0*bmag[1]*wmuR; 
  hamil0R[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 
  hamil0R[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamil0R[5] = (1.154700538379252*bmag[1])/rdmu2R; 

  double hamil1R[8]; 
  hamil1R[0] = 2.0*phi[0]*q_; 
  hamil1R[1] = 2.0*phi[1]*q_; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[8]; 
  BstarZdBmagEMR[0] = 1.732050807568877*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.732050807568877*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double alpha0R[4]; 

  double alpha0UpR[4]; 
  alpha0UpR[0] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alpha0UpR[1] = (0.25*hamil0R[1]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 
  alpha0UpR[2] = (0.25*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*hamil0R[5]*rdx2R)/m_; 
  alpha0UpR[3] = (0.25*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*hamil0R[5]*rdx2R)/m_; 

  double alpha1R[4]; 
  alpha1R[0] = -(1.414213562373095*dApardt[0]*q_)/m_; 
  alpha1R[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 

  double alpha1UpR[4]; 
  alpha1UpR[0] = (0.25*((3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*(BstarZdBmagR[0]*hamil1R[1]+BstarZdBmagEMR[0]*hamil0R[1]))*rdx2R-5.656854249492382*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = (0.25*((3.0*hamil1R[1]*BstarZdBmagR[4]-1.732050807568877*(BstarZdBmagR[1]*hamil1R[1]+BstarZdBmagEMR[1]*hamil0R[1]))*rdx2R-5.656854249492382*dApardtPrev[1]*q_))/m_; 
  alpha1UpR[2] = -(0.4330127018922193*BstarZdBmagEMR[0]*hamil0R[5]*rdx2R)/m_; 
  alpha1UpR[3] = -(0.4330127018922193*BstarZdBmagEMR[1]*hamil0R[5]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(((3.0*hamil1R[1]+3.0*hamil0R[1])*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamil1R[1]+((-1.732050807568877*BstarZdBmagR[0])-1.732050807568877*BstarZdBmagEMR[0])*hamil0R[1])*rdx2R-5.656854249492382*dApardtPrev[0]*q_))/m_; 

  double incr[8]; 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alpha1UpR[3]-0.5*(alpha1UpR[2]+alpha1UpR[1])+0.5*alpha1UpR[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7])-0.6123724356957944*(f1R[6]+f1L[6]+f0R[6]+f0L[6])-0.3535533905932737*(f1R[5]+f0R[5])+0.3535533905932737*(f1L[5]+f0L[5])-0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])+0.3535533905932737*(f1R[3]+f0R[3])-0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])+0.3535533905932737*(f1R[1]+f0R[1])-0.3535533905932737*(f1L[1]+f0L[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[0]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f0R[7])+0.6123724356957944*(f1L[7]+f0L[7]+f1R[6]+f0R[6])-0.6123724356957944*(f1L[6]+f0L[6])+0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])+0.6123724356957944*(f1R[4]+f0R[4])-0.6123724356957944*(f1L[4]+f0L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])-0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1])+0.3535533905932737*(f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[1]+alpha1UpR[0])-0.5*(alpha1UpR[3]+alpha1UpR[2]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]+f1R[6]+f1L[6]+f0R[6]+f0L[6]))+0.3535533905932737*(f1R[5]+f0R[5])-0.3535533905932737*(f1L[5]+f0L[5])+0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])+0.3535533905932737*(f1R[3]+f0R[3])-0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])-0.3535533905932737*(f1R[1]+f0R[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[1]+f0L[1]+f1L[0]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f0R[7]+f1R[6]+f0R[6])-0.6123724356957944*(f1L[7]+f0L[7]+f1L[6]+f0L[6])-0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])-0.6123724356957944*(f1R[4]+f0R[4])+0.6123724356957944*(f1L[4]+f0L[4])-0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])+0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1]+f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = (-0.5*alpha1UpR[3])+0.5*alpha1UpR[2]-0.5*alpha1UpR[1]+0.5*alpha1UpR[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]))+0.6123724356957944*(f1R[6]+f1L[6]+f0R[6]+f0L[6])+0.3535533905932737*(f1R[5]+f0R[5])-0.3535533905932737*(f1L[5]+f0L[5])-0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])-0.3535533905932737*(f1R[3]+f0R[3])+0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])+0.3535533905932737*(f1R[1]+f0R[1])-0.3535533905932737*(f1L[1]+f0L[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[0]+f0L[0]))*sgn(alphaOrdR)+0.6123724356957944*(f1R[7]+f0R[7])-0.6123724356957944*(f1L[7]+f0L[7]+f1R[6]+f0R[6])+0.6123724356957944*(f1L[6]+f0L[6])-0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])+0.6123724356957944*(f1R[4]+f0R[4])-0.6123724356957944*(f1L[4]+f0L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])-0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1])+0.3535533905932737*(f1R[0]+f1L[0]+f0R[0]+f0L[0])); 
  alphaOrdR = 0.5*(alpha1UpR[3]+alpha1UpR[2]+alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(f1R[7]+f1L[7]+f0R[7]+f0L[7]+f1R[6]+f1L[6]+f0R[6]+f0L[6])-0.3535533905932737*(f1R[5]+f0R[5])+0.3535533905932737*(f1L[5]+f0L[5])+0.6123724356957944*(f1R[4]+f1L[4]+f0R[4]+f0L[4])-0.3535533905932737*(f1R[3]+f0R[3])+0.3535533905932737*(f1L[3]+f0L[3])+0.6123724356957944*(f1R[2]+f1L[2]+f0R[2]+f0L[2])-0.3535533905932737*(f1R[1]+f0R[1]+f1R[0]+f0R[0])+0.3535533905932737*(f1L[1]+f0L[1]+f1L[0]+f0L[0]))*sgn(alphaOrdR)-0.6123724356957944*(f1R[7]+f0R[7]+f1R[6]+f0R[6])+0.6123724356957944*(f1L[7]+f0L[7]+f1L[6]+f0L[6])+0.3535533905932737*(f1R[5]+f1L[5]+f0R[5]+f0L[5])-0.6123724356957944*(f1R[4]+f0R[4])+0.6123724356957944*(f1L[4]+f0L[4])+0.3535533905932737*(f1R[3]+f1L[3]+f0R[3]+f0L[3])-0.6123724356957944*(f1R[2]+f0R[2])+0.6123724356957944*(f1L[2]+f0L[2])+0.3535533905932737*(f1R[1]+f1L[1]+f0R[1]+f0L[1]+f1R[0]+f1L[0]+f0R[0]+f0L[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.3535533905932737*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.6123724356957944*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = 0.3535533905932737*(alpha1R[1]*fUp[3]+alpha1R[0]*fUp[2]); 
  incr[4] = -0.6123724356957944*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[5] = 0.3535533905932737*(alpha1R[0]*fUp[3]+alpha1R[1]*fUp[2]); 
  incr[6] = -0.6123724356957944*(alpha1R[1]*fUp[3]+alpha1R[0]*fUp[2]); 
  incr[7] = -0.6123724356957944*(alpha1R[0]*fUp[3]+alpha1R[1]*fUp[2]); 


  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 
  outR[4] += incr[4]*rdvpar2R; 
  outR[5] += incr[5]*rdvpar2R; 
  outR[6] += incr[6]*rdvpar2R; 
  outR[7] += incr[7]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += -1.0*incr[3]*rdvpar2L; 
  outL[4] += incr[4]*rdvpar2L; 
  outL[5] += -1.0*incr[5]*rdvpar2L; 
  outL[6] += incr[6]*rdvpar2L; 
  outL[7] += incr[7]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
