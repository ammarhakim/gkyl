#include <DeltaFGyrokineticModDecl.h>
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSer_x_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double BstarZdBmagEML[4]; 
  BstarZdBmagEML[0] = -1.224744871391589*b_y[0]*jacobTotInv[0]*AparL[1]*rdx2R; 

  double alpha0R[2]; 
  alpha0R[0] = (0.6123724356957944*BstarZdBmagR[0]*hamil0R[2]*rdvpar2R)/m_; 

  double alpha1R[2]; 
  alpha1R[0] = (0.3061862178478973*(BstarZdBmagEMR[0]+BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(3.464101615137754*BstarZdBmagR[0]+1.732050807568877*BstarZdBmagEMR[0]+1.732050807568877*BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  double incr[4]; 
  // linear term alpha0*f1 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha0R[0]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f1R[3]+f1L[3]))+0.5*f1R[2]-0.5*f1L[2]+0.8660254037844386*(f1R[1]+f1L[1])-0.5*f1R[0]+0.5*f1L[0])*sgn(alphaOrdR)+0.8660254037844386*f1R[3]-0.8660254037844386*f1L[3]-0.5*(f1R[2]+f1L[2])-0.8660254037844386*f1R[1]+0.8660254037844386*f1L[1]+0.5*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.7071067811865475*alpha0R[0]; 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f1R[3]+f1L[3])-0.5*f1R[2]+0.5*f1L[2]+0.8660254037844386*(f1R[1]+f1L[1])-0.5*f1R[0]+0.5*f1L[0])*sgn(alphaOrdR)-0.8660254037844386*f1R[3]+0.8660254037844386*f1L[3]+0.5*(f1R[2]+f1L[2])-0.8660254037844386*f1R[1]+0.8660254037844386*f1L[1]+0.5*(f1R[0]+f1L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*alpha0R[0]*fUp[0]; 
  incr[1] = -0.8660254037844386*alpha0R[0]*fUp[0]; 
  incr[2] = 0.5*alpha0R[0]*fUp[1]; 
  incr[3] = -0.8660254037844386*alpha0R[0]*fUp[1]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += incr[3]*rdx2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.7071067811865475*alpha1R[0]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.5*f0R[2]-0.5*f0L[2]+0.8660254037844386*(f0R[1]+f0L[1])-0.5*f0R[0]+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*f0L[3]-0.5*(f0R[2]+f0L[2])-0.8660254037844386*f0R[1]+0.8660254037844386*f0L[1]+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*alpha1R[0]; 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3])-0.5*f0R[2]+0.5*f0L[2]+0.8660254037844386*(f0R[1]+f0L[1])-0.5*f0R[0]+0.5*f0L[0])*sgn(alphaOrdR)-0.8660254037844386*f0R[3]+0.8660254037844386*f0L[3]+0.5*(f0R[2]+f0L[2])-0.8660254037844386*f0R[1]+0.8660254037844386*f0L[1]+0.5*(f0R[0]+f0L[0])); 

  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*alpha1R[0]*fUp[0]; 
  incr[1] = -0.8660254037844386*alpha1R[0]*fUp[0]; 
  incr[2] = 0.5*alpha1R[0]*fUp[1]; 
  incr[3] = -0.8660254037844386*alpha1R[0]*fUp[1]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += incr[3]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSer_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double alpha0R[2]; 

  double alpha0UpR[2]; 

  double alpha1R[2]; 
  alpha1R[0] = -(0.6123724356957944*BstarZdBmagR[0]*hamil1R[1]*rdx2R)/m_; 

  double alpha1UpR[2]; 
  alpha1UpR[0] = -(0.3535533905932737*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = -(1.0*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.125*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double incrEmMod[4]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha1UpR[0]-0.7071067811865475*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.8660254037844386*(f0R[2]+f0L[2])+0.5*f0R[1]-0.5*(f0L[1]+f0R[0])+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*(f0L[3]+f0R[2])+0.8660254037844386*f0L[2]-0.5*(f0R[1]+f0L[1])+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.5*(f0R[1]+f0R[0])+0.5*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.8660254037844386*(f0R[3]+f0R[2])+0.8660254037844386*(f0L[3]+f0L[2])+0.5*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = -1.224744871391589*fUp[1]; 

  incr[0] = 0.5*alpha1R[0]*fUp[0]; 
  incr[1] = 0.5*alpha1R[0]*fUp[1]; 
  incr[2] = -0.8660254037844386*alpha1R[0]*fUp[0]; 
  incr[3] = -0.8660254037844386*alpha1R[0]*fUp[1]; 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += incrEmMod[3]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSerStep2_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double alpha0R[2]; 

  double alpha0UpR[2]; 

  double alpha1R[2]; 
  alpha1R[0] = -(1.0*dApardt[0]*q_)/m_; 
  alpha1R[1] = -(1.0*dApardt[1]*q_)/m_; 

  double alpha1UpR[2]; 
  alpha1UpR[0] = -(0.3535533905932737*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = -(1.0*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.125*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha1UpR[0]-0.7071067811865475*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.8660254037844386*(f0R[2]+f0L[2])+0.5*f0R[1]-0.5*(f0L[1]+f0R[0])+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*(f0L[3]+f0R[2])+0.8660254037844386*f0L[2]-0.5*(f0R[1]+f0L[1])+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.5*(f0R[1]+f0R[0])+0.5*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.8660254037844386*(f0R[3]+f0R[2])+0.8660254037844386*(f0L[3]+f0L[2])+0.5*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.5*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.8660254037844386*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 


  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSer_x_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmagR[1] = (1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_)/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.224744871391589*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double BstarZdBmagEML[4]; 
  BstarZdBmagEML[0] = -1.224744871391589*(2.0*AparL[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(b_y[0]*AparL[1]-1.0*AparL[0]*b_y[1]))*rdx2R; 
  BstarZdBmagEML[1] = 1.224744871391589*((AparL[0]*b_y[1]-1.0*b_y[0]*AparL[1])*jacobTotInv[1]-2.0*jacobTotInv[0]*AparL[1]*b_y[1])*rdx2R; 

  double alpha0R[2]; 
  alpha0R[0] = -(0.3535533905932737*(3.0*BstarZdBmagR[1]-1.732050807568877*BstarZdBmagR[0])*hamil0R[2]*rdvpar2R)/m_; 
  alpha0R[1] = -(0.3535533905932737*hamil0R[2]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[2])*rdvpar2R)/m_; 

  double alpha1R[2]; 
  alpha1R[0] = -(0.125*(4.242640687119286*(BstarZdBmagEMR[1]+BstarZdBmagEML[1])-2.449489742783178*(BstarZdBmagEMR[0]+BstarZdBmagEML[0]))*hamil0R[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(6.0*BstarZdBmagR[1]+3.0*BstarZdBmagEMR[1]+3.0*BstarZdBmagEML[1]-3.464101615137754*BstarZdBmagR[0]-1.732050807568877*BstarZdBmagEMR[0]-1.732050807568877*BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  double incr[4]; 
  // linear term alpha0*f1 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha0R[0]-0.7071067811865475*alpha0R[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f1R[3]+f1L[3]))+0.5*f1R[2]-0.5*f1L[2]+0.8660254037844386*(f1R[1]+f1L[1])-0.5*f1R[0]+0.5*f1L[0])*sgn(alphaOrdR)+0.8660254037844386*f1R[3]-0.8660254037844386*f1L[3]-0.5*(f1R[2]+f1L[2])-0.8660254037844386*f1R[1]+0.8660254037844386*f1L[1]+0.5*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.7071067811865475*(alpha0R[1]+alpha0R[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f1R[3]+f1L[3])-0.5*f1R[2]+0.5*f1L[2]+0.8660254037844386*(f1R[1]+f1L[1])-0.5*f1R[0]+0.5*f1L[0])*sgn(alphaOrdR)-0.8660254037844386*f1R[3]+0.8660254037844386*f1L[3]+0.5*(f1R[2]+f1L[2])-0.8660254037844386*f1R[1]+0.8660254037844386*f1L[1]+0.5*(f1R[0]+f1L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[1] = -0.8660254037844386*(alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[2] = 0.5*(alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1]); 
  incr[3] = -0.8660254037844386*(alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1]); 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += incr[3]*rdx2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.7071067811865475*alpha1R[0]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.5*f0R[2]-0.5*f0L[2]+0.8660254037844386*(f0R[1]+f0L[1])-0.5*f0R[0]+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*f0L[3]-0.5*(f0R[2]+f0L[2])-0.8660254037844386*f0R[1]+0.8660254037844386*f0L[1]+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*alpha1R[0]; 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3])-0.5*f0R[2]+0.5*f0L[2]+0.8660254037844386*(f0R[1]+f0L[1])-0.5*f0R[0]+0.5*f0L[0])*sgn(alphaOrdR)-0.8660254037844386*f0R[3]+0.8660254037844386*f0L[3]+0.5*(f0R[2]+f0L[2])-0.8660254037844386*f0R[1]+0.8660254037844386*f0L[1]+0.5*(f0R[0]+f0L[0])); 

  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*alpha1R[0]*fUp[0]; 
  incr[1] = -0.8660254037844386*alpha1R[0]*fUp[0]; 
  incr[2] = 0.5*alpha1R[0]*fUp[1]; 
  incr[3] = -0.8660254037844386*alpha1R[0]*fUp[1]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += incr[3]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSer_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmagR[1] = (1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_)/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.224744871391589*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double alpha0R[2]; 

  double alpha0UpR[2]; 

  double alpha1R[2]; 
  alpha1R[0] = (0.3535533905932737*hamil1R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alpha1R[1] = (0.3535533905932737*hamil1R[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 

  double alpha1UpR[2]; 
  alpha1UpR[0] = (0.3535533905932737*(hamil1R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = (0.3535533905932737*(hamil1R[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R-2.828427124746191*dApardtPrev[1]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.125*((3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamil1R[1])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double incrEmMod[4]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha1UpR[0]-0.7071067811865475*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.8660254037844386*(f0R[2]+f0L[2])+0.5*f0R[1]-0.5*(f0L[1]+f0R[0])+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*(f0L[3]+f0R[2])+0.8660254037844386*f0L[2]-0.5*(f0R[1]+f0L[1])+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.5*(f0R[1]+f0R[0])+0.5*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.8660254037844386*(f0R[3]+f0R[2])+0.8660254037844386*(f0L[3]+f0L[2])+0.5*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = -1.224744871391589*fUp[1]; 

  incr[0] = 0.5*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.5*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.8660254037844386*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += incrEmMod[3]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSerStep2_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmagR[1] = (1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_)/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.224744871391589*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double alpha0R[2]; 

  double alpha0UpR[2]; 

  double alpha1R[2]; 
  alpha1R[0] = -(1.0*dApardt[0]*q_)/m_; 
  alpha1R[1] = -(1.0*dApardt[1]*q_)/m_; 

  double alpha1UpR[2]; 
  alpha1UpR[0] = (0.3535533905932737*(hamil1R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = (0.3535533905932737*(hamil1R[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R-2.828427124746191*dApardtPrev[1]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.125*((3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamil1R[1])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha1UpR[0]-0.7071067811865475*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.8660254037844386*(f0R[2]+f0L[2])+0.5*f0R[1]-0.5*(f0L[1]+f0R[0])+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*(f0L[3]+f0R[2])+0.8660254037844386*f0L[2]-0.5*(f0R[1]+f0L[1])+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.5*(f0R[1]+f0R[0])+0.5*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.8660254037844386*(f0R[3]+f0R[2])+0.8660254037844386*(f0L[3]+f0L[2])+0.5*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.5*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.8660254037844386*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 


  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSer_x_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double BstarZdBmagEML[4]; 
  BstarZdBmagEML[0] = -1.224744871391589*b_y[0]*jacobTotInv[0]*AparL[1]*rdx2R; 

  double alpha0R[2]; 
  alpha0R[0] = (0.6123724356957944*BstarZdBmagR[0]*hamil0R[2]*rdvpar2R)/m_; 

  double alpha1R[2]; 
  alpha1R[0] = (0.3061862178478973*(BstarZdBmagEMR[0]+BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(3.464101615137754*BstarZdBmagR[0]+1.732050807568877*BstarZdBmagEMR[0]+1.732050807568877*BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  double incr[4]; 
  // linear term alpha0*f1 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha0R[0]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f1R[3]+f1L[3]))+0.5*f1R[2]-0.5*f1L[2]+0.8660254037844386*(f1R[1]+f1L[1])-0.5*f1R[0]+0.5*f1L[0])*sgn(alphaOrdR)+0.8660254037844386*f1R[3]-0.8660254037844386*f1L[3]-0.5*(f1R[2]+f1L[2])-0.8660254037844386*f1R[1]+0.8660254037844386*f1L[1]+0.5*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.7071067811865475*alpha0R[0]; 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f1R[3]+f1L[3])-0.5*f1R[2]+0.5*f1L[2]+0.8660254037844386*(f1R[1]+f1L[1])-0.5*f1R[0]+0.5*f1L[0])*sgn(alphaOrdR)-0.8660254037844386*f1R[3]+0.8660254037844386*f1L[3]+0.5*(f1R[2]+f1L[2])-0.8660254037844386*f1R[1]+0.8660254037844386*f1L[1]+0.5*(f1R[0]+f1L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*alpha0R[0]*fUp[0]; 
  incr[1] = -0.8660254037844386*alpha0R[0]*fUp[0]; 
  incr[2] = 0.5*alpha0R[0]*fUp[1]; 
  incr[3] = -0.8660254037844386*alpha0R[0]*fUp[1]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += incr[3]*rdx2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.7071067811865475*alpha1R[0]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.5*f0R[2]-0.5*f0L[2]+0.8660254037844386*(f0R[1]+f0L[1])-0.5*f0R[0]+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*f0L[3]-0.5*(f0R[2]+f0L[2])-0.8660254037844386*f0R[1]+0.8660254037844386*f0L[1]+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*alpha1R[0]; 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3])-0.5*f0R[2]+0.5*f0L[2]+0.8660254037844386*(f0R[1]+f0L[1])-0.5*f0R[0]+0.5*f0L[0])*sgn(alphaOrdR)-0.8660254037844386*f0R[3]+0.8660254037844386*f0L[3]+0.5*(f0R[2]+f0L[2])-0.8660254037844386*f0R[1]+0.8660254037844386*f0L[1]+0.5*(f0R[0]+f0L[0])); 

  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*alpha1R[0]*fUp[0]; 
  incr[1] = -0.8660254037844386*alpha1R[0]*fUp[0]; 
  incr[2] = 0.5*alpha1R[0]*fUp[1]; 
  incr[3] = -0.8660254037844386*alpha1R[0]*fUp[1]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += incr[3]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSer_vpar_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double alpha0R[2]; 

  double alpha0UpR[2]; 

  double alpha1R[2]; 
  alpha1R[0] = -(0.6123724356957944*BstarZdBmagR[0]*hamil1R[1]*rdx2R)/m_; 

  double alpha1UpR[2]; 
  alpha1UpR[0] = -(0.3535533905932737*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = -(1.0*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.125*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double incrEmMod[4]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha1UpR[0]-0.7071067811865475*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.8660254037844386*(f0R[2]+f0L[2])+0.5*f0R[1]-0.5*(f0L[1]+f0R[0])+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*(f0L[3]+f0R[2])+0.8660254037844386*f0L[2]-0.5*(f0R[1]+f0L[1])+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.5*(f0R[1]+f0R[0])+0.5*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.8660254037844386*(f0R[3]+f0R[2])+0.8660254037844386*(f0L[3]+f0L[2])+0.5*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = -1.224744871391589*fUp[1]; 

  incr[0] = 0.5*alpha1R[0]*fUp[0]; 
  incr[1] = 0.5*alpha1R[0]*fUp[1]; 
  incr[2] = -0.8660254037844386*alpha1R[0]*fUp[0]; 
  incr[3] = -0.8660254037844386*alpha1R[0]*fUp[1]; 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += incrEmMod[3]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSerStep2_vpar_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = cmag[0]*jacobTotInv[0]; 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*b_y[0]*jacobTotInv[0]*Apar[1]*rdx2R; 

  double alpha0R[2]; 

  double alpha0UpR[2]; 

  double alpha1R[2]; 
  alpha1R[0] = -(1.0*dApardt[0]*q_)/m_; 
  alpha1R[1] = -(1.0*dApardt[1]*q_)/m_; 

  double alpha1UpR[2]; 
  alpha1UpR[0] = -(0.3535533905932737*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = -(1.0*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.125*(1.732050807568877*BstarZdBmagR[0]*hamil1R[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha1UpR[0]-0.7071067811865475*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.8660254037844386*(f0R[2]+f0L[2])+0.5*f0R[1]-0.5*(f0L[1]+f0R[0])+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*(f0L[3]+f0R[2])+0.8660254037844386*f0L[2]-0.5*(f0R[1]+f0L[1])+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.5*(f0R[1]+f0R[0])+0.5*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.8660254037844386*(f0R[3]+f0R[2])+0.8660254037844386*(f0L[3]+f0L[2])+0.5*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.5*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.8660254037844386*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 


  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSer_x_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmagR[1] = (1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_)/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.224744871391589*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double BstarZdBmagEML[4]; 
  BstarZdBmagEML[0] = -1.224744871391589*(2.0*AparL[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(b_y[0]*AparL[1]-1.0*AparL[0]*b_y[1]))*rdx2R; 
  BstarZdBmagEML[1] = 1.224744871391589*((AparL[0]*b_y[1]-1.0*b_y[0]*AparL[1])*jacobTotInv[1]-2.0*jacobTotInv[0]*AparL[1]*b_y[1])*rdx2R; 

  double alpha0R[2]; 
  alpha0R[0] = -(0.3535533905932737*(3.0*BstarZdBmagR[1]-1.732050807568877*BstarZdBmagR[0])*hamil0R[2]*rdvpar2R)/m_; 
  alpha0R[1] = -(0.3535533905932737*hamil0R[2]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[2])*rdvpar2R)/m_; 

  double alpha1R[2]; 
  alpha1R[0] = -(0.125*(4.242640687119286*(BstarZdBmagEMR[1]+BstarZdBmagEML[1])-2.449489742783178*(BstarZdBmagEMR[0]+BstarZdBmagEML[0]))*hamil0R[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(6.0*BstarZdBmagR[1]+3.0*BstarZdBmagEMR[1]+3.0*BstarZdBmagEML[1]-3.464101615137754*BstarZdBmagR[0]-1.732050807568877*BstarZdBmagEMR[0]-1.732050807568877*BstarZdBmagEML[0])*hamil0R[2]*rdvpar2R)/m_; 

  double incr[4]; 
  // linear term alpha0*f1 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha0R[0]-0.7071067811865475*alpha0R[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f1R[3]+f1L[3]))+0.5*f1R[2]-0.5*f1L[2]+0.8660254037844386*(f1R[1]+f1L[1])-0.5*f1R[0]+0.5*f1L[0])*sgn(alphaOrdR)+0.8660254037844386*f1R[3]-0.8660254037844386*f1L[3]-0.5*(f1R[2]+f1L[2])-0.8660254037844386*f1R[1]+0.8660254037844386*f1L[1]+0.5*(f1R[0]+f1L[0])); 
  alphaOrdR = 0.7071067811865475*(alpha0R[1]+alpha0R[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f1R[3]+f1L[3])-0.5*f1R[2]+0.5*f1L[2]+0.8660254037844386*(f1R[1]+f1L[1])-0.5*f1R[0]+0.5*f1L[0])*sgn(alphaOrdR)-0.8660254037844386*f1R[3]+0.8660254037844386*f1L[3]+0.5*(f1R[2]+f1L[2])-0.8660254037844386*f1R[1]+0.8660254037844386*f1L[1]+0.5*(f1R[0]+f1L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[1] = -0.8660254037844386*(alpha0R[1]*fUp[1]+alpha0R[0]*fUp[0]); 
  incr[2] = 0.5*(alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1]); 
  incr[3] = -0.8660254037844386*(alpha0R[0]*fUp[1]+fUp[0]*alpha0R[1]); 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += incr[3]*rdx2L; 

  // linear term alpha1*f0 
  alphaOrdR = 0.7071067811865475*alpha1R[0]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.5*f0R[2]-0.5*f0L[2]+0.8660254037844386*(f0R[1]+f0L[1])-0.5*f0R[0]+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*f0L[3]-0.5*(f0R[2]+f0L[2])-0.8660254037844386*f0R[1]+0.8660254037844386*f0L[1]+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*alpha1R[0]; 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3])-0.5*f0R[2]+0.5*f0L[2]+0.8660254037844386*(f0R[1]+f0L[1])-0.5*f0R[0]+0.5*f0L[0])*sgn(alphaOrdR)-0.8660254037844386*f0R[3]+0.8660254037844386*f0L[3]+0.5*(f0R[2]+f0L[2])-0.8660254037844386*f0R[1]+0.8660254037844386*f0L[1]+0.5*(f0R[0]+f0L[0])); 

  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*alpha1R[0]*fUp[0]; 
  incr[1] = -0.8660254037844386*alpha1R[0]*fUp[0]; 
  incr[2] = 0.5*alpha1R[0]*fUp[1]; 
  incr[3] = -0.8660254037844386*alpha1R[0]*fUp[1]; 

  outR[0] += incr[0]*rdx2R; 
  outR[1] += incr[1]*rdx2R; 
  outR[2] += incr[2]*rdx2R; 
  outR[3] += incr[3]*rdx2R; 

  outL[0] += -1.0*incr[0]*rdx2L; 
  outL[1] += incr[1]*rdx2L; 
  outL[2] += -1.0*incr[2]*rdx2L; 
  outL[3] += incr[3]*rdx2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSer_vpar_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmagR[1] = (1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_)/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.224744871391589*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double alpha0R[2]; 

  double alpha0UpR[2]; 

  double alpha1R[2]; 
  alpha1R[0] = (0.3535533905932737*hamil1R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alpha1R[1] = (0.3535533905932737*hamil1R[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 

  double alpha1UpR[2]; 
  alpha1UpR[0] = (0.3535533905932737*(hamil1R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = (0.3535533905932737*(hamil1R[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R-2.828427124746191*dApardtPrev[1]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.125*((3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamil1R[1])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double incrEmMod[4]; 
  // alpha0 == 0, so nothing to do for alpha0*f1 term.
  // linear term alpha1*f0 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha1UpR[0]-0.7071067811865475*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.8660254037844386*(f0R[2]+f0L[2])+0.5*f0R[1]-0.5*(f0L[1]+f0R[0])+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*(f0L[3]+f0R[2])+0.8660254037844386*f0L[2]-0.5*(f0R[1]+f0L[1])+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.5*(f0R[1]+f0R[0])+0.5*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.8660254037844386*(f0R[3]+f0R[2])+0.8660254037844386*(f0L[3]+f0L[2])+0.5*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = -1.224744871391589*fUp[1]; 

  incr[0] = 0.5*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.5*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.8660254037844386*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 

  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += incrEmMod[2]*rdvpar2L; 
  emModL[3] += incrEmMod[3]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
double EmLinearDeltaFGyrokineticGenGeoSurf1x1vSerStep2_vpar_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *f0L, const double *f0R, const double *f1L, const double *f1R, double *outL, double *outR, double *emModL, double *emModR) 
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

  double wxSqL = wL[0]*wL[0];
  double wxSqR = wR[0]*wR[0];
  double rdx2SqL = rdx2L*rdx2L;
  double rdx2SqR = rdx2R*rdx2R;
  double wvparSqL = wL[1]*wL[1];
  double wvparSqR = wR[1]*wR[1];
  double rdvpar2SqL = rdvpar2L*rdvpar2L;
  double rdvpar2SqR = rdvpar2R*rdvpar2R;

  double hamil0R[4]; 
  hamil0R[0] = (0.3333333333333333*m_*(3.0*rdvpar2SqR*wvparSqR+1.0))/rdvpar2SqR; 
  hamil0R[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double hamil1R[4]; 
  hamil1R[0] = 1.414213562373095*phi[0]*q_; 
  hamil1R[1] = 1.414213562373095*phi[1]*q_; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_)/q_; 
  BstarZdBmagR[1] = (1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_)/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagEMR[4]; 
  BstarZdBmagEMR[0] = 1.224744871391589*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R; 
  BstarZdBmagEMR[1] = 1.224744871391589*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R; 

  double alpha0R[2]; 

  double alpha0UpR[2]; 

  double alpha1R[2]; 
  alpha1R[0] = -(1.0*dApardt[0]*q_)/m_; 
  alpha1R[1] = -(1.0*dApardt[1]*q_)/m_; 

  double alpha1UpR[2]; 
  alpha1UpR[0] = (0.3535533905932737*(hamil1R[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alpha1UpR[1] = (0.3535533905932737*(hamil1R[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R-2.828427124746191*dApardtPrev[1]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.125*((3.0*hamil1R[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamil1R[1])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alpha1UpR[0]-0.7071067811865475*alpha1UpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(f0R[3]+f0L[3]))+0.8660254037844386*(f0R[2]+f0L[2])+0.5*f0R[1]-0.5*(f0L[1]+f0R[0])+0.5*f0L[0])*sgn(alphaOrdR)+0.8660254037844386*f0R[3]-0.8660254037844386*(f0L[3]+f0R[2])+0.8660254037844386*f0L[2]-0.5*(f0R[1]+f0L[1])+0.5*(f0R[0]+f0L[0])); 
  alphaOrdR = 0.7071067811865475*(alpha1UpR[1]+alpha1UpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(f0R[3]+f0L[3]+f0R[2]+f0L[2])-0.5*(f0R[1]+f0R[0])+0.5*(f0L[1]+f0L[0]))*sgn(alphaOrdR)-0.8660254037844386*(f0R[3]+f0R[2])+0.8660254037844386*(f0L[3]+f0L[2])+0.5*(f0R[1]+f0L[1]+f0R[0]+f0L[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[1] = 0.5*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 
  incr[2] = -0.8660254037844386*(alpha1R[1]*fUp[1]+alpha1R[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alpha1R[0]*fUp[1]+fUp[0]*alpha1R[1]); 


  outR[0] += incr[0]*rdvpar2R; 
  outR[1] += incr[1]*rdvpar2R; 
  outR[2] += incr[2]*rdvpar2R; 
  outR[3] += incr[3]*rdvpar2R; 

  outL[0] += -1.0*incr[0]*rdvpar2L; 
  outL[1] += -1.0*incr[1]*rdvpar2L; 
  outL[2] += incr[2]*rdvpar2L; 
  outL[3] += incr[3]*rdvpar2L; 

  return std::abs(alphaSurfAvgR); 
} 
