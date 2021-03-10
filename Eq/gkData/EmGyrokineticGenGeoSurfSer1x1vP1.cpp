#include <GyrokineticModDecl.h>
double EmGyrokineticGenGeoSurf1x1vSer_x_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = 0.7071067811865475*jacobTotInv[0]*(1.732050807568877*b_y[0]*Apar[1]*rdx2R+1.414213562373095*cmag[0]); 

  double BstarZdBmagL[4]; 
  BstarZdBmagL[0] = -0.7071067811865475*jacobTotInv[0]*(1.732050807568877*b_y[0]*AparL[1]*rdx2R-1.414213562373095*cmag[0]); 

  double alphaR[2]; 
  alphaR[0] = (0.3061862178478973*(BstarZdBmagR[0]+BstarZdBmagL[0])*hamilR[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(1.732050807568877*BstarZdBmagR[0]+1.732050807568877*BstarZdBmagL[0])*hamilR[2]*rdvpar2R)/m_; 

  double incr[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.3535533905932737*alphaR[0]*(1.732050807568877*fL[1]+fL[0]); 
  incr[1] = -0.3535533905932737*alphaR[0]*(3.0*fL[1]+1.732050807568877*fL[0]); 
  incr[2] = 0.3535533905932737*alphaR[0]*(1.732050807568877*fL[3]+fL[2]); 
  incr[3] = -0.3535533905932737*alphaR[0]*(3.0*fL[3]+1.732050807568877*fL[2]); 
  } else { 
  incr[0] = -0.3535533905932737*alphaR[0]*(1.732050807568877*fR[1]-1.0*fR[0]); 
  incr[1] = 0.3535533905932737*alphaR[0]*(3.0*fR[1]-1.732050807568877*fR[0]); 
  incr[2] = -0.3535533905932737*alphaR[0]*(1.732050807568877*fR[3]-1.0*fR[2]); 
  incr[3] = 0.3535533905932737*alphaR[0]*(3.0*fR[3]-1.732050807568877*fR[2]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaR[0]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.5*fR[2]-0.5*fL[2]+0.8660254037844386*(fR[1]+fL[1])-0.5*fR[0]+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*fL[3]-0.5*(fR[2]+fL[2])-0.8660254037844386*fR[1]+0.8660254037844386*fL[1]+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*alphaR[0]; 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3])-0.5*fR[2]+0.5*fL[2]+0.8660254037844386*(fR[1]+fL[1])-0.5*fR[0]+0.5*fL[0])*sgn(alphaOrdR)-0.8660254037844386*fR[3]+0.8660254037844386*fL[3]+0.5*(fR[2]+fL[2])-0.8660254037844386*fR[1]+0.8660254037844386*fL[1]+0.5*(fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*alphaR[0]*fUp[0]; 
  incr[1] = -0.8660254037844386*alphaR[0]*fUp[0]; 
  incr[2] = 0.5*alphaR[0]*fUp[1]; 
  incr[3] = -0.8660254037844386*alphaR[0]*fUp[1]; 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSer_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = 0.7071067811865475*jacobTotInv[0]*(1.732050807568877*b_y[0]*Apar[1]*rdx2R+1.414213562373095*cmag[0]); 

  double alphaR[2]; 
  alphaR[0] = -(0.6123724356957944*BstarZdBmagR[0]*hamilR[1]*rdx2R)/m_; 

  double alphaUpR[2]; 
  alphaUpR[0] = -(0.3535533905932737*(1.732050807568877*BstarZdBmagR[0]*hamilR[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = -(1.0*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.125*(1.732050807568877*BstarZdBmagR[0]*hamilR[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double incrEmMod[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.3535533905932737*alphaR[0]*(1.732050807568877*fL[2]+fL[0]); 
  incr[1] = 0.3535533905932737*alphaR[0]*(1.732050807568877*fL[3]+fL[1]); 
  incr[2] = -0.3535533905932737*alphaR[0]*(3.0*fL[2]+1.732050807568877*fL[0]); 
  incr[3] = -0.3535533905932737*alphaR[0]*(3.0*fL[3]+1.732050807568877*fL[1]); 
  incrEmMod[0] = 0.5*(1.732050807568877*fL[2]+fL[0]); 
  incrEmMod[1] = 0.5*(1.732050807568877*fL[3]+fL[1]); 
  incrEmMod[2] = -0.5*(3.0*fL[2]+1.732050807568877*fL[0]); 
  incrEmMod[3] = -0.5*(3.0*fL[3]+1.732050807568877*fL[1]); 
  } else { 
  incrEmMod[0] = -0.5*(1.732050807568877*fR[2]-1.0*fR[0]); 
  incrEmMod[1] = -0.5*(1.732050807568877*fR[3]-1.0*fR[1]); 
  incrEmMod[2] = 0.5*(3.0*fR[2]-1.732050807568877*fR[0]); 
  incrEmMod[3] = 0.5*(3.0*fR[3]-1.732050807568877*fR[1]); 
  incr[0] = -0.3535533905932737*alphaR[0]*(1.732050807568877*fR[2]-1.0*fR[0]); 
  incr[1] = -0.3535533905932737*alphaR[0]*(1.732050807568877*fR[3]-1.0*fR[1]); 
  incr[2] = 0.3535533905932737*alphaR[0]*(3.0*fR[2]-1.732050807568877*fR[0]); 
  incr[3] = 0.3535533905932737*alphaR[0]*(3.0*fR[3]-1.732050807568877*fR[1]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaUpR[0]-0.7071067811865475*alphaUpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.8660254037844386*(fR[2]+fL[2])+0.5*fR[1]-0.5*(fL[1]+fR[0])+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*(fL[3]+fR[2])+0.8660254037844386*fL[2]-0.5*(fR[1]+fL[1])+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*(alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3]+fR[2]+fL[2])-0.5*(fR[1]+fR[0])+0.5*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.8660254037844386*(fR[3]+fR[2])+0.8660254037844386*(fL[3]+fL[2])+0.5*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = -1.224744871391589*fUp[1]; 

  incr[0] = 0.5*alphaR[0]*fUp[0]; 
  incr[1] = 0.5*alphaR[0]*fUp[1]; 
  incr[2] = -0.8660254037844386*alphaR[0]*fUp[0]; 
  incr[3] = -0.8660254037844386*alphaR[0]*fUp[1]; 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSerStep2_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = 0.7071067811865475*jacobTotInv[0]*(1.732050807568877*b_y[0]*Apar[1]*rdx2R+1.414213562373095*cmag[0]); 

  double alphaR[2]; 
  alphaR[0] = -(1.0*dApardt[0]*q_)/m_; 
  alphaR[1] = -(1.0*dApardt[1]*q_)/m_; 

  double alphaUpR[2]; 
  alphaUpR[0] = -(0.3535533905932737*(1.732050807568877*BstarZdBmagR[0]*hamilR[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = -(1.0*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.125*(1.732050807568877*BstarZdBmagR[0]*hamilR[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(2.449489742783178*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+1.414213562373095*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.25*(2.449489742783178*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = -0.25*(4.242640687119286*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+2.449489742783178*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[3] = -0.25*(4.242640687119286*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  } else { 
  incr[0] = -0.25*(2.449489742783178*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-1.414213562373095*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.25*(2.449489742783178*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = 0.25*(4.242640687119286*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-2.449489742783178*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[3] = 0.25*(4.242640687119286*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaUpR[0]-0.7071067811865475*alphaUpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.8660254037844386*(fR[2]+fL[2])+0.5*fR[1]-0.5*(fL[1]+fR[0])+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*(fL[3]+fR[2])+0.8660254037844386*fL[2]-0.5*(fR[1]+fL[1])+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*(alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3]+fR[2]+fL[2])-0.5*(fR[1]+fR[0])+0.5*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.8660254037844386*(fR[3]+fR[2])+0.8660254037844386*(fL[3]+fL[2])+0.5*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.5*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = -0.8660254037844386*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSer_x_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (0.5*(3.464101615137754*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R+2.0*(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0]))))/q_; 
  BstarZdBmagR[1] = (0.5*(3.464101615137754*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R+2.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))))/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagL[4]; 
  BstarZdBmagL[0] = (0.5*(3.464101615137754*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+q_*(2.0*(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])-2.449489742783178*(2.0*AparL[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(b_y[0]*AparL[1]-1.0*AparL[0]*b_y[1]))*rdx2R)))/q_; 
  BstarZdBmagL[1] = (0.5*(3.464101615137754*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*((AparL[0]*b_y[1]-1.0*b_y[0]*AparL[1])*jacobTotInv[1]-2.0*jacobTotInv[0]*AparL[1]*b_y[1])*rdx2R+2.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))))/q_; 
  BstarZdBmagL[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagL[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[2]; 
  alphaR[0] = -(0.1767766952966368*(3.0*(BstarZdBmagR[1]+BstarZdBmagL[1])-1.732050807568877*(BstarZdBmagR[0]+BstarZdBmagL[0]))*hamilR[2]*rdvpar2R)/m_; 
  alphaR[1] = -(0.1767766952966368*hamilR[2]*(3.0*(BstarZdBmagR[3]+BstarZdBmagL[3])-1.732050807568877*(BstarZdBmagR[2]+BstarZdBmagL[2]))*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(3.0*BstarZdBmagR[1]+3.0*BstarZdBmagL[1]-1.732050807568877*BstarZdBmagR[0]-1.732050807568877*BstarZdBmagL[0])*hamilR[2]*rdvpar2R)/m_; 

  double incr[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(alphaR[1]*(2.449489742783178*fL[3]+1.414213562373095*fL[2])+alphaR[0]*(2.449489742783178*fL[1]+1.414213562373095*fL[0])); 
  incr[1] = -0.25*(alphaR[1]*(4.242640687119286*fL[3]+2.449489742783178*fL[2])+alphaR[0]*(4.242640687119286*fL[1]+2.449489742783178*fL[0])); 
  incr[2] = 0.25*(alphaR[0]*(2.449489742783178*fL[3]+1.414213562373095*fL[2])+alphaR[1]*(2.449489742783178*fL[1]+1.414213562373095*fL[0])); 
  incr[3] = -0.25*(alphaR[0]*(4.242640687119286*fL[3]+2.449489742783178*fL[2])+alphaR[1]*(4.242640687119286*fL[1]+2.449489742783178*fL[0])); 
  } else { 
  incr[0] = -0.25*(alphaR[1]*(2.449489742783178*fR[3]-1.414213562373095*fR[2])+alphaR[0]*(2.449489742783178*fR[1]-1.414213562373095*fR[0])); 
  incr[1] = 0.25*(alphaR[1]*(4.242640687119286*fR[3]-2.449489742783178*fR[2])+alphaR[0]*(4.242640687119286*fR[1]-2.449489742783178*fR[0])); 
  incr[2] = -0.25*(alphaR[0]*(2.449489742783178*fR[3]-1.414213562373095*fR[2])+alphaR[1]*(2.449489742783178*fR[1]-1.414213562373095*fR[0])); 
  incr[3] = 0.25*(alphaR[0]*(4.242640687119286*fR[3]-2.449489742783178*fR[2])+alphaR[1]*(4.242640687119286*fR[1]-2.449489742783178*fR[0])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaR[0]-0.7071067811865475*alphaR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.5*fR[2]-0.5*fL[2]+0.8660254037844386*(fR[1]+fL[1])-0.5*fR[0]+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*fL[3]-0.5*(fR[2]+fL[2])-0.8660254037844386*fR[1]+0.8660254037844386*fL[1]+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*(alphaR[1]+alphaR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3])-0.5*fR[2]+0.5*fL[2]+0.8660254037844386*(fR[1]+fL[1])-0.5*fR[0]+0.5*fL[0])*sgn(alphaOrdR)-0.8660254037844386*fR[3]+0.8660254037844386*fL[3]+0.5*(fR[2]+fL[2])-0.8660254037844386*fR[1]+0.8660254037844386*fL[1]+0.5*(fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = -0.8660254037844386*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[2] = 0.5*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[3] = -0.8660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSer_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (0.5*(3.464101615137754*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R+2.0*(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0]))))/q_; 
  BstarZdBmagR[1] = (0.5*(3.464101615137754*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R+2.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))))/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[2]; 
  alphaR[0] = (0.3535533905932737*hamilR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alphaR[1] = (0.3535533905932737*hamilR[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 

  double alphaUpR[2]; 
  alphaUpR[0] = (0.3535533905932737*(hamilR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = (0.3535533905932737*(hamilR[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R-2.828427124746191*dApardtPrev[1]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.125*((3.0*hamilR[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamilR[1])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double incrEmMod[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(2.449489742783178*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+1.414213562373095*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.25*(2.449489742783178*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = -0.25*(4.242640687119286*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+2.449489742783178*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[3] = -0.25*(4.242640687119286*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incrEmMod[0] = 0.5*(1.732050807568877*fL[2]+fL[0]); 
  incrEmMod[1] = 0.5*(1.732050807568877*fL[3]+fL[1]); 
  incrEmMod[2] = -0.5*(3.0*fL[2]+1.732050807568877*fL[0]); 
  incrEmMod[3] = -0.5*(3.0*fL[3]+1.732050807568877*fL[1]); 
  } else { 
  incrEmMod[0] = -0.5*(1.732050807568877*fR[2]-1.0*fR[0]); 
  incrEmMod[1] = -0.5*(1.732050807568877*fR[3]-1.0*fR[1]); 
  incrEmMod[2] = 0.5*(3.0*fR[2]-1.732050807568877*fR[0]); 
  incrEmMod[3] = 0.5*(3.0*fR[3]-1.732050807568877*fR[1]); 
  incr[0] = -0.25*(2.449489742783178*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-1.414213562373095*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.25*(2.449489742783178*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = 0.25*(4.242640687119286*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-2.449489742783178*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[3] = 0.25*(4.242640687119286*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaUpR[0]-0.7071067811865475*alphaUpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.8660254037844386*(fR[2]+fL[2])+0.5*fR[1]-0.5*(fL[1]+fR[0])+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*(fL[3]+fR[2])+0.8660254037844386*fL[2]-0.5*(fR[1]+fL[1])+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*(alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3]+fR[2]+fL[2])-0.5*(fR[1]+fR[0])+0.5*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.8660254037844386*(fR[3]+fR[2])+0.8660254037844386*(fL[3]+fL[2])+0.5*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = -1.224744871391589*fUp[1]; 

  incr[0] = 0.5*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.5*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = -0.8660254037844386*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSerStep2_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (0.5*(3.464101615137754*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R+2.0*(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0]))))/q_; 
  BstarZdBmagR[1] = (0.5*(3.464101615137754*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R+2.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))))/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[2]; 
  alphaR[0] = -(1.0*dApardt[0]*q_)/m_; 
  alphaR[1] = -(1.0*dApardt[1]*q_)/m_; 

  double alphaUpR[2]; 
  alphaUpR[0] = (0.3535533905932737*(hamilR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = (0.3535533905932737*(hamilR[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R-2.828427124746191*dApardtPrev[1]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.125*((3.0*hamilR[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamilR[1])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(2.449489742783178*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+1.414213562373095*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.25*(2.449489742783178*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = -0.25*(4.242640687119286*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+2.449489742783178*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[3] = -0.25*(4.242640687119286*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  } else { 
  incr[0] = -0.25*(2.449489742783178*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-1.414213562373095*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.25*(2.449489742783178*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = 0.25*(4.242640687119286*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-2.449489742783178*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[3] = 0.25*(4.242640687119286*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaUpR[0]-0.7071067811865475*alphaUpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.8660254037844386*(fR[2]+fL[2])+0.5*fR[1]-0.5*(fL[1]+fR[0])+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*(fL[3]+fR[2])+0.8660254037844386*fL[2]-0.5*(fR[1]+fL[1])+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*(alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3]+fR[2]+fL[2])-0.5*(fR[1]+fR[0])+0.5*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.8660254037844386*(fR[3]+fR[2])+0.8660254037844386*(fL[3]+fL[2])+0.5*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.5*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = -0.8660254037844386*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSer_x_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = 0.7071067811865475*jacobTotInv[0]*(1.732050807568877*b_y[0]*Apar[1]*rdx2R+1.414213562373095*cmag[0]); 

  double BstarZdBmagL[4]; 
  BstarZdBmagL[0] = -0.7071067811865475*jacobTotInv[0]*(1.732050807568877*b_y[0]*AparL[1]*rdx2R-1.414213562373095*cmag[0]); 

  double alphaR[2]; 
  alphaR[0] = (0.3061862178478973*(BstarZdBmagR[0]+BstarZdBmagL[0])*hamilR[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.0625*(1.732050807568877*BstarZdBmagR[0]+1.732050807568877*BstarZdBmagL[0])*hamilR[2]*rdvpar2R)/m_; 

  double incr[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.3535533905932737*alphaR[0]*(1.732050807568877*fL[1]+fL[0]); 
  incr[1] = -0.3535533905932737*alphaR[0]*(3.0*fL[1]+1.732050807568877*fL[0]); 
  incr[2] = 0.3535533905932737*alphaR[0]*(1.732050807568877*fL[3]+fL[2]); 
  incr[3] = -0.3535533905932737*alphaR[0]*(3.0*fL[3]+1.732050807568877*fL[2]); 
  } else { 
  incr[0] = -0.3535533905932737*alphaR[0]*(1.732050807568877*fR[1]-1.0*fR[0]); 
  incr[1] = 0.3535533905932737*alphaR[0]*(3.0*fR[1]-1.732050807568877*fR[0]); 
  incr[2] = -0.3535533905932737*alphaR[0]*(1.732050807568877*fR[3]-1.0*fR[2]); 
  incr[3] = 0.3535533905932737*alphaR[0]*(3.0*fR[3]-1.732050807568877*fR[2]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaR[0]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.5*fR[2]-0.5*fL[2]+0.8660254037844386*(fR[1]+fL[1])-0.5*fR[0]+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*fL[3]-0.5*(fR[2]+fL[2])-0.8660254037844386*fR[1]+0.8660254037844386*fL[1]+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*alphaR[0]; 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3])-0.5*fR[2]+0.5*fL[2]+0.8660254037844386*(fR[1]+fL[1])-0.5*fR[0]+0.5*fL[0])*sgn(alphaOrdR)-0.8660254037844386*fR[3]+0.8660254037844386*fL[3]+0.5*(fR[2]+fL[2])-0.8660254037844386*fR[1]+0.8660254037844386*fL[1]+0.5*(fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*alphaR[0]*fUp[0]; 
  incr[1] = -0.8660254037844386*alphaR[0]*fUp[0]; 
  incr[2] = 0.5*alphaR[0]*fUp[1]; 
  incr[3] = -0.8660254037844386*alphaR[0]*fUp[1]; 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSer_vpar_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = 0.7071067811865475*jacobTotInv[0]*(1.732050807568877*b_y[0]*Apar[1]*rdx2R+1.414213562373095*cmag[0]); 

  double alphaR[2]; 
  alphaR[0] = -(0.6123724356957944*BstarZdBmagR[0]*hamilR[1]*rdx2R)/m_; 

  double alphaUpR[2]; 
  alphaUpR[0] = -(0.3535533905932737*(1.732050807568877*BstarZdBmagR[0]*hamilR[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = -(1.0*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.125*(1.732050807568877*BstarZdBmagR[0]*hamilR[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double incrEmMod[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.3535533905932737*alphaR[0]*(1.732050807568877*fL[2]+fL[0]); 
  incr[1] = 0.3535533905932737*alphaR[0]*(1.732050807568877*fL[3]+fL[1]); 
  incr[2] = -0.3535533905932737*alphaR[0]*(3.0*fL[2]+1.732050807568877*fL[0]); 
  incr[3] = -0.3535533905932737*alphaR[0]*(3.0*fL[3]+1.732050807568877*fL[1]); 
  incrEmMod[0] = 0.5*(1.732050807568877*fL[2]+fL[0]); 
  incrEmMod[1] = 0.5*(1.732050807568877*fL[3]+fL[1]); 
  incrEmMod[2] = -0.5*(3.0*fL[2]+1.732050807568877*fL[0]); 
  incrEmMod[3] = -0.5*(3.0*fL[3]+1.732050807568877*fL[1]); 
  } else { 
  incrEmMod[0] = -0.5*(1.732050807568877*fR[2]-1.0*fR[0]); 
  incrEmMod[1] = -0.5*(1.732050807568877*fR[3]-1.0*fR[1]); 
  incrEmMod[2] = 0.5*(3.0*fR[2]-1.732050807568877*fR[0]); 
  incrEmMod[3] = 0.5*(3.0*fR[3]-1.732050807568877*fR[1]); 
  incr[0] = -0.3535533905932737*alphaR[0]*(1.732050807568877*fR[2]-1.0*fR[0]); 
  incr[1] = -0.3535533905932737*alphaR[0]*(1.732050807568877*fR[3]-1.0*fR[1]); 
  incr[2] = 0.3535533905932737*alphaR[0]*(3.0*fR[2]-1.732050807568877*fR[0]); 
  incr[3] = 0.3535533905932737*alphaR[0]*(3.0*fR[3]-1.732050807568877*fR[1]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaUpR[0]-0.7071067811865475*alphaUpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.8660254037844386*(fR[2]+fL[2])+0.5*fR[1]-0.5*(fL[1]+fR[0])+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*(fL[3]+fR[2])+0.8660254037844386*fL[2]-0.5*(fR[1]+fL[1])+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*(alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3]+fR[2]+fL[2])-0.5*(fR[1]+fR[0])+0.5*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.8660254037844386*(fR[3]+fR[2])+0.8660254037844386*(fL[3]+fL[2])+0.5*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = -1.224744871391589*fUp[1]; 

  incr[0] = 0.5*alphaR[0]*fUp[0]; 
  incr[1] = 0.5*alphaR[0]*fUp[1]; 
  incr[2] = -0.8660254037844386*alphaR[0]*fUp[0]; 
  incr[3] = -0.8660254037844386*alphaR[0]*fUp[1]; 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSerStep2_vpar_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = 0.7071067811865475*jacobTotInv[0]*(1.732050807568877*b_y[0]*Apar[1]*rdx2R+1.414213562373095*cmag[0]); 

  double alphaR[2]; 
  alphaR[0] = -(1.0*dApardt[0]*q_)/m_; 
  alphaR[1] = -(1.0*dApardt[1]*q_)/m_; 

  double alphaUpR[2]; 
  alphaUpR[0] = -(0.3535533905932737*(1.732050807568877*BstarZdBmagR[0]*hamilR[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = -(1.0*dApardtPrev[1]*q_)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.125*(1.732050807568877*BstarZdBmagR[0]*hamilR[1]*rdx2R+2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(2.449489742783178*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+1.414213562373095*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.25*(2.449489742783178*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = -0.25*(4.242640687119286*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+2.449489742783178*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[3] = -0.25*(4.242640687119286*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  } else { 
  incr[0] = -0.25*(2.449489742783178*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-1.414213562373095*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.25*(2.449489742783178*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = 0.25*(4.242640687119286*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-2.449489742783178*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[3] = 0.25*(4.242640687119286*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaUpR[0]-0.7071067811865475*alphaUpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.8660254037844386*(fR[2]+fL[2])+0.5*fR[1]-0.5*(fL[1]+fR[0])+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*(fL[3]+fR[2])+0.8660254037844386*fL[2]-0.5*(fR[1]+fL[1])+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*(alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3]+fR[2]+fL[2])-0.5*(fR[1]+fR[0])+0.5*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.8660254037844386*(fR[3]+fR[2])+0.8660254037844386*(fL[3]+fL[2])+0.5*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.5*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = -0.8660254037844386*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSer_x_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (0.5*(3.464101615137754*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R+2.0*(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0]))))/q_; 
  BstarZdBmagR[1] = (0.5*(3.464101615137754*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R+2.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))))/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double BstarZdBmagL[4]; 
  BstarZdBmagL[0] = (0.5*(3.464101615137754*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+q_*(2.0*(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])-2.449489742783178*(2.0*AparL[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(b_y[0]*AparL[1]-1.0*AparL[0]*b_y[1]))*rdx2R)))/q_; 
  BstarZdBmagL[1] = (0.5*(3.464101615137754*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*((AparL[0]*b_y[1]-1.0*b_y[0]*AparL[1])*jacobTotInv[1]-2.0*jacobTotInv[0]*AparL[1]*b_y[1])*rdx2R+2.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))))/q_; 
  BstarZdBmagL[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagL[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[2]; 
  alphaR[0] = -(0.1767766952966368*(3.0*(BstarZdBmagR[1]+BstarZdBmagL[1])-1.732050807568877*(BstarZdBmagR[0]+BstarZdBmagL[0]))*hamilR[2]*rdvpar2R)/m_; 
  alphaR[1] = -(0.1767766952966368*hamilR[2]*(3.0*(BstarZdBmagR[3]+BstarZdBmagL[3])-1.732050807568877*(BstarZdBmagR[2]+BstarZdBmagL[2]))*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(3.0*BstarZdBmagR[1]+3.0*BstarZdBmagL[1]-1.732050807568877*BstarZdBmagR[0]-1.732050807568877*BstarZdBmagL[0])*hamilR[2]*rdvpar2R)/m_; 

  double incr[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(alphaR[1]*(2.449489742783178*fL[3]+1.414213562373095*fL[2])+alphaR[0]*(2.449489742783178*fL[1]+1.414213562373095*fL[0])); 
  incr[1] = -0.25*(alphaR[1]*(4.242640687119286*fL[3]+2.449489742783178*fL[2])+alphaR[0]*(4.242640687119286*fL[1]+2.449489742783178*fL[0])); 
  incr[2] = 0.25*(alphaR[0]*(2.449489742783178*fL[3]+1.414213562373095*fL[2])+alphaR[1]*(2.449489742783178*fL[1]+1.414213562373095*fL[0])); 
  incr[3] = -0.25*(alphaR[0]*(4.242640687119286*fL[3]+2.449489742783178*fL[2])+alphaR[1]*(4.242640687119286*fL[1]+2.449489742783178*fL[0])); 
  } else { 
  incr[0] = -0.25*(alphaR[1]*(2.449489742783178*fR[3]-1.414213562373095*fR[2])+alphaR[0]*(2.449489742783178*fR[1]-1.414213562373095*fR[0])); 
  incr[1] = 0.25*(alphaR[1]*(4.242640687119286*fR[3]-2.449489742783178*fR[2])+alphaR[0]*(4.242640687119286*fR[1]-2.449489742783178*fR[0])); 
  incr[2] = -0.25*(alphaR[0]*(2.449489742783178*fR[3]-1.414213562373095*fR[2])+alphaR[1]*(2.449489742783178*fR[1]-1.414213562373095*fR[0])); 
  incr[3] = 0.25*(alphaR[0]*(4.242640687119286*fR[3]-2.449489742783178*fR[2])+alphaR[1]*(4.242640687119286*fR[1]-2.449489742783178*fR[0])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaR[0]-0.7071067811865475*alphaR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.5*fR[2]-0.5*fL[2]+0.8660254037844386*(fR[1]+fL[1])-0.5*fR[0]+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*fL[3]-0.5*(fR[2]+fL[2])-0.8660254037844386*fR[1]+0.8660254037844386*fL[1]+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*(alphaR[1]+alphaR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3])-0.5*fR[2]+0.5*fL[2]+0.8660254037844386*(fR[1]+fL[1])-0.5*fR[0]+0.5*fL[0])*sgn(alphaOrdR)-0.8660254037844386*fR[3]+0.8660254037844386*fL[3]+0.5*(fR[2]+fL[2])-0.8660254037844386*fR[1]+0.8660254037844386*fL[1]+0.5*(fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = -0.8660254037844386*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[2] = 0.5*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[3] = -0.8660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSer_vpar_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (0.5*(3.464101615137754*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R+2.0*(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0]))))/q_; 
  BstarZdBmagR[1] = (0.5*(3.464101615137754*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R+2.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))))/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[2]; 
  alphaR[0] = (0.3535533905932737*hamilR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alphaR[1] = (0.3535533905932737*hamilR[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 

  double alphaUpR[2]; 
  alphaUpR[0] = (0.3535533905932737*(hamilR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = (0.3535533905932737*(hamilR[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R-2.828427124746191*dApardtPrev[1]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.125*((3.0*hamilR[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamilR[1])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
  double incrEmMod[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(2.449489742783178*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+1.414213562373095*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.25*(2.449489742783178*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = -0.25*(4.242640687119286*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+2.449489742783178*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[3] = -0.25*(4.242640687119286*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incrEmMod[0] = 0.5*(1.732050807568877*fL[2]+fL[0]); 
  incrEmMod[1] = 0.5*(1.732050807568877*fL[3]+fL[1]); 
  incrEmMod[2] = -0.5*(3.0*fL[2]+1.732050807568877*fL[0]); 
  incrEmMod[3] = -0.5*(3.0*fL[3]+1.732050807568877*fL[1]); 
  } else { 
  incrEmMod[0] = -0.5*(1.732050807568877*fR[2]-1.0*fR[0]); 
  incrEmMod[1] = -0.5*(1.732050807568877*fR[3]-1.0*fR[1]); 
  incrEmMod[2] = 0.5*(3.0*fR[2]-1.732050807568877*fR[0]); 
  incrEmMod[3] = 0.5*(3.0*fR[3]-1.732050807568877*fR[1]); 
  incr[0] = -0.25*(2.449489742783178*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-1.414213562373095*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.25*(2.449489742783178*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = 0.25*(4.242640687119286*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-2.449489742783178*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[3] = 0.25*(4.242640687119286*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaUpR[0]-0.7071067811865475*alphaUpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.8660254037844386*(fR[2]+fL[2])+0.5*fR[1]-0.5*(fL[1]+fR[0])+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*(fL[3]+fR[2])+0.8660254037844386*fL[2]-0.5*(fR[1]+fL[1])+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*(alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3]+fR[2]+fL[2])-0.5*(fR[1]+fR[0])+0.5*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.8660254037844386*(fR[3]+fR[2])+0.8660254037844386*(fL[3]+fL[2])+0.5*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = -1.224744871391589*fUp[0]; 
  incrEmMod[3] = -1.224744871391589*fUp[1]; 

  incr[0] = 0.5*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.5*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = -0.8660254037844386*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 

#endif 
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
double EmGyrokineticGenGeoSurf1x1vSerStep2_vpar_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *AparL, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[4]; 
  hamilR[0] = (0.3333333333333333*(3.0*rdvpar2SqR*(m_*wvparSqR+1.414213562373095*phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 1.414213562373095*phi[1]*q_; 
  hamilR[2] = (1.154700538379252*m_*wvparR)/rdvpar2R; 

  double BstarZdBmagR[4]; 
  BstarZdBmagR[0] = (0.5*(3.464101615137754*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*(2.0*Apar[1]*b_y[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_y[1]+b_y[0]*Apar[1]))*rdx2R+2.0*(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0]))))/q_; 
  BstarZdBmagR[1] = (0.5*(3.464101615137754*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+q_*(2.449489742783178*((Apar[0]*b_y[1]+b_y[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_y[1])*rdx2R+2.0*(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1]))))/q_; 
  BstarZdBmagR[2] = (jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[3] = (b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[2]; 
  alphaR[0] = -(1.0*dApardt[0]*q_)/m_; 
  alphaR[1] = -(1.0*dApardt[1]*q_)/m_; 

  double alphaUpR[2]; 
  alphaUpR[0] = (0.3535533905932737*(hamilR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = (0.3535533905932737*(hamilR[1]*(3.0*BstarZdBmagR[3]-1.732050807568877*BstarZdBmagR[1])*rdx2R-2.828427124746191*dApardtPrev[1]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.125*((3.0*hamilR[1]*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0]*hamilR[1])*rdx2R-2.828427124746191*dApardtPrev[0]*q_))/m_; 

  double incr[4]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(2.449489742783178*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+1.414213562373095*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.25*(2.449489742783178*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = -0.25*(4.242640687119286*(alphaR[1]*fL[3]+alphaR[0]*fL[2])+2.449489742783178*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[3] = -0.25*(4.242640687119286*(alphaR[0]*fL[3]+alphaR[1]*fL[2])+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  } else { 
  incr[0] = -0.25*(2.449489742783178*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-1.414213562373095*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.25*(2.449489742783178*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = 0.25*(4.242640687119286*(alphaR[1]*fR[3]+alphaR[0]*fR[2])-2.449489742783178*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[3] = 0.25*(4.242640687119286*(alphaR[0]*fR[3]+alphaR[1]*fR[2])-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[2];
  alphaOrdR = 0.7071067811865475*alphaUpR[0]-0.7071067811865475*alphaUpR[1]; 
  fUpOrd[0] = 0.5*(((-0.8660254037844386*(fR[3]+fL[3]))+0.8660254037844386*(fR[2]+fL[2])+0.5*fR[1]-0.5*(fL[1]+fR[0])+0.5*fL[0])*sgn(alphaOrdR)+0.8660254037844386*fR[3]-0.8660254037844386*(fL[3]+fR[2])+0.8660254037844386*fL[2]-0.5*(fR[1]+fL[1])+0.5*(fR[0]+fL[0])); 
  alphaOrdR = 0.7071067811865475*(alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[1] = 0.5*((0.8660254037844386*(fR[3]+fL[3]+fR[2]+fL[2])-0.5*(fR[1]+fR[0])+0.5*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.8660254037844386*(fR[3]+fR[2])+0.8660254037844386*(fL[3]+fL[2])+0.5*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[2];
  fUp[0] = 0.7071067811865475*(fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.7071067811865475*(fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.5*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.5*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = -0.8660254037844386*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[3] = -0.8660254037844386*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 

#endif 
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
