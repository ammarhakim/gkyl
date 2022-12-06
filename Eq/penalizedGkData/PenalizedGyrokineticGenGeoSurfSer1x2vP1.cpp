#include <PenalizedGyrokineticModDecl.h>
double PenalizedGyrokineticGenGeoSurf1x2vSer_x_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *fL, const double *fR, double *outL, double *outR) 
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
  // penaltyChi: penalization factor multiplying the mirror force term.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuBR[8]; 
  hamilNoMuBR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2SqR; 
  hamilNoMuBR[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuBR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 

  double hamilMuBR[8]; 
  hamilMuBR[0] = 2.0*bmag[0]*wmuR; 
  hamilMuBR[3] = (1.154700538379252*bmag[0])/rdmu2R; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double alphaR[4]; 
  alphaR[0] = (0.4330127018922193*BstarZdBmagR[0]*hamilNoMuBR[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.1082531754730548*BstarZdBmagR[0]*hamilNoMuBR[2]*rdvpar2R)/m_; 

  double incr[8]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*alphaR[0]*(1.732050807568877*fL[1]+fL[0]); 
  incr[1] = -0.25*alphaR[0]*(3.0*fL[1]+1.732050807568877*fL[0]); 
  incr[2] = 0.25*alphaR[0]*(1.732050807568877*fL[4]+fL[2]); 
  incr[3] = 0.25*alphaR[0]*(1.732050807568877*fL[5]+fL[3]); 
  incr[4] = -0.25*alphaR[0]*(3.0*fL[4]+1.732050807568877*fL[2]); 
  incr[5] = -0.25*alphaR[0]*(3.0*fL[5]+1.732050807568877*fL[3]); 
  incr[6] = 0.25*alphaR[0]*(1.732050807568877*fL[7]+fL[6]); 
  incr[7] = -0.25*alphaR[0]*(3.0*fL[7]+1.732050807568877*fL[6]); 
  } else { 
  incr[0] = -0.25*alphaR[0]*(1.732050807568877*fR[1]-1.0*fR[0]); 
  incr[1] = 0.25*alphaR[0]*(3.0*fR[1]-1.732050807568877*fR[0]); 
  incr[2] = -0.25*alphaR[0]*(1.732050807568877*fR[4]-1.0*fR[2]); 
  incr[3] = -0.25*alphaR[0]*(1.732050807568877*fR[5]-1.0*fR[3]); 
  incr[4] = 0.25*alphaR[0]*(3.0*fR[4]-1.732050807568877*fR[2]); 
  incr[5] = 0.25*alphaR[0]*(3.0*fR[5]-1.732050807568877*fR[3]); 
  incr[6] = -0.25*alphaR[0]*(1.732050807568877*fR[7]-1.0*fR[6]); 
  incr[7] = 0.25*alphaR[0]*(3.0*fR[7]-1.732050807568877*fR[6]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.3535533905932737*fR[6]+0.3535533905932737*fL[6]-0.6123724356957944*(fR[5]+fL[5]+fR[4]+fL[4])+0.3535533905932737*(fR[3]+fR[2])-0.3535533905932737*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*fL[7]+0.3535533905932737*(fR[6]+fL[6])+0.6123724356957944*(fR[5]+fR[4])-0.6123724356957944*(fL[5]+fL[4])-0.3535533905932737*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.3535533905932737*fR[6]-0.3535533905932737*fL[6]-0.6123724356957944*(fR[5]+fL[5])+0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*(fL[3]+fR[2])+0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*fL[7]-0.3535533905932737*(fR[6]+fL[6])+0.6123724356957944*fR[5]-0.6123724356957944*(fL[5]+fR[4])+0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])+0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.3535533905932737*fR[6]-0.3535533905932737*fL[6]+0.6123724356957944*(fR[5]+fL[5])-0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*(fL[3]+fR[2])-0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*fL[7]-0.3535533905932737*(fR[6]+fL[6])-0.6123724356957944*fR[5]+0.6123724356957944*(fL[5]+fR[4])-0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[3] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.3535533905932737*fR[6]+0.3535533905932737*fL[6]+0.6123724356957944*(fR[5]+fL[5]+fR[4]+fL[4])-0.3535533905932737*(fR[3]+fR[2])+0.3535533905932737*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*fL[7]+0.3535533905932737*(fR[6]+fL[6])-0.6123724356957944*(fR[5]+fR[4])+0.6123724356957944*(fL[5]+fL[4])+0.3535533905932737*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*alphaR[0]*fUp[0]; 
  incr[1] = -0.6123724356957944*alphaR[0]*fUp[0]; 
  incr[2] = 0.3535533905932737*alphaR[0]*fUp[1]; 
  incr[3] = 0.3535533905932737*alphaR[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alphaR[0]*fUp[1]; 
  incr[5] = -0.6123724356957944*alphaR[0]*fUp[2]; 
  incr[6] = 0.3535533905932737*alphaR[0]*fUp[3]; 
  incr[7] = -0.6123724356957944*alphaR[0]*fUp[3]; 

#endif 
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
double PenalizedGyrokineticGenGeoSurf1x2vSer_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *fL, const double *fR, double *outL, double *outR) 
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
  // penaltyChi: penalization factor multiplying the mirror force term.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuBR[8]; 
  hamilNoMuBR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2SqR; 
  hamilNoMuBR[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuBR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 

  double hamilMuBR[8]; 
  hamilMuBR[0] = 2.0*bmag[0]*wmuR; 
  hamilMuBR[3] = (1.154700538379252*bmag[0])/rdmu2R; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double alphaRNoMuB[4]; 
  alphaRNoMuB[0] = -(0.4330127018922193*BstarZdBmagR[0]*hamilNoMuBR[1]*rdx2R)/m_; 

  double alphaRMuB[4]; 

  double alphaR[4]; 
  alphaR[0] = alphaRNoMuB[0]; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.1082531754730548*BstarZdBmagR[0]*hamilNoMuBR[1]*rdx2R)/m_; 

  double incr[8]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*alphaR[0]*(1.732050807568877*fL[2]+fL[0]); 
  incr[1] = 0.25*alphaR[0]*(1.732050807568877*fL[4]+fL[1]); 
  incr[2] = -0.25*alphaR[0]*(3.0*fL[2]+1.732050807568877*fL[0]); 
  incr[3] = 0.25*alphaR[0]*(1.732050807568877*fL[6]+fL[3]); 
  incr[4] = -0.25*alphaR[0]*(3.0*fL[4]+1.732050807568877*fL[1]); 
  incr[5] = 0.25*alphaR[0]*(1.732050807568877*fL[7]+fL[5]); 
  incr[6] = -0.25*alphaR[0]*(3.0*fL[6]+1.732050807568877*fL[3]); 
  incr[7] = -0.25*alphaR[0]*(3.0*fL[7]+1.732050807568877*fL[5]); 
  } else { 
  incr[0] = -0.25*alphaR[0]*(1.732050807568877*fR[2]-1.0*fR[0]); 
  incr[1] = -0.25*alphaR[0]*(1.732050807568877*fR[4]-1.0*fR[1]); 
  incr[2] = 0.25*alphaR[0]*(3.0*fR[2]-1.732050807568877*fR[0]); 
  incr[3] = -0.25*alphaR[0]*(1.732050807568877*fR[6]-1.0*fR[3]); 
  incr[4] = 0.25*alphaR[0]*(3.0*fR[4]-1.732050807568877*fR[1]); 
  incr[5] = -0.25*alphaR[0]*(1.732050807568877*fR[7]-1.0*fR[5]); 
  incr[6] = 0.25*alphaR[0]*(3.0*fR[6]-1.732050807568877*fR[3]); 
  incr[7] = 0.25*alphaR[0]*(3.0*fR[7]-1.732050807568877*fR[5]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.6123724356957944*(fR[6]+fL[6])-0.3535533905932737*fR[5]+0.3535533905932737*fL[5]-0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])+0.3535533905932737*fR[1]-0.3535533905932737*(fL[1]+fR[0])+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*(fL[7]+fR[6])-0.6123724356957944*fL[6]+0.3535533905932737*(fR[5]+fL[5])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]-0.3535533905932737*(fR[1]+fL[1])+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]+fR[6]+fL[6]))+0.3535533905932737*fR[5]-0.3535533905932737*fL[5]+0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])-0.3535533905932737*(fR[1]+fR[0])+0.3535533905932737*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.6123724356957944*(fR[7]+fR[6])-0.6123724356957944*(fL[7]+fL[6])-0.3535533905932737*(fR[5]+fL[5])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]+0.3535533905932737*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.6123724356957944*(fR[6]+fL[6])+0.3535533905932737*fR[5]-0.3535533905932737*fL[5]-0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])+0.3535533905932737*fR[1]-0.3535533905932737*(fL[1]+fR[0])+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*(fL[7]+fR[6])+0.6123724356957944*fL[6]-0.3535533905932737*(fR[5]+fL[5])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]-0.3535533905932737*(fR[1]+fL[1])+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[3] = 0.5*((0.6123724356957944*(fR[7]+fL[7]+fR[6]+fL[6])-0.3535533905932737*fR[5]+0.3535533905932737*fL[5]+0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])-0.3535533905932737*(fR[1]+fR[0])+0.3535533905932737*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.6123724356957944*(fR[7]+fR[6])+0.6123724356957944*(fL[7]+fL[6])+0.3535533905932737*(fR[5]+fL[5])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]+0.3535533905932737*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*alphaR[0]*fUp[0]; 
  incr[1] = 0.3535533905932737*alphaR[0]*fUp[1]; 
  incr[2] = -0.6123724356957944*alphaR[0]*fUp[0]; 
  incr[3] = 0.3535533905932737*alphaR[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alphaR[0]*fUp[1]; 
  incr[5] = 0.3535533905932737*alphaR[0]*fUp[3]; 
  incr[6] = -0.6123724356957944*alphaR[0]*fUp[2]; 
  incr[7] = -0.6123724356957944*alphaR[0]*fUp[3]; 

#endif 
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
double PenalizedGyrokineticGenGeoSurf1x2vSer_x_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *fL, const double *fR, double *outL, double *outR) 
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
  // penaltyChi: penalization factor multiplying the mirror force term.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuBR[8]; 
  hamilNoMuBR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2SqR; 
  hamilNoMuBR[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuBR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 

  double hamilMuBR[8]; 
  hamilMuBR[0] = 2.0*bmag[0]*wmuR; 
  hamilMuBR[1] = 2.0*bmag[1]*wmuR; 
  hamilMuBR[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilMuBR[5] = (1.154700538379252*bmag[1])/rdmu2R; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[4]; 
  alphaR[0] = -(0.25*(3.0*BstarZdBmagR[1]-1.732050807568877*BstarZdBmagR[0])*hamilNoMuBR[2]*rdvpar2R)/m_; 
  alphaR[1] = -(0.25*hamilNoMuBR[2]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[2])*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(3.0*BstarZdBmagR[1]-1.732050807568877*BstarZdBmagR[0])*hamilNoMuBR[2]*rdvpar2R)/m_; 

  double incr[8]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(alphaR[1]*(1.732050807568877*fL[4]+fL[2])+alphaR[0]*(1.732050807568877*fL[1]+fL[0])); 
  incr[1] = -0.25*(alphaR[1]*(3.0*fL[4]+1.732050807568877*fL[2])+alphaR[0]*(3.0*fL[1]+1.732050807568877*fL[0])); 
  incr[2] = 0.25*(alphaR[0]*(1.732050807568877*fL[4]+fL[2])+alphaR[1]*(1.732050807568877*fL[1]+fL[0])); 
  incr[3] = 0.25*(alphaR[1]*(1.732050807568877*fL[7]+fL[6])+alphaR[0]*(1.732050807568877*fL[5]+fL[3])); 
  incr[4] = -0.25*(alphaR[0]*(3.0*fL[4]+1.732050807568877*fL[2])+alphaR[1]*(3.0*fL[1]+1.732050807568877*fL[0])); 
  incr[5] = -0.25*(alphaR[1]*(3.0*fL[7]+1.732050807568877*fL[6])+alphaR[0]*(3.0*fL[5]+1.732050807568877*fL[3])); 
  incr[6] = 0.25*(alphaR[0]*(1.732050807568877*fL[7]+fL[6])+alphaR[1]*(1.732050807568877*fL[5]+fL[3])); 
  incr[7] = -0.25*(alphaR[0]*(3.0*fL[7]+1.732050807568877*fL[6])+alphaR[1]*(3.0*fL[5]+1.732050807568877*fL[3])); 
  } else { 
  incr[0] = -0.25*(alphaR[1]*(1.732050807568877*fR[4]-1.0*fR[2])+alphaR[0]*(1.732050807568877*fR[1]-1.0*fR[0])); 
  incr[1] = 0.25*(alphaR[1]*(3.0*fR[4]-1.732050807568877*fR[2])+alphaR[0]*(3.0*fR[1]-1.732050807568877*fR[0])); 
  incr[2] = -0.25*(alphaR[0]*(1.732050807568877*fR[4]-1.0*fR[2])+alphaR[1]*(1.732050807568877*fR[1]-1.0*fR[0])); 
  incr[3] = -0.25*(alphaR[1]*(1.732050807568877*fR[7]-1.0*fR[6])+alphaR[0]*(1.732050807568877*fR[5]-1.0*fR[3])); 
  incr[4] = 0.25*(alphaR[0]*(3.0*fR[4]-1.732050807568877*fR[2])+alphaR[1]*(3.0*fR[1]-1.732050807568877*fR[0])); 
  incr[5] = 0.25*(alphaR[1]*(3.0*fR[7]-1.732050807568877*fR[6])+alphaR[0]*(3.0*fR[5]-1.732050807568877*fR[3])); 
  incr[6] = -0.25*(alphaR[0]*(1.732050807568877*fR[7]-1.0*fR[6])+alphaR[1]*(1.732050807568877*fR[5]-1.0*fR[3])); 
  incr[7] = 0.25*(alphaR[0]*(3.0*fR[7]-1.732050807568877*fR[6])+alphaR[1]*(3.0*fR[5]-1.732050807568877*fR[3])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alphaR[0]-0.5*alphaR[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.3535533905932737*fR[6]+0.3535533905932737*fL[6]-0.6123724356957944*(fR[5]+fL[5]+fR[4]+fL[4])+0.3535533905932737*(fR[3]+fR[2])-0.3535533905932737*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*fL[7]+0.3535533905932737*(fR[6]+fL[6])+0.6123724356957944*(fR[5]+fR[4])-0.6123724356957944*(fL[5]+fL[4])-0.3535533905932737*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[1]+alphaR[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.3535533905932737*fR[6]-0.3535533905932737*fL[6]-0.6123724356957944*(fR[5]+fL[5])+0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*(fL[3]+fR[2])+0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*fL[7]-0.3535533905932737*(fR[6]+fL[6])+0.6123724356957944*fR[5]-0.6123724356957944*(fL[5]+fR[4])+0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])+0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.5*alphaR[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.3535533905932737*fR[6]-0.3535533905932737*fL[6]+0.6123724356957944*(fR[5]+fL[5])-0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*(fL[3]+fR[2])-0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*fL[7]-0.3535533905932737*(fR[6]+fL[6])-0.6123724356957944*fR[5]+0.6123724356957944*(fL[5]+fR[4])-0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[1]+alphaR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.3535533905932737*fR[6]+0.3535533905932737*fL[6]+0.6123724356957944*(fR[5]+fL[5]+fR[4]+fL[4])-0.3535533905932737*(fR[3]+fR[2])+0.3535533905932737*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*fL[7]+0.3535533905932737*(fR[6]+fL[6])-0.6123724356957944*(fR[5]+fR[4])+0.6123724356957944*(fL[5]+fL[4])+0.3535533905932737*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = -0.6123724356957944*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[2] = 0.3535533905932737*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[3] = 0.3535533905932737*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2]); 
  incr[4] = -0.6123724356957944*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[5] = -0.6123724356957944*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2]); 
  incr[6] = 0.3535533905932737*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2]); 
  incr[7] = -0.6123724356957944*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2]); 

#endif 
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
double PenalizedGyrokineticGenGeoSurf1x2vSer_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *fL, const double *fR, double *outL, double *outR) 
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
  // penaltyChi: penalization factor multiplying the mirror force term.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuBR[8]; 
  hamilNoMuBR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2SqR; 
  hamilNoMuBR[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuBR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 

  double hamilMuBR[8]; 
  hamilMuBR[0] = 2.0*bmag[0]*wmuR; 
  hamilMuBR[1] = 2.0*bmag[1]*wmuR; 
  hamilMuBR[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilMuBR[5] = (1.154700538379252*bmag[1])/rdmu2R; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaRNoMuB[4]; 
  alphaRNoMuB[0] = (0.25*hamilNoMuBR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alphaRNoMuB[1] = (0.25*hamilNoMuBR[1]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 

  double alphaRMuB[4]; 
  alphaRMuB[0] = (0.25*hamilMuBR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alphaRMuB[1] = (0.25*hamilMuBR[1]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 
  alphaRMuB[2] = (0.25*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*hamilMuBR[5]*rdx2R)/m_; 
  alphaRMuB[3] = (0.25*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*hamilMuBR[5]*rdx2R)/m_; 

  double alphaR[4]; 
  alphaR[0] = 0.25*(1.414213562373095*(alphaRMuB[3]*penaltyChi[5]+alphaRMuB[2]*penaltyChi[3]+alphaRMuB[1]*penaltyChi[1]+alphaRMuB[0]*penaltyChi[0])+4.0*alphaRNoMuB[0]); 
  alphaR[1] = 0.25*(1.414213562373095*(alphaRMuB[2]*penaltyChi[5]+alphaRMuB[3]*penaltyChi[3]+alphaRMuB[0]*penaltyChi[1])+4.0*alphaRNoMuB[1]+1.414213562373095*penaltyChi[0]*alphaRMuB[1]); 
  alphaR[2] = 0.3535533905932738*(alphaRMuB[1]*penaltyChi[5]+alphaRMuB[0]*penaltyChi[3]+penaltyChi[1]*alphaRMuB[3]+penaltyChi[0]*alphaRMuB[2]); 
  alphaR[3] = 0.3535533905932738*(alphaRMuB[0]*penaltyChi[5]+alphaRMuB[1]*penaltyChi[3]+penaltyChi[0]*alphaRMuB[3]+penaltyChi[1]*alphaRMuB[2]); 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.015625*((4.242640687119286*BstarZdBmagR[4]-2.449489742783178*BstarZdBmagR[1])*hamilMuBR[5]*penaltyChi[5]+(4.242640687119286*BstarZdBmagR[2]-2.449489742783178*BstarZdBmagR[0])*penaltyChi[3]*hamilMuBR[5]+4.242640687119286*hamilMuBR[1]*penaltyChi[1]*BstarZdBmagR[4]+(12.0*hamilNoMuBR[1]+4.242640687119286*penaltyChi[0]*hamilMuBR[1])*BstarZdBmagR[2]-2.449489742783178*BstarZdBmagR[1]*hamilMuBR[1]*penaltyChi[1]-6.928203230275509*BstarZdBmagR[0]*hamilNoMuBR[1]-2.449489742783178*BstarZdBmagR[0]*penaltyChi[0]*hamilMuBR[1])*rdx2R)/m_; 

  double incr[8]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(1.732050807568877*(alphaR[3]*fL[7]+alphaR[2]*fL[6])+alphaR[3]*fL[5]+1.732050807568877*alphaR[1]*fL[4]+alphaR[2]*fL[3]+1.732050807568877*alphaR[0]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0]); 
  incr[1] = 0.25*(1.732050807568877*(alphaR[2]*fL[7]+alphaR[3]*fL[6])+alphaR[2]*fL[5]+1.732050807568877*alphaR[0]*fL[4]+alphaR[3]*fL[3]+1.732050807568877*alphaR[1]*fL[2]+alphaR[0]*fL[1]+fL[0]*alphaR[1]); 
  incr[2] = -0.25*(3.0*(alphaR[3]*fL[7]+alphaR[2]*fL[6])+1.732050807568877*alphaR[3]*fL[5]+3.0*alphaR[1]*fL[4]+1.732050807568877*alphaR[2]*fL[3]+3.0*alphaR[0]*fL[2]+1.732050807568877*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[3] = 0.25*(1.732050807568877*(alphaR[1]*fL[7]+alphaR[0]*fL[6])+alphaR[1]*fL[5]+1.732050807568877*alphaR[3]*fL[4]+alphaR[0]*fL[3]+fL[1]*alphaR[3]+alphaR[2]*(1.732050807568877*fL[2]+fL[0])); 
  incr[4] = -0.25*(3.0*(alphaR[2]*fL[7]+alphaR[3]*fL[6])+1.732050807568877*alphaR[2]*fL[5]+3.0*alphaR[0]*fL[4]+1.732050807568877*alphaR[3]*fL[3]+3.0*alphaR[1]*fL[2]+1.732050807568877*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[5] = 0.25*(1.732050807568877*(alphaR[0]*fL[7]+alphaR[1]*fL[6])+alphaR[0]*fL[5]+1.732050807568877*alphaR[2]*fL[4]+alphaR[1]*fL[3]+(1.732050807568877*fL[2]+fL[0])*alphaR[3]+fL[1]*alphaR[2]); 
  incr[6] = -0.25*(3.0*(alphaR[1]*fL[7]+alphaR[0]*fL[6])+1.732050807568877*alphaR[1]*fL[5]+3.0*alphaR[3]*fL[4]+1.732050807568877*(alphaR[0]*fL[3]+fL[1]*alphaR[3])+alphaR[2]*(3.0*fL[2]+1.732050807568877*fL[0])); 
  incr[7] = -0.25*(3.0*(alphaR[0]*fL[7]+alphaR[1]*fL[6])+1.732050807568877*alphaR[0]*fL[5]+3.0*alphaR[2]*fL[4]+1.732050807568877*alphaR[1]*fL[3]+(3.0*fL[2]+1.732050807568877*fL[0])*alphaR[3]+1.732050807568877*fL[1]*alphaR[2]); 
  } else { 
  incr[0] = -0.25*(1.732050807568877*(alphaR[3]*fR[7]+alphaR[2]*fR[6])-1.0*alphaR[3]*fR[5]+1.732050807568877*alphaR[1]*fR[4]-1.0*alphaR[2]*fR[3]+1.732050807568877*alphaR[0]*fR[2]-1.0*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.25*(1.732050807568877*(alphaR[2]*fR[7]+alphaR[3]*fR[6])-1.0*alphaR[2]*fR[5]+1.732050807568877*alphaR[0]*fR[4]-1.0*alphaR[3]*fR[3]+1.732050807568877*alphaR[1]*fR[2]-1.0*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = 0.25*(3.0*(alphaR[3]*fR[7]+alphaR[2]*fR[6])-1.732050807568877*alphaR[3]*fR[5]+3.0*alphaR[1]*fR[4]-1.732050807568877*alphaR[2]*fR[3]+3.0*alphaR[0]*fR[2]-1.732050807568877*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[3] = -0.25*(1.732050807568877*(alphaR[1]*fR[7]+alphaR[0]*fR[6])-1.0*alphaR[1]*fR[5]+1.732050807568877*alphaR[3]*fR[4]-1.0*(alphaR[0]*fR[3]+fR[1]*alphaR[3])+alphaR[2]*(1.732050807568877*fR[2]-1.0*fR[0])); 
  incr[4] = 0.25*(3.0*(alphaR[2]*fR[7]+alphaR[3]*fR[6])-1.732050807568877*alphaR[2]*fR[5]+3.0*alphaR[0]*fR[4]-1.732050807568877*alphaR[3]*fR[3]+3.0*alphaR[1]*fR[2]-1.732050807568877*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[5] = -0.25*(1.732050807568877*(alphaR[0]*fR[7]+alphaR[1]*fR[6])-1.0*alphaR[0]*fR[5]+1.732050807568877*alphaR[2]*fR[4]-1.0*alphaR[1]*fR[3]+(1.732050807568877*fR[2]-1.0*fR[0])*alphaR[3]-1.0*fR[1]*alphaR[2]); 
  incr[6] = 0.25*(3.0*(alphaR[1]*fR[7]+alphaR[0]*fR[6])-1.732050807568877*alphaR[1]*fR[5]+3.0*alphaR[3]*fR[4]-1.732050807568877*(alphaR[0]*fR[3]+fR[1]*alphaR[3])+alphaR[2]*(3.0*fR[2]-1.732050807568877*fR[0])); 
  incr[7] = 0.25*(3.0*(alphaR[0]*fR[7]+alphaR[1]*fR[6])-1.732050807568877*alphaR[0]*fR[5]+3.0*alphaR[2]*fR[4]-1.732050807568877*alphaR[1]*fR[3]+(3.0*fR[2]-1.732050807568877*fR[0])*alphaR[3]-1.732050807568877*fR[1]*alphaR[2]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alphaR[3]-0.5*(alphaR[2]+alphaR[1])+0.5*alphaR[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.6123724356957944*(fR[6]+fL[6])-0.3535533905932737*fR[5]+0.3535533905932737*fL[5]-0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])+0.3535533905932737*fR[1]-0.3535533905932737*(fL[1]+fR[0])+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*(fL[7]+fR[6])-0.6123724356957944*fL[6]+0.3535533905932737*(fR[5]+fL[5])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]-0.3535533905932737*(fR[1]+fL[1])+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[1]+alphaR[0])-0.5*(alphaR[3]+alphaR[2]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]+fR[6]+fL[6]))+0.3535533905932737*fR[5]-0.3535533905932737*fL[5]+0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])-0.3535533905932737*(fR[1]+fR[0])+0.3535533905932737*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.6123724356957944*(fR[7]+fR[6])-0.6123724356957944*(fL[7]+fL[6])-0.3535533905932737*(fR[5]+fL[5])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]+0.3535533905932737*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.5*alphaR[3])+0.5*alphaR[2]-0.5*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.6123724356957944*(fR[6]+fL[6])+0.3535533905932737*fR[5]-0.3535533905932737*fL[5]-0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])+0.3535533905932737*fR[1]-0.3535533905932737*(fL[1]+fR[0])+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*(fL[7]+fR[6])+0.6123724356957944*fL[6]-0.3535533905932737*(fR[5]+fL[5])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]-0.3535533905932737*(fR[1]+fL[1])+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(fR[7]+fL[7]+fR[6]+fL[6])-0.3535533905932737*fR[5]+0.3535533905932737*fL[5]+0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])-0.3535533905932737*(fR[1]+fR[0])+0.3535533905932737*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.6123724356957944*(fR[7]+fR[6])+0.6123724356957944*(fL[7]+fL[6])+0.3535533905932737*(fR[5]+fL[5])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]+0.3535533905932737*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alphaR[3]*fUp[3]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.3535533905932737*(alphaR[2]*fUp[3]+fUp[2]*alphaR[3]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = -0.6123724356957944*(alphaR[3]*fUp[3]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[3] = 0.3535533905932737*(alphaR[1]*fUp[3]+fUp[1]*alphaR[3]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[4] = -0.6123724356957944*(alphaR[2]*fUp[3]+fUp[2]*alphaR[3]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[5] = 0.3535533905932737*(alphaR[0]*fUp[3]+fUp[0]*alphaR[3]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[6] = -0.6123724356957944*(alphaR[1]*fUp[3]+fUp[1]*alphaR[3]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[7] = -0.6123724356957944*(alphaR[0]*fUp[3]+fUp[0]*alphaR[3]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 

#endif 
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
double PenalizedGyrokineticGenGeoSurf1x2vSer_x_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *fL, const double *fR, double *outL, double *outR) 
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
  // penaltyChi: penalization factor multiplying the mirror force term.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuBR[8]; 
  hamilNoMuBR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2SqR; 
  hamilNoMuBR[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuBR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 

  double hamilMuBR[8]; 
  hamilMuBR[0] = 2.0*bmag[0]*wmuR; 
  hamilMuBR[3] = (1.154700538379252*bmag[0])/rdmu2R; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double alphaR[4]; 
  alphaR[0] = (0.4330127018922193*BstarZdBmagR[0]*hamilNoMuBR[2]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.1082531754730548*BstarZdBmagR[0]*hamilNoMuBR[2]*rdvpar2R)/m_; 

  double incr[8]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*alphaR[0]*(1.732050807568877*fL[1]+fL[0]); 
  incr[1] = -0.25*alphaR[0]*(3.0*fL[1]+1.732050807568877*fL[0]); 
  incr[2] = 0.25*alphaR[0]*(1.732050807568877*fL[4]+fL[2]); 
  incr[3] = 0.25*alphaR[0]*(1.732050807568877*fL[5]+fL[3]); 
  incr[4] = -0.25*alphaR[0]*(3.0*fL[4]+1.732050807568877*fL[2]); 
  incr[5] = -0.25*alphaR[0]*(3.0*fL[5]+1.732050807568877*fL[3]); 
  incr[6] = 0.25*alphaR[0]*(1.732050807568877*fL[7]+fL[6]); 
  incr[7] = -0.25*alphaR[0]*(3.0*fL[7]+1.732050807568877*fL[6]); 
  } else { 
  incr[0] = -0.25*alphaR[0]*(1.732050807568877*fR[1]-1.0*fR[0]); 
  incr[1] = 0.25*alphaR[0]*(3.0*fR[1]-1.732050807568877*fR[0]); 
  incr[2] = -0.25*alphaR[0]*(1.732050807568877*fR[4]-1.0*fR[2]); 
  incr[3] = -0.25*alphaR[0]*(1.732050807568877*fR[5]-1.0*fR[3]); 
  incr[4] = 0.25*alphaR[0]*(3.0*fR[4]-1.732050807568877*fR[2]); 
  incr[5] = 0.25*alphaR[0]*(3.0*fR[5]-1.732050807568877*fR[3]); 
  incr[6] = -0.25*alphaR[0]*(1.732050807568877*fR[7]-1.0*fR[6]); 
  incr[7] = 0.25*alphaR[0]*(3.0*fR[7]-1.732050807568877*fR[6]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.3535533905932737*fR[6]+0.3535533905932737*fL[6]-0.6123724356957944*(fR[5]+fL[5]+fR[4]+fL[4])+0.3535533905932737*(fR[3]+fR[2])-0.3535533905932737*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*fL[7]+0.3535533905932737*(fR[6]+fL[6])+0.6123724356957944*(fR[5]+fR[4])-0.6123724356957944*(fL[5]+fL[4])-0.3535533905932737*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.3535533905932737*fR[6]-0.3535533905932737*fL[6]-0.6123724356957944*(fR[5]+fL[5])+0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*(fL[3]+fR[2])+0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*fL[7]-0.3535533905932737*(fR[6]+fL[6])+0.6123724356957944*fR[5]-0.6123724356957944*(fL[5]+fR[4])+0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])+0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.3535533905932737*fR[6]-0.3535533905932737*fL[6]+0.6123724356957944*(fR[5]+fL[5])-0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*(fL[3]+fR[2])-0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*fL[7]-0.3535533905932737*(fR[6]+fL[6])-0.6123724356957944*fR[5]+0.6123724356957944*(fL[5]+fR[4])-0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[3] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.3535533905932737*fR[6]+0.3535533905932737*fL[6]+0.6123724356957944*(fR[5]+fL[5]+fR[4]+fL[4])-0.3535533905932737*(fR[3]+fR[2])+0.3535533905932737*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*fL[7]+0.3535533905932737*(fR[6]+fL[6])-0.6123724356957944*(fR[5]+fR[4])+0.6123724356957944*(fL[5]+fL[4])+0.3535533905932737*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*alphaR[0]*fUp[0]; 
  incr[1] = -0.6123724356957944*alphaR[0]*fUp[0]; 
  incr[2] = 0.3535533905932737*alphaR[0]*fUp[1]; 
  incr[3] = 0.3535533905932737*alphaR[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alphaR[0]*fUp[1]; 
  incr[5] = -0.6123724356957944*alphaR[0]*fUp[2]; 
  incr[6] = 0.3535533905932737*alphaR[0]*fUp[3]; 
  incr[7] = -0.6123724356957944*alphaR[0]*fUp[3]; 

#endif 
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
double PenalizedGyrokineticGenGeoSurf1x2vSer_vpar_P1_Bvarsz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *fL, const double *fR, double *outL, double *outR) 
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
  // penaltyChi: penalization factor multiplying the mirror force term.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuBR[8]; 
  hamilNoMuBR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2SqR; 
  hamilNoMuBR[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuBR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 

  double hamilMuBR[8]; 
  hamilMuBR[0] = 2.0*bmag[0]*wmuR; 
  hamilMuBR[3] = (1.154700538379252*bmag[0])/rdmu2R; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = 1.414213562373095*cmag[0]*jacobTotInv[0]; 

  double alphaRNoMuB[4]; 
  alphaRNoMuB[0] = -(0.4330127018922193*BstarZdBmagR[0]*hamilNoMuBR[1]*rdx2R)/m_; 

  double alphaRMuB[4]; 

  double alphaR[4]; 
  alphaR[0] = alphaRNoMuB[0]; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.1082531754730548*BstarZdBmagR[0]*hamilNoMuBR[1]*rdx2R)/m_; 

  double incr[8]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*alphaR[0]*(1.732050807568877*fL[2]+fL[0]); 
  incr[1] = 0.25*alphaR[0]*(1.732050807568877*fL[4]+fL[1]); 
  incr[2] = -0.25*alphaR[0]*(3.0*fL[2]+1.732050807568877*fL[0]); 
  incr[3] = 0.25*alphaR[0]*(1.732050807568877*fL[6]+fL[3]); 
  incr[4] = -0.25*alphaR[0]*(3.0*fL[4]+1.732050807568877*fL[1]); 
  incr[5] = 0.25*alphaR[0]*(1.732050807568877*fL[7]+fL[5]); 
  incr[6] = -0.25*alphaR[0]*(3.0*fL[6]+1.732050807568877*fL[3]); 
  incr[7] = -0.25*alphaR[0]*(3.0*fL[7]+1.732050807568877*fL[5]); 
  } else { 
  incr[0] = -0.25*alphaR[0]*(1.732050807568877*fR[2]-1.0*fR[0]); 
  incr[1] = -0.25*alphaR[0]*(1.732050807568877*fR[4]-1.0*fR[1]); 
  incr[2] = 0.25*alphaR[0]*(3.0*fR[2]-1.732050807568877*fR[0]); 
  incr[3] = -0.25*alphaR[0]*(1.732050807568877*fR[6]-1.0*fR[3]); 
  incr[4] = 0.25*alphaR[0]*(3.0*fR[4]-1.732050807568877*fR[1]); 
  incr[5] = -0.25*alphaR[0]*(1.732050807568877*fR[7]-1.0*fR[5]); 
  incr[6] = 0.25*alphaR[0]*(3.0*fR[6]-1.732050807568877*fR[3]); 
  incr[7] = 0.25*alphaR[0]*(3.0*fR[7]-1.732050807568877*fR[5]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.6123724356957944*(fR[6]+fL[6])-0.3535533905932737*fR[5]+0.3535533905932737*fL[5]-0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])+0.3535533905932737*fR[1]-0.3535533905932737*(fL[1]+fR[0])+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*(fL[7]+fR[6])-0.6123724356957944*fL[6]+0.3535533905932737*(fR[5]+fL[5])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]-0.3535533905932737*(fR[1]+fL[1])+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]+fR[6]+fL[6]))+0.3535533905932737*fR[5]-0.3535533905932737*fL[5]+0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])-0.3535533905932737*(fR[1]+fR[0])+0.3535533905932737*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.6123724356957944*(fR[7]+fR[6])-0.6123724356957944*(fL[7]+fL[6])-0.3535533905932737*(fR[5]+fL[5])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]+0.3535533905932737*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.6123724356957944*(fR[6]+fL[6])+0.3535533905932737*fR[5]-0.3535533905932737*fL[5]-0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])+0.3535533905932737*fR[1]-0.3535533905932737*(fL[1]+fR[0])+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*(fL[7]+fR[6])+0.6123724356957944*fL[6]-0.3535533905932737*(fR[5]+fL[5])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]-0.3535533905932737*(fR[1]+fL[1])+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]; 
  fUpOrd[3] = 0.5*((0.6123724356957944*(fR[7]+fL[7]+fR[6]+fL[6])-0.3535533905932737*fR[5]+0.3535533905932737*fL[5]+0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])-0.3535533905932737*(fR[1]+fR[0])+0.3535533905932737*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.6123724356957944*(fR[7]+fR[6])+0.6123724356957944*(fL[7]+fL[6])+0.3535533905932737*(fR[5]+fL[5])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]+0.3535533905932737*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*alphaR[0]*fUp[0]; 
  incr[1] = 0.3535533905932737*alphaR[0]*fUp[1]; 
  incr[2] = -0.6123724356957944*alphaR[0]*fUp[0]; 
  incr[3] = 0.3535533905932737*alphaR[0]*fUp[2]; 
  incr[4] = -0.6123724356957944*alphaR[0]*fUp[1]; 
  incr[5] = 0.3535533905932737*alphaR[0]*fUp[3]; 
  incr[6] = -0.6123724356957944*alphaR[0]*fUp[2]; 
  incr[7] = -0.6123724356957944*alphaR[0]*fUp[3]; 

#endif 
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
double PenalizedGyrokineticGenGeoSurf1x2vSer_x_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *fL, const double *fR, double *outL, double *outR) 
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
  // penaltyChi: penalization factor multiplying the mirror force term.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuBR[8]; 
  hamilNoMuBR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2SqR; 
  hamilNoMuBR[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuBR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 

  double hamilMuBR[8]; 
  hamilMuBR[0] = 2.0*bmag[0]*wmuR; 
  hamilMuBR[1] = 2.0*bmag[1]*wmuR; 
  hamilMuBR[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilMuBR[5] = (1.154700538379252*bmag[1])/rdmu2R; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[4]; 
  alphaR[0] = -(0.25*(3.0*BstarZdBmagR[1]-1.732050807568877*BstarZdBmagR[0])*hamilNoMuBR[2]*rdvpar2R)/m_; 
  alphaR[1] = -(0.25*hamilNoMuBR[2]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[2])*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.0625*(3.0*BstarZdBmagR[1]-1.732050807568877*BstarZdBmagR[0])*hamilNoMuBR[2]*rdvpar2R)/m_; 

  double incr[8]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(alphaR[1]*(1.732050807568877*fL[4]+fL[2])+alphaR[0]*(1.732050807568877*fL[1]+fL[0])); 
  incr[1] = -0.25*(alphaR[1]*(3.0*fL[4]+1.732050807568877*fL[2])+alphaR[0]*(3.0*fL[1]+1.732050807568877*fL[0])); 
  incr[2] = 0.25*(alphaR[0]*(1.732050807568877*fL[4]+fL[2])+alphaR[1]*(1.732050807568877*fL[1]+fL[0])); 
  incr[3] = 0.25*(alphaR[1]*(1.732050807568877*fL[7]+fL[6])+alphaR[0]*(1.732050807568877*fL[5]+fL[3])); 
  incr[4] = -0.25*(alphaR[0]*(3.0*fL[4]+1.732050807568877*fL[2])+alphaR[1]*(3.0*fL[1]+1.732050807568877*fL[0])); 
  incr[5] = -0.25*(alphaR[1]*(3.0*fL[7]+1.732050807568877*fL[6])+alphaR[0]*(3.0*fL[5]+1.732050807568877*fL[3])); 
  incr[6] = 0.25*(alphaR[0]*(1.732050807568877*fL[7]+fL[6])+alphaR[1]*(1.732050807568877*fL[5]+fL[3])); 
  incr[7] = -0.25*(alphaR[0]*(3.0*fL[7]+1.732050807568877*fL[6])+alphaR[1]*(3.0*fL[5]+1.732050807568877*fL[3])); 
  } else { 
  incr[0] = -0.25*(alphaR[1]*(1.732050807568877*fR[4]-1.0*fR[2])+alphaR[0]*(1.732050807568877*fR[1]-1.0*fR[0])); 
  incr[1] = 0.25*(alphaR[1]*(3.0*fR[4]-1.732050807568877*fR[2])+alphaR[0]*(3.0*fR[1]-1.732050807568877*fR[0])); 
  incr[2] = -0.25*(alphaR[0]*(1.732050807568877*fR[4]-1.0*fR[2])+alphaR[1]*(1.732050807568877*fR[1]-1.0*fR[0])); 
  incr[3] = -0.25*(alphaR[1]*(1.732050807568877*fR[7]-1.0*fR[6])+alphaR[0]*(1.732050807568877*fR[5]-1.0*fR[3])); 
  incr[4] = 0.25*(alphaR[0]*(3.0*fR[4]-1.732050807568877*fR[2])+alphaR[1]*(3.0*fR[1]-1.732050807568877*fR[0])); 
  incr[5] = 0.25*(alphaR[1]*(3.0*fR[7]-1.732050807568877*fR[6])+alphaR[0]*(3.0*fR[5]-1.732050807568877*fR[3])); 
  incr[6] = -0.25*(alphaR[0]*(1.732050807568877*fR[7]-1.0*fR[6])+alphaR[1]*(1.732050807568877*fR[5]-1.0*fR[3])); 
  incr[7] = 0.25*(alphaR[0]*(3.0*fR[7]-1.732050807568877*fR[6])+alphaR[1]*(3.0*fR[5]-1.732050807568877*fR[3])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alphaR[0]-0.5*alphaR[1]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.3535533905932737*fR[6]+0.3535533905932737*fL[6]-0.6123724356957944*(fR[5]+fL[5]+fR[4]+fL[4])+0.3535533905932737*(fR[3]+fR[2])-0.3535533905932737*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*fL[7]+0.3535533905932737*(fR[6]+fL[6])+0.6123724356957944*(fR[5]+fR[4])-0.6123724356957944*(fL[5]+fL[4])-0.3535533905932737*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[1]+alphaR[0]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.3535533905932737*fR[6]-0.3535533905932737*fL[6]-0.6123724356957944*(fR[5]+fL[5])+0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*(fL[3]+fR[2])+0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*fL[7]-0.3535533905932737*(fR[6]+fL[6])+0.6123724356957944*fR[5]-0.6123724356957944*(fL[5]+fR[4])+0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])+0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*alphaR[0]-0.5*alphaR[1]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.3535533905932737*fR[6]-0.3535533905932737*fL[6]+0.6123724356957944*(fR[5]+fL[5])-0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*(fL[3]+fR[2])-0.3535533905932737*fL[2]+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*fL[7]-0.3535533905932737*(fR[6]+fL[6])-0.6123724356957944*fR[5]+0.6123724356957944*(fL[5]+fR[4])-0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.3535533905932737*(fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[1]+alphaR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.3535533905932737*fR[6]+0.3535533905932737*fL[6]+0.6123724356957944*(fR[5]+fL[5]+fR[4]+fL[4])-0.3535533905932737*(fR[3]+fR[2])+0.3535533905932737*(fL[3]+fL[2])+0.6123724356957944*(fR[1]+fL[1])-0.3535533905932737*fR[0]+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*fL[7]+0.3535533905932737*(fR[6]+fL[6])-0.6123724356957944*(fR[5]+fR[4])+0.6123724356957944*(fL[5]+fL[4])+0.3535533905932737*(fR[3]+fL[3]+fR[2]+fL[2])-0.6123724356957944*fR[1]+0.6123724356957944*fL[1]+0.3535533905932737*(fR[0]+fL[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = -0.6123724356957944*(alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[2] = 0.3535533905932737*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[3] = 0.3535533905932737*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2]); 
  incr[4] = -0.6123724356957944*(alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[5] = -0.6123724356957944*(alphaR[1]*fUp[3]+alphaR[0]*fUp[2]); 
  incr[6] = 0.3535533905932737*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2]); 
  incr[7] = -0.6123724356957944*(alphaR[0]*fUp[3]+alphaR[1]*fUp[2]); 

#endif 
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
double PenalizedGyrokineticGenGeoSurf1x2vSer_vpar_P1_Bvarsxz(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *penaltyChi, const double *phi, const double *exbEnergy, const double *fL, const double *fR, double *outL, double *outR) 
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
  // penaltyChi: penalization factor multiplying the mirror force term.
  // phi: electrostatic potential .
  // exbEnergy: 0.5*m*v_E^2 term in the Hamiltonian.
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

  double hamilNoMuBR[8]; 
  hamilNoMuBR[0] = (0.2357022603955158*(3.0*rdvpar2SqR*(2.0*m_*wvparSqR+2.828427124746191*(phi[0]*q_-1.0*exbEnergy[0]))+2.0*m_))/rdvpar2SqR; 
  hamilNoMuBR[1] = 2.0*(phi[1]*q_-1.0*exbEnergy[1]); 
  hamilNoMuBR[2] = (1.632993161855453*m_*wvparR)/rdvpar2R; 

  double hamilMuBR[8]; 
  hamilMuBR[0] = 2.0*bmag[0]*wmuR; 
  hamilMuBR[1] = 2.0*bmag[1]*wmuR; 
  hamilMuBR[3] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilMuBR[5] = (1.154700538379252*bmag[1])/rdmu2R; 

  double BstarZdBmagR[8]; 
  BstarZdBmagR[0] = (1.414213562373095*(1.732050807568877*jacobTotInv[0]*b_y[1]*m_*rdx2R*wvparR+(cmag[1]*jacobTotInv[1]+cmag[0]*jacobTotInv[0])*q_))/q_; 
  BstarZdBmagR[1] = (1.414213562373095*(1.732050807568877*b_y[1]*jacobTotInv[1]*m_*rdx2R*wvparR+(cmag[0]*jacobTotInv[1]+jacobTotInv[0]*cmag[1])*q_))/q_; 
  BstarZdBmagR[2] = (1.414213562373095*jacobTotInv[0]*b_y[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarZdBmagR[4] = (1.414213562373095*b_y[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaRNoMuB[4]; 
  alphaRNoMuB[0] = (0.25*hamilNoMuBR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alphaRNoMuB[1] = (0.25*hamilNoMuBR[1]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 

  double alphaRMuB[4]; 
  alphaRMuB[0] = (0.25*hamilMuBR[1]*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*rdx2R)/m_; 
  alphaRMuB[1] = (0.25*hamilMuBR[1]*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*rdx2R)/m_; 
  alphaRMuB[2] = (0.25*(3.0*BstarZdBmagR[2]-1.732050807568877*BstarZdBmagR[0])*hamilMuBR[5]*rdx2R)/m_; 
  alphaRMuB[3] = (0.25*(3.0*BstarZdBmagR[4]-1.732050807568877*BstarZdBmagR[1])*hamilMuBR[5]*rdx2R)/m_; 

  double alphaR[4]; 
  alphaR[0] = 0.25*(1.414213562373095*(alphaRMuB[3]*penaltyChi[5]+alphaRMuB[2]*penaltyChi[3]+alphaRMuB[1]*penaltyChi[1]+alphaRMuB[0]*penaltyChi[0])+4.0*alphaRNoMuB[0]); 
  alphaR[1] = 0.25*(1.414213562373095*(alphaRMuB[2]*penaltyChi[5]+alphaRMuB[3]*penaltyChi[3]+alphaRMuB[0]*penaltyChi[1])+4.0*alphaRNoMuB[1]+1.414213562373095*penaltyChi[0]*alphaRMuB[1]); 
  alphaR[2] = 0.3535533905932738*(alphaRMuB[1]*penaltyChi[5]+alphaRMuB[0]*penaltyChi[3]+penaltyChi[1]*alphaRMuB[3]+penaltyChi[0]*alphaRMuB[2]); 
  alphaR[3] = 0.3535533905932738*(alphaRMuB[0]*penaltyChi[5]+alphaRMuB[1]*penaltyChi[3]+penaltyChi[0]*alphaRMuB[3]+penaltyChi[1]*alphaRMuB[2]); 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.015625*((4.242640687119286*BstarZdBmagR[4]-2.449489742783178*BstarZdBmagR[1])*hamilMuBR[5]*penaltyChi[5]+(4.242640687119286*BstarZdBmagR[2]-2.449489742783178*BstarZdBmagR[0])*penaltyChi[3]*hamilMuBR[5]+4.242640687119286*hamilMuBR[1]*penaltyChi[1]*BstarZdBmagR[4]+(12.0*hamilNoMuBR[1]+4.242640687119286*penaltyChi[0]*hamilMuBR[1])*BstarZdBmagR[2]-2.449489742783178*BstarZdBmagR[1]*hamilMuBR[1]*penaltyChi[1]-6.928203230275509*BstarZdBmagR[0]*hamilNoMuBR[1]-2.449489742783178*BstarZdBmagR[0]*penaltyChi[0]*hamilMuBR[1])*rdx2R)/m_; 

  double incr[8]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.25*(1.732050807568877*(alphaR[3]*fL[7]+alphaR[2]*fL[6])+alphaR[3]*fL[5]+1.732050807568877*alphaR[1]*fL[4]+alphaR[2]*fL[3]+1.732050807568877*alphaR[0]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0]); 
  incr[1] = 0.25*(1.732050807568877*(alphaR[2]*fL[7]+alphaR[3]*fL[6])+alphaR[2]*fL[5]+1.732050807568877*alphaR[0]*fL[4]+alphaR[3]*fL[3]+1.732050807568877*alphaR[1]*fL[2]+alphaR[0]*fL[1]+fL[0]*alphaR[1]); 
  incr[2] = -0.25*(3.0*(alphaR[3]*fL[7]+alphaR[2]*fL[6])+1.732050807568877*alphaR[3]*fL[5]+3.0*alphaR[1]*fL[4]+1.732050807568877*alphaR[2]*fL[3]+3.0*alphaR[0]*fL[2]+1.732050807568877*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[3] = 0.25*(1.732050807568877*(alphaR[1]*fL[7]+alphaR[0]*fL[6])+alphaR[1]*fL[5]+1.732050807568877*alphaR[3]*fL[4]+alphaR[0]*fL[3]+fL[1]*alphaR[3]+alphaR[2]*(1.732050807568877*fL[2]+fL[0])); 
  incr[4] = -0.25*(3.0*(alphaR[2]*fL[7]+alphaR[3]*fL[6])+1.732050807568877*alphaR[2]*fL[5]+3.0*alphaR[0]*fL[4]+1.732050807568877*alphaR[3]*fL[3]+3.0*alphaR[1]*fL[2]+1.732050807568877*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[5] = 0.25*(1.732050807568877*(alphaR[0]*fL[7]+alphaR[1]*fL[6])+alphaR[0]*fL[5]+1.732050807568877*alphaR[2]*fL[4]+alphaR[1]*fL[3]+(1.732050807568877*fL[2]+fL[0])*alphaR[3]+fL[1]*alphaR[2]); 
  incr[6] = -0.25*(3.0*(alphaR[1]*fL[7]+alphaR[0]*fL[6])+1.732050807568877*alphaR[1]*fL[5]+3.0*alphaR[3]*fL[4]+1.732050807568877*(alphaR[0]*fL[3]+fL[1]*alphaR[3])+alphaR[2]*(3.0*fL[2]+1.732050807568877*fL[0])); 
  incr[7] = -0.25*(3.0*(alphaR[0]*fL[7]+alphaR[1]*fL[6])+1.732050807568877*alphaR[0]*fL[5]+3.0*alphaR[2]*fL[4]+1.732050807568877*alphaR[1]*fL[3]+(3.0*fL[2]+1.732050807568877*fL[0])*alphaR[3]+1.732050807568877*fL[1]*alphaR[2]); 
  } else { 
  incr[0] = -0.25*(1.732050807568877*(alphaR[3]*fR[7]+alphaR[2]*fR[6])-1.0*alphaR[3]*fR[5]+1.732050807568877*alphaR[1]*fR[4]-1.0*alphaR[2]*fR[3]+1.732050807568877*alphaR[0]*fR[2]-1.0*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.25*(1.732050807568877*(alphaR[2]*fR[7]+alphaR[3]*fR[6])-1.0*alphaR[2]*fR[5]+1.732050807568877*alphaR[0]*fR[4]-1.0*alphaR[3]*fR[3]+1.732050807568877*alphaR[1]*fR[2]-1.0*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = 0.25*(3.0*(alphaR[3]*fR[7]+alphaR[2]*fR[6])-1.732050807568877*alphaR[3]*fR[5]+3.0*alphaR[1]*fR[4]-1.732050807568877*alphaR[2]*fR[3]+3.0*alphaR[0]*fR[2]-1.732050807568877*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[3] = -0.25*(1.732050807568877*(alphaR[1]*fR[7]+alphaR[0]*fR[6])-1.0*alphaR[1]*fR[5]+1.732050807568877*alphaR[3]*fR[4]-1.0*(alphaR[0]*fR[3]+fR[1]*alphaR[3])+alphaR[2]*(1.732050807568877*fR[2]-1.0*fR[0])); 
  incr[4] = 0.25*(3.0*(alphaR[2]*fR[7]+alphaR[3]*fR[6])-1.732050807568877*alphaR[2]*fR[5]+3.0*alphaR[0]*fR[4]-1.732050807568877*alphaR[3]*fR[3]+3.0*alphaR[1]*fR[2]-1.732050807568877*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[5] = -0.25*(1.732050807568877*(alphaR[0]*fR[7]+alphaR[1]*fR[6])-1.0*alphaR[0]*fR[5]+1.732050807568877*alphaR[2]*fR[4]-1.0*alphaR[1]*fR[3]+(1.732050807568877*fR[2]-1.0*fR[0])*alphaR[3]-1.0*fR[1]*alphaR[2]); 
  incr[6] = 0.25*(3.0*(alphaR[1]*fR[7]+alphaR[0]*fR[6])-1.732050807568877*alphaR[1]*fR[5]+3.0*alphaR[3]*fR[4]-1.732050807568877*(alphaR[0]*fR[3]+fR[1]*alphaR[3])+alphaR[2]*(3.0*fR[2]-1.732050807568877*fR[0])); 
  incr[7] = 0.25*(3.0*(alphaR[0]*fR[7]+alphaR[1]*fR[6])-1.732050807568877*alphaR[0]*fR[5]+3.0*alphaR[2]*fR[4]-1.732050807568877*alphaR[1]*fR[3]+(3.0*fR[2]-1.732050807568877*fR[0])*alphaR[3]-1.732050807568877*fR[1]*alphaR[2]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[4];
  alphaOrdR = 0.5*alphaR[3]-0.5*(alphaR[2]+alphaR[1])+0.5*alphaR[0]; 
  fUpOrd[0] = 0.5*((0.6123724356957944*(fR[7]+fL[7])-0.6123724356957944*(fR[6]+fL[6])-0.3535533905932737*fR[5]+0.3535533905932737*fL[5]-0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])+0.3535533905932737*fR[1]-0.3535533905932737*(fL[1]+fR[0])+0.3535533905932737*fL[0])*sgn(alphaOrdR)-0.6123724356957944*fR[7]+0.6123724356957944*(fL[7]+fR[6])-0.6123724356957944*fL[6]+0.3535533905932737*(fR[5]+fL[5])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]-0.3535533905932737*(fR[1]+fL[1])+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[1]+alphaR[0])-0.5*(alphaR[3]+alphaR[2]); 
  fUpOrd[1] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]+fR[6]+fL[6]))+0.3535533905932737*fR[5]-0.3535533905932737*fL[5]+0.6123724356957944*(fR[4]+fL[4])+0.3535533905932737*fR[3]-0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])-0.3535533905932737*(fR[1]+fR[0])+0.3535533905932737*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.6123724356957944*(fR[7]+fR[6])-0.6123724356957944*(fL[7]+fL[6])-0.3535533905932737*(fR[5]+fL[5])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]-0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]+0.3535533905932737*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.5*alphaR[3])+0.5*alphaR[2]-0.5*alphaR[1]+0.5*alphaR[0]; 
  fUpOrd[2] = 0.5*(((-0.6123724356957944*(fR[7]+fL[7]))+0.6123724356957944*(fR[6]+fL[6])+0.3535533905932737*fR[5]-0.3535533905932737*fL[5]-0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])+0.3535533905932737*fR[1]-0.3535533905932737*(fL[1]+fR[0])+0.3535533905932737*fL[0])*sgn(alphaOrdR)+0.6123724356957944*fR[7]-0.6123724356957944*(fL[7]+fR[6])+0.6123724356957944*fL[6]-0.3535533905932737*(fR[5]+fL[5])+0.6123724356957944*fR[4]-0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]-0.3535533905932737*(fR[1]+fL[1])+0.3535533905932737*(fR[0]+fL[0])); 
  alphaOrdR = 0.5*(alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0]); 
  fUpOrd[3] = 0.5*((0.6123724356957944*(fR[7]+fL[7]+fR[6]+fL[6])-0.3535533905932737*fR[5]+0.3535533905932737*fL[5]+0.6123724356957944*(fR[4]+fL[4])-0.3535533905932737*fR[3]+0.3535533905932737*fL[3]+0.6123724356957944*(fR[2]+fL[2])-0.3535533905932737*(fR[1]+fR[0])+0.3535533905932737*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.6123724356957944*(fR[7]+fR[6])+0.6123724356957944*(fL[7]+fL[6])+0.3535533905932737*(fR[5]+fL[5])-0.6123724356957944*fR[4]+0.6123724356957944*fL[4]+0.3535533905932737*(fR[3]+fL[3])-0.6123724356957944*fR[2]+0.6123724356957944*fL[2]+0.3535533905932737*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[4];
  fUp[0] = 0.5*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.5*(fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.5*(fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.5*(fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 

  incr[0] = 0.3535533905932737*(alphaR[3]*fUp[3]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.3535533905932737*(alphaR[2]*fUp[3]+fUp[2]*alphaR[3]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = -0.6123724356957944*(alphaR[3]*fUp[3]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[3] = 0.3535533905932737*(alphaR[1]*fUp[3]+fUp[1]*alphaR[3]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[4] = -0.6123724356957944*(alphaR[2]*fUp[3]+fUp[2]*alphaR[3]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[5] = 0.3535533905932737*(alphaR[0]*fUp[3]+fUp[0]*alphaR[3]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[6] = -0.6123724356957944*(alphaR[1]*fUp[3]+fUp[1]*alphaR[3]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[7] = -0.6123724356957944*(alphaR[0]*fUp[3]+fUp[0]*alphaR[3]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 

#endif 
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
