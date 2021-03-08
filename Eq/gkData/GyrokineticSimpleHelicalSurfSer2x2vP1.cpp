#include <GyrokineticModDecl.h>
double GyrokineticSimpleHelicalSurf2x2vSer_x_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
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

  double hamilR[16]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagR[16]; 
  BstarXdBmagR[0] = (2.0*BdriftX[0]*m_*wvparR)/q_; 
  BstarXdBmagR[3] = (1.154700538379252*BdriftX[0]*m_)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = (0.125*(bmagInv[0]*(8.485281374238571*hamilR[5]-4.898979485566357*hamilR[2])*m_*rdy2R+2.449489742783178*BstarXdBmagR[0]*hamilR[3]*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = (0.3061862178478971*BstarXdBmagR[3]*hamilR[3]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*((6.0*bmagInv[0]*hamilR[5]-3.464101615137754*bmagInv[0]*hamilR[2])*m_*rdy2R+1.732050807568877*BstarXdBmagR[0]*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.125*(alphaR[2]*(2.449489742783178*fL[6]+1.414213562373095*fL[3])+alphaR[0]*(2.449489742783178*fL[1]+1.414213562373095*fL[0])); 
  incr[1] = -0.125*(alphaR[2]*(4.242640687119286*fL[6]+2.449489742783178*fL[3])+alphaR[0]*(4.242640687119286*fL[1]+2.449489742783178*fL[0])); 
  incr[2] = 0.125*(alphaR[2]*(2.449489742783178*fL[11]+1.414213562373095*fL[7])+alphaR[0]*(2.449489742783178*fL[5]+1.414213562373095*fL[2])); 
  incr[3] = 0.125*(alphaR[0]*(2.449489742783178*fL[6]+1.414213562373095*fL[3])+(2.449489742783178*fL[1]+1.414213562373095*fL[0])*alphaR[2]); 
  incr[4] = 0.125*(alphaR[2]*(2.449489742783178*fL[13]+1.414213562373095*fL[10])+alphaR[0]*(2.449489742783178*fL[8]+1.414213562373095*fL[4])); 
  incr[5] = -0.125*(alphaR[2]*(4.242640687119286*fL[11]+2.449489742783178*fL[7])+alphaR[0]*(4.242640687119286*fL[5]+2.449489742783178*fL[2])); 
  incr[6] = -0.125*(alphaR[0]*(4.242640687119286*fL[6]+2.449489742783178*fL[3])+(4.242640687119286*fL[1]+2.449489742783178*fL[0])*alphaR[2]); 
  incr[7] = 0.125*(alphaR[0]*(2.449489742783178*fL[11]+1.414213562373095*fL[7])+alphaR[2]*(2.449489742783178*fL[5]+1.414213562373095*fL[2])); 
  incr[8] = -0.125*(alphaR[2]*(4.242640687119286*fL[13]+2.449489742783178*fL[10])+alphaR[0]*(4.242640687119286*fL[8]+2.449489742783178*fL[4])); 
  incr[9] = 0.125*(alphaR[2]*(2.449489742783178*fL[15]+1.414213562373095*fL[14])+alphaR[0]*(2.449489742783178*fL[12]+1.414213562373095*fL[9])); 
  incr[10] = 0.125*(alphaR[0]*(2.449489742783178*fL[13]+1.414213562373095*fL[10])+alphaR[2]*(2.449489742783178*fL[8]+1.414213562373095*fL[4])); 
  incr[11] = -0.125*(alphaR[0]*(4.242640687119286*fL[11]+2.449489742783178*fL[7])+alphaR[2]*(4.242640687119286*fL[5]+2.449489742783178*fL[2])); 
  incr[12] = -0.125*(alphaR[2]*(4.242640687119286*fL[15]+2.449489742783178*fL[14])+alphaR[0]*(4.242640687119286*fL[12]+2.449489742783178*fL[9])); 
  incr[13] = -0.125*(alphaR[0]*(4.242640687119286*fL[13]+2.449489742783178*fL[10])+alphaR[2]*(4.242640687119286*fL[8]+2.449489742783178*fL[4])); 
  incr[14] = 0.125*(alphaR[0]*(2.449489742783178*fL[15]+1.414213562373095*fL[14])+alphaR[2]*(2.449489742783178*fL[12]+1.414213562373095*fL[9])); 
  incr[15] = -0.125*(alphaR[0]*(4.242640687119286*fL[15]+2.449489742783178*fL[14])+alphaR[2]*(4.242640687119286*fL[12]+2.449489742783178*fL[9])); 
  } else { 
  incr[0] = -0.125*(alphaR[2]*(2.449489742783178*fR[6]-1.414213562373095*fR[3])+alphaR[0]*(2.449489742783178*fR[1]-1.414213562373095*fR[0])); 
  incr[1] = 0.125*(alphaR[2]*(4.242640687119286*fR[6]-2.449489742783178*fR[3])+alphaR[0]*(4.242640687119286*fR[1]-2.449489742783178*fR[0])); 
  incr[2] = -0.125*(alphaR[2]*(2.449489742783178*fR[11]-1.414213562373095*fR[7])+alphaR[0]*(2.449489742783178*fR[5]-1.414213562373095*fR[2])); 
  incr[3] = -0.125*(alphaR[0]*(2.449489742783178*fR[6]-1.414213562373095*fR[3])+(2.449489742783178*fR[1]-1.414213562373095*fR[0])*alphaR[2]); 
  incr[4] = -0.125*(alphaR[2]*(2.449489742783178*fR[13]-1.414213562373095*fR[10])+alphaR[0]*(2.449489742783178*fR[8]-1.414213562373095*fR[4])); 
  incr[5] = 0.125*(alphaR[2]*(4.242640687119286*fR[11]-2.449489742783178*fR[7])+alphaR[0]*(4.242640687119286*fR[5]-2.449489742783178*fR[2])); 
  incr[6] = 0.125*(alphaR[0]*(4.242640687119286*fR[6]-2.449489742783178*fR[3])+(4.242640687119286*fR[1]-2.449489742783178*fR[0])*alphaR[2]); 
  incr[7] = -0.125*(alphaR[0]*(2.449489742783178*fR[11]-1.414213562373095*fR[7])+alphaR[2]*(2.449489742783178*fR[5]-1.414213562373095*fR[2])); 
  incr[8] = 0.125*(alphaR[2]*(4.242640687119286*fR[13]-2.449489742783178*fR[10])+alphaR[0]*(4.242640687119286*fR[8]-2.449489742783178*fR[4])); 
  incr[9] = -0.125*(alphaR[2]*(2.449489742783178*fR[15]-1.414213562373095*fR[14])+alphaR[0]*(2.449489742783178*fR[12]-1.414213562373095*fR[9])); 
  incr[10] = -0.125*(alphaR[0]*(2.449489742783178*fR[13]-1.414213562373095*fR[10])+alphaR[2]*(2.449489742783178*fR[8]-1.414213562373095*fR[4])); 
  incr[11] = 0.125*(alphaR[0]*(4.242640687119286*fR[11]-2.449489742783178*fR[7])+alphaR[2]*(4.242640687119286*fR[5]-2.449489742783178*fR[2])); 
  incr[12] = 0.125*(alphaR[2]*(4.242640687119286*fR[15]-2.449489742783178*fR[14])+alphaR[0]*(4.242640687119286*fR[12]-2.449489742783178*fR[9])); 
  incr[13] = 0.125*(alphaR[0]*(4.242640687119286*fR[13]-2.449489742783178*fR[10])+alphaR[2]*(4.242640687119286*fR[8]-2.449489742783178*fR[4])); 
  incr[14] = -0.125*(alphaR[0]*(2.449489742783178*fR[15]-1.414213562373095*fR[14])+alphaR[2]*(2.449489742783178*fR[12]-1.414213562373095*fR[9])); 
  incr[15] = 0.125*(alphaR[0]*(4.242640687119286*fR[15]-2.449489742783178*fR[14])+alphaR[2]*(4.242640687119286*fR[12]-2.449489742783178*fR[9])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12]+fR[11]+fL[11])-0.25*(fR[10]+fR[9])+0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*(fR[4]+fR[3]+fR[2])-0.25*(fL[4]+fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fR[12]+fR[11])+0.4330127018922193*(fL[13]+fL[12]+fL[11])+0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4]+fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*(fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3]+fR[2])+0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fR[5])+0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])+0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3]+fR[2])-0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*(fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12]+fR[11]+fL[11])-0.25*(fR[10]+fR[9])+0.25*(fL[10]+fL[9])+0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])-0.25*(fR[4]+fR[3]+fR[2])+0.25*(fL[4]+fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fR[12]+fR[11])+0.4330127018922193*(fL[13]+fL[12]+fL[11])+0.25*(fR[10]+fL[10]+fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fR[5])+0.4330127018922193*(fL[6]+fL[5])+0.25*(fR[4]+fL[4]+fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alphaR[2]*fUp[2]+alphaR[0]*fUp[0]); 
  incr[1] = -0.4330127018922193*(alphaR[2]*fUp[2]+alphaR[0]*fUp[0]); 
  incr[2] = 0.25*(alphaR[2]*fUp[4]+alphaR[0]*fUp[1]); 
  incr[3] = 0.25*(alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[4] = 0.25*(alphaR[2]*fUp[6]+alphaR[0]*fUp[3]); 
  incr[5] = -0.4330127018922193*(alphaR[2]*fUp[4]+alphaR[0]*fUp[1]); 
  incr[6] = -0.4330127018922193*(alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[7] = 0.25*(alphaR[0]*fUp[4]+fUp[1]*alphaR[2]); 
  incr[8] = -0.4330127018922193*(alphaR[2]*fUp[6]+alphaR[0]*fUp[3]); 
  incr[9] = 0.25*(alphaR[2]*fUp[7]+alphaR[0]*fUp[5]); 
  incr[10] = 0.25*(alphaR[0]*fUp[6]+alphaR[2]*fUp[3]); 
  incr[11] = -0.4330127018922193*(alphaR[0]*fUp[4]+fUp[1]*alphaR[2]); 
  incr[12] = -0.4330127018922193*(alphaR[2]*fUp[7]+alphaR[0]*fUp[5]); 
  incr[13] = -0.4330127018922193*(alphaR[0]*fUp[6]+alphaR[2]*fUp[3]); 
  incr[14] = 0.25*(alphaR[0]*fUp[7]+alphaR[2]*fUp[5]); 
  incr[15] = -0.4330127018922193*(alphaR[0]*fUp[7]+alphaR[2]*fUp[5]); 

#endif 
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
double GyrokineticSimpleHelicalSurf2x2vSer_y_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
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

  double hamilR[16]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = (2.0*BdriftY[0]*m_*wvparR)/q_; 
  BstarYdBmagR[3] = (1.154700538379252*BdriftY[0]*m_)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = -(0.125*(bmagInv[0]*(8.485281374238571*hamilR[5]-4.898979485566357*hamilR[1])*m_*rdx2R-2.449489742783178*BstarYdBmagR[0]*hamilR[3]*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = (0.3061862178478971*BstarYdBmagR[3]*hamilR[3]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*((6.0*bmagInv[0]*hamilR[5]-3.464101615137754*bmagInv[0]*hamilR[1])*m_*rdx2R-1.732050807568877*BstarYdBmagR[0]*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.125*(alphaR[2]*(2.449489742783178*fL[7]+1.414213562373095*fL[3])+alphaR[0]*(2.449489742783178*fL[2]+1.414213562373095*fL[0])); 
  incr[1] = 0.125*(alphaR[2]*(2.449489742783178*fL[11]+1.414213562373095*fL[6])+alphaR[0]*(2.449489742783178*fL[5]+1.414213562373095*fL[1])); 
  incr[2] = -0.125*(alphaR[2]*(4.242640687119286*fL[7]+2.449489742783178*fL[3])+alphaR[0]*(4.242640687119286*fL[2]+2.449489742783178*fL[0])); 
  incr[3] = 0.125*(alphaR[0]*(2.449489742783178*fL[7]+1.414213562373095*fL[3])+alphaR[2]*(2.449489742783178*fL[2]+1.414213562373095*fL[0])); 
  incr[4] = 0.125*(alphaR[2]*(2.449489742783178*fL[14]+1.414213562373095*fL[10])+alphaR[0]*(2.449489742783178*fL[9]+1.414213562373095*fL[4])); 
  incr[5] = -0.125*(alphaR[2]*(4.242640687119286*fL[11]+2.449489742783178*fL[6])+alphaR[0]*(4.242640687119286*fL[5]+2.449489742783178*fL[1])); 
  incr[6] = 0.125*(alphaR[0]*(2.449489742783178*fL[11]+1.414213562373095*fL[6])+alphaR[2]*(2.449489742783178*fL[5]+1.414213562373095*fL[1])); 
  incr[7] = -0.125*(alphaR[0]*(4.242640687119286*fL[7]+2.449489742783178*fL[3])+alphaR[2]*(4.242640687119286*fL[2]+2.449489742783178*fL[0])); 
  incr[8] = 0.125*(alphaR[2]*(2.449489742783178*fL[15]+1.414213562373095*fL[13])+alphaR[0]*(2.449489742783178*fL[12]+1.414213562373095*fL[8])); 
  incr[9] = -0.125*(alphaR[2]*(4.242640687119286*fL[14]+2.449489742783178*fL[10])+alphaR[0]*(4.242640687119286*fL[9]+2.449489742783178*fL[4])); 
  incr[10] = 0.125*(alphaR[0]*(2.449489742783178*fL[14]+1.414213562373095*fL[10])+alphaR[2]*(2.449489742783178*fL[9]+1.414213562373095*fL[4])); 
  incr[11] = -0.125*(alphaR[0]*(4.242640687119286*fL[11]+2.449489742783178*fL[6])+alphaR[2]*(4.242640687119286*fL[5]+2.449489742783178*fL[1])); 
  incr[12] = -0.125*(alphaR[2]*(4.242640687119286*fL[15]+2.449489742783178*fL[13])+alphaR[0]*(4.242640687119286*fL[12]+2.449489742783178*fL[8])); 
  incr[13] = 0.125*(alphaR[0]*(2.449489742783178*fL[15]+1.414213562373095*fL[13])+alphaR[2]*(2.449489742783178*fL[12]+1.414213562373095*fL[8])); 
  incr[14] = -0.125*(alphaR[0]*(4.242640687119286*fL[14]+2.449489742783178*fL[10])+alphaR[2]*(4.242640687119286*fL[9]+2.449489742783178*fL[4])); 
  incr[15] = -0.125*(alphaR[0]*(4.242640687119286*fL[15]+2.449489742783178*fL[13])+alphaR[2]*(4.242640687119286*fL[12]+2.449489742783178*fL[8])); 
  } else { 
  incr[0] = -0.125*(alphaR[2]*(2.449489742783178*fR[7]-1.414213562373095*fR[3])+alphaR[0]*(2.449489742783178*fR[2]-1.414213562373095*fR[0])); 
  incr[1] = -0.125*(alphaR[2]*(2.449489742783178*fR[11]-1.414213562373095*fR[6])+alphaR[0]*(2.449489742783178*fR[5]-1.414213562373095*fR[1])); 
  incr[2] = 0.125*(alphaR[2]*(4.242640687119286*fR[7]-2.449489742783178*fR[3])+alphaR[0]*(4.242640687119286*fR[2]-2.449489742783178*fR[0])); 
  incr[3] = -0.125*(alphaR[0]*(2.449489742783178*fR[7]-1.414213562373095*fR[3])+alphaR[2]*(2.449489742783178*fR[2]-1.414213562373095*fR[0])); 
  incr[4] = -0.125*(alphaR[2]*(2.449489742783178*fR[14]-1.414213562373095*fR[10])+alphaR[0]*(2.449489742783178*fR[9]-1.414213562373095*fR[4])); 
  incr[5] = 0.125*(alphaR[2]*(4.242640687119286*fR[11]-2.449489742783178*fR[6])+alphaR[0]*(4.242640687119286*fR[5]-2.449489742783178*fR[1])); 
  incr[6] = -0.125*(alphaR[0]*(2.449489742783178*fR[11]-1.414213562373095*fR[6])+alphaR[2]*(2.449489742783178*fR[5]-1.414213562373095*fR[1])); 
  incr[7] = 0.125*(alphaR[0]*(4.242640687119286*fR[7]-2.449489742783178*fR[3])+alphaR[2]*(4.242640687119286*fR[2]-2.449489742783178*fR[0])); 
  incr[8] = -0.125*(alphaR[2]*(2.449489742783178*fR[15]-1.414213562373095*fR[13])+alphaR[0]*(2.449489742783178*fR[12]-1.414213562373095*fR[8])); 
  incr[9] = 0.125*(alphaR[2]*(4.242640687119286*fR[14]-2.449489742783178*fR[10])+alphaR[0]*(4.242640687119286*fR[9]-2.449489742783178*fR[4])); 
  incr[10] = -0.125*(alphaR[0]*(2.449489742783178*fR[14]-1.414213562373095*fR[10])+alphaR[2]*(2.449489742783178*fR[9]-1.414213562373095*fR[4])); 
  incr[11] = 0.125*(alphaR[0]*(4.242640687119286*fR[11]-2.449489742783178*fR[6])+alphaR[2]*(4.242640687119286*fR[5]-2.449489742783178*fR[1])); 
  incr[12] = 0.125*(alphaR[2]*(4.242640687119286*fR[15]-2.449489742783178*fR[13])+alphaR[0]*(4.242640687119286*fR[12]-2.449489742783178*fR[8])); 
  incr[13] = -0.125*(alphaR[0]*(2.449489742783178*fR[15]-1.414213562373095*fR[13])+alphaR[2]*(2.449489742783178*fR[12]-1.414213562373095*fR[8])); 
  incr[14] = 0.125*(alphaR[0]*(4.242640687119286*fR[14]-2.449489742783178*fR[10])+alphaR[2]*(4.242640687119286*fR[9]-2.449489742783178*fR[4])); 
  incr[15] = 0.125*(alphaR[0]*(4.242640687119286*fR[15]-2.449489742783178*fR[13])+alphaR[2]*(4.242640687119286*fR[12]-2.449489742783178*fR[8])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])+0.25*fR[13]-0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*fL[14]-0.25*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fR[11])+0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14])+0.25*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*fL[14]+0.25*(fR[13]+fL[13])-0.4330127018922193*fR[12]+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.25*fR[13]-0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14])-0.25*(fR[13]+fL[13])+0.4330127018922193*fR[12]-0.4330127018922193*(fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*fL[14]+0.25*(fR[13]+fL[13])+0.4330127018922193*fR[12]-0.4330127018922193*(fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.25*fR[13]-0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14])-0.25*(fR[13]+fL[13])-0.4330127018922193*fR[12]+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])+0.25*fR[13]-0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*fL[14]-0.25*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14])+0.25*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fR[11])+0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alphaR[2]*fUp[2]+alphaR[0]*fUp[0]); 
  incr[1] = 0.25*(alphaR[2]*fUp[4]+alphaR[0]*fUp[1]); 
  incr[2] = -0.4330127018922193*(alphaR[2]*fUp[2]+alphaR[0]*fUp[0]); 
  incr[3] = 0.25*(alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[4] = 0.25*(alphaR[2]*fUp[6]+alphaR[0]*fUp[3]); 
  incr[5] = -0.4330127018922193*(alphaR[2]*fUp[4]+alphaR[0]*fUp[1]); 
  incr[6] = 0.25*(alphaR[0]*fUp[4]+fUp[1]*alphaR[2]); 
  incr[7] = -0.4330127018922193*(alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[8] = 0.25*(alphaR[2]*fUp[7]+alphaR[0]*fUp[5]); 
  incr[9] = -0.4330127018922193*(alphaR[2]*fUp[6]+alphaR[0]*fUp[3]); 
  incr[10] = 0.25*(alphaR[0]*fUp[6]+alphaR[2]*fUp[3]); 
  incr[11] = -0.4330127018922193*(alphaR[0]*fUp[4]+fUp[1]*alphaR[2]); 
  incr[12] = -0.4330127018922193*(alphaR[2]*fUp[7]+alphaR[0]*fUp[5]); 
  incr[13] = 0.25*(alphaR[0]*fUp[7]+alphaR[2]*fUp[5]); 
  incr[14] = -0.4330127018922193*(alphaR[0]*fUp[6]+alphaR[2]*fUp[3]); 
  incr[15] = -0.4330127018922193*(alphaR[0]*fUp[7]+alphaR[2]*fUp[5]); 

#endif 
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
double GyrokineticSimpleHelicalSurf2x2vSer_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
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

  double hamilR[16]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagR[16]; 
  BstarXdBmagR[0] = (2.0*BdriftX[0]*m_*wvparR)/q_; 
  BstarXdBmagR[3] = (1.154700538379252*BdriftX[0]*m_)/(q_*rdvpar2R); 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = (2.0*BdriftY[0]*m_*wvparR)/q_; 
  BstarYdBmagR[3] = (1.154700538379252*BdriftY[0]*m_)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = (0.1767766952966368*(hamilR[2]*(3.0*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0])*rdy2R+hamilR[1]*(3.0*BstarXdBmagR[3]-1.732050807568877*BstarXdBmagR[0])*rdx2R))/m_; 
  alphaR[1] = (0.1767766952966368*(3.0*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0])*hamilR[5]*rdy2R)/m_; 
  alphaR[2] = (0.1767766952966368*(3.0*BstarXdBmagR[3]-1.732050807568877*BstarXdBmagR[0])*hamilR[5]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*((3.0*hamilR[2]*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0]*hamilR[2])*rdy2R+(3.0*hamilR[1]*BstarXdBmagR[3]-1.732050807568877*BstarXdBmagR[0]*hamilR[1])*rdx2R))/m_; 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.125*(2.449489742783178*(alphaR[2]*fL[7]+alphaR[1]*fL[6]+alphaR[0]*fL[3])+1.414213562373095*(alphaR[2]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.125*(2.449489742783178*(alphaR[2]*fL[11]+alphaR[0]*fL[6])+1.414213562373095*alphaR[2]*fL[5]+2.449489742783178*alphaR[1]*fL[3]+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = 0.125*(2.449489742783178*(alphaR[1]*fL[11]+alphaR[0]*fL[7])+1.414213562373095*alphaR[1]*fL[5]+2.449489742783178*alphaR[2]*fL[3]+1.414213562373095*(alphaR[0]*fL[2]+fL[0]*alphaR[2])); 
  incr[3] = -0.125*(4.242640687119286*(alphaR[2]*fL[7]+alphaR[1]*fL[6]+alphaR[0]*fL[3])+2.449489742783178*(alphaR[2]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[4] = 0.125*(2.449489742783178*(alphaR[2]*fL[14]+alphaR[1]*fL[13]+alphaR[0]*fL[10])+1.414213562373095*(alphaR[2]*fL[9]+alphaR[1]*fL[8]+alphaR[0]*fL[4])); 
  incr[5] = 0.125*(2.449489742783178*(alphaR[0]*fL[11]+alphaR[1]*fL[7]+alphaR[2]*fL[6])+1.414213562373095*(alphaR[0]*fL[5]+alphaR[1]*fL[2]+fL[1]*alphaR[2])); 
  incr[6] = -0.125*(4.242640687119286*(alphaR[2]*fL[11]+alphaR[0]*fL[6])+2.449489742783178*alphaR[2]*fL[5]+4.242640687119286*alphaR[1]*fL[3]+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[7] = -0.125*(4.242640687119286*(alphaR[1]*fL[11]+alphaR[0]*fL[7])+2.449489742783178*alphaR[1]*fL[5]+4.242640687119286*alphaR[2]*fL[3]+2.449489742783178*(alphaR[0]*fL[2]+fL[0]*alphaR[2])); 
  incr[8] = 0.125*(2.449489742783178*(alphaR[2]*fL[15]+alphaR[0]*fL[13])+1.414213562373095*alphaR[2]*fL[12]+2.449489742783178*alphaR[1]*fL[10]+1.414213562373095*(alphaR[0]*fL[8]+alphaR[1]*fL[4])); 
  incr[9] = 0.125*(2.449489742783178*(alphaR[1]*fL[15]+alphaR[0]*fL[14])+1.414213562373095*alphaR[1]*fL[12]+2.449489742783178*alphaR[2]*fL[10]+1.414213562373095*(alphaR[0]*fL[9]+alphaR[2]*fL[4])); 
  incr[10] = -0.125*(4.242640687119286*(alphaR[2]*fL[14]+alphaR[1]*fL[13]+alphaR[0]*fL[10])+2.449489742783178*(alphaR[2]*fL[9]+alphaR[1]*fL[8]+alphaR[0]*fL[4])); 
  incr[11] = -0.125*(4.242640687119286*(alphaR[0]*fL[11]+alphaR[1]*fL[7]+alphaR[2]*fL[6])+2.449489742783178*(alphaR[0]*fL[5]+alphaR[1]*fL[2]+fL[1]*alphaR[2])); 
  incr[12] = 0.125*(2.449489742783178*(alphaR[0]*fL[15]+alphaR[1]*fL[14]+alphaR[2]*fL[13])+1.414213562373095*(alphaR[0]*fL[12]+alphaR[1]*fL[9]+alphaR[2]*fL[8])); 
  incr[13] = -0.125*(4.242640687119286*(alphaR[2]*fL[15]+alphaR[0]*fL[13])+2.449489742783178*alphaR[2]*fL[12]+4.242640687119286*alphaR[1]*fL[10]+2.449489742783178*(alphaR[0]*fL[8]+alphaR[1]*fL[4])); 
  incr[14] = -0.125*(4.242640687119286*(alphaR[1]*fL[15]+alphaR[0]*fL[14])+2.449489742783178*alphaR[1]*fL[12]+4.242640687119286*alphaR[2]*fL[10]+2.449489742783178*(alphaR[0]*fL[9]+alphaR[2]*fL[4])); 
  incr[15] = -0.125*(4.242640687119286*(alphaR[0]*fL[15]+alphaR[1]*fL[14]+alphaR[2]*fL[13])+2.449489742783178*(alphaR[0]*fL[12]+alphaR[1]*fL[9]+alphaR[2]*fL[8])); 
  } else { 
  incr[0] = -0.125*(2.449489742783178*(alphaR[2]*fR[7]+alphaR[1]*fR[6]+alphaR[0]*fR[3])-1.414213562373095*(alphaR[2]*fR[2]+alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.125*(2.449489742783178*(alphaR[2]*fR[11]+alphaR[0]*fR[6])-1.414213562373095*alphaR[2]*fR[5]+2.449489742783178*alphaR[1]*fR[3]-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = -0.125*(2.449489742783178*(alphaR[1]*fR[11]+alphaR[0]*fR[7])-1.414213562373095*alphaR[1]*fR[5]+2.449489742783178*alphaR[2]*fR[3]-1.414213562373095*(alphaR[0]*fR[2]+fR[0]*alphaR[2])); 
  incr[3] = 0.125*(4.242640687119286*(alphaR[2]*fR[7]+alphaR[1]*fR[6]+alphaR[0]*fR[3])-2.449489742783178*(alphaR[2]*fR[2]+alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[4] = -0.125*(2.449489742783178*(alphaR[2]*fR[14]+alphaR[1]*fR[13]+alphaR[0]*fR[10])-1.414213562373095*(alphaR[2]*fR[9]+alphaR[1]*fR[8]+alphaR[0]*fR[4])); 
  incr[5] = -0.125*(2.449489742783178*(alphaR[0]*fR[11]+alphaR[1]*fR[7]+alphaR[2]*fR[6])-1.414213562373095*(alphaR[0]*fR[5]+alphaR[1]*fR[2]+fR[1]*alphaR[2])); 
  incr[6] = 0.125*(4.242640687119286*(alphaR[2]*fR[11]+alphaR[0]*fR[6])-2.449489742783178*alphaR[2]*fR[5]+4.242640687119286*alphaR[1]*fR[3]-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[7] = 0.125*(4.242640687119286*(alphaR[1]*fR[11]+alphaR[0]*fR[7])-2.449489742783178*alphaR[1]*fR[5]+4.242640687119286*alphaR[2]*fR[3]-2.449489742783178*(alphaR[0]*fR[2]+fR[0]*alphaR[2])); 
  incr[8] = -0.125*(2.449489742783178*(alphaR[2]*fR[15]+alphaR[0]*fR[13])-1.414213562373095*alphaR[2]*fR[12]+2.449489742783178*alphaR[1]*fR[10]-1.414213562373095*(alphaR[0]*fR[8]+alphaR[1]*fR[4])); 
  incr[9] = -0.125*(2.449489742783178*(alphaR[1]*fR[15]+alphaR[0]*fR[14])-1.414213562373095*alphaR[1]*fR[12]+2.449489742783178*alphaR[2]*fR[10]-1.414213562373095*(alphaR[0]*fR[9]+alphaR[2]*fR[4])); 
  incr[10] = 0.125*(4.242640687119286*(alphaR[2]*fR[14]+alphaR[1]*fR[13]+alphaR[0]*fR[10])-2.449489742783178*(alphaR[2]*fR[9]+alphaR[1]*fR[8]+alphaR[0]*fR[4])); 
  incr[11] = 0.125*(4.242640687119286*(alphaR[0]*fR[11]+alphaR[1]*fR[7]+alphaR[2]*fR[6])-2.449489742783178*(alphaR[0]*fR[5]+alphaR[1]*fR[2]+fR[1]*alphaR[2])); 
  incr[12] = -0.125*(2.449489742783178*(alphaR[0]*fR[15]+alphaR[1]*fR[14]+alphaR[2]*fR[13])-1.414213562373095*(alphaR[0]*fR[12]+alphaR[1]*fR[9]+alphaR[2]*fR[8])); 
  incr[13] = 0.125*(4.242640687119286*(alphaR[2]*fR[15]+alphaR[0]*fR[13])-2.449489742783178*alphaR[2]*fR[12]+4.242640687119286*alphaR[1]*fR[10]-2.449489742783178*(alphaR[0]*fR[8]+alphaR[1]*fR[4])); 
  incr[14] = 0.125*(4.242640687119286*(alphaR[1]*fR[15]+alphaR[0]*fR[14])-2.449489742783178*alphaR[1]*fR[12]+4.242640687119286*alphaR[2]*fR[10]-2.449489742783178*(alphaR[0]*fR[9]+alphaR[2]*fR[4])); 
  incr[15] = 0.125*(4.242640687119286*(alphaR[0]*fR[15]+alphaR[1]*fR[14]+alphaR[2]*fR[13])-2.449489742783178*(alphaR[0]*fR[12]+alphaR[1]*fR[9]+alphaR[2]*fR[8])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*(alphaR[2]+alphaR[1]); 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])-0.25*(fR[9]+fR[8])+0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14]+fR[13])+0.4330127018922193*(fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[1]+alphaR[0])-0.3535533905932737*alphaR[2]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14]+fR[13])-0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[2]-0.3535533905932737*alphaR[1]+0.3535533905932737*alphaR[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*(fL[14]+fR[13])+0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[1]+alphaR[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]+fR[13]+fL[13]))+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])+0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*(fR[2]+fR[1]+fR[0])+0.25*(fL[2]+fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14]+fR[13])-0.4330127018922193*(fL[15]+fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9]+fR[8]+fL[8])-0.4330127018922193*(fR[7]+fR[6])+0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2]+fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*(alphaR[2]+alphaR[1]); 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*(fR[5]+fR[4])+0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14]+fR[13])-0.4330127018922193*(fL[14]+fL[13])+0.25*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fR[10])+0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[1]+alphaR[0])-0.3535533905932737*alphaR[2]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14]+fR[13])+0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[2]-0.3535533905932737*alphaR[1]+0.3535533905932737*alphaR[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*(fL[14]+fR[13])-0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[1]+alphaR[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]+fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])-0.25*(fR[9]+fR[8])+0.25*(fL[9]+fL[8])+0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*(fR[5]+fR[4])+0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])-0.25*(fR[2]+fR[1]+fR[0])+0.25*(fL[2]+fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14]+fR[13])+0.4330127018922193*(fL[15]+fL[14]+fL[13])+0.25*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fR[10])+0.4330127018922193*(fL[11]+fL[10])+0.25*(fR[9]+fL[9]+fR[8]+fL[8])-0.4330127018922193*(fR[7]+fR[6])+0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2]+fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.25*(alphaR[2]*fUp[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = 0.25*(alphaR[1]*fUp[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[3] = -0.4330127018922193*(alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[4] = 0.25*(alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+alphaR[0]*fUp[3]); 
  incr[5] = 0.25*(alphaR[0]*fUp[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[6] = -0.4330127018922193*(alphaR[2]*fUp[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[7] = -0.4330127018922193*(alphaR[1]*fUp[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[8] = 0.25*(alphaR[2]*fUp[7]+alphaR[0]*fUp[5]+alphaR[1]*fUp[3]); 
  incr[9] = 0.25*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[2]*fUp[3]); 
  incr[10] = -0.4330127018922193*(alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+alphaR[0]*fUp[3]); 
  incr[11] = -0.4330127018922193*(alphaR[0]*fUp[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[12] = 0.25*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]); 
  incr[13] = -0.4330127018922193*(alphaR[2]*fUp[7]+alphaR[0]*fUp[5]+alphaR[1]*fUp[3]); 
  incr[14] = -0.4330127018922193*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[2]*fUp[3]); 
  incr[15] = -0.4330127018922193*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]); 

#endif 
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
double GyrokineticSimpleHelicalSurf2x2vSer_x_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
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

  double hamilR[16]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 

  double BstarXdBmagR[16]; 
  BstarXdBmagR[0] = (2.0*BdriftX[0]*m_*wvparR)/q_; 
  BstarXdBmagR[1] = (2.0*BdriftX[1]*m_*wvparR)/q_; 
  BstarXdBmagR[3] = (1.154700538379252*BdriftX[0]*m_)/(q_*rdvpar2R); 
  BstarXdBmagR[6] = (1.154700538379252*BdriftX[1]*m_)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = -(0.125*(((14.69693845669907*bmagInv[1]-8.485281374238571*bmagInv[0])*hamilR[5]+(4.898979485566357*bmagInv[0]-8.485281374238571*bmagInv[1])*hamilR[2])*m_*rdy2R+(4.242640687119286*BstarXdBmagR[1]-2.449489742783178*BstarXdBmagR[0])*hamilR[3]*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = -(0.125*hamilR[3]*(4.242640687119286*BstarXdBmagR[6]-2.449489742783178*BstarXdBmagR[3])*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(((10.39230484541326*bmagInv[1]-6.0*bmagInv[0])*hamilR[5]+(3.464101615137754*bmagInv[0]-6.0*bmagInv[1])*hamilR[2])*m_*rdy2R+(3.0*BstarXdBmagR[1]-1.732050807568877*BstarXdBmagR[0])*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.125*(alphaR[2]*(2.449489742783178*fL[6]+1.414213562373095*fL[3])+alphaR[0]*(2.449489742783178*fL[1]+1.414213562373095*fL[0])); 
  incr[1] = -0.125*(alphaR[2]*(4.242640687119286*fL[6]+2.449489742783178*fL[3])+alphaR[0]*(4.242640687119286*fL[1]+2.449489742783178*fL[0])); 
  incr[2] = 0.125*(alphaR[2]*(2.449489742783178*fL[11]+1.414213562373095*fL[7])+alphaR[0]*(2.449489742783178*fL[5]+1.414213562373095*fL[2])); 
  incr[3] = 0.125*(alphaR[0]*(2.449489742783178*fL[6]+1.414213562373095*fL[3])+(2.449489742783178*fL[1]+1.414213562373095*fL[0])*alphaR[2]); 
  incr[4] = 0.125*(alphaR[2]*(2.449489742783178*fL[13]+1.414213562373095*fL[10])+alphaR[0]*(2.449489742783178*fL[8]+1.414213562373095*fL[4])); 
  incr[5] = -0.125*(alphaR[2]*(4.242640687119286*fL[11]+2.449489742783178*fL[7])+alphaR[0]*(4.242640687119286*fL[5]+2.449489742783178*fL[2])); 
  incr[6] = -0.125*(alphaR[0]*(4.242640687119286*fL[6]+2.449489742783178*fL[3])+(4.242640687119286*fL[1]+2.449489742783178*fL[0])*alphaR[2]); 
  incr[7] = 0.125*(alphaR[0]*(2.449489742783178*fL[11]+1.414213562373095*fL[7])+alphaR[2]*(2.449489742783178*fL[5]+1.414213562373095*fL[2])); 
  incr[8] = -0.125*(alphaR[2]*(4.242640687119286*fL[13]+2.449489742783178*fL[10])+alphaR[0]*(4.242640687119286*fL[8]+2.449489742783178*fL[4])); 
  incr[9] = 0.125*(alphaR[2]*(2.449489742783178*fL[15]+1.414213562373095*fL[14])+alphaR[0]*(2.449489742783178*fL[12]+1.414213562373095*fL[9])); 
  incr[10] = 0.125*(alphaR[0]*(2.449489742783178*fL[13]+1.414213562373095*fL[10])+alphaR[2]*(2.449489742783178*fL[8]+1.414213562373095*fL[4])); 
  incr[11] = -0.125*(alphaR[0]*(4.242640687119286*fL[11]+2.449489742783178*fL[7])+alphaR[2]*(4.242640687119286*fL[5]+2.449489742783178*fL[2])); 
  incr[12] = -0.125*(alphaR[2]*(4.242640687119286*fL[15]+2.449489742783178*fL[14])+alphaR[0]*(4.242640687119286*fL[12]+2.449489742783178*fL[9])); 
  incr[13] = -0.125*(alphaR[0]*(4.242640687119286*fL[13]+2.449489742783178*fL[10])+alphaR[2]*(4.242640687119286*fL[8]+2.449489742783178*fL[4])); 
  incr[14] = 0.125*(alphaR[0]*(2.449489742783178*fL[15]+1.414213562373095*fL[14])+alphaR[2]*(2.449489742783178*fL[12]+1.414213562373095*fL[9])); 
  incr[15] = -0.125*(alphaR[0]*(4.242640687119286*fL[15]+2.449489742783178*fL[14])+alphaR[2]*(4.242640687119286*fL[12]+2.449489742783178*fL[9])); 
  } else { 
  incr[0] = -0.125*(alphaR[2]*(2.449489742783178*fR[6]-1.414213562373095*fR[3])+alphaR[0]*(2.449489742783178*fR[1]-1.414213562373095*fR[0])); 
  incr[1] = 0.125*(alphaR[2]*(4.242640687119286*fR[6]-2.449489742783178*fR[3])+alphaR[0]*(4.242640687119286*fR[1]-2.449489742783178*fR[0])); 
  incr[2] = -0.125*(alphaR[2]*(2.449489742783178*fR[11]-1.414213562373095*fR[7])+alphaR[0]*(2.449489742783178*fR[5]-1.414213562373095*fR[2])); 
  incr[3] = -0.125*(alphaR[0]*(2.449489742783178*fR[6]-1.414213562373095*fR[3])+(2.449489742783178*fR[1]-1.414213562373095*fR[0])*alphaR[2]); 
  incr[4] = -0.125*(alphaR[2]*(2.449489742783178*fR[13]-1.414213562373095*fR[10])+alphaR[0]*(2.449489742783178*fR[8]-1.414213562373095*fR[4])); 
  incr[5] = 0.125*(alphaR[2]*(4.242640687119286*fR[11]-2.449489742783178*fR[7])+alphaR[0]*(4.242640687119286*fR[5]-2.449489742783178*fR[2])); 
  incr[6] = 0.125*(alphaR[0]*(4.242640687119286*fR[6]-2.449489742783178*fR[3])+(4.242640687119286*fR[1]-2.449489742783178*fR[0])*alphaR[2]); 
  incr[7] = -0.125*(alphaR[0]*(2.449489742783178*fR[11]-1.414213562373095*fR[7])+alphaR[2]*(2.449489742783178*fR[5]-1.414213562373095*fR[2])); 
  incr[8] = 0.125*(alphaR[2]*(4.242640687119286*fR[13]-2.449489742783178*fR[10])+alphaR[0]*(4.242640687119286*fR[8]-2.449489742783178*fR[4])); 
  incr[9] = -0.125*(alphaR[2]*(2.449489742783178*fR[15]-1.414213562373095*fR[14])+alphaR[0]*(2.449489742783178*fR[12]-1.414213562373095*fR[9])); 
  incr[10] = -0.125*(alphaR[0]*(2.449489742783178*fR[13]-1.414213562373095*fR[10])+alphaR[2]*(2.449489742783178*fR[8]-1.414213562373095*fR[4])); 
  incr[11] = 0.125*(alphaR[0]*(4.242640687119286*fR[11]-2.449489742783178*fR[7])+alphaR[2]*(4.242640687119286*fR[5]-2.449489742783178*fR[2])); 
  incr[12] = 0.125*(alphaR[2]*(4.242640687119286*fR[15]-2.449489742783178*fR[14])+alphaR[0]*(4.242640687119286*fR[12]-2.449489742783178*fR[9])); 
  incr[13] = 0.125*(alphaR[0]*(4.242640687119286*fR[13]-2.449489742783178*fR[10])+alphaR[2]*(4.242640687119286*fR[8]-2.449489742783178*fR[4])); 
  incr[14] = -0.125*(alphaR[0]*(2.449489742783178*fR[15]-1.414213562373095*fR[14])+alphaR[2]*(2.449489742783178*fR[12]-1.414213562373095*fR[9])); 
  incr[15] = 0.125*(alphaR[0]*(4.242640687119286*fR[15]-2.449489742783178*fR[14])+alphaR[2]*(4.242640687119286*fR[12]-2.449489742783178*fR[9])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12]+fR[11]+fL[11])-0.25*(fR[10]+fR[9])+0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*(fR[4]+fR[3]+fR[2])-0.25*(fL[4]+fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fR[12]+fR[11])+0.4330127018922193*(fL[13]+fL[12]+fL[11])+0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4]+fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*(fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3]+fR[2])+0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fR[5])+0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])+0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3]+fR[2])-0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]-0.3535533905932737*alphaR[2]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*(fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[2]+alphaR[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12]+fR[11]+fL[11])-0.25*(fR[10]+fR[9])+0.25*(fL[10]+fL[9])+0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])-0.25*(fR[4]+fR[3]+fR[2])+0.25*(fL[4]+fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fR[12]+fR[11])+0.4330127018922193*(fL[13]+fL[12]+fL[11])+0.25*(fR[10]+fL[10]+fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fR[5])+0.4330127018922193*(fL[6]+fL[5])+0.25*(fR[4]+fL[4]+fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alphaR[2]*fUp[2]+alphaR[0]*fUp[0]); 
  incr[1] = -0.4330127018922193*(alphaR[2]*fUp[2]+alphaR[0]*fUp[0]); 
  incr[2] = 0.25*(alphaR[2]*fUp[4]+alphaR[0]*fUp[1]); 
  incr[3] = 0.25*(alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[4] = 0.25*(alphaR[2]*fUp[6]+alphaR[0]*fUp[3]); 
  incr[5] = -0.4330127018922193*(alphaR[2]*fUp[4]+alphaR[0]*fUp[1]); 
  incr[6] = -0.4330127018922193*(alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[7] = 0.25*(alphaR[0]*fUp[4]+fUp[1]*alphaR[2]); 
  incr[8] = -0.4330127018922193*(alphaR[2]*fUp[6]+alphaR[0]*fUp[3]); 
  incr[9] = 0.25*(alphaR[2]*fUp[7]+alphaR[0]*fUp[5]); 
  incr[10] = 0.25*(alphaR[0]*fUp[6]+alphaR[2]*fUp[3]); 
  incr[11] = -0.4330127018922193*(alphaR[0]*fUp[4]+fUp[1]*alphaR[2]); 
  incr[12] = -0.4330127018922193*(alphaR[2]*fUp[7]+alphaR[0]*fUp[5]); 
  incr[13] = -0.4330127018922193*(alphaR[0]*fUp[6]+alphaR[2]*fUp[3]); 
  incr[14] = 0.25*(alphaR[0]*fUp[7]+alphaR[2]*fUp[5]); 
  incr[15] = -0.4330127018922193*(alphaR[0]*fUp[7]+alphaR[2]*fUp[5]); 

#endif 
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
double GyrokineticSimpleHelicalSurf2x2vSer_y_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
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

  double hamilR[16]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = (2.0*BdriftY[0]*m_*wvparR)/q_; 
  BstarYdBmagR[1] = (2.0*BdriftY[1]*m_*wvparR)/q_; 
  BstarYdBmagR[3] = (1.154700538379252*BdriftY[0]*m_)/(q_*rdvpar2R); 
  BstarYdBmagR[6] = (1.154700538379252*BdriftY[1]*m_)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = -(0.1767766952966368*(bmagInv[0]*(6.0*hamilR[5]-3.464101615137754*hamilR[1])*m_*rdx2R-1.732050807568877*BstarYdBmagR[0]*hamilR[3]*q_*rdvpar2R))/(m_*q_); 
  alphaR[1] = -(0.1767766952966368*(bmagInv[1]*(6.0*hamilR[5]-3.464101615137754*hamilR[1])*m_*rdx2R-1.732050807568877*BstarYdBmagR[1]*hamilR[3]*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = (0.3061862178478971*BstarYdBmagR[3]*hamilR[3]*rdvpar2R)/m_; 
  alphaR[3] = (0.6123724356957944*bmagInv[0]*hamilR[8]*rdx2R)/q_; 
  alphaR[4] = (0.3061862178478971*hamilR[3]*BstarYdBmagR[6]*rdvpar2R)/m_; 
  alphaR[5] = (0.6123724356957944*bmagInv[1]*hamilR[8]*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*((6.0*bmagInv[0]*hamilR[5]-3.464101615137754*bmagInv[0]*hamilR[1])*m_*rdx2R-1.732050807568877*BstarYdBmagR[0]*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.125*(2.449489742783178*(alphaR[5]*fL[12]+alphaR[4]*fL[11]+alphaR[3]*fL[9])+1.414213562373095*alphaR[5]*fL[8]+2.449489742783178*alphaR[2]*fL[7]+1.414213562373095*alphaR[4]*fL[6]+2.449489742783178*alphaR[1]*fL[5]+1.414213562373095*(alphaR[3]*fL[4]+alphaR[2]*fL[3])+2.449489742783178*alphaR[0]*fL[2]+1.414213562373095*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.125*(2.449489742783178*(alphaR[3]*fL[12]+alphaR[2]*fL[11]+alphaR[5]*fL[9])+1.414213562373095*alphaR[3]*fL[8]+2.449489742783178*alphaR[4]*fL[7]+1.414213562373095*alphaR[2]*fL[6]+2.449489742783178*alphaR[0]*fL[5]+1.414213562373095*(fL[4]*alphaR[5]+fL[3]*alphaR[4])+2.449489742783178*alphaR[1]*fL[2]+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = -0.125*(4.242640687119286*(alphaR[5]*fL[12]+alphaR[4]*fL[11]+alphaR[3]*fL[9])+2.449489742783178*alphaR[5]*fL[8]+4.242640687119286*alphaR[2]*fL[7]+2.449489742783178*alphaR[4]*fL[6]+4.242640687119286*alphaR[1]*fL[5]+2.449489742783178*(alphaR[3]*fL[4]+alphaR[2]*fL[3])+4.242640687119286*alphaR[0]*fL[2]+2.449489742783178*(alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[3] = 0.125*(2.449489742783178*(alphaR[5]*fL[15]+alphaR[3]*fL[14])+1.414213562373095*alphaR[5]*fL[13]+2.449489742783178*alphaR[1]*fL[11]+1.414213562373095*alphaR[3]*fL[10]+2.449489742783178*alphaR[0]*fL[7]+1.414213562373095*alphaR[1]*fL[6]+2.449489742783178*alphaR[4]*fL[5]+1.414213562373095*(fL[1]*alphaR[4]+alphaR[0]*fL[3])+alphaR[2]*(2.449489742783178*fL[2]+1.414213562373095*fL[0])); 
  incr[4] = 0.125*(2.449489742783178*(alphaR[4]*fL[15]+alphaR[2]*fL[14])+1.414213562373095*alphaR[4]*fL[13]+2.449489742783178*alphaR[1]*fL[12]+1.414213562373095*alphaR[2]*fL[10]+2.449489742783178*alphaR[0]*fL[9]+1.414213562373095*alphaR[1]*fL[8]+2.449489742783178*alphaR[5]*fL[5]+1.414213562373095*(fL[1]*alphaR[5]+alphaR[0]*fL[4])+(2.449489742783178*fL[2]+1.414213562373095*fL[0])*alphaR[3]); 
  incr[5] = -0.125*(4.242640687119286*(alphaR[3]*fL[12]+alphaR[2]*fL[11]+alphaR[5]*fL[9])+2.449489742783178*alphaR[3]*fL[8]+4.242640687119286*alphaR[4]*fL[7]+2.449489742783178*alphaR[2]*fL[6]+4.242640687119286*alphaR[0]*fL[5]+2.449489742783178*(fL[4]*alphaR[5]+fL[3]*alphaR[4])+4.242640687119286*alphaR[1]*fL[2]+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[6] = 0.125*(2.449489742783178*(alphaR[3]*fL[15]+alphaR[5]*fL[14])+1.414213562373095*alphaR[3]*fL[13]+2.449489742783178*alphaR[0]*fL[11]+1.414213562373095*alphaR[5]*fL[10]+2.449489742783178*alphaR[1]*fL[7]+1.414213562373095*alphaR[0]*fL[6]+2.449489742783178*alphaR[2]*fL[5]+(2.449489742783178*fL[2]+1.414213562373095*fL[0])*alphaR[4]+1.414213562373095*(alphaR[1]*fL[3]+fL[1]*alphaR[2])); 
  incr[7] = -0.125*(4.242640687119286*(alphaR[5]*fL[15]+alphaR[3]*fL[14])+2.449489742783178*alphaR[5]*fL[13]+4.242640687119286*alphaR[1]*fL[11]+2.449489742783178*alphaR[3]*fL[10]+4.242640687119286*alphaR[0]*fL[7]+2.449489742783178*alphaR[1]*fL[6]+4.242640687119286*alphaR[4]*fL[5]+2.449489742783178*(fL[1]*alphaR[4]+alphaR[0]*fL[3])+alphaR[2]*(4.242640687119286*fL[2]+2.449489742783178*fL[0])); 
  incr[8] = 0.125*(2.449489742783178*(alphaR[2]*fL[15]+alphaR[4]*fL[14])+1.414213562373095*alphaR[2]*fL[13]+2.449489742783178*alphaR[0]*fL[12]+1.414213562373095*alphaR[4]*fL[10]+2.449489742783178*alphaR[1]*fL[9]+1.414213562373095*alphaR[0]*fL[8]+2.449489742783178*alphaR[3]*fL[5]+(2.449489742783178*fL[2]+1.414213562373095*fL[0])*alphaR[5]+1.414213562373095*(alphaR[1]*fL[4]+fL[1]*alphaR[3])); 
  incr[9] = -0.125*(4.242640687119286*(alphaR[4]*fL[15]+alphaR[2]*fL[14])+2.449489742783178*alphaR[4]*fL[13]+4.242640687119286*alphaR[1]*fL[12]+2.449489742783178*alphaR[2]*fL[10]+4.242640687119286*alphaR[0]*fL[9]+2.449489742783178*alphaR[1]*fL[8]+4.242640687119286*alphaR[5]*fL[5]+2.449489742783178*(fL[1]*alphaR[5]+alphaR[0]*fL[4])+(4.242640687119286*fL[2]+2.449489742783178*fL[0])*alphaR[3]); 
  incr[10] = 0.125*(2.449489742783178*(alphaR[1]*fL[15]+alphaR[0]*fL[14])+1.414213562373095*alphaR[1]*fL[13]+2.449489742783178*(alphaR[4]*fL[12]+alphaR[5]*fL[11])+1.414213562373095*alphaR[0]*fL[10]+2.449489742783178*alphaR[2]*fL[9]+1.414213562373095*alphaR[4]*fL[8]+2.449489742783178*alphaR[3]*fL[7]+1.414213562373095*(alphaR[5]*fL[6]+alphaR[2]*fL[4]+alphaR[3]*fL[3])); 
  incr[11] = -0.125*(4.242640687119286*(alphaR[3]*fL[15]+alphaR[5]*fL[14])+2.449489742783178*alphaR[3]*fL[13]+4.242640687119286*alphaR[0]*fL[11]+2.449489742783178*alphaR[5]*fL[10]+4.242640687119286*alphaR[1]*fL[7]+2.449489742783178*alphaR[0]*fL[6]+4.242640687119286*alphaR[2]*fL[5]+(4.242640687119286*fL[2]+2.449489742783178*fL[0])*alphaR[4]+2.449489742783178*(alphaR[1]*fL[3]+fL[1]*alphaR[2])); 
  incr[12] = -0.125*(4.242640687119286*(alphaR[2]*fL[15]+alphaR[4]*fL[14])+2.449489742783178*alphaR[2]*fL[13]+4.242640687119286*alphaR[0]*fL[12]+2.449489742783178*alphaR[4]*fL[10]+4.242640687119286*alphaR[1]*fL[9]+2.449489742783178*alphaR[0]*fL[8]+4.242640687119286*alphaR[3]*fL[5]+(4.242640687119286*fL[2]+2.449489742783178*fL[0])*alphaR[5]+2.449489742783178*(alphaR[1]*fL[4]+fL[1]*alphaR[3])); 
  incr[13] = 0.125*(2.449489742783178*(alphaR[0]*fL[15]+alphaR[1]*fL[14])+1.414213562373095*alphaR[0]*fL[13]+2.449489742783178*(alphaR[2]*fL[12]+alphaR[3]*fL[11])+1.414213562373095*alphaR[1]*fL[10]+2.449489742783178*alphaR[4]*fL[9]+1.414213562373095*alphaR[2]*fL[8]+2.449489742783178*alphaR[5]*fL[7]+1.414213562373095*(alphaR[3]*fL[6]+fL[3]*alphaR[5]+alphaR[4]*fL[4])); 
  incr[14] = -0.125*(4.242640687119286*(alphaR[1]*fL[15]+alphaR[0]*fL[14])+2.449489742783178*alphaR[1]*fL[13]+4.242640687119286*(alphaR[4]*fL[12]+alphaR[5]*fL[11])+2.449489742783178*alphaR[0]*fL[10]+4.242640687119286*alphaR[2]*fL[9]+2.449489742783178*alphaR[4]*fL[8]+4.242640687119286*alphaR[3]*fL[7]+2.449489742783178*(alphaR[5]*fL[6]+alphaR[2]*fL[4]+alphaR[3]*fL[3])); 
  incr[15] = -0.125*(4.242640687119286*(alphaR[0]*fL[15]+alphaR[1]*fL[14])+2.449489742783178*alphaR[0]*fL[13]+4.242640687119286*(alphaR[2]*fL[12]+alphaR[3]*fL[11])+2.449489742783178*alphaR[1]*fL[10]+4.242640687119286*alphaR[4]*fL[9]+2.449489742783178*alphaR[2]*fL[8]+4.242640687119286*alphaR[5]*fL[7]+2.449489742783178*(alphaR[3]*fL[6]+fL[3]*alphaR[5]+alphaR[4]*fL[4])); 
  } else { 
  incr[0] = -0.125*(2.449489742783178*(alphaR[5]*fR[12]+alphaR[4]*fR[11]+alphaR[3]*fR[9])-1.414213562373095*alphaR[5]*fR[8]+2.449489742783178*alphaR[2]*fR[7]-1.414213562373095*alphaR[4]*fR[6]+2.449489742783178*alphaR[1]*fR[5]-1.414213562373095*(alphaR[3]*fR[4]+alphaR[2]*fR[3])+2.449489742783178*alphaR[0]*fR[2]-1.414213562373095*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.125*(2.449489742783178*(alphaR[3]*fR[12]+alphaR[2]*fR[11]+alphaR[5]*fR[9])-1.414213562373095*alphaR[3]*fR[8]+2.449489742783178*alphaR[4]*fR[7]-1.414213562373095*alphaR[2]*fR[6]+2.449489742783178*alphaR[0]*fR[5]-1.414213562373095*(fR[4]*alphaR[5]+fR[3]*alphaR[4])+2.449489742783178*alphaR[1]*fR[2]-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = 0.125*(4.242640687119286*(alphaR[5]*fR[12]+alphaR[4]*fR[11]+alphaR[3]*fR[9])-2.449489742783178*alphaR[5]*fR[8]+4.242640687119286*alphaR[2]*fR[7]-2.449489742783178*alphaR[4]*fR[6]+4.242640687119286*alphaR[1]*fR[5]-2.449489742783178*(alphaR[3]*fR[4]+alphaR[2]*fR[3])+4.242640687119286*alphaR[0]*fR[2]-2.449489742783178*(alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[3] = -0.125*(2.449489742783178*(alphaR[5]*fR[15]+alphaR[3]*fR[14])-1.414213562373095*alphaR[5]*fR[13]+2.449489742783178*alphaR[1]*fR[11]-1.414213562373095*alphaR[3]*fR[10]+2.449489742783178*alphaR[0]*fR[7]-1.414213562373095*alphaR[1]*fR[6]+2.449489742783178*alphaR[4]*fR[5]-1.414213562373095*(fR[1]*alphaR[4]+alphaR[0]*fR[3])+alphaR[2]*(2.449489742783178*fR[2]-1.414213562373095*fR[0])); 
  incr[4] = -0.125*(2.449489742783178*(alphaR[4]*fR[15]+alphaR[2]*fR[14])-1.414213562373095*alphaR[4]*fR[13]+2.449489742783178*alphaR[1]*fR[12]-1.414213562373095*alphaR[2]*fR[10]+2.449489742783178*alphaR[0]*fR[9]-1.414213562373095*alphaR[1]*fR[8]+2.449489742783178*alphaR[5]*fR[5]-1.414213562373095*(fR[1]*alphaR[5]+alphaR[0]*fR[4])+(2.449489742783178*fR[2]-1.414213562373095*fR[0])*alphaR[3]); 
  incr[5] = 0.125*(4.242640687119286*(alphaR[3]*fR[12]+alphaR[2]*fR[11]+alphaR[5]*fR[9])-2.449489742783178*alphaR[3]*fR[8]+4.242640687119286*alphaR[4]*fR[7]-2.449489742783178*alphaR[2]*fR[6]+4.242640687119286*alphaR[0]*fR[5]-2.449489742783178*(fR[4]*alphaR[5]+fR[3]*alphaR[4])+4.242640687119286*alphaR[1]*fR[2]-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[6] = -0.125*(2.449489742783178*(alphaR[3]*fR[15]+alphaR[5]*fR[14])-1.414213562373095*alphaR[3]*fR[13]+2.449489742783178*alphaR[0]*fR[11]-1.414213562373095*alphaR[5]*fR[10]+2.449489742783178*alphaR[1]*fR[7]-1.414213562373095*alphaR[0]*fR[6]+2.449489742783178*alphaR[2]*fR[5]+(2.449489742783178*fR[2]-1.414213562373095*fR[0])*alphaR[4]-1.414213562373095*(alphaR[1]*fR[3]+fR[1]*alphaR[2])); 
  incr[7] = 0.125*(4.242640687119286*(alphaR[5]*fR[15]+alphaR[3]*fR[14])-2.449489742783178*alphaR[5]*fR[13]+4.242640687119286*alphaR[1]*fR[11]-2.449489742783178*alphaR[3]*fR[10]+4.242640687119286*alphaR[0]*fR[7]-2.449489742783178*alphaR[1]*fR[6]+4.242640687119286*alphaR[4]*fR[5]-2.449489742783178*(fR[1]*alphaR[4]+alphaR[0]*fR[3])+alphaR[2]*(4.242640687119286*fR[2]-2.449489742783178*fR[0])); 
  incr[8] = -0.125*(2.449489742783178*(alphaR[2]*fR[15]+alphaR[4]*fR[14])-1.414213562373095*alphaR[2]*fR[13]+2.449489742783178*alphaR[0]*fR[12]-1.414213562373095*alphaR[4]*fR[10]+2.449489742783178*alphaR[1]*fR[9]-1.414213562373095*alphaR[0]*fR[8]+2.449489742783178*alphaR[3]*fR[5]+(2.449489742783178*fR[2]-1.414213562373095*fR[0])*alphaR[5]-1.414213562373095*(alphaR[1]*fR[4]+fR[1]*alphaR[3])); 
  incr[9] = 0.125*(4.242640687119286*(alphaR[4]*fR[15]+alphaR[2]*fR[14])-2.449489742783178*alphaR[4]*fR[13]+4.242640687119286*alphaR[1]*fR[12]-2.449489742783178*alphaR[2]*fR[10]+4.242640687119286*alphaR[0]*fR[9]-2.449489742783178*alphaR[1]*fR[8]+4.242640687119286*alphaR[5]*fR[5]-2.449489742783178*(fR[1]*alphaR[5]+alphaR[0]*fR[4])+(4.242640687119286*fR[2]-2.449489742783178*fR[0])*alphaR[3]); 
  incr[10] = -0.125*(2.449489742783178*(alphaR[1]*fR[15]+alphaR[0]*fR[14])-1.414213562373095*alphaR[1]*fR[13]+2.449489742783178*(alphaR[4]*fR[12]+alphaR[5]*fR[11])-1.414213562373095*alphaR[0]*fR[10]+2.449489742783178*alphaR[2]*fR[9]-1.414213562373095*alphaR[4]*fR[8]+2.449489742783178*alphaR[3]*fR[7]-1.414213562373095*(alphaR[5]*fR[6]+alphaR[2]*fR[4]+alphaR[3]*fR[3])); 
  incr[11] = 0.125*(4.242640687119286*(alphaR[3]*fR[15]+alphaR[5]*fR[14])-2.449489742783178*alphaR[3]*fR[13]+4.242640687119286*alphaR[0]*fR[11]-2.449489742783178*alphaR[5]*fR[10]+4.242640687119286*alphaR[1]*fR[7]-2.449489742783178*alphaR[0]*fR[6]+4.242640687119286*alphaR[2]*fR[5]+(4.242640687119286*fR[2]-2.449489742783178*fR[0])*alphaR[4]-2.449489742783178*(alphaR[1]*fR[3]+fR[1]*alphaR[2])); 
  incr[12] = 0.125*(4.242640687119286*(alphaR[2]*fR[15]+alphaR[4]*fR[14])-2.449489742783178*alphaR[2]*fR[13]+4.242640687119286*alphaR[0]*fR[12]-2.449489742783178*alphaR[4]*fR[10]+4.242640687119286*alphaR[1]*fR[9]-2.449489742783178*alphaR[0]*fR[8]+4.242640687119286*alphaR[3]*fR[5]+(4.242640687119286*fR[2]-2.449489742783178*fR[0])*alphaR[5]-2.449489742783178*(alphaR[1]*fR[4]+fR[1]*alphaR[3])); 
  incr[13] = -0.125*(2.449489742783178*(alphaR[0]*fR[15]+alphaR[1]*fR[14])-1.414213562373095*alphaR[0]*fR[13]+2.449489742783178*(alphaR[2]*fR[12]+alphaR[3]*fR[11])-1.414213562373095*alphaR[1]*fR[10]+2.449489742783178*alphaR[4]*fR[9]-1.414213562373095*alphaR[2]*fR[8]+2.449489742783178*alphaR[5]*fR[7]-1.414213562373095*(alphaR[3]*fR[6]+fR[3]*alphaR[5]+alphaR[4]*fR[4])); 
  incr[14] = 0.125*(4.242640687119286*(alphaR[1]*fR[15]+alphaR[0]*fR[14])-2.449489742783178*alphaR[1]*fR[13]+4.242640687119286*(alphaR[4]*fR[12]+alphaR[5]*fR[11])-2.449489742783178*alphaR[0]*fR[10]+4.242640687119286*alphaR[2]*fR[9]-2.449489742783178*alphaR[4]*fR[8]+4.242640687119286*alphaR[3]*fR[7]-2.449489742783178*(alphaR[5]*fR[6]+alphaR[2]*fR[4]+alphaR[3]*fR[3])); 
  incr[15] = 0.125*(4.242640687119286*(alphaR[0]*fR[15]+alphaR[1]*fR[14])-2.449489742783178*alphaR[0]*fR[13]+4.242640687119286*(alphaR[2]*fR[12]+alphaR[3]*fR[11])-2.449489742783178*alphaR[1]*fR[10]+4.242640687119286*alphaR[4]*fR[9]-2.449489742783178*alphaR[2]*fR[8]+4.242640687119286*alphaR[5]*fR[7]-2.449489742783178*(alphaR[3]*fR[6]+fR[3]*alphaR[5]+alphaR[4]*fR[4])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*(alphaR[5]+alphaR[4])-0.3535533905932737*(alphaR[3]+alphaR[2]+alphaR[1])+0.3535533905932737*alphaR[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])+0.25*fR[13]-0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*fL[14]-0.25*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fR[11])+0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[1]+alphaR[0])-0.3535533905932737*(alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14])+0.25*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[5]-0.3535533905932737*(alphaR[4]+alphaR[3])+0.3535533905932737*alphaR[2]-0.3535533905932737*alphaR[1]+0.3535533905932737*alphaR[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*fL[14]+0.25*(fR[13]+fL[13])-0.4330127018922193*fR[12]+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaR[5])+0.3535533905932737*alphaR[4]-0.3535533905932737*alphaR[3]+0.3535533905932737*(alphaR[2]+alphaR[1]+alphaR[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.25*fR[13]-0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14])-0.25*(fR[13]+fL[13])+0.4330127018922193*fR[12]-0.4330127018922193*(fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaR[5])+0.3535533905932737*(alphaR[4]+alphaR[3])-0.3535533905932737*(alphaR[2]+alphaR[1])+0.3535533905932737*alphaR[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*fL[14]+0.25*(fR[13]+fL[13])+0.4330127018922193*fR[12]-0.4330127018922193*(fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[5]-0.3535533905932737*alphaR[4]+0.3535533905932737*alphaR[3]-0.3535533905932737*alphaR[2]+0.3535533905932737*(alphaR[1]+alphaR[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.25*fR[13]-0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14])-0.25*(fR[13]+fL[13])-0.4330127018922193*fR[12]+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*(alphaR[5]+alphaR[4]))+0.3535533905932737*(alphaR[3]+alphaR[2])-0.3535533905932737*alphaR[1]+0.3535533905932737*alphaR[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])+0.25*fR[13]-0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*fL[14]-0.25*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14])+0.25*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fR[11])+0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alphaR[5]*fUp[5]+alphaR[4]*fUp[4]+alphaR[3]*fUp[3]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.25*(alphaR[3]*fUp[5]+fUp[3]*alphaR[5]+alphaR[2]*fUp[4]+fUp[2]*alphaR[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = -0.4330127018922193*(alphaR[5]*fUp[5]+alphaR[4]*fUp[4]+alphaR[3]*fUp[3]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[3] = 0.25*(alphaR[5]*fUp[7]+alphaR[3]*fUp[6]+alphaR[1]*fUp[4]+fUp[1]*alphaR[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[4] = 0.25*(alphaR[4]*fUp[7]+alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+fUp[1]*alphaR[5]+alphaR[0]*fUp[3]+fUp[0]*alphaR[3]); 
  incr[5] = -0.4330127018922193*(alphaR[3]*fUp[5]+fUp[3]*alphaR[5]+alphaR[2]*fUp[4]+fUp[2]*alphaR[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[6] = 0.25*(alphaR[3]*fUp[7]+alphaR[5]*fUp[6]+alphaR[0]*fUp[4]+fUp[0]*alphaR[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[7] = -0.4330127018922193*(alphaR[5]*fUp[7]+alphaR[3]*fUp[6]+alphaR[1]*fUp[4]+fUp[1]*alphaR[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[8] = 0.25*(alphaR[2]*fUp[7]+alphaR[4]*fUp[6]+alphaR[0]*fUp[5]+fUp[0]*alphaR[5]+alphaR[1]*fUp[3]+fUp[1]*alphaR[3]); 
  incr[9] = -0.4330127018922193*(alphaR[4]*fUp[7]+alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+fUp[1]*alphaR[5]+alphaR[0]*fUp[3]+fUp[0]*alphaR[3]); 
  incr[10] = 0.25*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[4]*fUp[5]+fUp[4]*alphaR[5]+alphaR[2]*fUp[3]+fUp[2]*alphaR[3]); 
  incr[11] = -0.4330127018922193*(alphaR[3]*fUp[7]+alphaR[5]*fUp[6]+alphaR[0]*fUp[4]+fUp[0]*alphaR[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[12] = -0.4330127018922193*(alphaR[2]*fUp[7]+alphaR[4]*fUp[6]+alphaR[0]*fUp[5]+fUp[0]*alphaR[5]+alphaR[1]*fUp[3]+fUp[1]*alphaR[3]); 
  incr[13] = 0.25*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]+fUp[2]*alphaR[5]+alphaR[3]*fUp[4]+fUp[3]*alphaR[4]); 
  incr[14] = -0.4330127018922193*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[4]*fUp[5]+fUp[4]*alphaR[5]+alphaR[2]*fUp[3]+fUp[2]*alphaR[3]); 
  incr[15] = -0.4330127018922193*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]+fUp[2]*alphaR[5]+alphaR[3]*fUp[4]+fUp[3]*alphaR[4]); 

#endif 
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
double GyrokineticSimpleHelicalSurf2x2vSer_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *bmagInv, const double *cmag, const double *BdriftX, const double *BdriftY, const double *phi, const double *fL, const double *fR, double *outL, double *outR) 
{ 
  // bmagInv: 1/bmag.
  // BdriftX,BdriftY: x,y components of gradB/curvature drift.
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

  double hamilR[16]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 

  double BstarXdBmagR[16]; 
  BstarXdBmagR[0] = (2.0*BdriftX[0]*m_*wvparR)/q_; 
  BstarXdBmagR[1] = (2.0*BdriftX[1]*m_*wvparR)/q_; 
  BstarXdBmagR[3] = (1.154700538379252*BdriftX[0]*m_)/(q_*rdvpar2R); 
  BstarXdBmagR[6] = (1.154700538379252*BdriftX[1]*m_)/(q_*rdvpar2R); 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = (2.0*BdriftY[0]*m_*wvparR)/q_; 
  BstarYdBmagR[1] = (2.0*BdriftY[1]*m_*wvparR)/q_; 
  BstarYdBmagR[3] = (1.154700538379252*BdriftY[0]*m_)/(q_*rdvpar2R); 
  BstarYdBmagR[6] = (1.154700538379252*BdriftY[1]*m_)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = (0.1767766952966368*((hamilR[5]*(3.0*BstarYdBmagR[6]-1.732050807568877*BstarYdBmagR[1])+hamilR[2]*(3.0*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0]))*rdy2R+hamilR[1]*(3.0*BstarXdBmagR[3]-1.732050807568877*BstarXdBmagR[0])*rdx2R))/m_; 
  alphaR[1] = (0.1767766952966368*((3.0*hamilR[2]*BstarYdBmagR[6]+(3.0*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0])*hamilR[5]-1.732050807568877*BstarYdBmagR[1]*hamilR[2])*rdy2R+hamilR[1]*(3.0*BstarXdBmagR[6]-1.732050807568877*BstarXdBmagR[1])*rdx2R))/m_; 
  alphaR[2] = (0.1767766952966368*(3.0*BstarXdBmagR[3]-1.732050807568877*BstarXdBmagR[0])*hamilR[5]*rdx2R)/m_; 
  alphaR[3] = (0.1767766952966368*(3.0*BstarXdBmagR[3]-1.732050807568877*BstarXdBmagR[0])*hamilR[8]*rdx2R)/m_; 
  alphaR[4] = (0.1767766952966368*hamilR[5]*(3.0*BstarXdBmagR[6]-1.732050807568877*BstarXdBmagR[1])*rdx2R)/m_; 
  alphaR[5] = (0.1767766952966368*(3.0*BstarXdBmagR[6]-1.732050807568877*BstarXdBmagR[1])*hamilR[8]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*((3.0*hamilR[5]*BstarYdBmagR[6]-1.732050807568877*BstarYdBmagR[1]*hamilR[5]+3.0*hamilR[2]*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0]*hamilR[2])*rdy2R+(3.0*hamilR[1]*BstarXdBmagR[3]-1.732050807568877*BstarXdBmagR[0]*hamilR[1])*rdx2R))/m_; 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.125*(2.449489742783178*(alphaR[5]*fL[13]+alphaR[4]*fL[11]+alphaR[3]*fL[10])+1.414213562373095*alphaR[5]*fL[8]+2.449489742783178*(alphaR[2]*fL[7]+alphaR[1]*fL[6])+1.414213562373095*(alphaR[4]*fL[5]+alphaR[3]*fL[4])+2.449489742783178*alphaR[0]*fL[3]+1.414213562373095*(alphaR[2]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.125*(2.449489742783178*(alphaR[3]*fL[13]+alphaR[2]*fL[11]+alphaR[5]*fL[10])+1.414213562373095*alphaR[3]*fL[8]+2.449489742783178*(alphaR[4]*fL[7]+alphaR[0]*fL[6])+1.414213562373095*(alphaR[2]*fL[5]+fL[4]*alphaR[5]+fL[2]*alphaR[4])+2.449489742783178*alphaR[1]*fL[3]+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = 0.125*(2.449489742783178*(alphaR[5]*fL[15]+alphaR[3]*fL[14])+1.414213562373095*alphaR[5]*fL[12]+2.449489742783178*alphaR[1]*fL[11]+1.414213562373095*alphaR[3]*fL[9]+2.449489742783178*(alphaR[0]*fL[7]+alphaR[4]*fL[6])+1.414213562373095*(alphaR[1]*fL[5]+fL[1]*alphaR[4])+2.449489742783178*alphaR[2]*fL[3]+1.414213562373095*(alphaR[0]*fL[2]+fL[0]*alphaR[2])); 
  incr[3] = -0.125*(4.242640687119286*(alphaR[5]*fL[13]+alphaR[4]*fL[11]+alphaR[3]*fL[10])+2.449489742783178*alphaR[5]*fL[8]+4.242640687119286*(alphaR[2]*fL[7]+alphaR[1]*fL[6])+2.449489742783178*(alphaR[4]*fL[5]+alphaR[3]*fL[4])+4.242640687119286*alphaR[0]*fL[3]+2.449489742783178*(alphaR[2]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[4] = 0.125*(2.449489742783178*(alphaR[4]*fL[15]+alphaR[2]*fL[14]+alphaR[1]*fL[13])+1.414213562373095*alphaR[4]*fL[12]+2.449489742783178*alphaR[0]*fL[10]+1.414213562373095*(alphaR[2]*fL[9]+alphaR[1]*fL[8])+2.449489742783178*alphaR[5]*fL[6]+1.414213562373095*(fL[1]*alphaR[5]+alphaR[0]*fL[4])+alphaR[3]*(2.449489742783178*fL[3]+1.414213562373095*fL[0])); 
  incr[5] = 0.125*(2.449489742783178*(alphaR[3]*fL[15]+alphaR[5]*fL[14])+1.414213562373095*alphaR[3]*fL[12]+2.449489742783178*alphaR[0]*fL[11]+1.414213562373095*alphaR[5]*fL[9]+2.449489742783178*(alphaR[1]*fL[7]+alphaR[2]*fL[6])+1.414213562373095*alphaR[0]*fL[5]+(2.449489742783178*fL[3]+1.414213562373095*fL[0])*alphaR[4]+1.414213562373095*(alphaR[1]*fL[2]+fL[1]*alphaR[2])); 
  incr[6] = -0.125*(4.242640687119286*(alphaR[3]*fL[13]+alphaR[2]*fL[11]+alphaR[5]*fL[10])+2.449489742783178*alphaR[3]*fL[8]+4.242640687119286*(alphaR[4]*fL[7]+alphaR[0]*fL[6])+2.449489742783178*(alphaR[2]*fL[5]+fL[4]*alphaR[5]+fL[2]*alphaR[4])+4.242640687119286*alphaR[1]*fL[3]+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[7] = -0.125*(4.242640687119286*(alphaR[5]*fL[15]+alphaR[3]*fL[14])+2.449489742783178*alphaR[5]*fL[12]+4.242640687119286*alphaR[1]*fL[11]+2.449489742783178*alphaR[3]*fL[9]+4.242640687119286*(alphaR[0]*fL[7]+alphaR[4]*fL[6])+2.449489742783178*(alphaR[1]*fL[5]+fL[1]*alphaR[4])+4.242640687119286*alphaR[2]*fL[3]+2.449489742783178*(alphaR[0]*fL[2]+fL[0]*alphaR[2])); 
  incr[8] = 0.125*(2.449489742783178*(alphaR[2]*fL[15]+alphaR[4]*fL[14]+alphaR[0]*fL[13])+1.414213562373095*alphaR[2]*fL[12]+2.449489742783178*alphaR[1]*fL[10]+1.414213562373095*(alphaR[4]*fL[9]+alphaR[0]*fL[8])+2.449489742783178*alphaR[3]*fL[6]+(2.449489742783178*fL[3]+1.414213562373095*fL[0])*alphaR[5]+1.414213562373095*(alphaR[1]*fL[4]+fL[1]*alphaR[3])); 
  incr[9] = 0.125*(2.449489742783178*(alphaR[1]*fL[15]+alphaR[0]*fL[14]+alphaR[4]*fL[13])+1.414213562373095*alphaR[1]*fL[12]+2.449489742783178*(alphaR[5]*fL[11]+alphaR[2]*fL[10])+1.414213562373095*(alphaR[0]*fL[9]+alphaR[4]*fL[8])+2.449489742783178*alphaR[3]*fL[7]+1.414213562373095*(alphaR[5]*fL[5]+alphaR[2]*fL[4]+fL[2]*alphaR[3])); 
  incr[10] = -0.125*(4.242640687119286*(alphaR[4]*fL[15]+alphaR[2]*fL[14]+alphaR[1]*fL[13])+2.449489742783178*alphaR[4]*fL[12]+4.242640687119286*alphaR[0]*fL[10]+2.449489742783178*(alphaR[2]*fL[9]+alphaR[1]*fL[8])+4.242640687119286*alphaR[5]*fL[6]+2.449489742783178*(fL[1]*alphaR[5]+alphaR[0]*fL[4])+alphaR[3]*(4.242640687119286*fL[3]+2.449489742783178*fL[0])); 
  incr[11] = -0.125*(4.242640687119286*(alphaR[3]*fL[15]+alphaR[5]*fL[14])+2.449489742783178*alphaR[3]*fL[12]+4.242640687119286*alphaR[0]*fL[11]+2.449489742783178*alphaR[5]*fL[9]+4.242640687119286*(alphaR[1]*fL[7]+alphaR[2]*fL[6])+2.449489742783178*alphaR[0]*fL[5]+(4.242640687119286*fL[3]+2.449489742783178*fL[0])*alphaR[4]+2.449489742783178*(alphaR[1]*fL[2]+fL[1]*alphaR[2])); 
  incr[12] = 0.125*(2.449489742783178*(alphaR[0]*fL[15]+alphaR[1]*fL[14]+alphaR[2]*fL[13])+1.414213562373095*alphaR[0]*fL[12]+2.449489742783178*(alphaR[3]*fL[11]+alphaR[4]*fL[10])+1.414213562373095*(alphaR[1]*fL[9]+alphaR[2]*fL[8])+2.449489742783178*alphaR[5]*fL[7]+1.414213562373095*(alphaR[3]*fL[5]+fL[2]*alphaR[5]+alphaR[4]*fL[4])); 
  incr[13] = -0.125*(4.242640687119286*(alphaR[2]*fL[15]+alphaR[4]*fL[14]+alphaR[0]*fL[13])+2.449489742783178*alphaR[2]*fL[12]+4.242640687119286*alphaR[1]*fL[10]+2.449489742783178*(alphaR[4]*fL[9]+alphaR[0]*fL[8])+4.242640687119286*alphaR[3]*fL[6]+(4.242640687119286*fL[3]+2.449489742783178*fL[0])*alphaR[5]+2.449489742783178*(alphaR[1]*fL[4]+fL[1]*alphaR[3])); 
  incr[14] = -0.125*(4.242640687119286*(alphaR[1]*fL[15]+alphaR[0]*fL[14]+alphaR[4]*fL[13])+2.449489742783178*alphaR[1]*fL[12]+4.242640687119286*(alphaR[5]*fL[11]+alphaR[2]*fL[10])+2.449489742783178*(alphaR[0]*fL[9]+alphaR[4]*fL[8])+4.242640687119286*alphaR[3]*fL[7]+2.449489742783178*(alphaR[5]*fL[5]+alphaR[2]*fL[4]+fL[2]*alphaR[3])); 
  incr[15] = -0.125*(4.242640687119286*(alphaR[0]*fL[15]+alphaR[1]*fL[14]+alphaR[2]*fL[13])+2.449489742783178*alphaR[0]*fL[12]+4.242640687119286*(alphaR[3]*fL[11]+alphaR[4]*fL[10])+2.449489742783178*(alphaR[1]*fL[9]+alphaR[2]*fL[8])+4.242640687119286*alphaR[5]*fL[7]+2.449489742783178*(alphaR[3]*fL[5]+fL[2]*alphaR[5]+alphaR[4]*fL[4])); 
  } else { 
  incr[0] = -0.125*(2.449489742783178*(alphaR[5]*fR[13]+alphaR[4]*fR[11]+alphaR[3]*fR[10])-1.414213562373095*alphaR[5]*fR[8]+2.449489742783178*(alphaR[2]*fR[7]+alphaR[1]*fR[6])-1.414213562373095*(alphaR[4]*fR[5]+alphaR[3]*fR[4])+2.449489742783178*alphaR[0]*fR[3]-1.414213562373095*(alphaR[2]*fR[2]+alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.125*(2.449489742783178*(alphaR[3]*fR[13]+alphaR[2]*fR[11]+alphaR[5]*fR[10])-1.414213562373095*alphaR[3]*fR[8]+2.449489742783178*(alphaR[4]*fR[7]+alphaR[0]*fR[6])-1.414213562373095*(alphaR[2]*fR[5]+fR[4]*alphaR[5]+fR[2]*alphaR[4])+2.449489742783178*alphaR[1]*fR[3]-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = -0.125*(2.449489742783178*(alphaR[5]*fR[15]+alphaR[3]*fR[14])-1.414213562373095*alphaR[5]*fR[12]+2.449489742783178*alphaR[1]*fR[11]-1.414213562373095*alphaR[3]*fR[9]+2.449489742783178*(alphaR[0]*fR[7]+alphaR[4]*fR[6])-1.414213562373095*(alphaR[1]*fR[5]+fR[1]*alphaR[4])+2.449489742783178*alphaR[2]*fR[3]-1.414213562373095*(alphaR[0]*fR[2]+fR[0]*alphaR[2])); 
  incr[3] = 0.125*(4.242640687119286*(alphaR[5]*fR[13]+alphaR[4]*fR[11]+alphaR[3]*fR[10])-2.449489742783178*alphaR[5]*fR[8]+4.242640687119286*(alphaR[2]*fR[7]+alphaR[1]*fR[6])-2.449489742783178*(alphaR[4]*fR[5]+alphaR[3]*fR[4])+4.242640687119286*alphaR[0]*fR[3]-2.449489742783178*(alphaR[2]*fR[2]+alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[4] = -0.125*(2.449489742783178*(alphaR[4]*fR[15]+alphaR[2]*fR[14]+alphaR[1]*fR[13])-1.414213562373095*alphaR[4]*fR[12]+2.449489742783178*alphaR[0]*fR[10]-1.414213562373095*(alphaR[2]*fR[9]+alphaR[1]*fR[8])+2.449489742783178*alphaR[5]*fR[6]-1.414213562373095*(fR[1]*alphaR[5]+alphaR[0]*fR[4])+alphaR[3]*(2.449489742783178*fR[3]-1.414213562373095*fR[0])); 
  incr[5] = -0.125*(2.449489742783178*(alphaR[3]*fR[15]+alphaR[5]*fR[14])-1.414213562373095*alphaR[3]*fR[12]+2.449489742783178*alphaR[0]*fR[11]-1.414213562373095*alphaR[5]*fR[9]+2.449489742783178*(alphaR[1]*fR[7]+alphaR[2]*fR[6])-1.414213562373095*alphaR[0]*fR[5]+(2.449489742783178*fR[3]-1.414213562373095*fR[0])*alphaR[4]-1.414213562373095*(alphaR[1]*fR[2]+fR[1]*alphaR[2])); 
  incr[6] = 0.125*(4.242640687119286*(alphaR[3]*fR[13]+alphaR[2]*fR[11]+alphaR[5]*fR[10])-2.449489742783178*alphaR[3]*fR[8]+4.242640687119286*(alphaR[4]*fR[7]+alphaR[0]*fR[6])-2.449489742783178*(alphaR[2]*fR[5]+fR[4]*alphaR[5]+fR[2]*alphaR[4])+4.242640687119286*alphaR[1]*fR[3]-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[7] = 0.125*(4.242640687119286*(alphaR[5]*fR[15]+alphaR[3]*fR[14])-2.449489742783178*alphaR[5]*fR[12]+4.242640687119286*alphaR[1]*fR[11]-2.449489742783178*alphaR[3]*fR[9]+4.242640687119286*(alphaR[0]*fR[7]+alphaR[4]*fR[6])-2.449489742783178*(alphaR[1]*fR[5]+fR[1]*alphaR[4])+4.242640687119286*alphaR[2]*fR[3]-2.449489742783178*(alphaR[0]*fR[2]+fR[0]*alphaR[2])); 
  incr[8] = -0.125*(2.449489742783178*(alphaR[2]*fR[15]+alphaR[4]*fR[14]+alphaR[0]*fR[13])-1.414213562373095*alphaR[2]*fR[12]+2.449489742783178*alphaR[1]*fR[10]-1.414213562373095*(alphaR[4]*fR[9]+alphaR[0]*fR[8])+2.449489742783178*alphaR[3]*fR[6]+(2.449489742783178*fR[3]-1.414213562373095*fR[0])*alphaR[5]-1.414213562373095*(alphaR[1]*fR[4]+fR[1]*alphaR[3])); 
  incr[9] = -0.125*(2.449489742783178*(alphaR[1]*fR[15]+alphaR[0]*fR[14]+alphaR[4]*fR[13])-1.414213562373095*alphaR[1]*fR[12]+2.449489742783178*(alphaR[5]*fR[11]+alphaR[2]*fR[10])-1.414213562373095*(alphaR[0]*fR[9]+alphaR[4]*fR[8])+2.449489742783178*alphaR[3]*fR[7]-1.414213562373095*(alphaR[5]*fR[5]+alphaR[2]*fR[4]+fR[2]*alphaR[3])); 
  incr[10] = 0.125*(4.242640687119286*(alphaR[4]*fR[15]+alphaR[2]*fR[14]+alphaR[1]*fR[13])-2.449489742783178*alphaR[4]*fR[12]+4.242640687119286*alphaR[0]*fR[10]-2.449489742783178*(alphaR[2]*fR[9]+alphaR[1]*fR[8])+4.242640687119286*alphaR[5]*fR[6]-2.449489742783178*(fR[1]*alphaR[5]+alphaR[0]*fR[4])+alphaR[3]*(4.242640687119286*fR[3]-2.449489742783178*fR[0])); 
  incr[11] = 0.125*(4.242640687119286*(alphaR[3]*fR[15]+alphaR[5]*fR[14])-2.449489742783178*alphaR[3]*fR[12]+4.242640687119286*alphaR[0]*fR[11]-2.449489742783178*alphaR[5]*fR[9]+4.242640687119286*(alphaR[1]*fR[7]+alphaR[2]*fR[6])-2.449489742783178*alphaR[0]*fR[5]+(4.242640687119286*fR[3]-2.449489742783178*fR[0])*alphaR[4]-2.449489742783178*(alphaR[1]*fR[2]+fR[1]*alphaR[2])); 
  incr[12] = -0.125*(2.449489742783178*(alphaR[0]*fR[15]+alphaR[1]*fR[14]+alphaR[2]*fR[13])-1.414213562373095*alphaR[0]*fR[12]+2.449489742783178*(alphaR[3]*fR[11]+alphaR[4]*fR[10])-1.414213562373095*(alphaR[1]*fR[9]+alphaR[2]*fR[8])+2.449489742783178*alphaR[5]*fR[7]-1.414213562373095*(alphaR[3]*fR[5]+fR[2]*alphaR[5]+alphaR[4]*fR[4])); 
  incr[13] = 0.125*(4.242640687119286*(alphaR[2]*fR[15]+alphaR[4]*fR[14]+alphaR[0]*fR[13])-2.449489742783178*alphaR[2]*fR[12]+4.242640687119286*alphaR[1]*fR[10]-2.449489742783178*(alphaR[4]*fR[9]+alphaR[0]*fR[8])+4.242640687119286*alphaR[3]*fR[6]+(4.242640687119286*fR[3]-2.449489742783178*fR[0])*alphaR[5]-2.449489742783178*(alphaR[1]*fR[4]+fR[1]*alphaR[3])); 
  incr[14] = 0.125*(4.242640687119286*(alphaR[1]*fR[15]+alphaR[0]*fR[14]+alphaR[4]*fR[13])-2.449489742783178*alphaR[1]*fR[12]+4.242640687119286*(alphaR[5]*fR[11]+alphaR[2]*fR[10])-2.449489742783178*(alphaR[0]*fR[9]+alphaR[4]*fR[8])+4.242640687119286*alphaR[3]*fR[7]-2.449489742783178*(alphaR[5]*fR[5]+alphaR[2]*fR[4]+fR[2]*alphaR[3])); 
  incr[15] = 0.125*(4.242640687119286*(alphaR[0]*fR[15]+alphaR[1]*fR[14]+alphaR[2]*fR[13])-2.449489742783178*alphaR[0]*fR[12]+4.242640687119286*(alphaR[3]*fR[11]+alphaR[4]*fR[10])-2.449489742783178*(alphaR[1]*fR[9]+alphaR[2]*fR[8])+4.242640687119286*alphaR[5]*fR[7]-2.449489742783178*(alphaR[3]*fR[5]+fR[2]*alphaR[5]+alphaR[4]*fR[4])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*(alphaR[5]+alphaR[4])-0.3535533905932737*(alphaR[3]+alphaR[2]+alphaR[1])+0.3535533905932737*alphaR[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])-0.25*(fR[9]+fR[8])+0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14]+fR[13])+0.4330127018922193*(fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[1]+alphaR[0])-0.3535533905932737*(alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14]+fR[13])-0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[5]-0.3535533905932737*(alphaR[4]+alphaR[3])+0.3535533905932737*alphaR[2]-0.3535533905932737*alphaR[1]+0.3535533905932737*alphaR[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*(fL[14]+fR[13])+0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaR[5])+0.3535533905932737*alphaR[4]-0.3535533905932737*alphaR[3]+0.3535533905932737*(alphaR[2]+alphaR[1]+alphaR[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]+fR[13]+fL[13]))+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])+0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*(fR[2]+fR[1]+fR[0])+0.25*(fL[2]+fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14]+fR[13])-0.4330127018922193*(fL[15]+fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9]+fR[8]+fL[8])-0.4330127018922193*(fR[7]+fR[6])+0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2]+fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaR[5])+0.3535533905932737*(alphaR[4]+alphaR[3])-0.3535533905932737*(alphaR[2]+alphaR[1])+0.3535533905932737*alphaR[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*(fR[5]+fR[4])+0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14]+fR[13])-0.4330127018922193*(fL[14]+fL[13])+0.25*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fR[10])+0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[5]-0.3535533905932737*alphaR[4]+0.3535533905932737*alphaR[3]-0.3535533905932737*alphaR[2]+0.3535533905932737*(alphaR[1]+alphaR[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14]+fR[13])+0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*(alphaR[5]+alphaR[4]))+0.3535533905932737*(alphaR[3]+alphaR[2])-0.3535533905932737*alphaR[1]+0.3535533905932737*alphaR[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*(fL[14]+fR[13])-0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaR[5]+alphaR[4]+alphaR[3]+alphaR[2]+alphaR[1]+alphaR[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]+fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])-0.25*(fR[9]+fR[8])+0.25*(fL[9]+fL[8])+0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*(fR[5]+fR[4])+0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])-0.25*(fR[2]+fR[1]+fR[0])+0.25*(fL[2]+fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14]+fR[13])+0.4330127018922193*(fL[15]+fL[14]+fL[13])+0.25*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fR[10])+0.4330127018922193*(fL[11]+fL[10])+0.25*(fR[9]+fL[9]+fR[8]+fL[8])-0.4330127018922193*(fR[7]+fR[6])+0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2]+fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alphaR[5]*fUp[5]+alphaR[4]*fUp[4]+alphaR[3]*fUp[3]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.25*(alphaR[3]*fUp[5]+fUp[3]*alphaR[5]+alphaR[2]*fUp[4]+fUp[2]*alphaR[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = 0.25*(alphaR[5]*fUp[7]+alphaR[3]*fUp[6]+alphaR[1]*fUp[4]+fUp[1]*alphaR[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[3] = -0.4330127018922193*(alphaR[5]*fUp[5]+alphaR[4]*fUp[4]+alphaR[3]*fUp[3]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[4] = 0.25*(alphaR[4]*fUp[7]+alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+fUp[1]*alphaR[5]+alphaR[0]*fUp[3]+fUp[0]*alphaR[3]); 
  incr[5] = 0.25*(alphaR[3]*fUp[7]+alphaR[5]*fUp[6]+alphaR[0]*fUp[4]+fUp[0]*alphaR[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[6] = -0.4330127018922193*(alphaR[3]*fUp[5]+fUp[3]*alphaR[5]+alphaR[2]*fUp[4]+fUp[2]*alphaR[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[7] = -0.4330127018922193*(alphaR[5]*fUp[7]+alphaR[3]*fUp[6]+alphaR[1]*fUp[4]+fUp[1]*alphaR[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[8] = 0.25*(alphaR[2]*fUp[7]+alphaR[4]*fUp[6]+alphaR[0]*fUp[5]+fUp[0]*alphaR[5]+alphaR[1]*fUp[3]+fUp[1]*alphaR[3]); 
  incr[9] = 0.25*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[4]*fUp[5]+fUp[4]*alphaR[5]+alphaR[2]*fUp[3]+fUp[2]*alphaR[3]); 
  incr[10] = -0.4330127018922193*(alphaR[4]*fUp[7]+alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+fUp[1]*alphaR[5]+alphaR[0]*fUp[3]+fUp[0]*alphaR[3]); 
  incr[11] = -0.4330127018922193*(alphaR[3]*fUp[7]+alphaR[5]*fUp[6]+alphaR[0]*fUp[4]+fUp[0]*alphaR[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[12] = 0.25*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]+fUp[2]*alphaR[5]+alphaR[3]*fUp[4]+fUp[3]*alphaR[4]); 
  incr[13] = -0.4330127018922193*(alphaR[2]*fUp[7]+alphaR[4]*fUp[6]+alphaR[0]*fUp[5]+fUp[0]*alphaR[5]+alphaR[1]*fUp[3]+fUp[1]*alphaR[3]); 
  incr[14] = -0.4330127018922193*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[4]*fUp[5]+fUp[4]*alphaR[5]+alphaR[2]*fUp[3]+fUp[2]*alphaR[3]); 
  incr[15] = -0.4330127018922193*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]+fUp[2]*alphaR[5]+alphaR[3]*fUp[4]+fUp[3]*alphaR[4]); 

#endif 
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
