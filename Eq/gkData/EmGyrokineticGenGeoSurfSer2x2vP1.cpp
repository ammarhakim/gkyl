#include <GyrokineticModDecl.h>
double EmGyrokineticGenGeoSurf2x2vSer_xL_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
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

  double hamilL[16]; 
  hamilL[0] = (0.6666666666666666*(3.0*rdvpar2SqL*(m_*wvparSqL+bmag[0]*wmuL+phi[0]*q_)+m_))/rdvpar2SqL; 
  hamilL[1] = 2.0*phi[1]*q_; 
  hamilL[2] = 2.0*phi[2]*q_; 
  hamilL[3] = (2.309401076758503*m_*wvparL)/rdvpar2L; 
  hamilL[4] = (1.154700538379252*bmag[0])/rdmu2L; 
  hamilL[5] = 2.0*phi[3]*q_; 

  double BstarXdBmagL[16]; 
  BstarXdBmagL[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2L; 
  BstarXdBmagL[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2L; 

  double alphaL[8]; 
  alphaL[0] = -(0.125*(b_z[0]*jacobTotInv[0]*(4.242640687119286*hamilL[5]+2.449489742783178*hamilL[2])*m_*rdy2L+((-4.242640687119286*BstarXdBmagL[1])-2.449489742783178*BstarXdBmagL[0])*hamilL[3]*q_*rdvpar2L))/(m_*q_); 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgL = -(0.03125*((3.0*b_z[0]*jacobTotInv[0]*hamilL[5]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilL[2])*m_*rdy2L+((-3.0*BstarXdBmagL[1])-1.732050807568877*BstarXdBmagL[0])*hamilL[3]*q_*rdvpar2L))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgL>0) { 
  incr[0] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[1]+fL[0]); 
  incr[1] = -0.1767766952966368*alphaL[0]*(3.0*fL[1]+1.732050807568877*fL[0]); 
  incr[2] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[5]+fL[2]); 
  incr[3] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[6]+fL[3]); 
  incr[4] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[8]+fL[4]); 
  incr[5] = -0.1767766952966368*alphaL[0]*(3.0*fL[5]+1.732050807568877*fL[2]); 
  incr[6] = -0.1767766952966368*alphaL[0]*(3.0*fL[6]+1.732050807568877*fL[3]); 
  incr[7] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[11]+fL[7]); 
  incr[8] = -0.1767766952966368*alphaL[0]*(3.0*fL[8]+1.732050807568877*fL[4]); 
  incr[9] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[12]+fL[9]); 
  incr[10] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[13]+fL[10]); 
  incr[11] = -0.1767766952966368*alphaL[0]*(3.0*fL[11]+1.732050807568877*fL[7]); 
  incr[12] = -0.1767766952966368*alphaL[0]*(3.0*fL[12]+1.732050807568877*fL[9]); 
  incr[13] = -0.1767766952966368*alphaL[0]*(3.0*fL[13]+1.732050807568877*fL[10]); 
  incr[14] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[15]+fL[14]); 
  incr[15] = -0.1767766952966368*alphaL[0]*(3.0*fL[15]+1.732050807568877*fL[14]); 
  } else { 
  incr[0] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[1]-1.0*fR[0]); 
  incr[1] = 0.1767766952966368*alphaL[0]*(3.0*fR[1]-1.732050807568877*fR[0]); 
  incr[2] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[5]-1.0*fR[2]); 
  incr[3] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[6]-1.0*fR[3]); 
  incr[4] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[8]-1.0*fR[4]); 
  incr[5] = 0.1767766952966368*alphaL[0]*(3.0*fR[5]-1.732050807568877*fR[2]); 
  incr[6] = 0.1767766952966368*alphaL[0]*(3.0*fR[6]-1.732050807568877*fR[3]); 
  incr[7] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[11]-1.0*fR[7]); 
  incr[8] = 0.1767766952966368*alphaL[0]*(3.0*fR[8]-1.732050807568877*fR[4]); 
  incr[9] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[12]-1.0*fR[9]); 
  incr[10] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[13]-1.0*fR[10]); 
  incr[11] = 0.1767766952966368*alphaL[0]*(3.0*fR[11]-1.732050807568877*fR[7]); 
  incr[12] = 0.1767766952966368*alphaL[0]*(3.0*fR[12]-1.732050807568877*fR[9]); 
  incr[13] = 0.1767766952966368*alphaL[0]*(3.0*fR[13]-1.732050807568877*fR[10]); 
  incr[14] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[15]-1.0*fR[14]); 
  incr[15] = 0.1767766952966368*alphaL[0]*(3.0*fR[15]-1.732050807568877*fR[14]); 
  }
#elif upwindType == QUAD 
  double alphaOrdL;
  double fUpOrd[8];
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12]+fR[11]+fL[11])-0.25*(fR[10]+fR[9])+0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*(fR[4]+fR[3]+fR[2])-0.25*(fL[4]+fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fR[12]+fR[11])+0.4330127018922193*(fL[13]+fL[12]+fL[11])+0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4]+fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*(fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3]+fR[2])+0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fR[5])+0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])+0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3]+fR[2])-0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*(fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12]+fR[11]+fL[11])-0.25*(fR[10]+fR[9])+0.25*(fL[10]+fL[9])+0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])-0.25*(fR[4]+fR[3]+fR[2])+0.25*(fL[4]+fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fR[12]+fR[11])+0.4330127018922193*(fL[13]+fL[12]+fL[11])+0.25*(fR[10]+fL[10]+fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fR[5])+0.4330127018922193*(fL[6]+fL[5])+0.25*(fR[4]+fL[4]+fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*alphaL[0]*fUp[0]; 
  incr[1] = -0.4330127018922193*alphaL[0]*fUp[0]; 
  incr[2] = 0.25*alphaL[0]*fUp[1]; 
  incr[3] = 0.25*alphaL[0]*fUp[2]; 
  incr[4] = 0.25*alphaL[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alphaL[0]*fUp[1]; 
  incr[6] = -0.4330127018922193*alphaL[0]*fUp[2]; 
  incr[7] = 0.25*alphaL[0]*fUp[4]; 
  incr[8] = -0.4330127018922193*alphaL[0]*fUp[3]; 
  incr[9] = 0.25*alphaL[0]*fUp[5]; 
  incr[10] = 0.25*alphaL[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alphaL[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alphaL[0]*fUp[5]; 
  incr[13] = -0.4330127018922193*alphaL[0]*fUp[6]; 
  incr[14] = 0.25*alphaL[0]*fUp[7]; 
  incr[15] = -0.4330127018922193*alphaL[0]*fUp[7]; 

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

  return std::abs(alphaSurfAvgL);
} 
double EmGyrokineticGenGeoSurf2x2vSer_xR_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
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
  BstarXdBmagR[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2R; 
  BstarXdBmagR[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2R; 

  double alphaR[8]; 
  alphaR[0] = (0.125*(b_z[0]*jacobTotInv[0]*(4.242640687119286*hamilR[5]-2.449489742783178*hamilR[2])*m_*rdy2R+(2.449489742783178*BstarXdBmagR[0]-4.242640687119286*BstarXdBmagR[1])*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*((3.0*b_z[0]*jacobTotInv[0]*hamilR[5]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamilR[2])*m_*rdy2R+(1.732050807568877*BstarXdBmagR[0]-3.0*BstarXdBmagR[1])*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[1]+fL[0]); 
  incr[1] = -0.1767766952966368*alphaR[0]*(3.0*fL[1]+1.732050807568877*fL[0]); 
  incr[2] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[5]+fL[2]); 
  incr[3] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[6]+fL[3]); 
  incr[4] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[8]+fL[4]); 
  incr[5] = -0.1767766952966368*alphaR[0]*(3.0*fL[5]+1.732050807568877*fL[2]); 
  incr[6] = -0.1767766952966368*alphaR[0]*(3.0*fL[6]+1.732050807568877*fL[3]); 
  incr[7] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[11]+fL[7]); 
  incr[8] = -0.1767766952966368*alphaR[0]*(3.0*fL[8]+1.732050807568877*fL[4]); 
  incr[9] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[12]+fL[9]); 
  incr[10] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[13]+fL[10]); 
  incr[11] = -0.1767766952966368*alphaR[0]*(3.0*fL[11]+1.732050807568877*fL[7]); 
  incr[12] = -0.1767766952966368*alphaR[0]*(3.0*fL[12]+1.732050807568877*fL[9]); 
  incr[13] = -0.1767766952966368*alphaR[0]*(3.0*fL[13]+1.732050807568877*fL[10]); 
  incr[14] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[15]+fL[14]); 
  incr[15] = -0.1767766952966368*alphaR[0]*(3.0*fL[15]+1.732050807568877*fL[14]); 
  } else { 
  incr[0] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[1]-1.0*fR[0]); 
  incr[1] = 0.1767766952966368*alphaR[0]*(3.0*fR[1]-1.732050807568877*fR[0]); 
  incr[2] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[5]-1.0*fR[2]); 
  incr[3] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[6]-1.0*fR[3]); 
  incr[4] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[8]-1.0*fR[4]); 
  incr[5] = 0.1767766952966368*alphaR[0]*(3.0*fR[5]-1.732050807568877*fR[2]); 
  incr[6] = 0.1767766952966368*alphaR[0]*(3.0*fR[6]-1.732050807568877*fR[3]); 
  incr[7] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[11]-1.0*fR[7]); 
  incr[8] = 0.1767766952966368*alphaR[0]*(3.0*fR[8]-1.732050807568877*fR[4]); 
  incr[9] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[12]-1.0*fR[9]); 
  incr[10] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[13]-1.0*fR[10]); 
  incr[11] = 0.1767766952966368*alphaR[0]*(3.0*fR[11]-1.732050807568877*fR[7]); 
  incr[12] = 0.1767766952966368*alphaR[0]*(3.0*fR[12]-1.732050807568877*fR[9]); 
  incr[13] = 0.1767766952966368*alphaR[0]*(3.0*fR[13]-1.732050807568877*fR[10]); 
  incr[14] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[15]-1.0*fR[14]); 
  incr[15] = 0.1767766952966368*alphaR[0]*(3.0*fR[15]-1.732050807568877*fR[14]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12]+fR[11]+fL[11])-0.25*(fR[10]+fR[9])+0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*(fR[4]+fR[3]+fR[2])-0.25*(fL[4]+fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fR[12]+fR[11])+0.4330127018922193*(fL[13]+fL[12]+fL[11])+0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4]+fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*(fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3]+fR[2])+0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fR[5])+0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])+0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3]+fR[2])-0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*(fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
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

  incr[0] = 0.25*alphaR[0]*fUp[0]; 
  incr[1] = -0.4330127018922193*alphaR[0]*fUp[0]; 
  incr[2] = 0.25*alphaR[0]*fUp[1]; 
  incr[3] = 0.25*alphaR[0]*fUp[2]; 
  incr[4] = 0.25*alphaR[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alphaR[0]*fUp[1]; 
  incr[6] = -0.4330127018922193*alphaR[0]*fUp[2]; 
  incr[7] = 0.25*alphaR[0]*fUp[4]; 
  incr[8] = -0.4330127018922193*alphaR[0]*fUp[3]; 
  incr[9] = 0.25*alphaR[0]*fUp[5]; 
  incr[10] = 0.25*alphaR[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alphaR[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alphaR[0]*fUp[5]; 
  incr[13] = -0.4330127018922193*alphaR[0]*fUp[6]; 
  incr[14] = 0.25*alphaR[0]*fUp[7]; 
  incr[15] = -0.4330127018922193*alphaR[0]*fUp[7]; 

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
double EmGyrokineticGenGeoSurf2x2vSer_yL_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
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

  double hamilL[16]; 
  hamilL[0] = (0.6666666666666666*(3.0*rdvpar2SqL*(m_*wvparSqL+bmag[0]*wmuL+phi[0]*q_)+m_))/rdvpar2SqL; 
  hamilL[1] = 2.0*phi[1]*q_; 
  hamilL[2] = 2.0*phi[2]*q_; 
  hamilL[3] = (2.309401076758503*m_*wvparL)/rdvpar2L; 
  hamilL[4] = (1.154700538379252*bmag[0])/rdmu2L; 
  hamilL[5] = 2.0*phi[3]*q_; 

  double BstarYdBmagL[16]; 
  BstarYdBmagL[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2L; 
  BstarYdBmagL[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2L; 

  double alphaL[8]; 
  alphaL[0] = (0.125*(b_z[0]*jacobTotInv[0]*(4.242640687119286*hamilL[5]+2.449489742783178*hamilL[1])*m_*rdx2L+(4.242640687119286*BstarYdBmagL[2]+2.449489742783178*BstarYdBmagL[0])*hamilL[3]*q_*rdvpar2L))/(m_*q_); 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgL = (0.03125*((3.0*b_z[0]*jacobTotInv[0]*hamilL[5]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilL[1])*m_*rdx2L+(3.0*BstarYdBmagL[2]+1.732050807568877*BstarYdBmagL[0])*hamilL[3]*q_*rdvpar2L))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgL>0) { 
  incr[0] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[2]+fL[0]); 
  incr[1] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[5]+fL[1]); 
  incr[2] = -0.1767766952966368*alphaL[0]*(3.0*fL[2]+1.732050807568877*fL[0]); 
  incr[3] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[7]+fL[3]); 
  incr[4] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[9]+fL[4]); 
  incr[5] = -0.1767766952966368*alphaL[0]*(3.0*fL[5]+1.732050807568877*fL[1]); 
  incr[6] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[11]+fL[6]); 
  incr[7] = -0.1767766952966368*alphaL[0]*(3.0*fL[7]+1.732050807568877*fL[3]); 
  incr[8] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[12]+fL[8]); 
  incr[9] = -0.1767766952966368*alphaL[0]*(3.0*fL[9]+1.732050807568877*fL[4]); 
  incr[10] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[14]+fL[10]); 
  incr[11] = -0.1767766952966368*alphaL[0]*(3.0*fL[11]+1.732050807568877*fL[6]); 
  incr[12] = -0.1767766952966368*alphaL[0]*(3.0*fL[12]+1.732050807568877*fL[8]); 
  incr[13] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[15]+fL[13]); 
  incr[14] = -0.1767766952966368*alphaL[0]*(3.0*fL[14]+1.732050807568877*fL[10]); 
  incr[15] = -0.1767766952966368*alphaL[0]*(3.0*fL[15]+1.732050807568877*fL[13]); 
  } else { 
  incr[0] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[2]-1.0*fR[0]); 
  incr[1] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[5]-1.0*fR[1]); 
  incr[2] = 0.1767766952966368*alphaL[0]*(3.0*fR[2]-1.732050807568877*fR[0]); 
  incr[3] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[7]-1.0*fR[3]); 
  incr[4] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[9]-1.0*fR[4]); 
  incr[5] = 0.1767766952966368*alphaL[0]*(3.0*fR[5]-1.732050807568877*fR[1]); 
  incr[6] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[11]-1.0*fR[6]); 
  incr[7] = 0.1767766952966368*alphaL[0]*(3.0*fR[7]-1.732050807568877*fR[3]); 
  incr[8] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[12]-1.0*fR[8]); 
  incr[9] = 0.1767766952966368*alphaL[0]*(3.0*fR[9]-1.732050807568877*fR[4]); 
  incr[10] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[14]-1.0*fR[10]); 
  incr[11] = 0.1767766952966368*alphaL[0]*(3.0*fR[11]-1.732050807568877*fR[6]); 
  incr[12] = 0.1767766952966368*alphaL[0]*(3.0*fR[12]-1.732050807568877*fR[8]); 
  incr[13] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[15]-1.0*fR[13]); 
  incr[14] = 0.1767766952966368*alphaL[0]*(3.0*fR[14]-1.732050807568877*fR[10]); 
  incr[15] = 0.1767766952966368*alphaL[0]*(3.0*fR[15]-1.732050807568877*fR[13]); 
  }
#elif upwindType == QUAD 
  double alphaOrdL;
  double fUpOrd[8];
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])+0.25*fR[13]-0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*fL[14]-0.25*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fR[11])+0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdL)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14])+0.25*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*fL[14]+0.25*(fR[13]+fL[13])-0.4330127018922193*fR[12]+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.25*fR[13]-0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdL)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14])-0.25*(fR[13]+fL[13])+0.4330127018922193*fR[12]-0.4330127018922193*(fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*fL[14]+0.25*(fR[13]+fL[13])+0.4330127018922193*fR[12]-0.4330127018922193*(fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.25*fR[13]-0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdL)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14])-0.25*(fR[13]+fL[13])-0.4330127018922193*fR[12]+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])+0.25*fR[13]-0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*fL[14]-0.25*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdL)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14])+0.25*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fR[11])+0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*alphaL[0]*fUp[0]; 
  incr[1] = 0.25*alphaL[0]*fUp[1]; 
  incr[2] = -0.4330127018922193*alphaL[0]*fUp[0]; 
  incr[3] = 0.25*alphaL[0]*fUp[2]; 
  incr[4] = 0.25*alphaL[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alphaL[0]*fUp[1]; 
  incr[6] = 0.25*alphaL[0]*fUp[4]; 
  incr[7] = -0.4330127018922193*alphaL[0]*fUp[2]; 
  incr[8] = 0.25*alphaL[0]*fUp[5]; 
  incr[9] = -0.4330127018922193*alphaL[0]*fUp[3]; 
  incr[10] = 0.25*alphaL[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alphaL[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alphaL[0]*fUp[5]; 
  incr[13] = 0.25*alphaL[0]*fUp[7]; 
  incr[14] = -0.4330127018922193*alphaL[0]*fUp[6]; 
  incr[15] = -0.4330127018922193*alphaL[0]*fUp[7]; 

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

  return std::abs(alphaSurfAvgL);
} 
double EmGyrokineticGenGeoSurf2x2vSer_yR_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
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
  BstarYdBmagR[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2R; 
  BstarYdBmagR[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2R; 

  double alphaR[8]; 
  alphaR[0] = -(0.125*(b_z[0]*jacobTotInv[0]*(4.242640687119286*hamilR[5]-2.449489742783178*hamilR[1])*m_*rdx2R+(4.242640687119286*BstarYdBmagR[2]-2.449489742783178*BstarYdBmagR[0])*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*((3.0*b_z[0]*jacobTotInv[0]*hamilR[5]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamilR[1])*m_*rdx2R+(3.0*BstarYdBmagR[2]-1.732050807568877*BstarYdBmagR[0])*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[2]+fL[0]); 
  incr[1] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[5]+fL[1]); 
  incr[2] = -0.1767766952966368*alphaR[0]*(3.0*fL[2]+1.732050807568877*fL[0]); 
  incr[3] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[7]+fL[3]); 
  incr[4] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[9]+fL[4]); 
  incr[5] = -0.1767766952966368*alphaR[0]*(3.0*fL[5]+1.732050807568877*fL[1]); 
  incr[6] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[11]+fL[6]); 
  incr[7] = -0.1767766952966368*alphaR[0]*(3.0*fL[7]+1.732050807568877*fL[3]); 
  incr[8] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[12]+fL[8]); 
  incr[9] = -0.1767766952966368*alphaR[0]*(3.0*fL[9]+1.732050807568877*fL[4]); 
  incr[10] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[14]+fL[10]); 
  incr[11] = -0.1767766952966368*alphaR[0]*(3.0*fL[11]+1.732050807568877*fL[6]); 
  incr[12] = -0.1767766952966368*alphaR[0]*(3.0*fL[12]+1.732050807568877*fL[8]); 
  incr[13] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[15]+fL[13]); 
  incr[14] = -0.1767766952966368*alphaR[0]*(3.0*fL[14]+1.732050807568877*fL[10]); 
  incr[15] = -0.1767766952966368*alphaR[0]*(3.0*fL[15]+1.732050807568877*fL[13]); 
  } else { 
  incr[0] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[2]-1.0*fR[0]); 
  incr[1] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[5]-1.0*fR[1]); 
  incr[2] = 0.1767766952966368*alphaR[0]*(3.0*fR[2]-1.732050807568877*fR[0]); 
  incr[3] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[7]-1.0*fR[3]); 
  incr[4] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[9]-1.0*fR[4]); 
  incr[5] = 0.1767766952966368*alphaR[0]*(3.0*fR[5]-1.732050807568877*fR[1]); 
  incr[6] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[11]-1.0*fR[6]); 
  incr[7] = 0.1767766952966368*alphaR[0]*(3.0*fR[7]-1.732050807568877*fR[3]); 
  incr[8] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[12]-1.0*fR[8]); 
  incr[9] = 0.1767766952966368*alphaR[0]*(3.0*fR[9]-1.732050807568877*fR[4]); 
  incr[10] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[14]-1.0*fR[10]); 
  incr[11] = 0.1767766952966368*alphaR[0]*(3.0*fR[11]-1.732050807568877*fR[6]); 
  incr[12] = 0.1767766952966368*alphaR[0]*(3.0*fR[12]-1.732050807568877*fR[8]); 
  incr[13] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[15]-1.0*fR[13]); 
  incr[14] = 0.1767766952966368*alphaR[0]*(3.0*fR[14]-1.732050807568877*fR[10]); 
  incr[15] = 0.1767766952966368*alphaR[0]*(3.0*fR[15]-1.732050807568877*fR[13]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])+0.25*fR[13]-0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*fL[14]-0.25*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fR[11])+0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14])+0.25*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*fL[14]+0.25*(fR[13]+fL[13])-0.4330127018922193*fR[12]+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.25*fR[13]-0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14])-0.25*(fR[13]+fL[13])+0.4330127018922193*fR[12]-0.4330127018922193*(fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*fL[14]+0.25*(fR[13]+fL[13])+0.4330127018922193*fR[12]-0.4330127018922193*(fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.25*fR[13]-0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14])-0.25*(fR[13]+fL[13])-0.4330127018922193*fR[12]+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])+0.25*fR[13]-0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*fL[14]-0.25*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
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

  incr[0] = 0.25*alphaR[0]*fUp[0]; 
  incr[1] = 0.25*alphaR[0]*fUp[1]; 
  incr[2] = -0.4330127018922193*alphaR[0]*fUp[0]; 
  incr[3] = 0.25*alphaR[0]*fUp[2]; 
  incr[4] = 0.25*alphaR[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alphaR[0]*fUp[1]; 
  incr[6] = 0.25*alphaR[0]*fUp[4]; 
  incr[7] = -0.4330127018922193*alphaR[0]*fUp[2]; 
  incr[8] = 0.25*alphaR[0]*fUp[5]; 
  incr[9] = -0.4330127018922193*alphaR[0]*fUp[3]; 
  incr[10] = 0.25*alphaR[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alphaR[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alphaR[0]*fUp[5]; 
  incr[13] = 0.25*alphaR[0]*fUp[7]; 
  incr[14] = -0.4330127018922193*alphaR[0]*fUp[6]; 
  incr[15] = -0.4330127018922193*alphaR[0]*fUp[7]; 

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
double EmGyrokineticGenGeoSurf2x2vSer_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
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
  BstarXdBmagR[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2R; 
  BstarXdBmagR[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2R; 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2R; 
  BstarYdBmagR[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2R; 

  double alphaR[8]; 
  alphaR[0] = -(0.3061862178478971*(BstarYdBmagR[0]*hamilR[2]*rdy2R+BstarXdBmagR[0]*hamilR[1]*rdx2R))/m_; 
  alphaR[1] = -(0.3061862178478971*(BstarYdBmagR[0]*hamilR[5]*rdy2R+BstarXdBmagR[1]*hamilR[1]*rdx2R))/m_; 
  alphaR[2] = -(0.3061862178478971*(BstarYdBmagR[2]*hamilR[2]*rdy2R+BstarXdBmagR[0]*hamilR[5]*rdx2R))/m_; 
  alphaR[4] = -(0.3061862178478971*hamilR[5]*(BstarYdBmagR[2]*rdy2R+BstarXdBmagR[1]*rdx2R))/m_; 

  double alphaUpR[8]; 
  alphaUpR[0] = -(0.1767766952966368*(1.732050807568877*(BstarYdBmagR[0]*hamilR[2]*rdy2R+BstarXdBmagR[0]*hamilR[1]*rdx2R)+8.0*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = -(0.1767766952966368*(1.732050807568877*(BstarYdBmagR[0]*hamilR[5]*rdy2R+BstarXdBmagR[1]*hamilR[1]*rdx2R)+8.0*dApardtPrev[1]*q_))/m_; 
  alphaUpR[2] = -(0.1767766952966368*(1.732050807568877*(BstarYdBmagR[2]*hamilR[2]*rdy2R+BstarXdBmagR[0]*hamilR[5]*rdx2R)+8.0*dApardtPrev[2]*q_))/m_; 
  alphaUpR[4] = -(0.1767766952966368*(1.732050807568877*hamilR[5]*(BstarYdBmagR[2]*rdy2R+BstarXdBmagR[1]*rdx2R)+8.0*dApardtPrev[3]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(1.732050807568877*BstarYdBmagR[0]*hamilR[2]*rdy2R+1.732050807568877*BstarXdBmagR[0]*hamilR[1]*rdx2R+8.0*dApardtPrev[0]*q_))/m_; 

  double incr[16]; 
  double incrEmMod[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.125*(2.449489742783178*(alphaR[4]*fL[11]+alphaR[2]*fL[7]+alphaR[1]*fL[6])+1.414213562373095*alphaR[4]*fL[5]+2.449489742783178*alphaR[0]*fL[3]+1.414213562373095*(alphaR[2]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.125*(2.449489742783178*(alphaR[2]*fL[11]+alphaR[4]*fL[7]+alphaR[0]*fL[6])+1.414213562373095*(alphaR[2]*fL[5]+fL[2]*alphaR[4])+2.449489742783178*alphaR[1]*fL[3]+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = 0.125*(2.449489742783178*(alphaR[1]*fL[11]+alphaR[0]*fL[7]+alphaR[4]*fL[6])+1.414213562373095*(alphaR[1]*fL[5]+fL[1]*alphaR[4])+2.449489742783178*alphaR[2]*fL[3]+1.414213562373095*(alphaR[0]*fL[2]+fL[0]*alphaR[2])); 
  incr[3] = -0.125*(4.242640687119286*(alphaR[4]*fL[11]+alphaR[2]*fL[7]+alphaR[1]*fL[6])+2.449489742783178*alphaR[4]*fL[5]+4.242640687119286*alphaR[0]*fL[3]+2.449489742783178*(alphaR[2]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[4] = 0.125*(2.449489742783178*(alphaR[4]*fL[15]+alphaR[2]*fL[14]+alphaR[1]*fL[13])+1.414213562373095*alphaR[4]*fL[12]+2.449489742783178*alphaR[0]*fL[10]+1.414213562373095*(alphaR[2]*fL[9]+alphaR[1]*fL[8]+alphaR[0]*fL[4])); 
  incr[5] = 0.125*(2.449489742783178*(alphaR[0]*fL[11]+alphaR[1]*fL[7]+alphaR[2]*fL[6])+1.414213562373095*alphaR[0]*fL[5]+(2.449489742783178*fL[3]+1.414213562373095*fL[0])*alphaR[4]+1.414213562373095*(alphaR[1]*fL[2]+fL[1]*alphaR[2])); 
  incr[6] = -0.125*(4.242640687119286*(alphaR[2]*fL[11]+alphaR[4]*fL[7]+alphaR[0]*fL[6])+2.449489742783178*(alphaR[2]*fL[5]+fL[2]*alphaR[4])+4.242640687119286*alphaR[1]*fL[3]+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[7] = -0.125*(4.242640687119286*(alphaR[1]*fL[11]+alphaR[0]*fL[7]+alphaR[4]*fL[6])+2.449489742783178*(alphaR[1]*fL[5]+fL[1]*alphaR[4])+4.242640687119286*alphaR[2]*fL[3]+2.449489742783178*(alphaR[0]*fL[2]+fL[0]*alphaR[2])); 
  incr[8] = 0.125*(2.449489742783178*(alphaR[2]*fL[15]+alphaR[4]*fL[14]+alphaR[0]*fL[13])+1.414213562373095*alphaR[2]*fL[12]+2.449489742783178*alphaR[1]*fL[10]+1.414213562373095*(alphaR[4]*fL[9]+alphaR[0]*fL[8]+alphaR[1]*fL[4])); 
  incr[9] = 0.125*(2.449489742783178*(alphaR[1]*fL[15]+alphaR[0]*fL[14]+alphaR[4]*fL[13])+1.414213562373095*alphaR[1]*fL[12]+2.449489742783178*alphaR[2]*fL[10]+1.414213562373095*(alphaR[0]*fL[9]+alphaR[4]*fL[8]+alphaR[2]*fL[4])); 
  incr[10] = -0.125*(4.242640687119286*(alphaR[4]*fL[15]+alphaR[2]*fL[14]+alphaR[1]*fL[13])+2.449489742783178*alphaR[4]*fL[12]+4.242640687119286*alphaR[0]*fL[10]+2.449489742783178*(alphaR[2]*fL[9]+alphaR[1]*fL[8]+alphaR[0]*fL[4])); 
  incr[11] = -0.125*(4.242640687119286*(alphaR[0]*fL[11]+alphaR[1]*fL[7]+alphaR[2]*fL[6])+2.449489742783178*alphaR[0]*fL[5]+(4.242640687119286*fL[3]+2.449489742783178*fL[0])*alphaR[4]+2.449489742783178*(alphaR[1]*fL[2]+fL[1]*alphaR[2])); 
  incr[12] = 0.125*(2.449489742783178*(alphaR[0]*fL[15]+alphaR[1]*fL[14]+alphaR[2]*fL[13])+1.414213562373095*alphaR[0]*fL[12]+2.449489742783178*alphaR[4]*fL[10]+1.414213562373095*(alphaR[1]*fL[9]+alphaR[2]*fL[8]+alphaR[4]*fL[4])); 
  incr[13] = -0.125*(4.242640687119286*(alphaR[2]*fL[15]+alphaR[4]*fL[14]+alphaR[0]*fL[13])+2.449489742783178*alphaR[2]*fL[12]+4.242640687119286*alphaR[1]*fL[10]+2.449489742783178*(alphaR[4]*fL[9]+alphaR[0]*fL[8]+alphaR[1]*fL[4])); 
  incr[14] = -0.125*(4.242640687119286*(alphaR[1]*fL[15]+alphaR[0]*fL[14]+alphaR[4]*fL[13])+2.449489742783178*alphaR[1]*fL[12]+4.242640687119286*alphaR[2]*fL[10]+2.449489742783178*(alphaR[0]*fL[9]+alphaR[4]*fL[8]+alphaR[2]*fL[4])); 
  incr[15] = -0.125*(4.242640687119286*(alphaR[0]*fL[15]+alphaR[1]*fL[14]+alphaR[2]*fL[13])+2.449489742783178*alphaR[0]*fL[12]+4.242640687119286*alphaR[4]*fL[10]+2.449489742783178*(alphaR[1]*fL[9]+alphaR[2]*fL[8]+alphaR[4]*fL[4])); 
  incrEmMod[0] = 0.5*(1.732050807568877*fL[3]+fL[0]); 
  incrEmMod[1] = 0.5*(1.732050807568877*fL[6]+fL[1]); 
  incrEmMod[2] = 0.5*(1.732050807568877*fL[7]+fL[2]); 
  incrEmMod[3] = -0.5*(3.0*fL[3]+1.732050807568877*fL[0]); 
  incrEmMod[4] = 0.5*(1.732050807568877*fL[10]+fL[4]); 
  incrEmMod[5] = 0.5*(1.732050807568877*fL[11]+fL[5]); 
  incrEmMod[6] = -0.5*(3.0*fL[6]+1.732050807568877*fL[1]); 
  incrEmMod[7] = -0.5*(3.0*fL[7]+1.732050807568877*fL[2]); 
  incrEmMod[8] = 0.5*(1.732050807568877*fL[13]+fL[8]); 
  incrEmMod[9] = 0.5*(1.732050807568877*fL[14]+fL[9]); 
  incrEmMod[10] = -0.5*(3.0*fL[10]+1.732050807568877*fL[4]); 
  incrEmMod[11] = -0.5*(3.0*fL[11]+1.732050807568877*fL[5]); 
  incrEmMod[12] = 0.5*(1.732050807568877*fL[15]+fL[12]); 
  incrEmMod[13] = -0.5*(3.0*fL[13]+1.732050807568877*fL[8]); 
  incrEmMod[14] = -0.5*(3.0*fL[14]+1.732050807568877*fL[9]); 
  incrEmMod[15] = -0.5*(3.0*fL[15]+1.732050807568877*fL[12]); 
  } else { 
  incrEmMod[0] = -0.5*(1.732050807568877*fR[3]-1.0*fR[0]); 
  incrEmMod[1] = -0.5*(1.732050807568877*fR[6]-1.0*fR[1]); 
  incrEmMod[2] = -0.5*(1.732050807568877*fR[7]-1.0*fR[2]); 
  incrEmMod[3] = 0.5*(3.0*fR[3]-1.732050807568877*fR[0]); 
  incrEmMod[4] = -0.5*(1.732050807568877*fR[10]-1.0*fR[4]); 
  incrEmMod[5] = -0.5*(1.732050807568877*fR[11]-1.0*fR[5]); 
  incrEmMod[6] = 0.5*(3.0*fR[6]-1.732050807568877*fR[1]); 
  incrEmMod[7] = 0.5*(3.0*fR[7]-1.732050807568877*fR[2]); 
  incrEmMod[8] = -0.5*(1.732050807568877*fR[13]-1.0*fR[8]); 
  incrEmMod[9] = -0.5*(1.732050807568877*fR[14]-1.0*fR[9]); 
  incrEmMod[10] = 0.5*(3.0*fR[10]-1.732050807568877*fR[4]); 
  incrEmMod[11] = 0.5*(3.0*fR[11]-1.732050807568877*fR[5]); 
  incrEmMod[12] = -0.5*(1.732050807568877*fR[15]-1.0*fR[12]); 
  incrEmMod[13] = 0.5*(3.0*fR[13]-1.732050807568877*fR[8]); 
  incrEmMod[14] = 0.5*(3.0*fR[14]-1.732050807568877*fR[9]); 
  incrEmMod[15] = 0.5*(3.0*fR[15]-1.732050807568877*fR[12]); 
  incr[0] = -0.125*(2.449489742783178*(alphaR[4]*fR[11]+alphaR[2]*fR[7]+alphaR[1]*fR[6])-1.414213562373095*alphaR[4]*fR[5]+2.449489742783178*alphaR[0]*fR[3]-1.414213562373095*(alphaR[2]*fR[2]+alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.125*(2.449489742783178*(alphaR[2]*fR[11]+alphaR[4]*fR[7]+alphaR[0]*fR[6])-1.414213562373095*(alphaR[2]*fR[5]+fR[2]*alphaR[4])+2.449489742783178*alphaR[1]*fR[3]-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = -0.125*(2.449489742783178*(alphaR[1]*fR[11]+alphaR[0]*fR[7]+alphaR[4]*fR[6])-1.414213562373095*(alphaR[1]*fR[5]+fR[1]*alphaR[4])+2.449489742783178*alphaR[2]*fR[3]-1.414213562373095*(alphaR[0]*fR[2]+fR[0]*alphaR[2])); 
  incr[3] = 0.125*(4.242640687119286*(alphaR[4]*fR[11]+alphaR[2]*fR[7]+alphaR[1]*fR[6])-2.449489742783178*alphaR[4]*fR[5]+4.242640687119286*alphaR[0]*fR[3]-2.449489742783178*(alphaR[2]*fR[2]+alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[4] = -0.125*(2.449489742783178*(alphaR[4]*fR[15]+alphaR[2]*fR[14]+alphaR[1]*fR[13])-1.414213562373095*alphaR[4]*fR[12]+2.449489742783178*alphaR[0]*fR[10]-1.414213562373095*(alphaR[2]*fR[9]+alphaR[1]*fR[8]+alphaR[0]*fR[4])); 
  incr[5] = -0.125*(2.449489742783178*(alphaR[0]*fR[11]+alphaR[1]*fR[7]+alphaR[2]*fR[6])-1.414213562373095*alphaR[0]*fR[5]+(2.449489742783178*fR[3]-1.414213562373095*fR[0])*alphaR[4]-1.414213562373095*(alphaR[1]*fR[2]+fR[1]*alphaR[2])); 
  incr[6] = 0.125*(4.242640687119286*(alphaR[2]*fR[11]+alphaR[4]*fR[7]+alphaR[0]*fR[6])-2.449489742783178*(alphaR[2]*fR[5]+fR[2]*alphaR[4])+4.242640687119286*alphaR[1]*fR[3]-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[7] = 0.125*(4.242640687119286*(alphaR[1]*fR[11]+alphaR[0]*fR[7]+alphaR[4]*fR[6])-2.449489742783178*(alphaR[1]*fR[5]+fR[1]*alphaR[4])+4.242640687119286*alphaR[2]*fR[3]-2.449489742783178*(alphaR[0]*fR[2]+fR[0]*alphaR[2])); 
  incr[8] = -0.125*(2.449489742783178*(alphaR[2]*fR[15]+alphaR[4]*fR[14]+alphaR[0]*fR[13])-1.414213562373095*alphaR[2]*fR[12]+2.449489742783178*alphaR[1]*fR[10]-1.414213562373095*(alphaR[4]*fR[9]+alphaR[0]*fR[8]+alphaR[1]*fR[4])); 
  incr[9] = -0.125*(2.449489742783178*(alphaR[1]*fR[15]+alphaR[0]*fR[14]+alphaR[4]*fR[13])-1.414213562373095*alphaR[1]*fR[12]+2.449489742783178*alphaR[2]*fR[10]-1.414213562373095*(alphaR[0]*fR[9]+alphaR[4]*fR[8]+alphaR[2]*fR[4])); 
  incr[10] = 0.125*(4.242640687119286*(alphaR[4]*fR[15]+alphaR[2]*fR[14]+alphaR[1]*fR[13])-2.449489742783178*alphaR[4]*fR[12]+4.242640687119286*alphaR[0]*fR[10]-2.449489742783178*(alphaR[2]*fR[9]+alphaR[1]*fR[8]+alphaR[0]*fR[4])); 
  incr[11] = 0.125*(4.242640687119286*(alphaR[0]*fR[11]+alphaR[1]*fR[7]+alphaR[2]*fR[6])-2.449489742783178*alphaR[0]*fR[5]+(4.242640687119286*fR[3]-2.449489742783178*fR[0])*alphaR[4]-2.449489742783178*(alphaR[1]*fR[2]+fR[1]*alphaR[2])); 
  incr[12] = -0.125*(2.449489742783178*(alphaR[0]*fR[15]+alphaR[1]*fR[14]+alphaR[2]*fR[13])-1.414213562373095*alphaR[0]*fR[12]+2.449489742783178*alphaR[4]*fR[10]-1.414213562373095*(alphaR[1]*fR[9]+alphaR[2]*fR[8]+alphaR[4]*fR[4])); 
  incr[13] = 0.125*(4.242640687119286*(alphaR[2]*fR[15]+alphaR[4]*fR[14]+alphaR[0]*fR[13])-2.449489742783178*alphaR[2]*fR[12]+4.242640687119286*alphaR[1]*fR[10]-2.449489742783178*(alphaR[4]*fR[9]+alphaR[0]*fR[8]+alphaR[1]*fR[4])); 
  incr[14] = 0.125*(4.242640687119286*(alphaR[1]*fR[15]+alphaR[0]*fR[14]+alphaR[4]*fR[13])-2.449489742783178*alphaR[1]*fR[12]+4.242640687119286*alphaR[2]*fR[10]-2.449489742783178*(alphaR[0]*fR[9]+alphaR[4]*fR[8]+alphaR[2]*fR[4])); 
  incr[15] = 0.125*(4.242640687119286*(alphaR[0]*fR[15]+alphaR[1]*fR[14]+alphaR[2]*fR[13])-2.449489742783178*alphaR[0]*fR[12]+4.242640687119286*alphaR[4]*fR[10]-2.449489742783178*(alphaR[1]*fR[9]+alphaR[2]*fR[8]+alphaR[4]*fR[4])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alphaUpR[4]-0.3535533905932737*(alphaUpR[2]+alphaUpR[1])+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])-0.25*(fR[9]+fR[8])+0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14]+fR[13])+0.4330127018922193*(fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[1]+alphaUpR[0])-0.3535533905932737*(alphaUpR[4]+alphaUpR[2]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14]+fR[13])-0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaUpR[4])+0.3535533905932737*alphaUpR[2]-0.3535533905932737*alphaUpR[1]+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*(fL[14]+fR[13])+0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[4]+alphaUpR[2]+alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]+fR[13]+fL[13]))+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])+0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*(fR[2]+fR[1]+fR[0])+0.25*(fL[2]+fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14]+fR[13])-0.4330127018922193*(fL[15]+fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9]+fR[8]+fL[8])-0.4330127018922193*(fR[7]+fR[6])+0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2]+fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaUpR[4]-0.3535533905932737*(alphaUpR[2]+alphaUpR[1])+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*(fR[5]+fR[4])+0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14]+fR[13])-0.4330127018922193*(fL[14]+fL[13])+0.25*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fR[10])+0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[1]+alphaUpR[0])-0.3535533905932737*(alphaUpR[4]+alphaUpR[2]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14]+fR[13])+0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaUpR[4])+0.3535533905932737*alphaUpR[2]-0.3535533905932737*alphaUpR[1]+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*(fL[14]+fR[13])-0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[4]+alphaUpR[2]+alphaUpR[1]+alphaUpR[0]); 
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

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = 0.7071067811865475*fUp[2]; 
  incrEmMod[3] = -1.224744871391589*fUp[0]; 
  incrEmMod[4] = 0.7071067811865475*fUp[3]; 
  incrEmMod[5] = 0.7071067811865475*fUp[4]; 
  incrEmMod[6] = -1.224744871391589*fUp[1]; 
  incrEmMod[7] = -1.224744871391589*fUp[2]; 
  incrEmMod[8] = 0.7071067811865475*fUp[5]; 
  incrEmMod[9] = 0.7071067811865475*fUp[6]; 
  incrEmMod[10] = -1.224744871391589*fUp[3]; 
  incrEmMod[11] = -1.224744871391589*fUp[4]; 
  incrEmMod[12] = 0.7071067811865475*fUp[7]; 
  incrEmMod[13] = -1.224744871391589*fUp[5]; 
  incrEmMod[14] = -1.224744871391589*fUp[6]; 
  incrEmMod[15] = -1.224744871391589*fUp[7]; 

  incr[0] = 0.25*(alphaR[4]*fUp[4]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.25*(alphaR[2]*fUp[4]+fUp[2]*alphaR[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = 0.25*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[3] = -0.4330127018922193*(alphaR[4]*fUp[4]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[4] = 0.25*(alphaR[4]*fUp[7]+alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+alphaR[0]*fUp[3]); 
  incr[5] = 0.25*(alphaR[0]*fUp[4]+fUp[0]*alphaR[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[6] = -0.4330127018922193*(alphaR[2]*fUp[4]+fUp[2]*alphaR[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[7] = -0.4330127018922193*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[8] = 0.25*(alphaR[2]*fUp[7]+alphaR[4]*fUp[6]+alphaR[0]*fUp[5]+alphaR[1]*fUp[3]); 
  incr[9] = 0.25*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[4]*fUp[5]+alphaR[2]*fUp[3]); 
  incr[10] = -0.4330127018922193*(alphaR[4]*fUp[7]+alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+alphaR[0]*fUp[3]); 
  incr[11] = -0.4330127018922193*(alphaR[0]*fUp[4]+fUp[0]*alphaR[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[12] = 0.25*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]+fUp[3]*alphaR[4]); 
  incr[13] = -0.4330127018922193*(alphaR[2]*fUp[7]+alphaR[4]*fUp[6]+alphaR[0]*fUp[5]+alphaR[1]*fUp[3]); 
  incr[14] = -0.4330127018922193*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[4]*fUp[5]+alphaR[2]*fUp[3]); 
  incr[15] = -0.4330127018922193*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]+fUp[3]*alphaR[4]); 

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

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModR[4] += incrEmMod[4]*rdvpar2R; 
  emModR[5] += incrEmMod[5]*rdvpar2R; 
  emModR[6] += incrEmMod[6]*rdvpar2R; 
  emModR[7] += incrEmMod[7]*rdvpar2R; 
  emModR[8] += incrEmMod[8]*rdvpar2R; 
  emModR[9] += incrEmMod[9]*rdvpar2R; 
  emModR[10] += incrEmMod[10]*rdvpar2R; 
  emModR[11] += incrEmMod[11]*rdvpar2R; 
  emModR[12] += incrEmMod[12]*rdvpar2R; 
  emModR[13] += incrEmMod[13]*rdvpar2R; 
  emModR[14] += incrEmMod[14]*rdvpar2R; 
  emModR[15] += incrEmMod[15]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += -1.0*incrEmMod[2]*rdvpar2L; 
  emModL[3] += incrEmMod[3]*rdvpar2L; 
  emModL[4] += -1.0*incrEmMod[4]*rdvpar2L; 
  emModL[5] += -1.0*incrEmMod[5]*rdvpar2L; 
  emModL[6] += incrEmMod[6]*rdvpar2L; 
  emModL[7] += incrEmMod[7]*rdvpar2L; 
  emModL[8] += -1.0*incrEmMod[8]*rdvpar2L; 
  emModL[9] += -1.0*incrEmMod[9]*rdvpar2L; 
  emModL[10] += incrEmMod[10]*rdvpar2L; 
  emModL[11] += incrEmMod[11]*rdvpar2L; 
  emModL[12] += -1.0*incrEmMod[12]*rdvpar2L; 
  emModL[13] += incrEmMod[13]*rdvpar2L; 
  emModL[14] += incrEmMod[14]*rdvpar2L; 
  emModL[15] += incrEmMod[15]*rdvpar2L; 

  return std::abs(alphaSurfAvgR);
} 
double EmGyrokineticGenGeoSurf2x2vSerStep2_vpar_P1_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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
  BstarXdBmagR[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2R; 
  BstarXdBmagR[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2R; 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2R; 
  BstarYdBmagR[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2R; 

  double alphaR[8]; 
  alphaR[0] = -(1.414213562373095*dApardt[0]*q_)/m_; 
  alphaR[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  alphaR[2] = -(1.414213562373095*dApardt[2]*q_)/m_; 
  alphaR[4] = -(1.414213562373095*dApardt[3]*q_)/m_; 

  double alphaUpR[8]; 
  alphaUpR[0] = -(0.1767766952966368*(1.732050807568877*(BstarYdBmagR[0]*hamilR[2]*rdy2R+BstarXdBmagR[0]*hamilR[1]*rdx2R)+8.0*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = -(0.1767766952966368*(1.732050807568877*(BstarYdBmagR[0]*hamilR[5]*rdy2R+BstarXdBmagR[1]*hamilR[1]*rdx2R)+8.0*dApardtPrev[1]*q_))/m_; 
  alphaUpR[2] = -(0.1767766952966368*(1.732050807568877*(BstarYdBmagR[2]*hamilR[2]*rdy2R+BstarXdBmagR[0]*hamilR[5]*rdx2R)+8.0*dApardtPrev[2]*q_))/m_; 
  alphaUpR[4] = -(0.1767766952966368*(1.732050807568877*hamilR[5]*(BstarYdBmagR[2]*rdy2R+BstarXdBmagR[1]*rdx2R)+8.0*dApardtPrev[3]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(1.732050807568877*BstarYdBmagR[0]*hamilR[2]*rdy2R+1.732050807568877*BstarXdBmagR[0]*hamilR[1]*rdx2R+8.0*dApardtPrev[0]*q_))/m_; 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.125*(2.449489742783178*(alphaR[4]*fL[11]+alphaR[2]*fL[7]+alphaR[1]*fL[6])+1.414213562373095*alphaR[4]*fL[5]+2.449489742783178*alphaR[0]*fL[3]+1.414213562373095*(alphaR[2]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.125*(2.449489742783178*(alphaR[2]*fL[11]+alphaR[4]*fL[7]+alphaR[0]*fL[6])+1.414213562373095*(alphaR[2]*fL[5]+fL[2]*alphaR[4])+2.449489742783178*alphaR[1]*fL[3]+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = 0.125*(2.449489742783178*(alphaR[1]*fL[11]+alphaR[0]*fL[7]+alphaR[4]*fL[6])+1.414213562373095*(alphaR[1]*fL[5]+fL[1]*alphaR[4])+2.449489742783178*alphaR[2]*fL[3]+1.414213562373095*(alphaR[0]*fL[2]+fL[0]*alphaR[2])); 
  incr[3] = -0.125*(4.242640687119286*(alphaR[4]*fL[11]+alphaR[2]*fL[7]+alphaR[1]*fL[6])+2.449489742783178*alphaR[4]*fL[5]+4.242640687119286*alphaR[0]*fL[3]+2.449489742783178*(alphaR[2]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[4] = 0.125*(2.449489742783178*(alphaR[4]*fL[15]+alphaR[2]*fL[14]+alphaR[1]*fL[13])+1.414213562373095*alphaR[4]*fL[12]+2.449489742783178*alphaR[0]*fL[10]+1.414213562373095*(alphaR[2]*fL[9]+alphaR[1]*fL[8]+alphaR[0]*fL[4])); 
  incr[5] = 0.125*(2.449489742783178*(alphaR[0]*fL[11]+alphaR[1]*fL[7]+alphaR[2]*fL[6])+1.414213562373095*alphaR[0]*fL[5]+(2.449489742783178*fL[3]+1.414213562373095*fL[0])*alphaR[4]+1.414213562373095*(alphaR[1]*fL[2]+fL[1]*alphaR[2])); 
  incr[6] = -0.125*(4.242640687119286*(alphaR[2]*fL[11]+alphaR[4]*fL[7]+alphaR[0]*fL[6])+2.449489742783178*(alphaR[2]*fL[5]+fL[2]*alphaR[4])+4.242640687119286*alphaR[1]*fL[3]+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[7] = -0.125*(4.242640687119286*(alphaR[1]*fL[11]+alphaR[0]*fL[7]+alphaR[4]*fL[6])+2.449489742783178*(alphaR[1]*fL[5]+fL[1]*alphaR[4])+4.242640687119286*alphaR[2]*fL[3]+2.449489742783178*(alphaR[0]*fL[2]+fL[0]*alphaR[2])); 
  incr[8] = 0.125*(2.449489742783178*(alphaR[2]*fL[15]+alphaR[4]*fL[14]+alphaR[0]*fL[13])+1.414213562373095*alphaR[2]*fL[12]+2.449489742783178*alphaR[1]*fL[10]+1.414213562373095*(alphaR[4]*fL[9]+alphaR[0]*fL[8]+alphaR[1]*fL[4])); 
  incr[9] = 0.125*(2.449489742783178*(alphaR[1]*fL[15]+alphaR[0]*fL[14]+alphaR[4]*fL[13])+1.414213562373095*alphaR[1]*fL[12]+2.449489742783178*alphaR[2]*fL[10]+1.414213562373095*(alphaR[0]*fL[9]+alphaR[4]*fL[8]+alphaR[2]*fL[4])); 
  incr[10] = -0.125*(4.242640687119286*(alphaR[4]*fL[15]+alphaR[2]*fL[14]+alphaR[1]*fL[13])+2.449489742783178*alphaR[4]*fL[12]+4.242640687119286*alphaR[0]*fL[10]+2.449489742783178*(alphaR[2]*fL[9]+alphaR[1]*fL[8]+alphaR[0]*fL[4])); 
  incr[11] = -0.125*(4.242640687119286*(alphaR[0]*fL[11]+alphaR[1]*fL[7]+alphaR[2]*fL[6])+2.449489742783178*alphaR[0]*fL[5]+(4.242640687119286*fL[3]+2.449489742783178*fL[0])*alphaR[4]+2.449489742783178*(alphaR[1]*fL[2]+fL[1]*alphaR[2])); 
  incr[12] = 0.125*(2.449489742783178*(alphaR[0]*fL[15]+alphaR[1]*fL[14]+alphaR[2]*fL[13])+1.414213562373095*alphaR[0]*fL[12]+2.449489742783178*alphaR[4]*fL[10]+1.414213562373095*(alphaR[1]*fL[9]+alphaR[2]*fL[8]+alphaR[4]*fL[4])); 
  incr[13] = -0.125*(4.242640687119286*(alphaR[2]*fL[15]+alphaR[4]*fL[14]+alphaR[0]*fL[13])+2.449489742783178*alphaR[2]*fL[12]+4.242640687119286*alphaR[1]*fL[10]+2.449489742783178*(alphaR[4]*fL[9]+alphaR[0]*fL[8]+alphaR[1]*fL[4])); 
  incr[14] = -0.125*(4.242640687119286*(alphaR[1]*fL[15]+alphaR[0]*fL[14]+alphaR[4]*fL[13])+2.449489742783178*alphaR[1]*fL[12]+4.242640687119286*alphaR[2]*fL[10]+2.449489742783178*(alphaR[0]*fL[9]+alphaR[4]*fL[8]+alphaR[2]*fL[4])); 
  incr[15] = -0.125*(4.242640687119286*(alphaR[0]*fL[15]+alphaR[1]*fL[14]+alphaR[2]*fL[13])+2.449489742783178*alphaR[0]*fL[12]+4.242640687119286*alphaR[4]*fL[10]+2.449489742783178*(alphaR[1]*fL[9]+alphaR[2]*fL[8]+alphaR[4]*fL[4])); 
  } else { 
  incr[0] = -0.125*(2.449489742783178*(alphaR[4]*fR[11]+alphaR[2]*fR[7]+alphaR[1]*fR[6])-1.414213562373095*alphaR[4]*fR[5]+2.449489742783178*alphaR[0]*fR[3]-1.414213562373095*(alphaR[2]*fR[2]+alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.125*(2.449489742783178*(alphaR[2]*fR[11]+alphaR[4]*fR[7]+alphaR[0]*fR[6])-1.414213562373095*(alphaR[2]*fR[5]+fR[2]*alphaR[4])+2.449489742783178*alphaR[1]*fR[3]-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = -0.125*(2.449489742783178*(alphaR[1]*fR[11]+alphaR[0]*fR[7]+alphaR[4]*fR[6])-1.414213562373095*(alphaR[1]*fR[5]+fR[1]*alphaR[4])+2.449489742783178*alphaR[2]*fR[3]-1.414213562373095*(alphaR[0]*fR[2]+fR[0]*alphaR[2])); 
  incr[3] = 0.125*(4.242640687119286*(alphaR[4]*fR[11]+alphaR[2]*fR[7]+alphaR[1]*fR[6])-2.449489742783178*alphaR[4]*fR[5]+4.242640687119286*alphaR[0]*fR[3]-2.449489742783178*(alphaR[2]*fR[2]+alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[4] = -0.125*(2.449489742783178*(alphaR[4]*fR[15]+alphaR[2]*fR[14]+alphaR[1]*fR[13])-1.414213562373095*alphaR[4]*fR[12]+2.449489742783178*alphaR[0]*fR[10]-1.414213562373095*(alphaR[2]*fR[9]+alphaR[1]*fR[8]+alphaR[0]*fR[4])); 
  incr[5] = -0.125*(2.449489742783178*(alphaR[0]*fR[11]+alphaR[1]*fR[7]+alphaR[2]*fR[6])-1.414213562373095*alphaR[0]*fR[5]+(2.449489742783178*fR[3]-1.414213562373095*fR[0])*alphaR[4]-1.414213562373095*(alphaR[1]*fR[2]+fR[1]*alphaR[2])); 
  incr[6] = 0.125*(4.242640687119286*(alphaR[2]*fR[11]+alphaR[4]*fR[7]+alphaR[0]*fR[6])-2.449489742783178*(alphaR[2]*fR[5]+fR[2]*alphaR[4])+4.242640687119286*alphaR[1]*fR[3]-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[7] = 0.125*(4.242640687119286*(alphaR[1]*fR[11]+alphaR[0]*fR[7]+alphaR[4]*fR[6])-2.449489742783178*(alphaR[1]*fR[5]+fR[1]*alphaR[4])+4.242640687119286*alphaR[2]*fR[3]-2.449489742783178*(alphaR[0]*fR[2]+fR[0]*alphaR[2])); 
  incr[8] = -0.125*(2.449489742783178*(alphaR[2]*fR[15]+alphaR[4]*fR[14]+alphaR[0]*fR[13])-1.414213562373095*alphaR[2]*fR[12]+2.449489742783178*alphaR[1]*fR[10]-1.414213562373095*(alphaR[4]*fR[9]+alphaR[0]*fR[8]+alphaR[1]*fR[4])); 
  incr[9] = -0.125*(2.449489742783178*(alphaR[1]*fR[15]+alphaR[0]*fR[14]+alphaR[4]*fR[13])-1.414213562373095*alphaR[1]*fR[12]+2.449489742783178*alphaR[2]*fR[10]-1.414213562373095*(alphaR[0]*fR[9]+alphaR[4]*fR[8]+alphaR[2]*fR[4])); 
  incr[10] = 0.125*(4.242640687119286*(alphaR[4]*fR[15]+alphaR[2]*fR[14]+alphaR[1]*fR[13])-2.449489742783178*alphaR[4]*fR[12]+4.242640687119286*alphaR[0]*fR[10]-2.449489742783178*(alphaR[2]*fR[9]+alphaR[1]*fR[8]+alphaR[0]*fR[4])); 
  incr[11] = 0.125*(4.242640687119286*(alphaR[0]*fR[11]+alphaR[1]*fR[7]+alphaR[2]*fR[6])-2.449489742783178*alphaR[0]*fR[5]+(4.242640687119286*fR[3]-2.449489742783178*fR[0])*alphaR[4]-2.449489742783178*(alphaR[1]*fR[2]+fR[1]*alphaR[2])); 
  incr[12] = -0.125*(2.449489742783178*(alphaR[0]*fR[15]+alphaR[1]*fR[14]+alphaR[2]*fR[13])-1.414213562373095*alphaR[0]*fR[12]+2.449489742783178*alphaR[4]*fR[10]-1.414213562373095*(alphaR[1]*fR[9]+alphaR[2]*fR[8]+alphaR[4]*fR[4])); 
  incr[13] = 0.125*(4.242640687119286*(alphaR[2]*fR[15]+alphaR[4]*fR[14]+alphaR[0]*fR[13])-2.449489742783178*alphaR[2]*fR[12]+4.242640687119286*alphaR[1]*fR[10]-2.449489742783178*(alphaR[4]*fR[9]+alphaR[0]*fR[8]+alphaR[1]*fR[4])); 
  incr[14] = 0.125*(4.242640687119286*(alphaR[1]*fR[15]+alphaR[0]*fR[14]+alphaR[4]*fR[13])-2.449489742783178*alphaR[1]*fR[12]+4.242640687119286*alphaR[2]*fR[10]-2.449489742783178*(alphaR[0]*fR[9]+alphaR[4]*fR[8]+alphaR[2]*fR[4])); 
  incr[15] = 0.125*(4.242640687119286*(alphaR[0]*fR[15]+alphaR[1]*fR[14]+alphaR[2]*fR[13])-2.449489742783178*alphaR[0]*fR[12]+4.242640687119286*alphaR[4]*fR[10]-2.449489742783178*(alphaR[1]*fR[9]+alphaR[2]*fR[8]+alphaR[4]*fR[4])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alphaUpR[4]-0.3535533905932737*(alphaUpR[2]+alphaUpR[1])+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])-0.25*(fR[9]+fR[8])+0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14]+fR[13])+0.4330127018922193*(fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[1]+alphaUpR[0])-0.3535533905932737*(alphaUpR[4]+alphaUpR[2]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14]+fR[13])-0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaUpR[4])+0.3535533905932737*alphaUpR[2]-0.3535533905932737*alphaUpR[1]+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*(fL[14]+fR[13])+0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[4]+alphaUpR[2]+alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]+fR[13]+fL[13]))+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])+0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*(fR[2]+fR[1]+fR[0])+0.25*(fL[2]+fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14]+fR[13])-0.4330127018922193*(fL[15]+fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9]+fR[8]+fL[8])-0.4330127018922193*(fR[7]+fR[6])+0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2]+fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaUpR[4]-0.3535533905932737*(alphaUpR[2]+alphaUpR[1])+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*(fR[5]+fR[4])+0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14]+fR[13])-0.4330127018922193*(fL[14]+fL[13])+0.25*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fR[10])+0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[1]+alphaUpR[0])-0.3535533905932737*(alphaUpR[4]+alphaUpR[2]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14]+fR[13])+0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaUpR[4])+0.3535533905932737*alphaUpR[2]-0.3535533905932737*alphaUpR[1]+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*(fL[14]+fR[13])-0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[4]+alphaUpR[2]+alphaUpR[1]+alphaUpR[0]); 
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

  incr[0] = 0.25*(alphaR[4]*fUp[4]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.25*(alphaR[2]*fUp[4]+fUp[2]*alphaR[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = 0.25*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[3] = -0.4330127018922193*(alphaR[4]*fUp[4]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[4] = 0.25*(alphaR[4]*fUp[7]+alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+alphaR[0]*fUp[3]); 
  incr[5] = 0.25*(alphaR[0]*fUp[4]+fUp[0]*alphaR[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[6] = -0.4330127018922193*(alphaR[2]*fUp[4]+fUp[2]*alphaR[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[7] = -0.4330127018922193*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[8] = 0.25*(alphaR[2]*fUp[7]+alphaR[4]*fUp[6]+alphaR[0]*fUp[5]+alphaR[1]*fUp[3]); 
  incr[9] = 0.25*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[4]*fUp[5]+alphaR[2]*fUp[3]); 
  incr[10] = -0.4330127018922193*(alphaR[4]*fUp[7]+alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+alphaR[0]*fUp[3]); 
  incr[11] = -0.4330127018922193*(alphaR[0]*fUp[4]+fUp[0]*alphaR[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[12] = 0.25*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]+fUp[3]*alphaR[4]); 
  incr[13] = -0.4330127018922193*(alphaR[2]*fUp[7]+alphaR[4]*fUp[6]+alphaR[0]*fUp[5]+alphaR[1]*fUp[3]); 
  incr[14] = -0.4330127018922193*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[4]*fUp[5]+alphaR[2]*fUp[3]); 
  incr[15] = -0.4330127018922193*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]+fUp[3]*alphaR[4]); 

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
double EmGyrokineticGenGeoSurf2x2vSer_xL_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
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

  double hamilL[16]; 
  hamilL[0] = (0.6666666666666666*(3.0*rdvpar2SqL*(m_*wvparSqL+bmag[0]*wmuL+phi[0]*q_)+m_))/rdvpar2SqL; 
  hamilL[1] = 2.0*(bmag[1]*wmuL+phi[1]*q_); 
  hamilL[2] = 2.0*phi[2]*q_; 
  hamilL[3] = (2.309401076758503*m_*wvparL)/rdvpar2L; 
  hamilL[4] = (1.154700538379252*bmag[0])/rdmu2L; 
  hamilL[5] = 2.0*phi[3]*q_; 
  hamilL[8] = (1.154700538379252*bmag[1])/rdmu2L; 

  double BstarXdBmagL[16]; 
  BstarXdBmagL[0] = 0.8660254037844386*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])*rdy2L; 
  BstarXdBmagL[1] = 0.1732050807568877*((9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])*rdy2L; 

  double alphaL[8]; 
  alphaL[0] = -(0.125*((((12.72792206135786*b_z[1]+7.348469228349534*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(7.348469228349534*b_z[1]+4.242640687119286*b_z[0]))*hamilL[5]+((7.348469228349534*b_z[1]+4.242640687119286*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(4.242640687119286*b_z[1]+2.449489742783178*b_z[0]))*hamilL[2])*m_*rdy2L+((-4.242640687119286*BstarXdBmagL[1])-2.449489742783178*BstarXdBmagL[0])*hamilL[3]*q_*rdvpar2L))/(m_*q_); 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgL = -(0.03125*((((9.0*b_z[1]+5.196152422706631*b_z[0])*jacobTotInv[1]+5.196152422706631*jacobTotInv[0]*b_z[1]+3.0*b_z[0]*jacobTotInv[0])*hamilL[5]+((5.196152422706631*b_z[1]+3.0*b_z[0])*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1]+1.732050807568877*b_z[0]*jacobTotInv[0])*hamilL[2])*m_*rdy2L+((-3.0*BstarXdBmagL[1])-1.732050807568877*BstarXdBmagL[0])*hamilL[3]*q_*rdvpar2L))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgL>0) { 
  incr[0] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[1]+fL[0]); 
  incr[1] = -0.1767766952966368*alphaL[0]*(3.0*fL[1]+1.732050807568877*fL[0]); 
  incr[2] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[5]+fL[2]); 
  incr[3] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[6]+fL[3]); 
  incr[4] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[8]+fL[4]); 
  incr[5] = -0.1767766952966368*alphaL[0]*(3.0*fL[5]+1.732050807568877*fL[2]); 
  incr[6] = -0.1767766952966368*alphaL[0]*(3.0*fL[6]+1.732050807568877*fL[3]); 
  incr[7] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[11]+fL[7]); 
  incr[8] = -0.1767766952966368*alphaL[0]*(3.0*fL[8]+1.732050807568877*fL[4]); 
  incr[9] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[12]+fL[9]); 
  incr[10] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[13]+fL[10]); 
  incr[11] = -0.1767766952966368*alphaL[0]*(3.0*fL[11]+1.732050807568877*fL[7]); 
  incr[12] = -0.1767766952966368*alphaL[0]*(3.0*fL[12]+1.732050807568877*fL[9]); 
  incr[13] = -0.1767766952966368*alphaL[0]*(3.0*fL[13]+1.732050807568877*fL[10]); 
  incr[14] = 0.1767766952966368*alphaL[0]*(1.732050807568877*fL[15]+fL[14]); 
  incr[15] = -0.1767766952966368*alphaL[0]*(3.0*fL[15]+1.732050807568877*fL[14]); 
  } else { 
  incr[0] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[1]-1.0*fR[0]); 
  incr[1] = 0.1767766952966368*alphaL[0]*(3.0*fR[1]-1.732050807568877*fR[0]); 
  incr[2] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[5]-1.0*fR[2]); 
  incr[3] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[6]-1.0*fR[3]); 
  incr[4] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[8]-1.0*fR[4]); 
  incr[5] = 0.1767766952966368*alphaL[0]*(3.0*fR[5]-1.732050807568877*fR[2]); 
  incr[6] = 0.1767766952966368*alphaL[0]*(3.0*fR[6]-1.732050807568877*fR[3]); 
  incr[7] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[11]-1.0*fR[7]); 
  incr[8] = 0.1767766952966368*alphaL[0]*(3.0*fR[8]-1.732050807568877*fR[4]); 
  incr[9] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[12]-1.0*fR[9]); 
  incr[10] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[13]-1.0*fR[10]); 
  incr[11] = 0.1767766952966368*alphaL[0]*(3.0*fR[11]-1.732050807568877*fR[7]); 
  incr[12] = 0.1767766952966368*alphaL[0]*(3.0*fR[12]-1.732050807568877*fR[9]); 
  incr[13] = 0.1767766952966368*alphaL[0]*(3.0*fR[13]-1.732050807568877*fR[10]); 
  incr[14] = -0.1767766952966368*alphaL[0]*(1.732050807568877*fR[15]-1.0*fR[14]); 
  incr[15] = 0.1767766952966368*alphaL[0]*(3.0*fR[15]-1.732050807568877*fR[14]); 
  }
#elif upwindType == QUAD 
  double alphaOrdL;
  double fUpOrd[8];
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12]+fR[11]+fL[11])-0.25*(fR[10]+fR[9])+0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*(fR[4]+fR[3]+fR[2])-0.25*(fL[4]+fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fR[12]+fR[11])+0.4330127018922193*(fL[13]+fL[12]+fL[11])+0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4]+fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*(fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3]+fR[2])+0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fR[5])+0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])+0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3]+fR[2])-0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*(fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[0]; 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12]+fR[11]+fL[11])-0.25*(fR[10]+fR[9])+0.25*(fL[10]+fL[9])+0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])-0.25*(fR[4]+fR[3]+fR[2])+0.25*(fL[4]+fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fR[12]+fR[11])+0.4330127018922193*(fL[13]+fL[12]+fL[11])+0.25*(fR[10]+fL[10]+fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fR[5])+0.4330127018922193*(fL[6]+fL[5])+0.25*(fR[4]+fL[4]+fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*alphaL[0]*fUp[0]; 
  incr[1] = -0.4330127018922193*alphaL[0]*fUp[0]; 
  incr[2] = 0.25*alphaL[0]*fUp[1]; 
  incr[3] = 0.25*alphaL[0]*fUp[2]; 
  incr[4] = 0.25*alphaL[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alphaL[0]*fUp[1]; 
  incr[6] = -0.4330127018922193*alphaL[0]*fUp[2]; 
  incr[7] = 0.25*alphaL[0]*fUp[4]; 
  incr[8] = -0.4330127018922193*alphaL[0]*fUp[3]; 
  incr[9] = 0.25*alphaL[0]*fUp[5]; 
  incr[10] = 0.25*alphaL[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alphaL[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alphaL[0]*fUp[5]; 
  incr[13] = -0.4330127018922193*alphaL[0]*fUp[6]; 
  incr[14] = 0.25*alphaL[0]*fUp[7]; 
  incr[15] = -0.4330127018922193*alphaL[0]*fUp[7]; 

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

  return std::abs(alphaSurfAvgL);
} 
double EmGyrokineticGenGeoSurf2x2vSer_xR_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
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
  BstarXdBmagR[0] = 0.8660254037844386*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])*rdy2R; 
  BstarXdBmagR[1] = 0.1732050807568877*((9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])*rdy2R; 

  double alphaR[8]; 
  alphaR[0] = (0.125*((((12.72792206135786*b_z[1]-7.348469228349534*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(4.242640687119286*b_z[0]-7.348469228349534*b_z[1]))*hamilR[5]+((4.242640687119286*b_z[0]-7.348469228349534*b_z[1])*jacobTotInv[1]+jacobTotInv[0]*(4.242640687119286*b_z[1]-2.449489742783178*b_z[0]))*hamilR[2])*m_*rdy2R+(2.449489742783178*BstarXdBmagR[0]-4.242640687119286*BstarXdBmagR[1])*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*((((9.0*b_z[1]-5.196152422706631*b_z[0])*jacobTotInv[1]-5.196152422706631*jacobTotInv[0]*b_z[1]+3.0*b_z[0]*jacobTotInv[0])*hamilR[5]+((3.0*b_z[0]-5.196152422706631*b_z[1])*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1]-1.732050807568877*b_z[0]*jacobTotInv[0])*hamilR[2])*m_*rdy2R+(1.732050807568877*BstarXdBmagR[0]-3.0*BstarXdBmagR[1])*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[1]+fL[0]); 
  incr[1] = -0.1767766952966368*alphaR[0]*(3.0*fL[1]+1.732050807568877*fL[0]); 
  incr[2] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[5]+fL[2]); 
  incr[3] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[6]+fL[3]); 
  incr[4] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[8]+fL[4]); 
  incr[5] = -0.1767766952966368*alphaR[0]*(3.0*fL[5]+1.732050807568877*fL[2]); 
  incr[6] = -0.1767766952966368*alphaR[0]*(3.0*fL[6]+1.732050807568877*fL[3]); 
  incr[7] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[11]+fL[7]); 
  incr[8] = -0.1767766952966368*alphaR[0]*(3.0*fL[8]+1.732050807568877*fL[4]); 
  incr[9] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[12]+fL[9]); 
  incr[10] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[13]+fL[10]); 
  incr[11] = -0.1767766952966368*alphaR[0]*(3.0*fL[11]+1.732050807568877*fL[7]); 
  incr[12] = -0.1767766952966368*alphaR[0]*(3.0*fL[12]+1.732050807568877*fL[9]); 
  incr[13] = -0.1767766952966368*alphaR[0]*(3.0*fL[13]+1.732050807568877*fL[10]); 
  incr[14] = 0.1767766952966368*alphaR[0]*(1.732050807568877*fL[15]+fL[14]); 
  incr[15] = -0.1767766952966368*alphaR[0]*(3.0*fL[15]+1.732050807568877*fL[14]); 
  } else { 
  incr[0] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[1]-1.0*fR[0]); 
  incr[1] = 0.1767766952966368*alphaR[0]*(3.0*fR[1]-1.732050807568877*fR[0]); 
  incr[2] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[5]-1.0*fR[2]); 
  incr[3] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[6]-1.0*fR[3]); 
  incr[4] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[8]-1.0*fR[4]); 
  incr[5] = 0.1767766952966368*alphaR[0]*(3.0*fR[5]-1.732050807568877*fR[2]); 
  incr[6] = 0.1767766952966368*alphaR[0]*(3.0*fR[6]-1.732050807568877*fR[3]); 
  incr[7] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[11]-1.0*fR[7]); 
  incr[8] = 0.1767766952966368*alphaR[0]*(3.0*fR[8]-1.732050807568877*fR[4]); 
  incr[9] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[12]-1.0*fR[9]); 
  incr[10] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[13]-1.0*fR[10]); 
  incr[11] = 0.1767766952966368*alphaR[0]*(3.0*fR[11]-1.732050807568877*fR[7]); 
  incr[12] = 0.1767766952966368*alphaR[0]*(3.0*fR[12]-1.732050807568877*fR[9]); 
  incr[13] = 0.1767766952966368*alphaR[0]*(3.0*fR[13]-1.732050807568877*fR[10]); 
  incr[14] = -0.1767766952966368*alphaR[0]*(1.732050807568877*fR[15]-1.0*fR[14]); 
  incr[15] = 0.1767766952966368*alphaR[0]*(3.0*fR[15]-1.732050807568877*fR[14]); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12]+fR[11]+fL[11])-0.25*(fR[10]+fR[9])+0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*(fR[4]+fR[3]+fR[2])-0.25*(fL[4]+fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fR[12]+fR[11])+0.4330127018922193*(fL[13]+fL[12]+fL[11])+0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4]+fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]-0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*(fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])-0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3]+fR[2])+0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])+0.4330127018922193*fR[8]-0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fR[5])+0.4330127018922193*(fL[6]+fL[5])-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.25*fR[14]+0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13]+fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*(fR[10]+fR[9])-0.25*(fL[10]+fL[9])+0.4330127018922193*(fR[8]+fL[8])-0.25*fR[7]+0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6]+fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3]+fR[2])-0.25*(fL[3]+fL[2])+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*fL[15]+0.25*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fR[12])-0.4330127018922193*(fL[13]+fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10]+fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]+0.25*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fR[5])-0.4330127018922193*(fL[6]+fL[5])+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3]+fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]-0.4330127018922193*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*(fL[10]+fR[9])+0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]-0.4330127018922193*(fR[6]+fL[6])+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*(fL[3]+fR[2])+0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])+0.4330127018922193*fR[13]-0.4330127018922193*(fL[13]+fR[12])+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])+0.4330127018922193*fR[6]-0.4330127018922193*(fL[6]+fR[5])+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])+0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.25*fR[14]-0.25*fL[14]+0.4330127018922193*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*(fL[10]+fR[9])-0.25*fL[9]+0.4330127018922193*(fR[8]+fL[8])+0.25*fR[7]-0.25*fL[7]+0.4330127018922193*(fR[6]+fL[6])-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3]+fR[2])-0.25*fL[2]+0.4330127018922193*(fR[1]+fL[1])-0.25*fR[0]+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*fL[15]-0.25*(fR[14]+fL[14])-0.4330127018922193*fR[13]+0.4330127018922193*(fL[13]+fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.25*(fR[9]+fL[9])-0.4330127018922193*fR[8]+0.4330127018922193*fL[8]-0.25*(fR[7]+fL[7])-0.4330127018922193*fR[6]+0.4330127018922193*(fL[6]+fR[5])-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.25*(fR[2]+fL[2])-0.4330127018922193*fR[1]+0.4330127018922193*fL[1]+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaR[0]; 
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

  incr[0] = 0.25*alphaR[0]*fUp[0]; 
  incr[1] = -0.4330127018922193*alphaR[0]*fUp[0]; 
  incr[2] = 0.25*alphaR[0]*fUp[1]; 
  incr[3] = 0.25*alphaR[0]*fUp[2]; 
  incr[4] = 0.25*alphaR[0]*fUp[3]; 
  incr[5] = -0.4330127018922193*alphaR[0]*fUp[1]; 
  incr[6] = -0.4330127018922193*alphaR[0]*fUp[2]; 
  incr[7] = 0.25*alphaR[0]*fUp[4]; 
  incr[8] = -0.4330127018922193*alphaR[0]*fUp[3]; 
  incr[9] = 0.25*alphaR[0]*fUp[5]; 
  incr[10] = 0.25*alphaR[0]*fUp[6]; 
  incr[11] = -0.4330127018922193*alphaR[0]*fUp[4]; 
  incr[12] = -0.4330127018922193*alphaR[0]*fUp[5]; 
  incr[13] = -0.4330127018922193*alphaR[0]*fUp[6]; 
  incr[14] = 0.25*alphaR[0]*fUp[7]; 
  incr[15] = -0.4330127018922193*alphaR[0]*fUp[7]; 

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
double EmGyrokineticGenGeoSurf2x2vSer_yL_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
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

  double hamilL[16]; 
  hamilL[0] = (0.6666666666666666*(3.0*rdvpar2SqL*(m_*wvparSqL+bmag[0]*wmuL+phi[0]*q_)+m_))/rdvpar2SqL; 
  hamilL[1] = 2.0*(bmag[1]*wmuL+phi[1]*q_); 
  hamilL[2] = 2.0*phi[2]*q_; 
  hamilL[3] = (2.309401076758503*m_*wvparL)/rdvpar2L; 
  hamilL[4] = (1.154700538379252*bmag[0])/rdmu2L; 
  hamilL[5] = 2.0*phi[3]*q_; 
  hamilL[8] = (1.154700538379252*bmag[1])/rdmu2L; 

  double BstarYdBmagL[16]; 
  BstarYdBmagL[0] = -(0.8660254037844386*rdx2L*(2.0*jacobTotInv[0]*b_z[1]*m_*wvparL+(2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*q_))/q_; 
  BstarYdBmagL[1] = -(0.8660254037844386*rdx2L*(2.0*b_z[1]*jacobTotInv[1]*m_*wvparL+((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])*q_))/q_; 
  BstarYdBmagL[2] = -0.8660254037844386*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])*rdx2L; 
  BstarYdBmagL[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2L)/(q_*rdvpar2L); 
  BstarYdBmagL[5] = -0.8660254037844386*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])*rdx2L; 
  BstarYdBmagL[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2L)/(q_*rdvpar2L); 

  double alphaL[8]; 
  alphaL[0] = (0.1767766952966368*((b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*(3.0*hamilL[5]+1.732050807568877*hamilL[1])*m_*rdx2L+(3.0*BstarYdBmagL[2]+1.732050807568877*BstarYdBmagL[0])*hamilL[3]*q_*rdvpar2L))/(m_*q_); 
  alphaL[1] = (0.1767766952966368*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*(3.0*hamilL[5]+1.732050807568877*hamilL[1])*m_*rdx2L+hamilL[3]*(3.0*BstarYdBmagL[5]+1.732050807568877*BstarYdBmagL[1])*q_*rdvpar2L))/(m_*q_); 
  alphaL[2] = (0.3061862178478971*BstarYdBmagL[3]*hamilL[3]*rdvpar2L)/m_; 
  alphaL[3] = (0.3061862178478971*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamilL[8]*rdx2L)/q_; 
  alphaL[4] = (0.3061862178478971*hamilL[3]*BstarYdBmagL[6]*rdvpar2L)/m_; 
  alphaL[5] = (0.3061862178478971*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamilL[8]*rdx2L)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgL = (0.03125*(((3.0*b_z[1]*jacobTotInv[1]+3.0*b_z[0]*jacobTotInv[0])*hamilL[5]+1.732050807568877*b_z[1]*hamilL[1]*jacobTotInv[1]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilL[1])*m_*rdx2L+(3.0*BstarYdBmagL[2]+1.732050807568877*BstarYdBmagL[0])*hamilL[3]*q_*rdvpar2L))/(m_*q_); 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgL>0) { 
  incr[0] = 0.125*(2.449489742783178*(alphaL[5]*fL[12]+alphaL[4]*fL[11]+alphaL[3]*fL[9])+1.414213562373095*alphaL[5]*fL[8]+2.449489742783178*alphaL[2]*fL[7]+1.414213562373095*alphaL[4]*fL[6]+2.449489742783178*alphaL[1]*fL[5]+1.414213562373095*(alphaL[3]*fL[4]+alphaL[2]*fL[3])+2.449489742783178*alphaL[0]*fL[2]+1.414213562373095*(alphaL[1]*fL[1]+alphaL[0]*fL[0])); 
  incr[1] = 0.125*(2.449489742783178*(alphaL[3]*fL[12]+alphaL[2]*fL[11]+alphaL[5]*fL[9])+1.414213562373095*alphaL[3]*fL[8]+2.449489742783178*alphaL[4]*fL[7]+1.414213562373095*alphaL[2]*fL[6]+2.449489742783178*alphaL[0]*fL[5]+1.414213562373095*(fL[4]*alphaL[5]+fL[3]*alphaL[4])+2.449489742783178*alphaL[1]*fL[2]+1.414213562373095*(alphaL[0]*fL[1]+fL[0]*alphaL[1])); 
  incr[2] = -0.125*(4.242640687119286*(alphaL[5]*fL[12]+alphaL[4]*fL[11]+alphaL[3]*fL[9])+2.449489742783178*alphaL[5]*fL[8]+4.242640687119286*alphaL[2]*fL[7]+2.449489742783178*alphaL[4]*fL[6]+4.242640687119286*alphaL[1]*fL[5]+2.449489742783178*(alphaL[3]*fL[4]+alphaL[2]*fL[3])+4.242640687119286*alphaL[0]*fL[2]+2.449489742783178*(alphaL[1]*fL[1]+alphaL[0]*fL[0])); 
  incr[3] = 0.125*(2.449489742783178*(alphaL[5]*fL[15]+alphaL[3]*fL[14])+1.414213562373095*alphaL[5]*fL[13]+2.449489742783178*alphaL[1]*fL[11]+1.414213562373095*alphaL[3]*fL[10]+2.449489742783178*alphaL[0]*fL[7]+1.414213562373095*alphaL[1]*fL[6]+2.449489742783178*alphaL[4]*fL[5]+1.414213562373095*(fL[1]*alphaL[4]+alphaL[0]*fL[3])+alphaL[2]*(2.449489742783178*fL[2]+1.414213562373095*fL[0])); 
  incr[4] = 0.125*(2.449489742783178*(alphaL[4]*fL[15]+alphaL[2]*fL[14])+1.414213562373095*alphaL[4]*fL[13]+2.449489742783178*alphaL[1]*fL[12]+1.414213562373095*alphaL[2]*fL[10]+2.449489742783178*alphaL[0]*fL[9]+1.414213562373095*alphaL[1]*fL[8]+2.449489742783178*alphaL[5]*fL[5]+1.414213562373095*(fL[1]*alphaL[5]+alphaL[0]*fL[4])+(2.449489742783178*fL[2]+1.414213562373095*fL[0])*alphaL[3]); 
  incr[5] = -0.125*(4.242640687119286*(alphaL[3]*fL[12]+alphaL[2]*fL[11]+alphaL[5]*fL[9])+2.449489742783178*alphaL[3]*fL[8]+4.242640687119286*alphaL[4]*fL[7]+2.449489742783178*alphaL[2]*fL[6]+4.242640687119286*alphaL[0]*fL[5]+2.449489742783178*(fL[4]*alphaL[5]+fL[3]*alphaL[4])+4.242640687119286*alphaL[1]*fL[2]+2.449489742783178*(alphaL[0]*fL[1]+fL[0]*alphaL[1])); 
  incr[6] = 0.125*(2.449489742783178*(alphaL[3]*fL[15]+alphaL[5]*fL[14])+1.414213562373095*alphaL[3]*fL[13]+2.449489742783178*alphaL[0]*fL[11]+1.414213562373095*alphaL[5]*fL[10]+2.449489742783178*alphaL[1]*fL[7]+1.414213562373095*alphaL[0]*fL[6]+2.449489742783178*alphaL[2]*fL[5]+(2.449489742783178*fL[2]+1.414213562373095*fL[0])*alphaL[4]+1.414213562373095*(alphaL[1]*fL[3]+fL[1]*alphaL[2])); 
  incr[7] = -0.125*(4.242640687119286*(alphaL[5]*fL[15]+alphaL[3]*fL[14])+2.449489742783178*alphaL[5]*fL[13]+4.242640687119286*alphaL[1]*fL[11]+2.449489742783178*alphaL[3]*fL[10]+4.242640687119286*alphaL[0]*fL[7]+2.449489742783178*alphaL[1]*fL[6]+4.242640687119286*alphaL[4]*fL[5]+2.449489742783178*(fL[1]*alphaL[4]+alphaL[0]*fL[3])+alphaL[2]*(4.242640687119286*fL[2]+2.449489742783178*fL[0])); 
  incr[8] = 0.125*(2.449489742783178*(alphaL[2]*fL[15]+alphaL[4]*fL[14])+1.414213562373095*alphaL[2]*fL[13]+2.449489742783178*alphaL[0]*fL[12]+1.414213562373095*alphaL[4]*fL[10]+2.449489742783178*alphaL[1]*fL[9]+1.414213562373095*alphaL[0]*fL[8]+2.449489742783178*alphaL[3]*fL[5]+(2.449489742783178*fL[2]+1.414213562373095*fL[0])*alphaL[5]+1.414213562373095*(alphaL[1]*fL[4]+fL[1]*alphaL[3])); 
  incr[9] = -0.125*(4.242640687119286*(alphaL[4]*fL[15]+alphaL[2]*fL[14])+2.449489742783178*alphaL[4]*fL[13]+4.242640687119286*alphaL[1]*fL[12]+2.449489742783178*alphaL[2]*fL[10]+4.242640687119286*alphaL[0]*fL[9]+2.449489742783178*alphaL[1]*fL[8]+4.242640687119286*alphaL[5]*fL[5]+2.449489742783178*(fL[1]*alphaL[5]+alphaL[0]*fL[4])+(4.242640687119286*fL[2]+2.449489742783178*fL[0])*alphaL[3]); 
  incr[10] = 0.125*(2.449489742783178*(alphaL[1]*fL[15]+alphaL[0]*fL[14])+1.414213562373095*alphaL[1]*fL[13]+2.449489742783178*(alphaL[4]*fL[12]+alphaL[5]*fL[11])+1.414213562373095*alphaL[0]*fL[10]+2.449489742783178*alphaL[2]*fL[9]+1.414213562373095*alphaL[4]*fL[8]+2.449489742783178*alphaL[3]*fL[7]+1.414213562373095*(alphaL[5]*fL[6]+alphaL[2]*fL[4]+alphaL[3]*fL[3])); 
  incr[11] = -0.125*(4.242640687119286*(alphaL[3]*fL[15]+alphaL[5]*fL[14])+2.449489742783178*alphaL[3]*fL[13]+4.242640687119286*alphaL[0]*fL[11]+2.449489742783178*alphaL[5]*fL[10]+4.242640687119286*alphaL[1]*fL[7]+2.449489742783178*alphaL[0]*fL[6]+4.242640687119286*alphaL[2]*fL[5]+(4.242640687119286*fL[2]+2.449489742783178*fL[0])*alphaL[4]+2.449489742783178*(alphaL[1]*fL[3]+fL[1]*alphaL[2])); 
  incr[12] = -0.125*(4.242640687119286*(alphaL[2]*fL[15]+alphaL[4]*fL[14])+2.449489742783178*alphaL[2]*fL[13]+4.242640687119286*alphaL[0]*fL[12]+2.449489742783178*alphaL[4]*fL[10]+4.242640687119286*alphaL[1]*fL[9]+2.449489742783178*alphaL[0]*fL[8]+4.242640687119286*alphaL[3]*fL[5]+(4.242640687119286*fL[2]+2.449489742783178*fL[0])*alphaL[5]+2.449489742783178*(alphaL[1]*fL[4]+fL[1]*alphaL[3])); 
  incr[13] = 0.125*(2.449489742783178*(alphaL[0]*fL[15]+alphaL[1]*fL[14])+1.414213562373095*alphaL[0]*fL[13]+2.449489742783178*(alphaL[2]*fL[12]+alphaL[3]*fL[11])+1.414213562373095*alphaL[1]*fL[10]+2.449489742783178*alphaL[4]*fL[9]+1.414213562373095*alphaL[2]*fL[8]+2.449489742783178*alphaL[5]*fL[7]+1.414213562373095*(alphaL[3]*fL[6]+fL[3]*alphaL[5]+alphaL[4]*fL[4])); 
  incr[14] = -0.125*(4.242640687119286*(alphaL[1]*fL[15]+alphaL[0]*fL[14])+2.449489742783178*alphaL[1]*fL[13]+4.242640687119286*(alphaL[4]*fL[12]+alphaL[5]*fL[11])+2.449489742783178*alphaL[0]*fL[10]+4.242640687119286*alphaL[2]*fL[9]+2.449489742783178*alphaL[4]*fL[8]+4.242640687119286*alphaL[3]*fL[7]+2.449489742783178*(alphaL[5]*fL[6]+alphaL[2]*fL[4]+alphaL[3]*fL[3])); 
  incr[15] = -0.125*(4.242640687119286*(alphaL[0]*fL[15]+alphaL[1]*fL[14])+2.449489742783178*alphaL[0]*fL[13]+4.242640687119286*(alphaL[2]*fL[12]+alphaL[3]*fL[11])+2.449489742783178*alphaL[1]*fL[10]+4.242640687119286*alphaL[4]*fL[9]+2.449489742783178*alphaL[2]*fL[8]+4.242640687119286*alphaL[5]*fL[7]+2.449489742783178*(alphaL[3]*fL[6]+fL[3]*alphaL[5]+alphaL[4]*fL[4])); 
  } else { 
  incr[0] = -0.125*(2.449489742783178*(alphaL[5]*fR[12]+alphaL[4]*fR[11]+alphaL[3]*fR[9])-1.414213562373095*alphaL[5]*fR[8]+2.449489742783178*alphaL[2]*fR[7]-1.414213562373095*alphaL[4]*fR[6]+2.449489742783178*alphaL[1]*fR[5]-1.414213562373095*(alphaL[3]*fR[4]+alphaL[2]*fR[3])+2.449489742783178*alphaL[0]*fR[2]-1.414213562373095*(alphaL[1]*fR[1]+alphaL[0]*fR[0])); 
  incr[1] = -0.125*(2.449489742783178*(alphaL[3]*fR[12]+alphaL[2]*fR[11]+alphaL[5]*fR[9])-1.414213562373095*alphaL[3]*fR[8]+2.449489742783178*alphaL[4]*fR[7]-1.414213562373095*alphaL[2]*fR[6]+2.449489742783178*alphaL[0]*fR[5]-1.414213562373095*(fR[4]*alphaL[5]+fR[3]*alphaL[4])+2.449489742783178*alphaL[1]*fR[2]-1.414213562373095*(alphaL[0]*fR[1]+fR[0]*alphaL[1])); 
  incr[2] = 0.125*(4.242640687119286*(alphaL[5]*fR[12]+alphaL[4]*fR[11]+alphaL[3]*fR[9])-2.449489742783178*alphaL[5]*fR[8]+4.242640687119286*alphaL[2]*fR[7]-2.449489742783178*alphaL[4]*fR[6]+4.242640687119286*alphaL[1]*fR[5]-2.449489742783178*(alphaL[3]*fR[4]+alphaL[2]*fR[3])+4.242640687119286*alphaL[0]*fR[2]-2.449489742783178*(alphaL[1]*fR[1]+alphaL[0]*fR[0])); 
  incr[3] = -0.125*(2.449489742783178*(alphaL[5]*fR[15]+alphaL[3]*fR[14])-1.414213562373095*alphaL[5]*fR[13]+2.449489742783178*alphaL[1]*fR[11]-1.414213562373095*alphaL[3]*fR[10]+2.449489742783178*alphaL[0]*fR[7]-1.414213562373095*alphaL[1]*fR[6]+2.449489742783178*alphaL[4]*fR[5]-1.414213562373095*(fR[1]*alphaL[4]+alphaL[0]*fR[3])+alphaL[2]*(2.449489742783178*fR[2]-1.414213562373095*fR[0])); 
  incr[4] = -0.125*(2.449489742783178*(alphaL[4]*fR[15]+alphaL[2]*fR[14])-1.414213562373095*alphaL[4]*fR[13]+2.449489742783178*alphaL[1]*fR[12]-1.414213562373095*alphaL[2]*fR[10]+2.449489742783178*alphaL[0]*fR[9]-1.414213562373095*alphaL[1]*fR[8]+2.449489742783178*alphaL[5]*fR[5]-1.414213562373095*(fR[1]*alphaL[5]+alphaL[0]*fR[4])+(2.449489742783178*fR[2]-1.414213562373095*fR[0])*alphaL[3]); 
  incr[5] = 0.125*(4.242640687119286*(alphaL[3]*fR[12]+alphaL[2]*fR[11]+alphaL[5]*fR[9])-2.449489742783178*alphaL[3]*fR[8]+4.242640687119286*alphaL[4]*fR[7]-2.449489742783178*alphaL[2]*fR[6]+4.242640687119286*alphaL[0]*fR[5]-2.449489742783178*(fR[4]*alphaL[5]+fR[3]*alphaL[4])+4.242640687119286*alphaL[1]*fR[2]-2.449489742783178*(alphaL[0]*fR[1]+fR[0]*alphaL[1])); 
  incr[6] = -0.125*(2.449489742783178*(alphaL[3]*fR[15]+alphaL[5]*fR[14])-1.414213562373095*alphaL[3]*fR[13]+2.449489742783178*alphaL[0]*fR[11]-1.414213562373095*alphaL[5]*fR[10]+2.449489742783178*alphaL[1]*fR[7]-1.414213562373095*alphaL[0]*fR[6]+2.449489742783178*alphaL[2]*fR[5]+(2.449489742783178*fR[2]-1.414213562373095*fR[0])*alphaL[4]-1.414213562373095*(alphaL[1]*fR[3]+fR[1]*alphaL[2])); 
  incr[7] = 0.125*(4.242640687119286*(alphaL[5]*fR[15]+alphaL[3]*fR[14])-2.449489742783178*alphaL[5]*fR[13]+4.242640687119286*alphaL[1]*fR[11]-2.449489742783178*alphaL[3]*fR[10]+4.242640687119286*alphaL[0]*fR[7]-2.449489742783178*alphaL[1]*fR[6]+4.242640687119286*alphaL[4]*fR[5]-2.449489742783178*(fR[1]*alphaL[4]+alphaL[0]*fR[3])+alphaL[2]*(4.242640687119286*fR[2]-2.449489742783178*fR[0])); 
  incr[8] = -0.125*(2.449489742783178*(alphaL[2]*fR[15]+alphaL[4]*fR[14])-1.414213562373095*alphaL[2]*fR[13]+2.449489742783178*alphaL[0]*fR[12]-1.414213562373095*alphaL[4]*fR[10]+2.449489742783178*alphaL[1]*fR[9]-1.414213562373095*alphaL[0]*fR[8]+2.449489742783178*alphaL[3]*fR[5]+(2.449489742783178*fR[2]-1.414213562373095*fR[0])*alphaL[5]-1.414213562373095*(alphaL[1]*fR[4]+fR[1]*alphaL[3])); 
  incr[9] = 0.125*(4.242640687119286*(alphaL[4]*fR[15]+alphaL[2]*fR[14])-2.449489742783178*alphaL[4]*fR[13]+4.242640687119286*alphaL[1]*fR[12]-2.449489742783178*alphaL[2]*fR[10]+4.242640687119286*alphaL[0]*fR[9]-2.449489742783178*alphaL[1]*fR[8]+4.242640687119286*alphaL[5]*fR[5]-2.449489742783178*(fR[1]*alphaL[5]+alphaL[0]*fR[4])+(4.242640687119286*fR[2]-2.449489742783178*fR[0])*alphaL[3]); 
  incr[10] = -0.125*(2.449489742783178*(alphaL[1]*fR[15]+alphaL[0]*fR[14])-1.414213562373095*alphaL[1]*fR[13]+2.449489742783178*(alphaL[4]*fR[12]+alphaL[5]*fR[11])-1.414213562373095*alphaL[0]*fR[10]+2.449489742783178*alphaL[2]*fR[9]-1.414213562373095*alphaL[4]*fR[8]+2.449489742783178*alphaL[3]*fR[7]-1.414213562373095*(alphaL[5]*fR[6]+alphaL[2]*fR[4]+alphaL[3]*fR[3])); 
  incr[11] = 0.125*(4.242640687119286*(alphaL[3]*fR[15]+alphaL[5]*fR[14])-2.449489742783178*alphaL[3]*fR[13]+4.242640687119286*alphaL[0]*fR[11]-2.449489742783178*alphaL[5]*fR[10]+4.242640687119286*alphaL[1]*fR[7]-2.449489742783178*alphaL[0]*fR[6]+4.242640687119286*alphaL[2]*fR[5]+(4.242640687119286*fR[2]-2.449489742783178*fR[0])*alphaL[4]-2.449489742783178*(alphaL[1]*fR[3]+fR[1]*alphaL[2])); 
  incr[12] = 0.125*(4.242640687119286*(alphaL[2]*fR[15]+alphaL[4]*fR[14])-2.449489742783178*alphaL[2]*fR[13]+4.242640687119286*alphaL[0]*fR[12]-2.449489742783178*alphaL[4]*fR[10]+4.242640687119286*alphaL[1]*fR[9]-2.449489742783178*alphaL[0]*fR[8]+4.242640687119286*alphaL[3]*fR[5]+(4.242640687119286*fR[2]-2.449489742783178*fR[0])*alphaL[5]-2.449489742783178*(alphaL[1]*fR[4]+fR[1]*alphaL[3])); 
  incr[13] = -0.125*(2.449489742783178*(alphaL[0]*fR[15]+alphaL[1]*fR[14])-1.414213562373095*alphaL[0]*fR[13]+2.449489742783178*(alphaL[2]*fR[12]+alphaL[3]*fR[11])-1.414213562373095*alphaL[1]*fR[10]+2.449489742783178*alphaL[4]*fR[9]-1.414213562373095*alphaL[2]*fR[8]+2.449489742783178*alphaL[5]*fR[7]-1.414213562373095*(alphaL[3]*fR[6]+fR[3]*alphaL[5]+alphaL[4]*fR[4])); 
  incr[14] = 0.125*(4.242640687119286*(alphaL[1]*fR[15]+alphaL[0]*fR[14])-2.449489742783178*alphaL[1]*fR[13]+4.242640687119286*(alphaL[4]*fR[12]+alphaL[5]*fR[11])-2.449489742783178*alphaL[0]*fR[10]+4.242640687119286*alphaL[2]*fR[9]-2.449489742783178*alphaL[4]*fR[8]+4.242640687119286*alphaL[3]*fR[7]-2.449489742783178*(alphaL[5]*fR[6]+alphaL[2]*fR[4]+alphaL[3]*fR[3])); 
  incr[15] = 0.125*(4.242640687119286*(alphaL[0]*fR[15]+alphaL[1]*fR[14])-2.449489742783178*alphaL[0]*fR[13]+4.242640687119286*(alphaL[2]*fR[12]+alphaL[3]*fR[11])-2.449489742783178*alphaL[1]*fR[10]+4.242640687119286*alphaL[4]*fR[9]-2.449489742783178*alphaL[2]*fR[8]+4.242640687119286*alphaL[5]*fR[7]-2.449489742783178*(alphaL[3]*fR[6]+fR[3]*alphaL[5]+alphaL[4]*fR[4])); 
  }
#elif upwindType == QUAD 
  double alphaOrdL;
  double fUpOrd[8];
  alphaOrdL = 0.3535533905932737*(alphaL[5]+alphaL[4])-0.3535533905932737*(alphaL[3]+alphaL[2]+alphaL[1])+0.3535533905932737*alphaL[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])+0.25*fR[13]-0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*fL[14]-0.25*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fR[11])+0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*(alphaL[1]+alphaL[0])-0.3535533905932737*(alphaL[5]+alphaL[4]+alphaL[3]+alphaL[2]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])+0.25*(fR[4]+fR[3])-0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdL)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14])+0.25*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[5]-0.3535533905932737*(alphaL[4]+alphaL[3])+0.3535533905932737*alphaL[2]-0.3535533905932737*alphaL[1]+0.3535533905932737*alphaL[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*fL[14]+0.25*(fR[13]+fL[13])-0.4330127018922193*fR[12]+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdL = (-0.3535533905932737*alphaL[5])+0.3535533905932737*alphaL[4]-0.3535533905932737*alphaL[3]+0.3535533905932737*(alphaL[2]+alphaL[1]+alphaL[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.25*fR[13]-0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]-0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])+0.25*fR[4]-0.25*(fL[4]+fR[3])+0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdL)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14])-0.25*(fR[13]+fL[13])+0.4330127018922193*fR[12]-0.4330127018922193*(fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])+0.4330127018922193*fR[9]-0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]-0.25*(fR[4]+fL[4])+0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdL = (-0.3535533905932737*alphaL[5])+0.3535533905932737*(alphaL[4]+alphaL[3])-0.3535533905932737*(alphaL[2]+alphaL[1])+0.3535533905932737*alphaL[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdL)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*fL[14]+0.25*(fR[13]+fL[13])+0.4330127018922193*fR[12]-0.4330127018922193*(fL[12]+fR[11])+0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*alphaL[5]-0.3535533905932737*alphaL[4]+0.3535533905932737*alphaL[3]-0.3535533905932737*alphaL[2]+0.3535533905932737*(alphaL[1]+alphaL[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.25*fR[13]-0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fL[11])+0.25*fR[10]-0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])-0.25*fR[4]+0.25*(fL[4]+fR[3])-0.25*fL[3]+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdL)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14])-0.25*(fR[13]+fL[13])-0.4330127018922193*fR[12]+0.4330127018922193*(fL[12]+fR[11])-0.4330127018922193*fL[11]-0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4])-0.25*(fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdL = (-0.3535533905932737*(alphaL[5]+alphaL[4]))+0.3535533905932737*(alphaL[3]+alphaL[2])-0.3535533905932737*alphaL[1]+0.3535533905932737*alphaL[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])+0.25*fR[13]-0.25*fL[13]-0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])+0.25*fR[8]-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])+0.25*fR[6]-0.25*fL[6]-0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])+0.25*fR[1]-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdL)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*fL[14]-0.25*(fR[13]+fL[13])+0.4330127018922193*(fR[12]+fR[11])-0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]-0.25*(fR[6]+fL[6])+0.4330127018922193*fR[5]-0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdL = 0.3535533905932737*(alphaL[5]+alphaL[4]+alphaL[3]+alphaL[2]+alphaL[1]+alphaL[0]); 
  fUpOrd[7] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.25*fR[13]+0.25*fL[13]+0.4330127018922193*(fR[12]+fL[12]+fR[11]+fL[11])-0.25*fR[10]+0.25*fL[10]+0.4330127018922193*(fR[9]+fL[9])-0.25*fR[8]+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.25*fR[6]+0.25*fL[6]+0.4330127018922193*(fR[5]+fL[5])-0.25*(fR[4]+fR[3])+0.25*(fL[4]+fL[3])+0.4330127018922193*(fR[2]+fL[2])-0.25*(fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdL)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14])+0.25*(fR[13]+fL[13])-0.4330127018922193*(fR[12]+fR[11])+0.4330127018922193*(fL[12]+fL[11])+0.25*(fR[10]+fL[10])-0.4330127018922193*fR[9]+0.4330127018922193*fL[9]+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*fL[7]+0.25*(fR[6]+fL[6])-0.4330127018922193*fR[5]+0.4330127018922193*fL[5]+0.25*(fR[4]+fL[4]+fR[3]+fL[3])-0.4330127018922193*fR[2]+0.4330127018922193*fL[2]+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 

  double fUp[8];
  fUp[0] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0]); 
  fUp[1] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*fUpOrd[4]+fUpOrd[3]-1.0*fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 
  fUp[2] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4])+fUpOrd[3]+fUpOrd[2]-1.0*(fUpOrd[1]+fUpOrd[0])); 
  fUp[3] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]+fUpOrd[5]+fUpOrd[4]-1.0*(fUpOrd[3]+fUpOrd[2]+fUpOrd[1]+fUpOrd[0])); 
  fUp[4] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]+fUpOrd[3]-1.0*(fUpOrd[2]+fUpOrd[1])+fUpOrd[0]); 
  fUp[5] = 0.3535533905932737*(fUpOrd[7]-1.0*fUpOrd[6]+fUpOrd[5]-1.0*(fUpOrd[4]+fUpOrd[3])+fUpOrd[2]-1.0*fUpOrd[1]+fUpOrd[0]); 
  fUp[6] = 0.3535533905932737*(fUpOrd[7]+fUpOrd[6]-1.0*(fUpOrd[5]+fUpOrd[4]+fUpOrd[3]+fUpOrd[2])+fUpOrd[1]+fUpOrd[0]); 
  fUp[7] = 0.3535533905932737*(fUpOrd[7]-1.0*(fUpOrd[6]+fUpOrd[5])+fUpOrd[4]-1.0*fUpOrd[3]+fUpOrd[2]+fUpOrd[1]-1.0*fUpOrd[0]); 

  incr[0] = 0.25*(alphaL[5]*fUp[5]+alphaL[4]*fUp[4]+alphaL[3]*fUp[3]+alphaL[2]*fUp[2]+alphaL[1]*fUp[1]+alphaL[0]*fUp[0]); 
  incr[1] = 0.25*(alphaL[3]*fUp[5]+fUp[3]*alphaL[5]+alphaL[2]*fUp[4]+fUp[2]*alphaL[4]+alphaL[0]*fUp[1]+fUp[0]*alphaL[1]); 
  incr[2] = -0.4330127018922193*(alphaL[5]*fUp[5]+alphaL[4]*fUp[4]+alphaL[3]*fUp[3]+alphaL[2]*fUp[2]+alphaL[1]*fUp[1]+alphaL[0]*fUp[0]); 
  incr[3] = 0.25*(alphaL[5]*fUp[7]+alphaL[3]*fUp[6]+alphaL[1]*fUp[4]+fUp[1]*alphaL[4]+alphaL[0]*fUp[2]+fUp[0]*alphaL[2]); 
  incr[4] = 0.25*(alphaL[4]*fUp[7]+alphaL[2]*fUp[6]+alphaL[1]*fUp[5]+fUp[1]*alphaL[5]+alphaL[0]*fUp[3]+fUp[0]*alphaL[3]); 
  incr[5] = -0.4330127018922193*(alphaL[3]*fUp[5]+fUp[3]*alphaL[5]+alphaL[2]*fUp[4]+fUp[2]*alphaL[4]+alphaL[0]*fUp[1]+fUp[0]*alphaL[1]); 
  incr[6] = 0.25*(alphaL[3]*fUp[7]+alphaL[5]*fUp[6]+alphaL[0]*fUp[4]+fUp[0]*alphaL[4]+alphaL[1]*fUp[2]+fUp[1]*alphaL[2]); 
  incr[7] = -0.4330127018922193*(alphaL[5]*fUp[7]+alphaL[3]*fUp[6]+alphaL[1]*fUp[4]+fUp[1]*alphaL[4]+alphaL[0]*fUp[2]+fUp[0]*alphaL[2]); 
  incr[8] = 0.25*(alphaL[2]*fUp[7]+alphaL[4]*fUp[6]+alphaL[0]*fUp[5]+fUp[0]*alphaL[5]+alphaL[1]*fUp[3]+fUp[1]*alphaL[3]); 
  incr[9] = -0.4330127018922193*(alphaL[4]*fUp[7]+alphaL[2]*fUp[6]+alphaL[1]*fUp[5]+fUp[1]*alphaL[5]+alphaL[0]*fUp[3]+fUp[0]*alphaL[3]); 
  incr[10] = 0.25*(alphaL[1]*fUp[7]+alphaL[0]*fUp[6]+alphaL[4]*fUp[5]+fUp[4]*alphaL[5]+alphaL[2]*fUp[3]+fUp[2]*alphaL[3]); 
  incr[11] = -0.4330127018922193*(alphaL[3]*fUp[7]+alphaL[5]*fUp[6]+alphaL[0]*fUp[4]+fUp[0]*alphaL[4]+alphaL[1]*fUp[2]+fUp[1]*alphaL[2]); 
  incr[12] = -0.4330127018922193*(alphaL[2]*fUp[7]+alphaL[4]*fUp[6]+alphaL[0]*fUp[5]+fUp[0]*alphaL[5]+alphaL[1]*fUp[3]+fUp[1]*alphaL[3]); 
  incr[13] = 0.25*(alphaL[0]*fUp[7]+alphaL[1]*fUp[6]+alphaL[2]*fUp[5]+fUp[2]*alphaL[5]+alphaL[3]*fUp[4]+fUp[3]*alphaL[4]); 
  incr[14] = -0.4330127018922193*(alphaL[1]*fUp[7]+alphaL[0]*fUp[6]+alphaL[4]*fUp[5]+fUp[4]*alphaL[5]+alphaL[2]*fUp[3]+fUp[2]*alphaL[3]); 
  incr[15] = -0.4330127018922193*(alphaL[0]*fUp[7]+alphaL[1]*fUp[6]+alphaL[2]*fUp[5]+fUp[2]*alphaL[5]+alphaL[3]*fUp[4]+fUp[3]*alphaL[4]); 

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

  return std::abs(alphaSurfAvgL);
} 
double EmGyrokineticGenGeoSurf2x2vSer_yR_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
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
  BstarYdBmagR[0] = -(0.8660254037844386*rdx2R*(2.0*jacobTotInv[0]*b_z[1]*m_*wvparR+(2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*q_))/q_; 
  BstarYdBmagR[1] = -(0.8660254037844386*rdx2R*(2.0*b_z[1]*jacobTotInv[1]*m_*wvparR+((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])*q_))/q_; 
  BstarYdBmagR[2] = -0.8660254037844386*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])*rdx2R; 
  BstarYdBmagR[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[5] = -0.8660254037844386*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])*rdx2R; 
  BstarYdBmagR[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = -(0.1767766952966368*((b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*(3.0*hamilR[5]-1.732050807568877*hamilR[1])*m_*rdx2R+(3.0*BstarYdBmagR[2]-1.732050807568877*BstarYdBmagR[0])*hamilR[3]*q_*rdvpar2R))/(m_*q_); 
  alphaR[1] = -(0.1767766952966368*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*(3.0*hamilR[5]-1.732050807568877*hamilR[1])*m_*rdx2R+hamilR[3]*(3.0*BstarYdBmagR[5]-1.732050807568877*BstarYdBmagR[1])*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = (0.3061862178478971*BstarYdBmagR[3]*hamilR[3]*rdvpar2R)/m_; 
  alphaR[3] = (0.3061862178478971*(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamilR[8]*rdx2R)/q_; 
  alphaR[4] = (0.3061862178478971*hamilR[3]*BstarYdBmagR[6]*rdvpar2R)/m_; 
  alphaR[5] = (0.3061862178478971*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*hamilR[8]*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*(((3.0*b_z[1]*jacobTotInv[1]+3.0*b_z[0]*jacobTotInv[0])*hamilR[5]-1.732050807568877*b_z[1]*hamilR[1]*jacobTotInv[1]-1.732050807568877*b_z[0]*jacobTotInv[0]*hamilR[1])*m_*rdx2R+(3.0*BstarYdBmagR[2]-1.732050807568877*BstarYdBmagR[0])*hamilR[3]*q_*rdvpar2R))/(m_*q_); 

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
double EmGyrokineticGenGeoSurf2x2vSer_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
{ 
  // Apar: parallel component of magnetic vector potential.
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
  BstarXdBmagR[0] = 0.8660254037844386*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])*rdy2R; 
  BstarXdBmagR[1] = 0.1732050807568877*((9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])*rdy2R; 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = -(0.8660254037844386*rdx2R*(2.0*jacobTotInv[0]*b_z[1]*m_*wvparR+(2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*q_))/q_; 
  BstarYdBmagR[1] = -(0.8660254037844386*rdx2R*(2.0*b_z[1]*jacobTotInv[1]*m_*wvparR+((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])*q_))/q_; 
  BstarYdBmagR[2] = -0.8660254037844386*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])*rdx2R; 
  BstarYdBmagR[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[5] = -0.8660254037844386*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])*rdx2R; 
  BstarYdBmagR[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = (0.1767766952966368*((hamilR[5]*(3.0*BstarYdBmagR[6]-1.732050807568877*BstarYdBmagR[1])+hamilR[2]*(3.0*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0]))*rdy2R-1.732050807568877*BstarXdBmagR[0]*hamilR[1]*rdx2R))/m_; 
  alphaR[1] = (0.1767766952966368*((3.0*hamilR[2]*BstarYdBmagR[6]+(3.0*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0])*hamilR[5]-1.732050807568877*BstarYdBmagR[1]*hamilR[2])*rdy2R-1.732050807568877*BstarXdBmagR[1]*hamilR[1]*rdx2R))/m_; 
  alphaR[2] = -(0.3061862178478971*((BstarYdBmagR[5]*hamilR[5]+BstarYdBmagR[2]*hamilR[2])*rdy2R+BstarXdBmagR[0]*hamilR[5]*rdx2R))/m_; 
  alphaR[3] = -(0.3061862178478971*BstarXdBmagR[0]*hamilR[8]*rdx2R)/m_; 
  alphaR[4] = -(0.3061862178478971*((BstarYdBmagR[2]*hamilR[5]+hamilR[2]*BstarYdBmagR[5])*rdy2R+BstarXdBmagR[1]*hamilR[5]*rdx2R))/m_; 
  alphaR[5] = -(0.3061862178478971*BstarXdBmagR[1]*hamilR[8]*rdx2R)/m_; 

  double alphaUpR[8]; 
  alphaUpR[0] = (0.1767766952966368*((hamilR[5]*(3.0*BstarYdBmagR[6]-1.732050807568877*BstarYdBmagR[1])+hamilR[2]*(3.0*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0]))*rdy2R-1.732050807568877*BstarXdBmagR[0]*hamilR[1]*rdx2R-8.0*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = (0.1767766952966368*((3.0*hamilR[2]*BstarYdBmagR[6]+(3.0*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0])*hamilR[5]-1.732050807568877*BstarYdBmagR[1]*hamilR[2])*rdy2R-1.732050807568877*BstarXdBmagR[1]*hamilR[1]*rdx2R-8.0*dApardtPrev[1]*q_))/m_; 
  alphaUpR[2] = -(0.1020620726159657*(3.0*((BstarYdBmagR[5]*hamilR[5]+BstarYdBmagR[2]*hamilR[2])*rdy2R+BstarXdBmagR[0]*hamilR[5]*rdx2R)+13.85640646055102*dApardtPrev[2]*q_))/m_; 
  alphaUpR[3] = -(0.3061862178478971*BstarXdBmagR[0]*hamilR[8]*rdx2R)/m_; 
  alphaUpR[4] = -(0.1767766952966368*(1.732050807568877*((BstarYdBmagR[2]*hamilR[5]+hamilR[2]*BstarYdBmagR[5])*rdy2R+BstarXdBmagR[1]*hamilR[5]*rdx2R)+8.0*dApardtPrev[3]*q_))/m_; 
  alphaUpR[5] = -(0.3061862178478971*BstarXdBmagR[1]*hamilR[8]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*((3.0*hamilR[5]*BstarYdBmagR[6]-1.732050807568877*BstarYdBmagR[1]*hamilR[5]+3.0*hamilR[2]*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0]*hamilR[2])*rdy2R-1.732050807568877*BstarXdBmagR[0]*hamilR[1]*rdx2R-8.0*dApardtPrev[0]*q_))/m_; 

  double incr[16]; 
  double incrEmMod[16]; 
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
  incrEmMod[0] = 0.5*(1.732050807568877*fL[3]+fL[0]); 
  incrEmMod[1] = 0.5*(1.732050807568877*fL[6]+fL[1]); 
  incrEmMod[2] = 0.5*(1.732050807568877*fL[7]+fL[2]); 
  incrEmMod[3] = -0.5*(3.0*fL[3]+1.732050807568877*fL[0]); 
  incrEmMod[4] = 0.5*(1.732050807568877*fL[10]+fL[4]); 
  incrEmMod[5] = 0.5*(1.732050807568877*fL[11]+fL[5]); 
  incrEmMod[6] = -0.5*(3.0*fL[6]+1.732050807568877*fL[1]); 
  incrEmMod[7] = -0.5*(3.0*fL[7]+1.732050807568877*fL[2]); 
  incrEmMod[8] = 0.5*(1.732050807568877*fL[13]+fL[8]); 
  incrEmMod[9] = 0.5*(1.732050807568877*fL[14]+fL[9]); 
  incrEmMod[10] = -0.5*(3.0*fL[10]+1.732050807568877*fL[4]); 
  incrEmMod[11] = -0.5*(3.0*fL[11]+1.732050807568877*fL[5]); 
  incrEmMod[12] = 0.5*(1.732050807568877*fL[15]+fL[12]); 
  incrEmMod[13] = -0.5*(3.0*fL[13]+1.732050807568877*fL[8]); 
  incrEmMod[14] = -0.5*(3.0*fL[14]+1.732050807568877*fL[9]); 
  incrEmMod[15] = -0.5*(3.0*fL[15]+1.732050807568877*fL[12]); 
  } else { 
  incrEmMod[0] = -0.5*(1.732050807568877*fR[3]-1.0*fR[0]); 
  incrEmMod[1] = -0.5*(1.732050807568877*fR[6]-1.0*fR[1]); 
  incrEmMod[2] = -0.5*(1.732050807568877*fR[7]-1.0*fR[2]); 
  incrEmMod[3] = 0.5*(3.0*fR[3]-1.732050807568877*fR[0]); 
  incrEmMod[4] = -0.5*(1.732050807568877*fR[10]-1.0*fR[4]); 
  incrEmMod[5] = -0.5*(1.732050807568877*fR[11]-1.0*fR[5]); 
  incrEmMod[6] = 0.5*(3.0*fR[6]-1.732050807568877*fR[1]); 
  incrEmMod[7] = 0.5*(3.0*fR[7]-1.732050807568877*fR[2]); 
  incrEmMod[8] = -0.5*(1.732050807568877*fR[13]-1.0*fR[8]); 
  incrEmMod[9] = -0.5*(1.732050807568877*fR[14]-1.0*fR[9]); 
  incrEmMod[10] = 0.5*(3.0*fR[10]-1.732050807568877*fR[4]); 
  incrEmMod[11] = 0.5*(3.0*fR[11]-1.732050807568877*fR[5]); 
  incrEmMod[12] = -0.5*(1.732050807568877*fR[15]-1.0*fR[12]); 
  incrEmMod[13] = 0.5*(3.0*fR[13]-1.732050807568877*fR[8]); 
  incrEmMod[14] = 0.5*(3.0*fR[14]-1.732050807568877*fR[9]); 
  incrEmMod[15] = 0.5*(3.0*fR[15]-1.732050807568877*fR[12]); 
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
  alphaOrdR = 0.3535533905932737*(alphaUpR[5]+alphaUpR[4])-0.3535533905932737*(alphaUpR[3]+alphaUpR[2]+alphaUpR[1])+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])-0.25*(fR[9]+fR[8])+0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14]+fR[13])+0.4330127018922193*(fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[1]+alphaUpR[0])-0.3535533905932737*(alphaUpR[5]+alphaUpR[4]+alphaUpR[3]+alphaUpR[2]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14]+fR[13])-0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaUpR[5]-0.3535533905932737*(alphaUpR[4]+alphaUpR[3])+0.3535533905932737*alphaUpR[2]-0.3535533905932737*alphaUpR[1]+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*(fL[14]+fR[13])+0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaUpR[5])+0.3535533905932737*alphaUpR[4]-0.3535533905932737*alphaUpR[3]+0.3535533905932737*(alphaUpR[2]+alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]+fR[13]+fL[13]))+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])+0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*(fR[2]+fR[1]+fR[0])+0.25*(fL[2]+fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14]+fR[13])-0.4330127018922193*(fL[15]+fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9]+fR[8]+fL[8])-0.4330127018922193*(fR[7]+fR[6])+0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2]+fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaUpR[5])+0.3535533905932737*(alphaUpR[4]+alphaUpR[3])-0.3535533905932737*(alphaUpR[2]+alphaUpR[1])+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*(fR[5]+fR[4])+0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14]+fR[13])-0.4330127018922193*(fL[14]+fL[13])+0.25*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fR[10])+0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaUpR[5]-0.3535533905932737*alphaUpR[4]+0.3535533905932737*alphaUpR[3]-0.3535533905932737*alphaUpR[2]+0.3535533905932737*(alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14]+fR[13])+0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*(alphaUpR[5]+alphaUpR[4]))+0.3535533905932737*(alphaUpR[3]+alphaUpR[2])-0.3535533905932737*alphaUpR[1]+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*(fL[14]+fR[13])-0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[5]+alphaUpR[4]+alphaUpR[3]+alphaUpR[2]+alphaUpR[1]+alphaUpR[0]); 
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

  incrEmMod[0] = 0.7071067811865475*fUp[0]; 
  incrEmMod[1] = 0.7071067811865475*fUp[1]; 
  incrEmMod[2] = 0.7071067811865475*fUp[2]; 
  incrEmMod[3] = -1.224744871391589*fUp[0]; 
  incrEmMod[4] = 0.7071067811865475*fUp[3]; 
  incrEmMod[5] = 0.7071067811865475*fUp[4]; 
  incrEmMod[6] = -1.224744871391589*fUp[1]; 
  incrEmMod[7] = -1.224744871391589*fUp[2]; 
  incrEmMod[8] = 0.7071067811865475*fUp[5]; 
  incrEmMod[9] = 0.7071067811865475*fUp[6]; 
  incrEmMod[10] = -1.224744871391589*fUp[3]; 
  incrEmMod[11] = -1.224744871391589*fUp[4]; 
  incrEmMod[12] = 0.7071067811865475*fUp[7]; 
  incrEmMod[13] = -1.224744871391589*fUp[5]; 
  incrEmMod[14] = -1.224744871391589*fUp[6]; 
  incrEmMod[15] = -1.224744871391589*fUp[7]; 

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

  emModR[0] += incrEmMod[0]*rdvpar2R; 
  emModR[1] += incrEmMod[1]*rdvpar2R; 
  emModR[2] += incrEmMod[2]*rdvpar2R; 
  emModR[3] += incrEmMod[3]*rdvpar2R; 
  emModR[4] += incrEmMod[4]*rdvpar2R; 
  emModR[5] += incrEmMod[5]*rdvpar2R; 
  emModR[6] += incrEmMod[6]*rdvpar2R; 
  emModR[7] += incrEmMod[7]*rdvpar2R; 
  emModR[8] += incrEmMod[8]*rdvpar2R; 
  emModR[9] += incrEmMod[9]*rdvpar2R; 
  emModR[10] += incrEmMod[10]*rdvpar2R; 
  emModR[11] += incrEmMod[11]*rdvpar2R; 
  emModR[12] += incrEmMod[12]*rdvpar2R; 
  emModR[13] += incrEmMod[13]*rdvpar2R; 
  emModR[14] += incrEmMod[14]*rdvpar2R; 
  emModR[15] += incrEmMod[15]*rdvpar2R; 
  emModL[0] += -1.0*incrEmMod[0]*rdvpar2L; 
  emModL[1] += -1.0*incrEmMod[1]*rdvpar2L; 
  emModL[2] += -1.0*incrEmMod[2]*rdvpar2L; 
  emModL[3] += incrEmMod[3]*rdvpar2L; 
  emModL[4] += -1.0*incrEmMod[4]*rdvpar2L; 
  emModL[5] += -1.0*incrEmMod[5]*rdvpar2L; 
  emModL[6] += incrEmMod[6]*rdvpar2L; 
  emModL[7] += incrEmMod[7]*rdvpar2L; 
  emModL[8] += -1.0*incrEmMod[8]*rdvpar2L; 
  emModL[9] += -1.0*incrEmMod[9]*rdvpar2L; 
  emModL[10] += incrEmMod[10]*rdvpar2L; 
  emModL[11] += incrEmMod[11]*rdvpar2L; 
  emModL[12] += -1.0*incrEmMod[12]*rdvpar2L; 
  emModL[13] += incrEmMod[13]*rdvpar2L; 
  emModL[14] += incrEmMod[14]*rdvpar2L; 
  emModL[15] += incrEmMod[15]*rdvpar2L; 

  return std::abs(alphaSurfAvgR);
} 
double EmGyrokineticGenGeoSurf2x2vSerStep2_vpar_P1_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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
  BstarXdBmagR[0] = 0.8660254037844386*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])*rdy2R; 
  BstarXdBmagR[1] = 0.1732050807568877*((9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])*rdy2R; 

  double BstarYdBmagR[16]; 
  BstarYdBmagR[0] = -(0.8660254037844386*rdx2R*(2.0*jacobTotInv[0]*b_z[1]*m_*wvparR+(2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*q_))/q_; 
  BstarYdBmagR[1] = -(0.8660254037844386*rdx2R*(2.0*b_z[1]*jacobTotInv[1]*m_*wvparR+((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])*q_))/q_; 
  BstarYdBmagR[2] = -0.8660254037844386*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])*rdx2R; 
  BstarYdBmagR[3] = -(1.0*jacobTotInv[0]*b_z[1]*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[5] = -0.8660254037844386*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])*rdx2R; 
  BstarYdBmagR[6] = -(1.0*b_z[1]*jacobTotInv[1]*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[8]; 
  alphaR[0] = -(1.414213562373095*dApardt[0]*q_)/m_; 
  alphaR[1] = -(1.414213562373095*dApardt[1]*q_)/m_; 
  alphaR[2] = -(1.414213562373095*dApardt[2]*q_)/m_; 
  alphaR[4] = -(1.414213562373095*dApardt[3]*q_)/m_; 

  double alphaUpR[8]; 
  alphaUpR[0] = (0.1767766952966368*((hamilR[5]*(3.0*BstarYdBmagR[6]-1.732050807568877*BstarYdBmagR[1])+hamilR[2]*(3.0*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0]))*rdy2R-1.732050807568877*BstarXdBmagR[0]*hamilR[1]*rdx2R-8.0*dApardtPrev[0]*q_))/m_; 
  alphaUpR[1] = (0.1767766952966368*((3.0*hamilR[2]*BstarYdBmagR[6]+(3.0*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0])*hamilR[5]-1.732050807568877*BstarYdBmagR[1]*hamilR[2])*rdy2R-1.732050807568877*BstarXdBmagR[1]*hamilR[1]*rdx2R-8.0*dApardtPrev[1]*q_))/m_; 
  alphaUpR[2] = -(0.1020620726159657*(3.0*((BstarYdBmagR[5]*hamilR[5]+BstarYdBmagR[2]*hamilR[2])*rdy2R+BstarXdBmagR[0]*hamilR[5]*rdx2R)+13.85640646055102*dApardtPrev[2]*q_))/m_; 
  alphaUpR[3] = -(0.3061862178478971*BstarXdBmagR[0]*hamilR[8]*rdx2R)/m_; 
  alphaUpR[4] = -(0.1767766952966368*(1.732050807568877*((BstarYdBmagR[2]*hamilR[5]+hamilR[2]*BstarYdBmagR[5])*rdy2R+BstarXdBmagR[1]*hamilR[5]*rdx2R)+8.0*dApardtPrev[3]*q_))/m_; 
  alphaUpR[5] = -(0.3061862178478971*BstarXdBmagR[1]*hamilR[8]*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*((3.0*hamilR[5]*BstarYdBmagR[6]-1.732050807568877*BstarYdBmagR[1]*hamilR[5]+3.0*hamilR[2]*BstarYdBmagR[3]-1.732050807568877*BstarYdBmagR[0]*hamilR[2])*rdy2R-1.732050807568877*BstarXdBmagR[0]*hamilR[1]*rdx2R-8.0*dApardtPrev[0]*q_))/m_; 

  double incr[16]; 
#if upwindType == SURFAVG 
  if (alphaSurfAvgR>0) { 
  incr[0] = 0.125*(2.449489742783178*(alphaR[4]*fL[11]+alphaR[2]*fL[7]+alphaR[1]*fL[6])+1.414213562373095*alphaR[4]*fL[5]+2.449489742783178*alphaR[0]*fL[3]+1.414213562373095*(alphaR[2]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[1] = 0.125*(2.449489742783178*(alphaR[2]*fL[11]+alphaR[4]*fL[7]+alphaR[0]*fL[6])+1.414213562373095*(alphaR[2]*fL[5]+fL[2]*alphaR[4])+2.449489742783178*alphaR[1]*fL[3]+1.414213562373095*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[2] = 0.125*(2.449489742783178*(alphaR[1]*fL[11]+alphaR[0]*fL[7]+alphaR[4]*fL[6])+1.414213562373095*(alphaR[1]*fL[5]+fL[1]*alphaR[4])+2.449489742783178*alphaR[2]*fL[3]+1.414213562373095*(alphaR[0]*fL[2]+fL[0]*alphaR[2])); 
  incr[3] = -0.125*(4.242640687119286*(alphaR[4]*fL[11]+alphaR[2]*fL[7]+alphaR[1]*fL[6])+2.449489742783178*alphaR[4]*fL[5]+4.242640687119286*alphaR[0]*fL[3]+2.449489742783178*(alphaR[2]*fL[2]+alphaR[1]*fL[1]+alphaR[0]*fL[0])); 
  incr[4] = 0.125*(2.449489742783178*(alphaR[4]*fL[15]+alphaR[2]*fL[14]+alphaR[1]*fL[13])+1.414213562373095*alphaR[4]*fL[12]+2.449489742783178*alphaR[0]*fL[10]+1.414213562373095*(alphaR[2]*fL[9]+alphaR[1]*fL[8]+alphaR[0]*fL[4])); 
  incr[5] = 0.125*(2.449489742783178*(alphaR[0]*fL[11]+alphaR[1]*fL[7]+alphaR[2]*fL[6])+1.414213562373095*alphaR[0]*fL[5]+(2.449489742783178*fL[3]+1.414213562373095*fL[0])*alphaR[4]+1.414213562373095*(alphaR[1]*fL[2]+fL[1]*alphaR[2])); 
  incr[6] = -0.125*(4.242640687119286*(alphaR[2]*fL[11]+alphaR[4]*fL[7]+alphaR[0]*fL[6])+2.449489742783178*(alphaR[2]*fL[5]+fL[2]*alphaR[4])+4.242640687119286*alphaR[1]*fL[3]+2.449489742783178*(alphaR[0]*fL[1]+fL[0]*alphaR[1])); 
  incr[7] = -0.125*(4.242640687119286*(alphaR[1]*fL[11]+alphaR[0]*fL[7]+alphaR[4]*fL[6])+2.449489742783178*(alphaR[1]*fL[5]+fL[1]*alphaR[4])+4.242640687119286*alphaR[2]*fL[3]+2.449489742783178*(alphaR[0]*fL[2]+fL[0]*alphaR[2])); 
  incr[8] = 0.125*(2.449489742783178*(alphaR[2]*fL[15]+alphaR[4]*fL[14]+alphaR[0]*fL[13])+1.414213562373095*alphaR[2]*fL[12]+2.449489742783178*alphaR[1]*fL[10]+1.414213562373095*(alphaR[4]*fL[9]+alphaR[0]*fL[8]+alphaR[1]*fL[4])); 
  incr[9] = 0.125*(2.449489742783178*(alphaR[1]*fL[15]+alphaR[0]*fL[14]+alphaR[4]*fL[13])+1.414213562373095*alphaR[1]*fL[12]+2.449489742783178*alphaR[2]*fL[10]+1.414213562373095*(alphaR[0]*fL[9]+alphaR[4]*fL[8]+alphaR[2]*fL[4])); 
  incr[10] = -0.125*(4.242640687119286*(alphaR[4]*fL[15]+alphaR[2]*fL[14]+alphaR[1]*fL[13])+2.449489742783178*alphaR[4]*fL[12]+4.242640687119286*alphaR[0]*fL[10]+2.449489742783178*(alphaR[2]*fL[9]+alphaR[1]*fL[8]+alphaR[0]*fL[4])); 
  incr[11] = -0.125*(4.242640687119286*(alphaR[0]*fL[11]+alphaR[1]*fL[7]+alphaR[2]*fL[6])+2.449489742783178*alphaR[0]*fL[5]+(4.242640687119286*fL[3]+2.449489742783178*fL[0])*alphaR[4]+2.449489742783178*(alphaR[1]*fL[2]+fL[1]*alphaR[2])); 
  incr[12] = 0.125*(2.449489742783178*(alphaR[0]*fL[15]+alphaR[1]*fL[14]+alphaR[2]*fL[13])+1.414213562373095*alphaR[0]*fL[12]+2.449489742783178*alphaR[4]*fL[10]+1.414213562373095*(alphaR[1]*fL[9]+alphaR[2]*fL[8]+alphaR[4]*fL[4])); 
  incr[13] = -0.125*(4.242640687119286*(alphaR[2]*fL[15]+alphaR[4]*fL[14]+alphaR[0]*fL[13])+2.449489742783178*alphaR[2]*fL[12]+4.242640687119286*alphaR[1]*fL[10]+2.449489742783178*(alphaR[4]*fL[9]+alphaR[0]*fL[8]+alphaR[1]*fL[4])); 
  incr[14] = -0.125*(4.242640687119286*(alphaR[1]*fL[15]+alphaR[0]*fL[14]+alphaR[4]*fL[13])+2.449489742783178*alphaR[1]*fL[12]+4.242640687119286*alphaR[2]*fL[10]+2.449489742783178*(alphaR[0]*fL[9]+alphaR[4]*fL[8]+alphaR[2]*fL[4])); 
  incr[15] = -0.125*(4.242640687119286*(alphaR[0]*fL[15]+alphaR[1]*fL[14]+alphaR[2]*fL[13])+2.449489742783178*alphaR[0]*fL[12]+4.242640687119286*alphaR[4]*fL[10]+2.449489742783178*(alphaR[1]*fL[9]+alphaR[2]*fL[8]+alphaR[4]*fL[4])); 
  } else { 
  incr[0] = -0.125*(2.449489742783178*(alphaR[4]*fR[11]+alphaR[2]*fR[7]+alphaR[1]*fR[6])-1.414213562373095*alphaR[4]*fR[5]+2.449489742783178*alphaR[0]*fR[3]-1.414213562373095*(alphaR[2]*fR[2]+alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[1] = -0.125*(2.449489742783178*(alphaR[2]*fR[11]+alphaR[4]*fR[7]+alphaR[0]*fR[6])-1.414213562373095*(alphaR[2]*fR[5]+fR[2]*alphaR[4])+2.449489742783178*alphaR[1]*fR[3]-1.414213562373095*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[2] = -0.125*(2.449489742783178*(alphaR[1]*fR[11]+alphaR[0]*fR[7]+alphaR[4]*fR[6])-1.414213562373095*(alphaR[1]*fR[5]+fR[1]*alphaR[4])+2.449489742783178*alphaR[2]*fR[3]-1.414213562373095*(alphaR[0]*fR[2]+fR[0]*alphaR[2])); 
  incr[3] = 0.125*(4.242640687119286*(alphaR[4]*fR[11]+alphaR[2]*fR[7]+alphaR[1]*fR[6])-2.449489742783178*alphaR[4]*fR[5]+4.242640687119286*alphaR[0]*fR[3]-2.449489742783178*(alphaR[2]*fR[2]+alphaR[1]*fR[1]+alphaR[0]*fR[0])); 
  incr[4] = -0.125*(2.449489742783178*(alphaR[4]*fR[15]+alphaR[2]*fR[14]+alphaR[1]*fR[13])-1.414213562373095*alphaR[4]*fR[12]+2.449489742783178*alphaR[0]*fR[10]-1.414213562373095*(alphaR[2]*fR[9]+alphaR[1]*fR[8]+alphaR[0]*fR[4])); 
  incr[5] = -0.125*(2.449489742783178*(alphaR[0]*fR[11]+alphaR[1]*fR[7]+alphaR[2]*fR[6])-1.414213562373095*alphaR[0]*fR[5]+(2.449489742783178*fR[3]-1.414213562373095*fR[0])*alphaR[4]-1.414213562373095*(alphaR[1]*fR[2]+fR[1]*alphaR[2])); 
  incr[6] = 0.125*(4.242640687119286*(alphaR[2]*fR[11]+alphaR[4]*fR[7]+alphaR[0]*fR[6])-2.449489742783178*(alphaR[2]*fR[5]+fR[2]*alphaR[4])+4.242640687119286*alphaR[1]*fR[3]-2.449489742783178*(alphaR[0]*fR[1]+fR[0]*alphaR[1])); 
  incr[7] = 0.125*(4.242640687119286*(alphaR[1]*fR[11]+alphaR[0]*fR[7]+alphaR[4]*fR[6])-2.449489742783178*(alphaR[1]*fR[5]+fR[1]*alphaR[4])+4.242640687119286*alphaR[2]*fR[3]-2.449489742783178*(alphaR[0]*fR[2]+fR[0]*alphaR[2])); 
  incr[8] = -0.125*(2.449489742783178*(alphaR[2]*fR[15]+alphaR[4]*fR[14]+alphaR[0]*fR[13])-1.414213562373095*alphaR[2]*fR[12]+2.449489742783178*alphaR[1]*fR[10]-1.414213562373095*(alphaR[4]*fR[9]+alphaR[0]*fR[8]+alphaR[1]*fR[4])); 
  incr[9] = -0.125*(2.449489742783178*(alphaR[1]*fR[15]+alphaR[0]*fR[14]+alphaR[4]*fR[13])-1.414213562373095*alphaR[1]*fR[12]+2.449489742783178*alphaR[2]*fR[10]-1.414213562373095*(alphaR[0]*fR[9]+alphaR[4]*fR[8]+alphaR[2]*fR[4])); 
  incr[10] = 0.125*(4.242640687119286*(alphaR[4]*fR[15]+alphaR[2]*fR[14]+alphaR[1]*fR[13])-2.449489742783178*alphaR[4]*fR[12]+4.242640687119286*alphaR[0]*fR[10]-2.449489742783178*(alphaR[2]*fR[9]+alphaR[1]*fR[8]+alphaR[0]*fR[4])); 
  incr[11] = 0.125*(4.242640687119286*(alphaR[0]*fR[11]+alphaR[1]*fR[7]+alphaR[2]*fR[6])-2.449489742783178*alphaR[0]*fR[5]+(4.242640687119286*fR[3]-2.449489742783178*fR[0])*alphaR[4]-2.449489742783178*(alphaR[1]*fR[2]+fR[1]*alphaR[2])); 
  incr[12] = -0.125*(2.449489742783178*(alphaR[0]*fR[15]+alphaR[1]*fR[14]+alphaR[2]*fR[13])-1.414213562373095*alphaR[0]*fR[12]+2.449489742783178*alphaR[4]*fR[10]-1.414213562373095*(alphaR[1]*fR[9]+alphaR[2]*fR[8]+alphaR[4]*fR[4])); 
  incr[13] = 0.125*(4.242640687119286*(alphaR[2]*fR[15]+alphaR[4]*fR[14]+alphaR[0]*fR[13])-2.449489742783178*alphaR[2]*fR[12]+4.242640687119286*alphaR[1]*fR[10]-2.449489742783178*(alphaR[4]*fR[9]+alphaR[0]*fR[8]+alphaR[1]*fR[4])); 
  incr[14] = 0.125*(4.242640687119286*(alphaR[1]*fR[15]+alphaR[0]*fR[14]+alphaR[4]*fR[13])-2.449489742783178*alphaR[1]*fR[12]+4.242640687119286*alphaR[2]*fR[10]-2.449489742783178*(alphaR[0]*fR[9]+alphaR[4]*fR[8]+alphaR[2]*fR[4])); 
  incr[15] = 0.125*(4.242640687119286*(alphaR[0]*fR[15]+alphaR[1]*fR[14]+alphaR[2]*fR[13])-2.449489742783178*alphaR[0]*fR[12]+4.242640687119286*alphaR[4]*fR[10]-2.449489742783178*(alphaR[1]*fR[9]+alphaR[2]*fR[8]+alphaR[4]*fR[4])); 
  }
#elif upwindType == QUAD 
  double alphaOrdR;
  double fUpOrd[8];
  alphaOrdR = 0.3535533905932737*(alphaUpR[5]+alphaUpR[4])-0.3535533905932737*(alphaUpR[3]+alphaUpR[2]+alphaUpR[1])+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[0] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])-0.25*(fR[9]+fR[8])+0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14]+fR[13])+0.4330127018922193*(fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[1]+alphaUpR[0])-0.3535533905932737*(alphaUpR[5]+alphaUpR[4]+alphaUpR[3]+alphaUpR[2]); 
  fUpOrd[1] = 0.5*((0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)-0.4330127018922193*(fR[15]+fR[14])+0.4330127018922193*(fL[15]+fL[14]+fR[13])-0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaUpR[5]-0.3535533905932737*(alphaUpR[4]+alphaUpR[3])+0.3535533905932737*alphaUpR[2]-0.3535533905932737*alphaUpR[1]+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[2] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14])+0.4330127018922193*(fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*(fR[5]+fR[4])-0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14])-0.4330127018922193*(fL[14]+fR[13])+0.4330127018922193*fL[13]+0.25*(fR[12]+fL[12])+0.4330127018922193*(fR[11]+fR[10])-0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaUpR[5])+0.3535533905932737*alphaUpR[4]-0.3535533905932737*alphaUpR[3]+0.3535533905932737*(alphaUpR[2]+alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[3] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]+fR[13]+fL[13]))+0.25*fR[12]-0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11])-0.4330127018922193*(fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])+0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*fR[5]+0.25*(fL[5]+fR[4])-0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*(fR[2]+fR[1]+fR[0])+0.25*(fL[2]+fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14]+fR[13])-0.4330127018922193*(fL[15]+fL[14]+fL[13])-0.25*(fR[12]+fL[12])-0.4330127018922193*fR[11]+0.4330127018922193*(fL[11]+fR[10])-0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9]+fR[8]+fL[8])-0.4330127018922193*(fR[7]+fR[6])+0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5])-0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2]+fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*alphaUpR[5])+0.3535533905932737*(alphaUpR[4]+alphaUpR[3])-0.3535533905932737*(alphaUpR[2]+alphaUpR[1])+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[4] = 0.5*((0.4330127018922193*(fR[15]+fL[15])-0.4330127018922193*(fR[14]+fL[14]+fR[13]+fL[13])-0.25*fR[12]+0.25*fL[12]+0.4330127018922193*(fR[11]+fL[11]+fR[10]+fL[10])+0.25*(fR[9]+fR[8])-0.25*(fL[9]+fL[8])-0.4330127018922193*(fR[7]+fL[7]+fR[6]+fL[6])-0.25*(fR[5]+fR[4])+0.25*(fL[5]+fL[4])+0.4330127018922193*(fR[3]+fL[3])+0.25*(fR[2]+fR[1])-0.25*(fL[2]+fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)-0.4330127018922193*fR[15]+0.4330127018922193*(fL[15]+fR[14]+fR[13])-0.4330127018922193*(fL[14]+fL[13])+0.25*(fR[12]+fL[12])-0.4330127018922193*(fR[11]+fR[10])+0.4330127018922193*(fL[11]+fL[10])-0.25*(fR[9]+fL[9]+fR[8]+fL[8])+0.4330127018922193*(fR[7]+fR[6])-0.4330127018922193*(fL[7]+fL[6])+0.25*(fR[5]+fL[5]+fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2]+fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*alphaUpR[5]-0.3535533905932737*alphaUpR[4]+0.3535533905932737*alphaUpR[3]-0.3535533905932737*alphaUpR[2]+0.3535533905932737*(alphaUpR[1]+alphaUpR[0]); 
  fUpOrd[5] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]+fR[14]+fL[14]))+0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])+0.25*fR[9]-0.25*(fL[9]+fR[8])+0.25*fL[8]-0.4330127018922193*(fR[7]+fL[7])+0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])+0.25*fR[2]-0.25*(fL[2]+fR[1]+fR[0])+0.25*(fL[1]+fL[0]))*sgn(alphaOrdR)+0.4330127018922193*(fR[15]+fR[14])-0.4330127018922193*(fL[15]+fL[14]+fR[13])+0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]-0.25*(fR[9]+fL[9])+0.25*(fR[8]+fL[8])+0.4330127018922193*fR[7]-0.4330127018922193*(fL[7]+fR[6])+0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]-0.25*(fR[2]+fL[2])+0.25*(fR[1]+fL[1]+fR[0]+fL[0])); 
  alphaOrdR = (-0.3535533905932737*(alphaUpR[5]+alphaUpR[4]))+0.3535533905932737*(alphaUpR[3]+alphaUpR[2])-0.3535533905932737*alphaUpR[1]+0.3535533905932737*alphaUpR[0]; 
  fUpOrd[6] = 0.5*(((-0.4330127018922193*(fR[15]+fL[15]))+0.4330127018922193*(fR[14]+fL[14])-0.4330127018922193*(fR[13]+fL[13])+0.25*fR[12]-0.25*fL[12]-0.4330127018922193*(fR[11]+fL[11])+0.4330127018922193*(fR[10]+fL[10])-0.25*fR[9]+0.25*(fL[9]+fR[8])-0.25*fL[8]+0.4330127018922193*(fR[7]+fL[7])-0.4330127018922193*(fR[6]+fL[6])+0.25*fR[5]-0.25*(fL[5]+fR[4])+0.25*fL[4]+0.4330127018922193*(fR[3]+fL[3])-0.25*fR[2]+0.25*(fL[2]+fR[1])-0.25*(fL[1]+fR[0])+0.25*fL[0])*sgn(alphaOrdR)+0.4330127018922193*fR[15]-0.4330127018922193*(fL[15]+fR[14])+0.4330127018922193*(fL[14]+fR[13])-0.4330127018922193*fL[13]-0.25*(fR[12]+fL[12])+0.4330127018922193*fR[11]-0.4330127018922193*(fL[11]+fR[10])+0.4330127018922193*fL[10]+0.25*(fR[9]+fL[9])-0.25*(fR[8]+fL[8])-0.4330127018922193*fR[7]+0.4330127018922193*(fL[7]+fR[6])-0.4330127018922193*fL[6]-0.25*(fR[5]+fL[5])+0.25*(fR[4]+fL[4])-0.4330127018922193*fR[3]+0.4330127018922193*fL[3]+0.25*(fR[2]+fL[2])-0.25*(fR[1]+fL[1])+0.25*(fR[0]+fL[0])); 
  alphaOrdR = 0.3535533905932737*(alphaUpR[5]+alphaUpR[4]+alphaUpR[3]+alphaUpR[2]+alphaUpR[1]+alphaUpR[0]); 
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

  incr[0] = 0.25*(alphaR[4]*fUp[4]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[1] = 0.25*(alphaR[2]*fUp[4]+fUp[2]*alphaR[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[2] = 0.25*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[3] = -0.4330127018922193*(alphaR[4]*fUp[4]+alphaR[2]*fUp[2]+alphaR[1]*fUp[1]+alphaR[0]*fUp[0]); 
  incr[4] = 0.25*(alphaR[4]*fUp[7]+alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+alphaR[0]*fUp[3]); 
  incr[5] = 0.25*(alphaR[0]*fUp[4]+fUp[0]*alphaR[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[6] = -0.4330127018922193*(alphaR[2]*fUp[4]+fUp[2]*alphaR[4]+alphaR[0]*fUp[1]+fUp[0]*alphaR[1]); 
  incr[7] = -0.4330127018922193*(alphaR[1]*fUp[4]+fUp[1]*alphaR[4]+alphaR[0]*fUp[2]+fUp[0]*alphaR[2]); 
  incr[8] = 0.25*(alphaR[2]*fUp[7]+alphaR[4]*fUp[6]+alphaR[0]*fUp[5]+alphaR[1]*fUp[3]); 
  incr[9] = 0.25*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[4]*fUp[5]+alphaR[2]*fUp[3]); 
  incr[10] = -0.4330127018922193*(alphaR[4]*fUp[7]+alphaR[2]*fUp[6]+alphaR[1]*fUp[5]+alphaR[0]*fUp[3]); 
  incr[11] = -0.4330127018922193*(alphaR[0]*fUp[4]+fUp[0]*alphaR[4]+alphaR[1]*fUp[2]+fUp[1]*alphaR[2]); 
  incr[12] = 0.25*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]+fUp[3]*alphaR[4]); 
  incr[13] = -0.4330127018922193*(alphaR[2]*fUp[7]+alphaR[4]*fUp[6]+alphaR[0]*fUp[5]+alphaR[1]*fUp[3]); 
  incr[14] = -0.4330127018922193*(alphaR[1]*fUp[7]+alphaR[0]*fUp[6]+alphaR[4]*fUp[5]+alphaR[2]*fUp[3]); 
  incr[15] = -0.4330127018922193*(alphaR[0]*fUp[7]+alphaR[1]*fUp[6]+alphaR[2]*fUp[5]+fUp[3]*alphaR[4]); 

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
