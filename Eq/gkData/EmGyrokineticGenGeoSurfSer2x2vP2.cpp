#include <GyrokineticModDecl.h>
double EmGyrokineticGenGeoSurf2x2vSer_xL_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilL[48]; 
  hamilL[0] = (0.6666666666666666*(3.0*rdvpar2SqL*(m_*wvparSqL+bmag[0]*wmuL+phi[0]*q_)+m_))/rdvpar2SqL; 
  hamilL[1] = 2.0*phi[1]*q_; 
  hamilL[2] = 2.0*phi[2]*q_; 
  hamilL[3] = (2.309401076758503*m_*wvparL)/rdvpar2L; 
  hamilL[4] = (1.154700538379252*bmag[0])/rdmu2L; 
  hamilL[5] = 2.0*phi[3]*q_; 
  hamilL[11] = 2.0*phi[4]*q_; 
  hamilL[12] = 2.0*phi[5]*q_; 
  hamilL[13] = (0.5962847939999438*m_)/rdvpar2SqL; 
  hamilL[19] = 2.0*phi[6]*q_; 
  hamilL[20] = 2.0*phi[7]*q_; 

  double BstarXdBmagL[48]; 
  BstarXdBmagL[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2L; 
  BstarXdBmagL[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2L; 
  BstarXdBmagL[2] = 1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[5]*rdy2L; 
  BstarXdBmagL[5] = 1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[7]*rdy2L; 
  BstarXdBmagL[11] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[6]*rdy2L; 

  double alphaL[20]; 
  alphaL[0] = -(0.1767766952966368*(b_z[0]*jacobTotInv[0]*(3.872983346207417*hamilL[19]+3.0*hamilL[5]+1.732050807568877*hamilL[2])*m_*rdy2L+hamilL[3]*((-3.872983346207417*BstarXdBmagL[11])-3.0*BstarXdBmagL[1]-1.732050807568877*BstarXdBmagL[0])*q_*rdvpar2L))/(m_*q_); 
  alphaL[1] = -(0.1767766952966368*(b_z[0]*jacobTotInv[0]*(6.708203932499369*hamilL[20]+3.872983346207417*hamilL[12])*m_*rdy2L+hamilL[3]*((-3.0*BstarXdBmagL[5])-1.732050807568877*BstarXdBmagL[2])*q_*rdvpar2L))/(m_*q_); 
  alphaL[2] = (0.1767766952966368*(8.660254037844386*BstarXdBmagL[11]+6.708203932499369*BstarXdBmagL[1]+3.872983346207417*BstarXdBmagL[0])*hamilL[13]*rdvpar2L)/m_; 
  alphaL[4] = (0.1767766952966368*(6.708203932499369*BstarXdBmagL[5]+3.872983346207417*BstarXdBmagL[2])*hamilL[13]*rdvpar2L)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgL = -(0.03125*((3.872983346207417*b_z[0]*jacobTotInv[0]*hamilL[19]+3.0*b_z[0]*jacobTotInv[0]*hamilL[5]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilL[2])*m_*rdy2L+(((-3.0*BstarXdBmagL[1])-1.732050807568877*BstarXdBmagL[0])*hamilL[3]-3.872983346207417*hamilL[3]*BstarXdBmagL[11])*q_*rdvpar2L))/(m_*q_); 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[11]+fL[11])+1.732050807568877*(fL[1]-1.0*fR[1])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[19]+fL[19])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[2]+fL[2])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[21]+fL[21])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[25]+fL[25])+3.0*(fL[8]-1.0*fR[8]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[32]+fL[32])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[7]+fL[7]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[35]+fL[35])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[9]+fL[9]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[37]+fL[37])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[20]-1.0*(8.660254037844387*fL[20]+5.0*(fR[12]+fL[12]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[23]-1.0*(8.660254037844387*fL[23]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[28]-1.0*(8.660254037844387*fL[28]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[44]+fL[44])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[18]+fL[18])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[22]+fL[22]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[24]+fL[24]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[36]-1.0*(8.660254037844387*fL[36]+5.0*(fR[26]+fL[26]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[39]-1.0*(8.660254037844387*fL[39]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[29]+fL[29]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[38]+fL[38]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[40]+fL[40]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[43]+fL[43]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[11]-6.324555320336761*fL[11]-4.898979485566357*(fR[1]+fL[1])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaL[4]*fAvg[4]+alphaL[2]*fAvg[2]+alphaL[1]*fAvg[1]+alphaL[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[19]-94.8683298050514*fL[19]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*(alphaL[4]*fAvg[11]+alphaL[1]*fAvg[7])-21.21320343559643*(alphaL[2]*fAvg[4]+fAvg[2]*alphaL[4]+alphaL[0]*fAvg[1]+fAvg[0]*alphaL[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[21]-94.8683298050514*fL[21]-73.48469228349535*(fR[6]+fL[6])+42.42640687119286*fR[3]-42.42640687119286*fL[3])*amax-18.97366596101028*(alphaL[4]*fAvg[12]+alphaL[2]*fAvg[8])-21.21320343559643*(alphaL[1]*fAvg[4]+fAvg[1]*alphaL[4]+alphaL[0]*fAvg[2]+fAvg[0]*alphaL[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[25]-18.97366596101028*fL[25]-14.69693845669907*(fR[8]+fL[8])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaL[4]*fAvg[10]+alphaL[2]*fAvg[6]+alphaL[1]*fAvg[5]+alphaL[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[32]-94.86832980505142*fL[32]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[7]-42.42640687119286*fL[7])*amax-18.97366596101028*(alphaL[2]*fAvg[12]+alphaL[1]*fAvg[11]+alphaL[4]*(fAvg[8]+fAvg[7]))-21.21320343559643*(alphaL[0]*fAvg[4]+fAvg[0]*alphaL[4]+alphaL[1]*fAvg[2]+fAvg[1]*alphaL[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[35]-94.86832980505142*fL[35]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*(alphaL[4]*fAvg[17]+alphaL[1]*fAvg[13])-21.21320343559643*(alphaL[2]*fAvg[10]+alphaL[4]*fAvg[6]+alphaL[0]*fAvg[5]+alphaL[1]*fAvg[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[37]-94.86832980505142*fL[37]-73.48469228349535*(fR[17]+fL[17])+42.42640687119286*fR[10]-42.42640687119286*fL[10])*amax-18.97366596101028*(alphaL[4]*fAvg[18]+alphaL[2]*fAvg[14])-21.21320343559643*(alphaL[1]*fAvg[10]+alphaL[0]*fAvg[6]+alphaL[4]*fAvg[5]+alphaL[2]*fAvg[3])); 
  Ghat[7] = 0.008333333333333333*((73.48469228349536*(fR[20]+fL[20])-42.42640687119286*fR[12]+42.42640687119286*fL[12])*amax+21.21320343559643*alphaL[2]*fAvg[11]+21.21320343559643*alphaL[0]*fAvg[7]+18.97366596101028*(alphaL[4]*fAvg[4]+alphaL[1]*fAvg[1])); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[23]+fL[23])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaL[1]*fAvg[12]+21.21320343559643*alphaL[0]*fAvg[8]+18.97366596101028*(alphaL[4]*fAvg[4]+alphaL[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[28]+fL[28])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaL[4]*fAvg[19]+21.21320343559643*(alphaL[2]*fAvg[16]+alphaL[1]*fAvg[15])+21.21320343559643*alphaL[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[44]-94.8683298050514*fL[44]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[18]-42.42640687119286*fL[18])*amax-18.97366596101028*(alphaL[2]*fAvg[18]+alphaL[1]*fAvg[17]+alphaL[4]*(fAvg[14]+fAvg[13]))-21.21320343559643*(alphaL[0]*fAvg[10]+alphaL[1]*fAvg[6]+alphaL[2]*fAvg[5]+fAvg[3]*alphaL[4])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[33]+fL[33])-42.42640687119286*fR[22]+42.42640687119286*fL[22])*amax+16.97056274847715*alphaL[4]*fAvg[12]+21.21320343559643*alphaL[0]*fAvg[11]+21.21320343559643*alphaL[2]*fAvg[7]+18.97366596101028*(alphaL[1]*fAvg[4]+fAvg[1]*alphaL[4])); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[24]+42.42640687119286*fL[24])*amax+21.21320343559643*alphaL[0]*fAvg[12]+16.97056274847715*alphaL[4]*fAvg[11]+21.21320343559643*alphaL[1]*fAvg[8]+18.97366596101028*(alphaL[2]*fAvg[4]+fAvg[2]*alphaL[4])); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[36]+fL[36])-42.42640687119286*fR[26]+42.42640687119286*fL[26])*amax+21.21320343559643*alphaL[2]*fAvg[17]+21.21320343559643*alphaL[0]*fAvg[13]+18.97366596101028*(alphaL[4]*fAvg[10]+alphaL[1]*fAvg[5])); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[39]+fL[39])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaL[1]*fAvg[18]+21.21320343559643*alphaL[0]*fAvg[14]+18.97366596101028*(alphaL[4]*fAvg[10]+alphaL[2]*fAvg[6])); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+21.21320343559643*alphaL[2]*fAvg[19]+21.21320343559643*(alphaL[4]*fAvg[16]+alphaL[0]*fAvg[15])+21.21320343559643*alphaL[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaL[1]*fAvg[19]+21.21320343559643*(alphaL[0]*fAvg[16]+alphaL[4]*fAvg[15])+21.21320343559643*alphaL[2]*fAvg[9]); 
  Ghat[17] = 0.008333333333333333*((73.48469228349536*(fR[45]+fL[45])-42.42640687119286*fR[38]+42.42640687119286*fL[38])*amax+16.97056274847715*alphaL[4]*fAvg[18]+21.21320343559643*alphaL[0]*fAvg[17]+21.21320343559643*alphaL[2]*fAvg[13]+18.97366596101028*(alphaL[1]*fAvg[10]+alphaL[4]*fAvg[5])); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[40]+42.42640687119286*fL[40])*amax+21.21320343559643*alphaL[0]*fAvg[18]+16.97056274847715*alphaL[4]*fAvg[17]+21.21320343559643*alphaL[1]*fAvg[14]+18.97366596101028*(alphaL[2]*fAvg[10]+alphaL[4]*fAvg[6])); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[43]+42.42640687119286*fL[43])*amax+21.21320343559643*alphaL[0]*fAvg[19]+21.21320343559643*(alphaL[1]*fAvg[16]+alphaL[2]*fAvg[15])+21.21320343559643*alphaL[4]*fAvg[9]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = -1.224744871391589*Ghat[0]; 
  incr[2] = 0.7071067811865475*Ghat[1]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = -1.224744871391589*Ghat[2]; 
  incr[7] = 0.7071067811865475*Ghat[4]; 
  incr[8] = -1.224744871391589*Ghat[3]; 
  incr[9] = 0.7071067811865475*Ghat[5]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 1.58113883008419*Ghat[0]; 
  incr[12] = 0.7071067811865475*Ghat[7]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = -1.224744871391589*Ghat[6]; 
  incr[18] = 0.7071067811865475*Ghat[10]; 
  incr[19] = 1.58113883008419*Ghat[1]; 
  incr[20] = -1.224744871391589*Ghat[7]; 
  incr[21] = 1.58113883008419*Ghat[2]; 
  incr[22] = 0.7071067811865475*Ghat[11]; 
  incr[23] = -1.224744871391589*Ghat[8]; 
  incr[24] = 0.7071067811865475*Ghat[12]; 
  incr[25] = 1.58113883008419*Ghat[3]; 
  incr[26] = 0.7071067811865475*Ghat[13]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = -1.224744871391589*Ghat[9]; 
  incr[29] = 0.7071067811865475*Ghat[15]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = 1.58113883008419*Ghat[4]; 
  incr[33] = -1.224744871391589*Ghat[11]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = 1.58113883008419*Ghat[5]; 
  incr[36] = -1.224744871391589*Ghat[13]; 
  incr[37] = 1.58113883008419*Ghat[6]; 
  incr[38] = 0.7071067811865475*Ghat[17]; 
  incr[39] = -1.224744871391589*Ghat[14]; 
  incr[40] = 0.7071067811865475*Ghat[18]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = -1.224744871391589*Ghat[16]; 
  incr[43] = 0.7071067811865475*Ghat[19]; 
  incr[44] = 1.58113883008419*Ghat[10]; 
  incr[45] = -1.224744871391589*Ghat[17]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

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
  outR[16] += incr[16]*rdx2R; 
  outR[17] += incr[17]*rdx2R; 
  outR[18] += incr[18]*rdx2R; 
  outR[19] += incr[19]*rdx2R; 
  outR[20] += incr[20]*rdx2R; 
  outR[21] += incr[21]*rdx2R; 
  outR[22] += incr[22]*rdx2R; 
  outR[23] += incr[23]*rdx2R; 
  outR[24] += incr[24]*rdx2R; 
  outR[25] += incr[25]*rdx2R; 
  outR[26] += incr[26]*rdx2R; 
  outR[27] += incr[27]*rdx2R; 
  outR[28] += incr[28]*rdx2R; 
  outR[29] += incr[29]*rdx2R; 
  outR[30] += incr[30]*rdx2R; 
  outR[31] += incr[31]*rdx2R; 
  outR[32] += incr[32]*rdx2R; 
  outR[33] += incr[33]*rdx2R; 
  outR[34] += incr[34]*rdx2R; 
  outR[35] += incr[35]*rdx2R; 
  outR[36] += incr[36]*rdx2R; 
  outR[37] += incr[37]*rdx2R; 
  outR[38] += incr[38]*rdx2R; 
  outR[39] += incr[39]*rdx2R; 
  outR[40] += incr[40]*rdx2R; 
  outR[41] += incr[41]*rdx2R; 
  outR[42] += incr[42]*rdx2R; 
  outR[43] += incr[43]*rdx2R; 
  outR[44] += incr[44]*rdx2R; 
  outR[45] += incr[45]*rdx2R; 
  outR[46] += incr[46]*rdx2R; 
  outR[47] += incr[47]*rdx2R; 

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
  outL[11] += -1.0*incr[11]*rdx2L; 
  outL[12] += -1.0*incr[12]*rdx2L; 
  outL[13] += -1.0*incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 
  outL[16] += incr[16]*rdx2L; 
  outL[17] += incr[17]*rdx2L; 
  outL[18] += -1.0*incr[18]*rdx2L; 
  outL[19] += -1.0*incr[19]*rdx2L; 
  outL[20] += incr[20]*rdx2L; 
  outL[21] += -1.0*incr[21]*rdx2L; 
  outL[22] += -1.0*incr[22]*rdx2L; 
  outL[23] += incr[23]*rdx2L; 
  outL[24] += -1.0*incr[24]*rdx2L; 
  outL[25] += -1.0*incr[25]*rdx2L; 
  outL[26] += -1.0*incr[26]*rdx2L; 
  outL[27] += -1.0*incr[27]*rdx2L; 
  outL[28] += incr[28]*rdx2L; 
  outL[29] += -1.0*incr[29]*rdx2L; 
  outL[30] += -1.0*incr[30]*rdx2L; 
  outL[31] += incr[31]*rdx2L; 
  outL[32] += -1.0*incr[32]*rdx2L; 
  outL[33] += incr[33]*rdx2L; 
  outL[34] += incr[34]*rdx2L; 
  outL[35] += -1.0*incr[35]*rdx2L; 
  outL[36] += incr[36]*rdx2L; 
  outL[37] += -1.0*incr[37]*rdx2L; 
  outL[38] += -1.0*incr[38]*rdx2L; 
  outL[39] += incr[39]*rdx2L; 
  outL[40] += -1.0*incr[40]*rdx2L; 
  outL[41] += incr[41]*rdx2L; 
  outL[42] += incr[42]*rdx2L; 
  outL[43] += -1.0*incr[43]*rdx2L; 
  outL[44] += -1.0*incr[44]*rdx2L; 
  outL[45] += incr[45]*rdx2L; 
  outL[46] += incr[46]*rdx2L; 
  outL[47] += incr[47]*rdx2L; 

  return std::abs(alphaSurfAvgL);
} 
double EmGyrokineticGenGeoSurf2x2vSer_xR_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[11] = 2.0*phi[4]*q_; 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 

  double BstarXdBmagR[48]; 
  BstarXdBmagR[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2R; 
  BstarXdBmagR[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2R; 
  BstarXdBmagR[2] = 1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[5]*rdy2R; 
  BstarXdBmagR[5] = 1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[7]*rdy2R; 
  BstarXdBmagR[11] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[6]*rdy2R; 

  double alphaR[20]; 
  alphaR[0] = -(0.1767766952966368*(b_z[0]*jacobTotInv[0]*(3.872983346207417*hamilR[19]-3.0*hamilR[5]+1.732050807568877*hamilR[2])*m_*rdy2R+hamilR[3]*((-3.872983346207417*BstarXdBmagR[11])+3.0*BstarXdBmagR[1]-1.732050807568877*BstarXdBmagR[0])*q_*rdvpar2R))/(m_*q_); 
  alphaR[1] = (0.1767766952966368*(b_z[0]*jacobTotInv[0]*(6.708203932499369*hamilR[20]-3.872983346207417*hamilR[12])*m_*rdy2R+hamilR[3]*(1.732050807568877*BstarXdBmagR[2]-3.0*BstarXdBmagR[5])*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = (0.1767766952966368*(8.660254037844386*BstarXdBmagR[11]-6.708203932499369*BstarXdBmagR[1]+3.872983346207417*BstarXdBmagR[0])*hamilR[13]*rdvpar2R)/m_; 
  alphaR[4] = -(0.1767766952966368*(6.708203932499369*BstarXdBmagR[5]-3.872983346207417*BstarXdBmagR[2])*hamilR[13]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*((3.872983346207417*b_z[0]*jacobTotInv[0]*hamilR[19]-3.0*b_z[0]*jacobTotInv[0]*hamilR[5]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilR[2])*m_*rdy2R+((3.0*BstarXdBmagR[1]-1.732050807568877*BstarXdBmagR[0])*hamilR[3]-3.872983346207417*hamilR[3]*BstarXdBmagR[11])*q_*rdvpar2R))/(m_*q_); 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[11]+fL[11])+1.732050807568877*(fL[1]-1.0*fR[1])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[19]+fL[19])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[2]+fL[2])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[21]+fL[21])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[25]+fL[25])+3.0*(fL[8]-1.0*fR[8]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[32]+fL[32])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[7]+fL[7]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[35]+fL[35])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[9]+fL[9]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[37]+fL[37])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[20]-1.0*(8.660254037844387*fL[20]+5.0*(fR[12]+fL[12]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[23]-1.0*(8.660254037844387*fL[23]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[28]-1.0*(8.660254037844387*fL[28]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[44]+fL[44])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[18]+fL[18])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[22]+fL[22]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[24]+fL[24]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[36]-1.0*(8.660254037844387*fL[36]+5.0*(fR[26]+fL[26]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[39]-1.0*(8.660254037844387*fL[39]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[29]+fL[29]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[38]+fL[38]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[40]+fL[40]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[43]+fL[43]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[11]-6.324555320336761*fL[11]-4.898979485566357*(fR[1]+fL[1])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[19]-94.8683298050514*fL[19]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*(alphaR[4]*fAvg[11]+alphaR[1]*fAvg[7])-21.21320343559643*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[21]-94.8683298050514*fL[21]-73.48469228349535*(fR[6]+fL[6])+42.42640687119286*fR[3]-42.42640687119286*fL[3])*amax-18.97366596101028*(alphaR[4]*fAvg[12]+alphaR[2]*fAvg[8])-21.21320343559643*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[25]-18.97366596101028*fL[25]-14.69693845669907*(fR[8]+fL[8])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6]+alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[32]-94.86832980505142*fL[32]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[7]-42.42640687119286*fL[7])*amax-18.97366596101028*(alphaR[2]*fAvg[12]+alphaR[1]*fAvg[11]+alphaR[4]*(fAvg[8]+fAvg[7]))-21.21320343559643*(alphaR[0]*fAvg[4]+fAvg[0]*alphaR[4]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[35]-94.86832980505142*fL[35]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*(alphaR[4]*fAvg[17]+alphaR[1]*fAvg[13])-21.21320343559643*(alphaR[2]*fAvg[10]+alphaR[4]*fAvg[6]+alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[37]-94.86832980505142*fL[37]-73.48469228349535*(fR[17]+fL[17])+42.42640687119286*fR[10]-42.42640687119286*fL[10])*amax-18.97366596101028*(alphaR[4]*fAvg[18]+alphaR[2]*fAvg[14])-21.21320343559643*(alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6]+alphaR[4]*fAvg[5]+alphaR[2]*fAvg[3])); 
  Ghat[7] = 0.008333333333333333*((73.48469228349536*(fR[20]+fL[20])-42.42640687119286*fR[12]+42.42640687119286*fL[12])*amax+21.21320343559643*alphaR[2]*fAvg[11]+21.21320343559643*alphaR[0]*fAvg[7]+18.97366596101028*(alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1])); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[23]+fL[23])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaR[1]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[8]+18.97366596101028*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[28]+fL[28])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[4]*fAvg[19]+21.21320343559643*(alphaR[2]*fAvg[16]+alphaR[1]*fAvg[15])+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[44]-94.8683298050514*fL[44]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[18]-42.42640687119286*fL[18])*amax-18.97366596101028*(alphaR[2]*fAvg[18]+alphaR[1]*fAvg[17]+alphaR[4]*(fAvg[14]+fAvg[13]))-21.21320343559643*(alphaR[0]*fAvg[10]+alphaR[1]*fAvg[6]+alphaR[2]*fAvg[5]+fAvg[3]*alphaR[4])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[33]+fL[33])-42.42640687119286*fR[22]+42.42640687119286*fL[22])*amax+16.97056274847715*alphaR[4]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[11]+21.21320343559643*alphaR[2]*fAvg[7]+18.97366596101028*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[24]+42.42640687119286*fL[24])*amax+21.21320343559643*alphaR[0]*fAvg[12]+16.97056274847715*alphaR[4]*fAvg[11]+21.21320343559643*alphaR[1]*fAvg[8]+18.97366596101028*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[36]+fL[36])-42.42640687119286*fR[26]+42.42640687119286*fL[26])*amax+21.21320343559643*alphaR[2]*fAvg[17]+21.21320343559643*alphaR[0]*fAvg[13]+18.97366596101028*(alphaR[4]*fAvg[10]+alphaR[1]*fAvg[5])); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[39]+fL[39])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaR[1]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[14]+18.97366596101028*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6])); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+21.21320343559643*alphaR[2]*fAvg[19]+21.21320343559643*(alphaR[4]*fAvg[16]+alphaR[0]*fAvg[15])+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*(alphaR[0]*fAvg[16]+alphaR[4]*fAvg[15])+21.21320343559643*alphaR[2]*fAvg[9]); 
  Ghat[17] = 0.008333333333333333*((73.48469228349536*(fR[45]+fL[45])-42.42640687119286*fR[38]+42.42640687119286*fL[38])*amax+16.97056274847715*alphaR[4]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[17]+21.21320343559643*alphaR[2]*fAvg[13]+18.97366596101028*(alphaR[1]*fAvg[10]+alphaR[4]*fAvg[5])); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[40]+42.42640687119286*fL[40])*amax+21.21320343559643*alphaR[0]*fAvg[18]+16.97056274847715*alphaR[4]*fAvg[17]+21.21320343559643*alphaR[1]*fAvg[14]+18.97366596101028*(alphaR[2]*fAvg[10]+alphaR[4]*fAvg[6])); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[43]+42.42640687119286*fL[43])*amax+21.21320343559643*alphaR[0]*fAvg[19]+21.21320343559643*(alphaR[1]*fAvg[16]+alphaR[2]*fAvg[15])+21.21320343559643*alphaR[4]*fAvg[9]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = -1.224744871391589*Ghat[0]; 
  incr[2] = 0.7071067811865475*Ghat[1]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = -1.224744871391589*Ghat[2]; 
  incr[7] = 0.7071067811865475*Ghat[4]; 
  incr[8] = -1.224744871391589*Ghat[3]; 
  incr[9] = 0.7071067811865475*Ghat[5]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 1.58113883008419*Ghat[0]; 
  incr[12] = 0.7071067811865475*Ghat[7]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = -1.224744871391589*Ghat[6]; 
  incr[18] = 0.7071067811865475*Ghat[10]; 
  incr[19] = 1.58113883008419*Ghat[1]; 
  incr[20] = -1.224744871391589*Ghat[7]; 
  incr[21] = 1.58113883008419*Ghat[2]; 
  incr[22] = 0.7071067811865475*Ghat[11]; 
  incr[23] = -1.224744871391589*Ghat[8]; 
  incr[24] = 0.7071067811865475*Ghat[12]; 
  incr[25] = 1.58113883008419*Ghat[3]; 
  incr[26] = 0.7071067811865475*Ghat[13]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = -1.224744871391589*Ghat[9]; 
  incr[29] = 0.7071067811865475*Ghat[15]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = 1.58113883008419*Ghat[4]; 
  incr[33] = -1.224744871391589*Ghat[11]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = 1.58113883008419*Ghat[5]; 
  incr[36] = -1.224744871391589*Ghat[13]; 
  incr[37] = 1.58113883008419*Ghat[6]; 
  incr[38] = 0.7071067811865475*Ghat[17]; 
  incr[39] = -1.224744871391589*Ghat[14]; 
  incr[40] = 0.7071067811865475*Ghat[18]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = -1.224744871391589*Ghat[16]; 
  incr[43] = 0.7071067811865475*Ghat[19]; 
  incr[44] = 1.58113883008419*Ghat[10]; 
  incr[45] = -1.224744871391589*Ghat[17]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

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
  outR[16] += incr[16]*rdx2R; 
  outR[17] += incr[17]*rdx2R; 
  outR[18] += incr[18]*rdx2R; 
  outR[19] += incr[19]*rdx2R; 
  outR[20] += incr[20]*rdx2R; 
  outR[21] += incr[21]*rdx2R; 
  outR[22] += incr[22]*rdx2R; 
  outR[23] += incr[23]*rdx2R; 
  outR[24] += incr[24]*rdx2R; 
  outR[25] += incr[25]*rdx2R; 
  outR[26] += incr[26]*rdx2R; 
  outR[27] += incr[27]*rdx2R; 
  outR[28] += incr[28]*rdx2R; 
  outR[29] += incr[29]*rdx2R; 
  outR[30] += incr[30]*rdx2R; 
  outR[31] += incr[31]*rdx2R; 
  outR[32] += incr[32]*rdx2R; 
  outR[33] += incr[33]*rdx2R; 
  outR[34] += incr[34]*rdx2R; 
  outR[35] += incr[35]*rdx2R; 
  outR[36] += incr[36]*rdx2R; 
  outR[37] += incr[37]*rdx2R; 
  outR[38] += incr[38]*rdx2R; 
  outR[39] += incr[39]*rdx2R; 
  outR[40] += incr[40]*rdx2R; 
  outR[41] += incr[41]*rdx2R; 
  outR[42] += incr[42]*rdx2R; 
  outR[43] += incr[43]*rdx2R; 
  outR[44] += incr[44]*rdx2R; 
  outR[45] += incr[45]*rdx2R; 
  outR[46] += incr[46]*rdx2R; 
  outR[47] += incr[47]*rdx2R; 

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
  outL[11] += -1.0*incr[11]*rdx2L; 
  outL[12] += -1.0*incr[12]*rdx2L; 
  outL[13] += -1.0*incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 
  outL[16] += incr[16]*rdx2L; 
  outL[17] += incr[17]*rdx2L; 
  outL[18] += -1.0*incr[18]*rdx2L; 
  outL[19] += -1.0*incr[19]*rdx2L; 
  outL[20] += incr[20]*rdx2L; 
  outL[21] += -1.0*incr[21]*rdx2L; 
  outL[22] += -1.0*incr[22]*rdx2L; 
  outL[23] += incr[23]*rdx2L; 
  outL[24] += -1.0*incr[24]*rdx2L; 
  outL[25] += -1.0*incr[25]*rdx2L; 
  outL[26] += -1.0*incr[26]*rdx2L; 
  outL[27] += -1.0*incr[27]*rdx2L; 
  outL[28] += incr[28]*rdx2L; 
  outL[29] += -1.0*incr[29]*rdx2L; 
  outL[30] += -1.0*incr[30]*rdx2L; 
  outL[31] += incr[31]*rdx2L; 
  outL[32] += -1.0*incr[32]*rdx2L; 
  outL[33] += incr[33]*rdx2L; 
  outL[34] += incr[34]*rdx2L; 
  outL[35] += -1.0*incr[35]*rdx2L; 
  outL[36] += incr[36]*rdx2L; 
  outL[37] += -1.0*incr[37]*rdx2L; 
  outL[38] += -1.0*incr[38]*rdx2L; 
  outL[39] += incr[39]*rdx2L; 
  outL[40] += -1.0*incr[40]*rdx2L; 
  outL[41] += incr[41]*rdx2L; 
  outL[42] += incr[42]*rdx2L; 
  outL[43] += -1.0*incr[43]*rdx2L; 
  outL[44] += -1.0*incr[44]*rdx2L; 
  outL[45] += incr[45]*rdx2L; 
  outL[46] += incr[46]*rdx2L; 
  outL[47] += incr[47]*rdx2L; 

  return std::abs(alphaSurfAvgR);
} 
double EmGyrokineticGenGeoSurf2x2vSer_yL_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilL[48]; 
  hamilL[0] = (0.6666666666666666*(3.0*rdvpar2SqL*(m_*wvparSqL+bmag[0]*wmuL+phi[0]*q_)+m_))/rdvpar2SqL; 
  hamilL[1] = 2.0*phi[1]*q_; 
  hamilL[2] = 2.0*phi[2]*q_; 
  hamilL[3] = (2.309401076758503*m_*wvparL)/rdvpar2L; 
  hamilL[4] = (1.154700538379252*bmag[0])/rdmu2L; 
  hamilL[5] = 2.0*phi[3]*q_; 
  hamilL[11] = 2.0*phi[4]*q_; 
  hamilL[12] = 2.0*phi[5]*q_; 
  hamilL[13] = (0.5962847939999438*m_)/rdvpar2SqL; 
  hamilL[19] = 2.0*phi[6]*q_; 
  hamilL[20] = 2.0*phi[7]*q_; 

  double BstarYdBmagL[48]; 
  BstarYdBmagL[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2L; 
  BstarYdBmagL[1] = -1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[4]*rdx2L; 
  BstarYdBmagL[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2L; 
  BstarYdBmagL[5] = -1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[6]*rdx2L; 
  BstarYdBmagL[12] = -0.8660254037844387*b_z[0]*jacobTotInv[0]*Apar[7]*rdx2L; 

  double alphaL[20]; 
  alphaL[0] = (0.1767766952966368*(b_z[0]*jacobTotInv[0]*(3.872983346207417*hamilL[20]+3.0*hamilL[5]+1.732050807568877*hamilL[1])*m_*rdx2L+hamilL[3]*(3.872983346207417*BstarYdBmagL[12]+3.0*BstarYdBmagL[2]+1.732050807568877*BstarYdBmagL[0])*q_*rdvpar2L))/(m_*q_); 
  alphaL[1] = (0.1767766952966368*(b_z[0]*jacobTotInv[0]*(6.708203932499369*hamilL[19]+3.872983346207417*hamilL[11])*m_*rdx2L+hamilL[3]*(3.0*BstarYdBmagL[5]+1.732050807568877*BstarYdBmagL[1])*q_*rdvpar2L))/(m_*q_); 
  alphaL[2] = (0.1767766952966368*(8.660254037844386*BstarYdBmagL[12]+6.708203932499369*BstarYdBmagL[2]+3.872983346207417*BstarYdBmagL[0])*hamilL[13]*rdvpar2L)/m_; 
  alphaL[4] = (0.1767766952966368*(6.708203932499369*BstarYdBmagL[5]+3.872983346207417*BstarYdBmagL[1])*hamilL[13]*rdvpar2L)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgL = (0.03125*((3.872983346207417*b_z[0]*jacobTotInv[0]*hamilL[20]+3.0*b_z[0]*jacobTotInv[0]*hamilL[5]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilL[1])*m_*rdx2L+(3.872983346207417*hamilL[3]*BstarYdBmagL[12]+(3.0*BstarYdBmagL[2]+1.732050807568877*BstarYdBmagL[0])*hamilL[3])*q_*rdvpar2L))/(m_*q_); 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[12]+fL[12])+1.732050807568877*(fL[2]-1.0*fR[2])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[20]+fL[20])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[22]+fL[22])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[26]+fL[26])+3.0*(fL[9]-1.0*fR[9]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[33]+fL[33])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[6]+fL[6]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[36]+fL[36])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[38]+fL[38])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[19]-1.0*(8.660254037844387*fL[19]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[24]-1.0*(8.660254037844387*fL[24]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[29]-1.0*(8.660254037844387*fL[29]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[45]+fL[45])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[17]+fL[17])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[21]+fL[21]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[23]+fL[23]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[35]-1.0*(8.660254037844387*fL[35]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[40]-1.0*(8.660254037844387*fL[40]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[37]+fL[37]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[39]+fL[39]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[42]+fL[42]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[12]-6.324555320336761*fL[12]-4.898979485566357*(fR[2]+fL[2])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaL[4]*fAvg[4]+alphaL[2]*fAvg[2]+alphaL[1]*fAvg[1]+alphaL[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[20]-94.8683298050514*fL[20]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*(alphaL[4]*fAvg[11]+alphaL[1]*fAvg[7])-21.21320343559643*(alphaL[2]*fAvg[4]+fAvg[2]*alphaL[4]+alphaL[0]*fAvg[1]+fAvg[0]*alphaL[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[22]-94.8683298050514*fL[22]-73.48469228349535*(fR[7]+fL[7])+42.42640687119286*fR[3]-42.42640687119286*fL[3])*amax-18.97366596101028*(alphaL[4]*fAvg[12]+alphaL[2]*fAvg[8])-21.21320343559643*(alphaL[1]*fAvg[4]+fAvg[1]*alphaL[4]+alphaL[0]*fAvg[2]+fAvg[0]*alphaL[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[26]-18.97366596101028*fL[26]-14.69693845669907*(fR[9]+fL[9])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaL[4]*fAvg[10]+alphaL[2]*fAvg[6]+alphaL[1]*fAvg[5]+alphaL[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[33]-94.86832980505142*fL[33]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[6]-42.42640687119286*fL[6])*amax-18.97366596101028*(alphaL[2]*fAvg[12]+alphaL[1]*fAvg[11]+alphaL[4]*(fAvg[8]+fAvg[7]))-21.21320343559643*(alphaL[0]*fAvg[4]+fAvg[0]*alphaL[4]+alphaL[1]*fAvg[2]+fAvg[1]*alphaL[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[36]-94.86832980505142*fL[36]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-18.97366596101028*(alphaL[4]*fAvg[17]+alphaL[1]*fAvg[13])-21.21320343559643*(alphaL[2]*fAvg[10]+alphaL[4]*fAvg[6]+alphaL[0]*fAvg[5]+alphaL[1]*fAvg[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[38]-94.86832980505142*fL[38]-73.48469228349535*(fR[18]+fL[18])+42.42640687119286*fR[10]-42.42640687119286*fL[10])*amax-18.97366596101028*(alphaL[4]*fAvg[18]+alphaL[2]*fAvg[14])-21.21320343559643*(alphaL[1]*fAvg[10]+alphaL[0]*fAvg[6]+alphaL[4]*fAvg[5]+alphaL[2]*fAvg[3])); 
  Ghat[7] = 0.008333333333333333*((73.48469228349536*(fR[19]+fL[19])-42.42640687119286*fR[11]+42.42640687119286*fL[11])*amax+21.21320343559643*alphaL[2]*fAvg[11]+21.21320343559643*alphaL[0]*fAvg[7]+18.97366596101028*(alphaL[4]*fAvg[4]+alphaL[1]*fAvg[1])); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[24]+fL[24])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaL[1]*fAvg[12]+21.21320343559643*alphaL[0]*fAvg[8]+18.97366596101028*(alphaL[4]*fAvg[4]+alphaL[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[29]+fL[29])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaL[4]*fAvg[19]+21.21320343559643*(alphaL[2]*fAvg[16]+alphaL[1]*fAvg[15])+21.21320343559643*alphaL[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[45]-94.8683298050514*fL[45]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[17]-42.42640687119286*fL[17])*amax-18.97366596101028*(alphaL[2]*fAvg[18]+alphaL[1]*fAvg[17]+alphaL[4]*(fAvg[14]+fAvg[13]))-21.21320343559643*(alphaL[0]*fAvg[10]+alphaL[1]*fAvg[6]+alphaL[2]*fAvg[5]+fAvg[3]*alphaL[4])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[32]+fL[32])-42.42640687119286*fR[21]+42.42640687119286*fL[21])*amax+16.97056274847715*alphaL[4]*fAvg[12]+21.21320343559643*alphaL[0]*fAvg[11]+21.21320343559643*alphaL[2]*fAvg[7]+18.97366596101028*(alphaL[1]*fAvg[4]+fAvg[1]*alphaL[4])); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[23]+42.42640687119286*fL[23])*amax+21.21320343559643*alphaL[0]*fAvg[12]+16.97056274847715*alphaL[4]*fAvg[11]+21.21320343559643*alphaL[1]*fAvg[8]+18.97366596101028*(alphaL[2]*fAvg[4]+fAvg[2]*alphaL[4])); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[35]+fL[35])-42.42640687119286*fR[25]+42.42640687119286*fL[25])*amax+21.21320343559643*alphaL[2]*fAvg[17]+21.21320343559643*alphaL[0]*fAvg[13]+18.97366596101028*(alphaL[4]*fAvg[10]+alphaL[1]*fAvg[5])); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[40]+fL[40])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaL[1]*fAvg[18]+21.21320343559643*alphaL[0]*fAvg[14]+18.97366596101028*(alphaL[4]*fAvg[10]+alphaL[2]*fAvg[6])); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+21.21320343559643*alphaL[2]*fAvg[19]+21.21320343559643*(alphaL[4]*fAvg[16]+alphaL[0]*fAvg[15])+21.21320343559643*alphaL[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaL[1]*fAvg[19]+21.21320343559643*(alphaL[0]*fAvg[16]+alphaL[4]*fAvg[15])+21.21320343559643*alphaL[2]*fAvg[9]); 
  Ghat[17] = 0.008333333333333333*((73.48469228349536*(fR[44]+fL[44])-42.42640687119286*fR[37]+42.42640687119286*fL[37])*amax+16.97056274847715*alphaL[4]*fAvg[18]+21.21320343559643*alphaL[0]*fAvg[17]+21.21320343559643*alphaL[2]*fAvg[13]+18.97366596101028*(alphaL[1]*fAvg[10]+alphaL[4]*fAvg[5])); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[39]+42.42640687119286*fL[39])*amax+21.21320343559643*alphaL[0]*fAvg[18]+16.97056274847715*alphaL[4]*fAvg[17]+21.21320343559643*alphaL[1]*fAvg[14]+18.97366596101028*(alphaL[2]*fAvg[10]+alphaL[4]*fAvg[6])); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[42]+42.42640687119286*fL[42])*amax+21.21320343559643*alphaL[0]*fAvg[19]+21.21320343559643*(alphaL[1]*fAvg[16]+alphaL[2]*fAvg[15])+21.21320343559643*alphaL[4]*fAvg[9]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = -1.224744871391589*Ghat[0]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = 0.7071067811865475*Ghat[4]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = -1.224744871391589*Ghat[3]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 1.58113883008419*Ghat[0]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = 0.7071067811865475*Ghat[10]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = -1.224744871391589*Ghat[7]; 
  incr[20] = 1.58113883008419*Ghat[1]; 
  incr[21] = 0.7071067811865475*Ghat[11]; 
  incr[22] = 1.58113883008419*Ghat[2]; 
  incr[23] = 0.7071067811865475*Ghat[12]; 
  incr[24] = -1.224744871391589*Ghat[8]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 1.58113883008419*Ghat[3]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = -1.224744871391589*Ghat[9]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = 1.58113883008419*Ghat[4]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = -1.224744871391589*Ghat[13]; 
  incr[36] = 1.58113883008419*Ghat[5]; 
  incr[37] = 0.7071067811865475*Ghat[17]; 
  incr[38] = 1.58113883008419*Ghat[6]; 
  incr[39] = 0.7071067811865475*Ghat[18]; 
  incr[40] = -1.224744871391589*Ghat[14]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = 0.7071067811865475*Ghat[19]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = 1.58113883008419*Ghat[10]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

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
  outR[16] += incr[16]*rdy2R; 
  outR[17] += incr[17]*rdy2R; 
  outR[18] += incr[18]*rdy2R; 
  outR[19] += incr[19]*rdy2R; 
  outR[20] += incr[20]*rdy2R; 
  outR[21] += incr[21]*rdy2R; 
  outR[22] += incr[22]*rdy2R; 
  outR[23] += incr[23]*rdy2R; 
  outR[24] += incr[24]*rdy2R; 
  outR[25] += incr[25]*rdy2R; 
  outR[26] += incr[26]*rdy2R; 
  outR[27] += incr[27]*rdy2R; 
  outR[28] += incr[28]*rdy2R; 
  outR[29] += incr[29]*rdy2R; 
  outR[30] += incr[30]*rdy2R; 
  outR[31] += incr[31]*rdy2R; 
  outR[32] += incr[32]*rdy2R; 
  outR[33] += incr[33]*rdy2R; 
  outR[34] += incr[34]*rdy2R; 
  outR[35] += incr[35]*rdy2R; 
  outR[36] += incr[36]*rdy2R; 
  outR[37] += incr[37]*rdy2R; 
  outR[38] += incr[38]*rdy2R; 
  outR[39] += incr[39]*rdy2R; 
  outR[40] += incr[40]*rdy2R; 
  outR[41] += incr[41]*rdy2R; 
  outR[42] += incr[42]*rdy2R; 
  outR[43] += incr[43]*rdy2R; 
  outR[44] += incr[44]*rdy2R; 
  outR[45] += incr[45]*rdy2R; 
  outR[46] += incr[46]*rdy2R; 
  outR[47] += incr[47]*rdy2R; 

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
  outL[11] += -1.0*incr[11]*rdy2L; 
  outL[12] += -1.0*incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += -1.0*incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 
  outL[16] += incr[16]*rdy2L; 
  outL[17] += -1.0*incr[17]*rdy2L; 
  outL[18] += incr[18]*rdy2L; 
  outL[19] += incr[19]*rdy2L; 
  outL[20] += -1.0*incr[20]*rdy2L; 
  outL[21] += -1.0*incr[21]*rdy2L; 
  outL[22] += -1.0*incr[22]*rdy2L; 
  outL[23] += -1.0*incr[23]*rdy2L; 
  outL[24] += incr[24]*rdy2L; 
  outL[25] += -1.0*incr[25]*rdy2L; 
  outL[26] += -1.0*incr[26]*rdy2L; 
  outL[27] += -1.0*incr[27]*rdy2L; 
  outL[28] += -1.0*incr[28]*rdy2L; 
  outL[29] += incr[29]*rdy2L; 
  outL[30] += -1.0*incr[30]*rdy2L; 
  outL[31] += incr[31]*rdy2L; 
  outL[32] += incr[32]*rdy2L; 
  outL[33] += -1.0*incr[33]*rdy2L; 
  outL[34] += incr[34]*rdy2L; 
  outL[35] += incr[35]*rdy2L; 
  outL[36] += -1.0*incr[36]*rdy2L; 
  outL[37] += -1.0*incr[37]*rdy2L; 
  outL[38] += -1.0*incr[38]*rdy2L; 
  outL[39] += -1.0*incr[39]*rdy2L; 
  outL[40] += incr[40]*rdy2L; 
  outL[41] += incr[41]*rdy2L; 
  outL[42] += -1.0*incr[42]*rdy2L; 
  outL[43] += incr[43]*rdy2L; 
  outL[44] += incr[44]*rdy2L; 
  outL[45] += -1.0*incr[45]*rdy2L; 
  outL[46] += incr[46]*rdy2L; 
  outL[47] += incr[47]*rdy2L; 

  return std::abs(alphaSurfAvgL);
} 
double EmGyrokineticGenGeoSurf2x2vSer_yR_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[11] = 2.0*phi[4]*q_; 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 

  double BstarYdBmagR[48]; 
  BstarYdBmagR[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2R; 
  BstarYdBmagR[1] = -1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[4]*rdx2R; 
  BstarYdBmagR[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2R; 
  BstarYdBmagR[5] = -1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[6]*rdx2R; 
  BstarYdBmagR[12] = -0.8660254037844387*b_z[0]*jacobTotInv[0]*Apar[7]*rdx2R; 

  double alphaR[20]; 
  alphaR[0] = (0.1767766952966368*(b_z[0]*jacobTotInv[0]*(3.872983346207417*hamilR[20]-3.0*hamilR[5]+1.732050807568877*hamilR[1])*m_*rdx2R+hamilR[3]*(3.872983346207417*BstarYdBmagR[12]-3.0*BstarYdBmagR[2]+1.732050807568877*BstarYdBmagR[0])*q_*rdvpar2R))/(m_*q_); 
  alphaR[1] = -(0.1767766952966368*(b_z[0]*jacobTotInv[0]*(6.708203932499369*hamilR[19]-3.872983346207417*hamilR[11])*m_*rdx2R+hamilR[3]*(3.0*BstarYdBmagR[5]-1.732050807568877*BstarYdBmagR[1])*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = (0.1767766952966368*(8.660254037844386*BstarYdBmagR[12]-6.708203932499369*BstarYdBmagR[2]+3.872983346207417*BstarYdBmagR[0])*hamilR[13]*rdvpar2R)/m_; 
  alphaR[4] = -(0.1767766952966368*(6.708203932499369*BstarYdBmagR[5]-3.872983346207417*BstarYdBmagR[1])*hamilR[13]*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.03125*((3.872983346207417*b_z[0]*jacobTotInv[0]*hamilR[20]-3.0*b_z[0]*jacobTotInv[0]*hamilR[5]+1.732050807568877*b_z[0]*jacobTotInv[0]*hamilR[1])*m_*rdx2R+(3.872983346207417*hamilR[3]*BstarYdBmagR[12]+(1.732050807568877*BstarYdBmagR[0]-3.0*BstarYdBmagR[2])*hamilR[3])*q_*rdvpar2R))/(m_*q_); 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[12]+fL[12])+1.732050807568877*(fL[2]-1.0*fR[2])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[20]+fL[20])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[22]+fL[22])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[26]+fL[26])+3.0*(fL[9]-1.0*fR[9]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[33]+fL[33])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[6]+fL[6]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[36]+fL[36])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[38]+fL[38])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[19]-1.0*(8.660254037844387*fL[19]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[24]-1.0*(8.660254037844387*fL[24]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[29]-1.0*(8.660254037844387*fL[29]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[45]+fL[45])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[17]+fL[17])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[21]+fL[21]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[23]+fL[23]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[35]-1.0*(8.660254037844387*fL[35]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[40]-1.0*(8.660254037844387*fL[40]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[37]+fL[37]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[39]+fL[39]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[42]+fL[42]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[12]-6.324555320336761*fL[12]-4.898979485566357*(fR[2]+fL[2])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[20]-94.8683298050514*fL[20]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*(alphaR[4]*fAvg[11]+alphaR[1]*fAvg[7])-21.21320343559643*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[22]-94.8683298050514*fL[22]-73.48469228349535*(fR[7]+fL[7])+42.42640687119286*fR[3]-42.42640687119286*fL[3])*amax-18.97366596101028*(alphaR[4]*fAvg[12]+alphaR[2]*fAvg[8])-21.21320343559643*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[26]-18.97366596101028*fL[26]-14.69693845669907*(fR[9]+fL[9])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6]+alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[33]-94.86832980505142*fL[33]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[6]-42.42640687119286*fL[6])*amax-18.97366596101028*(alphaR[2]*fAvg[12]+alphaR[1]*fAvg[11]+alphaR[4]*(fAvg[8]+fAvg[7]))-21.21320343559643*(alphaR[0]*fAvg[4]+fAvg[0]*alphaR[4]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[36]-94.86832980505142*fL[36]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-18.97366596101028*(alphaR[4]*fAvg[17]+alphaR[1]*fAvg[13])-21.21320343559643*(alphaR[2]*fAvg[10]+alphaR[4]*fAvg[6]+alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[38]-94.86832980505142*fL[38]-73.48469228349535*(fR[18]+fL[18])+42.42640687119286*fR[10]-42.42640687119286*fL[10])*amax-18.97366596101028*(alphaR[4]*fAvg[18]+alphaR[2]*fAvg[14])-21.21320343559643*(alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6]+alphaR[4]*fAvg[5]+alphaR[2]*fAvg[3])); 
  Ghat[7] = 0.008333333333333333*((73.48469228349536*(fR[19]+fL[19])-42.42640687119286*fR[11]+42.42640687119286*fL[11])*amax+21.21320343559643*alphaR[2]*fAvg[11]+21.21320343559643*alphaR[0]*fAvg[7]+18.97366596101028*(alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1])); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[24]+fL[24])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaR[1]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[8]+18.97366596101028*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[29]+fL[29])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[4]*fAvg[19]+21.21320343559643*(alphaR[2]*fAvg[16]+alphaR[1]*fAvg[15])+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[45]-94.8683298050514*fL[45]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[17]-42.42640687119286*fL[17])*amax-18.97366596101028*(alphaR[2]*fAvg[18]+alphaR[1]*fAvg[17]+alphaR[4]*(fAvg[14]+fAvg[13]))-21.21320343559643*(alphaR[0]*fAvg[10]+alphaR[1]*fAvg[6]+alphaR[2]*fAvg[5]+fAvg[3]*alphaR[4])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[32]+fL[32])-42.42640687119286*fR[21]+42.42640687119286*fL[21])*amax+16.97056274847715*alphaR[4]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[11]+21.21320343559643*alphaR[2]*fAvg[7]+18.97366596101028*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[23]+42.42640687119286*fL[23])*amax+21.21320343559643*alphaR[0]*fAvg[12]+16.97056274847715*alphaR[4]*fAvg[11]+21.21320343559643*alphaR[1]*fAvg[8]+18.97366596101028*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[35]+fL[35])-42.42640687119286*fR[25]+42.42640687119286*fL[25])*amax+21.21320343559643*alphaR[2]*fAvg[17]+21.21320343559643*alphaR[0]*fAvg[13]+18.97366596101028*(alphaR[4]*fAvg[10]+alphaR[1]*fAvg[5])); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[40]+fL[40])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaR[1]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[14]+18.97366596101028*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6])); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+21.21320343559643*alphaR[2]*fAvg[19]+21.21320343559643*(alphaR[4]*fAvg[16]+alphaR[0]*fAvg[15])+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*(alphaR[0]*fAvg[16]+alphaR[4]*fAvg[15])+21.21320343559643*alphaR[2]*fAvg[9]); 
  Ghat[17] = 0.008333333333333333*((73.48469228349536*(fR[44]+fL[44])-42.42640687119286*fR[37]+42.42640687119286*fL[37])*amax+16.97056274847715*alphaR[4]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[17]+21.21320343559643*alphaR[2]*fAvg[13]+18.97366596101028*(alphaR[1]*fAvg[10]+alphaR[4]*fAvg[5])); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[39]+42.42640687119286*fL[39])*amax+21.21320343559643*alphaR[0]*fAvg[18]+16.97056274847715*alphaR[4]*fAvg[17]+21.21320343559643*alphaR[1]*fAvg[14]+18.97366596101028*(alphaR[2]*fAvg[10]+alphaR[4]*fAvg[6])); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[42]+42.42640687119286*fL[42])*amax+21.21320343559643*alphaR[0]*fAvg[19]+21.21320343559643*(alphaR[1]*fAvg[16]+alphaR[2]*fAvg[15])+21.21320343559643*alphaR[4]*fAvg[9]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = -1.224744871391589*Ghat[0]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = 0.7071067811865475*Ghat[4]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = -1.224744871391589*Ghat[3]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 1.58113883008419*Ghat[0]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = 0.7071067811865475*Ghat[10]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = -1.224744871391589*Ghat[7]; 
  incr[20] = 1.58113883008419*Ghat[1]; 
  incr[21] = 0.7071067811865475*Ghat[11]; 
  incr[22] = 1.58113883008419*Ghat[2]; 
  incr[23] = 0.7071067811865475*Ghat[12]; 
  incr[24] = -1.224744871391589*Ghat[8]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 1.58113883008419*Ghat[3]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = -1.224744871391589*Ghat[9]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = 1.58113883008419*Ghat[4]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = -1.224744871391589*Ghat[13]; 
  incr[36] = 1.58113883008419*Ghat[5]; 
  incr[37] = 0.7071067811865475*Ghat[17]; 
  incr[38] = 1.58113883008419*Ghat[6]; 
  incr[39] = 0.7071067811865475*Ghat[18]; 
  incr[40] = -1.224744871391589*Ghat[14]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = 0.7071067811865475*Ghat[19]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = 1.58113883008419*Ghat[10]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

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
  outR[16] += incr[16]*rdy2R; 
  outR[17] += incr[17]*rdy2R; 
  outR[18] += incr[18]*rdy2R; 
  outR[19] += incr[19]*rdy2R; 
  outR[20] += incr[20]*rdy2R; 
  outR[21] += incr[21]*rdy2R; 
  outR[22] += incr[22]*rdy2R; 
  outR[23] += incr[23]*rdy2R; 
  outR[24] += incr[24]*rdy2R; 
  outR[25] += incr[25]*rdy2R; 
  outR[26] += incr[26]*rdy2R; 
  outR[27] += incr[27]*rdy2R; 
  outR[28] += incr[28]*rdy2R; 
  outR[29] += incr[29]*rdy2R; 
  outR[30] += incr[30]*rdy2R; 
  outR[31] += incr[31]*rdy2R; 
  outR[32] += incr[32]*rdy2R; 
  outR[33] += incr[33]*rdy2R; 
  outR[34] += incr[34]*rdy2R; 
  outR[35] += incr[35]*rdy2R; 
  outR[36] += incr[36]*rdy2R; 
  outR[37] += incr[37]*rdy2R; 
  outR[38] += incr[38]*rdy2R; 
  outR[39] += incr[39]*rdy2R; 
  outR[40] += incr[40]*rdy2R; 
  outR[41] += incr[41]*rdy2R; 
  outR[42] += incr[42]*rdy2R; 
  outR[43] += incr[43]*rdy2R; 
  outR[44] += incr[44]*rdy2R; 
  outR[45] += incr[45]*rdy2R; 
  outR[46] += incr[46]*rdy2R; 
  outR[47] += incr[47]*rdy2R; 

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
  outL[11] += -1.0*incr[11]*rdy2L; 
  outL[12] += -1.0*incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += -1.0*incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 
  outL[16] += incr[16]*rdy2L; 
  outL[17] += -1.0*incr[17]*rdy2L; 
  outL[18] += incr[18]*rdy2L; 
  outL[19] += incr[19]*rdy2L; 
  outL[20] += -1.0*incr[20]*rdy2L; 
  outL[21] += -1.0*incr[21]*rdy2L; 
  outL[22] += -1.0*incr[22]*rdy2L; 
  outL[23] += -1.0*incr[23]*rdy2L; 
  outL[24] += incr[24]*rdy2L; 
  outL[25] += -1.0*incr[25]*rdy2L; 
  outL[26] += -1.0*incr[26]*rdy2L; 
  outL[27] += -1.0*incr[27]*rdy2L; 
  outL[28] += -1.0*incr[28]*rdy2L; 
  outL[29] += incr[29]*rdy2L; 
  outL[30] += -1.0*incr[30]*rdy2L; 
  outL[31] += incr[31]*rdy2L; 
  outL[32] += incr[32]*rdy2L; 
  outL[33] += -1.0*incr[33]*rdy2L; 
  outL[34] += incr[34]*rdy2L; 
  outL[35] += incr[35]*rdy2L; 
  outL[36] += -1.0*incr[36]*rdy2L; 
  outL[37] += -1.0*incr[37]*rdy2L; 
  outL[38] += -1.0*incr[38]*rdy2L; 
  outL[39] += -1.0*incr[39]*rdy2L; 
  outL[40] += incr[40]*rdy2L; 
  outL[41] += incr[41]*rdy2L; 
  outL[42] += -1.0*incr[42]*rdy2L; 
  outL[43] += incr[43]*rdy2L; 
  outL[44] += incr[44]*rdy2L; 
  outL[45] += -1.0*incr[45]*rdy2L; 
  outL[46] += incr[46]*rdy2L; 
  outL[47] += incr[47]*rdy2L; 

  return std::abs(alphaSurfAvgR);
} 
double EmGyrokineticGenGeoSurf2x2vSer_vpar_P2_Bvars(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*phi[1]*q_; 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[11] = 2.0*phi[4]*q_; 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 

  double BstarXdBmagR[48]; 
  BstarXdBmagR[0] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[2]*rdy2R; 
  BstarXdBmagR[1] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdy2R; 
  BstarXdBmagR[2] = 1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[5]*rdy2R; 
  BstarXdBmagR[5] = 1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[7]*rdy2R; 
  BstarXdBmagR[11] = 0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[6]*rdy2R; 

  double BstarYdBmagR[48]; 
  BstarYdBmagR[0] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[1]*rdx2R; 
  BstarYdBmagR[1] = -1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[4]*rdx2R; 
  BstarYdBmagR[2] = -0.8660254037844386*b_z[0]*jacobTotInv[0]*Apar[3]*rdx2R; 
  BstarYdBmagR[5] = -1.936491673103709*b_z[0]*jacobTotInv[0]*Apar[6]*rdx2R; 
  BstarYdBmagR[12] = -0.8660254037844387*b_z[0]*jacobTotInv[0]*Apar[7]*rdx2R; 

  double alphaR[20]; 
  alphaR[0] = -(0.1767766952966368*((3.872983346207417*(BstarYdBmagR[5]*hamilR[20]+BstarYdBmagR[2]*hamilR[12])+1.732050807568877*(BstarYdBmagR[1]*hamilR[5]+BstarYdBmagR[0]*hamilR[2]))*rdy2R+(3.872983346207417*(BstarXdBmagR[5]*hamilR[19]+BstarXdBmagR[1]*hamilR[11])+1.732050807568877*(BstarXdBmagR[2]*hamilR[5]+BstarXdBmagR[0]*hamilR[1]))*rdx2R+8.0*dApardt[0]*q_))/m_; 
  alphaR[1] = -(0.03535533905932736*((19.36491673103708*BstarYdBmagR[2]*hamilR[20]+7.745966692414834*BstarYdBmagR[1]*hamilR[19]+19.36491673103709*BstarYdBmagR[5]*hamilR[12]+8.660254037844386*(BstarYdBmagR[0]*hamilR[5]+BstarYdBmagR[1]*hamilR[2]))*rdy2R+(19.36491673103708*BstarXdBmagR[2]*hamilR[19]+(17.32050807568877*BstarXdBmagR[11]+19.36491673103709*BstarXdBmagR[0])*hamilR[11]+8.660254037844386*(BstarXdBmagR[5]*hamilR[5]+BstarXdBmagR[1]*hamilR[1]))*rdx2R+40.0*dApardt[1]*q_))/m_; 
  alphaR[2] = -(0.03535533905932736*((19.36491673103708*BstarYdBmagR[1]*hamilR[20]+(17.32050807568877*BstarYdBmagR[12]+19.36491673103709*BstarYdBmagR[0])*hamilR[12]+8.660254037844386*(BstarYdBmagR[5]*hamilR[5]+BstarYdBmagR[2]*hamilR[2]))*rdy2R+(7.745966692414834*BstarXdBmagR[2]*hamilR[20]+19.36491673103708*BstarXdBmagR[1]*hamilR[19]+19.36491673103709*BstarXdBmagR[5]*hamilR[11]+8.660254037844386*(BstarXdBmagR[0]*hamilR[5]+hamilR[1]*BstarXdBmagR[2]))*rdx2R+40.0*dApardt[2]*q_))/m_; 
  alphaR[4] = -(0.03535533905932736*(((17.32050807568877*BstarYdBmagR[12]+19.36491673103708*BstarYdBmagR[0])*hamilR[20]+7.745966692414834*BstarYdBmagR[5]*hamilR[19]+19.36491673103709*BstarYdBmagR[1]*hamilR[12]+8.660254037844386*(BstarYdBmagR[2]*hamilR[5]+hamilR[2]*BstarYdBmagR[5]))*rdy2R+(7.745966692414834*BstarXdBmagR[5]*hamilR[20]+(17.32050807568877*BstarXdBmagR[11]+19.36491673103708*BstarXdBmagR[0])*hamilR[19]+19.36491673103709*BstarXdBmagR[2]*hamilR[11]+8.660254037844386*(BstarXdBmagR[1]*hamilR[5]+hamilR[1]*BstarXdBmagR[5]))*rdx2R+40.0*dApardt[3]*q_))/m_; 
  alphaR[7] = -(0.03535533905932736*((17.32050807568877*BstarYdBmagR[5]*hamilR[20]+8.660254037844387*BstarYdBmagR[0]*hamilR[19]+7.745966692414834*BstarYdBmagR[1]*hamilR[5])*rdy2R+(17.32050807568877*BstarXdBmagR[5]*hamilR[19]+17.32050807568877*BstarXdBmagR[1]*hamilR[11]+8.660254037844386*hamilR[1]*BstarXdBmagR[11])*rdx2R+40.0*dApardt[4]*q_))/m_; 
  alphaR[8] = -(0.03535533905932736*((17.32050807568877*BstarYdBmagR[5]*hamilR[20]+17.32050807568877*BstarYdBmagR[2]*hamilR[12]+8.660254037844386*hamilR[2]*BstarYdBmagR[12])*rdy2R+(8.660254037844387*BstarXdBmagR[0]*hamilR[20]+17.32050807568877*BstarXdBmagR[5]*hamilR[19]+7.745966692414834*BstarXdBmagR[2]*hamilR[5])*rdx2R+40.0*dApardt[5]*q_))/m_; 
  alphaR[11] = -(0.03535533905932736*((17.32050807568877*BstarYdBmagR[1]*hamilR[20]+8.660254037844386*BstarYdBmagR[2]*hamilR[19]+7.745966692414834*BstarYdBmagR[5]*hamilR[5])*rdy2R+(17.32050807568877*BstarXdBmagR[1]*hamilR[19]+17.32050807568877*BstarXdBmagR[5]*hamilR[11]+8.660254037844387*hamilR[5]*BstarXdBmagR[11])*rdx2R+40.0*dApardt[6]*q_))/m_; 
  alphaR[12] = -(0.03535533905932736*((17.32050807568877*BstarYdBmagR[2]*hamilR[20]+17.32050807568877*BstarYdBmagR[5]*hamilR[12]+8.660254037844387*hamilR[5]*BstarYdBmagR[12])*rdy2R+(8.660254037844386*BstarXdBmagR[1]*hamilR[20]+17.32050807568877*BstarXdBmagR[2]*hamilR[19]+7.745966692414834*BstarXdBmagR[5]*hamilR[5])*rdx2R+40.0*dApardt[7]*q_))/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*((3.872983346207417*BstarYdBmagR[5]*hamilR[20]+3.872983346207417*BstarYdBmagR[2]*hamilR[12]+1.732050807568877*BstarYdBmagR[1]*hamilR[5]+1.732050807568877*BstarYdBmagR[0]*hamilR[2])*rdy2R+(3.872983346207417*BstarXdBmagR[5]*hamilR[19]+3.872983346207417*BstarXdBmagR[1]*hamilR[11]+1.732050807568877*BstarXdBmagR[2]*hamilR[5]+1.732050807568877*BstarXdBmagR[0]*hamilR[1])*rdx2R+8.0*dApardt[0]*q_))/m_; 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[13]+fL[13])+1.732050807568877*(fL[3]-1.0*fR[3])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[23]+fL[23])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[24]+fL[24])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[2]+fL[2])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[27]+fL[27])+3.0*(fL[10]-1.0*fR[10]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[34]+fL[34])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[5]+fL[5]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[39]+fL[39])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[40]+fL[40])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[9]+fL[9]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[21]-1.0*(8.660254037844387*fL[21]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[22]-1.0*(8.660254037844387*fL[22]+5.0*(fR[12]+fL[12]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[30]-1.0*(8.660254037844387*fL[30]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[46]+fL[46])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[16]+fL[16])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[19]+fL[19]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[20]+fL[20]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[37]-1.0*(8.660254037844387*fL[37]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[38]-1.0*(8.660254037844387*fL[38]+5.0*(fR[26]+fL[26]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[29]+fL[29]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[35]+fL[35]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[36]+fL[36]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[41]+fL[41]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[13]-6.324555320336761*fL[13]-4.898979485566357*(fR[3]+fL[3])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[12]*fAvg[12]+alphaR[11]*fAvg[11]+alphaR[8]*fAvg[8]+alphaR[7]*fAvg[7]+alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[23]-94.8683298050514*fL[23]-73.48469228349535*(fR[6]+fL[6])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-21.21320343559643*(alphaR[8]*fAvg[12]+fAvg[8]*alphaR[12])-18.97366596101028*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11]+alphaR[1]*fAvg[7]+fAvg[1]*alphaR[7])-21.21320343559643*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[24]-94.8683298050514*fL[24]-73.48469228349535*(fR[7]+fL[7])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*(alphaR[4]*fAvg[12]+fAvg[4]*alphaR[12])-21.21320343559643*(alphaR[7]*fAvg[11]+fAvg[7]*alphaR[11])-18.97366596101028*(alphaR[2]*fAvg[8]+fAvg[2]*alphaR[8])-21.21320343559643*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.008333333333333333*((94.8683298050514*fR[27]-94.8683298050514*fL[27]-73.48469228349535*(fR[10]+fL[10])+42.42640687119286*fR[4]-42.42640687119286*fL[4])*amax-21.21320343559643*(alphaR[12]*fAvg[18]+alphaR[11]*fAvg[17]+alphaR[8]*fAvg[14]+alphaR[7]*fAvg[13])-21.21320343559643*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6]+alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[34]-94.86832980505142*fL[34]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[5]-42.42640687119286*fL[5])*amax+((-16.97056274847715*alphaR[11])-18.97366596101028*alphaR[2])*fAvg[12]+((-16.97056274847715*fAvg[11])-18.97366596101028*fAvg[2])*alphaR[12]-18.97366596101028*(alphaR[1]*fAvg[11]+fAvg[1]*alphaR[11]+alphaR[4]*fAvg[8]+fAvg[4]*alphaR[8]+alphaR[4]*fAvg[7]+fAvg[4]*alphaR[7])-21.21320343559643*(alphaR[0]*fAvg[4]+fAvg[0]*alphaR[4]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[39]-94.86832980505142*fL[39]-73.48469228349535*(fR[17]+fL[17])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-21.21320343559643*alphaR[8]*fAvg[18]-18.97366596101028*alphaR[4]*fAvg[17]-21.21320343559643*alphaR[12]*fAvg[14]-18.97366596101028*alphaR[1]*fAvg[13]+fAvg[10]*((-18.97366596101028*alphaR[11])-21.21320343559643*alphaR[2])-18.97366596101028*fAvg[5]*alphaR[7]-21.21320343559643*(alphaR[4]*fAvg[6]+alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[40]-94.86832980505142*fL[40]-73.48469228349535*(fR[18]+fL[18])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*alphaR[4]*fAvg[18]-21.21320343559643*alphaR[7]*fAvg[17]-18.97366596101028*alphaR[2]*fAvg[14]-21.21320343559643*alphaR[11]*fAvg[13]+fAvg[10]*((-18.97366596101028*alphaR[12])-21.21320343559643*alphaR[1])-18.97366596101028*fAvg[6]*alphaR[8]-21.21320343559643*(alphaR[0]*fAvg[6]+alphaR[4]*fAvg[5]+alphaR[2]*fAvg[3])); 
  Ghat[7] = 0.001190476190476191*((514.3928459844675*(fR[21]+fL[21])-296.98484809835*fR[11]+296.98484809835*fL[11])*amax+132.815661727072*alphaR[12]*fAvg[12]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[11]+148.492424049175*fAvg[2]*alphaR[11]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[7]+148.492424049175*fAvg[0]*alphaR[7]+132.815661727072*(alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1])); 
  Ghat[8] = 0.001190476190476191*((514.3928459844675*(fR[22]+fL[22])-296.98484809835*fR[12]+296.98484809835*fL[12])*amax+(94.86832980505142*alphaR[12]+148.492424049175*alphaR[1])*fAvg[12]+148.492424049175*fAvg[1]*alphaR[12]+132.815661727072*alphaR[11]*fAvg[11]+(94.86832980505142*alphaR[8]+148.492424049175*alphaR[0])*fAvg[8]+148.492424049175*fAvg[0]*alphaR[8]+132.815661727072*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[30]+fL[30])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[4]*fAvg[19]+21.21320343559643*(alphaR[2]*fAvg[16]+alphaR[1]*fAvg[15])+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.001666666666666667*((474.341649025257*fR[46]-474.341649025257*fL[46]-367.4234614174767*(fR[31]+fL[31])+212.1320343559643*fR[16]-212.1320343559643*fL[16])*amax+((-84.85281374238573*alphaR[11])-94.86832980505142*alphaR[2])*fAvg[18]+((-84.85281374238573*alphaR[12])-94.86832980505142*alphaR[1])*fAvg[17]-94.8683298050514*(alphaR[4]*(fAvg[14]+fAvg[13])+fAvg[6]*alphaR[12]+fAvg[5]*alphaR[11])+((-94.86832980505142*(alphaR[8]+alphaR[7]))-106.0660171779821*alphaR[0])*fAvg[10]-106.0660171779821*(alphaR[1]*fAvg[6]+alphaR[2]*fAvg[5]+fAvg[3]*alphaR[4])); 
  Ghat[11] = 0.001190476190476191*((514.3928459844675*(fR[32]+fL[32])-296.98484809835*fR[19]+296.98484809835*fL[19])*amax+118.79393923934*(alphaR[4]*fAvg[12]+fAvg[4]*alphaR[12])+(132.815661727072*alphaR[8]+94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[11]+(132.815661727072*fAvg[8]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaR[11]+148.492424049175*(alphaR[2]*fAvg[7]+fAvg[2]*alphaR[7])+132.815661727072*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])); 
  Ghat[12] = 0.001190476190476191*((514.3928459844675*(fR[33]+fL[33])-296.98484809835*fR[20]+296.98484809835*fL[20])*amax+(94.86832980505142*alphaR[8]+132.815661727072*alphaR[7]+148.492424049175*alphaR[0])*fAvg[12]+(94.86832980505142*fAvg[8]+132.815661727072*fAvg[7]+148.492424049175*fAvg[0])*alphaR[12]+118.79393923934*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11])+148.492424049175*(alphaR[1]*fAvg[8]+fAvg[1]*alphaR[8])+132.815661727072*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])); 
  Ghat[13] = 0.001190476190476191*((514.3928459844675*(fR[37]+fL[37])-296.98484809835*fR[25]+296.98484809835*fL[25])*amax+132.815661727072*alphaR[12]*fAvg[18]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[17]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[13]+148.492424049175*fAvg[6]*alphaR[11]+132.815661727072*alphaR[4]*fAvg[10]+148.492424049175*fAvg[3]*alphaR[7]+132.815661727072*alphaR[1]*fAvg[5]); 
  Ghat[14] = 0.001190476190476191*((514.3928459844675*(fR[38]+fL[38])-296.98484809835*fR[26]+296.98484809835*fL[26])*amax+(94.86832980505142*alphaR[12]+148.492424049175*alphaR[1])*fAvg[18]+132.815661727072*alphaR[11]*fAvg[17]+(94.86832980505142*alphaR[8]+148.492424049175*alphaR[0])*fAvg[14]+148.492424049175*fAvg[5]*alphaR[12]+132.815661727072*alphaR[4]*fAvg[10]+148.492424049175*fAvg[3]*alphaR[8]+132.815661727072*alphaR[2]*fAvg[6]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+(18.97366596101028*alphaR[11]+21.21320343559643*alphaR[2])*fAvg[19]+21.21320343559643*alphaR[4]*fAvg[16]+(18.97366596101028*alphaR[7]+21.21320343559643*alphaR[0])*fAvg[15]+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+(18.97366596101028*alphaR[12]+21.21320343559643*alphaR[1])*fAvg[19]+(18.97366596101028*alphaR[8]+21.21320343559643*alphaR[0])*fAvg[16]+21.21320343559643*alphaR[4]*fAvg[15]+21.21320343559643*alphaR[2]*fAvg[9]); 
  Ghat[17] = 2.380952380952381e-4*((2571.964229922338*(fR[44]+fL[44])-1484.92424049175*fR[35]+1484.92424049175*fL[35])*amax+593.9696961967002*alphaR[4]*fAvg[18]+(664.0783086353599*alphaR[8]+474.3416490252571*alphaR[7]+742.462120245875*alphaR[0])*fAvg[17]+664.0783086353599*alphaR[11]*fAvg[14]+(474.3416490252571*alphaR[11]+742.462120245875*alphaR[2])*fAvg[13]+593.9696961967002*fAvg[10]*alphaR[12]+742.462120245875*fAvg[3]*alphaR[11]+664.0783086353599*alphaR[1]*fAvg[10]+742.462120245875*fAvg[6]*alphaR[7]+664.0783086353599*alphaR[4]*fAvg[5]); 
  Ghat[18] = 2.380952380952381e-4*((2571.964229922338*(fR[45]+fL[45])-1484.92424049175*fR[36]+1484.92424049175*fL[36])*amax+(474.3416490252571*alphaR[8]+664.0783086353599*alphaR[7]+742.462120245875*alphaR[0])*fAvg[18]+593.9696961967002*alphaR[4]*fAvg[17]+(474.3416490252571*alphaR[12]+742.462120245875*alphaR[1])*fAvg[14]+alphaR[12]*(664.0783086353599*fAvg[13]+742.462120245875*fAvg[3])+fAvg[10]*(593.9696961967002*alphaR[11]+664.0783086353599*alphaR[2])+742.462120245875*fAvg[5]*alphaR[8]+664.0783086353599*alphaR[4]*fAvg[6]); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[41]+42.42640687119286*fL[41])*amax+(18.97366596101028*(alphaR[8]+alphaR[7])+21.21320343559643*alphaR[0])*fAvg[19]+(18.97366596101028*alphaR[12]+21.21320343559643*alphaR[1])*fAvg[16]+(18.97366596101028*alphaR[11]+21.21320343559643*alphaR[2])*fAvg[15]+21.21320343559643*alphaR[4]*fAvg[9]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = 0.7071067811865475*Ghat[2]; 
  incr[3] = -1.224744871391589*Ghat[0]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = 0.7071067811865475*Ghat[4]; 
  incr[6] = -1.224744871391589*Ghat[1]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = 0.7071067811865475*Ghat[6]; 
  incr[10] = -1.224744871391589*Ghat[3]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 0.7071067811865475*Ghat[8]; 
  incr[13] = 1.58113883008419*Ghat[0]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = 0.7071067811865475*Ghat[10]; 
  incr[17] = -1.224744871391589*Ghat[5]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = 0.7071067811865475*Ghat[11]; 
  incr[20] = 0.7071067811865475*Ghat[12]; 
  incr[21] = -1.224744871391589*Ghat[7]; 
  incr[22] = -1.224744871391589*Ghat[8]; 
  incr[23] = 1.58113883008419*Ghat[1]; 
  incr[24] = 1.58113883008419*Ghat[2]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 0.7071067811865475*Ghat[14]; 
  incr[27] = 1.58113883008419*Ghat[3]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = 0.7071067811865475*Ghat[16]; 
  incr[30] = -1.224744871391589*Ghat[9]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = -1.224744871391589*Ghat[12]; 
  incr[34] = 1.58113883008419*Ghat[4]; 
  incr[35] = 0.7071067811865475*Ghat[17]; 
  incr[36] = 0.7071067811865475*Ghat[18]; 
  incr[37] = -1.224744871391589*Ghat[13]; 
  incr[38] = -1.224744871391589*Ghat[14]; 
  incr[39] = 1.58113883008419*Ghat[5]; 
  incr[40] = 1.58113883008419*Ghat[6]; 
  incr[41] = 0.7071067811865475*Ghat[19]; 
  incr[42] = -1.224744871391589*Ghat[15]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = -1.224744871391589*Ghat[18]; 
  incr[46] = 1.58113883008419*Ghat[10]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

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
  outR[16] += incr[16]*rdvpar2R; 
  outR[17] += incr[17]*rdvpar2R; 
  outR[18] += incr[18]*rdvpar2R; 
  outR[19] += incr[19]*rdvpar2R; 
  outR[20] += incr[20]*rdvpar2R; 
  outR[21] += incr[21]*rdvpar2R; 
  outR[22] += incr[22]*rdvpar2R; 
  outR[23] += incr[23]*rdvpar2R; 
  outR[24] += incr[24]*rdvpar2R; 
  outR[25] += incr[25]*rdvpar2R; 
  outR[26] += incr[26]*rdvpar2R; 
  outR[27] += incr[27]*rdvpar2R; 
  outR[28] += incr[28]*rdvpar2R; 
  outR[29] += incr[29]*rdvpar2R; 
  outR[30] += incr[30]*rdvpar2R; 
  outR[31] += incr[31]*rdvpar2R; 
  outR[32] += incr[32]*rdvpar2R; 
  outR[33] += incr[33]*rdvpar2R; 
  outR[34] += incr[34]*rdvpar2R; 
  outR[35] += incr[35]*rdvpar2R; 
  outR[36] += incr[36]*rdvpar2R; 
  outR[37] += incr[37]*rdvpar2R; 
  outR[38] += incr[38]*rdvpar2R; 
  outR[39] += incr[39]*rdvpar2R; 
  outR[40] += incr[40]*rdvpar2R; 
  outR[41] += incr[41]*rdvpar2R; 
  outR[42] += incr[42]*rdvpar2R; 
  outR[43] += incr[43]*rdvpar2R; 
  outR[44] += incr[44]*rdvpar2R; 
  outR[45] += incr[45]*rdvpar2R; 
  outR[46] += incr[46]*rdvpar2R; 
  outR[47] += incr[47]*rdvpar2R; 

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
  outL[11] += -1.0*incr[11]*rdvpar2L; 
  outL[12] += -1.0*incr[12]*rdvpar2L; 
  outL[13] += -1.0*incr[13]*rdvpar2L; 
  outL[14] += -1.0*incr[14]*rdvpar2L; 
  outL[15] += incr[15]*rdvpar2L; 
  outL[16] += -1.0*incr[16]*rdvpar2L; 
  outL[17] += incr[17]*rdvpar2L; 
  outL[18] += incr[18]*rdvpar2L; 
  outL[19] += -1.0*incr[19]*rdvpar2L; 
  outL[20] += -1.0*incr[20]*rdvpar2L; 
  outL[21] += incr[21]*rdvpar2L; 
  outL[22] += incr[22]*rdvpar2L; 
  outL[23] += -1.0*incr[23]*rdvpar2L; 
  outL[24] += -1.0*incr[24]*rdvpar2L; 
  outL[25] += -1.0*incr[25]*rdvpar2L; 
  outL[26] += -1.0*incr[26]*rdvpar2L; 
  outL[27] += -1.0*incr[27]*rdvpar2L; 
  outL[28] += -1.0*incr[28]*rdvpar2L; 
  outL[29] += -1.0*incr[29]*rdvpar2L; 
  outL[30] += incr[30]*rdvpar2L; 
  outL[31] += incr[31]*rdvpar2L; 
  outL[32] += incr[32]*rdvpar2L; 
  outL[33] += incr[33]*rdvpar2L; 
  outL[34] += -1.0*incr[34]*rdvpar2L; 
  outL[35] += -1.0*incr[35]*rdvpar2L; 
  outL[36] += -1.0*incr[36]*rdvpar2L; 
  outL[37] += incr[37]*rdvpar2L; 
  outL[38] += incr[38]*rdvpar2L; 
  outL[39] += -1.0*incr[39]*rdvpar2L; 
  outL[40] += -1.0*incr[40]*rdvpar2L; 
  outL[41] += -1.0*incr[41]*rdvpar2L; 
  outL[42] += incr[42]*rdvpar2L; 
  outL[43] += incr[43]*rdvpar2L; 
  outL[44] += incr[44]*rdvpar2L; 
  outL[45] += incr[45]*rdvpar2L; 
  outL[46] += -1.0*incr[46]*rdvpar2L; 
  outL[47] += incr[47]*rdvpar2L; 

  return std::abs(alphaSurfAvgR);
} 
double EmGyrokineticGenGeoSurf2x2vSer_xL_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilL[48]; 
  hamilL[0] = (0.6666666666666666*(3.0*rdvpar2SqL*(m_*wvparSqL+bmag[0]*wmuL+phi[0]*q_)+m_))/rdvpar2SqL; 
  hamilL[1] = 2.0*(bmag[1]*wmuL+phi[1]*q_); 
  hamilL[2] = 2.0*phi[2]*q_; 
  hamilL[3] = (2.309401076758503*m_*wvparL)/rdvpar2L; 
  hamilL[4] = (1.154700538379252*bmag[0])/rdmu2L; 
  hamilL[5] = 2.0*phi[3]*q_; 
  hamilL[8] = (1.154700538379252*bmag[1])/rdmu2L; 
  hamilL[11] = 2.0*(bmag[4]*wmuL+phi[4]*q_); 
  hamilL[12] = 2.0*phi[5]*q_; 
  hamilL[13] = (0.5962847939999438*m_)/rdvpar2SqL; 
  hamilL[19] = 2.0*phi[6]*q_; 
  hamilL[20] = 2.0*phi[7]*q_; 
  hamilL[25] = (1.154700538379251*bmag[4])/rdmu2L; 

  double BstarXdBmagL[48]; 
  BstarXdBmagL[0] = 0.01428571428571429*(3.872983346207417*((10.0*b_z[4]+15.65247584249853*b_z[0])*jacobTotInv[4]+7.0*(2.23606797749979*jacobTotInv[0]*b_z[4]+2.0*b_z[1]*jacobTotInv[1]))*Apar[6]+12.12435565298214*((5.0*Apar[2]*b_z[4]+4.47213595499958*b_z[1]*Apar[3])*jacobTotInv[4]+4.47213595499958*jacobTotInv[1]*Apar[3]*b_z[4]+5.0*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])))*rdy2L; 
  BstarXdBmagL[1] = 0.01428571428571429*(3.872983346207417*(24.59674775249769*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+14.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[6]+1.732050807568877*((55.0*Apar[3]*b_z[4]+31.30495168499706*(b_z[0]*Apar[3]+b_z[1]*Apar[2]))*jacobTotInv[4]+7.0*(4.47213595499958*(jacobTotInv[0]*Apar[3]+jacobTotInv[1]*Apar[2])*b_z[4]+(9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])))*rdy2L; 
  BstarXdBmagL[2] = 0.1*(3.872983346207417*(4.47213595499958*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[7]+19.36491673103709*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[5])*rdy2L; 
  BstarXdBmagL[5] = 0.01428571428571429*(3.872983346207417*((55.0*b_z[4]+31.30495168499706*b_z[0])*jacobTotInv[4]+7.0*(4.47213595499958*jacobTotInv[0]*b_z[4]+9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0]))*Apar[7]+12.12435565298214*(10.0*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+11.18033988749895*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[5])*rdy2L; 
  BstarXdBmagL[11] = 0.01428571428571429*(3.872983346207417*((33.54101966249685*b_z[4]+10.0*b_z[0])*jacobTotInv[4]+10.0*jacobTotInv[0]*b_z[4]+2.23606797749979*(11.0*b_z[1]*jacobTotInv[1]+7.0*b_z[0]*jacobTotInv[0]))*Apar[6]+1.732050807568877*((22.3606797749979*Apar[2]*b_z[4]+5.0*(11.0*b_z[1]*Apar[3]+7.0*b_z[0]*Apar[2]))*jacobTotInv[4]+5.0*(11.0*jacobTotInv[1]*Apar[3]+7.0*jacobTotInv[0]*Apar[2])*b_z[4]+31.30495168499706*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])))*rdy2L; 
  BstarXdBmagL[19] = 0.01428571428571429*(1.732050807568877*(122.9837387624885*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+70.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[7]+3.872983346207417*((22.3606797749979*b_z[4]+35.0*b_z[0])*jacobTotInv[4]+7.0*(5.0*jacobTotInv[0]*b_z[4]+4.47213595499958*b_z[1]*jacobTotInv[1]))*Apar[5])*rdy2L; 

  double alphaL[20]; 
  alphaL[0] = -(0.1767766952966368*((((19.36491673103708*b_z[4]+15.0*b_z[1]+8.660254037844387*b_z[0])*jacobTotInv[4]+(15.0*jacobTotInv[1]+8.660254037844387*jacobTotInv[0])*b_z[4]+(11.61895003862225*b_z[1]+6.708203932499369*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(6.708203932499369*b_z[1]+3.872983346207417*b_z[0]))*hamilL[19]+((15.0*b_z[4]+11.61895003862225*b_z[1]+6.708203932499369*b_z[0])*jacobTotInv[4]+(11.61895003862225*jacobTotInv[1]+6.708203932499369*jacobTotInv[0])*b_z[4]+(9.0*b_z[1]+5.196152422706631*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(5.196152422706631*b_z[1]+3.0*b_z[0]))*hamilL[5]+hamilL[2]*((8.660254037844386*b_z[4]+6.708203932499369*b_z[1]+3.872983346207417*b_z[0])*jacobTotInv[4]+(6.708203932499369*jacobTotInv[1]+3.872983346207417*jacobTotInv[0])*b_z[4]+(5.196152422706631*b_z[1]+3.0*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(3.0*b_z[1]+1.732050807568877*b_z[0])))*m_*rdy2L+hamilL[3]*((-3.872983346207417*BstarXdBmagL[11])-3.0*BstarXdBmagL[1]-1.732050807568877*BstarXdBmagL[0])*q_*rdvpar2L))/(m_*q_); 
  alphaL[1] = -(0.1767766952966368*((((33.54101966249684*b_z[4]+25.98076211353316*b_z[1]+15.0*b_z[0])*jacobTotInv[4]+(25.98076211353316*jacobTotInv[1]+15.0*jacobTotInv[0])*b_z[4]+(20.1246117974981*b_z[1]+11.61895003862225*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(11.61895003862225*b_z[1]+6.708203932499369*b_z[0]))*hamilL[20]+((19.36491673103709*b_z[4]+15.0*b_z[1]+8.660254037844386*b_z[0])*jacobTotInv[4]+(15.0*jacobTotInv[1]+8.660254037844386*jacobTotInv[0])*b_z[4]+(11.61895003862225*b_z[1]+6.708203932499369*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(6.708203932499369*b_z[1]+3.872983346207417*b_z[0]))*hamilL[12])*m_*rdy2L+hamilL[3]*((-3.872983346207417*BstarXdBmagL[19])-3.0*BstarXdBmagL[5]-1.732050807568877*BstarXdBmagL[2])*q_*rdvpar2L))/(m_*q_); 
  alphaL[2] = (0.1767766952966368*(8.660254037844386*BstarXdBmagL[11]+6.708203932499369*BstarXdBmagL[1]+3.872983346207417*BstarXdBmagL[0])*hamilL[13]*rdvpar2L)/m_; 
  alphaL[4] = (0.1767766952966368*hamilL[13]*(8.660254037844387*BstarXdBmagL[19]+6.708203932499369*BstarXdBmagL[5]+3.872983346207417*BstarXdBmagL[2])*rdvpar2L)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgL = -(0.03125*((((19.36491673103708*b_z[4]+15.0*b_z[1]+8.660254037844387*b_z[0])*jacobTotInv[4]+(15.0*jacobTotInv[1]+8.660254037844387*jacobTotInv[0])*b_z[4]+(11.61895003862225*b_z[1]+6.708203932499369*b_z[0])*jacobTotInv[1]+6.708203932499369*jacobTotInv[0]*b_z[1]+3.872983346207417*b_z[0]*jacobTotInv[0])*hamilL[19]+((15.0*b_z[4]+11.61895003862225*b_z[1]+6.708203932499369*b_z[0])*jacobTotInv[4]+(11.61895003862225*jacobTotInv[1]+6.708203932499369*jacobTotInv[0])*b_z[4]+(9.0*b_z[1]+5.196152422706631*b_z[0])*jacobTotInv[1]+5.196152422706631*jacobTotInv[0]*b_z[1]+3.0*b_z[0]*jacobTotInv[0])*hamilL[5]+(8.660254037844386*hamilL[2]*b_z[4]+(6.708203932499369*b_z[1]+3.872983346207417*b_z[0])*hamilL[2])*jacobTotInv[4]+(6.708203932499369*jacobTotInv[1]+3.872983346207417*jacobTotInv[0])*hamilL[2]*b_z[4]+((5.196152422706631*b_z[1]+3.0*b_z[0])*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1]+1.732050807568877*b_z[0]*jacobTotInv[0])*hamilL[2])*m_*rdy2L+(((-3.0*BstarXdBmagL[1])-1.732050807568877*BstarXdBmagL[0])*hamilL[3]-3.872983346207417*hamilL[3]*BstarXdBmagL[11])*q_*rdvpar2L))/(m_*q_); 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[11]+fL[11])+1.732050807568877*(fL[1]-1.0*fR[1])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[19]+fL[19])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[2]+fL[2])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[21]+fL[21])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[25]+fL[25])+3.0*(fL[8]-1.0*fR[8]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[32]+fL[32])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[7]+fL[7]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[35]+fL[35])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[9]+fL[9]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[37]+fL[37])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[20]-1.0*(8.660254037844387*fL[20]+5.0*(fR[12]+fL[12]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[23]-1.0*(8.660254037844387*fL[23]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[28]-1.0*(8.660254037844387*fL[28]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[44]+fL[44])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[18]+fL[18])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[22]+fL[22]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[24]+fL[24]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[36]-1.0*(8.660254037844387*fL[36]+5.0*(fR[26]+fL[26]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[39]-1.0*(8.660254037844387*fL[39]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[29]+fL[29]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[38]+fL[38]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[40]+fL[40]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[43]+fL[43]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[11]-6.324555320336761*fL[11]-4.898979485566357*(fR[1]+fL[1])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaL[4]*fAvg[4]+alphaL[2]*fAvg[2]+alphaL[1]*fAvg[1]+alphaL[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[19]-94.8683298050514*fL[19]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*(alphaL[4]*fAvg[11]+alphaL[1]*fAvg[7])-21.21320343559643*(alphaL[2]*fAvg[4]+fAvg[2]*alphaL[4]+alphaL[0]*fAvg[1]+fAvg[0]*alphaL[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[21]-94.8683298050514*fL[21]-73.48469228349535*(fR[6]+fL[6])+42.42640687119286*fR[3]-42.42640687119286*fL[3])*amax-18.97366596101028*(alphaL[4]*fAvg[12]+alphaL[2]*fAvg[8])-21.21320343559643*(alphaL[1]*fAvg[4]+fAvg[1]*alphaL[4]+alphaL[0]*fAvg[2]+fAvg[0]*alphaL[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[25]-18.97366596101028*fL[25]-14.69693845669907*(fR[8]+fL[8])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaL[4]*fAvg[10]+alphaL[2]*fAvg[6]+alphaL[1]*fAvg[5]+alphaL[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[32]-94.86832980505142*fL[32]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[7]-42.42640687119286*fL[7])*amax-18.97366596101028*(alphaL[2]*fAvg[12]+alphaL[1]*fAvg[11]+alphaL[4]*(fAvg[8]+fAvg[7]))-21.21320343559643*(alphaL[0]*fAvg[4]+fAvg[0]*alphaL[4]+alphaL[1]*fAvg[2]+fAvg[1]*alphaL[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[35]-94.86832980505142*fL[35]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*(alphaL[4]*fAvg[17]+alphaL[1]*fAvg[13])-21.21320343559643*(alphaL[2]*fAvg[10]+alphaL[4]*fAvg[6]+alphaL[0]*fAvg[5]+alphaL[1]*fAvg[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[37]-94.86832980505142*fL[37]-73.48469228349535*(fR[17]+fL[17])+42.42640687119286*fR[10]-42.42640687119286*fL[10])*amax-18.97366596101028*(alphaL[4]*fAvg[18]+alphaL[2]*fAvg[14])-21.21320343559643*(alphaL[1]*fAvg[10]+alphaL[0]*fAvg[6]+alphaL[4]*fAvg[5]+alphaL[2]*fAvg[3])); 
  Ghat[7] = 0.008333333333333333*((73.48469228349536*(fR[20]+fL[20])-42.42640687119286*fR[12]+42.42640687119286*fL[12])*amax+21.21320343559643*alphaL[2]*fAvg[11]+21.21320343559643*alphaL[0]*fAvg[7]+18.97366596101028*(alphaL[4]*fAvg[4]+alphaL[1]*fAvg[1])); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[23]+fL[23])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaL[1]*fAvg[12]+21.21320343559643*alphaL[0]*fAvg[8]+18.97366596101028*(alphaL[4]*fAvg[4]+alphaL[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[28]+fL[28])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaL[4]*fAvg[19]+21.21320343559643*(alphaL[2]*fAvg[16]+alphaL[1]*fAvg[15])+21.21320343559643*alphaL[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[44]-94.8683298050514*fL[44]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[18]-42.42640687119286*fL[18])*amax-18.97366596101028*(alphaL[2]*fAvg[18]+alphaL[1]*fAvg[17]+alphaL[4]*(fAvg[14]+fAvg[13]))-21.21320343559643*(alphaL[0]*fAvg[10]+alphaL[1]*fAvg[6]+alphaL[2]*fAvg[5]+fAvg[3]*alphaL[4])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[33]+fL[33])-42.42640687119286*fR[22]+42.42640687119286*fL[22])*amax+16.97056274847715*alphaL[4]*fAvg[12]+21.21320343559643*alphaL[0]*fAvg[11]+21.21320343559643*alphaL[2]*fAvg[7]+18.97366596101028*(alphaL[1]*fAvg[4]+fAvg[1]*alphaL[4])); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[24]+42.42640687119286*fL[24])*amax+21.21320343559643*alphaL[0]*fAvg[12]+16.97056274847715*alphaL[4]*fAvg[11]+21.21320343559643*alphaL[1]*fAvg[8]+18.97366596101028*(alphaL[2]*fAvg[4]+fAvg[2]*alphaL[4])); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[36]+fL[36])-42.42640687119286*fR[26]+42.42640687119286*fL[26])*amax+21.21320343559643*alphaL[2]*fAvg[17]+21.21320343559643*alphaL[0]*fAvg[13]+18.97366596101028*(alphaL[4]*fAvg[10]+alphaL[1]*fAvg[5])); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[39]+fL[39])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaL[1]*fAvg[18]+21.21320343559643*alphaL[0]*fAvg[14]+18.97366596101028*(alphaL[4]*fAvg[10]+alphaL[2]*fAvg[6])); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+21.21320343559643*alphaL[2]*fAvg[19]+21.21320343559643*(alphaL[4]*fAvg[16]+alphaL[0]*fAvg[15])+21.21320343559643*alphaL[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaL[1]*fAvg[19]+21.21320343559643*(alphaL[0]*fAvg[16]+alphaL[4]*fAvg[15])+21.21320343559643*alphaL[2]*fAvg[9]); 
  Ghat[17] = 0.008333333333333333*((73.48469228349536*(fR[45]+fL[45])-42.42640687119286*fR[38]+42.42640687119286*fL[38])*amax+16.97056274847715*alphaL[4]*fAvg[18]+21.21320343559643*alphaL[0]*fAvg[17]+21.21320343559643*alphaL[2]*fAvg[13]+18.97366596101028*(alphaL[1]*fAvg[10]+alphaL[4]*fAvg[5])); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[40]+42.42640687119286*fL[40])*amax+21.21320343559643*alphaL[0]*fAvg[18]+16.97056274847715*alphaL[4]*fAvg[17]+21.21320343559643*alphaL[1]*fAvg[14]+18.97366596101028*(alphaL[2]*fAvg[10]+alphaL[4]*fAvg[6])); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[43]+42.42640687119286*fL[43])*amax+21.21320343559643*alphaL[0]*fAvg[19]+21.21320343559643*(alphaL[1]*fAvg[16]+alphaL[2]*fAvg[15])+21.21320343559643*alphaL[4]*fAvg[9]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = -1.224744871391589*Ghat[0]; 
  incr[2] = 0.7071067811865475*Ghat[1]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = -1.224744871391589*Ghat[2]; 
  incr[7] = 0.7071067811865475*Ghat[4]; 
  incr[8] = -1.224744871391589*Ghat[3]; 
  incr[9] = 0.7071067811865475*Ghat[5]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 1.58113883008419*Ghat[0]; 
  incr[12] = 0.7071067811865475*Ghat[7]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = -1.224744871391589*Ghat[6]; 
  incr[18] = 0.7071067811865475*Ghat[10]; 
  incr[19] = 1.58113883008419*Ghat[1]; 
  incr[20] = -1.224744871391589*Ghat[7]; 
  incr[21] = 1.58113883008419*Ghat[2]; 
  incr[22] = 0.7071067811865475*Ghat[11]; 
  incr[23] = -1.224744871391589*Ghat[8]; 
  incr[24] = 0.7071067811865475*Ghat[12]; 
  incr[25] = 1.58113883008419*Ghat[3]; 
  incr[26] = 0.7071067811865475*Ghat[13]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = -1.224744871391589*Ghat[9]; 
  incr[29] = 0.7071067811865475*Ghat[15]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = 1.58113883008419*Ghat[4]; 
  incr[33] = -1.224744871391589*Ghat[11]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = 1.58113883008419*Ghat[5]; 
  incr[36] = -1.224744871391589*Ghat[13]; 
  incr[37] = 1.58113883008419*Ghat[6]; 
  incr[38] = 0.7071067811865475*Ghat[17]; 
  incr[39] = -1.224744871391589*Ghat[14]; 
  incr[40] = 0.7071067811865475*Ghat[18]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = -1.224744871391589*Ghat[16]; 
  incr[43] = 0.7071067811865475*Ghat[19]; 
  incr[44] = 1.58113883008419*Ghat[10]; 
  incr[45] = -1.224744871391589*Ghat[17]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

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
  outR[16] += incr[16]*rdx2R; 
  outR[17] += incr[17]*rdx2R; 
  outR[18] += incr[18]*rdx2R; 
  outR[19] += incr[19]*rdx2R; 
  outR[20] += incr[20]*rdx2R; 
  outR[21] += incr[21]*rdx2R; 
  outR[22] += incr[22]*rdx2R; 
  outR[23] += incr[23]*rdx2R; 
  outR[24] += incr[24]*rdx2R; 
  outR[25] += incr[25]*rdx2R; 
  outR[26] += incr[26]*rdx2R; 
  outR[27] += incr[27]*rdx2R; 
  outR[28] += incr[28]*rdx2R; 
  outR[29] += incr[29]*rdx2R; 
  outR[30] += incr[30]*rdx2R; 
  outR[31] += incr[31]*rdx2R; 
  outR[32] += incr[32]*rdx2R; 
  outR[33] += incr[33]*rdx2R; 
  outR[34] += incr[34]*rdx2R; 
  outR[35] += incr[35]*rdx2R; 
  outR[36] += incr[36]*rdx2R; 
  outR[37] += incr[37]*rdx2R; 
  outR[38] += incr[38]*rdx2R; 
  outR[39] += incr[39]*rdx2R; 
  outR[40] += incr[40]*rdx2R; 
  outR[41] += incr[41]*rdx2R; 
  outR[42] += incr[42]*rdx2R; 
  outR[43] += incr[43]*rdx2R; 
  outR[44] += incr[44]*rdx2R; 
  outR[45] += incr[45]*rdx2R; 
  outR[46] += incr[46]*rdx2R; 
  outR[47] += incr[47]*rdx2R; 

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
  outL[11] += -1.0*incr[11]*rdx2L; 
  outL[12] += -1.0*incr[12]*rdx2L; 
  outL[13] += -1.0*incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 
  outL[16] += incr[16]*rdx2L; 
  outL[17] += incr[17]*rdx2L; 
  outL[18] += -1.0*incr[18]*rdx2L; 
  outL[19] += -1.0*incr[19]*rdx2L; 
  outL[20] += incr[20]*rdx2L; 
  outL[21] += -1.0*incr[21]*rdx2L; 
  outL[22] += -1.0*incr[22]*rdx2L; 
  outL[23] += incr[23]*rdx2L; 
  outL[24] += -1.0*incr[24]*rdx2L; 
  outL[25] += -1.0*incr[25]*rdx2L; 
  outL[26] += -1.0*incr[26]*rdx2L; 
  outL[27] += -1.0*incr[27]*rdx2L; 
  outL[28] += incr[28]*rdx2L; 
  outL[29] += -1.0*incr[29]*rdx2L; 
  outL[30] += -1.0*incr[30]*rdx2L; 
  outL[31] += incr[31]*rdx2L; 
  outL[32] += -1.0*incr[32]*rdx2L; 
  outL[33] += incr[33]*rdx2L; 
  outL[34] += incr[34]*rdx2L; 
  outL[35] += -1.0*incr[35]*rdx2L; 
  outL[36] += incr[36]*rdx2L; 
  outL[37] += -1.0*incr[37]*rdx2L; 
  outL[38] += -1.0*incr[38]*rdx2L; 
  outL[39] += incr[39]*rdx2L; 
  outL[40] += -1.0*incr[40]*rdx2L; 
  outL[41] += incr[41]*rdx2L; 
  outL[42] += incr[42]*rdx2L; 
  outL[43] += -1.0*incr[43]*rdx2L; 
  outL[44] += -1.0*incr[44]*rdx2L; 
  outL[45] += incr[45]*rdx2L; 
  outL[46] += incr[46]*rdx2L; 
  outL[47] += incr[47]*rdx2L; 

  return std::abs(alphaSurfAvgL);
} 
double EmGyrokineticGenGeoSurf2x2vSer_xR_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[11] = 2.0*(bmag[4]*wmuR+phi[4]*q_); 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 
  hamilR[25] = (1.154700538379251*bmag[4])/rdmu2R; 

  double BstarXdBmagR[48]; 
  BstarXdBmagR[0] = 0.01428571428571429*(3.872983346207417*((10.0*b_z[4]+15.65247584249853*b_z[0])*jacobTotInv[4]+7.0*(2.23606797749979*jacobTotInv[0]*b_z[4]+2.0*b_z[1]*jacobTotInv[1]))*Apar[6]+12.12435565298214*((5.0*Apar[2]*b_z[4]+4.47213595499958*b_z[1]*Apar[3])*jacobTotInv[4]+4.47213595499958*jacobTotInv[1]*Apar[3]*b_z[4]+5.0*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])))*rdy2R; 
  BstarXdBmagR[1] = 0.01428571428571429*(3.872983346207417*(24.59674775249769*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+14.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[6]+1.732050807568877*((55.0*Apar[3]*b_z[4]+31.30495168499706*(b_z[0]*Apar[3]+b_z[1]*Apar[2]))*jacobTotInv[4]+7.0*(4.47213595499958*(jacobTotInv[0]*Apar[3]+jacobTotInv[1]*Apar[2])*b_z[4]+(9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])))*rdy2R; 
  BstarXdBmagR[2] = 0.1*(3.872983346207417*(4.47213595499958*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[7]+19.36491673103709*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[5])*rdy2R; 
  BstarXdBmagR[5] = 0.01428571428571429*(3.872983346207417*((55.0*b_z[4]+31.30495168499706*b_z[0])*jacobTotInv[4]+7.0*(4.47213595499958*jacobTotInv[0]*b_z[4]+9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0]))*Apar[7]+12.12435565298214*(10.0*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+11.18033988749895*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[5])*rdy2R; 
  BstarXdBmagR[11] = 0.01428571428571429*(3.872983346207417*((33.54101966249685*b_z[4]+10.0*b_z[0])*jacobTotInv[4]+10.0*jacobTotInv[0]*b_z[4]+2.23606797749979*(11.0*b_z[1]*jacobTotInv[1]+7.0*b_z[0]*jacobTotInv[0]))*Apar[6]+1.732050807568877*((22.3606797749979*Apar[2]*b_z[4]+5.0*(11.0*b_z[1]*Apar[3]+7.0*b_z[0]*Apar[2]))*jacobTotInv[4]+5.0*(11.0*jacobTotInv[1]*Apar[3]+7.0*jacobTotInv[0]*Apar[2])*b_z[4]+31.30495168499706*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])))*rdy2R; 
  BstarXdBmagR[19] = 0.01428571428571429*(1.732050807568877*(122.9837387624885*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+70.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[7]+3.872983346207417*((22.3606797749979*b_z[4]+35.0*b_z[0])*jacobTotInv[4]+7.0*(5.0*jacobTotInv[0]*b_z[4]+4.47213595499958*b_z[1]*jacobTotInv[1]))*Apar[5])*rdy2R; 

  double alphaR[20]; 
  alphaR[0] = -(0.1767766952966368*((((19.36491673103708*b_z[4]-15.0*b_z[1]+8.660254037844387*b_z[0])*jacobTotInv[4]+(8.660254037844387*jacobTotInv[0]-15.0*jacobTotInv[1])*b_z[4]+(11.61895003862225*b_z[1]-6.708203932499369*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(3.872983346207417*b_z[0]-6.708203932499369*b_z[1]))*hamilR[19]+(((-15.0*b_z[4])+11.61895003862225*b_z[1]-6.708203932499369*b_z[0])*jacobTotInv[4]+(11.61895003862225*jacobTotInv[1]-6.708203932499369*jacobTotInv[0])*b_z[4]+(5.196152422706631*b_z[0]-9.0*b_z[1])*jacobTotInv[1]+jacobTotInv[0]*(5.196152422706631*b_z[1]-3.0*b_z[0]))*hamilR[5]+hamilR[2]*((8.660254037844386*b_z[4]-6.708203932499369*b_z[1]+3.872983346207417*b_z[0])*jacobTotInv[4]+(3.872983346207417*jacobTotInv[0]-6.708203932499369*jacobTotInv[1])*b_z[4]+(5.196152422706631*b_z[1]-3.0*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(1.732050807568877*b_z[0]-3.0*b_z[1])))*m_*rdy2R+hamilR[3]*((-3.872983346207417*BstarXdBmagR[11])+3.0*BstarXdBmagR[1]-1.732050807568877*BstarXdBmagR[0])*q_*rdvpar2R))/(m_*q_); 
  alphaR[1] = (0.1767766952966368*((((33.54101966249684*b_z[4]-25.98076211353316*b_z[1]+15.0*b_z[0])*jacobTotInv[4]+(15.0*jacobTotInv[0]-25.98076211353316*jacobTotInv[1])*b_z[4]+(20.1246117974981*b_z[1]-11.61895003862225*b_z[0])*jacobTotInv[1]+jacobTotInv[0]*(6.708203932499369*b_z[0]-11.61895003862225*b_z[1]))*hamilR[20]+(((-19.36491673103709*b_z[4])+15.0*b_z[1]-8.660254037844386*b_z[0])*jacobTotInv[4]+(15.0*jacobTotInv[1]-8.660254037844386*jacobTotInv[0])*b_z[4]+(6.708203932499369*b_z[0]-11.61895003862225*b_z[1])*jacobTotInv[1]+jacobTotInv[0]*(6.708203932499369*b_z[1]-3.872983346207417*b_z[0]))*hamilR[12])*m_*rdy2R+hamilR[3]*(3.872983346207417*BstarXdBmagR[19]-3.0*BstarXdBmagR[5]+1.732050807568877*BstarXdBmagR[2])*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = (0.1767766952966368*(8.660254037844386*BstarXdBmagR[11]-6.708203932499369*BstarXdBmagR[1]+3.872983346207417*BstarXdBmagR[0])*hamilR[13]*rdvpar2R)/m_; 
  alphaR[4] = (0.1767766952966368*hamilR[13]*(8.660254037844387*BstarXdBmagR[19]-6.708203932499369*BstarXdBmagR[5]+3.872983346207417*BstarXdBmagR[2])*rdvpar2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = -(0.03125*((((19.36491673103708*b_z[4]-15.0*b_z[1]+8.660254037844387*b_z[0])*jacobTotInv[4]+(8.660254037844387*jacobTotInv[0]-15.0*jacobTotInv[1])*b_z[4]+(11.61895003862225*b_z[1]-6.708203932499369*b_z[0])*jacobTotInv[1]-6.708203932499369*jacobTotInv[0]*b_z[1]+3.872983346207417*b_z[0]*jacobTotInv[0])*hamilR[19]+(((-15.0*b_z[4])+11.61895003862225*b_z[1]-6.708203932499369*b_z[0])*jacobTotInv[4]+(11.61895003862225*jacobTotInv[1]-6.708203932499369*jacobTotInv[0])*b_z[4]+(5.196152422706631*b_z[0]-9.0*b_z[1])*jacobTotInv[1]+5.196152422706631*jacobTotInv[0]*b_z[1]-3.0*b_z[0]*jacobTotInv[0])*hamilR[5]+(8.660254037844386*hamilR[2]*b_z[4]+(3.872983346207417*b_z[0]-6.708203932499369*b_z[1])*hamilR[2])*jacobTotInv[4]+(3.872983346207417*jacobTotInv[0]-6.708203932499369*jacobTotInv[1])*hamilR[2]*b_z[4]+((5.196152422706631*b_z[1]-3.0*b_z[0])*jacobTotInv[1]-3.0*jacobTotInv[0]*b_z[1]+1.732050807568877*b_z[0]*jacobTotInv[0])*hamilR[2])*m_*rdy2R+((3.0*BstarXdBmagR[1]-1.732050807568877*BstarXdBmagR[0])*hamilR[3]-3.872983346207417*hamilR[3]*BstarXdBmagR[11])*q_*rdvpar2R))/(m_*q_); 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[11]+fL[11])+1.732050807568877*(fL[1]-1.0*fR[1])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[19]+fL[19])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[2]+fL[2])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[21]+fL[21])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[25]+fL[25])+3.0*(fL[8]-1.0*fR[8]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[32]+fL[32])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[7]+fL[7]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[35]+fL[35])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[9]+fL[9]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[37]+fL[37])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[20]-1.0*(8.660254037844387*fL[20]+5.0*(fR[12]+fL[12]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[23]-1.0*(8.660254037844387*fL[23]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[28]-1.0*(8.660254037844387*fL[28]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[44]+fL[44])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[18]+fL[18])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[22]+fL[22]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[24]+fL[24]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[36]-1.0*(8.660254037844387*fL[36]+5.0*(fR[26]+fL[26]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[39]-1.0*(8.660254037844387*fL[39]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[29]+fL[29]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[38]+fL[38]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[40]+fL[40]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[43]+fL[43]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[11]-6.324555320336761*fL[11]-4.898979485566357*(fR[1]+fL[1])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[19]-94.8683298050514*fL[19]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*(alphaR[4]*fAvg[11]+alphaR[1]*fAvg[7])-21.21320343559643*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[21]-94.8683298050514*fL[21]-73.48469228349535*(fR[6]+fL[6])+42.42640687119286*fR[3]-42.42640687119286*fL[3])*amax-18.97366596101028*(alphaR[4]*fAvg[12]+alphaR[2]*fAvg[8])-21.21320343559643*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.04166666666666666*((18.97366596101028*fR[25]-18.97366596101028*fL[25]-14.69693845669907*(fR[8]+fL[8])+8.485281374238571*fR[4]-8.485281374238571*fL[4])*amax-4.242640687119286*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6]+alphaR[1]*fAvg[5]+alphaR[0]*fAvg[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[32]-94.86832980505142*fL[32]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[7]-42.42640687119286*fL[7])*amax-18.97366596101028*(alphaR[2]*fAvg[12]+alphaR[1]*fAvg[11]+alphaR[4]*(fAvg[8]+fAvg[7]))-21.21320343559643*(alphaR[0]*fAvg[4]+fAvg[0]*alphaR[4]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[35]-94.86832980505142*fL[35]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*(alphaR[4]*fAvg[17]+alphaR[1]*fAvg[13])-21.21320343559643*(alphaR[2]*fAvg[10]+alphaR[4]*fAvg[6]+alphaR[0]*fAvg[5]+alphaR[1]*fAvg[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[37]-94.86832980505142*fL[37]-73.48469228349535*(fR[17]+fL[17])+42.42640687119286*fR[10]-42.42640687119286*fL[10])*amax-18.97366596101028*(alphaR[4]*fAvg[18]+alphaR[2]*fAvg[14])-21.21320343559643*(alphaR[1]*fAvg[10]+alphaR[0]*fAvg[6]+alphaR[4]*fAvg[5]+alphaR[2]*fAvg[3])); 
  Ghat[7] = 0.008333333333333333*((73.48469228349536*(fR[20]+fL[20])-42.42640687119286*fR[12]+42.42640687119286*fL[12])*amax+21.21320343559643*alphaR[2]*fAvg[11]+21.21320343559643*alphaR[0]*fAvg[7]+18.97366596101028*(alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1])); 
  Ghat[8] = 0.008333333333333333*((73.48469228349536*(fR[23]+fL[23])-42.42640687119286*fR[13]+42.42640687119286*fL[13])*amax+21.21320343559643*alphaR[1]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[8]+18.97366596101028*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[28]+fL[28])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[4]*fAvg[19]+21.21320343559643*(alphaR[2]*fAvg[16]+alphaR[1]*fAvg[15])+21.21320343559643*alphaR[0]*fAvg[9]); 
  Ghat[10] = -0.008333333333333333*((94.8683298050514*fR[44]-94.8683298050514*fL[44]-73.48469228349535*(fR[31]+fL[31])+42.42640687119286*fR[18]-42.42640687119286*fL[18])*amax-18.97366596101028*(alphaR[2]*fAvg[18]+alphaR[1]*fAvg[17]+alphaR[4]*(fAvg[14]+fAvg[13]))-21.21320343559643*(alphaR[0]*fAvg[10]+alphaR[1]*fAvg[6]+alphaR[2]*fAvg[5]+fAvg[3]*alphaR[4])); 
  Ghat[11] = 0.008333333333333333*((73.48469228349536*(fR[33]+fL[33])-42.42640687119286*fR[22]+42.42640687119286*fL[22])*amax+16.97056274847715*alphaR[4]*fAvg[12]+21.21320343559643*alphaR[0]*fAvg[11]+21.21320343559643*alphaR[2]*fAvg[7]+18.97366596101028*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])); 
  Ghat[12] = 0.008333333333333333*((73.48469228349536*(fR[34]+fL[34])-42.42640687119286*fR[24]+42.42640687119286*fL[24])*amax+21.21320343559643*alphaR[0]*fAvg[12]+16.97056274847715*alphaR[4]*fAvg[11]+21.21320343559643*alphaR[1]*fAvg[8]+18.97366596101028*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])); 
  Ghat[13] = 0.008333333333333333*((73.48469228349536*(fR[36]+fL[36])-42.42640687119286*fR[26]+42.42640687119286*fL[26])*amax+21.21320343559643*alphaR[2]*fAvg[17]+21.21320343559643*alphaR[0]*fAvg[13]+18.97366596101028*(alphaR[4]*fAvg[10]+alphaR[1]*fAvg[5])); 
  Ghat[14] = 0.008333333333333333*((73.48469228349536*(fR[39]+fL[39])-42.42640687119286*fR[27]+42.42640687119286*fL[27])*amax+21.21320343559643*alphaR[1]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[14]+18.97366596101028*(alphaR[4]*fAvg[10]+alphaR[2]*fAvg[6])); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[29]+42.42640687119286*fL[29])*amax+21.21320343559643*alphaR[2]*fAvg[19]+21.21320343559643*(alphaR[4]*fAvg[16]+alphaR[0]*fAvg[15])+21.21320343559643*alphaR[1]*fAvg[9]); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[42]+fL[42])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+21.21320343559643*alphaR[1]*fAvg[19]+21.21320343559643*(alphaR[0]*fAvg[16]+alphaR[4]*fAvg[15])+21.21320343559643*alphaR[2]*fAvg[9]); 
  Ghat[17] = 0.008333333333333333*((73.48469228349536*(fR[45]+fL[45])-42.42640687119286*fR[38]+42.42640687119286*fL[38])*amax+16.97056274847715*alphaR[4]*fAvg[18]+21.21320343559643*alphaR[0]*fAvg[17]+21.21320343559643*alphaR[2]*fAvg[13]+18.97366596101028*(alphaR[1]*fAvg[10]+alphaR[4]*fAvg[5])); 
  Ghat[18] = 0.008333333333333333*((73.48469228349536*(fR[46]+fL[46])-42.42640687119286*fR[40]+42.42640687119286*fL[40])*amax+21.21320343559643*alphaR[0]*fAvg[18]+16.97056274847715*alphaR[4]*fAvg[17]+21.21320343559643*alphaR[1]*fAvg[14]+18.97366596101028*(alphaR[2]*fAvg[10]+alphaR[4]*fAvg[6])); 
  Ghat[19] = 0.008333333333333333*((73.48469228349536*(fR[47]+fL[47])-42.42640687119286*fR[43]+42.42640687119286*fL[43])*amax+21.21320343559643*alphaR[0]*fAvg[19]+21.21320343559643*(alphaR[1]*fAvg[16]+alphaR[2]*fAvg[15])+21.21320343559643*alphaR[4]*fAvg[9]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = -1.224744871391589*Ghat[0]; 
  incr[2] = 0.7071067811865475*Ghat[1]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = -1.224744871391589*Ghat[2]; 
  incr[7] = 0.7071067811865475*Ghat[4]; 
  incr[8] = -1.224744871391589*Ghat[3]; 
  incr[9] = 0.7071067811865475*Ghat[5]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 1.58113883008419*Ghat[0]; 
  incr[12] = 0.7071067811865475*Ghat[7]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = -1.224744871391589*Ghat[6]; 
  incr[18] = 0.7071067811865475*Ghat[10]; 
  incr[19] = 1.58113883008419*Ghat[1]; 
  incr[20] = -1.224744871391589*Ghat[7]; 
  incr[21] = 1.58113883008419*Ghat[2]; 
  incr[22] = 0.7071067811865475*Ghat[11]; 
  incr[23] = -1.224744871391589*Ghat[8]; 
  incr[24] = 0.7071067811865475*Ghat[12]; 
  incr[25] = 1.58113883008419*Ghat[3]; 
  incr[26] = 0.7071067811865475*Ghat[13]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = -1.224744871391589*Ghat[9]; 
  incr[29] = 0.7071067811865475*Ghat[15]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = 1.58113883008419*Ghat[4]; 
  incr[33] = -1.224744871391589*Ghat[11]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = 1.58113883008419*Ghat[5]; 
  incr[36] = -1.224744871391589*Ghat[13]; 
  incr[37] = 1.58113883008419*Ghat[6]; 
  incr[38] = 0.7071067811865475*Ghat[17]; 
  incr[39] = -1.224744871391589*Ghat[14]; 
  incr[40] = 0.7071067811865475*Ghat[18]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = -1.224744871391589*Ghat[16]; 
  incr[43] = 0.7071067811865475*Ghat[19]; 
  incr[44] = 1.58113883008419*Ghat[10]; 
  incr[45] = -1.224744871391589*Ghat[17]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

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
  outR[16] += incr[16]*rdx2R; 
  outR[17] += incr[17]*rdx2R; 
  outR[18] += incr[18]*rdx2R; 
  outR[19] += incr[19]*rdx2R; 
  outR[20] += incr[20]*rdx2R; 
  outR[21] += incr[21]*rdx2R; 
  outR[22] += incr[22]*rdx2R; 
  outR[23] += incr[23]*rdx2R; 
  outR[24] += incr[24]*rdx2R; 
  outR[25] += incr[25]*rdx2R; 
  outR[26] += incr[26]*rdx2R; 
  outR[27] += incr[27]*rdx2R; 
  outR[28] += incr[28]*rdx2R; 
  outR[29] += incr[29]*rdx2R; 
  outR[30] += incr[30]*rdx2R; 
  outR[31] += incr[31]*rdx2R; 
  outR[32] += incr[32]*rdx2R; 
  outR[33] += incr[33]*rdx2R; 
  outR[34] += incr[34]*rdx2R; 
  outR[35] += incr[35]*rdx2R; 
  outR[36] += incr[36]*rdx2R; 
  outR[37] += incr[37]*rdx2R; 
  outR[38] += incr[38]*rdx2R; 
  outR[39] += incr[39]*rdx2R; 
  outR[40] += incr[40]*rdx2R; 
  outR[41] += incr[41]*rdx2R; 
  outR[42] += incr[42]*rdx2R; 
  outR[43] += incr[43]*rdx2R; 
  outR[44] += incr[44]*rdx2R; 
  outR[45] += incr[45]*rdx2R; 
  outR[46] += incr[46]*rdx2R; 
  outR[47] += incr[47]*rdx2R; 

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
  outL[11] += -1.0*incr[11]*rdx2L; 
  outL[12] += -1.0*incr[12]*rdx2L; 
  outL[13] += -1.0*incr[13]*rdx2L; 
  outL[14] += -1.0*incr[14]*rdx2L; 
  outL[15] += incr[15]*rdx2L; 
  outL[16] += incr[16]*rdx2L; 
  outL[17] += incr[17]*rdx2L; 
  outL[18] += -1.0*incr[18]*rdx2L; 
  outL[19] += -1.0*incr[19]*rdx2L; 
  outL[20] += incr[20]*rdx2L; 
  outL[21] += -1.0*incr[21]*rdx2L; 
  outL[22] += -1.0*incr[22]*rdx2L; 
  outL[23] += incr[23]*rdx2L; 
  outL[24] += -1.0*incr[24]*rdx2L; 
  outL[25] += -1.0*incr[25]*rdx2L; 
  outL[26] += -1.0*incr[26]*rdx2L; 
  outL[27] += -1.0*incr[27]*rdx2L; 
  outL[28] += incr[28]*rdx2L; 
  outL[29] += -1.0*incr[29]*rdx2L; 
  outL[30] += -1.0*incr[30]*rdx2L; 
  outL[31] += incr[31]*rdx2L; 
  outL[32] += -1.0*incr[32]*rdx2L; 
  outL[33] += incr[33]*rdx2L; 
  outL[34] += incr[34]*rdx2L; 
  outL[35] += -1.0*incr[35]*rdx2L; 
  outL[36] += incr[36]*rdx2L; 
  outL[37] += -1.0*incr[37]*rdx2L; 
  outL[38] += -1.0*incr[38]*rdx2L; 
  outL[39] += incr[39]*rdx2L; 
  outL[40] += -1.0*incr[40]*rdx2L; 
  outL[41] += incr[41]*rdx2L; 
  outL[42] += incr[42]*rdx2L; 
  outL[43] += -1.0*incr[43]*rdx2L; 
  outL[44] += -1.0*incr[44]*rdx2L; 
  outL[45] += incr[45]*rdx2L; 
  outL[46] += incr[46]*rdx2L; 
  outL[47] += incr[47]*rdx2L; 

  return std::abs(alphaSurfAvgR);
} 
double EmGyrokineticGenGeoSurf2x2vSer_yL_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilL[48]; 
  hamilL[0] = (0.6666666666666666*(3.0*rdvpar2SqL*(m_*wvparSqL+bmag[0]*wmuL+phi[0]*q_)+m_))/rdvpar2SqL; 
  hamilL[1] = 2.0*(bmag[1]*wmuL+phi[1]*q_); 
  hamilL[2] = 2.0*phi[2]*q_; 
  hamilL[3] = (2.309401076758503*m_*wvparL)/rdvpar2L; 
  hamilL[4] = (1.154700538379252*bmag[0])/rdmu2L; 
  hamilL[5] = 2.0*phi[3]*q_; 
  hamilL[8] = (1.154700538379252*bmag[1])/rdmu2L; 
  hamilL[11] = 2.0*(bmag[4]*wmuL+phi[4]*q_); 
  hamilL[12] = 2.0*phi[5]*q_; 
  hamilL[13] = (0.5962847939999438*m_)/rdvpar2SqL; 
  hamilL[19] = 2.0*phi[6]*q_; 
  hamilL[20] = 2.0*phi[7]*q_; 
  hamilL[25] = (1.154700538379251*bmag[4])/rdmu2L; 

  double BstarYdBmagL[48]; 
  BstarYdBmagL[0] = -(0.8660254037844386*rdx2L*(2.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*wvparL+(3.0*(Apar[1]*b_z[4]+b_z[1]*Apar[4])*jacobTotInv[4]+(4.0*jacobTotInv[1]*Apar[4]+2.23606797749979*(Apar[0]*jacobTotInv[1]+jacobTotInv[0]*Apar[1]))*b_z[4]+2.23606797749979*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[4]+2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*q_))/q_; 
  BstarYdBmagL[1] = -(0.02474358296526967*rdx2L*(14.0*(b_z[4]*(10.0*jacobTotInv[4]+11.18033988749895*jacobTotInv[0])+5.0*b_z[1]*jacobTotInv[1])*m_*wvparL+(2.0*((122.9837387624885*Apar[4]+35.0*Apar[0])*b_z[4]+7.0*(5.0*b_z[0]*Apar[4]+4.47213595499958*Apar[1]*b_z[1]))*jacobTotInv[4]+7.0*((20.0*jacobTotInv[0]*Apar[4]+2.23606797749979*(11.0*Apar[1]*jacobTotInv[1]+5.0*Apar[0]*jacobTotInv[0]))*b_z[4]+2.23606797749979*(11.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[4]+5.0*((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])))*q_))/q_; 
  BstarYdBmagL[2] = -0.05773502691896259*(6.708203932499369*(2.23606797749979*(3.0*b_z[1]*jacobTotInv[4]+4.0*jacobTotInv[1]*b_z[4])+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[6]+3.0*(b_z[4]*(15.0*Apar[3]*jacobTotInv[4]+11.18033988749895*(jacobTotInv[0]*Apar[3]+jacobTotInv[1]*Apar[2]))+5.0*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])))*rdx2L; 
  BstarYdBmagL[3] = -(1.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2L)/(q_*rdvpar2L); 
  BstarYdBmagL[5] = -0.01428571428571429*(3.872983346207417*(2.0*(55.0*b_z[4]+15.65247584249853*b_z[0])*jacobTotInv[4]+7.0*(8.94427190999916*jacobTotInv[0]*b_z[4]+11.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0]))*Apar[6]+12.12435565298214*(2.0*(5.0*Apar[2]*b_z[4]+4.47213595499958*b_z[1]*Apar[3])*jacobTotInv[4]+2.23606797749979*(11.0*jacobTotInv[1]*Apar[3]+5.0*jacobTotInv[0]*Apar[2])*b_z[4]+5.0*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])))*rdx2L; 
  BstarYdBmagL[6] = -(1.0*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2L)/(q_*rdvpar2L); 
  BstarYdBmagL[11] = -(0.02474358296526967*rdx2L*(70.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*wvparL+((145.3444185374863*(Apar[1]*b_z[4]+b_z[1]*Apar[4])+35.0*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*jacobTotInv[4]+(245.9674775249769*jacobTotInv[1]*Apar[4]+35.0*(2.0*Apar[0]*jacobTotInv[1]+3.0*jacobTotInv[0]*Apar[1]))*b_z[4]+7.0*(5.0*(2.0*b_z[0]*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1])*Apar[4]+8.94427190999916*Apar[1]*b_z[1]*jacobTotInv[1]))*q_))/q_; 
  BstarYdBmagL[12] = -0.1*(3.872983346207417*(b_z[4]*(6.708203932499369*jacobTotInv[4]+5.0*jacobTotInv[0])+2.23606797749979*(2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*Apar[7]+1.732050807568877*(11.18033988749895*jacobTotInv[1]*b_z[4]+5.0*jacobTotInv[0]*b_z[1])*Apar[5])*rdx2L; 
  BstarYdBmagL[19] = -0.01844277783908294*(6.708203932499369*(2.23606797749979*(13.0*b_z[1]*jacobTotInv[4]+22.0*jacobTotInv[1]*b_z[4])+7.0*(2.0*b_z[0]*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1]))*Apar[6]+3.0*((65.0*Apar[3]*b_z[4]+15.65247584249853*(b_z[0]*Apar[3]+b_z[1]*Apar[2]))*jacobTotInv[4]+7.0*(2.23606797749979*(3.0*jacobTotInv[0]*Apar[3]+2.0*jacobTotInv[1]*Apar[2])*b_z[4]+4.0*b_z[1]*jacobTotInv[1]*Apar[3])))*rdx2L; 
  BstarYdBmagL[20] = -0.1*(1.732050807568877*(2.23606797749979*(4.0*b_z[1]*jacobTotInv[4]+11.0*jacobTotInv[1]*b_z[4])+5.0*(b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1]))*Apar[7]+3.872983346207417*(b_z[4]*(4.47213595499958*jacobTotInv[4]+5.0*jacobTotInv[0])+2.23606797749979*b_z[1]*jacobTotInv[1])*Apar[5])*rdx2L; 
  BstarYdBmagL[21] = -(1.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2L)/(q_*rdvpar2L); 

  double alphaL[20]; 
  alphaL[0] = (0.03535533905932736*((19.36491673103708*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamilL[20]+(30.0*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+33.54101966249684*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilL[19]+(17.32050807568877*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+19.36491673103709*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilL[11]+(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*(15.0*hamilL[5]+8.660254037844386*hamilL[1]))*m_*rdx2L+(19.36491673103709*BstarYdBmagL[3]*hamilL[13]+hamilL[3]*(19.36491673103709*BstarYdBmagL[12]+15.0*BstarYdBmagL[2]+8.660254037844386*BstarYdBmagL[0]))*q_*rdvpar2L))/(m_*q_); 
  alphaL[1] = (0.005050762722761052*(((121.2435565298214*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+135.5544171172596*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilL[20]+((368.9512162874653*b_z[4]+210.0*b_z[0])*jacobTotInv[4]+210.0*jacobTotInv[0]*b_z[4]+422.6168477474601*b_z[1]*jacobTotInv[1]+234.7871376374779*b_z[0]*jacobTotInv[0])*hamilL[19]+((213.014084041408*b_z[4]+121.2435565298214*b_z[0])*jacobTotInv[4]+121.2435565298214*jacobTotInv[0]*b_z[4]+243.9979508110672*b_z[1]*jacobTotInv[1]+135.5544171172596*b_z[0]*jacobTotInv[0])*hamilL[11]+(93.91485505499116*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+105.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilL[5]+hamilL[1]*(54.22176684690384*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+60.6217782649107*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])))*m_*rdx2L+(135.5544171172596*hamilL[3]*BstarYdBmagL[20]+135.5544171172596*BstarYdBmagL[6]*hamilL[13]+hamilL[3]*(105.0*BstarYdBmagL[5]+60.6217782649107*BstarYdBmagL[1]))*q_*rdvpar2L))/(m_*q_); 
  alphaL[2] = (0.1767766952966368*((8.660254037844386*BstarYdBmagL[12]+6.708203932499369*BstarYdBmagL[2]+3.872983346207417*BstarYdBmagL[0])*hamilL[13]+1.732050807568877*BstarYdBmagL[3]*hamilL[3])*rdvpar2L)/m_; 
  alphaL[3] = (0.03535533905932736*((17.32050807568877*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+19.36491673103708*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilL[25]+8.660254037844386*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamilL[8])*rdx2L)/q_; 
  alphaL[4] = (0.1767766952966368*(hamilL[13]*(8.660254037844387*BstarYdBmagL[20]+6.708203932499369*BstarYdBmagL[5]+3.872983346207417*BstarYdBmagL[1])+1.732050807568877*hamilL[3]*BstarYdBmagL[6])*rdvpar2L)/m_; 
  alphaL[5] = (0.005050762722761052*(((213.0140840414079*b_z[4]+121.2435565298214*b_z[0])*jacobTotInv[4]+121.2435565298214*jacobTotInv[0]*b_z[4]+243.9979508110673*b_z[1]*jacobTotInv[1]+135.5544171172596*b_z[0]*jacobTotInv[0])*hamilL[25]+(54.22176684690384*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+60.6217782649107*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilL[8])*rdx2L)/q_; 
  alphaL[7] = (0.005050762722761052*((((86.60254037844389*b_z[4]+135.5544171172596*b_z[0])*jacobTotInv[4]+135.5544171172596*jacobTotInv[0]*b_z[4]+121.2435565298214*b_z[1]*jacobTotInv[1])*hamilL[20]+(368.9512162874653*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+210.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilL[19]+(213.014084041408*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+121.2435565298214*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilL[11]+((67.0820393249937*b_z[4]+105.0*b_z[0])*jacobTotInv[4]+105.0*jacobTotInv[0]*b_z[4]+93.91485505499116*b_z[1]*jacobTotInv[1])*hamilL[5]+hamilL[1]*((38.72983346207418*b_z[4]+60.6217782649107*b_z[0])*jacobTotInv[4]+60.6217782649107*jacobTotInv[0]*b_z[4]+54.22176684690384*b_z[1]*jacobTotInv[1]))*m_*rdx2L+(135.5544171172596*hamilL[13]*BstarYdBmagL[21]+hamilL[3]*(105.0*BstarYdBmagL[19]+60.6217782649107*BstarYdBmagL[11]))*q_*rdvpar2L))/(m_*q_); 
  alphaL[8] = (0.6123724356957944*BstarYdBmagL[3]*hamilL[13]*rdvpar2L)/m_; 
  alphaL[11] = (0.04564354645876382*(6.708203932499369*hamilL[3]*BstarYdBmagL[21]+hamilL[13]*(25.98076211353316*BstarYdBmagL[19]+15.0*BstarYdBmagL[11]))*rdvpar2L)/m_; 
  alphaL[12] = (0.6123724356957944*BstarYdBmagL[6]*hamilL[13]*rdvpar2L)/m_; 
  alphaL[13] = (0.005050762722761052*((213.014084041408*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+121.2435565298214*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilL[25]+((38.72983346207417*b_z[4]+60.62177826491071*b_z[0])*jacobTotInv[4]+60.62177826491071*jacobTotInv[0]*b_z[4]+54.22176684690384*b_z[1]*jacobTotInv[1])*hamilL[8])*rdx2L)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgL = (0.00625*(((19.36491673103708*b_z[4]*jacobTotInv[4]+19.36491673103708*b_z[1]*jacobTotInv[1]+19.36491673103708*b_z[0]*jacobTotInv[0])*hamilL[20]+(30.0*b_z[1]*jacobTotInv[4]+30.0*jacobTotInv[1]*b_z[4]+33.54101966249684*b_z[0]*jacobTotInv[1]+33.54101966249684*jacobTotInv[0]*b_z[1])*hamilL[19]+(17.32050807568877*b_z[1]*jacobTotInv[4]+17.32050807568877*jacobTotInv[1]*b_z[4]+19.36491673103709*b_z[0]*jacobTotInv[1]+19.36491673103709*jacobTotInv[0]*b_z[1])*hamilL[11]+(15.0*b_z[4]*jacobTotInv[4]+15.0*b_z[1]*jacobTotInv[1]+15.0*b_z[0]*jacobTotInv[0])*hamilL[5]+8.660254037844386*hamilL[1]*b_z[4]*jacobTotInv[4]+8.660254037844386*b_z[1]*hamilL[1]*jacobTotInv[1]+8.660254037844386*b_z[0]*jacobTotInv[0]*hamilL[1])*m_*rdx2L+(19.36491673103709*BstarYdBmagL[3]*hamilL[13]+19.36491673103709*hamilL[3]*BstarYdBmagL[12]+(15.0*BstarYdBmagL[2]+8.660254037844386*BstarYdBmagL[0])*hamilL[3])*q_*rdvpar2L))/(m_*q_); 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[12]+fL[12])+1.732050807568877*(fL[2]-1.0*fR[2])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[20]+fL[20])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[22]+fL[22])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[26]+fL[26])+3.0*(fL[9]-1.0*fR[9]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[33]+fL[33])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[6]+fL[6]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[36]+fL[36])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[38]+fL[38])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[19]-1.0*(8.660254037844387*fL[19]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[24]-1.0*(8.660254037844387*fL[24]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[29]-1.0*(8.660254037844387*fL[29]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[45]+fL[45])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[17]+fL[17])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[21]+fL[21]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[23]+fL[23]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[35]-1.0*(8.660254037844387*fL[35]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[40]-1.0*(8.660254037844387*fL[40]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[37]+fL[37]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[39]+fL[39]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[42]+fL[42]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[12]-6.324555320336761*fL[12]-4.898979485566357*(fR[2]+fL[2])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaL[13]*fAvg[13]+alphaL[12]*fAvg[12]+alphaL[11]*fAvg[11]+alphaL[8]*fAvg[8]+alphaL[7]*fAvg[7]+alphaL[5]*fAvg[5]+alphaL[4]*fAvg[4]+alphaL[3]*fAvg[3]+alphaL[2]*fAvg[2]+alphaL[1]*fAvg[1]+alphaL[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[20]-94.8683298050514*fL[20]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*(alphaL[5]*fAvg[13]+fAvg[5]*alphaL[13])-21.21320343559643*(alphaL[8]*fAvg[12]+fAvg[8]*alphaL[12])-18.97366596101028*(alphaL[4]*fAvg[11]+fAvg[4]*alphaL[11]+alphaL[1]*fAvg[7]+fAvg[1]*alphaL[7])-21.21320343559643*(alphaL[3]*fAvg[5]+fAvg[3]*alphaL[5]+alphaL[2]*fAvg[4]+fAvg[2]*alphaL[4]+alphaL[0]*fAvg[1]+fAvg[0]*alphaL[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[22]-94.8683298050514*fL[22]-73.48469228349535*(fR[7]+fL[7])+42.42640687119286*fR[3]-42.42640687119286*fL[3])*amax-21.21320343559643*alphaL[13]*fAvg[17]-18.97366596101028*(alphaL[4]*fAvg[12]+fAvg[4]*alphaL[12])-21.21320343559643*(alphaL[7]*fAvg[11]+fAvg[7]*alphaL[11])-21.21320343559643*alphaL[5]*fAvg[10]-18.97366596101028*(alphaL[2]*fAvg[8]+fAvg[2]*alphaL[8])-21.21320343559643*(alphaL[3]*fAvg[6]+alphaL[1]*fAvg[4]+fAvg[1]*alphaL[4]+alphaL[0]*fAvg[2]+fAvg[0]*alphaL[2])); 
  Ghat[3] = -0.008333333333333333*((94.8683298050514*fR[26]-94.8683298050514*fL[26]-73.48469228349535*(fR[9]+fL[9])+42.42640687119286*fR[4]-42.42640687119286*fL[4])*amax-21.21320343559643*(alphaL[12]*fAvg[18]+alphaL[11]*fAvg[17])-18.97366596101028*alphaL[5]*fAvg[15]-21.21320343559643*(alphaL[8]*fAvg[14]+alphaL[7]*fAvg[13]+fAvg[7]*alphaL[13])-21.21320343559643*alphaL[4]*fAvg[10]-18.97366596101028*alphaL[3]*fAvg[9]-21.21320343559643*(alphaL[2]*fAvg[6]+alphaL[1]*fAvg[5]+fAvg[1]*alphaL[5]+alphaL[0]*fAvg[3]+fAvg[0]*alphaL[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[33]-94.86832980505142*fL[33]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[6]-42.42640687119286*fL[6])*amax-18.97366596101028*(alphaL[5]*fAvg[17]+fAvg[10]*alphaL[13])+((-16.97056274847715*alphaL[11])-18.97366596101028*alphaL[2])*fAvg[12]+((-16.97056274847715*fAvg[11])-18.97366596101028*fAvg[2])*alphaL[12]-18.97366596101028*(alphaL[1]*fAvg[11]+fAvg[1]*alphaL[11])-21.21320343559643*alphaL[3]*fAvg[10]-18.97366596101028*(alphaL[4]*fAvg[8]+fAvg[4]*alphaL[8]+alphaL[4]*fAvg[7]+fAvg[4]*alphaL[7])-21.21320343559643*(alphaL[5]*fAvg[6]+alphaL[0]*fAvg[4]+fAvg[0]*alphaL[4]+alphaL[1]*fAvg[2]+fAvg[1]*alphaL[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[36]-94.86832980505142*fL[36]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-21.21320343559643*alphaL[8]*fAvg[18]-18.97366596101028*alphaL[4]*fAvg[17]+((-16.97056274847715*alphaL[13])-18.97366596101028*alphaL[3])*fAvg[15]-21.21320343559643*alphaL[12]*fAvg[14]-18.97366596101028*(alphaL[1]*fAvg[13]+fAvg[1]*alphaL[13])+fAvg[10]*((-18.97366596101028*alphaL[11])-21.21320343559643*alphaL[2])-18.97366596101028*(alphaL[5]*(fAvg[9]+fAvg[7])+fAvg[5]*alphaL[7])-21.21320343559643*(alphaL[4]*fAvg[6]+alphaL[0]*fAvg[5]+fAvg[0]*alphaL[5]+alphaL[1]*fAvg[3]+fAvg[1]*alphaL[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[38]-94.86832980505142*fL[38]-73.48469228349535*(fR[18]+fL[18])+42.42640687119286*fR[10]-42.42640687119286*fL[10])*amax-18.97366596101028*(alphaL[5]*fAvg[19]+alphaL[4]*fAvg[18])-21.21320343559643*alphaL[7]*fAvg[17]-18.97366596101028*(alphaL[3]*fAvg[16]+alphaL[2]*fAvg[14])-21.21320343559643*(alphaL[11]*fAvg[13]+fAvg[11]*alphaL[13])+fAvg[10]*((-18.97366596101028*alphaL[12])-21.21320343559643*alphaL[1])-18.97366596101028*fAvg[6]*alphaL[8]-21.21320343559643*(alphaL[0]*fAvg[6]+alphaL[4]*fAvg[5]+fAvg[4]*alphaL[5]+alphaL[2]*fAvg[3]+fAvg[2]*alphaL[3])); 
  Ghat[7] = 0.001190476190476191*((514.3928459844675*(fR[19]+fL[19])-296.98484809835*fR[11]+296.98484809835*fL[11])*amax+(94.86832980505142*alphaL[13]+148.492424049175*alphaL[3])*fAvg[13]+148.492424049175*fAvg[3]*alphaL[13]+132.815661727072*alphaL[12]*fAvg[12]+(94.86832980505142*alphaL[11]+148.492424049175*alphaL[2])*fAvg[11]+148.492424049175*fAvg[2]*alphaL[11]+(94.86832980505142*alphaL[7]+148.492424049175*alphaL[0])*fAvg[7]+148.492424049175*fAvg[0]*alphaL[7]+132.815661727072*(alphaL[5]*fAvg[5]+alphaL[4]*fAvg[4]+alphaL[1]*fAvg[1])); 
  Ghat[8] = 0.001190476190476191*((514.3928459844675*(fR[24]+fL[24])-296.98484809835*fR[13]+296.98484809835*fL[13])*amax+148.492424049175*(alphaL[5]*fAvg[18]+alphaL[3]*fAvg[14])+(94.86832980505142*alphaL[12]+148.492424049175*alphaL[1])*fAvg[12]+148.492424049175*fAvg[1]*alphaL[12]+132.815661727072*alphaL[11]*fAvg[11]+(94.86832980505142*alphaL[8]+148.492424049175*alphaL[0])*fAvg[8]+148.492424049175*fAvg[0]*alphaL[8]+132.815661727072*(alphaL[4]*fAvg[4]+alphaL[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[29]+fL[29])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaL[4]*fAvg[19]+21.21320343559643*(alphaL[2]*fAvg[16]+alphaL[1]*fAvg[15])+18.97366596101028*alphaL[13]*fAvg[13]+21.21320343559643*alphaL[0]*fAvg[9]+18.97366596101028*(alphaL[5]*fAvg[5]+alphaL[3]*fAvg[3])); 
  Ghat[10] = -0.001666666666666667*((474.341649025257*fR[45]-474.341649025257*fL[45]-367.4234614174767*(fR[31]+fL[31])+212.1320343559643*fR[17]-212.1320343559643*fL[17])*amax+((-84.85281374238573*alphaL[13])-94.86832980505142*alphaL[3])*fAvg[19]+((-84.85281374238573*alphaL[11])-94.86832980505142*alphaL[2])*fAvg[18]+((-84.85281374238573*alphaL[12])-94.86832980505142*alphaL[1])*fAvg[17]-94.8683298050514*(alphaL[5]*fAvg[16]+alphaL[4]*(fAvg[14]+fAvg[13])+fAvg[4]*alphaL[13]+fAvg[6]*alphaL[12]+alphaL[5]*fAvg[11]+fAvg[5]*alphaL[11])+((-94.86832980505142*(alphaL[8]+alphaL[7]))-106.0660171779821*alphaL[0])*fAvg[10]-106.0660171779821*(alphaL[1]*fAvg[6]+alphaL[2]*fAvg[5]+fAvg[2]*alphaL[5]+alphaL[3]*fAvg[4]+fAvg[3]*alphaL[4])); 
  Ghat[11] = 0.001190476190476191*((514.3928459844675*(fR[32]+fL[32])-296.98484809835*fR[21]+296.98484809835*fL[21])*amax+(94.86832980505142*alphaL[13]+148.492424049175*alphaL[3])*fAvg[17]+148.492424049175*fAvg[6]*alphaL[13]+118.79393923934*(alphaL[4]*fAvg[12]+fAvg[4]*alphaL[12])+(132.815661727072*alphaL[8]+94.86832980505142*alphaL[7]+148.492424049175*alphaL[0])*fAvg[11]+(132.815661727072*fAvg[8]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaL[11]+132.815661727072*alphaL[5]*fAvg[10]+148.492424049175*(alphaL[2]*fAvg[7]+fAvg[2]*alphaL[7])+132.815661727072*(alphaL[1]*fAvg[4]+fAvg[1]*alphaL[4])); 
  Ghat[12] = 0.001190476190476191*((514.3928459844675*(fR[34]+fL[34])-296.98484809835*fR[23]+296.98484809835*fL[23])*amax+(132.815661727072*alphaL[13]+148.492424049175*alphaL[3])*fAvg[18]+148.492424049175*alphaL[5]*fAvg[14]+(94.86832980505142*alphaL[8]+132.815661727072*alphaL[7]+148.492424049175*alphaL[0])*fAvg[12]+(94.86832980505142*fAvg[8]+132.815661727072*fAvg[7]+148.492424049175*fAvg[0])*alphaL[12]+118.79393923934*(alphaL[4]*fAvg[11]+fAvg[4]*alphaL[11])+148.492424049175*(alphaL[1]*fAvg[8]+fAvg[1]*alphaL[8])+132.815661727072*(alphaL[2]*fAvg[4]+fAvg[2]*alphaL[4])); 
  Ghat[13] = 0.001190476190476191*((514.3928459844675*(fR[35]+fL[35])-296.98484809835*fR[25]+296.98484809835*fL[25])*amax+132.815661727072*alphaL[12]*fAvg[18]+(94.86832980505142*alphaL[11]+148.492424049175*alphaL[2])*fAvg[17]+118.79393923934*alphaL[5]*fAvg[15]+(94.86832980505142*alphaL[7]+148.492424049175*alphaL[0])*fAvg[13]+(132.815661727072*fAvg[9]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaL[13]+148.492424049175*fAvg[6]*alphaL[11]+132.815661727072*alphaL[4]*fAvg[10]+148.492424049175*(alphaL[3]*fAvg[7]+fAvg[3]*alphaL[7])+132.815661727072*(alphaL[1]*fAvg[5]+fAvg[1]*alphaL[5])); 
  Ghat[14] = 0.001190476190476191*((514.3928459844675*(fR[40]+fL[40])-296.98484809835*fR[27]+296.98484809835*fL[27])*amax+(94.86832980505142*alphaL[12]+148.492424049175*alphaL[1])*fAvg[18]+132.815661727072*alphaL[11]*fAvg[17]+(94.86832980505142*alphaL[8]+148.492424049175*alphaL[0])*fAvg[14]+148.492424049175*(alphaL[5]*fAvg[12]+fAvg[5]*alphaL[12])+132.815661727072*alphaL[4]*fAvg[10]+148.492424049175*(alphaL[3]*fAvg[8]+fAvg[3]*alphaL[8])+132.815661727072*alphaL[2]*fAvg[6]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+(18.97366596101028*alphaL[11]+21.21320343559643*alphaL[2])*fAvg[19]+21.21320343559643*alphaL[4]*fAvg[16]+(18.97366596101028*alphaL[7]+21.21320343559643*alphaL[0])*fAvg[15]+16.97056274847715*(alphaL[5]*fAvg[13]+fAvg[5]*alphaL[13])+21.21320343559643*alphaL[1]*fAvg[9]+18.97366596101028*(alphaL[3]*fAvg[5]+fAvg[3]*alphaL[5])); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+(18.97366596101028*alphaL[12]+21.21320343559643*alphaL[1])*fAvg[19]+18.97366596101028*alphaL[13]*fAvg[17]+(18.97366596101028*alphaL[8]+21.21320343559643*alphaL[0])*fAvg[16]+21.21320343559643*alphaL[4]*fAvg[15]+18.97366596101028*alphaL[5]*fAvg[10]+21.21320343559643*alphaL[2]*fAvg[9]+18.97366596101028*alphaL[3]*fAvg[6]); 
  Ghat[17] = 2.380952380952381e-4*((2571.964229922338*(fR[44]+fL[44])-1484.92424049175*fR[37]+1484.92424049175*fL[37])*amax+593.9696961967002*(alphaL[5]*fAvg[19]+alphaL[4]*fAvg[18])+(664.0783086353599*alphaL[8]+474.3416490252571*alphaL[7]+742.462120245875*alphaL[0])*fAvg[17]+664.0783086353599*(alphaL[13]*fAvg[16]+alphaL[11]*fAvg[14])+(474.3416490252571*alphaL[11]+742.462120245875*alphaL[2])*fAvg[13]+(474.3416490252571*fAvg[11]+742.462120245875*fAvg[2])*alphaL[13]+593.9696961967002*fAvg[10]*alphaL[12]+742.462120245875*(alphaL[3]*fAvg[11]+fAvg[3]*alphaL[11])+664.0783086353599*alphaL[1]*fAvg[10]+742.462120245875*fAvg[6]*alphaL[7]+664.0783086353599*(alphaL[4]*fAvg[5]+fAvg[4]*alphaL[5])); 
  Ghat[18] = 2.380952380952381e-4*((2571.964229922338*(fR[46]+fL[46])-1484.92424049175*fR[39]+1484.92424049175*fL[39])*amax+(474.3416490252571*alphaL[8]+664.0783086353599*alphaL[7]+742.462120245875*alphaL[0])*fAvg[18]+593.9696961967002*alphaL[4]*fAvg[17]+(474.3416490252571*alphaL[12]+742.462120245875*alphaL[1])*fAvg[14]+664.0783086353599*(alphaL[12]*fAvg[13]+fAvg[12]*alphaL[13])+742.462120245875*(alphaL[3]*fAvg[12]+fAvg[3]*alphaL[12])+fAvg[10]*(593.9696961967002*alphaL[11]+664.0783086353599*alphaL[2])+742.462120245875*(alphaL[5]*fAvg[8]+fAvg[5]*alphaL[8])+664.0783086353599*alphaL[4]*fAvg[6]); 
  Ghat[19] = 0.001666666666666667*((367.4234614174769*(fR[47]+fL[47])-212.1320343559643*fR[42]+212.1320343559643*fL[42])*amax+(94.86832980505142*(alphaL[8]+alphaL[7])+106.0660171779821*alphaL[0])*fAvg[19]+84.85281374238573*alphaL[5]*fAvg[17]+(94.86832980505142*alphaL[12]+106.0660171779822*alphaL[1])*fAvg[16]+(94.86832980505142*alphaL[11]+106.0660171779822*alphaL[2])*fAvg[15]+fAvg[10]*(84.85281374238573*alphaL[13]+94.86832980505142*alphaL[3])+106.0660171779821*alphaL[4]*fAvg[9]+94.86832980505142*alphaL[5]*fAvg[6]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = -1.224744871391589*Ghat[0]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = 0.7071067811865475*Ghat[4]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = -1.224744871391589*Ghat[3]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 1.58113883008419*Ghat[0]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = 0.7071067811865475*Ghat[10]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = -1.224744871391589*Ghat[7]; 
  incr[20] = 1.58113883008419*Ghat[1]; 
  incr[21] = 0.7071067811865475*Ghat[11]; 
  incr[22] = 1.58113883008419*Ghat[2]; 
  incr[23] = 0.7071067811865475*Ghat[12]; 
  incr[24] = -1.224744871391589*Ghat[8]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 1.58113883008419*Ghat[3]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = -1.224744871391589*Ghat[9]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = 1.58113883008419*Ghat[4]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = -1.224744871391589*Ghat[13]; 
  incr[36] = 1.58113883008419*Ghat[5]; 
  incr[37] = 0.7071067811865475*Ghat[17]; 
  incr[38] = 1.58113883008419*Ghat[6]; 
  incr[39] = 0.7071067811865475*Ghat[18]; 
  incr[40] = -1.224744871391589*Ghat[14]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = 0.7071067811865475*Ghat[19]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = 1.58113883008419*Ghat[10]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

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
  outR[16] += incr[16]*rdy2R; 
  outR[17] += incr[17]*rdy2R; 
  outR[18] += incr[18]*rdy2R; 
  outR[19] += incr[19]*rdy2R; 
  outR[20] += incr[20]*rdy2R; 
  outR[21] += incr[21]*rdy2R; 
  outR[22] += incr[22]*rdy2R; 
  outR[23] += incr[23]*rdy2R; 
  outR[24] += incr[24]*rdy2R; 
  outR[25] += incr[25]*rdy2R; 
  outR[26] += incr[26]*rdy2R; 
  outR[27] += incr[27]*rdy2R; 
  outR[28] += incr[28]*rdy2R; 
  outR[29] += incr[29]*rdy2R; 
  outR[30] += incr[30]*rdy2R; 
  outR[31] += incr[31]*rdy2R; 
  outR[32] += incr[32]*rdy2R; 
  outR[33] += incr[33]*rdy2R; 
  outR[34] += incr[34]*rdy2R; 
  outR[35] += incr[35]*rdy2R; 
  outR[36] += incr[36]*rdy2R; 
  outR[37] += incr[37]*rdy2R; 
  outR[38] += incr[38]*rdy2R; 
  outR[39] += incr[39]*rdy2R; 
  outR[40] += incr[40]*rdy2R; 
  outR[41] += incr[41]*rdy2R; 
  outR[42] += incr[42]*rdy2R; 
  outR[43] += incr[43]*rdy2R; 
  outR[44] += incr[44]*rdy2R; 
  outR[45] += incr[45]*rdy2R; 
  outR[46] += incr[46]*rdy2R; 
  outR[47] += incr[47]*rdy2R; 

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
  outL[11] += -1.0*incr[11]*rdy2L; 
  outL[12] += -1.0*incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += -1.0*incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 
  outL[16] += incr[16]*rdy2L; 
  outL[17] += -1.0*incr[17]*rdy2L; 
  outL[18] += incr[18]*rdy2L; 
  outL[19] += incr[19]*rdy2L; 
  outL[20] += -1.0*incr[20]*rdy2L; 
  outL[21] += -1.0*incr[21]*rdy2L; 
  outL[22] += -1.0*incr[22]*rdy2L; 
  outL[23] += -1.0*incr[23]*rdy2L; 
  outL[24] += incr[24]*rdy2L; 
  outL[25] += -1.0*incr[25]*rdy2L; 
  outL[26] += -1.0*incr[26]*rdy2L; 
  outL[27] += -1.0*incr[27]*rdy2L; 
  outL[28] += -1.0*incr[28]*rdy2L; 
  outL[29] += incr[29]*rdy2L; 
  outL[30] += -1.0*incr[30]*rdy2L; 
  outL[31] += incr[31]*rdy2L; 
  outL[32] += incr[32]*rdy2L; 
  outL[33] += -1.0*incr[33]*rdy2L; 
  outL[34] += incr[34]*rdy2L; 
  outL[35] += incr[35]*rdy2L; 
  outL[36] += -1.0*incr[36]*rdy2L; 
  outL[37] += -1.0*incr[37]*rdy2L; 
  outL[38] += -1.0*incr[38]*rdy2L; 
  outL[39] += -1.0*incr[39]*rdy2L; 
  outL[40] += incr[40]*rdy2L; 
  outL[41] += incr[41]*rdy2L; 
  outL[42] += -1.0*incr[42]*rdy2L; 
  outL[43] += incr[43]*rdy2L; 
  outL[44] += incr[44]*rdy2L; 
  outL[45] += -1.0*incr[45]*rdy2L; 
  outL[46] += incr[46]*rdy2L; 
  outL[47] += incr[47]*rdy2L; 

  return std::abs(alphaSurfAvgL);
} 
double EmGyrokineticGenGeoSurf2x2vSer_yR_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[11] = 2.0*(bmag[4]*wmuR+phi[4]*q_); 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 
  hamilR[25] = (1.154700538379251*bmag[4])/rdmu2R; 

  double BstarYdBmagR[48]; 
  BstarYdBmagR[0] = -(0.8660254037844386*rdx2R*(2.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*wvparR+(3.0*(Apar[1]*b_z[4]+b_z[1]*Apar[4])*jacobTotInv[4]+(4.0*jacobTotInv[1]*Apar[4]+2.23606797749979*(Apar[0]*jacobTotInv[1]+jacobTotInv[0]*Apar[1]))*b_z[4]+2.23606797749979*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[4]+2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*q_))/q_; 
  BstarYdBmagR[1] = -(0.02474358296526967*rdx2R*(14.0*(b_z[4]*(10.0*jacobTotInv[4]+11.18033988749895*jacobTotInv[0])+5.0*b_z[1]*jacobTotInv[1])*m_*wvparR+(2.0*((122.9837387624885*Apar[4]+35.0*Apar[0])*b_z[4]+7.0*(5.0*b_z[0]*Apar[4]+4.47213595499958*Apar[1]*b_z[1]))*jacobTotInv[4]+7.0*((20.0*jacobTotInv[0]*Apar[4]+2.23606797749979*(11.0*Apar[1]*jacobTotInv[1]+5.0*Apar[0]*jacobTotInv[0]))*b_z[4]+2.23606797749979*(11.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[4]+5.0*((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])))*q_))/q_; 
  BstarYdBmagR[2] = -0.05773502691896259*(6.708203932499369*(2.23606797749979*(3.0*b_z[1]*jacobTotInv[4]+4.0*jacobTotInv[1]*b_z[4])+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[6]+3.0*(b_z[4]*(15.0*Apar[3]*jacobTotInv[4]+11.18033988749895*(jacobTotInv[0]*Apar[3]+jacobTotInv[1]*Apar[2]))+5.0*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])))*rdx2R; 
  BstarYdBmagR[3] = -(1.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[5] = -0.01428571428571429*(3.872983346207417*(2.0*(55.0*b_z[4]+15.65247584249853*b_z[0])*jacobTotInv[4]+7.0*(8.94427190999916*jacobTotInv[0]*b_z[4]+11.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0]))*Apar[6]+12.12435565298214*(2.0*(5.0*Apar[2]*b_z[4]+4.47213595499958*b_z[1]*Apar[3])*jacobTotInv[4]+2.23606797749979*(11.0*jacobTotInv[1]*Apar[3]+5.0*jacobTotInv[0]*Apar[2])*b_z[4]+5.0*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])))*rdx2R; 
  BstarYdBmagR[6] = -(1.0*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[11] = -(0.02474358296526967*rdx2R*(70.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*wvparR+((145.3444185374863*(Apar[1]*b_z[4]+b_z[1]*Apar[4])+35.0*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*jacobTotInv[4]+(245.9674775249769*jacobTotInv[1]*Apar[4]+35.0*(2.0*Apar[0]*jacobTotInv[1]+3.0*jacobTotInv[0]*Apar[1]))*b_z[4]+7.0*(5.0*(2.0*b_z[0]*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1])*Apar[4]+8.94427190999916*Apar[1]*b_z[1]*jacobTotInv[1]))*q_))/q_; 
  BstarYdBmagR[12] = -0.1*(3.872983346207417*(b_z[4]*(6.708203932499369*jacobTotInv[4]+5.0*jacobTotInv[0])+2.23606797749979*(2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*Apar[7]+1.732050807568877*(11.18033988749895*jacobTotInv[1]*b_z[4]+5.0*jacobTotInv[0]*b_z[1])*Apar[5])*rdx2R; 
  BstarYdBmagR[19] = -0.01844277783908294*(6.708203932499369*(2.23606797749979*(13.0*b_z[1]*jacobTotInv[4]+22.0*jacobTotInv[1]*b_z[4])+7.0*(2.0*b_z[0]*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1]))*Apar[6]+3.0*((65.0*Apar[3]*b_z[4]+15.65247584249853*(b_z[0]*Apar[3]+b_z[1]*Apar[2]))*jacobTotInv[4]+7.0*(2.23606797749979*(3.0*jacobTotInv[0]*Apar[3]+2.0*jacobTotInv[1]*Apar[2])*b_z[4]+4.0*b_z[1]*jacobTotInv[1]*Apar[3])))*rdx2R; 
  BstarYdBmagR[20] = -0.1*(1.732050807568877*(2.23606797749979*(4.0*b_z[1]*jacobTotInv[4]+11.0*jacobTotInv[1]*b_z[4])+5.0*(b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1]))*Apar[7]+3.872983346207417*(b_z[4]*(4.47213595499958*jacobTotInv[4]+5.0*jacobTotInv[0])+2.23606797749979*b_z[1]*jacobTotInv[1])*Apar[5])*rdx2R; 
  BstarYdBmagR[21] = -(1.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[20]; 
  alphaR[0] = (0.03535533905932736*((19.36491673103708*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamilR[20]+((-30.0*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-33.54101966249684*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[19]+(17.32050807568877*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+19.36491673103709*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[11]+(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*(8.660254037844386*hamilR[1]-15.0*hamilR[5]))*m_*rdx2R+(19.36491673103709*BstarYdBmagR[3]*hamilR[13]+hamilR[3]*(19.36491673103709*BstarYdBmagR[12]-15.0*BstarYdBmagR[2]+8.660254037844386*BstarYdBmagR[0]))*q_*rdvpar2R))/(m_*q_); 
  alphaR[1] = (0.005050762722761052*(((121.2435565298214*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+135.5544171172596*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[20]+(((-368.9512162874653*b_z[4])-210.0*b_z[0])*jacobTotInv[4]-210.0*jacobTotInv[0]*b_z[4]-422.6168477474601*b_z[1]*jacobTotInv[1]-234.7871376374779*b_z[0]*jacobTotInv[0])*hamilR[19]+((213.014084041408*b_z[4]+121.2435565298214*b_z[0])*jacobTotInv[4]+121.2435565298214*jacobTotInv[0]*b_z[4]+243.9979508110672*b_z[1]*jacobTotInv[1]+135.5544171172596*b_z[0]*jacobTotInv[0])*hamilR[11]+((-93.91485505499116*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-105.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[5]+hamilR[1]*(54.22176684690384*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+60.6217782649107*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])))*m_*rdx2R+(135.5544171172596*hamilR[3]*BstarYdBmagR[20]+135.5544171172596*BstarYdBmagR[6]*hamilR[13]+hamilR[3]*(60.6217782649107*BstarYdBmagR[1]-105.0*BstarYdBmagR[5]))*q_*rdvpar2R))/(m_*q_); 
  alphaR[2] = (0.1767766952966368*((8.660254037844386*BstarYdBmagR[12]-6.708203932499369*BstarYdBmagR[2]+3.872983346207417*BstarYdBmagR[0])*hamilR[13]+1.732050807568877*BstarYdBmagR[3]*hamilR[3])*rdvpar2R)/m_; 
  alphaR[3] = (0.03535533905932736*((17.32050807568877*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+19.36491673103708*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[25]+8.660254037844386*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*hamilR[8])*rdx2R)/q_; 
  alphaR[4] = (0.1767766952966368*(hamilR[13]*(8.660254037844387*BstarYdBmagR[20]-6.708203932499369*BstarYdBmagR[5]+3.872983346207417*BstarYdBmagR[1])+1.732050807568877*hamilR[3]*BstarYdBmagR[6])*rdvpar2R)/m_; 
  alphaR[5] = (0.005050762722761052*(((213.0140840414079*b_z[4]+121.2435565298214*b_z[0])*jacobTotInv[4]+121.2435565298214*jacobTotInv[0]*b_z[4]+243.9979508110673*b_z[1]*jacobTotInv[1]+135.5544171172596*b_z[0]*jacobTotInv[0])*hamilR[25]+(54.22176684690384*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+60.6217782649107*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[8])*rdx2R)/q_; 
  alphaR[7] = (0.005050762722761052*((((86.60254037844389*b_z[4]+135.5544171172596*b_z[0])*jacobTotInv[4]+135.5544171172596*jacobTotInv[0]*b_z[4]+121.2435565298214*b_z[1]*jacobTotInv[1])*hamilR[20]+((-368.9512162874653*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4]))-210.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[19]+(213.014084041408*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+121.2435565298214*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[11]+(((-67.0820393249937*b_z[4])-105.0*b_z[0])*jacobTotInv[4]-105.0*jacobTotInv[0]*b_z[4]-93.91485505499116*b_z[1]*jacobTotInv[1])*hamilR[5]+hamilR[1]*((38.72983346207418*b_z[4]+60.6217782649107*b_z[0])*jacobTotInv[4]+60.6217782649107*jacobTotInv[0]*b_z[4]+54.22176684690384*b_z[1]*jacobTotInv[1]))*m_*rdx2R+(135.5544171172596*hamilR[13]*BstarYdBmagR[21]+hamilR[3]*(60.6217782649107*BstarYdBmagR[11]-105.0*BstarYdBmagR[19]))*q_*rdvpar2R))/(m_*q_); 
  alphaR[8] = (0.6123724356957944*BstarYdBmagR[3]*hamilR[13]*rdvpar2R)/m_; 
  alphaR[11] = (0.04564354645876382*(6.708203932499369*hamilR[3]*BstarYdBmagR[21]+hamilR[13]*(15.0*BstarYdBmagR[11]-25.98076211353316*BstarYdBmagR[19]))*rdvpar2R)/m_; 
  alphaR[12] = (0.6123724356957944*BstarYdBmagR[6]*hamilR[13]*rdvpar2R)/m_; 
  alphaR[13] = (0.005050762722761052*((213.014084041408*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+121.2435565298214*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*hamilR[25]+((38.72983346207417*b_z[4]+60.62177826491071*b_z[0])*jacobTotInv[4]+60.62177826491071*jacobTotInv[0]*b_z[4]+54.22176684690384*b_z[1]*jacobTotInv[1])*hamilR[8])*rdx2R)/q_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.00625*(((19.36491673103708*b_z[4]*jacobTotInv[4]+19.36491673103708*b_z[1]*jacobTotInv[1]+19.36491673103708*b_z[0]*jacobTotInv[0])*hamilR[20]+((-30.0*b_z[1]*jacobTotInv[4])-30.0*jacobTotInv[1]*b_z[4]-33.54101966249684*b_z[0]*jacobTotInv[1]-33.54101966249684*jacobTotInv[0]*b_z[1])*hamilR[19]+(17.32050807568877*b_z[1]*jacobTotInv[4]+17.32050807568877*jacobTotInv[1]*b_z[4]+19.36491673103709*b_z[0]*jacobTotInv[1]+19.36491673103709*jacobTotInv[0]*b_z[1])*hamilR[11]+((-15.0*b_z[4]*jacobTotInv[4])-15.0*b_z[1]*jacobTotInv[1]-15.0*b_z[0]*jacobTotInv[0])*hamilR[5]+8.660254037844386*hamilR[1]*b_z[4]*jacobTotInv[4]+8.660254037844386*b_z[1]*hamilR[1]*jacobTotInv[1]+8.660254037844386*b_z[0]*jacobTotInv[0]*hamilR[1])*m_*rdx2R+(19.36491673103709*BstarYdBmagR[3]*hamilR[13]+19.36491673103709*hamilR[3]*BstarYdBmagR[12]+(8.660254037844386*BstarYdBmagR[0]-15.0*BstarYdBmagR[2])*hamilR[3])*q_*rdvpar2R))/(m_*q_); 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[12]+fL[12])+1.732050807568877*(fL[2]-1.0*fR[2])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[20]+fL[20])+3.0*(fL[5]-1.0*fR[5]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[22]+fL[22])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[3]+fL[3])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[26]+fL[26])+3.0*(fL[9]-1.0*fR[9]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[33]+fL[33])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[6]+fL[6]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[36]+fL[36])+1.732050807568877*(fL[16]-1.0*fR[16])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[38]+fL[38])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[10]+fL[10]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[19]-1.0*(8.660254037844387*fL[19]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[24]-1.0*(8.660254037844387*fL[24]+5.0*(fR[13]+fL[13]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[29]-1.0*(8.660254037844387*fL[29]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[45]+fL[45])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[17]+fL[17])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[21]+fL[21]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[34]-1.0*(8.660254037844387*fL[34]+5.0*(fR[23]+fL[23]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[35]-1.0*(8.660254037844387*fL[35]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[40]-1.0*(8.660254037844387*fL[40]+5.0*(fR[27]+fL[27]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[41]-1.0*(8.660254037844387*fL[41]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[30]+fL[30]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[37]+fL[37]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[46]-1.0*(8.660254037844387*fL[46]+5.0*(fR[39]+fL[39]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[42]+fL[42]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[12]-6.324555320336761*fL[12]-4.898979485566357*(fR[2]+fL[2])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[13]*fAvg[13]+alphaR[12]*fAvg[12]+alphaR[11]*fAvg[11]+alphaR[8]*fAvg[8]+alphaR[7]*fAvg[7]+alphaR[5]*fAvg[5]+alphaR[4]*fAvg[4]+alphaR[3]*fAvg[3]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[20]-94.8683298050514*fL[20]-73.48469228349535*(fR[5]+fL[5])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*(alphaR[5]*fAvg[13]+fAvg[5]*alphaR[13])-21.21320343559643*(alphaR[8]*fAvg[12]+fAvg[8]*alphaR[12])-18.97366596101028*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11]+alphaR[1]*fAvg[7]+fAvg[1]*alphaR[7])-21.21320343559643*(alphaR[3]*fAvg[5]+fAvg[3]*alphaR[5]+alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[22]-94.8683298050514*fL[22]-73.48469228349535*(fR[7]+fL[7])+42.42640687119286*fR[3]-42.42640687119286*fL[3])*amax-21.21320343559643*alphaR[13]*fAvg[17]-18.97366596101028*(alphaR[4]*fAvg[12]+fAvg[4]*alphaR[12])-21.21320343559643*(alphaR[7]*fAvg[11]+fAvg[7]*alphaR[11])-21.21320343559643*alphaR[5]*fAvg[10]-18.97366596101028*(alphaR[2]*fAvg[8]+fAvg[2]*alphaR[8])-21.21320343559643*(alphaR[3]*fAvg[6]+alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.008333333333333333*((94.8683298050514*fR[26]-94.8683298050514*fL[26]-73.48469228349535*(fR[9]+fL[9])+42.42640687119286*fR[4]-42.42640687119286*fL[4])*amax-21.21320343559643*(alphaR[12]*fAvg[18]+alphaR[11]*fAvg[17])-18.97366596101028*alphaR[5]*fAvg[15]-21.21320343559643*(alphaR[8]*fAvg[14]+alphaR[7]*fAvg[13]+fAvg[7]*alphaR[13])-21.21320343559643*alphaR[4]*fAvg[10]-18.97366596101028*alphaR[3]*fAvg[9]-21.21320343559643*(alphaR[2]*fAvg[6]+alphaR[1]*fAvg[5]+fAvg[1]*alphaR[5]+alphaR[0]*fAvg[3]+fAvg[0]*alphaR[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[33]-94.86832980505142*fL[33]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[6]-42.42640687119286*fL[6])*amax-18.97366596101028*(alphaR[5]*fAvg[17]+fAvg[10]*alphaR[13])+((-16.97056274847715*alphaR[11])-18.97366596101028*alphaR[2])*fAvg[12]+((-16.97056274847715*fAvg[11])-18.97366596101028*fAvg[2])*alphaR[12]-18.97366596101028*(alphaR[1]*fAvg[11]+fAvg[1]*alphaR[11])-21.21320343559643*alphaR[3]*fAvg[10]-18.97366596101028*(alphaR[4]*fAvg[8]+fAvg[4]*alphaR[8]+alphaR[4]*fAvg[7]+fAvg[4]*alphaR[7])-21.21320343559643*(alphaR[5]*fAvg[6]+alphaR[0]*fAvg[4]+fAvg[0]*alphaR[4]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[36]-94.86832980505142*fL[36]-73.48469228349535*(fR[16]+fL[16])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax-21.21320343559643*alphaR[8]*fAvg[18]-18.97366596101028*alphaR[4]*fAvg[17]+((-16.97056274847715*alphaR[13])-18.97366596101028*alphaR[3])*fAvg[15]-21.21320343559643*alphaR[12]*fAvg[14]-18.97366596101028*(alphaR[1]*fAvg[13]+fAvg[1]*alphaR[13])+fAvg[10]*((-18.97366596101028*alphaR[11])-21.21320343559643*alphaR[2])-18.97366596101028*(alphaR[5]*(fAvg[9]+fAvg[7])+fAvg[5]*alphaR[7])-21.21320343559643*(alphaR[4]*fAvg[6]+alphaR[0]*fAvg[5]+fAvg[0]*alphaR[5]+alphaR[1]*fAvg[3]+fAvg[1]*alphaR[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[38]-94.86832980505142*fL[38]-73.48469228349535*(fR[18]+fL[18])+42.42640687119286*fR[10]-42.42640687119286*fL[10])*amax-18.97366596101028*(alphaR[5]*fAvg[19]+alphaR[4]*fAvg[18])-21.21320343559643*alphaR[7]*fAvg[17]-18.97366596101028*(alphaR[3]*fAvg[16]+alphaR[2]*fAvg[14])-21.21320343559643*(alphaR[11]*fAvg[13]+fAvg[11]*alphaR[13])+fAvg[10]*((-18.97366596101028*alphaR[12])-21.21320343559643*alphaR[1])-18.97366596101028*fAvg[6]*alphaR[8]-21.21320343559643*(alphaR[0]*fAvg[6]+alphaR[4]*fAvg[5]+fAvg[4]*alphaR[5]+alphaR[2]*fAvg[3]+fAvg[2]*alphaR[3])); 
  Ghat[7] = 0.001190476190476191*((514.3928459844675*(fR[19]+fL[19])-296.98484809835*fR[11]+296.98484809835*fL[11])*amax+(94.86832980505142*alphaR[13]+148.492424049175*alphaR[3])*fAvg[13]+148.492424049175*fAvg[3]*alphaR[13]+132.815661727072*alphaR[12]*fAvg[12]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[11]+148.492424049175*fAvg[2]*alphaR[11]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[7]+148.492424049175*fAvg[0]*alphaR[7]+132.815661727072*(alphaR[5]*fAvg[5]+alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1])); 
  Ghat[8] = 0.001190476190476191*((514.3928459844675*(fR[24]+fL[24])-296.98484809835*fR[13]+296.98484809835*fL[13])*amax+148.492424049175*(alphaR[5]*fAvg[18]+alphaR[3]*fAvg[14])+(94.86832980505142*alphaR[12]+148.492424049175*alphaR[1])*fAvg[12]+148.492424049175*fAvg[1]*alphaR[12]+132.815661727072*alphaR[11]*fAvg[11]+(94.86832980505142*alphaR[8]+148.492424049175*alphaR[0])*fAvg[8]+148.492424049175*fAvg[0]*alphaR[8]+132.815661727072*(alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[29]+fL[29])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[4]*fAvg[19]+21.21320343559643*(alphaR[2]*fAvg[16]+alphaR[1]*fAvg[15])+18.97366596101028*alphaR[13]*fAvg[13]+21.21320343559643*alphaR[0]*fAvg[9]+18.97366596101028*(alphaR[5]*fAvg[5]+alphaR[3]*fAvg[3])); 
  Ghat[10] = -0.001666666666666667*((474.341649025257*fR[45]-474.341649025257*fL[45]-367.4234614174767*(fR[31]+fL[31])+212.1320343559643*fR[17]-212.1320343559643*fL[17])*amax+((-84.85281374238573*alphaR[13])-94.86832980505142*alphaR[3])*fAvg[19]+((-84.85281374238573*alphaR[11])-94.86832980505142*alphaR[2])*fAvg[18]+((-84.85281374238573*alphaR[12])-94.86832980505142*alphaR[1])*fAvg[17]-94.8683298050514*(alphaR[5]*fAvg[16]+alphaR[4]*(fAvg[14]+fAvg[13])+fAvg[4]*alphaR[13]+fAvg[6]*alphaR[12]+alphaR[5]*fAvg[11]+fAvg[5]*alphaR[11])+((-94.86832980505142*(alphaR[8]+alphaR[7]))-106.0660171779821*alphaR[0])*fAvg[10]-106.0660171779821*(alphaR[1]*fAvg[6]+alphaR[2]*fAvg[5]+fAvg[2]*alphaR[5]+alphaR[3]*fAvg[4]+fAvg[3]*alphaR[4])); 
  Ghat[11] = 0.001190476190476191*((514.3928459844675*(fR[32]+fL[32])-296.98484809835*fR[21]+296.98484809835*fL[21])*amax+(94.86832980505142*alphaR[13]+148.492424049175*alphaR[3])*fAvg[17]+148.492424049175*fAvg[6]*alphaR[13]+118.79393923934*(alphaR[4]*fAvg[12]+fAvg[4]*alphaR[12])+(132.815661727072*alphaR[8]+94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[11]+(132.815661727072*fAvg[8]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaR[11]+132.815661727072*alphaR[5]*fAvg[10]+148.492424049175*(alphaR[2]*fAvg[7]+fAvg[2]*alphaR[7])+132.815661727072*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])); 
  Ghat[12] = 0.001190476190476191*((514.3928459844675*(fR[34]+fL[34])-296.98484809835*fR[23]+296.98484809835*fL[23])*amax+(132.815661727072*alphaR[13]+148.492424049175*alphaR[3])*fAvg[18]+148.492424049175*alphaR[5]*fAvg[14]+(94.86832980505142*alphaR[8]+132.815661727072*alphaR[7]+148.492424049175*alphaR[0])*fAvg[12]+(94.86832980505142*fAvg[8]+132.815661727072*fAvg[7]+148.492424049175*fAvg[0])*alphaR[12]+118.79393923934*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11])+148.492424049175*(alphaR[1]*fAvg[8]+fAvg[1]*alphaR[8])+132.815661727072*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])); 
  Ghat[13] = 0.001190476190476191*((514.3928459844675*(fR[35]+fL[35])-296.98484809835*fR[25]+296.98484809835*fL[25])*amax+132.815661727072*alphaR[12]*fAvg[18]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[17]+118.79393923934*alphaR[5]*fAvg[15]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[13]+(132.815661727072*fAvg[9]+94.86832980505142*fAvg[7]+148.492424049175*fAvg[0])*alphaR[13]+148.492424049175*fAvg[6]*alphaR[11]+132.815661727072*alphaR[4]*fAvg[10]+148.492424049175*(alphaR[3]*fAvg[7]+fAvg[3]*alphaR[7])+132.815661727072*(alphaR[1]*fAvg[5]+fAvg[1]*alphaR[5])); 
  Ghat[14] = 0.001190476190476191*((514.3928459844675*(fR[40]+fL[40])-296.98484809835*fR[27]+296.98484809835*fL[27])*amax+(94.86832980505142*alphaR[12]+148.492424049175*alphaR[1])*fAvg[18]+132.815661727072*alphaR[11]*fAvg[17]+(94.86832980505142*alphaR[8]+148.492424049175*alphaR[0])*fAvg[14]+148.492424049175*(alphaR[5]*fAvg[12]+fAvg[5]*alphaR[12])+132.815661727072*alphaR[4]*fAvg[10]+148.492424049175*(alphaR[3]*fAvg[8]+fAvg[3]*alphaR[8])+132.815661727072*alphaR[2]*fAvg[6]); 
  Ghat[15] = 0.008333333333333333*((73.48469228349536*(fR[41]+fL[41])-42.42640687119286*fR[28]+42.42640687119286*fL[28])*amax+(18.97366596101028*alphaR[11]+21.21320343559643*alphaR[2])*fAvg[19]+21.21320343559643*alphaR[4]*fAvg[16]+(18.97366596101028*alphaR[7]+21.21320343559643*alphaR[0])*fAvg[15]+16.97056274847715*(alphaR[5]*fAvg[13]+fAvg[5]*alphaR[13])+21.21320343559643*alphaR[1]*fAvg[9]+18.97366596101028*(alphaR[3]*fAvg[5]+fAvg[3]*alphaR[5])); 
  Ghat[16] = 0.008333333333333333*((73.48469228349536*(fR[43]+fL[43])-42.42640687119286*fR[30]+42.42640687119286*fL[30])*amax+(18.97366596101028*alphaR[12]+21.21320343559643*alphaR[1])*fAvg[19]+18.97366596101028*alphaR[13]*fAvg[17]+(18.97366596101028*alphaR[8]+21.21320343559643*alphaR[0])*fAvg[16]+21.21320343559643*alphaR[4]*fAvg[15]+18.97366596101028*alphaR[5]*fAvg[10]+21.21320343559643*alphaR[2]*fAvg[9]+18.97366596101028*alphaR[3]*fAvg[6]); 
  Ghat[17] = 2.380952380952381e-4*((2571.964229922338*(fR[44]+fL[44])-1484.92424049175*fR[37]+1484.92424049175*fL[37])*amax+593.9696961967002*(alphaR[5]*fAvg[19]+alphaR[4]*fAvg[18])+(664.0783086353599*alphaR[8]+474.3416490252571*alphaR[7]+742.462120245875*alphaR[0])*fAvg[17]+664.0783086353599*(alphaR[13]*fAvg[16]+alphaR[11]*fAvg[14])+(474.3416490252571*alphaR[11]+742.462120245875*alphaR[2])*fAvg[13]+(474.3416490252571*fAvg[11]+742.462120245875*fAvg[2])*alphaR[13]+593.9696961967002*fAvg[10]*alphaR[12]+742.462120245875*(alphaR[3]*fAvg[11]+fAvg[3]*alphaR[11])+664.0783086353599*alphaR[1]*fAvg[10]+742.462120245875*fAvg[6]*alphaR[7]+664.0783086353599*(alphaR[4]*fAvg[5]+fAvg[4]*alphaR[5])); 
  Ghat[18] = 2.380952380952381e-4*((2571.964229922338*(fR[46]+fL[46])-1484.92424049175*fR[39]+1484.92424049175*fL[39])*amax+(474.3416490252571*alphaR[8]+664.0783086353599*alphaR[7]+742.462120245875*alphaR[0])*fAvg[18]+593.9696961967002*alphaR[4]*fAvg[17]+(474.3416490252571*alphaR[12]+742.462120245875*alphaR[1])*fAvg[14]+664.0783086353599*(alphaR[12]*fAvg[13]+fAvg[12]*alphaR[13])+742.462120245875*(alphaR[3]*fAvg[12]+fAvg[3]*alphaR[12])+fAvg[10]*(593.9696961967002*alphaR[11]+664.0783086353599*alphaR[2])+742.462120245875*(alphaR[5]*fAvg[8]+fAvg[5]*alphaR[8])+664.0783086353599*alphaR[4]*fAvg[6]); 
  Ghat[19] = 0.001666666666666667*((367.4234614174769*(fR[47]+fL[47])-212.1320343559643*fR[42]+212.1320343559643*fL[42])*amax+(94.86832980505142*(alphaR[8]+alphaR[7])+106.0660171779821*alphaR[0])*fAvg[19]+84.85281374238573*alphaR[5]*fAvg[17]+(94.86832980505142*alphaR[12]+106.0660171779822*alphaR[1])*fAvg[16]+(94.86832980505142*alphaR[11]+106.0660171779822*alphaR[2])*fAvg[15]+fAvg[10]*(84.85281374238573*alphaR[13]+94.86832980505142*alphaR[3])+106.0660171779821*alphaR[4]*fAvg[9]+94.86832980505142*alphaR[5]*fAvg[6]); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = -1.224744871391589*Ghat[0]; 
  incr[3] = 0.7071067811865475*Ghat[2]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = -1.224744871391589*Ghat[1]; 
  incr[6] = 0.7071067811865475*Ghat[4]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = -1.224744871391589*Ghat[3]; 
  incr[10] = 0.7071067811865475*Ghat[6]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 1.58113883008419*Ghat[0]; 
  incr[13] = 0.7071067811865475*Ghat[8]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = -1.224744871391589*Ghat[5]; 
  incr[17] = 0.7071067811865475*Ghat[10]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = -1.224744871391589*Ghat[7]; 
  incr[20] = 1.58113883008419*Ghat[1]; 
  incr[21] = 0.7071067811865475*Ghat[11]; 
  incr[22] = 1.58113883008419*Ghat[2]; 
  incr[23] = 0.7071067811865475*Ghat[12]; 
  incr[24] = -1.224744871391589*Ghat[8]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 1.58113883008419*Ghat[3]; 
  incr[27] = 0.7071067811865475*Ghat[14]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = -1.224744871391589*Ghat[9]; 
  incr[30] = 0.7071067811865475*Ghat[16]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = 1.58113883008419*Ghat[4]; 
  incr[34] = -1.224744871391589*Ghat[12]; 
  incr[35] = -1.224744871391589*Ghat[13]; 
  incr[36] = 1.58113883008419*Ghat[5]; 
  incr[37] = 0.7071067811865475*Ghat[17]; 
  incr[38] = 1.58113883008419*Ghat[6]; 
  incr[39] = 0.7071067811865475*Ghat[18]; 
  incr[40] = -1.224744871391589*Ghat[14]; 
  incr[41] = -1.224744871391589*Ghat[15]; 
  incr[42] = 0.7071067811865475*Ghat[19]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = 1.58113883008419*Ghat[10]; 
  incr[46] = -1.224744871391589*Ghat[18]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

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
  outR[16] += incr[16]*rdy2R; 
  outR[17] += incr[17]*rdy2R; 
  outR[18] += incr[18]*rdy2R; 
  outR[19] += incr[19]*rdy2R; 
  outR[20] += incr[20]*rdy2R; 
  outR[21] += incr[21]*rdy2R; 
  outR[22] += incr[22]*rdy2R; 
  outR[23] += incr[23]*rdy2R; 
  outR[24] += incr[24]*rdy2R; 
  outR[25] += incr[25]*rdy2R; 
  outR[26] += incr[26]*rdy2R; 
  outR[27] += incr[27]*rdy2R; 
  outR[28] += incr[28]*rdy2R; 
  outR[29] += incr[29]*rdy2R; 
  outR[30] += incr[30]*rdy2R; 
  outR[31] += incr[31]*rdy2R; 
  outR[32] += incr[32]*rdy2R; 
  outR[33] += incr[33]*rdy2R; 
  outR[34] += incr[34]*rdy2R; 
  outR[35] += incr[35]*rdy2R; 
  outR[36] += incr[36]*rdy2R; 
  outR[37] += incr[37]*rdy2R; 
  outR[38] += incr[38]*rdy2R; 
  outR[39] += incr[39]*rdy2R; 
  outR[40] += incr[40]*rdy2R; 
  outR[41] += incr[41]*rdy2R; 
  outR[42] += incr[42]*rdy2R; 
  outR[43] += incr[43]*rdy2R; 
  outR[44] += incr[44]*rdy2R; 
  outR[45] += incr[45]*rdy2R; 
  outR[46] += incr[46]*rdy2R; 
  outR[47] += incr[47]*rdy2R; 

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
  outL[11] += -1.0*incr[11]*rdy2L; 
  outL[12] += -1.0*incr[12]*rdy2L; 
  outL[13] += -1.0*incr[13]*rdy2L; 
  outL[14] += -1.0*incr[14]*rdy2L; 
  outL[15] += incr[15]*rdy2L; 
  outL[16] += incr[16]*rdy2L; 
  outL[17] += -1.0*incr[17]*rdy2L; 
  outL[18] += incr[18]*rdy2L; 
  outL[19] += incr[19]*rdy2L; 
  outL[20] += -1.0*incr[20]*rdy2L; 
  outL[21] += -1.0*incr[21]*rdy2L; 
  outL[22] += -1.0*incr[22]*rdy2L; 
  outL[23] += -1.0*incr[23]*rdy2L; 
  outL[24] += incr[24]*rdy2L; 
  outL[25] += -1.0*incr[25]*rdy2L; 
  outL[26] += -1.0*incr[26]*rdy2L; 
  outL[27] += -1.0*incr[27]*rdy2L; 
  outL[28] += -1.0*incr[28]*rdy2L; 
  outL[29] += incr[29]*rdy2L; 
  outL[30] += -1.0*incr[30]*rdy2L; 
  outL[31] += incr[31]*rdy2L; 
  outL[32] += incr[32]*rdy2L; 
  outL[33] += -1.0*incr[33]*rdy2L; 
  outL[34] += incr[34]*rdy2L; 
  outL[35] += incr[35]*rdy2L; 
  outL[36] += -1.0*incr[36]*rdy2L; 
  outL[37] += -1.0*incr[37]*rdy2L; 
  outL[38] += -1.0*incr[38]*rdy2L; 
  outL[39] += -1.0*incr[39]*rdy2L; 
  outL[40] += incr[40]*rdy2L; 
  outL[41] += incr[41]*rdy2L; 
  outL[42] += -1.0*incr[42]*rdy2L; 
  outL[43] += incr[43]*rdy2L; 
  outL[44] += incr[44]*rdy2L; 
  outL[45] += -1.0*incr[45]*rdy2L; 
  outL[46] += incr[46]*rdy2L; 
  outL[47] += incr[47]*rdy2L; 

  return std::abs(alphaSurfAvgR);
} 
double EmGyrokineticGenGeoSurf2x2vSer_vpar_P2_Bvarsx(const double q_, const double m_, const double cflL, const double cflR, const double *wL, const double *dxvL, const double *wR, const double *dxvR, const double amax_in, const double *bmag, const double *jacobTotInv, const double *cmag, const double *b_x, const double *b_y, const double *b_z, const double *phi, const double *Apar, const double *dApardt, const double *dApardtPrev, const double *fL, const double *fR, double *outL, double *outR, double *emModL, double *emModR) 
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

  double hamilR[48]; 
  hamilR[0] = (0.6666666666666666*(3.0*rdvpar2SqR*(m_*wvparSqR+bmag[0]*wmuR+phi[0]*q_)+m_))/rdvpar2SqR; 
  hamilR[1] = 2.0*(bmag[1]*wmuR+phi[1]*q_); 
  hamilR[2] = 2.0*phi[2]*q_; 
  hamilR[3] = (2.309401076758503*m_*wvparR)/rdvpar2R; 
  hamilR[4] = (1.154700538379252*bmag[0])/rdmu2R; 
  hamilR[5] = 2.0*phi[3]*q_; 
  hamilR[8] = (1.154700538379252*bmag[1])/rdmu2R; 
  hamilR[11] = 2.0*(bmag[4]*wmuR+phi[4]*q_); 
  hamilR[12] = 2.0*phi[5]*q_; 
  hamilR[13] = (0.5962847939999438*m_)/rdvpar2SqR; 
  hamilR[19] = 2.0*phi[6]*q_; 
  hamilR[20] = 2.0*phi[7]*q_; 
  hamilR[25] = (1.154700538379251*bmag[4])/rdmu2R; 

  double BstarXdBmagR[48]; 
  BstarXdBmagR[0] = 0.01428571428571429*(3.872983346207417*((10.0*b_z[4]+15.65247584249853*b_z[0])*jacobTotInv[4]+7.0*(2.23606797749979*jacobTotInv[0]*b_z[4]+2.0*b_z[1]*jacobTotInv[1]))*Apar[6]+12.12435565298214*((5.0*Apar[2]*b_z[4]+4.47213595499958*b_z[1]*Apar[3])*jacobTotInv[4]+4.47213595499958*jacobTotInv[1]*Apar[3]*b_z[4]+5.0*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+(b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[2])))*rdy2R; 
  BstarXdBmagR[1] = 0.01428571428571429*(3.872983346207417*(24.59674775249769*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+14.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[6]+1.732050807568877*((55.0*Apar[3]*b_z[4]+31.30495168499706*(b_z[0]*Apar[3]+b_z[1]*Apar[2]))*jacobTotInv[4]+7.0*(4.47213595499958*(jacobTotInv[0]*Apar[3]+jacobTotInv[1]*Apar[2])*b_z[4]+(9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[3]+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[2])))*rdy2R; 
  BstarXdBmagR[2] = 0.1*(3.872983346207417*(4.47213595499958*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[7]+19.36491673103709*(b_z[4]*jacobTotInv[4]+b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[5])*rdy2R; 
  BstarXdBmagR[5] = 0.01428571428571429*(3.872983346207417*((55.0*b_z[4]+31.30495168499706*b_z[0])*jacobTotInv[4]+7.0*(4.47213595499958*jacobTotInv[0]*b_z[4]+9.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0]))*Apar[7]+12.12435565298214*(10.0*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+11.18033988749895*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[5])*rdy2R; 
  BstarXdBmagR[11] = 0.01428571428571429*(3.872983346207417*((33.54101966249685*b_z[4]+10.0*b_z[0])*jacobTotInv[4]+10.0*jacobTotInv[0]*b_z[4]+2.23606797749979*(11.0*b_z[1]*jacobTotInv[1]+7.0*b_z[0]*jacobTotInv[0]))*Apar[6]+1.732050807568877*((22.3606797749979*Apar[2]*b_z[4]+5.0*(11.0*b_z[1]*Apar[3]+7.0*b_z[0]*Apar[2]))*jacobTotInv[4]+5.0*(11.0*jacobTotInv[1]*Apar[3]+7.0*jacobTotInv[0]*Apar[2])*b_z[4]+31.30495168499706*((b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])))*rdy2R; 
  BstarXdBmagR[19] = 0.01428571428571429*(1.732050807568877*(122.9837387624885*(b_z[1]*jacobTotInv[4]+jacobTotInv[1]*b_z[4])+70.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[7]+3.872983346207417*((22.3606797749979*b_z[4]+35.0*b_z[0])*jacobTotInv[4]+7.0*(5.0*jacobTotInv[0]*b_z[4]+4.47213595499958*b_z[1]*jacobTotInv[1]))*Apar[5])*rdy2R; 

  double BstarYdBmagR[48]; 
  BstarYdBmagR[0] = -(0.8660254037844386*rdx2R*(2.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*wvparR+(3.0*(Apar[1]*b_z[4]+b_z[1]*Apar[4])*jacobTotInv[4]+(4.0*jacobTotInv[1]*Apar[4]+2.23606797749979*(Apar[0]*jacobTotInv[1]+jacobTotInv[0]*Apar[1]))*b_z[4]+2.23606797749979*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1])*Apar[4]+2.0*Apar[1]*b_z[1]*jacobTotInv[1]+jacobTotInv[0]*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*q_))/q_; 
  BstarYdBmagR[1] = -(0.02474358296526967*rdx2R*(14.0*(b_z[4]*(10.0*jacobTotInv[4]+11.18033988749895*jacobTotInv[0])+5.0*b_z[1]*jacobTotInv[1])*m_*wvparR+(2.0*((122.9837387624885*Apar[4]+35.0*Apar[0])*b_z[4]+7.0*(5.0*b_z[0]*Apar[4]+4.47213595499958*Apar[1]*b_z[1]))*jacobTotInv[4]+7.0*((20.0*jacobTotInv[0]*Apar[4]+2.23606797749979*(11.0*Apar[1]*jacobTotInv[1]+5.0*Apar[0]*jacobTotInv[0]))*b_z[4]+2.23606797749979*(11.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0])*Apar[4]+5.0*((Apar[0]*b_z[1]+b_z[0]*Apar[1])*jacobTotInv[1]+2.0*jacobTotInv[0]*Apar[1]*b_z[1])))*q_))/q_; 
  BstarYdBmagR[2] = -0.05773502691896259*(6.708203932499369*(2.23606797749979*(3.0*b_z[1]*jacobTotInv[4]+4.0*jacobTotInv[1]*b_z[4])+5.0*(b_z[0]*jacobTotInv[1]+jacobTotInv[0]*b_z[1]))*Apar[6]+3.0*(b_z[4]*(15.0*Apar[3]*jacobTotInv[4]+11.18033988749895*(jacobTotInv[0]*Apar[3]+jacobTotInv[1]*Apar[2]))+5.0*((2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0])*Apar[3]+jacobTotInv[0]*b_z[1]*Apar[2])))*rdx2R; 
  BstarYdBmagR[3] = -(1.0*(2.23606797749979*jacobTotInv[1]*b_z[4]+jacobTotInv[0]*b_z[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[5] = -0.01428571428571429*(3.872983346207417*(2.0*(55.0*b_z[4]+15.65247584249853*b_z[0])*jacobTotInv[4]+7.0*(8.94427190999916*jacobTotInv[0]*b_z[4]+11.0*b_z[1]*jacobTotInv[1]+5.0*b_z[0]*jacobTotInv[0]))*Apar[6]+12.12435565298214*(2.0*(5.0*Apar[2]*b_z[4]+4.47213595499958*b_z[1]*Apar[3])*jacobTotInv[4]+2.23606797749979*(11.0*jacobTotInv[1]*Apar[3]+5.0*jacobTotInv[0]*Apar[2])*b_z[4]+5.0*((b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1])*Apar[3]+b_z[1]*jacobTotInv[1]*Apar[2])))*rdx2R; 
  BstarYdBmagR[6] = -(1.0*(b_z[4]*(2.0*jacobTotInv[4]+2.23606797749979*jacobTotInv[0])+b_z[1]*jacobTotInv[1])*m_*rdx2R)/(q_*rdvpar2R); 
  BstarYdBmagR[11] = -(0.02474358296526967*rdx2R*(70.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*wvparR+((145.3444185374863*(Apar[1]*b_z[4]+b_z[1]*Apar[4])+35.0*(Apar[0]*b_z[1]+b_z[0]*Apar[1]))*jacobTotInv[4]+(245.9674775249769*jacobTotInv[1]*Apar[4]+35.0*(2.0*Apar[0]*jacobTotInv[1]+3.0*jacobTotInv[0]*Apar[1]))*b_z[4]+7.0*(5.0*(2.0*b_z[0]*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1])*Apar[4]+8.94427190999916*Apar[1]*b_z[1]*jacobTotInv[1]))*q_))/q_; 
  BstarYdBmagR[12] = -0.1*(3.872983346207417*(b_z[4]*(6.708203932499369*jacobTotInv[4]+5.0*jacobTotInv[0])+2.23606797749979*(2.0*b_z[1]*jacobTotInv[1]+b_z[0]*jacobTotInv[0]))*Apar[7]+1.732050807568877*(11.18033988749895*jacobTotInv[1]*b_z[4]+5.0*jacobTotInv[0]*b_z[1])*Apar[5])*rdx2R; 
  BstarYdBmagR[19] = -0.01844277783908294*(6.708203932499369*(2.23606797749979*(13.0*b_z[1]*jacobTotInv[4]+22.0*jacobTotInv[1]*b_z[4])+7.0*(2.0*b_z[0]*jacobTotInv[1]+3.0*jacobTotInv[0]*b_z[1]))*Apar[6]+3.0*((65.0*Apar[3]*b_z[4]+15.65247584249853*(b_z[0]*Apar[3]+b_z[1]*Apar[2]))*jacobTotInv[4]+7.0*(2.23606797749979*(3.0*jacobTotInv[0]*Apar[3]+2.0*jacobTotInv[1]*Apar[2])*b_z[4]+4.0*b_z[1]*jacobTotInv[1]*Apar[3])))*rdx2R; 
  BstarYdBmagR[20] = -0.1*(1.732050807568877*(2.23606797749979*(4.0*b_z[1]*jacobTotInv[4]+11.0*jacobTotInv[1]*b_z[4])+5.0*(b_z[0]*jacobTotInv[1]+2.0*jacobTotInv[0]*b_z[1]))*Apar[7]+3.872983346207417*(b_z[4]*(4.47213595499958*jacobTotInv[4]+5.0*jacobTotInv[0])+2.23606797749979*b_z[1]*jacobTotInv[1])*Apar[5])*rdx2R; 
  BstarYdBmagR[21] = -(1.0*(b_z[1]*jacobTotInv[4]+2.0*jacobTotInv[1]*b_z[4])*m_*rdx2R)/(q_*rdvpar2R); 

  double alphaR[20]; 
  alphaR[0] = (0.03535533905932736*((15.0*hamilR[19]*BstarYdBmagR[21]-19.36491673103708*BstarYdBmagR[5]*hamilR[20]-8.660254037844387*BstarYdBmagR[11]*hamilR[19]-19.36491673103709*BstarYdBmagR[2]*hamilR[12]+hamilR[5]*(15.0*BstarYdBmagR[6]-8.660254037844386*BstarYdBmagR[1])+hamilR[2]*(15.0*BstarYdBmagR[3]-8.660254037844386*BstarYdBmagR[0]))*rdy2R+((-19.36491673103708*BstarXdBmagR[5]*hamilR[19])-19.36491673103709*BstarXdBmagR[1]*hamilR[11]-8.660254037844386*(BstarXdBmagR[2]*hamilR[5]+BstarXdBmagR[0]*hamilR[1]))*rdx2R-40.0*dApardt[0]*q_))/m_; 
  alphaR[1] = (0.03535533905932736*((13.41640786499874*hamilR[5]*BstarYdBmagR[21]+((-17.32050807568877*BstarYdBmagR[19])-19.36491673103708*BstarYdBmagR[2])*hamilR[20]+(13.41640786499874*BstarYdBmagR[6]-7.745966692414834*BstarYdBmagR[1])*hamilR[19]-19.36491673103709*BstarYdBmagR[5]*hamilR[12]-7.745966692414834*hamilR[5]*BstarYdBmagR[11]+15.0*hamilR[2]*BstarYdBmagR[6]+(15.0*BstarYdBmagR[3]-8.660254037844386*BstarYdBmagR[0])*hamilR[5]-8.660254037844386*BstarYdBmagR[1]*hamilR[2])*rdy2R+(((-17.32050807568877*BstarXdBmagR[19])-19.36491673103708*BstarXdBmagR[2])*hamilR[19]+((-17.32050807568877*BstarXdBmagR[11])-19.36491673103709*BstarXdBmagR[0])*hamilR[11]-8.660254037844386*(BstarXdBmagR[5]*hamilR[5]+BstarXdBmagR[1]*hamilR[1]))*rdx2R-40.0*dApardt[1]*q_))/m_; 
  alphaR[2] = -(0.03535533905932736*(((17.32050807568877*BstarYdBmagR[20]-33.54101966249684*BstarYdBmagR[6]+19.36491673103708*BstarYdBmagR[1])*hamilR[20]+8.660254037844386*BstarYdBmagR[19]*hamilR[19]+(17.32050807568877*BstarYdBmagR[12]-33.54101966249685*BstarYdBmagR[3]+19.36491673103709*BstarYdBmagR[0])*hamilR[12]+8.660254037844386*(BstarYdBmagR[5]*hamilR[5]+BstarYdBmagR[2]*hamilR[2]))*rdy2R+(7.745966692414834*BstarXdBmagR[2]*hamilR[20]+19.36491673103708*BstarXdBmagR[1]*hamilR[19]+19.36491673103709*BstarXdBmagR[5]*hamilR[11]+8.660254037844386*(BstarXdBmagR[0]*hamilR[5]+hamilR[1]*BstarXdBmagR[2]))*rdx2R+40.0*dApardt[2]*q_))/m_; 
  alphaR[3] = -(0.1767766952966368*(3.872983346207417*BstarXdBmagR[1]*hamilR[25]+1.732050807568877*BstarXdBmagR[0]*hamilR[8])*rdx2R)/m_; 
  alphaR[4] = (0.03535533905932736*((hamilR[20]*(30.0*BstarYdBmagR[21]-17.32050807568877*(BstarYdBmagR[12]+BstarYdBmagR[11])+33.54101966249684*BstarYdBmagR[3]-19.36491673103708*BstarYdBmagR[0])-17.32050807568877*hamilR[12]*BstarYdBmagR[20]-7.745966692414834*(BstarYdBmagR[5]*hamilR[19]+hamilR[5]*BstarYdBmagR[19])+(33.54101966249685*BstarYdBmagR[6]-19.36491673103709*BstarYdBmagR[1])*hamilR[12]-8.660254037844386*(BstarYdBmagR[2]*hamilR[5]+hamilR[2]*BstarYdBmagR[5]))*rdy2R+((-7.745966692414834*BstarXdBmagR[5]*hamilR[20])+((-17.32050807568877*BstarXdBmagR[11])-19.36491673103708*BstarXdBmagR[0])*hamilR[19]+hamilR[11]*((-17.32050807568877*BstarXdBmagR[19])-19.36491673103709*BstarXdBmagR[2])-8.660254037844386*(BstarXdBmagR[1]*hamilR[5]+hamilR[1]*BstarXdBmagR[5]))*rdx2R-40.0*dApardt[3]*q_))/m_; 
  alphaR[5] = -(0.03535533905932736*((17.32050807568877*BstarXdBmagR[11]+19.36491673103708*BstarXdBmagR[0])*hamilR[25]+8.660254037844386*BstarXdBmagR[1]*hamilR[8])*rdx2R)/m_; 
  alphaR[6] = -(0.1767766952966368*(3.872983346207417*BstarXdBmagR[5]*hamilR[25]+1.732050807568877*BstarXdBmagR[2]*hamilR[8])*rdx2R)/m_; 
  alphaR[7] = (0.005050762722761052*(((67.0820393249937*hamilR[19]+105.0*hamilR[2])*BstarYdBmagR[21]-121.2435565298214*BstarYdBmagR[5]*hamilR[20]+((-38.72983346207417*BstarYdBmagR[11])+105.0*BstarYdBmagR[3]-60.62177826491071*BstarYdBmagR[0])*hamilR[19]-135.5544171172596*hamilR[12]*BstarYdBmagR[19]-60.6217782649107*hamilR[2]*BstarYdBmagR[11]+hamilR[5]*(93.91485505499116*BstarYdBmagR[6]-54.22176684690384*BstarYdBmagR[1]))*rdy2R+((-121.2435565298214*BstarXdBmagR[5]*hamilR[19])-60.62177826491071*hamilR[5]*BstarXdBmagR[19]-121.2435565298214*BstarXdBmagR[1]*hamilR[11]-60.6217782649107*hamilR[1]*BstarXdBmagR[11])*rdx2R-280.0*dApardt[4]*q_))/m_; 
  alphaR[8] = -(0.07905694150420944*((7.745966692414834*BstarYdBmagR[5]*hamilR[20]+3.872983346207417*hamilR[5]*BstarYdBmagR[20]+7.745966692414834*BstarYdBmagR[2]*hamilR[12]+3.872983346207417*hamilR[2]*BstarYdBmagR[12])*rdy2R+(3.872983346207417*BstarXdBmagR[0]*hamilR[20]+7.745966692414834*BstarXdBmagR[5]*hamilR[19]+3.464101615137754*BstarXdBmagR[2]*hamilR[5])*rdx2R+17.88854381999832*dApardt[5]*q_))/m_; 
  alphaR[10] = -(0.1767766952966368*((3.464101615137754*BstarXdBmagR[19]+3.872983346207417*BstarXdBmagR[2])*hamilR[25]+1.732050807568877*BstarXdBmagR[5]*hamilR[8])*rdx2R)/m_; 
  alphaR[11] = (0.005050762722761052*((234.787137637478*hamilR[12]*BstarYdBmagR[21]+((-108.4435336938077*BstarYdBmagR[20])+210.0*BstarYdBmagR[6]-121.2435565298214*BstarYdBmagR[1])*hamilR[20]+((-38.72983346207417*BstarYdBmagR[19])-60.6217782649107*BstarYdBmagR[2])*hamilR[19]-60.6217782649107*hamilR[2]*BstarYdBmagR[19]-135.5544171172596*BstarYdBmagR[11]*hamilR[12]-54.22176684690384*BstarYdBmagR[5]*hamilR[5])*rdy2R+((-54.22176684690384*BstarXdBmagR[19]*hamilR[20])-121.2435565298214*BstarXdBmagR[1]*hamilR[19]-60.6217782649107*hamilR[1]*BstarXdBmagR[19]-121.2435565298214*BstarXdBmagR[5]*hamilR[11]-60.62177826491071*hamilR[5]*BstarXdBmagR[11])*rdx2R-280.0*dApardt[6]*q_))/m_; 
  alphaR[12] = -(0.04564354645876382*(((12.0*BstarYdBmagR[19]+13.41640786499874*BstarYdBmagR[2])*hamilR[20]+(6.0*hamilR[19]+6.708203932499369*hamilR[2])*BstarYdBmagR[20]+13.41640786499874*BstarYdBmagR[5]*hamilR[12]+6.708203932499369*hamilR[5]*BstarYdBmagR[12])*rdy2R+(6.708203932499369*BstarXdBmagR[1]*hamilR[20]+(12.0*BstarXdBmagR[19]+13.41640786499874*BstarXdBmagR[2])*hamilR[19]+6.0*BstarXdBmagR[5]*hamilR[5])*rdx2R+30.98386676965934*dApardt[7]*q_))/m_; 
  alphaR[13] = -(0.03535533905932736*(17.32050807568877*BstarXdBmagR[1]*hamilR[25]+8.660254037844387*hamilR[8]*BstarXdBmagR[11])*rdx2R)/m_; 
  alphaR[17] = -(0.03535533905932736*(17.32050807568877*BstarXdBmagR[5]*hamilR[25]+8.660254037844387*hamilR[8]*BstarXdBmagR[19])*rdx2R)/m_; 

  // Surface-averaged phase velocity in this direction.
  double alphaSurfAvgR = (0.00625*((15.0*hamilR[19]*BstarYdBmagR[21]-19.36491673103708*BstarYdBmagR[5]*hamilR[20]-8.660254037844387*BstarYdBmagR[11]*hamilR[19]-19.36491673103709*BstarYdBmagR[2]*hamilR[12]+15.0*hamilR[5]*BstarYdBmagR[6]-8.660254037844386*BstarYdBmagR[1]*hamilR[5]+15.0*hamilR[2]*BstarYdBmagR[3]-8.660254037844386*BstarYdBmagR[0]*hamilR[2])*rdy2R+((-19.36491673103708*BstarXdBmagR[5]*hamilR[19])-19.36491673103709*BstarXdBmagR[1]*hamilR[11]-8.660254037844386*BstarXdBmagR[2]*hamilR[5]-8.660254037844386*BstarXdBmagR[0]*hamilR[1])*rdx2R-40.0*dApardt[0]*q_))/m_; 

  double incr[48]; 
  double amax = amax_in; 

  double fAvg[20]; 
  fAvg[0] = 0.7071067811865475*(2.23606797749979*(fR[13]+fL[13])+1.732050807568877*(fL[3]-1.0*fR[3])+fR[0]+fL[0]); 
  fAvg[1] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[23]+fL[23])+3.0*(fL[6]-1.0*fR[6]))+3.0*(fR[1]+fL[1])); 
  fAvg[2] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[24]+fL[24])+3.0*(fL[7]-1.0*fR[7]))+3.0*(fR[2]+fL[2])); 
  fAvg[3] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[27]+fL[27])+3.0*(fL[10]-1.0*fR[10]))+3.0*(fR[4]+fL[4])); 
  fAvg[4] = 0.7071067811865475*(2.23606797749979*(fR[34]+fL[34])+1.732050807568877*(fL[15]-1.0*fR[15])+fR[5]+fL[5]); 
  fAvg[5] = 0.7071067811865475*(2.23606797749979*(fR[39]+fL[39])+1.732050807568877*(fL[17]-1.0*fR[17])+fR[8]+fL[8]); 
  fAvg[6] = 0.7071067811865475*(2.23606797749979*(fR[40]+fL[40])+1.732050807568877*(fL[18]-1.0*fR[18])+fR[9]+fL[9]); 
  fAvg[7] = -0.1414213562373095*(8.660254037844387*fR[21]-1.0*(8.660254037844387*fL[21]+5.0*(fR[11]+fL[11]))); 
  fAvg[8] = -0.1414213562373095*(8.660254037844387*fR[22]-1.0*(8.660254037844387*fL[22]+5.0*(fR[12]+fL[12]))); 
  fAvg[9] = -0.1414213562373095*(8.660254037844387*fR[30]-1.0*(8.660254037844387*fL[30]+5.0*(fR[14]+fL[14]))); 
  fAvg[10] = 0.2357022603955158*(1.732050807568877*(3.872983346207417*(fR[46]+fL[46])+3.0*(fL[31]-1.0*fR[31]))+3.0*(fR[16]+fL[16])); 
  fAvg[11] = -0.1414213562373095*(8.660254037844387*fR[32]-1.0*(8.660254037844387*fL[32]+5.0*(fR[19]+fL[19]))); 
  fAvg[12] = -0.1414213562373095*(8.660254037844387*fR[33]-1.0*(8.660254037844387*fL[33]+5.0*(fR[20]+fL[20]))); 
  fAvg[13] = -0.1414213562373095*(8.660254037844387*fR[37]-1.0*(8.660254037844387*fL[37]+5.0*(fR[25]+fL[25]))); 
  fAvg[14] = -0.1414213562373095*(8.660254037844387*fR[38]-1.0*(8.660254037844387*fL[38]+5.0*(fR[26]+fL[26]))); 
  fAvg[15] = -0.1414213562373095*(8.660254037844387*fR[42]-1.0*(8.660254037844387*fL[42]+5.0*(fR[28]+fL[28]))); 
  fAvg[16] = -0.1414213562373095*(8.660254037844387*fR[43]-1.0*(8.660254037844387*fL[43]+5.0*(fR[29]+fL[29]))); 
  fAvg[17] = -0.1414213562373095*(8.660254037844387*fR[44]-1.0*(8.660254037844387*fL[44]+5.0*(fR[35]+fL[35]))); 
  fAvg[18] = -0.1414213562373095*(8.660254037844387*fR[45]-1.0*(8.660254037844387*fL[45]+5.0*(fR[36]+fL[36]))); 
  fAvg[19] = -0.1414213562373095*(8.660254037844387*fR[47]-1.0*(8.660254037844387*fL[47]+5.0*(fR[41]+fL[41]))); 

  double Ghat[20]; 
  Ghat[0] = -0.125*((6.324555320336761*fR[13]-6.324555320336761*fL[13]-4.898979485566357*(fR[3]+fL[3])+2.828427124746191*fR[0]-2.828427124746191*fL[0])*amax-1.414213562373095*(alphaR[17]*fAvg[17]+alphaR[13]*fAvg[13]+alphaR[12]*fAvg[12]+alphaR[11]*fAvg[11]+alphaR[10]*fAvg[10]+alphaR[8]*fAvg[8]+alphaR[7]*fAvg[7]+alphaR[6]*fAvg[6]+alphaR[5]*fAvg[5]+alphaR[4]*fAvg[4]+alphaR[3]*fAvg[3]+alphaR[2]*fAvg[2]+alphaR[1]*fAvg[1]+alphaR[0]*fAvg[0])); 
  Ghat[1] = -0.008333333333333333*((94.8683298050514*fR[23]-94.8683298050514*fL[23]-73.48469228349535*(fR[6]+fL[6])+42.42640687119286*fR[1]-42.42640687119286*fL[1])*amax-18.97366596101028*(alphaR[10]*fAvg[17]+fAvg[10]*alphaR[17]+alphaR[5]*fAvg[13]+fAvg[5]*alphaR[13])-21.21320343559643*(alphaR[8]*fAvg[12]+fAvg[8]*alphaR[12])-18.97366596101028*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11])-21.21320343559643*(alphaR[6]*fAvg[10]+fAvg[6]*alphaR[10])-18.97366596101028*(alphaR[1]*fAvg[7]+fAvg[1]*alphaR[7])-21.21320343559643*(alphaR[3]*fAvg[5]+fAvg[3]*alphaR[5]+alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4]+alphaR[0]*fAvg[1]+fAvg[0]*alphaR[1])); 
  Ghat[2] = -0.008333333333333333*((94.8683298050514*fR[24]-94.8683298050514*fL[24]-73.48469228349535*(fR[7]+fL[7])+42.42640687119286*fR[2]-42.42640687119286*fL[2])*amax-18.97366596101028*alphaR[10]*fAvg[18]-21.21320343559643*(alphaR[13]*fAvg[17]+fAvg[13]*alphaR[17])-18.97366596101028*(alphaR[6]*fAvg[14]+alphaR[4]*fAvg[12]+fAvg[4]*alphaR[12])-21.21320343559643*(alphaR[7]*fAvg[11]+fAvg[7]*alphaR[11])-21.21320343559643*(alphaR[5]*fAvg[10]+fAvg[5]*alphaR[10])-18.97366596101028*(alphaR[2]*fAvg[8]+fAvg[2]*alphaR[8])-21.21320343559643*(alphaR[3]*fAvg[6]+fAvg[3]*alphaR[6]+alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4]+alphaR[0]*fAvg[2]+fAvg[0]*alphaR[2])); 
  Ghat[3] = -0.008333333333333333*((94.8683298050514*fR[27]-94.8683298050514*fL[27]-73.48469228349535*(fR[10]+fL[10])+42.42640687119286*fR[4]-42.42640687119286*fL[4])*amax-18.97366596101028*alphaR[10]*fAvg[19]-21.21320343559643*(alphaR[12]*fAvg[18]+alphaR[11]*fAvg[17]+fAvg[11]*alphaR[17])-18.97366596101028*(alphaR[6]*fAvg[16]+alphaR[5]*fAvg[15])-21.21320343559643*(alphaR[8]*fAvg[14]+alphaR[7]*fAvg[13]+fAvg[7]*alphaR[13])-21.21320343559643*(alphaR[4]*fAvg[10]+fAvg[4]*alphaR[10])-18.97366596101028*alphaR[3]*fAvg[9]-21.21320343559643*(alphaR[2]*fAvg[6]+fAvg[2]*alphaR[6]+alphaR[1]*fAvg[5]+fAvg[1]*alphaR[5]+alphaR[0]*fAvg[3]+fAvg[0]*alphaR[3])); 
  Ghat[4] = -0.008333333333333333*((94.86832980505142*fR[34]-94.86832980505142*fL[34]-73.48469228349535*(fR[15]+fL[15])+42.42640687119286*fR[5]-42.42640687119286*fL[5])*amax+((-16.97056274847715*alphaR[17])-18.97366596101028*alphaR[6])*fAvg[18]-18.97366596101028*(alphaR[5]*fAvg[17]+fAvg[5]*alphaR[17]+alphaR[10]*(fAvg[14]+fAvg[13])+fAvg[10]*alphaR[13])+((-16.97056274847715*alphaR[11])-18.97366596101028*alphaR[2])*fAvg[12]+((-16.97056274847715*fAvg[11])-18.97366596101028*fAvg[2])*alphaR[12]-18.97366596101028*(alphaR[1]*fAvg[11]+fAvg[1]*alphaR[11])-21.21320343559643*(alphaR[3]*fAvg[10]+fAvg[3]*alphaR[10])-18.97366596101028*(alphaR[4]*fAvg[8]+fAvg[4]*alphaR[8]+alphaR[4]*fAvg[7]+fAvg[4]*alphaR[7])-21.21320343559643*(alphaR[5]*fAvg[6]+fAvg[5]*alphaR[6]+alphaR[0]*fAvg[4]+fAvg[0]*alphaR[4]+alphaR[1]*fAvg[2]+fAvg[1]*alphaR[2])); 
  Ghat[5] = -0.008333333333333333*((94.86832980505142*fR[39]-94.86832980505142*fL[39]-73.48469228349535*(fR[17]+fL[17])+42.42640687119286*fR[8]-42.42640687119286*fL[8])*amax+((-16.97056274847715*alphaR[17])-18.97366596101028*alphaR[6])*fAvg[19]-21.21320343559643*alphaR[8]*fAvg[18]-18.97366596101028*(alphaR[4]*fAvg[17]+fAvg[4]*alphaR[17]+alphaR[10]*fAvg[16])+((-16.97056274847715*alphaR[13])-18.97366596101028*alphaR[3])*fAvg[15]-21.21320343559643*alphaR[12]*fAvg[14]-18.97366596101028*(alphaR[1]*fAvg[13]+fAvg[1]*alphaR[13]+alphaR[10]*fAvg[11]+fAvg[10]*alphaR[11])-21.21320343559643*(alphaR[2]*fAvg[10]+fAvg[2]*alphaR[10])-18.97366596101028*(alphaR[5]*(fAvg[9]+fAvg[7])+fAvg[5]*alphaR[7])-21.21320343559643*(alphaR[4]*fAvg[6]+fAvg[4]*alphaR[6]+alphaR[0]*fAvg[5]+fAvg[0]*alphaR[5]+alphaR[1]*fAvg[3]+fAvg[1]*alphaR[3])); 
  Ghat[6] = -0.008333333333333333*((94.86832980505142*fR[40]-94.86832980505142*fL[40]-73.48469228349535*(fR[18]+fL[18])+42.42640687119286*fR[9]-42.42640687119286*fL[9])*amax-18.97366596101028*(alphaR[5]*fAvg[19]+alphaR[4]*fAvg[18])-21.21320343559643*(alphaR[7]*fAvg[17]+fAvg[7]*alphaR[17])-18.97366596101028*(alphaR[3]*fAvg[16]+alphaR[10]*fAvg[15]+alphaR[2]*fAvg[14])-21.21320343559643*(alphaR[11]*fAvg[13]+fAvg[11]*alphaR[13])-18.97366596101028*(alphaR[10]*fAvg[12]+fAvg[10]*alphaR[12])-21.21320343559643*(alphaR[1]*fAvg[10]+fAvg[1]*alphaR[10])-18.97366596101028*(alphaR[6]*(fAvg[9]+fAvg[8])+fAvg[6]*alphaR[8])-21.21320343559643*(alphaR[0]*fAvg[6]+fAvg[0]*alphaR[6]+alphaR[4]*fAvg[5]+fAvg[4]*alphaR[5]+alphaR[2]*fAvg[3]+fAvg[2]*alphaR[3])); 
  Ghat[7] = 0.001190476190476191*((514.3928459844675*(fR[21]+fL[21])-296.98484809835*fR[11]+296.98484809835*fL[11])*amax+(94.86832980505142*alphaR[17]+148.492424049175*alphaR[6])*fAvg[17]+148.492424049175*fAvg[6]*alphaR[17]+(94.86832980505142*alphaR[13]+148.492424049175*alphaR[3])*fAvg[13]+148.492424049175*fAvg[3]*alphaR[13]+132.815661727072*alphaR[12]*fAvg[12]+(94.86832980505142*alphaR[11]+148.492424049175*alphaR[2])*fAvg[11]+148.492424049175*fAvg[2]*alphaR[11]+132.815661727072*alphaR[10]*fAvg[10]+(94.86832980505142*alphaR[7]+148.492424049175*alphaR[0])*fAvg[7]+148.492424049175*fAvg[0]*alphaR[7]+132.815661727072*(alphaR[5]*fAvg[5]+alphaR[4]*fAvg[4]+alphaR[1]*fAvg[1])); 
  Ghat[8] = 0.001190476190476191*((514.3928459844675*(fR[22]+fL[22])-296.98484809835*fR[12]+296.98484809835*fL[12])*amax+148.492424049175*alphaR[5]*fAvg[18]+132.815661727072*alphaR[17]*fAvg[17]+148.492424049175*alphaR[3]*fAvg[14]+(94.86832980505142*alphaR[12]+148.492424049175*alphaR[1])*fAvg[12]+148.492424049175*fAvg[1]*alphaR[12]+132.815661727072*(alphaR[11]*fAvg[11]+alphaR[10]*fAvg[10])+(94.86832980505142*alphaR[8]+148.492424049175*alphaR[0])*fAvg[8]+148.492424049175*fAvg[0]*alphaR[8]+132.815661727072*(alphaR[6]*fAvg[6]+alphaR[4]*fAvg[4]+alphaR[2]*fAvg[2])); 
  Ghat[9] = 0.008333333333333333*((73.48469228349536*(fR[30]+fL[30])-42.42640687119286*fR[14]+42.42640687119286*fL[14])*amax+21.21320343559643*alphaR[4]*fAvg[19]+18.97366596101028*alphaR[17]*fAvg[17]+21.21320343559643*(alphaR[2]*fAvg[16]+alphaR[1]*fAvg[15])+18.97366596101028*(alphaR[13]*fAvg[13]+alphaR[10]*fAvg[10])+21.21320343559643*alphaR[0]*fAvg[9]+18.97366596101028*(alphaR[6]*fAvg[6]+alphaR[5]*fAvg[5]+alphaR[3]*fAvg[3])); 
  Ghat[10] = -0.001666666666666667*((474.341649025257*fR[46]-474.341649025257*fL[46]-367.4234614174767*(fR[31]+fL[31])+212.1320343559643*fR[16]-212.1320343559643*fL[16])*amax+((-84.85281374238573*alphaR[13])-94.86832980505142*alphaR[3])*fAvg[19]+((-84.85281374238573*alphaR[11])-94.86832980505142*alphaR[2])*fAvg[18]+((-84.85281374238573*alphaR[12])-94.86832980505142*alphaR[1])*fAvg[17]+((-84.85281374238573*(fAvg[15]+fAvg[12]))-94.86832980505142*fAvg[1])*alphaR[17]-94.8683298050514*(alphaR[5]*fAvg[16]+alphaR[6]*fAvg[15]+alphaR[4]*(fAvg[14]+fAvg[13])+fAvg[4]*alphaR[13]+alphaR[6]*fAvg[12]+fAvg[6]*alphaR[12]+alphaR[5]*fAvg[11]+fAvg[5]*alphaR[11])+((-94.86832980505142*(alphaR[8]+alphaR[7]))-106.0660171779821*alphaR[0])*fAvg[10]+((-94.86832980505142*(fAvg[9]+fAvg[8]+fAvg[7]))-106.0660171779821*fAvg[0])*alphaR[10]-106.0660171779821*(alphaR[1]*fAvg[6]+fAvg[1]*alphaR[6]+alphaR[2]*fAvg[5]+fAvg[2]*alphaR[5]+alphaR[3]*fAvg[4]+fAvg[3]*alphaR[4])); 
  Ghat[11] = 2.380952380952381e-4*((2571.964229922338*(fR[32]+fL[32])-1484.92424049175*fR[19]+1484.92424049175*fL[19])*amax+593.9696961967002*alphaR[10]*fAvg[18]+(474.3416490252571*alphaR[13]+742.462120245875*alphaR[3])*fAvg[17]+(664.0783086353599*fAvg[14]+474.3416490252571*fAvg[13]+742.462120245875*fAvg[3])*alphaR[17]+742.462120245875*(alphaR[6]*fAvg[13]+fAvg[6]*alphaR[13])+593.9696961967002*(alphaR[4]*fAvg[12]+fAvg[4]*alphaR[12])+(664.0783086353599*alphaR[8]+474.3416490252571*alphaR[7]+742.462120245875*alphaR[0])*fAvg[11]+(664.0783086353599*fAvg[8]+474.3416490252571*fAvg[7]+742.462120245875*fAvg[0])*alphaR[11]+664.0783086353598*(alphaR[5]*fAvg[10]+fAvg[5]*alphaR[10])+742.462120245875*(alphaR[2]*fAvg[7]+fAvg[2]*alphaR[7])+664.0783086353598*(alphaR[1]*fAvg[4]+fAvg[1]*alphaR[4])); 
  Ghat[12] = 2.380952380952381e-4*((2571.964229922338*(fR[33]+fL[33])-1484.92424049175*fR[20]+1484.92424049175*fL[20])*amax+(664.0783086353599*alphaR[13]+742.462120245875*alphaR[3])*fAvg[18]+593.9696961967002*(alphaR[10]*fAvg[17]+fAvg[10]*alphaR[17])+742.462120245875*alphaR[5]*fAvg[14]+(474.3416490252571*alphaR[8]+664.0783086353599*alphaR[7]+742.462120245875*alphaR[0])*fAvg[12]+(474.3416490252571*fAvg[8]+664.0783086353599*fAvg[7]+742.462120245875*fAvg[0])*alphaR[12]+593.9696961967002*(alphaR[4]*fAvg[11]+fAvg[4]*alphaR[11])+664.0783086353598*(alphaR[6]*fAvg[10]+fAvg[6]*alphaR[10])+742.462120245875*(alphaR[1]*fAvg[8]+fAvg[1]*alphaR[8])+664.0783086353598*(alphaR[2]*fAvg[4]+fAvg[2]*alphaR[4])); 
  Ghat[13] = 2.380952380952381e-4*((2571.964229922338*(fR[37]+fL[37])-1484.92424049175*fR[25]+1484.92424049175*fL[25])*amax+593.9696961967002*alphaR[10]*fAvg[19]+664.0783086353599*alphaR[12]*fAvg[18]+(474.3416490252571*alphaR[11]+742.462120245875*alphaR[2])*fAvg[17]+(664.0783086353599*fAvg[16]+474.3416490252571*fAvg[11]+742.462120245875*fAvg[2])*alphaR[17]+593.9696961967002*alphaR[5]*fAvg[15]+(474.3416490252571*alphaR[7]+742.462120245875*alphaR[0])*fAvg[13]+(664.0783086353599*fAvg[9]+474.3416490252571*fAvg[7]+742.462120245875*fAvg[0])*alphaR[13]+742.462120245875*(alphaR[6]*fAvg[11]+fAvg[6]*alphaR[11])+664.0783086353598*(alphaR[4]*fAvg[10]+fAvg[4]*alphaR[10])+742.462120245875*(alphaR[3]*fAvg[7]+fAvg[3]*alphaR[7])+664.0783086353598*(alphaR[1]*fAvg[5]+fAvg[1]*alphaR[5])); 
  Ghat[14] = 2.380952380952381e-4*((2571.964229922338*(fR[38]+fL[38])-1484.92424049175*fR[26]+1484.92424049175*fL[26])*amax+593.9696961967002*alphaR[10]*fAvg[19]+(474.3416490252571*alphaR[12]+742.462120245875*alphaR[1])*fAvg[18]+664.0783086353599*(alphaR[11]*fAvg[17]+fAvg[11]*alphaR[17])+593.9696961967002*alphaR[6]*fAvg[16]+(474.3416490252571*alphaR[8]+742.462120245875*alphaR[0])*fAvg[14]+742.462120245875*(alphaR[5]*fAvg[12]+fAvg[5]*alphaR[12])+664.0783086353598*(alphaR[4]*fAvg[10]+fAvg[4]*alphaR[10])+742.462120245875*(alphaR[3]*fAvg[8]+fAvg[3]*alphaR[8])+664.0783086353598*(alphaR[2]*fAvg[6]+fAvg[2]*alphaR[6])); 
  Ghat[15] = 0.001666666666666667*((367.4234614174769*(fR[42]+fL[42])-212.1320343559643*fR[28]+212.1320343559643*fL[28])*amax+(94.86832980505142*alphaR[11]+106.0660171779822*alphaR[2])*fAvg[19]+84.85281374238573*(alphaR[10]*fAvg[17]+fAvg[10]*alphaR[17])+106.0660171779821*alphaR[4]*fAvg[16]+(94.86832980505142*alphaR[7]+106.0660171779821*alphaR[0])*fAvg[15]+84.85281374238573*(alphaR[5]*fAvg[13]+fAvg[5]*alphaR[13])+94.8683298050514*(alphaR[6]*fAvg[10]+fAvg[6]*alphaR[10])+106.0660171779822*alphaR[1]*fAvg[9]+94.8683298050514*(alphaR[3]*fAvg[5]+fAvg[3]*alphaR[5])); 
  Ghat[16] = 0.001666666666666667*((367.4234614174769*(fR[43]+fL[43])-212.1320343559643*fR[29]+212.1320343559643*fL[29])*amax+(94.86832980505142*alphaR[12]+106.0660171779822*alphaR[1])*fAvg[19]+84.85281374238573*alphaR[10]*fAvg[18]+94.86832980505142*(alphaR[13]*fAvg[17]+fAvg[13]*alphaR[17])+(94.86832980505142*alphaR[8]+106.0660171779821*alphaR[0])*fAvg[16]+106.0660171779821*alphaR[4]*fAvg[15]+84.85281374238573*alphaR[6]*fAvg[14]+94.8683298050514*(alphaR[5]*fAvg[10]+fAvg[5]*alphaR[10])+106.0660171779822*alphaR[2]*fAvg[9]+94.8683298050514*(alphaR[3]*fAvg[6]+fAvg[3]*alphaR[6])); 
  Ghat[17] = 2.380952380952381e-4*((2571.964229922338*(fR[44]+fL[44])-1484.92424049175*fR[35]+1484.92424049175*fL[35])*amax+593.9696961967002*(alphaR[5]*fAvg[19]+alphaR[4]*fAvg[18])+(664.0783086353599*alphaR[8]+474.3416490252571*alphaR[7]+742.462120245875*alphaR[0])*fAvg[17]+(664.0783086353599*(fAvg[9]+fAvg[8])+474.3416490252571*fAvg[7]+742.462120245875*fAvg[0])*alphaR[17]+664.0783086353599*alphaR[13]*fAvg[16]+593.9696961967002*alphaR[10]*fAvg[15]+664.0783086353599*alphaR[11]*fAvg[14]+(474.3416490252571*alphaR[11]+742.462120245875*alphaR[2])*fAvg[13]+(474.3416490252571*fAvg[11]+742.462120245875*fAvg[2])*alphaR[13]+593.9696961967002*(alphaR[10]*fAvg[12]+fAvg[10]*alphaR[12])+742.462120245875*(alphaR[3]*fAvg[11]+fAvg[3]*alphaR[11])+664.0783086353599*(alphaR[1]*fAvg[10]+fAvg[1]*alphaR[10])+742.462120245875*(alphaR[6]*fAvg[7]+fAvg[6]*alphaR[7])+664.0783086353599*(alphaR[4]*fAvg[5]+fAvg[4]*alphaR[5])); 
  Ghat[18] = 2.380952380952381e-4*((2571.964229922338*(fR[45]+fL[45])-1484.92424049175*fR[36]+1484.92424049175*fL[36])*amax+(531.262646908288*alphaR[17]+593.9696961967002*alphaR[6])*fAvg[19]+(474.3416490252571*alphaR[8]+664.0783086353599*alphaR[7]+742.462120245875*alphaR[0])*fAvg[18]+593.9696961967002*(alphaR[4]*fAvg[17]+fAvg[4]*alphaR[17]+alphaR[10]*fAvg[16])+(474.3416490252571*alphaR[12]+742.462120245875*alphaR[1])*fAvg[14]+664.0783086353599*(alphaR[12]*fAvg[13]+fAvg[12]*alphaR[13])+742.462120245875*(alphaR[3]*fAvg[12]+fAvg[3]*alphaR[12])+593.9696961967002*(alphaR[10]*fAvg[11]+fAvg[10]*alphaR[11])+664.0783086353599*(alphaR[2]*fAvg[10]+fAvg[2]*alphaR[10])+742.462120245875*(alphaR[5]*fAvg[8]+fAvg[5]*alphaR[8])+664.0783086353599*(alphaR[4]*fAvg[6]+fAvg[4]*alphaR[6])); 
  Ghat[19] = 0.001666666666666667*((367.4234614174769*(fR[47]+fL[47])-212.1320343559643*fR[41]+212.1320343559643*fL[41])*amax+(94.86832980505142*(alphaR[8]+alphaR[7])+106.0660171779821*alphaR[0])*fAvg[19]+(75.89466384404115*alphaR[17]+84.85281374238573*alphaR[6])*fAvg[18]+84.85281374238573*(alphaR[5]*fAvg[17]+fAvg[5]*alphaR[17])+(94.86832980505142*alphaR[12]+106.0660171779822*alphaR[1])*fAvg[16]+(94.86832980505142*alphaR[11]+106.0660171779822*alphaR[2])*fAvg[15]+84.85281374238573*(alphaR[10]*(fAvg[14]+fAvg[13])+fAvg[10]*alphaR[13])+94.86832980505142*(alphaR[3]*fAvg[10]+fAvg[3]*alphaR[10])+106.0660171779821*alphaR[4]*fAvg[9]+94.86832980505142*(alphaR[5]*fAvg[6]+fAvg[5]*alphaR[6])); 

  incr[0] = 0.7071067811865475*Ghat[0]; 
  incr[1] = 0.7071067811865475*Ghat[1]; 
  incr[2] = 0.7071067811865475*Ghat[2]; 
  incr[3] = -1.224744871391589*Ghat[0]; 
  incr[4] = 0.7071067811865475*Ghat[3]; 
  incr[5] = 0.7071067811865475*Ghat[4]; 
  incr[6] = -1.224744871391589*Ghat[1]; 
  incr[7] = -1.224744871391589*Ghat[2]; 
  incr[8] = 0.7071067811865475*Ghat[5]; 
  incr[9] = 0.7071067811865475*Ghat[6]; 
  incr[10] = -1.224744871391589*Ghat[3]; 
  incr[11] = 0.7071067811865475*Ghat[7]; 
  incr[12] = 0.7071067811865475*Ghat[8]; 
  incr[13] = 1.58113883008419*Ghat[0]; 
  incr[14] = 0.7071067811865475*Ghat[9]; 
  incr[15] = -1.224744871391589*Ghat[4]; 
  incr[16] = 0.7071067811865475*Ghat[10]; 
  incr[17] = -1.224744871391589*Ghat[5]; 
  incr[18] = -1.224744871391589*Ghat[6]; 
  incr[19] = 0.7071067811865475*Ghat[11]; 
  incr[20] = 0.7071067811865475*Ghat[12]; 
  incr[21] = -1.224744871391589*Ghat[7]; 
  incr[22] = -1.224744871391589*Ghat[8]; 
  incr[23] = 1.58113883008419*Ghat[1]; 
  incr[24] = 1.58113883008419*Ghat[2]; 
  incr[25] = 0.7071067811865475*Ghat[13]; 
  incr[26] = 0.7071067811865475*Ghat[14]; 
  incr[27] = 1.58113883008419*Ghat[3]; 
  incr[28] = 0.7071067811865475*Ghat[15]; 
  incr[29] = 0.7071067811865475*Ghat[16]; 
  incr[30] = -1.224744871391589*Ghat[9]; 
  incr[31] = -1.224744871391589*Ghat[10]; 
  incr[32] = -1.224744871391589*Ghat[11]; 
  incr[33] = -1.224744871391589*Ghat[12]; 
  incr[34] = 1.58113883008419*Ghat[4]; 
  incr[35] = 0.7071067811865475*Ghat[17]; 
  incr[36] = 0.7071067811865475*Ghat[18]; 
  incr[37] = -1.224744871391589*Ghat[13]; 
  incr[38] = -1.224744871391589*Ghat[14]; 
  incr[39] = 1.58113883008419*Ghat[5]; 
  incr[40] = 1.58113883008419*Ghat[6]; 
  incr[41] = 0.7071067811865475*Ghat[19]; 
  incr[42] = -1.224744871391589*Ghat[15]; 
  incr[43] = -1.224744871391589*Ghat[16]; 
  incr[44] = -1.224744871391589*Ghat[17]; 
  incr[45] = -1.224744871391589*Ghat[18]; 
  incr[46] = 1.58113883008419*Ghat[10]; 
  incr[47] = -1.224744871391589*Ghat[19]; 

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
  outR[16] += incr[16]*rdvpar2R; 
  outR[17] += incr[17]*rdvpar2R; 
  outR[18] += incr[18]*rdvpar2R; 
  outR[19] += incr[19]*rdvpar2R; 
  outR[20] += incr[20]*rdvpar2R; 
  outR[21] += incr[21]*rdvpar2R; 
  outR[22] += incr[22]*rdvpar2R; 
  outR[23] += incr[23]*rdvpar2R; 
  outR[24] += incr[24]*rdvpar2R; 
  outR[25] += incr[25]*rdvpar2R; 
  outR[26] += incr[26]*rdvpar2R; 
  outR[27] += incr[27]*rdvpar2R; 
  outR[28] += incr[28]*rdvpar2R; 
  outR[29] += incr[29]*rdvpar2R; 
  outR[30] += incr[30]*rdvpar2R; 
  outR[31] += incr[31]*rdvpar2R; 
  outR[32] += incr[32]*rdvpar2R; 
  outR[33] += incr[33]*rdvpar2R; 
  outR[34] += incr[34]*rdvpar2R; 
  outR[35] += incr[35]*rdvpar2R; 
  outR[36] += incr[36]*rdvpar2R; 
  outR[37] += incr[37]*rdvpar2R; 
  outR[38] += incr[38]*rdvpar2R; 
  outR[39] += incr[39]*rdvpar2R; 
  outR[40] += incr[40]*rdvpar2R; 
  outR[41] += incr[41]*rdvpar2R; 
  outR[42] += incr[42]*rdvpar2R; 
  outR[43] += incr[43]*rdvpar2R; 
  outR[44] += incr[44]*rdvpar2R; 
  outR[45] += incr[45]*rdvpar2R; 
  outR[46] += incr[46]*rdvpar2R; 
  outR[47] += incr[47]*rdvpar2R; 

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
  outL[11] += -1.0*incr[11]*rdvpar2L; 
  outL[12] += -1.0*incr[12]*rdvpar2L; 
  outL[13] += -1.0*incr[13]*rdvpar2L; 
  outL[14] += -1.0*incr[14]*rdvpar2L; 
  outL[15] += incr[15]*rdvpar2L; 
  outL[16] += -1.0*incr[16]*rdvpar2L; 
  outL[17] += incr[17]*rdvpar2L; 
  outL[18] += incr[18]*rdvpar2L; 
  outL[19] += -1.0*incr[19]*rdvpar2L; 
  outL[20] += -1.0*incr[20]*rdvpar2L; 
  outL[21] += incr[21]*rdvpar2L; 
  outL[22] += incr[22]*rdvpar2L; 
  outL[23] += -1.0*incr[23]*rdvpar2L; 
  outL[24] += -1.0*incr[24]*rdvpar2L; 
  outL[25] += -1.0*incr[25]*rdvpar2L; 
  outL[26] += -1.0*incr[26]*rdvpar2L; 
  outL[27] += -1.0*incr[27]*rdvpar2L; 
  outL[28] += -1.0*incr[28]*rdvpar2L; 
  outL[29] += -1.0*incr[29]*rdvpar2L; 
  outL[30] += incr[30]*rdvpar2L; 
  outL[31] += incr[31]*rdvpar2L; 
  outL[32] += incr[32]*rdvpar2L; 
  outL[33] += incr[33]*rdvpar2L; 
  outL[34] += -1.0*incr[34]*rdvpar2L; 
  outL[35] += -1.0*incr[35]*rdvpar2L; 
  outL[36] += -1.0*incr[36]*rdvpar2L; 
  outL[37] += incr[37]*rdvpar2L; 
  outL[38] += incr[38]*rdvpar2L; 
  outL[39] += -1.0*incr[39]*rdvpar2L; 
  outL[40] += -1.0*incr[40]*rdvpar2L; 
  outL[41] += -1.0*incr[41]*rdvpar2L; 
  outL[42] += incr[42]*rdvpar2L; 
  outL[43] += incr[43]*rdvpar2L; 
  outL[44] += incr[44]*rdvpar2L; 
  outL[45] += incr[45]*rdvpar2L; 
  outL[46] += -1.0*incr[46]*rdvpar2L; 
  outL[47] += incr[47]*rdvpar2L; 

  return std::abs(alphaSurfAvgR);
} 
