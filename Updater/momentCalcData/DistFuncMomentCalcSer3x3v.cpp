#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc3x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[3]*volFact; 
  out[4] += 2.828427124746191*f[7]*volFact; 
  out[5] += 2.828427124746191*f[8]*volFact; 
  out[6] += 2.828427124746191*f[9]*volFact; 
  out[7] += 2.828427124746191*f[22]*volFact; 
} 
void MomentCalc3x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  out[0] += (2.828427124746191*f[0]*wx1+0.8164965809277261*f[4]*dv1)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1+0.8164965809277261*f[10]*dv1)*volFact; 
  out[2] += (2.828427124746191*f[2]*wx1+0.8164965809277261*f[11]*dv1)*volFact; 
  out[3] += (2.828427124746191*f[3]*wx1+0.8164965809277261*f[12]*dv1)*volFact; 
  out[4] += (2.828427124746191*f[7]*wx1+0.8164965809277261*f[23]*dv1)*volFact; 
  out[5] += (2.828427124746191*f[8]*wx1+0.8164965809277261*f[24]*dv1)*volFact; 
  out[6] += (2.828427124746191*f[9]*wx1+0.8164965809277261*f[25]*dv1)*volFact; 
  out[7] += (2.828427124746191*f[22]*wx1+0.8164965809277261*f[42]*dv1)*volFact; 
  out[8] += (2.828427124746191*f[0]*wx2+0.8164965809277261*f[5]*dv2)*volFact; 
  out[9] += (2.828427124746191*f[1]*wx2+0.8164965809277261*f[13]*dv2)*volFact; 
  out[10] += (2.828427124746191*f[2]*wx2+0.8164965809277261*f[14]*dv2)*volFact; 
  out[11] += (2.828427124746191*f[3]*wx2+0.8164965809277261*f[15]*dv2)*volFact; 
  out[12] += (2.828427124746191*f[7]*wx2+0.8164965809277261*f[26]*dv2)*volFact; 
  out[13] += (2.828427124746191*f[8]*wx2+0.8164965809277261*f[27]*dv2)*volFact; 
  out[14] += (2.828427124746191*f[9]*wx2+0.8164965809277261*f[28]*dv2)*volFact; 
  out[15] += (2.828427124746191*f[22]*wx2+0.8164965809277261*f[43]*dv2)*volFact; 
  out[16] += (2.828427124746191*f[0]*wx3+0.8164965809277261*f[6]*dv3)*volFact; 
  out[17] += (2.828427124746191*f[1]*wx3+0.8164965809277261*f[17]*dv3)*volFact; 
  out[18] += (2.828427124746191*f[2]*wx3+0.8164965809277261*f[18]*dv3)*volFact; 
  out[19] += (2.828427124746191*f[3]*wx3+0.8164965809277261*f[19]*dv3)*volFact; 
  out[20] += (2.828427124746191*f[7]*wx3+0.8164965809277261*f[32]*dv3)*volFact; 
  out[21] += (2.828427124746191*f[8]*wx3+0.8164965809277261*f[33]*dv3)*volFact; 
  out[22] += (2.828427124746191*f[9]*wx3+0.8164965809277261*f[34]*dv3)*volFact; 
  out[23] += (2.828427124746191*f[22]*wx3+0.8164965809277261*f[47]*dv3)*volFact; 
} 
void MomentCalc3x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[4]*dv1*wx1+0.2357022603955158*f[0]*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[10]*dv1*wx1+0.2357022603955158*f[1]*dv1_sq)*volFact; 
  out[2] += (2.828427124746191*f[2]*wx1_sq+1.632993161855453*f[11]*dv1*wx1+0.2357022603955158*f[2]*dv1_sq)*volFact; 
  out[3] += (2.828427124746191*f[3]*wx1_sq+1.632993161855453*f[12]*dv1*wx1+0.2357022603955158*f[3]*dv1_sq)*volFact; 
  out[4] += (2.828427124746191*f[7]*wx1_sq+1.632993161855453*f[23]*dv1*wx1+0.2357022603955158*f[7]*dv1_sq)*volFact; 
  out[5] += (2.828427124746191*f[8]*wx1_sq+1.632993161855453*f[24]*dv1*wx1+0.2357022603955158*f[8]*dv1_sq)*volFact; 
  out[6] += (2.828427124746191*f[9]*wx1_sq+1.632993161855453*f[25]*dv1*wx1+0.2357022603955158*f[9]*dv1_sq)*volFact; 
  out[7] += (2.828427124746191*f[22]*wx1_sq+1.632993161855453*f[42]*dv1*wx1+0.2357022603955158*f[22]*dv1_sq)*volFact; 
  out[8] += (2.828427124746191*f[0]*wx1*wx2+0.8164965809277261*f[4]*dv1*wx2+0.8164965809277261*f[5]*dv2*wx1+0.2357022603955158*f[16]*dv1*dv2)*volFact; 
  out[9] += (2.828427124746191*f[1]*wx1*wx2+0.8164965809277261*f[10]*dv1*wx2+0.8164965809277261*f[13]*dv2*wx1+0.2357022603955158*f[29]*dv1*dv2)*volFact; 
  out[10] += (2.828427124746191*f[2]*wx1*wx2+0.8164965809277261*f[11]*dv1*wx2+0.8164965809277261*f[14]*dv2*wx1+0.2357022603955158*f[30]*dv1*dv2)*volFact; 
  out[11] += (2.828427124746191*f[3]*wx1*wx2+0.8164965809277261*f[12]*dv1*wx2+0.8164965809277261*f[15]*dv2*wx1+0.2357022603955158*f[31]*dv1*dv2)*volFact; 
  out[12] += (2.828427124746191*f[7]*wx1*wx2+0.8164965809277261*f[23]*dv1*wx2+0.8164965809277261*f[26]*dv2*wx1+0.2357022603955158*f[44]*dv1*dv2)*volFact; 
  out[13] += (2.828427124746191*f[8]*wx1*wx2+0.8164965809277261*f[24]*dv1*wx2+0.8164965809277261*f[27]*dv2*wx1+0.2357022603955158*f[45]*dv1*dv2)*volFact; 
  out[14] += (2.828427124746191*f[9]*wx1*wx2+0.8164965809277261*f[25]*dv1*wx2+0.8164965809277261*f[28]*dv2*wx1+0.2357022603955158*f[46]*dv1*dv2)*volFact; 
  out[15] += (2.828427124746191*f[22]*wx1*wx2+0.8164965809277261*f[42]*dv1*wx2+0.8164965809277261*f[43]*dv2*wx1+0.2357022603955158*f[57]*dv1*dv2)*volFact; 
  out[16] += (2.828427124746191*f[0]*wx1*wx3+0.8164965809277261*f[4]*dv1*wx3+0.8164965809277261*f[6]*dv3*wx1+0.2357022603955158*f[20]*dv1*dv3)*volFact; 
  out[17] += (2.828427124746191*f[1]*wx1*wx3+0.8164965809277261*f[10]*dv1*wx3+0.8164965809277261*f[17]*dv3*wx1+0.2357022603955158*f[35]*dv1*dv3)*volFact; 
  out[18] += (2.828427124746191*f[2]*wx1*wx3+0.8164965809277261*f[11]*dv1*wx3+0.8164965809277261*f[18]*dv3*wx1+0.2357022603955158*f[36]*dv1*dv3)*volFact; 
  out[19] += (2.828427124746191*f[3]*wx1*wx3+0.8164965809277261*f[12]*dv1*wx3+0.8164965809277261*f[19]*dv3*wx1+0.2357022603955158*f[37]*dv1*dv3)*volFact; 
  out[20] += (2.828427124746191*f[7]*wx1*wx3+0.8164965809277261*f[23]*dv1*wx3+0.8164965809277261*f[32]*dv3*wx1+0.2357022603955158*f[48]*dv1*dv3)*volFact; 
  out[21] += (2.828427124746191*f[8]*wx1*wx3+0.8164965809277261*f[24]*dv1*wx3+0.8164965809277261*f[33]*dv3*wx1+0.2357022603955158*f[49]*dv1*dv3)*volFact; 
  out[22] += (2.828427124746191*f[9]*wx1*wx3+0.8164965809277261*f[25]*dv1*wx3+0.8164965809277261*f[34]*dv3*wx1+0.2357022603955158*f[50]*dv1*dv3)*volFact; 
  out[23] += (2.828427124746191*f[22]*wx1*wx3+0.8164965809277261*f[42]*dv1*wx3+0.8164965809277261*f[47]*dv3*wx1+0.2357022603955158*f[58]*dv1*dv3)*volFact; 
  out[24] += (2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[5]*dv2*wx2+0.2357022603955158*f[0]*dv2_sq)*volFact; 
  out[25] += (2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[13]*dv2*wx2+0.2357022603955158*f[1]*dv2_sq)*volFact; 
  out[26] += (2.828427124746191*f[2]*wx2_sq+1.632993161855453*f[14]*dv2*wx2+0.2357022603955158*f[2]*dv2_sq)*volFact; 
  out[27] += (2.828427124746191*f[3]*wx2_sq+1.632993161855453*f[15]*dv2*wx2+0.2357022603955158*f[3]*dv2_sq)*volFact; 
  out[28] += (2.828427124746191*f[7]*wx2_sq+1.632993161855453*f[26]*dv2*wx2+0.2357022603955158*f[7]*dv2_sq)*volFact; 
  out[29] += (2.828427124746191*f[8]*wx2_sq+1.632993161855453*f[27]*dv2*wx2+0.2357022603955158*f[8]*dv2_sq)*volFact; 
  out[30] += (2.828427124746191*f[9]*wx2_sq+1.632993161855453*f[28]*dv2*wx2+0.2357022603955158*f[9]*dv2_sq)*volFact; 
  out[31] += (2.828427124746191*f[22]*wx2_sq+1.632993161855453*f[43]*dv2*wx2+0.2357022603955158*f[22]*dv2_sq)*volFact; 
  out[32] += (2.828427124746191*f[0]*wx2*wx3+0.8164965809277261*f[5]*dv2*wx3+0.8164965809277261*f[6]*dv3*wx2+0.2357022603955158*f[21]*dv2*dv3)*volFact; 
  out[33] += (2.828427124746191*f[1]*wx2*wx3+0.8164965809277261*f[13]*dv2*wx3+0.8164965809277261*f[17]*dv3*wx2+0.2357022603955158*f[38]*dv2*dv3)*volFact; 
  out[34] += (2.828427124746191*f[2]*wx2*wx3+0.8164965809277261*f[14]*dv2*wx3+0.8164965809277261*f[18]*dv3*wx2+0.2357022603955158*f[39]*dv2*dv3)*volFact; 
  out[35] += (2.828427124746191*f[3]*wx2*wx3+0.8164965809277261*f[15]*dv2*wx3+0.8164965809277261*f[19]*dv3*wx2+0.2357022603955158*f[40]*dv2*dv3)*volFact; 
  out[36] += (2.828427124746191*f[7]*wx2*wx3+0.8164965809277261*f[26]*dv2*wx3+0.8164965809277261*f[32]*dv3*wx2+0.2357022603955158*f[51]*dv2*dv3)*volFact; 
  out[37] += (2.828427124746191*f[8]*wx2*wx3+0.8164965809277261*f[27]*dv2*wx3+0.8164965809277261*f[33]*dv3*wx2+0.2357022603955158*f[52]*dv2*dv3)*volFact; 
  out[38] += (2.828427124746191*f[9]*wx2*wx3+0.8164965809277261*f[28]*dv2*wx3+0.8164965809277261*f[34]*dv3*wx2+0.2357022603955158*f[53]*dv2*dv3)*volFact; 
  out[39] += (2.828427124746191*f[22]*wx2*wx3+0.8164965809277261*f[43]*dv2*wx3+0.8164965809277261*f[47]*dv3*wx2+0.2357022603955158*f[59]*dv2*dv3)*volFact; 
  out[40] += (2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[6]*dv3*wx3+0.2357022603955158*f[0]*dv3_sq)*volFact; 
  out[41] += (2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[17]*dv3*wx3+0.2357022603955158*f[1]*dv3_sq)*volFact; 
  out[42] += (2.828427124746191*f[2]*wx3_sq+1.632993161855453*f[18]*dv3*wx3+0.2357022603955158*f[2]*dv3_sq)*volFact; 
  out[43] += (2.828427124746191*f[3]*wx3_sq+1.632993161855453*f[19]*dv3*wx3+0.2357022603955158*f[3]*dv3_sq)*volFact; 
  out[44] += (2.828427124746191*f[7]*wx3_sq+1.632993161855453*f[32]*dv3*wx3+0.2357022603955158*f[7]*dv3_sq)*volFact; 
  out[45] += (2.828427124746191*f[8]*wx3_sq+1.632993161855453*f[33]*dv3*wx3+0.2357022603955158*f[8]*dv3_sq)*volFact; 
  out[46] += (2.828427124746191*f[9]*wx3_sq+1.632993161855453*f[34]*dv3*wx3+0.2357022603955158*f[9]*dv3_sq)*volFact; 
  out[47] += (2.828427124746191*f[22]*wx3_sq+1.632993161855453*f[47]*dv3*wx3+0.2357022603955158*f[22]*dv3_sq)*volFact; 
} 
void MomentCalc3x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[6]*dv3*wx3+2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[5]*dv2*wx2+2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[4]*dv1*wx1+0.2357022603955158*f[0]*dv3_sq+0.2357022603955158*f[0]*dv2_sq+0.2357022603955158*f[0]*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[17]*dv3*wx3+2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[13]*dv2*wx2+2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[10]*dv1*wx1+0.2357022603955158*f[1]*dv3_sq+0.2357022603955158*f[1]*dv2_sq+0.2357022603955158*f[1]*dv1_sq)*volFact; 
  out[2] += (2.828427124746191*f[2]*wx3_sq+1.632993161855453*f[18]*dv3*wx3+2.828427124746191*f[2]*wx2_sq+1.632993161855453*f[14]*dv2*wx2+2.828427124746191*f[2]*wx1_sq+1.632993161855453*f[11]*dv1*wx1+0.2357022603955158*f[2]*dv3_sq+0.2357022603955158*f[2]*dv2_sq+0.2357022603955158*f[2]*dv1_sq)*volFact; 
  out[3] += (2.828427124746191*f[3]*wx3_sq+1.632993161855453*f[19]*dv3*wx3+2.828427124746191*f[3]*wx2_sq+1.632993161855453*f[15]*dv2*wx2+2.828427124746191*f[3]*wx1_sq+1.632993161855453*f[12]*dv1*wx1+0.2357022603955158*f[3]*dv3_sq+0.2357022603955158*f[3]*dv2_sq+0.2357022603955158*f[3]*dv1_sq)*volFact; 
  out[4] += (2.828427124746191*f[7]*wx3_sq+1.632993161855453*f[32]*dv3*wx3+2.828427124746191*f[7]*wx2_sq+1.632993161855453*f[26]*dv2*wx2+2.828427124746191*f[7]*wx1_sq+1.632993161855453*f[23]*dv1*wx1+0.2357022603955158*f[7]*dv3_sq+0.2357022603955158*f[7]*dv2_sq+0.2357022603955158*f[7]*dv1_sq)*volFact; 
  out[5] += (2.828427124746191*f[8]*wx3_sq+1.632993161855453*f[33]*dv3*wx3+2.828427124746191*f[8]*wx2_sq+1.632993161855453*f[27]*dv2*wx2+2.828427124746191*f[8]*wx1_sq+1.632993161855453*f[24]*dv1*wx1+0.2357022603955158*f[8]*dv3_sq+0.2357022603955158*f[8]*dv2_sq+0.2357022603955158*f[8]*dv1_sq)*volFact; 
  out[6] += (2.828427124746191*f[9]*wx3_sq+1.632993161855453*f[34]*dv3*wx3+2.828427124746191*f[9]*wx2_sq+1.632993161855453*f[28]*dv2*wx2+2.828427124746191*f[9]*wx1_sq+1.632993161855453*f[25]*dv1*wx1+0.2357022603955158*f[9]*dv3_sq+0.2357022603955158*f[9]*dv2_sq+0.2357022603955158*f[9]*dv1_sq)*volFact; 
  out[7] += (2.828427124746191*f[22]*wx3_sq+1.632993161855453*f[47]*dv3*wx3+2.828427124746191*f[22]*wx2_sq+1.632993161855453*f[43]*dv2*wx2+2.828427124746191*f[22]*wx1_sq+1.632993161855453*f[42]*dv1*wx1+0.2357022603955158*f[22]*dv3_sq+0.2357022603955158*f[22]*dv2_sq+0.2357022603955158*f[22]*dv1_sq)*volFact; 
} 
void MomentCalc3x3vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double wx3_cu = wx3*wx3*wx3, dv3_cu = dv3*dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx1*wx3_sq+0.8164965809277261*f[4]*dv1*wx3_sq+1.632993161855453*f[6]*dv3*wx1*wx3+0.4714045207910317*f[20]*dv1*dv3*wx3+2.828427124746191*f[0]*wx1*wx2_sq+0.8164965809277261*f[4]*dv1*wx2_sq+1.632993161855453*f[5]*dv2*wx1*wx2+0.4714045207910317*f[16]*dv1*dv2*wx2+2.828427124746191*f[0]*wx1*wx1_sq+2.449489742783178*f[4]*dv1*wx1_sq+0.2357022603955158*f[0]*dv3_sq*wx1+0.2357022603955158*f[0]*dv2_sq*wx1+0.7071067811865475*f[0]*dv1_sq*wx1+0.06804138174397717*f[4]*dv1*dv3_sq+0.06804138174397717*f[4]*dv1*dv2_sq+0.1224744871391589*f[4]*dv1*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1*wx3_sq+0.8164965809277261*f[10]*dv1*wx3_sq+1.632993161855453*f[17]*dv3*wx1*wx3+0.4714045207910317*f[35]*dv1*dv3*wx3+2.828427124746191*f[1]*wx1*wx2_sq+0.8164965809277261*f[10]*dv1*wx2_sq+1.632993161855453*f[13]*dv2*wx1*wx2+0.4714045207910317*f[29]*dv1*dv2*wx2+2.828427124746191*f[1]*wx1*wx1_sq+2.449489742783178*f[10]*dv1*wx1_sq+0.2357022603955158*f[1]*dv3_sq*wx1+0.2357022603955158*f[1]*dv2_sq*wx1+0.7071067811865475*f[1]*dv1_sq*wx1+0.06804138174397717*f[10]*dv1*dv3_sq+0.06804138174397717*f[10]*dv1*dv2_sq+0.1224744871391589*f[10]*dv1*dv1_sq)*volFact; 
  out[2] += (2.828427124746191*f[2]*wx1*wx3_sq+0.8164965809277261*f[11]*dv1*wx3_sq+1.632993161855453*f[18]*dv3*wx1*wx3+0.4714045207910317*f[36]*dv1*dv3*wx3+2.828427124746191*f[2]*wx1*wx2_sq+0.8164965809277261*f[11]*dv1*wx2_sq+1.632993161855453*f[14]*dv2*wx1*wx2+0.4714045207910317*f[30]*dv1*dv2*wx2+2.828427124746191*f[2]*wx1*wx1_sq+2.449489742783178*f[11]*dv1*wx1_sq+0.2357022603955158*f[2]*dv3_sq*wx1+0.2357022603955158*f[2]*dv2_sq*wx1+0.7071067811865475*f[2]*dv1_sq*wx1+0.06804138174397717*f[11]*dv1*dv3_sq+0.06804138174397717*f[11]*dv1*dv2_sq+0.1224744871391589*f[11]*dv1*dv1_sq)*volFact; 
  out[3] += (2.828427124746191*f[3]*wx1*wx3_sq+0.8164965809277261*f[12]*dv1*wx3_sq+1.632993161855453*f[19]*dv3*wx1*wx3+0.4714045207910317*f[37]*dv1*dv3*wx3+2.828427124746191*f[3]*wx1*wx2_sq+0.8164965809277261*f[12]*dv1*wx2_sq+1.632993161855453*f[15]*dv2*wx1*wx2+0.4714045207910317*f[31]*dv1*dv2*wx2+2.828427124746191*f[3]*wx1*wx1_sq+2.449489742783178*f[12]*dv1*wx1_sq+0.2357022603955158*f[3]*dv3_sq*wx1+0.2357022603955158*f[3]*dv2_sq*wx1+0.7071067811865475*f[3]*dv1_sq*wx1+0.06804138174397717*f[12]*dv1*dv3_sq+0.06804138174397717*f[12]*dv1*dv2_sq+0.1224744871391589*f[12]*dv1*dv1_sq)*volFact; 
  out[4] += (2.828427124746191*f[7]*wx1*wx3_sq+0.8164965809277261*f[23]*dv1*wx3_sq+1.632993161855453*f[32]*dv3*wx1*wx3+0.4714045207910317*f[48]*dv1*dv3*wx3+2.828427124746191*f[7]*wx1*wx2_sq+0.8164965809277261*f[23]*dv1*wx2_sq+1.632993161855453*f[26]*dv2*wx1*wx2+0.4714045207910317*f[44]*dv1*dv2*wx2+2.828427124746191*f[7]*wx1*wx1_sq+2.449489742783178*f[23]*dv1*wx1_sq+0.2357022603955158*f[7]*dv3_sq*wx1+0.2357022603955158*f[7]*dv2_sq*wx1+0.7071067811865475*f[7]*dv1_sq*wx1+0.06804138174397717*f[23]*dv1*dv3_sq+0.06804138174397717*f[23]*dv1*dv2_sq+0.1224744871391589*f[23]*dv1*dv1_sq)*volFact; 
  out[5] += (2.828427124746191*f[8]*wx1*wx3_sq+0.8164965809277261*f[24]*dv1*wx3_sq+1.632993161855453*f[33]*dv3*wx1*wx3+0.4714045207910317*f[49]*dv1*dv3*wx3+2.828427124746191*f[8]*wx1*wx2_sq+0.8164965809277261*f[24]*dv1*wx2_sq+1.632993161855453*f[27]*dv2*wx1*wx2+0.4714045207910317*f[45]*dv1*dv2*wx2+2.828427124746191*f[8]*wx1*wx1_sq+2.449489742783178*f[24]*dv1*wx1_sq+0.2357022603955158*f[8]*dv3_sq*wx1+0.2357022603955158*f[8]*dv2_sq*wx1+0.7071067811865475*f[8]*dv1_sq*wx1+0.06804138174397717*f[24]*dv1*dv3_sq+0.06804138174397717*f[24]*dv1*dv2_sq+0.1224744871391589*f[24]*dv1*dv1_sq)*volFact; 
  out[6] += (2.828427124746191*f[9]*wx1*wx3_sq+0.8164965809277261*f[25]*dv1*wx3_sq+1.632993161855453*f[34]*dv3*wx1*wx3+0.4714045207910317*f[50]*dv1*dv3*wx3+2.828427124746191*f[9]*wx1*wx2_sq+0.8164965809277261*f[25]*dv1*wx2_sq+1.632993161855453*f[28]*dv2*wx1*wx2+0.4714045207910317*f[46]*dv1*dv2*wx2+2.828427124746191*f[9]*wx1*wx1_sq+2.449489742783178*f[25]*dv1*wx1_sq+0.2357022603955158*f[9]*dv3_sq*wx1+0.2357022603955158*f[9]*dv2_sq*wx1+0.7071067811865475*f[9]*dv1_sq*wx1+0.06804138174397717*f[25]*dv1*dv3_sq+0.06804138174397717*f[25]*dv1*dv2_sq+0.1224744871391589*f[25]*dv1*dv1_sq)*volFact; 
  out[7] += (2.828427124746191*f[22]*wx1*wx3_sq+0.8164965809277261*f[42]*dv1*wx3_sq+1.632993161855453*f[47]*dv3*wx1*wx3+0.4714045207910317*f[58]*dv1*dv3*wx3+2.828427124746191*f[22]*wx1*wx2_sq+0.8164965809277261*f[42]*dv1*wx2_sq+1.632993161855453*f[43]*dv2*wx1*wx2+0.4714045207910317*f[57]*dv1*dv2*wx2+2.828427124746191*f[22]*wx1*wx1_sq+2.449489742783178*f[42]*dv1*wx1_sq+0.2357022603955158*f[22]*dv3_sq*wx1+0.2357022603955158*f[22]*dv2_sq*wx1+0.7071067811865475*f[22]*dv1_sq*wx1+0.06804138174397717*f[42]*dv1*dv3_sq+0.06804138174397717*f[42]*dv1*dv2_sq+0.1224744871391589*f[42]*dv1*dv1_sq)*volFact; 
  out[8] += (2.828427124746191*f[0]*wx2*wx3_sq+0.8164965809277261*f[5]*dv2*wx3_sq+1.632993161855453*f[6]*dv3*wx2*wx3+0.4714045207910317*f[21]*dv2*dv3*wx3+2.828427124746191*f[0]*wx2*wx2_sq+2.449489742783178*f[5]*dv2*wx2_sq+2.828427124746191*f[0]*wx1_sq*wx2+1.632993161855453*f[4]*dv1*wx1*wx2+0.2357022603955158*f[0]*dv3_sq*wx2+0.7071067811865475*f[0]*dv2_sq*wx2+0.2357022603955158*f[0]*dv1_sq*wx2+0.8164965809277261*f[5]*dv2*wx1_sq+0.4714045207910317*f[16]*dv1*dv2*wx1+0.06804138174397717*f[5]*dv2*dv3_sq+0.1224744871391589*f[5]*dv2*dv2_sq+0.06804138174397717*f[5]*dv1_sq*dv2)*volFact; 
  out[9] += (2.828427124746191*f[1]*wx2*wx3_sq+0.8164965809277261*f[13]*dv2*wx3_sq+1.632993161855453*f[17]*dv3*wx2*wx3+0.4714045207910317*f[38]*dv2*dv3*wx3+2.828427124746191*f[1]*wx2*wx2_sq+2.449489742783178*f[13]*dv2*wx2_sq+2.828427124746191*f[1]*wx1_sq*wx2+1.632993161855453*f[10]*dv1*wx1*wx2+0.2357022603955158*f[1]*dv3_sq*wx2+0.7071067811865475*f[1]*dv2_sq*wx2+0.2357022603955158*f[1]*dv1_sq*wx2+0.8164965809277261*f[13]*dv2*wx1_sq+0.4714045207910317*f[29]*dv1*dv2*wx1+0.06804138174397717*f[13]*dv2*dv3_sq+0.1224744871391589*f[13]*dv2*dv2_sq+0.06804138174397717*f[13]*dv1_sq*dv2)*volFact; 
  out[10] += (2.828427124746191*f[2]*wx2*wx3_sq+0.8164965809277261*f[14]*dv2*wx3_sq+1.632993161855453*f[18]*dv3*wx2*wx3+0.4714045207910317*f[39]*dv2*dv3*wx3+2.828427124746191*f[2]*wx2*wx2_sq+2.449489742783178*f[14]*dv2*wx2_sq+2.828427124746191*f[2]*wx1_sq*wx2+1.632993161855453*f[11]*dv1*wx1*wx2+0.2357022603955158*f[2]*dv3_sq*wx2+0.7071067811865475*f[2]*dv2_sq*wx2+0.2357022603955158*f[2]*dv1_sq*wx2+0.8164965809277261*f[14]*dv2*wx1_sq+0.4714045207910317*f[30]*dv1*dv2*wx1+0.06804138174397717*f[14]*dv2*dv3_sq+0.1224744871391589*f[14]*dv2*dv2_sq+0.06804138174397717*f[14]*dv1_sq*dv2)*volFact; 
  out[11] += (2.828427124746191*f[3]*wx2*wx3_sq+0.8164965809277261*f[15]*dv2*wx3_sq+1.632993161855453*f[19]*dv3*wx2*wx3+0.4714045207910317*f[40]*dv2*dv3*wx3+2.828427124746191*f[3]*wx2*wx2_sq+2.449489742783178*f[15]*dv2*wx2_sq+2.828427124746191*f[3]*wx1_sq*wx2+1.632993161855453*f[12]*dv1*wx1*wx2+0.2357022603955158*f[3]*dv3_sq*wx2+0.7071067811865475*f[3]*dv2_sq*wx2+0.2357022603955158*f[3]*dv1_sq*wx2+0.8164965809277261*f[15]*dv2*wx1_sq+0.4714045207910317*f[31]*dv1*dv2*wx1+0.06804138174397717*f[15]*dv2*dv3_sq+0.1224744871391589*f[15]*dv2*dv2_sq+0.06804138174397717*f[15]*dv1_sq*dv2)*volFact; 
  out[12] += (2.828427124746191*f[7]*wx2*wx3_sq+0.8164965809277261*f[26]*dv2*wx3_sq+1.632993161855453*f[32]*dv3*wx2*wx3+0.4714045207910317*f[51]*dv2*dv3*wx3+2.828427124746191*f[7]*wx2*wx2_sq+2.449489742783178*f[26]*dv2*wx2_sq+2.828427124746191*f[7]*wx1_sq*wx2+1.632993161855453*f[23]*dv1*wx1*wx2+0.2357022603955158*f[7]*dv3_sq*wx2+0.7071067811865475*f[7]*dv2_sq*wx2+0.2357022603955158*f[7]*dv1_sq*wx2+0.8164965809277261*f[26]*dv2*wx1_sq+0.4714045207910317*f[44]*dv1*dv2*wx1+0.06804138174397717*f[26]*dv2*dv3_sq+0.1224744871391589*f[26]*dv2*dv2_sq+0.06804138174397717*f[26]*dv1_sq*dv2)*volFact; 
  out[13] += (2.828427124746191*f[8]*wx2*wx3_sq+0.8164965809277261*f[27]*dv2*wx3_sq+1.632993161855453*f[33]*dv3*wx2*wx3+0.4714045207910317*f[52]*dv2*dv3*wx3+2.828427124746191*f[8]*wx2*wx2_sq+2.449489742783178*f[27]*dv2*wx2_sq+2.828427124746191*f[8]*wx1_sq*wx2+1.632993161855453*f[24]*dv1*wx1*wx2+0.2357022603955158*f[8]*dv3_sq*wx2+0.7071067811865475*f[8]*dv2_sq*wx2+0.2357022603955158*f[8]*dv1_sq*wx2+0.8164965809277261*f[27]*dv2*wx1_sq+0.4714045207910317*f[45]*dv1*dv2*wx1+0.06804138174397717*f[27]*dv2*dv3_sq+0.1224744871391589*f[27]*dv2*dv2_sq+0.06804138174397717*f[27]*dv1_sq*dv2)*volFact; 
  out[14] += (2.828427124746191*f[9]*wx2*wx3_sq+0.8164965809277261*f[28]*dv2*wx3_sq+1.632993161855453*f[34]*dv3*wx2*wx3+0.4714045207910317*f[53]*dv2*dv3*wx3+2.828427124746191*f[9]*wx2*wx2_sq+2.449489742783178*f[28]*dv2*wx2_sq+2.828427124746191*f[9]*wx1_sq*wx2+1.632993161855453*f[25]*dv1*wx1*wx2+0.2357022603955158*f[9]*dv3_sq*wx2+0.7071067811865475*f[9]*dv2_sq*wx2+0.2357022603955158*f[9]*dv1_sq*wx2+0.8164965809277261*f[28]*dv2*wx1_sq+0.4714045207910317*f[46]*dv1*dv2*wx1+0.06804138174397717*f[28]*dv2*dv3_sq+0.1224744871391589*f[28]*dv2*dv2_sq+0.06804138174397717*f[28]*dv1_sq*dv2)*volFact; 
  out[15] += (2.828427124746191*f[22]*wx2*wx3_sq+0.8164965809277261*f[43]*dv2*wx3_sq+1.632993161855453*f[47]*dv3*wx2*wx3+0.4714045207910317*f[59]*dv2*dv3*wx3+2.828427124746191*f[22]*wx2*wx2_sq+2.449489742783178*f[43]*dv2*wx2_sq+2.828427124746191*f[22]*wx1_sq*wx2+1.632993161855453*f[42]*dv1*wx1*wx2+0.2357022603955158*f[22]*dv3_sq*wx2+0.7071067811865475*f[22]*dv2_sq*wx2+0.2357022603955158*f[22]*dv1_sq*wx2+0.8164965809277261*f[43]*dv2*wx1_sq+0.4714045207910317*f[57]*dv1*dv2*wx1+0.06804138174397717*f[43]*dv2*dv3_sq+0.1224744871391589*f[43]*dv2*dv2_sq+0.06804138174397717*f[43]*dv1_sq*dv2)*volFact; 
  out[16] += (2.828427124746191*f[0]*wx3*wx3_sq+2.449489742783178*f[6]*dv3*wx3_sq+2.828427124746191*f[0]*wx2_sq*wx3+1.632993161855453*f[5]*dv2*wx2*wx3+2.828427124746191*f[0]*wx1_sq*wx3+1.632993161855453*f[4]*dv1*wx1*wx3+0.7071067811865475*f[0]*dv3_sq*wx3+0.2357022603955158*f[0]*dv2_sq*wx3+0.2357022603955158*f[0]*dv1_sq*wx3+0.8164965809277261*f[6]*dv3*wx2_sq+0.4714045207910317*f[21]*dv2*dv3*wx2+0.8164965809277261*f[6]*dv3*wx1_sq+0.4714045207910317*f[20]*dv1*dv3*wx1+0.1224744871391589*f[6]*dv3*dv3_sq+0.06804138174397717*f[6]*dv2_sq*dv3+0.06804138174397717*f[6]*dv1_sq*dv3)*volFact; 
  out[17] += (2.828427124746191*f[1]*wx3*wx3_sq+2.449489742783178*f[17]*dv3*wx3_sq+2.828427124746191*f[1]*wx2_sq*wx3+1.632993161855453*f[13]*dv2*wx2*wx3+2.828427124746191*f[1]*wx1_sq*wx3+1.632993161855453*f[10]*dv1*wx1*wx3+0.7071067811865475*f[1]*dv3_sq*wx3+0.2357022603955158*f[1]*dv2_sq*wx3+0.2357022603955158*f[1]*dv1_sq*wx3+0.8164965809277261*f[17]*dv3*wx2_sq+0.4714045207910317*f[38]*dv2*dv3*wx2+0.8164965809277261*f[17]*dv3*wx1_sq+0.4714045207910317*f[35]*dv1*dv3*wx1+0.1224744871391589*f[17]*dv3*dv3_sq+0.06804138174397717*f[17]*dv2_sq*dv3+0.06804138174397717*f[17]*dv1_sq*dv3)*volFact; 
  out[18] += (2.828427124746191*f[2]*wx3*wx3_sq+2.449489742783178*f[18]*dv3*wx3_sq+2.828427124746191*f[2]*wx2_sq*wx3+1.632993161855453*f[14]*dv2*wx2*wx3+2.828427124746191*f[2]*wx1_sq*wx3+1.632993161855453*f[11]*dv1*wx1*wx3+0.7071067811865475*f[2]*dv3_sq*wx3+0.2357022603955158*f[2]*dv2_sq*wx3+0.2357022603955158*f[2]*dv1_sq*wx3+0.8164965809277261*f[18]*dv3*wx2_sq+0.4714045207910317*f[39]*dv2*dv3*wx2+0.8164965809277261*f[18]*dv3*wx1_sq+0.4714045207910317*f[36]*dv1*dv3*wx1+0.1224744871391589*f[18]*dv3*dv3_sq+0.06804138174397717*f[18]*dv2_sq*dv3+0.06804138174397717*f[18]*dv1_sq*dv3)*volFact; 
  out[19] += (2.828427124746191*f[3]*wx3*wx3_sq+2.449489742783178*f[19]*dv3*wx3_sq+2.828427124746191*f[3]*wx2_sq*wx3+1.632993161855453*f[15]*dv2*wx2*wx3+2.828427124746191*f[3]*wx1_sq*wx3+1.632993161855453*f[12]*dv1*wx1*wx3+0.7071067811865475*f[3]*dv3_sq*wx3+0.2357022603955158*f[3]*dv2_sq*wx3+0.2357022603955158*f[3]*dv1_sq*wx3+0.8164965809277261*f[19]*dv3*wx2_sq+0.4714045207910317*f[40]*dv2*dv3*wx2+0.8164965809277261*f[19]*dv3*wx1_sq+0.4714045207910317*f[37]*dv1*dv3*wx1+0.1224744871391589*f[19]*dv3*dv3_sq+0.06804138174397717*f[19]*dv2_sq*dv3+0.06804138174397717*f[19]*dv1_sq*dv3)*volFact; 
  out[20] += (2.828427124746191*f[7]*wx3*wx3_sq+2.449489742783178*f[32]*dv3*wx3_sq+2.828427124746191*f[7]*wx2_sq*wx3+1.632993161855453*f[26]*dv2*wx2*wx3+2.828427124746191*f[7]*wx1_sq*wx3+1.632993161855453*f[23]*dv1*wx1*wx3+0.7071067811865475*f[7]*dv3_sq*wx3+0.2357022603955158*f[7]*dv2_sq*wx3+0.2357022603955158*f[7]*dv1_sq*wx3+0.8164965809277261*f[32]*dv3*wx2_sq+0.4714045207910317*f[51]*dv2*dv3*wx2+0.8164965809277261*f[32]*dv3*wx1_sq+0.4714045207910317*f[48]*dv1*dv3*wx1+0.1224744871391589*f[32]*dv3*dv3_sq+0.06804138174397717*f[32]*dv2_sq*dv3+0.06804138174397717*f[32]*dv1_sq*dv3)*volFact; 
  out[21] += (2.828427124746191*f[8]*wx3*wx3_sq+2.449489742783178*f[33]*dv3*wx3_sq+2.828427124746191*f[8]*wx2_sq*wx3+1.632993161855453*f[27]*dv2*wx2*wx3+2.828427124746191*f[8]*wx1_sq*wx3+1.632993161855453*f[24]*dv1*wx1*wx3+0.7071067811865475*f[8]*dv3_sq*wx3+0.2357022603955158*f[8]*dv2_sq*wx3+0.2357022603955158*f[8]*dv1_sq*wx3+0.8164965809277261*f[33]*dv3*wx2_sq+0.4714045207910317*f[52]*dv2*dv3*wx2+0.8164965809277261*f[33]*dv3*wx1_sq+0.4714045207910317*f[49]*dv1*dv3*wx1+0.1224744871391589*f[33]*dv3*dv3_sq+0.06804138174397717*f[33]*dv2_sq*dv3+0.06804138174397717*f[33]*dv1_sq*dv3)*volFact; 
  out[22] += (2.828427124746191*f[9]*wx3*wx3_sq+2.449489742783178*f[34]*dv3*wx3_sq+2.828427124746191*f[9]*wx2_sq*wx3+1.632993161855453*f[28]*dv2*wx2*wx3+2.828427124746191*f[9]*wx1_sq*wx3+1.632993161855453*f[25]*dv1*wx1*wx3+0.7071067811865475*f[9]*dv3_sq*wx3+0.2357022603955158*f[9]*dv2_sq*wx3+0.2357022603955158*f[9]*dv1_sq*wx3+0.8164965809277261*f[34]*dv3*wx2_sq+0.4714045207910317*f[53]*dv2*dv3*wx2+0.8164965809277261*f[34]*dv3*wx1_sq+0.4714045207910317*f[50]*dv1*dv3*wx1+0.1224744871391589*f[34]*dv3*dv3_sq+0.06804138174397717*f[34]*dv2_sq*dv3+0.06804138174397717*f[34]*dv1_sq*dv3)*volFact; 
  out[23] += (2.828427124746191*f[22]*wx3*wx3_sq+2.449489742783178*f[47]*dv3*wx3_sq+2.828427124746191*f[22]*wx2_sq*wx3+1.632993161855453*f[43]*dv2*wx2*wx3+2.828427124746191*f[22]*wx1_sq*wx3+1.632993161855453*f[42]*dv1*wx1*wx3+0.7071067811865475*f[22]*dv3_sq*wx3+0.2357022603955158*f[22]*dv2_sq*wx3+0.2357022603955158*f[22]*dv1_sq*wx3+0.8164965809277261*f[47]*dv3*wx2_sq+0.4714045207910317*f[59]*dv2*dv3*wx2+0.8164965809277261*f[47]*dv3*wx1_sq+0.4714045207910317*f[58]*dv1*dv3*wx1+0.1224744871391589*f[47]*dv3*dv3_sq+0.06804138174397717*f[47]*dv2_sq*dv3+0.06804138174397717*f[47]*dv1_sq*dv3)*volFact; 
} 
void MomentCalc3x3vSer_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double tempM0[8], tempM1i[24]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 
  tempM0[2] = 2.828427124746191*f[2]*volFact; 
  tempM0[3] = 2.828427124746191*f[3]*volFact; 
  tempM0[4] = 2.828427124746191*f[7]*volFact; 
  tempM0[5] = 2.828427124746191*f[8]*volFact; 
  tempM0[6] = 2.828427124746191*f[9]*volFact; 
  tempM0[7] = 2.828427124746191*f[22]*volFact; 

  tempM1i[0] = 0.8164965809277261*f[4]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.8164965809277261*f[10]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.8164965809277261*f[11]*dv1*volFact+tempM0[2]*wx1; 
  tempM1i[3] = 0.8164965809277261*f[12]*dv1*volFact+tempM0[3]*wx1; 
  tempM1i[4] = 0.8164965809277261*f[23]*dv1*volFact+tempM0[4]*wx1; 
  tempM1i[5] = 0.8164965809277261*f[24]*dv1*volFact+tempM0[5]*wx1; 
  tempM1i[6] = 0.8164965809277261*f[25]*dv1*volFact+tempM0[6]*wx1; 
  tempM1i[7] = 0.8164965809277261*f[42]*dv1*volFact+tempM0[7]*wx1; 
  tempM1i[8] = 0.8164965809277261*f[5]*dv2*volFact+tempM0[0]*wx2; 
  tempM1i[9] = 0.8164965809277261*f[13]*dv2*volFact+tempM0[1]*wx2; 
  tempM1i[10] = 0.8164965809277261*f[14]*dv2*volFact+tempM0[2]*wx2; 
  tempM1i[11] = 0.8164965809277261*f[15]*dv2*volFact+tempM0[3]*wx2; 
  tempM1i[12] = 0.8164965809277261*f[26]*dv2*volFact+tempM0[4]*wx2; 
  tempM1i[13] = 0.8164965809277261*f[27]*dv2*volFact+tempM0[5]*wx2; 
  tempM1i[14] = 0.8164965809277261*f[28]*dv2*volFact+tempM0[6]*wx2; 
  tempM1i[15] = 0.8164965809277261*f[43]*dv2*volFact+tempM0[7]*wx2; 
  tempM1i[16] = 0.8164965809277261*f[6]*dv3*volFact+tempM0[0]*wx3; 
  tempM1i[17] = 0.8164965809277261*f[17]*dv3*volFact+tempM0[1]*wx3; 
  tempM1i[18] = 0.8164965809277261*f[18]*dv3*volFact+tempM0[2]*wx3; 
  tempM1i[19] = 0.8164965809277261*f[19]*dv3*volFact+tempM0[3]*wx3; 
  tempM1i[20] = 0.8164965809277261*f[32]*dv3*volFact+tempM0[4]*wx3; 
  tempM1i[21] = 0.8164965809277261*f[33]*dv3*volFact+tempM0[5]*wx3; 
  tempM1i[22] = 0.8164965809277261*f[34]*dv3*volFact+tempM0[6]*wx3; 
  tempM1i[23] = 0.8164965809277261*f[47]*dv3*volFact+tempM0[7]*wx3; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM0[2] += tempM0[2]; 
  outM0[3] += tempM0[3]; 
  outM0[4] += tempM0[4]; 
  outM0[5] += tempM0[5]; 
  outM0[6] += tempM0[6]; 
  outM0[7] += tempM0[7]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM1i[2] += tempM1i[2]; 
  outM1i[3] += tempM1i[3]; 
  outM1i[4] += tempM1i[4]; 
  outM1i[5] += tempM1i[5]; 
  outM1i[6] += tempM1i[6]; 
  outM1i[7] += tempM1i[7]; 
  outM1i[8] += tempM1i[8]; 
  outM1i[9] += tempM1i[9]; 
  outM1i[10] += tempM1i[10]; 
  outM1i[11] += tempM1i[11]; 
  outM1i[12] += tempM1i[12]; 
  outM1i[13] += tempM1i[13]; 
  outM1i[14] += tempM1i[14]; 
  outM1i[15] += tempM1i[15]; 
  outM1i[16] += tempM1i[16]; 
  outM1i[17] += tempM1i[17]; 
  outM1i[18] += tempM1i[18]; 
  outM1i[19] += tempM1i[19]; 
  outM1i[20] += tempM1i[20]; 
  outM1i[21] += tempM1i[21]; 
  outM1i[22] += tempM1i[22]; 
  outM1i[23] += tempM1i[23]; 
  outM2[0] += (0.2357022603955158*f[0]*dv3_sq+0.2357022603955158*f[0]*dv2_sq+0.2357022603955158*f[0]*dv1_sq)*volFact+tempM0[0]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[16]*wx3+2.0*tempM1i[8]*wx2+2.0*tempM1i[0]*wx1; 
  outM2[1] += (0.2357022603955158*f[1]*dv3_sq+0.2357022603955158*f[1]*dv2_sq+0.2357022603955158*f[1]*dv1_sq)*volFact+tempM0[1]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[17]*wx3+2.0*tempM1i[9]*wx2+2.0*tempM1i[1]*wx1; 
  outM2[2] += (0.2357022603955158*f[2]*dv3_sq+0.2357022603955158*f[2]*dv2_sq+0.2357022603955158*f[2]*dv1_sq)*volFact+tempM0[2]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[18]*wx3+2.0*tempM1i[10]*wx2+2.0*tempM1i[2]*wx1; 
  outM2[3] += (0.2357022603955158*f[3]*dv3_sq+0.2357022603955158*f[3]*dv2_sq+0.2357022603955158*f[3]*dv1_sq)*volFact+tempM0[3]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[19]*wx3+2.0*tempM1i[11]*wx2+2.0*tempM1i[3]*wx1; 
  outM2[4] += (0.2357022603955158*f[7]*dv3_sq+0.2357022603955158*f[7]*dv2_sq+0.2357022603955158*f[7]*dv1_sq)*volFact+tempM0[4]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[20]*wx3+2.0*tempM1i[12]*wx2+2.0*tempM1i[4]*wx1; 
  outM2[5] += (0.2357022603955158*f[8]*dv3_sq+0.2357022603955158*f[8]*dv2_sq+0.2357022603955158*f[8]*dv1_sq)*volFact+tempM0[5]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[21]*wx3+2.0*tempM1i[13]*wx2+2.0*tempM1i[5]*wx1; 
  outM2[6] += (0.2357022603955158*f[9]*dv3_sq+0.2357022603955158*f[9]*dv2_sq+0.2357022603955158*f[9]*dv1_sq)*volFact+tempM0[6]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[22]*wx3+2.0*tempM1i[14]*wx2+2.0*tempM1i[6]*wx1; 
  outM2[7] += (0.2357022603955158*f[22]*dv3_sq+0.2357022603955158*f[22]*dv2_sq+0.2357022603955158*f[22]*dv1_sq)*volFact+tempM0[7]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[23]*wx3+2.0*tempM1i[15]*wx2+2.0*tempM1i[7]*wx1; 
} 
