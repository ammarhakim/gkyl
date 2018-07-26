#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc2x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[6]*volFact; 
} 
void MomentCalc2x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[6]*volFact; 
  out[4] += 2.828427124746191*f[16]*volFact; 
  out[5] += 2.828427124746191*f[17]*volFact; 
  out[6] += 2.828427124746191*f[31]*volFact; 
  out[7] += 2.828427124746191*f[32]*volFact; 
} 
void MomentCalc2x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  out[0] += (2.828427124746191*f[0]*wx1+0.8164965809277261*f[3]*dv1)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1+0.8164965809277261*f[7]*dv1)*volFact; 
  out[2] += (2.828427124746191*f[2]*wx1+0.8164965809277261*f[8]*dv1)*volFact; 
  out[3] += (2.828427124746191*f[6]*wx1+0.8164965809277261*f[16]*dv1)*volFact; 
  out[4] += (2.828427124746191*f[0]*wx2+0.8164965809277261*f[4]*dv2)*volFact; 
  out[5] += (2.828427124746191*f[1]*wx2+0.8164965809277261*f[9]*dv2)*volFact; 
  out[6] += (2.828427124746191*f[2]*wx2+0.8164965809277261*f[10]*dv2)*volFact; 
  out[7] += (2.828427124746191*f[6]*wx2+0.8164965809277261*f[17]*dv2)*volFact; 
  out[8] += (2.828427124746191*f[0]*wx3+0.8164965809277261*f[5]*dv3)*volFact; 
  out[9] += (2.828427124746191*f[1]*wx3+0.8164965809277261*f[12]*dv3)*volFact; 
  out[10] += (2.828427124746191*f[2]*wx3+0.8164965809277261*f[13]*dv3)*volFact; 
  out[11] += (2.828427124746191*f[6]*wx3+0.8164965809277261*f[20]*dv3)*volFact; 
} 
void MomentCalc2x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  out[0] += (2.828427124746191*f[0]*wx1+0.8164965809277261*f[3]*dv1)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1+0.8164965809277261*f[7]*dv1)*volFact; 
  out[2] += (2.828427124746191*f[2]*wx1+0.8164965809277261*f[8]*dv1)*volFact; 
  out[3] += (2.828427124746191*f[6]*wx1+0.8164965809277261*f[21]*dv1)*volFact; 
  out[4] += (2.828427124746191*f[16]*wx1+0.816496580927726*f[33]*dv1)*volFact; 
  out[5] += (2.828427124746191*f[17]*wx1+0.816496580927726*f[34]*dv1)*volFact; 
  out[6] += (2.828427124746191*f[31]*wx1+0.816496580927726*f[56]*dv1)*volFact; 
  out[7] += (2.828427124746191*f[32]*wx1+0.816496580927726*f[57]*dv1)*volFact; 
  out[8] += (2.828427124746191*f[0]*wx2+0.8164965809277261*f[4]*dv2)*volFact; 
  out[9] += (2.828427124746191*f[1]*wx2+0.8164965809277261*f[9]*dv2)*volFact; 
  out[10] += (2.828427124746191*f[2]*wx2+0.8164965809277261*f[10]*dv2)*volFact; 
  out[11] += (2.828427124746191*f[6]*wx2+0.8164965809277261*f[22]*dv2)*volFact; 
  out[12] += (2.828427124746191*f[16]*wx2+0.816496580927726*f[37]*dv2)*volFact; 
  out[13] += (2.828427124746191*f[17]*wx2+0.816496580927726*f[38]*dv2)*volFact; 
  out[14] += (2.828427124746191*f[31]*wx2+0.816496580927726*f[59]*dv2)*volFact; 
  out[15] += (2.828427124746191*f[32]*wx2+0.816496580927726*f[60]*dv2)*volFact; 
  out[16] += (2.828427124746191*f[0]*wx3+0.8164965809277261*f[5]*dv3)*volFact; 
  out[17] += (2.828427124746191*f[1]*wx3+0.8164965809277261*f[12]*dv3)*volFact; 
  out[18] += (2.828427124746191*f[2]*wx3+0.8164965809277261*f[13]*dv3)*volFact; 
  out[19] += (2.828427124746191*f[6]*wx3+0.8164965809277261*f[25]*dv3)*volFact; 
  out[20] += (2.828427124746191*f[16]*wx3+0.816496580927726*f[43]*dv3)*volFact; 
  out[21] += (2.828427124746191*f[17]*wx3+0.816496580927726*f[44]*dv3)*volFact; 
  out[22] += (2.828427124746191*f[31]*wx3+0.816496580927726*f[68]*dv3)*volFact; 
  out[23] += (2.828427124746191*f[32]*wx3+0.816496580927726*f[69]*dv3)*volFact; 
} 
void MomentCalc2x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double zeta1[32]; 

  zeta1[0] = 5.656854249492382*wx1_sq+0.4714045207910317*dv1_sq; 
  zeta1[3] = 3.265986323710906*dv1*wx1; 
  double zeta2[32]; 

  zeta2[0] = 5.656854249492382*wx1*wx2; 
  zeta2[3] = 1.632993161855453*dv1*wx2; 
  zeta2[4] = 1.632993161855453*dv2*wx1; 
  zeta2[11] = 0.4714045207910317*dv1*dv2; 
  double zeta3[32]; 

  zeta3[0] = 5.656854249492382*wx1*wx3; 
  zeta3[3] = 1.632993161855453*dv1*wx3; 
  zeta3[5] = 1.632993161855453*dv3*wx1; 
  zeta3[14] = 0.4714045207910317*dv1*dv3; 
  double zeta4[32]; 

  zeta4[0] = 5.656854249492382*wx2_sq+0.4714045207910317*dv2_sq; 
  zeta4[4] = 3.265986323710906*dv2*wx2; 
  double zeta5[32]; 

  zeta5[0] = 5.656854249492382*wx2*wx3; 
  zeta5[4] = 1.632993161855453*dv2*wx3; 
  zeta5[5] = 1.632993161855453*dv3*wx2; 
  zeta5[15] = 0.4714045207910317*dv2*dv3; 
  double zeta6[32]; 

  zeta6[0] = 5.656854249492382*wx3_sq+0.4714045207910317*dv3_sq; 
  zeta6[5] = 3.265986323710906*dv3*wx3; 
  out[0] += (0.5*f[3]*zeta1[3]+0.5*f[0]*zeta1[0])*volFact; 
  out[1] += (0.5*zeta1[3]*f[7]+0.5*zeta1[0]*f[1])*volFact; 
  out[2] += (0.5*zeta1[3]*f[8]+0.5*zeta1[0]*f[2])*volFact; 
  out[3] += (0.5*zeta1[3]*f[16]+0.5*zeta1[0]*f[6])*volFact; 
  out[4] += (0.5*f[11]*zeta2[11]+0.5*f[4]*zeta2[4]+0.5*f[3]*zeta2[3]+0.5*f[0]*zeta2[0])*volFact; 
  out[5] += (0.5*zeta2[11]*f[18]+0.5*zeta2[4]*f[9]+0.5*zeta2[3]*f[7]+0.5*zeta2[0]*f[1])*volFact; 
  out[6] += (0.5*zeta2[11]*f[19]+0.5*zeta2[4]*f[10]+0.5*zeta2[3]*f[8]+0.5*zeta2[0]*f[2])*volFact; 
  out[7] += (0.5*zeta2[11]*f[26]+0.5*zeta2[4]*f[17]+0.5*zeta2[3]*f[16]+0.5*zeta2[0]*f[6])*volFact; 
  out[8] += (0.5*f[14]*zeta3[14]+0.5*f[5]*zeta3[5]+0.5*f[3]*zeta3[3]+0.5*f[0]*zeta3[0])*volFact; 
  out[9] += (0.5*zeta3[14]*f[21]+0.5*zeta3[5]*f[12]+0.5*zeta3[3]*f[7]+0.5*zeta3[0]*f[1])*volFact; 
  out[10] += (0.5*zeta3[14]*f[22]+0.5*zeta3[5]*f[13]+0.5*zeta3[3]*f[8]+0.5*zeta3[0]*f[2])*volFact; 
  out[11] += (0.5*zeta3[14]*f[27]+0.5*zeta3[5]*f[20]+0.5*zeta3[3]*f[16]+0.5*zeta3[0]*f[6])*volFact; 
  out[12] += (0.5*f[4]*zeta4[4]+0.5*f[0]*zeta4[0])*volFact; 
  out[13] += (0.5*zeta4[4]*f[9]+0.5*zeta4[0]*f[1])*volFact; 
  out[14] += (0.5*zeta4[4]*f[10]+0.5*zeta4[0]*f[2])*volFact; 
  out[15] += (0.5*zeta4[4]*f[17]+0.5*zeta4[0]*f[6])*volFact; 
  out[16] += (0.5*f[15]*zeta5[15]+0.5*f[5]*zeta5[5]+0.5*f[4]*zeta5[4]+0.5*f[0]*zeta5[0])*volFact; 
  out[17] += (0.5*zeta5[15]*f[23]+0.5*zeta5[5]*f[12]+0.5*zeta5[4]*f[9]+0.5*zeta5[0]*f[1])*volFact; 
  out[18] += (0.5*zeta5[15]*f[24]+0.5*zeta5[5]*f[13]+0.5*zeta5[4]*f[10]+0.5*zeta5[0]*f[2])*volFact; 
  out[19] += (0.5*zeta5[15]*f[28]+0.5*zeta5[5]*f[20]+0.5*zeta5[4]*f[17]+0.5*zeta5[0]*f[6])*volFact; 
  out[20] += (0.5*f[5]*zeta6[5]+0.5*f[0]*zeta6[0])*volFact; 
  out[21] += (0.5*zeta6[5]*f[12]+0.5*zeta6[0]*f[1])*volFact; 
  out[22] += (0.5*zeta6[5]*f[13]+0.5*zeta6[0]*f[2])*volFact; 
  out[23] += (0.5*zeta6[5]*f[20]+0.5*zeta6[0]*f[6])*volFact; 
} 
void MomentCalc2x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[3]*dv1*wx1+0.210818510677892*f[18]*dv1_sq+0.2357022603955158*f[0]*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[7]*dv1*wx1+0.2108185106778921*f[35]*dv1_sq+0.2357022603955158*f[1]*dv1_sq)*volFact; 
  out[2] += (2.828427124746191*f[2]*wx1_sq+1.632993161855453*f[8]*dv1*wx1+0.2108185106778921*f[36]*dv1_sq+0.2357022603955158*f[2]*dv1_sq)*volFact; 
  out[3] += (2.828427124746191*f[6]*wx1_sq+1.632993161855453*f[21]*dv1*wx1+0.210818510677892*f[58]*dv1_sq+0.2357022603955158*f[6]*dv1_sq)*volFact; 
  out[4] += (2.828427124746191*f[16]*wx1_sq+1.632993161855453*f[33]*dv1*wx1+0.2357022603955158*f[16]*dv1_sq)*volFact; 
  out[5] += (2.828427124746191*f[17]*wx1_sq+1.632993161855453*f[34]*dv1*wx1+0.2357022603955158*f[17]*dv1_sq)*volFact; 
  out[6] += (2.828427124746191*f[31]*wx1_sq+1.632993161855453*f[56]*dv1*wx1+0.2357022603955158*f[31]*dv1_sq)*volFact; 
  out[7] += (2.828427124746191*f[32]*wx1_sq+1.632993161855453*f[57]*dv1*wx1+0.2357022603955158*f[32]*dv1_sq)*volFact; 
  out[8] += (2.828427124746191*f[0]*wx1*wx2+0.8164965809277261*f[3]*dv1*wx2+0.8164965809277261*f[4]*dv2*wx1+0.2357022603955158*f[11]*dv1*dv2)*volFact; 
  out[9] += (2.828427124746191*f[1]*wx1*wx2+0.8164965809277261*f[7]*dv1*wx2+0.8164965809277261*f[9]*dv2*wx1+0.2357022603955158*f[23]*dv1*dv2)*volFact; 
  out[10] += (2.828427124746191*f[2]*wx1*wx2+0.8164965809277261*f[8]*dv1*wx2+0.8164965809277261*f[10]*dv2*wx1+0.2357022603955158*f[24]*dv1*dv2)*volFact; 
  out[11] += (2.828427124746191*f[6]*wx1*wx2+0.8164965809277261*f[21]*dv1*wx2+0.8164965809277261*f[22]*dv2*wx1+0.2357022603955158*f[51]*dv1*dv2)*volFact; 
  out[12] += (2.828427124746191*f[16]*wx1*wx2+0.816496580927726*f[33]*dv1*wx2+0.816496580927726*f[37]*dv2*wx1+0.2357022603955158*f[61]*dv1*dv2)*volFact; 
  out[13] += (2.828427124746191*f[17]*wx1*wx2+0.816496580927726*f[34]*dv1*wx2+0.816496580927726*f[38]*dv2*wx1+0.2357022603955158*f[62]*dv1*dv2)*volFact; 
  out[14] += (2.828427124746191*f[31]*wx1*wx2+0.816496580927726*f[56]*dv1*wx2+0.816496580927726*f[59]*dv2*wx1+0.2357022603955158*f[87]*dv1*dv2)*volFact; 
  out[15] += (2.828427124746191*f[32]*wx1*wx2+0.816496580927726*f[57]*dv1*wx2+0.816496580927726*f[60]*dv2*wx1+0.2357022603955158*f[88]*dv1*dv2)*volFact; 
  out[16] += (2.828427124746191*f[0]*wx1*wx3+0.8164965809277261*f[3]*dv1*wx3+0.8164965809277261*f[5]*dv3*wx1+0.2357022603955158*f[14]*dv1*dv3)*volFact; 
  out[17] += (2.828427124746191*f[1]*wx1*wx3+0.8164965809277261*f[7]*dv1*wx3+0.8164965809277261*f[12]*dv3*wx1+0.2357022603955158*f[26]*dv1*dv3)*volFact; 
  out[18] += (2.828427124746191*f[2]*wx1*wx3+0.8164965809277261*f[8]*dv1*wx3+0.8164965809277261*f[13]*dv3*wx1+0.2357022603955158*f[27]*dv1*dv3)*volFact; 
  out[19] += (2.828427124746191*f[6]*wx1*wx3+0.8164965809277261*f[21]*dv1*wx3+0.8164965809277261*f[25]*dv3*wx1+0.2357022603955158*f[52]*dv1*dv3)*volFact; 
  out[20] += (2.828427124746191*f[16]*wx1*wx3+0.816496580927726*f[33]*dv1*wx3+0.816496580927726*f[43]*dv3*wx1+0.2357022603955158*f[70]*dv1*dv3)*volFact; 
  out[21] += (2.828427124746191*f[17]*wx1*wx3+0.816496580927726*f[34]*dv1*wx3+0.816496580927726*f[44]*dv3*wx1+0.2357022603955158*f[71]*dv1*dv3)*volFact; 
  out[22] += (2.828427124746191*f[31]*wx1*wx3+0.816496580927726*f[56]*dv1*wx3+0.816496580927726*f[68]*dv3*wx1+0.2357022603955158*f[91]*dv1*dv3)*volFact; 
  out[23] += (2.828427124746191*f[32]*wx1*wx3+0.816496580927726*f[57]*dv1*wx3+0.816496580927726*f[69]*dv3*wx1+0.2357022603955158*f[92]*dv1*dv3)*volFact; 
  out[24] += (2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[4]*dv2*wx2+0.210818510677892*f[19]*dv2_sq+0.2357022603955158*f[0]*dv2_sq)*volFact; 
  out[25] += (2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[9]*dv2*wx2+0.2108185106778921*f[40]*dv2_sq+0.2357022603955158*f[1]*dv2_sq)*volFact; 
  out[26] += (2.828427124746191*f[2]*wx2_sq+1.632993161855453*f[10]*dv2*wx2+0.2108185106778921*f[41]*dv2_sq+0.2357022603955158*f[2]*dv2_sq)*volFact; 
  out[27] += (2.828427124746191*f[6]*wx2_sq+1.632993161855453*f[22]*dv2*wx2+0.210818510677892*f[65]*dv2_sq+0.2357022603955158*f[6]*dv2_sq)*volFact; 
  out[28] += (2.828427124746191*f[16]*wx2_sq+1.632993161855453*f[37]*dv2*wx2+0.2357022603955158*f[16]*dv2_sq)*volFact; 
  out[29] += (2.828427124746191*f[17]*wx2_sq+1.632993161855453*f[38]*dv2*wx2+0.2357022603955158*f[17]*dv2_sq)*volFact; 
  out[30] += (2.828427124746191*f[31]*wx2_sq+1.632993161855453*f[59]*dv2*wx2+0.2357022603955158*f[31]*dv2_sq)*volFact; 
  out[31] += (2.828427124746191*f[32]*wx2_sq+1.632993161855453*f[60]*dv2*wx2+0.2357022603955158*f[32]*dv2_sq)*volFact; 
  out[32] += (2.828427124746191*f[0]*wx2*wx3+0.8164965809277261*f[4]*dv2*wx3+0.8164965809277261*f[5]*dv3*wx2+0.2357022603955158*f[15]*dv2*dv3)*volFact; 
  out[33] += (2.828427124746191*f[1]*wx2*wx3+0.8164965809277261*f[9]*dv2*wx3+0.8164965809277261*f[12]*dv3*wx2+0.2357022603955158*f[28]*dv2*dv3)*volFact; 
  out[34] += (2.828427124746191*f[2]*wx2*wx3+0.8164965809277261*f[10]*dv2*wx3+0.8164965809277261*f[13]*dv3*wx2+0.2357022603955158*f[29]*dv2*dv3)*volFact; 
  out[35] += (2.828427124746191*f[6]*wx2*wx3+0.8164965809277261*f[22]*dv2*wx3+0.8164965809277261*f[25]*dv3*wx2+0.2357022603955158*f[53]*dv2*dv3)*volFact; 
  out[36] += (2.828427124746191*f[16]*wx2*wx3+0.816496580927726*f[37]*dv2*wx3+0.816496580927726*f[43]*dv3*wx2+0.2357022603955158*f[74]*dv2*dv3)*volFact; 
  out[37] += (2.828427124746191*f[17]*wx2*wx3+0.816496580927726*f[38]*dv2*wx3+0.816496580927726*f[44]*dv3*wx2+0.2357022603955158*f[75]*dv2*dv3)*volFact; 
  out[38] += (2.828427124746191*f[31]*wx2*wx3+0.816496580927726*f[59]*dv2*wx3+0.816496580927726*f[68]*dv3*wx2+0.2357022603955158*f[94]*dv2*dv3)*volFact; 
  out[39] += (2.828427124746191*f[32]*wx2*wx3+0.816496580927726*f[60]*dv2*wx3+0.816496580927726*f[69]*dv3*wx2+0.2357022603955158*f[95]*dv2*dv3)*volFact; 
  out[40] += (2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[5]*dv3*wx3+0.210818510677892*f[20]*dv3_sq+0.2357022603955158*f[0]*dv3_sq)*volFact; 
  out[41] += (2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[12]*dv3*wx3+0.2108185106778921*f[47]*dv3_sq+0.2357022603955158*f[1]*dv3_sq)*volFact; 
  out[42] += (2.828427124746191*f[2]*wx3_sq+1.632993161855453*f[13]*dv3*wx3+0.2108185106778921*f[48]*dv3_sq+0.2357022603955158*f[2]*dv3_sq)*volFact; 
  out[43] += (2.828427124746191*f[6]*wx3_sq+1.632993161855453*f[25]*dv3*wx3+0.210818510677892*f[80]*dv3_sq+0.2357022603955158*f[6]*dv3_sq)*volFact; 
  out[44] += (2.828427124746191*f[16]*wx3_sq+1.632993161855453*f[43]*dv3*wx3+0.2357022603955158*f[16]*dv3_sq)*volFact; 
  out[45] += (2.828427124746191*f[17]*wx3_sq+1.632993161855453*f[44]*dv3*wx3+0.2357022603955158*f[17]*dv3_sq)*volFact; 
  out[46] += (2.828427124746191*f[31]*wx3_sq+1.632993161855453*f[68]*dv3*wx3+0.2357022603955158*f[31]*dv3_sq)*volFact; 
  out[47] += (2.828427124746191*f[32]*wx3_sq+1.632993161855453*f[69]*dv3*wx3+0.2357022603955158*f[32]*dv3_sq)*volFact; 
} 
void MomentCalc2x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double zeta[32]; 

  zeta[0] = 5.656854249492382*wx3_sq+5.656854249492382*wx2_sq+5.656854249492382*wx1_sq+0.4714045207910317*dv3_sq+0.4714045207910317*dv2_sq+0.4714045207910317*dv1_sq; 
  zeta[3] = 3.265986323710906*dv1*wx1; 
  zeta[4] = 3.265986323710906*dv2*wx2; 
  zeta[5] = 3.265986323710906*dv3*wx3; 
  out[0] += (0.5*f[5]*zeta[5]+0.5*f[4]*zeta[4]+0.5*f[3]*zeta[3]+0.5*f[0]*zeta[0])*volFact; 
  out[1] += (0.5*zeta[5]*f[12]+0.5*zeta[4]*f[9]+0.5*zeta[3]*f[7]+0.5*zeta[0]*f[1])*volFact; 
  out[2] += (0.5*zeta[5]*f[13]+0.5*zeta[4]*f[10]+0.5*zeta[3]*f[8]+0.5*zeta[0]*f[2])*volFact; 
  out[3] += (0.5*zeta[5]*f[20]+0.5*zeta[4]*f[17]+0.5*zeta[3]*f[16]+0.5*zeta[0]*f[6])*volFact; 
} 
void MomentCalc2x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[5]*dv3*wx3+2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[4]*dv2*wx2+2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[3]*dv1*wx1+0.210818510677892*f[20]*dv3_sq+0.2357022603955158*f[0]*dv3_sq+0.210818510677892*f[19]*dv2_sq+0.2357022603955158*f[0]*dv2_sq+0.210818510677892*f[18]*dv1_sq+0.2357022603955158*f[0]*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[12]*dv3*wx3+2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[9]*dv2*wx2+2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[7]*dv1*wx1+0.2108185106778921*f[47]*dv3_sq+0.2357022603955158*f[1]*dv3_sq+0.2108185106778921*f[40]*dv2_sq+0.2357022603955158*f[1]*dv2_sq+0.2108185106778921*f[35]*dv1_sq+0.2357022603955158*f[1]*dv1_sq)*volFact; 
  out[2] += (2.828427124746191*f[2]*wx3_sq+1.632993161855453*f[13]*dv3*wx3+2.828427124746191*f[2]*wx2_sq+1.632993161855453*f[10]*dv2*wx2+2.828427124746191*f[2]*wx1_sq+1.632993161855453*f[8]*dv1*wx1+0.2108185106778921*f[48]*dv3_sq+0.2357022603955158*f[2]*dv3_sq+0.2108185106778921*f[41]*dv2_sq+0.2357022603955158*f[2]*dv2_sq+0.2108185106778921*f[36]*dv1_sq+0.2357022603955158*f[2]*dv1_sq)*volFact; 
  out[3] += (2.828427124746191*f[6]*wx3_sq+1.632993161855453*f[25]*dv3*wx3+2.828427124746191*f[6]*wx2_sq+1.632993161855453*f[22]*dv2*wx2+2.828427124746191*f[6]*wx1_sq+1.632993161855453*f[21]*dv1*wx1+0.210818510677892*f[80]*dv3_sq+0.2357022603955158*f[6]*dv3_sq+0.210818510677892*f[65]*dv2_sq+0.2357022603955158*f[6]*dv2_sq+0.210818510677892*f[58]*dv1_sq+0.2357022603955158*f[6]*dv1_sq)*volFact; 
  out[4] += (2.828427124746191*f[16]*wx3_sq+1.632993161855453*f[43]*dv3*wx3+2.828427124746191*f[16]*wx2_sq+1.632993161855453*f[37]*dv2*wx2+2.828427124746191*f[16]*wx1_sq+1.632993161855453*f[33]*dv1*wx1+0.2357022603955158*f[16]*dv3_sq+0.2357022603955158*f[16]*dv2_sq+0.2357022603955158*f[16]*dv1_sq)*volFact; 
  out[5] += (2.828427124746191*f[17]*wx3_sq+1.632993161855453*f[44]*dv3*wx3+2.828427124746191*f[17]*wx2_sq+1.632993161855453*f[38]*dv2*wx2+2.828427124746191*f[17]*wx1_sq+1.632993161855453*f[34]*dv1*wx1+0.2357022603955158*f[17]*dv3_sq+0.2357022603955158*f[17]*dv2_sq+0.2357022603955158*f[17]*dv1_sq)*volFact; 
  out[6] += (2.828427124746191*f[31]*wx3_sq+1.632993161855453*f[68]*dv3*wx3+2.828427124746191*f[31]*wx2_sq+1.632993161855453*f[59]*dv2*wx2+2.828427124746191*f[31]*wx1_sq+1.632993161855453*f[56]*dv1*wx1+0.2357022603955158*f[31]*dv3_sq+0.2357022603955158*f[31]*dv2_sq+0.2357022603955158*f[31]*dv1_sq)*volFact; 
  out[7] += (2.828427124746191*f[32]*wx3_sq+1.632993161855453*f[69]*dv3*wx3+2.828427124746191*f[32]*wx2_sq+1.632993161855453*f[60]*dv2*wx2+2.828427124746191*f[32]*wx1_sq+1.632993161855453*f[57]*dv1*wx1+0.2357022603955158*f[32]*dv3_sq+0.2357022603955158*f[32]*dv2_sq+0.2357022603955158*f[32]*dv1_sq)*volFact; 
} 
void MomentCalc2x3vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double wx3_cu = wx3*wx3*wx3, dv3_cu = dv3*dv3*dv3; 
  double zeta1[32]; 

  zeta1[0] = 5.656854249492382*wx1*wx3_sq+5.656854249492382*wx1*wx2_sq+5.656854249492382*wx1_cu+0.4714045207910317*dv3_sq*wx1+0.4714045207910317*dv2_sq*wx1+1.414213562373095*dv1_sq*wx1; 
  zeta1[3] = 1.632993161855453*dv1*wx3_sq+1.632993161855453*dv1*wx2_sq+4.898979485566357*dv1*wx1_sq+0.1360827634879543*dv1*dv3_sq+0.1360827634879543*dv1*dv2_sq+0.2449489742783178*dv1_cu; 
  zeta1[4] = 3.265986323710906*dv2*wx1*wx2; 
  zeta1[5] = 3.265986323710906*dv3*wx1*wx3; 
  zeta1[11] = 0.9428090415820636*dv1*dv2*wx2; 
  zeta1[14] = 0.9428090415820636*dv1*dv3*wx3; 
  double zeta2[32]; 

  zeta2[0] = 5.656854249492382*wx2*wx3_sq+5.656854249492382*wx2_cu+5.656854249492382*wx1_sq*wx2+0.4714045207910317*dv3_sq*wx2+1.414213562373095*dv2_sq*wx2+0.4714045207910317*dv1_sq*wx2; 
  zeta2[3] = 3.265986323710906*dv1*wx1*wx2; 
  zeta2[4] = 1.632993161855453*dv2*wx3_sq+4.898979485566357*dv2*wx2_sq+1.632993161855453*dv2*wx1_sq+0.1360827634879543*dv2*dv3_sq+0.2449489742783178*dv2_cu+0.1360827634879543*dv1_sq*dv2; 
  zeta2[5] = 3.265986323710906*dv3*wx2*wx3; 
  zeta2[11] = 0.9428090415820636*dv1*dv2*wx1; 
  zeta2[15] = 0.9428090415820636*dv2*dv3*wx3; 
  double zeta3[32]; 

  zeta3[0] = 5.656854249492382*wx3_cu+5.656854249492382*wx2_sq*wx3+5.656854249492382*wx1_sq*wx3+1.414213562373095*dv3_sq*wx3+0.4714045207910317*dv2_sq*wx3+0.4714045207910317*dv1_sq*wx3; 
  zeta3[3] = 3.265986323710906*dv1*wx1*wx3; 
  zeta3[4] = 3.265986323710906*dv2*wx2*wx3; 
  zeta3[5] = 4.898979485566357*dv3*wx3_sq+1.632993161855453*dv3*wx2_sq+1.632993161855453*dv3*wx1_sq+0.2449489742783178*dv3_cu+0.1360827634879543*dv2_sq*dv3+0.1360827634879543*dv1_sq*dv3; 
  zeta3[14] = 0.9428090415820636*dv1*dv3*wx1; 
  zeta3[15] = 0.9428090415820636*dv2*dv3*wx2; 
  out[0] += (0.5*f[14]*zeta1[14]+0.5*f[11]*zeta1[11]+0.5*f[5]*zeta1[5]+0.5*f[4]*zeta1[4]+0.5*f[3]*zeta1[3]+0.5*f[0]*zeta1[0])*volFact; 
  out[1] += (0.5*zeta1[14]*f[21]+0.5*zeta1[11]*f[18]+0.5*zeta1[5]*f[12]+0.5*zeta1[4]*f[9]+0.5*zeta1[3]*f[7]+0.5*zeta1[0]*f[1])*volFact; 
  out[2] += (0.5*zeta1[14]*f[22]+0.5*zeta1[11]*f[19]+0.5*zeta1[5]*f[13]+0.5*zeta1[4]*f[10]+0.5*zeta1[3]*f[8]+0.5*zeta1[0]*f[2])*volFact; 
  out[3] += (0.5*zeta1[14]*f[27]+0.5*zeta1[11]*f[26]+0.5*zeta1[5]*f[20]+0.5*zeta1[4]*f[17]+0.5*zeta1[3]*f[16]+0.5*zeta1[0]*f[6])*volFact; 
  out[4] += (0.5*f[15]*zeta2[15]+0.5*f[11]*zeta2[11]+0.5*f[5]*zeta2[5]+0.5*f[4]*zeta2[4]+0.5*f[3]*zeta2[3]+0.5*f[0]*zeta2[0])*volFact; 
  out[5] += (0.5*zeta2[15]*f[23]+0.5*zeta2[11]*f[18]+0.5*zeta2[5]*f[12]+0.5*zeta2[4]*f[9]+0.5*zeta2[3]*f[7]+0.5*zeta2[0]*f[1])*volFact; 
  out[6] += (0.5*zeta2[15]*f[24]+0.5*zeta2[11]*f[19]+0.5*zeta2[5]*f[13]+0.5*zeta2[4]*f[10]+0.5*zeta2[3]*f[8]+0.5*zeta2[0]*f[2])*volFact; 
  out[7] += (0.5*zeta2[15]*f[28]+0.5*zeta2[11]*f[26]+0.5*zeta2[5]*f[20]+0.5*zeta2[4]*f[17]+0.5*zeta2[3]*f[16]+0.5*zeta2[0]*f[6])*volFact; 
  out[8] += (0.5*f[15]*zeta3[15]+0.5*f[14]*zeta3[14]+0.5*f[5]*zeta3[5]+0.5*f[4]*zeta3[4]+0.5*f[3]*zeta3[3]+0.5*f[0]*zeta3[0])*volFact; 
  out[9] += (0.5*zeta3[15]*f[23]+0.5*zeta3[14]*f[21]+0.5*zeta3[5]*f[12]+0.5*zeta3[4]*f[9]+0.5*zeta3[3]*f[7]+0.5*zeta3[0]*f[1])*volFact; 
  out[10] += (0.5*zeta3[15]*f[24]+0.5*zeta3[14]*f[22]+0.5*zeta3[5]*f[13]+0.5*zeta3[4]*f[10]+0.5*zeta3[3]*f[8]+0.5*zeta3[0]*f[2])*volFact; 
  out[11] += (0.5*zeta3[15]*f[28]+0.5*zeta3[14]*f[27]+0.5*zeta3[5]*f[20]+0.5*zeta3[4]*f[17]+0.5*zeta3[3]*f[16]+0.5*zeta3[0]*f[6])*volFact; 
} 
void MomentCalc2x3vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double wx3_cu = wx3*wx3*wx3, dv3_cu = dv3*dv3*dv3; 
  double zeta1[112]; 

  zeta1[0] = 5.656854249492382*wx1*wx3_sq+5.656854249492382*wx1*wx2_sq+5.656854249492382*wx1_cu+0.4714045207910317*dv3_sq*wx1+0.4714045207910317*dv2_sq*wx1+1.414213562373095*dv1_sq*wx1; 
  zeta1[3] = 1.632993161855453*dv1*wx3_sq+1.632993161855453*dv1*wx2_sq+4.898979485566357*dv1*wx1_sq+0.1360827634879543*dv1*dv3_sq+0.1360827634879543*dv1*dv2_sq+0.2449489742783178*dv1_cu; 
  zeta1[4] = 3.265986323710906*dv2*wx1*wx2; 
  zeta1[5] = 3.265986323710906*dv3*wx1*wx3; 
  zeta1[11] = 0.9428090415820636*dv1*dv2*wx2; 
  zeta1[14] = 0.9428090415820636*dv1*dv3*wx3; 
  zeta1[18] = 1.264911064067352*dv1_sq*wx1; 
  zeta1[19] = 0.421637021355784*dv2_sq*wx1; 
  zeta1[20] = 0.421637021355784*dv3_sq*wx1; 
  zeta1[42] = 0.1217161238900369*dv1*dv2_sq; 
  zeta1[49] = 0.1217161238900369*dv1*dv3_sq; 
  double zeta2[112]; 

  zeta2[0] = 5.656854249492382*wx2*wx3_sq+5.656854249492382*wx2_cu+5.656854249492382*wx1_sq*wx2+0.4714045207910317*dv3_sq*wx2+1.414213562373095*dv2_sq*wx2+0.4714045207910317*dv1_sq*wx2; 
  zeta2[3] = 3.265986323710906*dv1*wx1*wx2; 
  zeta2[4] = 1.632993161855453*dv2*wx3_sq+4.898979485566357*dv2*wx2_sq+1.632993161855453*dv2*wx1_sq+0.1360827634879543*dv2*dv3_sq+0.2449489742783178*dv2_cu+0.1360827634879543*dv1_sq*dv2; 
  zeta2[5] = 3.265986323710906*dv3*wx2*wx3; 
  zeta2[11] = 0.9428090415820636*dv1*dv2*wx1; 
  zeta2[15] = 0.9428090415820636*dv2*dv3*wx3; 
  zeta2[18] = 0.421637021355784*dv1_sq*wx2; 
  zeta2[19] = 1.264911064067352*dv2_sq*wx2; 
  zeta2[20] = 0.421637021355784*dv3_sq*wx2; 
  zeta2[39] = 0.1217161238900369*dv1_sq*dv2; 
  zeta2[50] = 0.1217161238900369*dv2*dv3_sq; 
  double zeta3[112]; 

  zeta3[0] = 5.656854249492382*wx3_cu+5.656854249492382*wx2_sq*wx3+5.656854249492382*wx1_sq*wx3+1.414213562373095*dv3_sq*wx3+0.4714045207910317*dv2_sq*wx3+0.4714045207910317*dv1_sq*wx3; 
  zeta3[3] = 3.265986323710906*dv1*wx1*wx3; 
  zeta3[4] = 3.265986323710906*dv2*wx2*wx3; 
  zeta3[5] = 4.898979485566357*dv3*wx3_sq+1.632993161855453*dv3*wx2_sq+1.632993161855453*dv3*wx1_sq+0.2449489742783178*dv3_cu+0.1360827634879543*dv2_sq*dv3+0.1360827634879543*dv1_sq*dv3; 
  zeta3[14] = 0.9428090415820636*dv1*dv3*wx1; 
  zeta3[15] = 0.9428090415820636*dv2*dv3*wx2; 
  zeta3[18] = 0.421637021355784*dv1_sq*wx3; 
  zeta3[19] = 0.421637021355784*dv2_sq*wx3; 
  zeta3[20] = 1.264911064067352*dv3_sq*wx3; 
  zeta3[45] = 0.1217161238900369*dv1_sq*dv3; 
  zeta3[46] = 0.1217161238900369*dv2_sq*dv3; 
  out[0] += (0.5*f[49]*zeta1[49]+0.5*f[42]*zeta1[42]+0.5*f[20]*zeta1[20]+0.5*f[19]*zeta1[19]+0.5*f[18]*zeta1[18]+0.5*f[14]*zeta1[14]+0.5*f[11]*zeta1[11]+0.5*f[5]*zeta1[5]+0.5*f[4]*zeta1[4]+0.5*f[3]*zeta1[3]+0.5*f[0]*zeta1[0])*volFact; 
  out[1] += (0.5000000000000001*zeta1[49]*f[81]+0.5000000000000001*zeta1[42]*f[66]+0.5000000000000001*zeta1[20]*f[47]+0.5000000000000001*zeta1[19]*f[40]+0.5000000000000001*zeta1[18]*f[35]+0.5*zeta1[14]*f[26]+0.5*zeta1[11]*f[23]+0.5*zeta1[5]*f[12]+0.5*zeta1[4]*f[9]+0.5*zeta1[3]*f[7]+0.5*zeta1[0]*f[1])*volFact; 
  out[2] += (0.5000000000000001*zeta1[49]*f[82]+0.5000000000000001*zeta1[42]*f[67]+0.5000000000000001*zeta1[20]*f[48]+0.5000000000000001*zeta1[19]*f[41]+0.5000000000000001*zeta1[18]*f[36]+0.5*zeta1[14]*f[27]+0.5*zeta1[11]*f[24]+0.5*zeta1[5]*f[13]+0.5*zeta1[4]*f[10]+0.5*zeta1[3]*f[8]+0.5*zeta1[0]*f[2])*volFact; 
  out[3] += (0.5*zeta1[49]*f[103]+0.5*zeta1[42]*f[90]+0.5*zeta1[20]*f[80]+0.5*zeta1[19]*f[65]+0.5*zeta1[18]*f[58]+0.5*zeta1[14]*f[52]+0.5*zeta1[11]*f[51]+0.5*zeta1[5]*f[25]+0.5*zeta1[4]*f[22]+0.5*zeta1[3]*f[21]+0.5*zeta1[0]*f[6])*volFact; 
  out[4] += (0.5*zeta1[14]*f[70]+0.5*zeta1[11]*f[61]+0.5000000000000001*zeta1[5]*f[43]+0.5000000000000001*zeta1[4]*f[37]+0.5000000000000001*zeta1[3]*f[33]+0.5*zeta1[0]*f[16])*volFact; 
  out[5] += (0.5*zeta1[14]*f[71]+0.5*zeta1[11]*f[62]+0.5000000000000001*zeta1[5]*f[44]+0.5000000000000001*zeta1[4]*f[38]+0.5000000000000001*zeta1[3]*f[34]+0.5*zeta1[0]*f[17])*volFact; 
  out[6] += (0.5*zeta1[14]*f[91]+0.5*zeta1[11]*f[87]+0.5000000000000001*zeta1[5]*f[68]+0.5000000000000001*zeta1[4]*f[59]+0.5000000000000001*zeta1[3]*f[56]+0.5*zeta1[0]*f[31])*volFact; 
  out[7] += (0.5*zeta1[14]*f[92]+0.5*zeta1[11]*f[88]+0.5000000000000001*zeta1[5]*f[69]+0.5000000000000001*zeta1[4]*f[60]+0.5000000000000001*zeta1[3]*f[57]+0.5*zeta1[0]*f[32])*volFact; 
  out[8] += (0.5*f[50]*zeta2[50]+0.5*f[39]*zeta2[39]+0.5*f[20]*zeta2[20]+0.5*f[19]*zeta2[19]+0.5*f[18]*zeta2[18]+0.5*f[15]*zeta2[15]+0.5*f[11]*zeta2[11]+0.5*f[5]*zeta2[5]+0.5*f[4]*zeta2[4]+0.5*f[3]*zeta2[3]+0.5*f[0]*zeta2[0])*volFact; 
  out[9] += (0.5000000000000001*zeta2[50]*f[83]+0.5000000000000001*zeta2[39]*f[63]+0.5000000000000001*zeta2[20]*f[47]+0.5000000000000001*zeta2[19]*f[40]+0.5000000000000001*zeta2[18]*f[35]+0.5*zeta2[15]*f[28]+0.5*zeta2[11]*f[23]+0.5*zeta2[5]*f[12]+0.5*zeta2[4]*f[9]+0.5*zeta2[3]*f[7]+0.5*zeta2[0]*f[1])*volFact; 
  out[10] += (0.5000000000000001*zeta2[50]*f[84]+0.5000000000000001*zeta2[39]*f[64]+0.5000000000000001*zeta2[20]*f[48]+0.5000000000000001*zeta2[19]*f[41]+0.5000000000000001*zeta2[18]*f[36]+0.5*zeta2[15]*f[29]+0.5*zeta2[11]*f[24]+0.5*zeta2[5]*f[13]+0.5*zeta2[4]*f[10]+0.5*zeta2[3]*f[8]+0.5*zeta2[0]*f[2])*volFact; 
  out[11] += (0.5*zeta2[50]*f[104]+0.5*zeta2[39]*f[89]+0.5*zeta2[20]*f[80]+0.5*zeta2[19]*f[65]+0.5*zeta2[18]*f[58]+0.5*zeta2[15]*f[53]+0.5*zeta2[11]*f[51]+0.5*zeta2[5]*f[25]+0.5*zeta2[4]*f[22]+0.5*zeta2[3]*f[21]+0.5*zeta2[0]*f[6])*volFact; 
  out[12] += (0.5*zeta2[15]*f[74]+0.5*zeta2[11]*f[61]+0.5000000000000001*zeta2[5]*f[43]+0.5000000000000001*zeta2[4]*f[37]+0.5000000000000001*zeta2[3]*f[33]+0.5*zeta2[0]*f[16])*volFact; 
  out[13] += (0.5*zeta2[15]*f[75]+0.5*zeta2[11]*f[62]+0.5000000000000001*zeta2[5]*f[44]+0.5000000000000001*zeta2[4]*f[38]+0.5000000000000001*zeta2[3]*f[34]+0.5*zeta2[0]*f[17])*volFact; 
  out[14] += (0.5*zeta2[15]*f[94]+0.5*zeta2[11]*f[87]+0.5000000000000001*zeta2[5]*f[68]+0.5000000000000001*zeta2[4]*f[59]+0.5000000000000001*zeta2[3]*f[56]+0.5*zeta2[0]*f[31])*volFact; 
  out[15] += (0.5*zeta2[15]*f[95]+0.5*zeta2[11]*f[88]+0.5000000000000001*zeta2[5]*f[69]+0.5000000000000001*zeta2[4]*f[60]+0.5000000000000001*zeta2[3]*f[57]+0.5*zeta2[0]*f[32])*volFact; 
  out[16] += (0.5*f[46]*zeta3[46]+0.5*f[45]*zeta3[45]+0.5*f[20]*zeta3[20]+0.5*f[19]*zeta3[19]+0.5*f[18]*zeta3[18]+0.5*f[15]*zeta3[15]+0.5*f[14]*zeta3[14]+0.5*f[5]*zeta3[5]+0.5*f[4]*zeta3[4]+0.5*f[3]*zeta3[3]+0.5*f[0]*zeta3[0])*volFact; 
  out[17] += (0.5000000000000001*zeta3[46]*f[77]+0.5000000000000001*zeta3[45]*f[72]+0.5000000000000001*zeta3[20]*f[47]+0.5000000000000001*zeta3[19]*f[40]+0.5000000000000001*zeta3[18]*f[35]+0.5*zeta3[15]*f[28]+0.5*zeta3[14]*f[26]+0.5*zeta3[5]*f[12]+0.5*zeta3[4]*f[9]+0.5*zeta3[3]*f[7]+0.5*zeta3[0]*f[1])*volFact; 
  out[18] += (0.5000000000000001*zeta3[46]*f[78]+0.5000000000000001*zeta3[45]*f[73]+0.5000000000000001*zeta3[20]*f[48]+0.5000000000000001*zeta3[19]*f[41]+0.5000000000000001*zeta3[18]*f[36]+0.5*zeta3[15]*f[29]+0.5*zeta3[14]*f[27]+0.5*zeta3[5]*f[13]+0.5*zeta3[4]*f[10]+0.5*zeta3[3]*f[8]+0.5*zeta3[0]*f[2])*volFact; 
  out[19] += (0.5*zeta3[46]*f[100]+0.5*zeta3[45]*f[93]+0.5*zeta3[20]*f[80]+0.5*zeta3[19]*f[65]+0.5*zeta3[18]*f[58]+0.5*zeta3[15]*f[53]+0.5*zeta3[14]*f[52]+0.5*zeta3[5]*f[25]+0.5*zeta3[4]*f[22]+0.5*zeta3[3]*f[21]+0.5*zeta3[0]*f[6])*volFact; 
  out[20] += (0.5*zeta3[15]*f[74]+0.5*zeta3[14]*f[70]+0.5000000000000001*zeta3[5]*f[43]+0.5000000000000001*zeta3[4]*f[37]+0.5000000000000001*zeta3[3]*f[33]+0.5*zeta3[0]*f[16])*volFact; 
  out[21] += (0.5*zeta3[15]*f[75]+0.5*zeta3[14]*f[71]+0.5000000000000001*zeta3[5]*f[44]+0.5000000000000001*zeta3[4]*f[38]+0.5000000000000001*zeta3[3]*f[34]+0.5*zeta3[0]*f[17])*volFact; 
  out[22] += (0.5*zeta3[15]*f[94]+0.5*zeta3[14]*f[91]+0.5000000000000001*zeta3[5]*f[68]+0.5000000000000001*zeta3[4]*f[59]+0.5000000000000001*zeta3[3]*f[56]+0.5*zeta3[0]*f[31])*volFact; 
  out[23] += (0.5*zeta3[15]*f[95]+0.5*zeta3[14]*f[92]+0.5000000000000001*zeta3[5]*f[69]+0.5000000000000001*zeta3[4]*f[60]+0.5000000000000001*zeta3[3]*f[57]+0.5*zeta3[0]*f[32])*volFact; 
} 
void MomentCalc2x3vSer_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double tempM0[4], tempM1i[12]; 

  double zeta[32]; 

  zeta[0] = 5.656854249492382*wx3_sq+5.656854249492382*wx2_sq+5.656854249492382*wx1_sq+0.4714045207910317*dv3_sq+0.4714045207910317*dv2_sq+0.4714045207910317*dv1_sq; 
  zeta[3] = 3.265986323710906*dv1*wx1; 
  zeta[4] = 3.265986323710906*dv2*wx2; 
  zeta[5] = 3.265986323710906*dv3*wx3; 
  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 
  tempM0[2] = 2.828427124746191*f[2]*volFact; 
  tempM0[3] = 2.828427124746191*f[6]*volFact; 

  tempM1i[0] = 0.8164965809277261*f[3]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.8164965809277261*f[7]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.8164965809277261*f[8]*dv1*volFact+tempM0[2]*wx1; 
  tempM1i[3] = 0.8164965809277261*f[16]*dv1*volFact+tempM0[3]*wx1; 
  tempM1i[4] = 0.8164965809277261*f[4]*dv2*volFact+tempM0[0]*wx2; 
  tempM1i[5] = 0.8164965809277261*f[9]*dv2*volFact+tempM0[1]*wx2; 
  tempM1i[6] = 0.8164965809277261*f[10]*dv2*volFact+tempM0[2]*wx2; 
  tempM1i[7] = 0.8164965809277261*f[17]*dv2*volFact+tempM0[3]*wx2; 
  tempM1i[8] = 0.8164965809277261*f[5]*dv3*volFact+tempM0[0]*wx3; 
  tempM1i[9] = 0.8164965809277261*f[12]*dv3*volFact+tempM0[1]*wx3; 
  tempM1i[10] = 0.8164965809277261*f[13]*dv3*volFact+tempM0[2]*wx3; 
  tempM1i[11] = 0.8164965809277261*f[20]*dv3*volFact+tempM0[3]*wx3; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM0[2] += tempM0[2]; 
  outM0[3] += tempM0[3]; 
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
  outM2[0] += (0.5*f[5]*zeta[5]+0.5*f[4]*zeta[4]+0.5*f[3]*zeta[3]+0.5*f[0]*zeta[0])*volFact; 
  outM2[1] += (0.5*zeta[5]*f[12]+0.5*zeta[4]*f[9]+0.5*zeta[3]*f[7]+0.5*zeta[0]*f[1])*volFact; 
  outM2[2] += (0.5*zeta[5]*f[13]+0.5*zeta[4]*f[10]+0.5*zeta[3]*f[8]+0.5*zeta[0]*f[2])*volFact; 
  outM2[3] += (0.5*zeta[5]*f[20]+0.5*zeta[4]*f[17]+0.5*zeta[3]*f[16]+0.5*zeta[0]*f[6])*volFact; 
} 
void MomentCalc2x3vSer_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double tempM0[8], tempM1i[24]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 
  tempM0[2] = 2.828427124746191*f[2]*volFact; 
  tempM0[3] = 2.828427124746191*f[6]*volFact; 
  tempM0[4] = 2.828427124746191*f[16]*volFact; 
  tempM0[5] = 2.828427124746191*f[17]*volFact; 
  tempM0[6] = 2.828427124746191*f[31]*volFact; 
  tempM0[7] = 2.828427124746191*f[32]*volFact; 

  tempM1i[0] = 0.8164965809277261*f[3]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.8164965809277261*f[7]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.8164965809277261*f[8]*dv1*volFact+tempM0[2]*wx1; 
  tempM1i[3] = 0.8164965809277261*f[21]*dv1*volFact+tempM0[3]*wx1; 
  tempM1i[4] = 0.816496580927726*f[33]*dv1*volFact+tempM0[4]*wx1; 
  tempM1i[5] = 0.816496580927726*f[34]*dv1*volFact+tempM0[5]*wx1; 
  tempM1i[6] = 0.816496580927726*f[56]*dv1*volFact+tempM0[6]*wx1; 
  tempM1i[7] = 0.816496580927726*f[57]*dv1*volFact+tempM0[7]*wx1; 
  tempM1i[8] = 0.8164965809277261*f[4]*dv2*volFact+tempM0[0]*wx2; 
  tempM1i[9] = 0.8164965809277261*f[9]*dv2*volFact+tempM0[1]*wx2; 
  tempM1i[10] = 0.8164965809277261*f[10]*dv2*volFact+tempM0[2]*wx2; 
  tempM1i[11] = 0.8164965809277261*f[22]*dv2*volFact+tempM0[3]*wx2; 
  tempM1i[12] = 0.816496580927726*f[37]*dv2*volFact+tempM0[4]*wx2; 
  tempM1i[13] = 0.816496580927726*f[38]*dv2*volFact+tempM0[5]*wx2; 
  tempM1i[14] = 0.816496580927726*f[59]*dv2*volFact+tempM0[6]*wx2; 
  tempM1i[15] = 0.816496580927726*f[60]*dv2*volFact+tempM0[7]*wx2; 
  tempM1i[16] = 0.8164965809277261*f[5]*dv3*volFact+tempM0[0]*wx3; 
  tempM1i[17] = 0.8164965809277261*f[12]*dv3*volFact+tempM0[1]*wx3; 
  tempM1i[18] = 0.8164965809277261*f[13]*dv3*volFact+tempM0[2]*wx3; 
  tempM1i[19] = 0.8164965809277261*f[25]*dv3*volFact+tempM0[3]*wx3; 
  tempM1i[20] = 0.816496580927726*f[43]*dv3*volFact+tempM0[4]*wx3; 
  tempM1i[21] = 0.816496580927726*f[44]*dv3*volFact+tempM0[5]*wx3; 
  tempM1i[22] = 0.816496580927726*f[68]*dv3*volFact+tempM0[6]*wx3; 
  tempM1i[23] = 0.816496580927726*f[69]*dv3*volFact+tempM0[7]*wx3; 

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
  outM2[0] += (0.210818510677892*f[20]*dv3_sq+0.2357022603955158*f[0]*dv3_sq+0.210818510677892*f[19]*dv2_sq+0.2357022603955158*f[0]*dv2_sq+0.210818510677892*f[18]*dv1_sq+0.2357022603955158*f[0]*dv1_sq)*volFact+tempM0[0]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[16]*wx3+2.0*tempM1i[8]*wx2+2.0*tempM1i[0]*wx1; 
  outM2[1] += (0.2108185106778921*f[47]*dv3_sq+0.2357022603955158*f[1]*dv3_sq+0.2108185106778921*f[40]*dv2_sq+0.2357022603955158*f[1]*dv2_sq+0.2108185106778921*f[35]*dv1_sq+0.2357022603955158*f[1]*dv1_sq)*volFact+tempM0[1]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[17]*wx3+2.0*tempM1i[9]*wx2+2.0*tempM1i[1]*wx1; 
  outM2[2] += (0.2108185106778921*f[48]*dv3_sq+0.2357022603955158*f[2]*dv3_sq+0.2108185106778921*f[41]*dv2_sq+0.2357022603955158*f[2]*dv2_sq+0.2108185106778921*f[36]*dv1_sq+0.2357022603955158*f[2]*dv1_sq)*volFact+tempM0[2]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[18]*wx3+2.0*tempM1i[10]*wx2+2.0*tempM1i[2]*wx1; 
  outM2[3] += (0.210818510677892*f[80]*dv3_sq+0.2357022603955158*f[6]*dv3_sq+0.210818510677892*f[65]*dv2_sq+0.2357022603955158*f[6]*dv2_sq+0.210818510677892*f[58]*dv1_sq+0.2357022603955158*f[6]*dv1_sq)*volFact+tempM0[3]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[19]*wx3+2.0*tempM1i[11]*wx2+2.0*tempM1i[3]*wx1; 
  outM2[4] += (0.2357022603955158*f[16]*dv3_sq+0.2357022603955158*f[16]*dv2_sq+0.2357022603955158*f[16]*dv1_sq)*volFact+tempM0[4]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[20]*wx3+2.0*tempM1i[12]*wx2+2.0*tempM1i[4]*wx1; 
  outM2[5] += (0.2357022603955158*f[17]*dv3_sq+0.2357022603955158*f[17]*dv2_sq+0.2357022603955158*f[17]*dv1_sq)*volFact+tempM0[5]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[21]*wx3+2.0*tempM1i[13]*wx2+2.0*tempM1i[5]*wx1; 
  outM2[6] += (0.2357022603955158*f[31]*dv3_sq+0.2357022603955158*f[31]*dv2_sq+0.2357022603955158*f[31]*dv1_sq)*volFact+tempM0[6]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[22]*wx3+2.0*tempM1i[14]*wx2+2.0*tempM1i[6]*wx1; 
  outM2[7] += (0.2357022603955158*f[32]*dv3_sq+0.2357022603955158*f[32]*dv2_sq+0.2357022603955158*f[32]*dv1_sq)*volFact+tempM0[7]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[23]*wx3+2.0*tempM1i[15]*wx2+2.0*tempM1i[7]*wx1; 
} 
