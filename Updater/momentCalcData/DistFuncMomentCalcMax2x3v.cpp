#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc2x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
} 
void MomentCalc2x3vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[6]*volFact; 
  out[4] += 2.828427124746191*f[16]*volFact; 
  out[5] += 2.828427124746191*f[17]*volFact; 
} 
void MomentCalc2x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  out[0] += 2.828427124746191*f[0]*volFact*wx1+0.8164965809277261*f[3]*dv1*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx1; 
  out[2] += 2.828427124746191*f[2]*volFact*wx1; 
  out[3] += 2.828427124746191*f[0]*volFact*wx2+0.8164965809277261*f[4]*dv2*volFact; 
  out[4] += 2.828427124746191*f[1]*volFact*wx2; 
  out[5] += 2.828427124746191*f[2]*volFact*wx2; 
  out[6] += 2.828427124746191*f[0]*volFact*wx3+0.8164965809277261*f[5]*dv3*volFact; 
  out[7] += 2.828427124746191*f[1]*volFact*wx3; 
  out[8] += 2.828427124746191*f[2]*volFact*wx3; 
} 
void MomentCalc2x3vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  out[0] += 2.828427124746191*f[0]*volFact*wx1+0.8164965809277261*f[3]*dv1*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx1+0.8164965809277261*f[7]*dv1*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact*wx1+0.8164965809277261*f[8]*dv1*volFact; 
  out[3] += 2.828427124746191*f[6]*volFact*wx1; 
  out[4] += 2.828427124746191*f[16]*volFact*wx1; 
  out[5] += 2.828427124746191*f[17]*volFact*wx1; 
  out[6] += 2.828427124746191*f[0]*volFact*wx2+0.8164965809277261*f[4]*dv2*volFact; 
  out[7] += 2.828427124746191*f[1]*volFact*wx2+0.8164965809277261*f[9]*dv2*volFact; 
  out[8] += 2.828427124746191*f[2]*volFact*wx2+0.8164965809277261*f[10]*dv2*volFact; 
  out[9] += 2.828427124746191*f[6]*volFact*wx2; 
  out[10] += 2.828427124746191*f[16]*volFact*wx2; 
  out[11] += 2.828427124746191*f[17]*volFact*wx2; 
  out[12] += 2.828427124746191*f[0]*volFact*wx3+0.8164965809277261*f[5]*dv3*volFact; 
  out[13] += 2.828427124746191*f[1]*volFact*wx3+0.8164965809277261*f[12]*dv3*volFact; 
  out[14] += 2.828427124746191*f[2]*volFact*wx3+0.8164965809277261*f[13]*dv3*volFact; 
  out[15] += 2.828427124746191*f[6]*volFact*wx3; 
  out[16] += 2.828427124746191*f[16]*volFact*wx3; 
  out[17] += 2.828427124746191*f[17]*volFact*wx3; 
} 
void MomentCalc2x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += 2.828427124746191*f[0]*volFact*wx1_sq+1.632993161855453*f[3]*dv1*volFact*wx1+0.2357022603955158*f[0]*dv1_sq*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx1_sq+0.2357022603955158*f[1]*dv1_sq*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact*wx1_sq+0.2357022603955158*f[2]*dv1_sq*volFact; 
  out[3] += 2.828427124746191*f[0]*volFact*wx1*wx2+0.8164965809277261*f[3]*dv1*volFact*wx2+0.8164965809277261*f[4]*dv2*volFact*wx1; 
  out[4] += 2.828427124746191*f[1]*volFact*wx1*wx2; 
  out[5] += 2.828427124746191*f[2]*volFact*wx1*wx2; 
  out[6] += 2.828427124746191*f[0]*volFact*wx1*wx3+0.8164965809277261*f[3]*dv1*volFact*wx3+0.8164965809277261*f[5]*dv3*volFact*wx1; 
  out[7] += 2.828427124746191*f[1]*volFact*wx1*wx3; 
  out[8] += 2.828427124746191*f[2]*volFact*wx1*wx3; 
  out[9] += 2.828427124746191*f[0]*volFact*wx2_sq+1.632993161855453*f[4]*dv2*volFact*wx2+0.2357022603955158*f[0]*dv2_sq*volFact; 
  out[10] += 2.828427124746191*f[1]*volFact*wx2_sq+0.2357022603955158*f[1]*dv2_sq*volFact; 
  out[11] += 2.828427124746191*f[2]*volFact*wx2_sq+0.2357022603955158*f[2]*dv2_sq*volFact; 
  out[12] += 2.828427124746191*f[0]*volFact*wx2*wx3+0.8164965809277261*f[4]*dv2*volFact*wx3+0.8164965809277261*f[5]*dv3*volFact*wx2; 
  out[13] += 2.828427124746191*f[1]*volFact*wx2*wx3; 
  out[14] += 2.828427124746191*f[2]*volFact*wx2*wx3; 
  out[15] += 2.828427124746191*f[0]*volFact*wx3_sq+1.632993161855453*f[5]*dv3*volFact*wx3+0.2357022603955158*f[0]*dv3_sq*volFact; 
  out[16] += 2.828427124746191*f[1]*volFact*wx3_sq+0.2357022603955158*f[1]*dv3_sq*volFact; 
  out[17] += 2.828427124746191*f[2]*volFact*wx3_sq+0.2357022603955158*f[2]*dv3_sq*volFact; 
} 
void MomentCalc2x3vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += 2.828427124746191*f[0]*volFact*wx1_sq+1.632993161855453*f[3]*dv1*volFact*wx1+0.210818510677892*f[18]*dv1_sq*volFact+0.2357022603955158*f[0]*dv1_sq*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx1_sq+1.632993161855453*f[7]*dv1*volFact*wx1+0.2357022603955158*f[1]*dv1_sq*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact*wx1_sq+1.632993161855453*f[8]*dv1*volFact*wx1+0.2357022603955158*f[2]*dv1_sq*volFact; 
  out[3] += 2.828427124746191*f[6]*volFact*wx1_sq+0.2357022603955158*f[6]*dv1_sq*volFact; 
  out[4] += 2.828427124746191*f[16]*volFact*wx1_sq+0.2357022603955158*f[16]*dv1_sq*volFact; 
  out[5] += 2.828427124746191*f[17]*volFact*wx1_sq+0.2357022603955158*f[17]*dv1_sq*volFact; 
  out[6] += 2.828427124746191*f[0]*volFact*wx1*wx2+0.8164965809277261*f[3]*dv1*volFact*wx2+0.8164965809277261*f[4]*dv2*volFact*wx1+0.2357022603955158*f[11]*dv1*dv2*volFact; 
  out[7] += 2.828427124746191*f[1]*volFact*wx1*wx2+0.8164965809277261*f[7]*dv1*volFact*wx2+0.8164965809277261*f[9]*dv2*volFact*wx1; 
  out[8] += 2.828427124746191*f[2]*volFact*wx1*wx2+0.8164965809277261*f[8]*dv1*volFact*wx2+0.8164965809277261*f[10]*dv2*volFact*wx1; 
  out[9] += 2.828427124746191*f[6]*volFact*wx1*wx2; 
  out[10] += 2.828427124746191*f[16]*volFact*wx1*wx2; 
  out[11] += 2.828427124746191*f[17]*volFact*wx1*wx2; 
  out[12] += 2.828427124746191*f[0]*volFact*wx1*wx3+0.8164965809277261*f[3]*dv1*volFact*wx3+0.8164965809277261*f[5]*dv3*volFact*wx1+0.2357022603955158*f[14]*dv1*dv3*volFact; 
  out[13] += 2.828427124746191*f[1]*volFact*wx1*wx3+0.8164965809277261*f[7]*dv1*volFact*wx3+0.8164965809277261*f[12]*dv3*volFact*wx1; 
  out[14] += 2.828427124746191*f[2]*volFact*wx1*wx3+0.8164965809277261*f[8]*dv1*volFact*wx3+0.8164965809277261*f[13]*dv3*volFact*wx1; 
  out[15] += 2.828427124746191*f[6]*volFact*wx1*wx3; 
  out[16] += 2.828427124746191*f[16]*volFact*wx1*wx3; 
  out[17] += 2.828427124746191*f[17]*volFact*wx1*wx3; 
  out[18] += 2.828427124746191*f[0]*volFact*wx2_sq+1.632993161855453*f[4]*dv2*volFact*wx2+0.210818510677892*f[19]*dv2_sq*volFact+0.2357022603955158*f[0]*dv2_sq*volFact; 
  out[19] += 2.828427124746191*f[1]*volFact*wx2_sq+1.632993161855453*f[9]*dv2*volFact*wx2+0.2357022603955158*f[1]*dv2_sq*volFact; 
  out[20] += 2.828427124746191*f[2]*volFact*wx2_sq+1.632993161855453*f[10]*dv2*volFact*wx2+0.2357022603955158*f[2]*dv2_sq*volFact; 
  out[21] += 2.828427124746191*f[6]*volFact*wx2_sq+0.2357022603955158*f[6]*dv2_sq*volFact; 
  out[22] += 2.828427124746191*f[16]*volFact*wx2_sq+0.2357022603955158*f[16]*dv2_sq*volFact; 
  out[23] += 2.828427124746191*f[17]*volFact*wx2_sq+0.2357022603955158*f[17]*dv2_sq*volFact; 
  out[24] += 2.828427124746191*f[0]*volFact*wx2*wx3+0.8164965809277261*f[4]*dv2*volFact*wx3+0.8164965809277261*f[5]*dv3*volFact*wx2+0.2357022603955158*f[15]*dv2*dv3*volFact; 
  out[25] += 2.828427124746191*f[1]*volFact*wx2*wx3+0.8164965809277261*f[9]*dv2*volFact*wx3+0.8164965809277261*f[12]*dv3*volFact*wx2; 
  out[26] += 2.828427124746191*f[2]*volFact*wx2*wx3+0.8164965809277261*f[10]*dv2*volFact*wx3+0.8164965809277261*f[13]*dv3*volFact*wx2; 
  out[27] += 2.828427124746191*f[6]*volFact*wx2*wx3; 
  out[28] += 2.828427124746191*f[16]*volFact*wx2*wx3; 
  out[29] += 2.828427124746191*f[17]*volFact*wx2*wx3; 
  out[30] += 2.828427124746191*f[0]*volFact*wx3_sq+1.632993161855453*f[5]*dv3*volFact*wx3+0.210818510677892*f[20]*dv3_sq*volFact+0.2357022603955158*f[0]*dv3_sq*volFact; 
  out[31] += 2.828427124746191*f[1]*volFact*wx3_sq+1.632993161855453*f[12]*dv3*volFact*wx3+0.2357022603955158*f[1]*dv3_sq*volFact; 
  out[32] += 2.828427124746191*f[2]*volFact*wx3_sq+1.632993161855453*f[13]*dv3*volFact*wx3+0.2357022603955158*f[2]*dv3_sq*volFact; 
  out[33] += 2.828427124746191*f[6]*volFact*wx3_sq+0.2357022603955158*f[6]*dv3_sq*volFact; 
  out[34] += 2.828427124746191*f[16]*volFact*wx3_sq+0.2357022603955158*f[16]*dv3_sq*volFact; 
  out[35] += 2.828427124746191*f[17]*volFact*wx3_sq+0.2357022603955158*f[17]*dv3_sq*volFact; 
} 
void MomentCalc2x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += 2.828427124746191*f[0]*volFact*wx3_sq+1.632993161855453*f[5]*dv3*volFact*wx3+2.828427124746191*f[0]*volFact*wx2_sq+1.632993161855453*f[4]*dv2*volFact*wx2+2.828427124746191*f[0]*volFact*wx1_sq+1.632993161855453*f[3]*dv1*volFact*wx1+0.2357022603955158*f[0]*dv3_sq*volFact+0.2357022603955158*f[0]*dv2_sq*volFact+0.2357022603955158*f[0]*dv1_sq*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx3_sq+2.828427124746191*f[1]*volFact*wx2_sq+2.828427124746191*f[1]*volFact*wx1_sq+0.2357022603955158*f[1]*dv3_sq*volFact+0.2357022603955158*f[1]*dv2_sq*volFact+0.2357022603955158*f[1]*dv1_sq*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact*wx3_sq+2.828427124746191*f[2]*volFact*wx2_sq+2.828427124746191*f[2]*volFact*wx1_sq+0.2357022603955158*f[2]*dv3_sq*volFact+0.2357022603955158*f[2]*dv2_sq*volFact+0.2357022603955158*f[2]*dv1_sq*volFact; 
} 
void MomentCalc2x3vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[2]*dxv[3]*dxv[4]/8; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += 2.828427124746191*f[0]*volFact*wx3_sq+1.632993161855453*f[5]*dv3*volFact*wx3+2.828427124746191*f[0]*volFact*wx2_sq+1.632993161855453*f[4]*dv2*volFact*wx2+2.828427124746191*f[0]*volFact*wx1_sq+1.632993161855453*f[3]*dv1*volFact*wx1+0.210818510677892*f[20]*dv3_sq*volFact+0.2357022603955158*f[0]*dv3_sq*volFact+0.210818510677892*f[19]*dv2_sq*volFact+0.2357022603955158*f[0]*dv2_sq*volFact+0.210818510677892*f[18]*dv1_sq*volFact+0.2357022603955158*f[0]*dv1_sq*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx3_sq+1.632993161855453*f[12]*dv3*volFact*wx3+2.828427124746191*f[1]*volFact*wx2_sq+1.632993161855453*f[9]*dv2*volFact*wx2+2.828427124746191*f[1]*volFact*wx1_sq+1.632993161855453*f[7]*dv1*volFact*wx1+0.2357022603955158*f[1]*dv3_sq*volFact+0.2357022603955158*f[1]*dv2_sq*volFact+0.2357022603955158*f[1]*dv1_sq*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact*wx3_sq+1.632993161855453*f[13]*dv3*volFact*wx3+2.828427124746191*f[2]*volFact*wx2_sq+1.632993161855453*f[10]*dv2*volFact*wx2+2.828427124746191*f[2]*volFact*wx1_sq+1.632993161855453*f[8]*dv1*volFact*wx1+0.2357022603955158*f[2]*dv3_sq*volFact+0.2357022603955158*f[2]*dv2_sq*volFact+0.2357022603955158*f[2]*dv1_sq*volFact; 
  out[3] += 2.828427124746191*f[6]*volFact*wx3_sq+2.828427124746191*f[6]*volFact*wx2_sq+2.828427124746191*f[6]*volFact*wx1_sq+0.2357022603955158*f[6]*dv3_sq*volFact+0.2357022603955158*f[6]*dv2_sq*volFact+0.2357022603955158*f[6]*dv1_sq*volFact; 
  out[4] += 2.828427124746191*f[16]*volFact*wx3_sq+2.828427124746191*f[16]*volFact*wx2_sq+2.828427124746191*f[16]*volFact*wx1_sq+0.2357022603955158*f[16]*dv3_sq*volFact+0.2357022603955158*f[16]*dv2_sq*volFact+0.2357022603955158*f[16]*dv1_sq*volFact; 
  out[5] += 2.828427124746191*f[17]*volFact*wx3_sq+2.828427124746191*f[17]*volFact*wx2_sq+2.828427124746191*f[17]*volFact*wx1_sq+0.2357022603955158*f[17]*dv3_sq*volFact+0.2357022603955158*f[17]*dv2_sq*volFact+0.2357022603955158*f[17]*dv1_sq*volFact; 
} 
