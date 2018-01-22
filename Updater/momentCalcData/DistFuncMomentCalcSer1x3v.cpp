#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc1x3vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
} 
void MomentCalc1x3vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[11]*volFact; 
} 
void MomentCalc1x3vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  out[0] += 2.828427124746191*f[0]*volFact*wx1+0.8164965809277261*f[2]*dv1*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx1+0.8164965809277261*f[5]*dv1*volFact; 
  out[2] += 2.828427124746191*f[0]*volFact*wx2+0.8164965809277261*f[3]*dv2*volFact; 
  out[3] += 2.828427124746191*f[1]*volFact*wx2+0.8164965809277261*f[6]*dv2*volFact; 
  out[4] += 2.828427124746191*f[0]*volFact*wx3+0.8164965809277261*f[4]*dv3*volFact; 
  out[5] += 2.828427124746191*f[1]*volFact*wx3+0.8164965809277261*f[8]*dv3*volFact; 
} 
void MomentCalc1x3vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  out[0] += 2.828427124746191*f[0]*volFact*wx1+0.8164965809277261*f[2]*dv1*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx1+0.8164965809277261*f[5]*dv1*volFact; 
  out[2] += 2.828427124746191*f[11]*volFact*wx1+0.816496580927726*f[19]*dv1*volFact; 
  out[3] += 2.828427124746191*f[0]*volFact*wx2+0.8164965809277261*f[3]*dv2*volFact; 
  out[4] += 2.828427124746191*f[1]*volFact*wx2+0.8164965809277261*f[6]*dv2*volFact; 
  out[5] += 2.828427124746191*f[11]*volFact*wx2+0.816496580927726*f[21]*dv2*volFact; 
  out[6] += 2.828427124746191*f[0]*volFact*wx3+0.8164965809277261*f[4]*dv3*volFact; 
  out[7] += 2.828427124746191*f[1]*volFact*wx3+0.8164965809277261*f[8]*dv3*volFact; 
  out[8] += 2.828427124746191*f[11]*volFact*wx3+0.816496580927726*f[25]*dv3*volFact; 
} 
void MomentCalc1x3vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += 2.828427124746191*f[0]*volFact*wx1_sq+1.632993161855453*f[2]*dv1*volFact*wx1+0.2357022603955158*f[0]*dv1_sq*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx1_sq+1.632993161855453*f[5]*dv1*volFact*wx1+0.2357022603955158*f[1]*dv1_sq*volFact; 
  out[2] += 2.828427124746191*f[0]*volFact*wx1*wx2+0.8164965809277261*f[2]*dv1*volFact*wx2+0.8164965809277261*f[3]*dv2*volFact*wx1+0.2357022603955158*f[7]*dv1*dv2*volFact; 
  out[3] += 2.828427124746191*f[1]*volFact*wx1*wx2+0.8164965809277261*f[5]*dv1*volFact*wx2+0.8164965809277261*f[6]*dv2*volFact*wx1+0.2357022603955158*f[11]*dv1*dv2*volFact; 
  out[4] += 2.828427124746191*f[0]*volFact*wx1*wx3+0.8164965809277261*f[2]*dv1*volFact*wx3+0.8164965809277261*f[4]*dv3*volFact*wx1+0.2357022603955158*f[9]*dv1*dv3*volFact; 
  out[5] += 2.828427124746191*f[1]*volFact*wx1*wx3+0.8164965809277261*f[5]*dv1*volFact*wx3+0.8164965809277261*f[8]*dv3*volFact*wx1+0.2357022603955158*f[12]*dv1*dv3*volFact; 
  out[6] += 2.828427124746191*f[0]*volFact*wx2_sq+1.632993161855453*f[3]*dv2*volFact*wx2+0.2357022603955158*f[0]*dv2_sq*volFact; 
  out[7] += 2.828427124746191*f[1]*volFact*wx2_sq+1.632993161855453*f[6]*dv2*volFact*wx2+0.2357022603955158*f[1]*dv2_sq*volFact; 
  out[8] += 2.828427124746191*f[0]*volFact*wx2*wx3+0.8164965809277261*f[3]*dv2*volFact*wx3+0.8164965809277261*f[4]*dv3*volFact*wx2+0.2357022603955158*f[10]*dv2*dv3*volFact; 
  out[9] += 2.828427124746191*f[1]*volFact*wx2*wx3+0.8164965809277261*f[6]*dv2*volFact*wx3+0.8164965809277261*f[8]*dv3*volFact*wx2+0.2357022603955158*f[13]*dv2*dv3*volFact; 
  out[10] += 2.828427124746191*f[0]*volFact*wx3_sq+1.632993161855453*f[4]*dv3*volFact*wx3+0.2357022603955158*f[0]*dv3_sq*volFact; 
  out[11] += 2.828427124746191*f[1]*volFact*wx3_sq+1.632993161855453*f[8]*dv3*volFact*wx3+0.2357022603955158*f[1]*dv3_sq*volFact; 
} 
void MomentCalc1x3vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += 2.828427124746191*f[0]*volFact*wx1_sq+1.632993161855453*f[2]*dv1*volFact*wx1+0.210818510677892*f[12]*dv1_sq*volFact+0.2357022603955158*f[0]*dv1_sq*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx1_sq+1.632993161855453*f[5]*dv1*volFact*wx1+0.2108185106778921*f[20]*dv1_sq*volFact+0.2357022603955158*f[1]*dv1_sq*volFact; 
  out[2] += 2.828427124746191*f[11]*volFact*wx1_sq+1.632993161855453*f[19]*dv1*volFact*wx1+0.2357022603955158*f[11]*dv1_sq*volFact; 
  out[3] += 2.828427124746191*f[0]*volFact*wx1*wx2+0.8164965809277261*f[2]*dv1*volFact*wx2+0.8164965809277261*f[3]*dv2*volFact*wx1+0.2357022603955158*f[7]*dv1*dv2*volFact; 
  out[4] += 2.828427124746191*f[1]*volFact*wx1*wx2+0.8164965809277261*f[5]*dv1*volFact*wx2+0.8164965809277261*f[6]*dv2*volFact*wx1+0.2357022603955158*f[15]*dv1*dv2*volFact; 
  out[5] += 2.828427124746191*f[11]*volFact*wx1*wx2+0.816496580927726*f[19]*dv1*volFact*wx2+0.816496580927726*f[21]*dv2*volFact*wx1+0.2357022603955158*f[32]*dv1*dv2*volFact; 
  out[6] += 2.828427124746191*f[0]*volFact*wx1*wx3+0.8164965809277261*f[2]*dv1*volFact*wx3+0.8164965809277261*f[4]*dv3*volFact*wx1+0.2357022603955158*f[9]*dv1*dv3*volFact; 
  out[7] += 2.828427124746191*f[1]*volFact*wx1*wx3+0.8164965809277261*f[5]*dv1*volFact*wx3+0.8164965809277261*f[8]*dv3*volFact*wx1+0.2357022603955158*f[16]*dv1*dv3*volFact; 
  out[8] += 2.828427124746191*f[11]*volFact*wx1*wx3+0.816496580927726*f[19]*dv1*volFact*wx3+0.816496580927726*f[25]*dv3*volFact*wx1+0.2357022603955158*f[35]*dv1*dv3*volFact; 
  out[9] += 2.828427124746191*f[0]*volFact*wx2_sq+1.632993161855453*f[3]*dv2*volFact*wx2+0.210818510677892*f[13]*dv2_sq*volFact+0.2357022603955158*f[0]*dv2_sq*volFact; 
  out[10] += 2.828427124746191*f[1]*volFact*wx2_sq+1.632993161855453*f[6]*dv2*volFact*wx2+0.2108185106778921*f[23]*dv2_sq*volFact+0.2357022603955158*f[1]*dv2_sq*volFact; 
  out[11] += 2.828427124746191*f[11]*volFact*wx2_sq+1.632993161855453*f[21]*dv2*volFact*wx2+0.2357022603955158*f[11]*dv2_sq*volFact; 
  out[12] += 2.828427124746191*f[0]*volFact*wx2*wx3+0.8164965809277261*f[3]*dv2*volFact*wx3+0.8164965809277261*f[4]*dv3*volFact*wx2+0.2357022603955158*f[10]*dv2*dv3*volFact; 
  out[13] += 2.828427124746191*f[1]*volFact*wx2*wx3+0.8164965809277261*f[6]*dv2*volFact*wx3+0.8164965809277261*f[8]*dv3*volFact*wx2+0.2357022603955158*f[17]*dv2*dv3*volFact; 
  out[14] += 2.828427124746191*f[11]*volFact*wx2*wx3+0.816496580927726*f[21]*dv2*volFact*wx3+0.816496580927726*f[25]*dv3*volFact*wx2+0.2357022603955158*f[37]*dv2*dv3*volFact; 
  out[15] += 2.828427124746191*f[0]*volFact*wx3_sq+1.632993161855453*f[4]*dv3*volFact*wx3+0.210818510677892*f[14]*dv3_sq*volFact+0.2357022603955158*f[0]*dv3_sq*volFact; 
  out[16] += 2.828427124746191*f[1]*volFact*wx3_sq+1.632993161855453*f[8]*dv3*volFact*wx3+0.2108185106778921*f[28]*dv3_sq*volFact+0.2357022603955158*f[1]*dv3_sq*volFact; 
  out[17] += 2.828427124746191*f[11]*volFact*wx3_sq+1.632993161855453*f[25]*dv3*volFact*wx3+0.2357022603955158*f[11]*dv3_sq*volFact; 
} 
void MomentCalc1x3vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += 2.828427124746191*f[0]*volFact*wx3_sq+1.632993161855453*f[4]*dv3*volFact*wx3+2.828427124746191*f[0]*volFact*wx2_sq+1.632993161855453*f[3]*dv2*volFact*wx2+2.828427124746191*f[0]*volFact*wx1_sq+1.632993161855453*f[2]*dv1*volFact*wx1+0.2357022603955158*f[0]*dv3_sq*volFact+0.2357022603955158*f[0]*dv2_sq*volFact+0.2357022603955158*f[0]*dv1_sq*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx3_sq+1.632993161855453*f[8]*dv3*volFact*wx3+2.828427124746191*f[1]*volFact*wx2_sq+1.632993161855453*f[6]*dv2*volFact*wx2+2.828427124746191*f[1]*volFact*wx1_sq+1.632993161855453*f[5]*dv1*volFact*wx1+0.2357022603955158*f[1]*dv3_sq*volFact+0.2357022603955158*f[1]*dv2_sq*volFact+0.2357022603955158*f[1]*dv1_sq*volFact; 
} 
void MomentCalc1x3vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += 2.828427124746191*f[0]*volFact*wx3_sq+1.632993161855453*f[4]*dv3*volFact*wx3+2.828427124746191*f[0]*volFact*wx2_sq+1.632993161855453*f[3]*dv2*volFact*wx2+2.828427124746191*f[0]*volFact*wx1_sq+1.632993161855453*f[2]*dv1*volFact*wx1+0.210818510677892*f[14]*dv3_sq*volFact+0.2357022603955158*f[0]*dv3_sq*volFact+0.210818510677892*f[13]*dv2_sq*volFact+0.2357022603955158*f[0]*dv2_sq*volFact+0.210818510677892*f[12]*dv1_sq*volFact+0.2357022603955158*f[0]*dv1_sq*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact*wx3_sq+1.632993161855453*f[8]*dv3*volFact*wx3+2.828427124746191*f[1]*volFact*wx2_sq+1.632993161855453*f[6]*dv2*volFact*wx2+2.828427124746191*f[1]*volFact*wx1_sq+1.632993161855453*f[5]*dv1*volFact*wx1+0.2108185106778921*f[28]*dv3_sq*volFact+0.2357022603955158*f[1]*dv3_sq*volFact+0.2108185106778921*f[23]*dv2_sq*volFact+0.2357022603955158*f[1]*dv2_sq*volFact+0.2108185106778921*f[20]*dv1_sq*volFact+0.2357022603955158*f[1]*dv1_sq*volFact; 
  out[2] += 2.828427124746191*f[11]*volFact*wx3_sq+1.632993161855453*f[25]*dv3*volFact*wx3+2.828427124746191*f[11]*volFact*wx2_sq+1.632993161855453*f[21]*dv2*volFact*wx2+2.828427124746191*f[11]*volFact*wx1_sq+1.632993161855453*f[19]*dv1*volFact*wx1+0.2357022603955158*f[11]*dv3_sq*volFact+0.2357022603955158*f[11]*dv2_sq*volFact+0.2357022603955158*f[11]*dv1_sq*volFact; 
} 
