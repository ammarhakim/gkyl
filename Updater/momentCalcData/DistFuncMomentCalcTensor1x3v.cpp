#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc1x3vTensor_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
} 
void MomentCalc1x3vTensor_M0_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[11]*volFact; 
} 
void MomentCalc1x3vTensor_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  out[0] += (2.828427124746191*f[0]*wx1+0.8164965809277261*f[2]*dv1)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1+0.8164965809277261*f[5]*dv1)*volFact; 
  out[2] += (2.828427124746191*f[0]*wx2+0.8164965809277261*f[3]*dv2)*volFact; 
  out[3] += (2.828427124746191*f[1]*wx2+0.8164965809277261*f[6]*dv2)*volFact; 
  out[4] += (2.828427124746191*f[0]*wx3+0.8164965809277261*f[4]*dv3)*volFact; 
  out[5] += (2.828427124746191*f[1]*wx3+0.8164965809277261*f[8]*dv3)*volFact; 
} 
void MomentCalc1x3vTensor_M1i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  out[0] += (2.828427124746191*f[0]*wx1+0.8164965809277261*f[2]*dv1)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1+0.8164965809277261*f[5]*dv1)*volFact; 
  out[2] += (2.828427124746191*f[11]*wx1+0.816496580927726*f[19]*dv1)*volFact; 
  out[3] += (2.828427124746191*f[0]*wx2+0.8164965809277261*f[3]*dv2)*volFact; 
  out[4] += (2.828427124746191*f[1]*wx2+0.8164965809277261*f[6]*dv2)*volFact; 
  out[5] += (2.828427124746191*f[11]*wx2+0.816496580927726*f[21]*dv2)*volFact; 
  out[6] += (2.828427124746191*f[0]*wx3+0.8164965809277261*f[4]*dv3)*volFact; 
  out[7] += (2.828427124746191*f[1]*wx3+0.8164965809277261*f[8]*dv3)*volFact; 
  out[8] += (2.828427124746191*f[11]*wx3+0.816496580927726*f[25]*dv3)*volFact; 
} 
void MomentCalc1x3vTensor_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[2]*dv1*wx1+0.2357022603955158*f[0]*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[5]*dv1*wx1+0.2357022603955158*f[1]*dv1_sq)*volFact; 
  out[2] += (2.828427124746191*f[0]*wx1*wx2+0.8164965809277261*f[2]*dv1*wx2+0.8164965809277261*f[3]*dv2*wx1+0.2357022603955158*f[7]*dv1*dv2)*volFact; 
  out[3] += (2.828427124746191*f[1]*wx1*wx2+0.8164965809277261*f[5]*dv1*wx2+0.8164965809277261*f[6]*dv2*wx1+0.2357022603955158*f[11]*dv1*dv2)*volFact; 
  out[4] += (2.828427124746191*f[0]*wx1*wx3+0.8164965809277261*f[2]*dv1*wx3+0.8164965809277261*f[4]*dv3*wx1+0.2357022603955158*f[9]*dv1*dv3)*volFact; 
  out[5] += (2.828427124746191*f[1]*wx1*wx3+0.8164965809277261*f[5]*dv1*wx3+0.8164965809277261*f[8]*dv3*wx1+0.2357022603955158*f[12]*dv1*dv3)*volFact; 
  out[6] += (2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[3]*dv2*wx2+0.2357022603955158*f[0]*dv2_sq)*volFact; 
  out[7] += (2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[6]*dv2*wx2+0.2357022603955158*f[1]*dv2_sq)*volFact; 
  out[8] += (2.828427124746191*f[0]*wx2*wx3+0.8164965809277261*f[3]*dv2*wx3+0.8164965809277261*f[4]*dv3*wx2+0.2357022603955158*f[10]*dv2*dv3)*volFact; 
  out[9] += (2.828427124746191*f[1]*wx2*wx3+0.8164965809277261*f[6]*dv2*wx3+0.8164965809277261*f[8]*dv3*wx2+0.2357022603955158*f[13]*dv2*dv3)*volFact; 
  out[10] += (2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[4]*dv3*wx3+0.2357022603955158*f[0]*dv3_sq)*volFact; 
  out[11] += (2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[8]*dv3*wx3+0.2357022603955158*f[1]*dv3_sq)*volFact; 
} 
void MomentCalc1x3vTensor_M2ij_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[2]*dv1*wx1+0.210818510677892*f[12]*dv1_sq+0.2357022603955158*f[0]*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[5]*dv1*wx1+0.2108185106778921*f[20]*dv1_sq+0.2357022603955158*f[1]*dv1_sq)*volFact; 
  out[2] += (2.828427124746191*f[11]*wx1_sq+1.632993161855453*f[19]*dv1*wx1+0.210818510677892*f[44]*dv1_sq+0.2357022603955158*f[11]*dv1_sq)*volFact; 
  out[3] += (2.828427124746191*f[0]*wx1*wx2+0.8164965809277261*f[2]*dv1*wx2+0.8164965809277261*f[3]*dv2*wx1+0.2357022603955158*f[7]*dv1*dv2)*volFact; 
  out[4] += (2.828427124746191*f[1]*wx1*wx2+0.8164965809277261*f[5]*dv1*wx2+0.8164965809277261*f[6]*dv2*wx1+0.2357022603955158*f[15]*dv1*dv2)*volFact; 
  out[5] += (2.828427124746191*f[11]*wx1*wx2+0.816496580927726*f[19]*dv1*wx2+0.816496580927726*f[21]*dv2*wx1+0.2357022603955158*f[32]*dv1*dv2)*volFact; 
  out[6] += (2.828427124746191*f[0]*wx1*wx3+0.8164965809277261*f[2]*dv1*wx3+0.8164965809277261*f[4]*dv3*wx1+0.2357022603955158*f[9]*dv1*dv3)*volFact; 
  out[7] += (2.828427124746191*f[1]*wx1*wx3+0.8164965809277261*f[5]*dv1*wx3+0.8164965809277261*f[8]*dv3*wx1+0.2357022603955158*f[16]*dv1*dv3)*volFact; 
  out[8] += (2.828427124746191*f[11]*wx1*wx3+0.816496580927726*f[19]*dv1*wx3+0.816496580927726*f[25]*dv3*wx1+0.2357022603955158*f[35]*dv1*dv3)*volFact; 
  out[9] += (2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[3]*dv2*wx2+0.210818510677892*f[13]*dv2_sq+0.2357022603955158*f[0]*dv2_sq)*volFact; 
  out[10] += (2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[6]*dv2*wx2+0.2108185106778921*f[23]*dv2_sq+0.2357022603955158*f[1]*dv2_sq)*volFact; 
  out[11] += (2.828427124746191*f[11]*wx2_sq+1.632993161855453*f[21]*dv2*wx2+0.210818510677892*f[45]*dv2_sq+0.2357022603955158*f[11]*dv2_sq)*volFact; 
  out[12] += (2.828427124746191*f[0]*wx2*wx3+0.8164965809277261*f[3]*dv2*wx3+0.8164965809277261*f[4]*dv3*wx2+0.2357022603955158*f[10]*dv2*dv3)*volFact; 
  out[13] += (2.828427124746191*f[1]*wx2*wx3+0.8164965809277261*f[6]*dv2*wx3+0.8164965809277261*f[8]*dv3*wx2+0.2357022603955158*f[17]*dv2*dv3)*volFact; 
  out[14] += (2.828427124746191*f[11]*wx2*wx3+0.816496580927726*f[21]*dv2*wx3+0.816496580927726*f[25]*dv3*wx2+0.2357022603955158*f[37]*dv2*dv3)*volFact; 
  out[15] += (2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[4]*dv3*wx3+0.210818510677892*f[14]*dv3_sq+0.2357022603955158*f[0]*dv3_sq)*volFact; 
  out[16] += (2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[8]*dv3*wx3+0.2108185106778921*f[28]*dv3_sq+0.2357022603955158*f[1]*dv3_sq)*volFact; 
  out[17] += (2.828427124746191*f[11]*wx3_sq+1.632993161855453*f[25]*dv3*wx3+0.210818510677892*f[47]*dv3_sq+0.2357022603955158*f[11]*dv3_sq)*volFact; 
} 
void MomentCalc1x3vTensor_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[4]*dv3*wx3+2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[3]*dv2*wx2+2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[2]*dv1*wx1+0.2357022603955158*f[0]*dv3_sq+0.2357022603955158*f[0]*dv2_sq+0.2357022603955158*f[0]*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[8]*dv3*wx3+2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[6]*dv2*wx2+2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[5]*dv1*wx1+0.2357022603955158*f[1]*dv3_sq+0.2357022603955158*f[1]*dv2_sq+0.2357022603955158*f[1]*dv1_sq)*volFact; 
} 
void MomentCalc1x3vTensor_M2_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[4]*dv3*wx3+2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[3]*dv2*wx2+2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[2]*dv1*wx1+0.210818510677892*f[14]*dv3_sq+0.2357022603955158*f[0]*dv3_sq+0.210818510677892*f[13]*dv2_sq+0.2357022603955158*f[0]*dv2_sq+0.210818510677892*f[12]*dv1_sq+0.2357022603955158*f[0]*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx3_sq+1.632993161855453*f[8]*dv3*wx3+2.828427124746191*f[1]*wx2_sq+1.632993161855453*f[6]*dv2*wx2+2.828427124746191*f[1]*wx1_sq+1.632993161855453*f[5]*dv1*wx1+0.2108185106778921*f[28]*dv3_sq+0.2357022603955158*f[1]*dv3_sq+0.2108185106778921*f[23]*dv2_sq+0.2357022603955158*f[1]*dv2_sq+0.2108185106778921*f[20]*dv1_sq+0.2357022603955158*f[1]*dv1_sq)*volFact; 
  out[2] += (2.828427124746191*f[11]*wx3_sq+1.632993161855453*f[25]*dv3*wx3+2.828427124746191*f[11]*wx2_sq+1.632993161855453*f[21]*dv2*wx2+2.828427124746191*f[11]*wx1_sq+1.632993161855453*f[19]*dv1*wx1+0.210818510677892*f[47]*dv3_sq+0.2357022603955158*f[11]*dv3_sq+0.210818510677892*f[45]*dv2_sq+0.2357022603955158*f[11]*dv2_sq+0.210818510677892*f[44]*dv1_sq+0.2357022603955158*f[11]*dv1_sq)*volFact; 
} 
void MomentCalc1x3vTensor_M3i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double wx3_cu = wx3*wx3*wx3, dv3_cu = dv3*dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx1*wx3_sq+0.8164965809277261*f[2]*dv1*wx3_sq+1.632993161855453*f[4]*dv3*wx1*wx3+0.4714045207910317*f[9]*dv1*dv3*wx3+2.828427124746191*f[0]*wx1*wx2_sq+0.8164965809277261*f[2]*dv1*wx2_sq+1.632993161855453*f[3]*dv2*wx1*wx2+0.4714045207910317*f[7]*dv1*dv2*wx2+2.828427124746191*f[0]*wx1*wx1_sq+2.449489742783178*f[2]*dv1*wx1_sq+0.2357022603955158*f[0]*dv3_sq*wx1+0.2357022603955158*f[0]*dv2_sq*wx1+0.7071067811865475*f[0]*dv1_sq*wx1+0.06804138174397717*f[2]*dv1*dv3_sq+0.06804138174397717*f[2]*dv1*dv2_sq+0.1224744871391589*f[2]*dv1*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1*wx3_sq+0.8164965809277261*f[5]*dv1*wx3_sq+1.632993161855453*f[8]*dv3*wx1*wx3+0.4714045207910317*f[12]*dv1*dv3*wx3+2.828427124746191*f[1]*wx1*wx2_sq+0.8164965809277261*f[5]*dv1*wx2_sq+1.632993161855453*f[6]*dv2*wx1*wx2+0.4714045207910317*f[11]*dv1*dv2*wx2+2.828427124746191*f[1]*wx1*wx1_sq+2.449489742783178*f[5]*dv1*wx1_sq+0.2357022603955158*f[1]*dv3_sq*wx1+0.2357022603955158*f[1]*dv2_sq*wx1+0.7071067811865475*f[1]*dv1_sq*wx1+0.06804138174397717*f[5]*dv1*dv3_sq+0.06804138174397717*f[5]*dv1*dv2_sq+0.1224744871391589*f[5]*dv1*dv1_sq)*volFact; 
  out[2] += (2.828427124746191*f[0]*wx2*wx3_sq+0.8164965809277261*f[3]*dv2*wx3_sq+1.632993161855453*f[4]*dv3*wx2*wx3+0.4714045207910317*f[10]*dv2*dv3*wx3+2.828427124746191*f[0]*wx2*wx2_sq+2.449489742783178*f[3]*dv2*wx2_sq+2.828427124746191*f[0]*wx1_sq*wx2+1.632993161855453*f[2]*dv1*wx1*wx2+0.2357022603955158*f[0]*dv3_sq*wx2+0.7071067811865475*f[0]*dv2_sq*wx2+0.2357022603955158*f[0]*dv1_sq*wx2+0.8164965809277261*f[3]*dv2*wx1_sq+0.4714045207910317*f[7]*dv1*dv2*wx1+0.06804138174397717*f[3]*dv2*dv3_sq+0.1224744871391589*f[3]*dv2*dv2_sq+0.06804138174397717*f[3]*dv1_sq*dv2)*volFact; 
  out[3] += (2.828427124746191*f[1]*wx2*wx3_sq+0.8164965809277261*f[6]*dv2*wx3_sq+1.632993161855453*f[8]*dv3*wx2*wx3+0.4714045207910317*f[13]*dv2*dv3*wx3+2.828427124746191*f[1]*wx2*wx2_sq+2.449489742783178*f[6]*dv2*wx2_sq+2.828427124746191*f[1]*wx1_sq*wx2+1.632993161855453*f[5]*dv1*wx1*wx2+0.2357022603955158*f[1]*dv3_sq*wx2+0.7071067811865475*f[1]*dv2_sq*wx2+0.2357022603955158*f[1]*dv1_sq*wx2+0.8164965809277261*f[6]*dv2*wx1_sq+0.4714045207910317*f[11]*dv1*dv2*wx1+0.06804138174397717*f[6]*dv2*dv3_sq+0.1224744871391589*f[6]*dv2*dv2_sq+0.06804138174397717*f[6]*dv1_sq*dv2)*volFact; 
  out[4] += (2.828427124746191*f[0]*wx3*wx3_sq+2.449489742783178*f[4]*dv3*wx3_sq+2.828427124746191*f[0]*wx2_sq*wx3+1.632993161855453*f[3]*dv2*wx2*wx3+2.828427124746191*f[0]*wx1_sq*wx3+1.632993161855453*f[2]*dv1*wx1*wx3+0.7071067811865475*f[0]*dv3_sq*wx3+0.2357022603955158*f[0]*dv2_sq*wx3+0.2357022603955158*f[0]*dv1_sq*wx3+0.8164965809277261*f[4]*dv3*wx2_sq+0.4714045207910317*f[10]*dv2*dv3*wx2+0.8164965809277261*f[4]*dv3*wx1_sq+0.4714045207910317*f[9]*dv1*dv3*wx1+0.1224744871391589*f[4]*dv3*dv3_sq+0.06804138174397717*f[4]*dv2_sq*dv3+0.06804138174397717*f[4]*dv1_sq*dv3)*volFact; 
  out[5] += (2.828427124746191*f[1]*wx3*wx3_sq+2.449489742783178*f[8]*dv3*wx3_sq+2.828427124746191*f[1]*wx2_sq*wx3+1.632993161855453*f[6]*dv2*wx2*wx3+2.828427124746191*f[1]*wx1_sq*wx3+1.632993161855453*f[5]*dv1*wx1*wx3+0.7071067811865475*f[1]*dv3_sq*wx3+0.2357022603955158*f[1]*dv2_sq*wx3+0.2357022603955158*f[1]*dv1_sq*wx3+0.8164965809277261*f[8]*dv3*wx2_sq+0.4714045207910317*f[13]*dv2*dv3*wx2+0.8164965809277261*f[8]*dv3*wx1_sq+0.4714045207910317*f[12]*dv1*dv3*wx1+0.1224744871391589*f[8]*dv3*dv3_sq+0.06804138174397717*f[8]*dv2_sq*dv3+0.06804138174397717*f[8]*dv1_sq*dv3)*volFact; 
} 
void MomentCalc1x3vTensor_M3i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  const double wx3_cu = wx3*wx3*wx3, dv3_cu = dv3*dv3*dv3; 
  out[0] += (2.828427124746191*f[0]*wx1*wx3_sq+0.8164965809277261*f[2]*dv1*wx3_sq+1.632993161855453*f[4]*dv3*wx1*wx3+0.4714045207910317*f[9]*dv1*dv3*wx3+2.828427124746191*f[0]*wx1*wx2_sq+0.8164965809277261*f[2]*dv1*wx2_sq+1.632993161855453*f[3]*dv2*wx1*wx2+0.4714045207910317*f[7]*dv1*dv2*wx2+2.828427124746191*f[0]*wx1*wx1_sq+2.449489742783178*f[2]*dv1*wx1_sq+0.210818510677892*f[14]*dv3_sq*wx1+0.2357022603955158*f[0]*dv3_sq*wx1+0.210818510677892*f[13]*dv2_sq*wx1+0.2357022603955158*f[0]*dv2_sq*wx1+0.6324555320336759*f[12]*dv1_sq*wx1+0.7071067811865475*f[0]*dv1_sq*wx1+0.06085806194501844*f[29]*dv1*dv3_sq+0.06804138174397717*f[2]*dv1*dv3_sq+0.06085806194501844*f[24]*dv1*dv2_sq+0.06804138174397717*f[2]*dv1*dv2_sq+0.1224744871391589*f[2]*dv1*dv1_sq)*volFact; 
  out[1] += (2.828427124746191*f[1]*wx1*wx3_sq+0.8164965809277261*f[5]*dv1*wx3_sq+1.632993161855453*f[8]*dv3*wx1*wx3+0.4714045207910317*f[16]*dv1*dv3*wx3+2.828427124746191*f[1]*wx1*wx2_sq+0.8164965809277261*f[5]*dv1*wx2_sq+1.632993161855453*f[6]*dv2*wx1*wx2+0.4714045207910317*f[15]*dv1*dv2*wx2+2.828427124746191*f[1]*wx1*wx1_sq+2.449489742783178*f[5]*dv1*wx1_sq+0.2108185106778921*f[28]*dv3_sq*wx1+0.2357022603955158*f[1]*dv3_sq*wx1+0.2108185106778921*f[23]*dv2_sq*wx1+0.2357022603955158*f[1]*dv2_sq*wx1+0.632455532033676*f[20]*dv1_sq*wx1+0.7071067811865475*f[1]*dv1_sq*wx1+0.06085806194501846*f[41]*dv1*dv3_sq+0.06804138174397717*f[5]*dv1*dv3_sq+0.06085806194501846*f[34]*dv1*dv2_sq+0.06804138174397717*f[5]*dv1*dv2_sq+0.1224744871391589*f[5]*dv1*dv1_sq)*volFact; 
  out[2] += (2.828427124746191*f[11]*wx1*wx3_sq+0.816496580927726*f[19]*dv1*wx3_sq+1.632993161855453*f[25]*dv3*wx1*wx3+0.4714045207910317*f[35]*dv1*dv3*wx3+2.828427124746191*f[11]*wx1*wx2_sq+0.816496580927726*f[19]*dv1*wx2_sq+1.632993161855453*f[21]*dv2*wx1*wx2+0.4714045207910317*f[32]*dv1*dv2*wx2+2.828427124746191*f[11]*wx1*wx1_sq+2.449489742783178*f[19]*dv1*wx1_sq+0.210818510677892*f[47]*dv3_sq*wx1+0.2357022603955158*f[11]*dv3_sq*wx1+0.210818510677892*f[45]*dv2_sq*wx1+0.2357022603955158*f[11]*dv2_sq*wx1+0.6324555320336759*f[44]*dv1_sq*wx1+0.7071067811865475*f[11]*dv1_sq*wx1+0.06085806194501846*f[60]*dv1*dv3_sq+0.06804138174397714*f[19]*dv1*dv3_sq+0.06085806194501846*f[55]*dv1*dv2_sq+0.06804138174397714*f[19]*dv1*dv2_sq+0.1224744871391588*f[19]*dv1*dv1_sq)*volFact; 
  out[3] += (2.828427124746191*f[0]*wx2*wx3_sq+0.8164965809277261*f[3]*dv2*wx3_sq+1.632993161855453*f[4]*dv3*wx2*wx3+0.4714045207910317*f[10]*dv2*dv3*wx3+2.828427124746191*f[0]*wx2*wx2_sq+2.449489742783178*f[3]*dv2*wx2_sq+2.828427124746191*f[0]*wx1_sq*wx2+1.632993161855453*f[2]*dv1*wx1*wx2+0.210818510677892*f[14]*dv3_sq*wx2+0.2357022603955158*f[0]*dv3_sq*wx2+0.6324555320336759*f[13]*dv2_sq*wx2+0.7071067811865475*f[0]*dv2_sq*wx2+0.210818510677892*f[12]*dv1_sq*wx2+0.2357022603955158*f[0]*dv1_sq*wx2+0.8164965809277261*f[3]*dv2*wx1_sq+0.4714045207910317*f[7]*dv1*dv2*wx1+0.06085806194501844*f[30]*dv2*dv3_sq+0.06804138174397717*f[3]*dv2*dv3_sq+0.1224744871391589*f[3]*dv2*dv2_sq+0.06085806194501844*f[22]*dv1_sq*dv2+0.06804138174397717*f[3]*dv1_sq*dv2)*volFact; 
  out[4] += (2.828427124746191*f[1]*wx2*wx3_sq+0.8164965809277261*f[6]*dv2*wx3_sq+1.632993161855453*f[8]*dv3*wx2*wx3+0.4714045207910317*f[17]*dv2*dv3*wx3+2.828427124746191*f[1]*wx2*wx2_sq+2.449489742783178*f[6]*dv2*wx2_sq+2.828427124746191*f[1]*wx1_sq*wx2+1.632993161855453*f[5]*dv1*wx1*wx2+0.2108185106778921*f[28]*dv3_sq*wx2+0.2357022603955158*f[1]*dv3_sq*wx2+0.632455532033676*f[23]*dv2_sq*wx2+0.7071067811865475*f[1]*dv2_sq*wx2+0.2108185106778921*f[20]*dv1_sq*wx2+0.2357022603955158*f[1]*dv1_sq*wx2+0.8164965809277261*f[6]*dv2*wx1_sq+0.4714045207910317*f[15]*dv1*dv2*wx1+0.06085806194501846*f[42]*dv2*dv3_sq+0.06804138174397717*f[6]*dv2*dv3_sq+0.1224744871391589*f[6]*dv2*dv2_sq+0.06085806194501846*f[33]*dv1_sq*dv2+0.06804138174397717*f[6]*dv1_sq*dv2)*volFact; 
  out[5] += (2.828427124746191*f[11]*wx2*wx3_sq+0.816496580927726*f[21]*dv2*wx3_sq+1.632993161855453*f[25]*dv3*wx2*wx3+0.4714045207910317*f[37]*dv2*dv3*wx3+2.828427124746191*f[11]*wx2*wx2_sq+2.449489742783178*f[21]*dv2*wx2_sq+2.828427124746191*f[11]*wx1_sq*wx2+1.632993161855453*f[19]*dv1*wx1*wx2+0.210818510677892*f[47]*dv3_sq*wx2+0.2357022603955158*f[11]*dv3_sq*wx2+0.6324555320336759*f[45]*dv2_sq*wx2+0.7071067811865475*f[11]*dv2_sq*wx2+0.210818510677892*f[44]*dv1_sq*wx2+0.2357022603955158*f[11]*dv1_sq*wx2+0.816496580927726*f[21]*dv2*wx1_sq+0.4714045207910317*f[32]*dv1*dv2*wx1+0.06085806194501846*f[62]*dv2*dv3_sq+0.06804138174397714*f[21]*dv2*dv3_sq+0.1224744871391588*f[21]*dv2*dv2_sq+0.06085806194501846*f[54]*dv1_sq*dv2+0.06804138174397714*f[21]*dv1_sq*dv2)*volFact; 
  out[6] += (2.828427124746191*f[0]*wx3*wx3_sq+2.449489742783178*f[4]*dv3*wx3_sq+2.828427124746191*f[0]*wx2_sq*wx3+1.632993161855453*f[3]*dv2*wx2*wx3+2.828427124746191*f[0]*wx1_sq*wx3+1.632993161855453*f[2]*dv1*wx1*wx3+0.6324555320336759*f[14]*dv3_sq*wx3+0.7071067811865475*f[0]*dv3_sq*wx3+0.210818510677892*f[13]*dv2_sq*wx3+0.2357022603955158*f[0]*dv2_sq*wx3+0.210818510677892*f[12]*dv1_sq*wx3+0.2357022603955158*f[0]*dv1_sq*wx3+0.8164965809277261*f[4]*dv3*wx2_sq+0.4714045207910317*f[10]*dv2*dv3*wx2+0.8164965809277261*f[4]*dv3*wx1_sq+0.4714045207910317*f[9]*dv1*dv3*wx1+0.1224744871391589*f[4]*dv3*dv3_sq+0.06085806194501844*f[27]*dv2_sq*dv3+0.06804138174397717*f[4]*dv2_sq*dv3+0.06085806194501844*f[26]*dv1_sq*dv3+0.06804138174397717*f[4]*dv1_sq*dv3)*volFact; 
  out[7] += (2.828427124746191*f[1]*wx3*wx3_sq+2.449489742783178*f[8]*dv3*wx3_sq+2.828427124746191*f[1]*wx2_sq*wx3+1.632993161855453*f[6]*dv2*wx2*wx3+2.828427124746191*f[1]*wx1_sq*wx3+1.632993161855453*f[5]*dv1*wx1*wx3+0.632455532033676*f[28]*dv3_sq*wx3+0.7071067811865475*f[1]*dv3_sq*wx3+0.2108185106778921*f[23]*dv2_sq*wx3+0.2357022603955158*f[1]*dv2_sq*wx3+0.2108185106778921*f[20]*dv1_sq*wx3+0.2357022603955158*f[1]*dv1_sq*wx3+0.8164965809277261*f[8]*dv3*wx2_sq+0.4714045207910317*f[17]*dv2*dv3*wx2+0.8164965809277261*f[8]*dv3*wx1_sq+0.4714045207910317*f[16]*dv1*dv3*wx1+0.1224744871391589*f[8]*dv3*dv3_sq+0.06085806194501846*f[39]*dv2_sq*dv3+0.06804138174397717*f[8]*dv2_sq*dv3+0.06085806194501846*f[36]*dv1_sq*dv3+0.06804138174397717*f[8]*dv1_sq*dv3)*volFact; 
  out[8] += (2.828427124746191*f[11]*wx3*wx3_sq+2.449489742783178*f[25]*dv3*wx3_sq+2.828427124746191*f[11]*wx2_sq*wx3+1.632993161855453*f[21]*dv2*wx2*wx3+2.828427124746191*f[11]*wx1_sq*wx3+1.632993161855453*f[19]*dv1*wx1*wx3+0.6324555320336759*f[47]*dv3_sq*wx3+0.7071067811865475*f[11]*dv3_sq*wx3+0.210818510677892*f[45]*dv2_sq*wx3+0.2357022603955158*f[11]*dv2_sq*wx3+0.210818510677892*f[44]*dv1_sq*wx3+0.2357022603955158*f[11]*dv1_sq*wx3+0.816496580927726*f[25]*dv3*wx2_sq+0.4714045207910317*f[37]*dv2*dv3*wx2+0.816496580927726*f[25]*dv3*wx1_sq+0.4714045207910317*f[35]*dv1*dv3*wx1+0.1224744871391588*f[25]*dv3*dv3_sq+0.06085806194501846*f[58]*dv2_sq*dv3+0.06804138174397714*f[25]*dv2_sq*dv3+0.06085806194501846*f[57]*dv1_sq*dv3+0.06804138174397714*f[25]*dv1_sq*dv3)*volFact; 
} 
void MomentCalc1x3vTensor_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double tempM0[2], tempM1i[6]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 

  tempM1i[0] = 0.8164965809277261*f[2]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.8164965809277261*f[5]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.8164965809277261*f[3]*dv2*volFact+tempM0[0]*wx2; 
  tempM1i[3] = 0.8164965809277261*f[6]*dv2*volFact+tempM0[1]*wx2; 
  tempM1i[4] = 0.8164965809277261*f[4]*dv3*volFact+tempM0[0]*wx3; 
  tempM1i[5] = 0.8164965809277261*f[8]*dv3*volFact+tempM0[1]*wx3; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM1i[2] += tempM1i[2]; 
  outM1i[3] += tempM1i[3]; 
  outM1i[4] += tempM1i[4]; 
  outM1i[5] += tempM1i[5]; 
  outM2[0] += (0.2357022603955158*f[0]*dv3_sq+0.2357022603955158*f[0]*dv2_sq+0.2357022603955158*f[0]*dv1_sq)*volFact+tempM0[0]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[4]*wx3+2.0*tempM1i[2]*wx2+2.0*tempM1i[0]*wx1; 
  outM2[1] += (0.2357022603955158*f[1]*dv3_sq+0.2357022603955158*f[1]*dv2_sq+0.2357022603955158*f[1]*dv1_sq)*volFact+tempM0[1]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[5]*wx3+2.0*tempM1i[3]*wx2+2.0*tempM1i[1]*wx1; 
} 
void MomentCalc1x3vTensor_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[1]*dxv[2]*dxv[3]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double tempM0[3], tempM1i[9]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 
  tempM0[2] = 2.828427124746191*f[11]*volFact; 

  tempM1i[0] = 0.8164965809277261*f[2]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.8164965809277261*f[5]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.816496580927726*f[19]*dv1*volFact+tempM0[2]*wx1; 
  tempM1i[3] = 0.8164965809277261*f[3]*dv2*volFact+tempM0[0]*wx2; 
  tempM1i[4] = 0.8164965809277261*f[6]*dv2*volFact+tempM0[1]*wx2; 
  tempM1i[5] = 0.816496580927726*f[21]*dv2*volFact+tempM0[2]*wx2; 
  tempM1i[6] = 0.8164965809277261*f[4]*dv3*volFact+tempM0[0]*wx3; 
  tempM1i[7] = 0.8164965809277261*f[8]*dv3*volFact+tempM0[1]*wx3; 
  tempM1i[8] = 0.816496580927726*f[25]*dv3*volFact+tempM0[2]*wx3; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM0[2] += tempM0[2]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM1i[2] += tempM1i[2]; 
  outM1i[3] += tempM1i[3]; 
  outM1i[4] += tempM1i[4]; 
  outM1i[5] += tempM1i[5]; 
  outM1i[6] += tempM1i[6]; 
  outM1i[7] += tempM1i[7]; 
  outM1i[8] += tempM1i[8]; 
  outM2[0] += (0.210818510677892*f[14]*dv3_sq+0.2357022603955158*f[0]*dv3_sq+0.210818510677892*f[13]*dv2_sq+0.2357022603955158*f[0]*dv2_sq+0.210818510677892*f[12]*dv1_sq+0.2357022603955158*f[0]*dv1_sq)*volFact+tempM0[0]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[6]*wx3+2.0*tempM1i[3]*wx2+2.0*tempM1i[0]*wx1; 
  outM2[1] += (0.2108185106778921*f[28]*dv3_sq+0.2357022603955158*f[1]*dv3_sq+0.2108185106778921*f[23]*dv2_sq+0.2357022603955158*f[1]*dv2_sq+0.2108185106778921*f[20]*dv1_sq+0.2357022603955158*f[1]*dv1_sq)*volFact+tempM0[1]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[7]*wx3+2.0*tempM1i[4]*wx2+2.0*tempM1i[1]*wx1; 
  outM2[2] += (0.210818510677892*f[47]*dv3_sq+0.2357022603955158*f[11]*dv3_sq+0.210818510677892*f[45]*dv2_sq+0.2357022603955158*f[11]*dv2_sq+0.210818510677892*f[44]*dv1_sq+0.2357022603955158*f[11]*dv1_sq)*volFact+tempM0[2]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[8]*wx3+2.0*tempM1i[5]*wx2+2.0*tempM1i[2]*wx1; 
} 
