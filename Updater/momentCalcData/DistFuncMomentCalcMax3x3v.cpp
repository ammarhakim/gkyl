#include <DistFuncMomentCalcModDecl.h> 
__host__ __device__ void MomentCalc3x3vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[3]*volFact; 
} 
__host__ __device__ void MomentCalc3x3vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[4]*dv1); 
  out[1] += 2.828427124746191*f[1]*volFact*wx1; 
  out[2] += 2.828427124746191*f[2]*volFact*wx1; 
  out[3] += 2.828427124746191*f[3]*volFact*wx1; 
  out[4] += volFact*(2.828427124746191*f[0]*wx2+0.8164965809277261*f[5]*dv2); 
  out[5] += 2.828427124746191*f[1]*volFact*wx2; 
  out[6] += 2.828427124746191*f[2]*volFact*wx2; 
  out[7] += 2.828427124746191*f[3]*volFact*wx2; 
  out[8] += volFact*(2.828427124746191*f[0]*wx3+0.8164965809277261*f[6]*dv3); 
  out[9] += 2.828427124746191*f[1]*volFact*wx3; 
  out[10] += 2.828427124746191*f[2]*volFact*wx3; 
  out[11] += 2.828427124746191*f[3]*volFact*wx3; 
} 
__host__ __device__ void MomentCalc3x3vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += volFact*(2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[4]*dv1*wx1+0.2357022603955158*f[0]*dv1_sq); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1_sq+0.2357022603955158*f[1]*dv1_sq); 
  out[2] += volFact*(2.828427124746191*f[2]*wx1_sq+0.2357022603955158*f[2]*dv1_sq); 
  out[3] += volFact*(2.828427124746191*f[3]*wx1_sq+0.2357022603955158*f[3]*dv1_sq); 
  out[4] += volFact*(2.828427124746191*f[0]*wx1*wx2+0.8164965809277261*f[4]*dv1*wx2+0.8164965809277261*f[5]*dv2*wx1); 
  out[5] += 2.828427124746191*f[1]*volFact*wx1*wx2; 
  out[6] += 2.828427124746191*f[2]*volFact*wx1*wx2; 
  out[7] += 2.828427124746191*f[3]*volFact*wx1*wx2; 
  out[8] += volFact*(2.828427124746191*f[0]*wx1*wx3+0.8164965809277261*f[4]*dv1*wx3+0.8164965809277261*f[6]*dv3*wx1); 
  out[9] += 2.828427124746191*f[1]*volFact*wx1*wx3; 
  out[10] += 2.828427124746191*f[2]*volFact*wx1*wx3; 
  out[11] += 2.828427124746191*f[3]*volFact*wx1*wx3; 
  out[12] += volFact*(2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[5]*dv2*wx2+0.2357022603955158*f[0]*dv2_sq); 
  out[13] += volFact*(2.828427124746191*f[1]*wx2_sq+0.2357022603955158*f[1]*dv2_sq); 
  out[14] += volFact*(2.828427124746191*f[2]*wx2_sq+0.2357022603955158*f[2]*dv2_sq); 
  out[15] += volFact*(2.828427124746191*f[3]*wx2_sq+0.2357022603955158*f[3]*dv2_sq); 
  out[16] += volFact*(2.828427124746191*f[0]*wx2*wx3+0.8164965809277261*f[5]*dv2*wx3+0.8164965809277261*f[6]*dv3*wx2); 
  out[17] += 2.828427124746191*f[1]*volFact*wx2*wx3; 
  out[18] += 2.828427124746191*f[2]*volFact*wx2*wx3; 
  out[19] += 2.828427124746191*f[3]*volFact*wx2*wx3; 
  out[20] += volFact*(2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[6]*dv3*wx3+0.2357022603955158*f[0]*dv3_sq); 
  out[21] += volFact*(2.828427124746191*f[1]*wx3_sq+0.2357022603955158*f[1]*dv3_sq); 
  out[22] += volFact*(2.828427124746191*f[2]*wx3_sq+0.2357022603955158*f[2]*dv3_sq); 
  out[23] += volFact*(2.828427124746191*f[3]*wx3_sq+0.2357022603955158*f[3]*dv3_sq); 
} 
__host__ __device__ void MomentCalc3x3vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  out[0] += volFact*(2.828427124746191*f[0]*wx3_sq+1.632993161855453*f[6]*dv3*wx3+2.828427124746191*f[0]*wx2_sq+1.632993161855453*f[5]*dv2*wx2+2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[4]*dv1*wx1+0.2357022603955158*f[0]*dv3_sq+0.2357022603955158*f[0]*dv2_sq+0.2357022603955158*f[0]*dv1_sq); 
  out[1] += volFact*(2.828427124746191*f[1]*wx3_sq+2.828427124746191*f[1]*wx2_sq+2.828427124746191*f[1]*wx1_sq+0.2357022603955158*f[1]*dv3_sq+0.2357022603955158*f[1]*dv2_sq+0.2357022603955158*f[1]*dv1_sq); 
  out[2] += volFact*(2.828427124746191*f[2]*wx3_sq+2.828427124746191*f[2]*wx2_sq+2.828427124746191*f[2]*wx1_sq+0.2357022603955158*f[2]*dv3_sq+0.2357022603955158*f[2]*dv2_sq+0.2357022603955158*f[2]*dv1_sq); 
  out[3] += volFact*(2.828427124746191*f[3]*wx3_sq+2.828427124746191*f[3]*wx2_sq+2.828427124746191*f[3]*wx1_sq+0.2357022603955158*f[3]*dv3_sq+0.2357022603955158*f[3]*dv2_sq+0.2357022603955158*f[3]*dv1_sq); 
} 
__host__ __device__ void MomentCalc3x3vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out) 
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
  out[0] += volFact*(2.828427124746191*f[0]*wx1*wx3_sq+0.8164965809277261*f[4]*dv1*wx3_sq+1.632993161855453*f[6]*dv3*wx1*wx3+2.828427124746191*f[0]*wx1*wx2_sq+0.8164965809277261*f[4]*dv1*wx2_sq+1.632993161855453*f[5]*dv2*wx1*wx2+2.828427124746191*f[0]*wx1*wx1_sq+2.449489742783178*f[4]*dv1*wx1_sq+0.2357022603955158*f[0]*dv3_sq*wx1+0.2357022603955158*f[0]*dv2_sq*wx1+0.7071067811865475*f[0]*dv1_sq*wx1+0.06804138174397717*f[4]*dv1*dv3_sq+0.06804138174397717*f[4]*dv1*dv2_sq+0.1224744871391589*f[4]*dv1*dv1_sq); 
  out[1] += volFact*(2.828427124746191*f[1]*wx1*wx3_sq+2.828427124746191*f[1]*wx1*wx2_sq+2.828427124746191*f[1]*wx1*wx1_sq+0.2357022603955158*f[1]*dv3_sq*wx1+0.2357022603955158*f[1]*dv2_sq*wx1+0.7071067811865475*f[1]*dv1_sq*wx1); 
  out[2] += volFact*(2.828427124746191*f[2]*wx1*wx3_sq+2.828427124746191*f[2]*wx1*wx2_sq+2.828427124746191*f[2]*wx1*wx1_sq+0.2357022603955158*f[2]*dv3_sq*wx1+0.2357022603955158*f[2]*dv2_sq*wx1+0.7071067811865475*f[2]*dv1_sq*wx1); 
  out[3] += volFact*(2.828427124746191*f[3]*wx1*wx3_sq+2.828427124746191*f[3]*wx1*wx2_sq+2.828427124746191*f[3]*wx1*wx1_sq+0.2357022603955158*f[3]*dv3_sq*wx1+0.2357022603955158*f[3]*dv2_sq*wx1+0.7071067811865475*f[3]*dv1_sq*wx1); 
  out[4] += volFact*(2.828427124746191*f[0]*wx2*wx3_sq+0.8164965809277261*f[5]*dv2*wx3_sq+1.632993161855453*f[6]*dv3*wx2*wx3+2.828427124746191*f[0]*wx2*wx2_sq+2.449489742783178*f[5]*dv2*wx2_sq+2.828427124746191*f[0]*wx1_sq*wx2+1.632993161855453*f[4]*dv1*wx1*wx2+0.2357022603955158*f[0]*dv3_sq*wx2+0.7071067811865475*f[0]*dv2_sq*wx2+0.2357022603955158*f[0]*dv1_sq*wx2+0.8164965809277261*f[5]*dv2*wx1_sq+0.06804138174397717*f[5]*dv2*dv3_sq+0.1224744871391589*f[5]*dv2*dv2_sq+0.06804138174397717*f[5]*dv1_sq*dv2); 
  out[5] += volFact*(2.828427124746191*f[1]*wx2*wx3_sq+2.828427124746191*f[1]*wx2*wx2_sq+2.828427124746191*f[1]*wx1_sq*wx2+0.2357022603955158*f[1]*dv3_sq*wx2+0.7071067811865475*f[1]*dv2_sq*wx2+0.2357022603955158*f[1]*dv1_sq*wx2); 
  out[6] += volFact*(2.828427124746191*f[2]*wx2*wx3_sq+2.828427124746191*f[2]*wx2*wx2_sq+2.828427124746191*f[2]*wx1_sq*wx2+0.2357022603955158*f[2]*dv3_sq*wx2+0.7071067811865475*f[2]*dv2_sq*wx2+0.2357022603955158*f[2]*dv1_sq*wx2); 
  out[7] += volFact*(2.828427124746191*f[3]*wx2*wx3_sq+2.828427124746191*f[3]*wx2*wx2_sq+2.828427124746191*f[3]*wx1_sq*wx2+0.2357022603955158*f[3]*dv3_sq*wx2+0.7071067811865475*f[3]*dv2_sq*wx2+0.2357022603955158*f[3]*dv1_sq*wx2); 
  out[8] += volFact*(2.828427124746191*f[0]*wx3*wx3_sq+2.449489742783178*f[6]*dv3*wx3_sq+2.828427124746191*f[0]*wx2_sq*wx3+1.632993161855453*f[5]*dv2*wx2*wx3+2.828427124746191*f[0]*wx1_sq*wx3+1.632993161855453*f[4]*dv1*wx1*wx3+0.7071067811865475*f[0]*dv3_sq*wx3+0.2357022603955158*f[0]*dv2_sq*wx3+0.2357022603955158*f[0]*dv1_sq*wx3+0.8164965809277261*f[6]*dv3*wx2_sq+0.8164965809277261*f[6]*dv3*wx1_sq+0.1224744871391589*f[6]*dv3*dv3_sq+0.06804138174397717*f[6]*dv2_sq*dv3+0.06804138174397717*f[6]*dv1_sq*dv3); 
  out[9] += volFact*(2.828427124746191*f[1]*wx3*wx3_sq+2.828427124746191*f[1]*wx2_sq*wx3+2.828427124746191*f[1]*wx1_sq*wx3+0.7071067811865475*f[1]*dv3_sq*wx3+0.2357022603955158*f[1]*dv2_sq*wx3+0.2357022603955158*f[1]*dv1_sq*wx3); 
  out[10] += volFact*(2.828427124746191*f[2]*wx3*wx3_sq+2.828427124746191*f[2]*wx2_sq*wx3+2.828427124746191*f[2]*wx1_sq*wx3+0.7071067811865475*f[2]*dv3_sq*wx3+0.2357022603955158*f[2]*dv2_sq*wx3+0.2357022603955158*f[2]*dv1_sq*wx3); 
  out[11] += volFact*(2.828427124746191*f[3]*wx3*wx3_sq+2.828427124746191*f[3]*wx2_sq*wx3+2.828427124746191*f[3]*wx1_sq*wx3+0.7071067811865475*f[3]*dv3_sq*wx3+0.2357022603955158*f[3]*dv2_sq*wx3+0.2357022603955158*f[3]*dv1_sq*wx3); 
} 
__host__ __device__ void MomentCalc3x3vMax_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[3]*dxv[4]*dxv[5]/8; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
  double tempM0[4], tempM1i[12]; 

  tempM0[0] = 2.828427124746191*f[0]*volFact; 
  tempM0[1] = 2.828427124746191*f[1]*volFact; 
  tempM0[2] = 2.828427124746191*f[2]*volFact; 
  tempM0[3] = 2.828427124746191*f[3]*volFact; 

  tempM1i[0] = tempM0[0]*wx1+0.8164965809277261*f[4]*dv1*volFact; 
  tempM1i[1] = tempM0[1]*wx1; 
  tempM1i[2] = tempM0[2]*wx1; 
  tempM1i[3] = tempM0[3]*wx1; 
  tempM1i[4] = tempM0[0]*wx2+0.8164965809277261*f[5]*dv2*volFact; 
  tempM1i[5] = tempM0[1]*wx2; 
  tempM1i[6] = tempM0[2]*wx2; 
  tempM1i[7] = tempM0[3]*wx2; 
  tempM1i[8] = tempM0[0]*wx3+0.8164965809277261*f[6]*dv3*volFact; 
  tempM1i[9] = tempM0[1]*wx3; 
  tempM1i[10] = tempM0[2]*wx3; 
  tempM1i[11] = tempM0[3]*wx3; 

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
  outM2[0] += tempM0[0]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[8]*wx3+2.0*tempM1i[4]*wx2+2.0*tempM1i[0]*wx1+(0.2357022603955158*f[0]*dv3_sq+0.2357022603955158*f[0]*dv2_sq+0.2357022603955158*f[0]*dv1_sq)*volFact; 
  outM2[1] += tempM0[1]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[9]*wx3+2.0*tempM1i[5]*wx2+2.0*tempM1i[1]*wx1+(0.2357022603955158*f[1]*dv3_sq+0.2357022603955158*f[1]*dv2_sq+0.2357022603955158*f[1]*dv1_sq)*volFact; 
  outM2[2] += tempM0[2]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[10]*wx3+2.0*tempM1i[6]*wx2+2.0*tempM1i[2]*wx1+(0.2357022603955158*f[2]*dv3_sq+0.2357022603955158*f[2]*dv2_sq+0.2357022603955158*f[2]*dv1_sq)*volFact; 
  outM2[3] += tempM0[3]*((-1.0*wx3_sq)-1.0*wx2_sq-1.0*wx1_sq)+2.0*tempM1i[11]*wx3+2.0*tempM1i[7]*wx2+2.0*tempM1i[3]*wx1+(0.2357022603955158*f[3]*dv3_sq+0.2357022603955158*f[3]*dv2_sq+0.2357022603955158*f[3]*dv1_sq)*volFact; 
} 
