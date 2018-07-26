#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc1x1vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
} 
void MomentCalc1x1vSer_M0_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact; 
} 
void MomentCalc1x1vSer_M0_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact; 
  out[3] += 1.414213562373095*f[8]*volFact; 
} 
void MomentCalc1x1vSer_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += (1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1)*volFact; 
} 
void MomentCalc1x1vSer_M1i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += (1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1+0.408248290463863*f[6]*dv1)*volFact; 
} 
void MomentCalc1x1vSer_M1i_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += (1.414213562373095*f[0]*wx1+0.408248290463863*f[2]*dv1)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1+0.408248290463863*f[3]*dv1)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1+0.408248290463863*f[6]*dv1)*volFact; 
  out[3] += (1.414213562373095*f[8]*wx1+0.4082482904638629*f[10]*dv1)*volFact; 
} 
void MomentCalc1x1vSer_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double zeta1[4]; 

  zeta1[0] = 2.0*wx1_sq+0.1666666666666667*dv1_sq; 
  zeta1[2] = 1.154700538379252*dv1*wx1; 
  out[0] += (0.7071067811865475*f[2]*zeta1[2]+0.7071067811865475*f[0]*zeta1[0])*volFact; 
  out[1] += (0.7071067811865475*zeta1[2]*f[3]+0.7071067811865475*zeta1[0]*f[1])*volFact; 
} 
void MomentCalc1x1vSer_M2ij_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq)*volFact; 
} 
void MomentCalc1x1vSer_M2ij_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq)*volFact; 
  out[3] += (1.414213562373095*f[8]*wx1_sq+0.8164965809277258*f[10]*dv1*wx1+0.1178511301977579*f[8]*dv1_sq)*volFact; 
} 
void MomentCalc1x1vSer_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double zeta[4]; 

  zeta[0] = 2.0*wx1_sq+0.1666666666666667*dv1_sq; 
  zeta[2] = 1.154700538379252*dv1*wx1; 
  out[0] += (0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  out[1] += (0.7071067811865475*zeta[2]*f[3]+0.7071067811865475*zeta[0]*f[1])*volFact; 
} 
void MomentCalc1x1vSer_M2_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq)*volFact; 
} 
void MomentCalc1x1vSer_M2_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += (1.414213562373095*f[0]*wx1_sq+0.8164965809277261*f[2]*dv1*wx1+0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact; 
  out[1] += (1.414213562373095*f[1]*wx1_sq+0.8164965809277261*f[3]*dv1*wx1+0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact; 
  out[2] += (1.414213562373095*f[4]*wx1_sq+0.816496580927726*f[6]*dv1*wx1+0.1178511301977579*f[4]*dv1_sq)*volFact; 
  out[3] += (1.414213562373095*f[8]*wx1_sq+0.8164965809277258*f[10]*dv1*wx1+0.1178511301977579*f[8]*dv1_sq)*volFact; 
} 
void MomentCalc1x1vSer_M3i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  double zeta1[4]; 

  zeta1[0] = 2.0*wx1_cu+0.5*dv1_sq*wx1; 
  zeta1[2] = 1.732050807568877*dv1*wx1_sq+0.08660254037844387*dv1_cu; 
  out[0] += (0.7071067811865475*f[2]*zeta1[2]+0.7071067811865475*f[0]*zeta1[0])*volFact; 
  out[1] += (0.7071067811865475*zeta1[2]*f[3]+0.7071067811865475*zeta1[0]*f[1])*volFact; 
} 
void MomentCalc1x1vSer_M3i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  double zeta1[8]; 

  zeta1[0] = 2.0*wx1_cu+0.5*dv1_sq*wx1; 
  zeta1[2] = 1.732050807568877*dv1*wx1_sq+0.08660254037844387*dv1_cu; 
  zeta1[5] = 0.4472135954999579*dv1_sq*wx1; 
  out[0] += (0.7071067811865475*f[5]*zeta1[5]+0.7071067811865475*f[2]*zeta1[2]+0.7071067811865475*f[0]*zeta1[0])*volFact; 
  out[1] += (0.7071067811865475*zeta1[5]*f[7]+0.7071067811865475*zeta1[2]*f[3]+0.7071067811865475*zeta1[0]*f[1])*volFact; 
  out[2] += (0.7071067811865475*zeta1[2]*f[6]+0.7071067811865475*zeta1[0]*f[4])*volFact; 
} 
void MomentCalc1x1vSer_M3i_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += (1.224744871391589*f[2]*dv1*wx1_sq+1.414213562373095*f[0]*wx1_cu+0.3162277660168379*f[5]*dv1_sq*wx1+0.3535533905932737*f[0]*dv1_sq*wx1+0.02672612419124243*f[9]*dv1_cu+0.06123724356957942*f[2]*dv1_cu)*volFact; 
  out[1] += (1.224744871391589*f[3]*dv1*wx1_sq+1.414213562373095*f[1]*wx1_cu+0.3162277660168379*f[7]*dv1_sq*wx1+0.3535533905932737*f[1]*dv1_sq*wx1+0.02672612419124244*f[11]*dv1_cu+0.06123724356957942*f[3]*dv1_cu)*volFact; 
  out[2] += (1.224744871391589*f[6]*dv1*wx1_sq+1.414213562373095*f[4]*wx1_cu+0.3535533905932737*f[4]*dv1_sq*wx1+0.06123724356957942*f[6]*dv1_cu)*volFact; 
  out[3] += (1.224744871391589*f[10]*dv1*wx1_sq+1.414213562373095*f[8]*wx1_cu+0.3535533905932737*f[8]*dv1_sq*wx1+0.06123724356957942*f[10]*dv1_cu)*volFact; 
} 
void MomentCalc1x1vSer_FiveMoments_P1(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double tempM0[2], tempM1i[2]; 

  double zeta[4]; 

  zeta[0] = 2.0*wx1_sq+0.1666666666666667*dv1_sq; 
  zeta[2] = 1.154700538379252*dv1*wx1; 
  tempM0[0] = 1.414213562373095*f[0]*volFact; 
  tempM0[1] = 1.414213562373095*f[1]*volFact; 

  tempM1i[0] = 0.408248290463863*f[2]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.408248290463863*f[3]*dv1*volFact+tempM0[1]*wx1; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM2[0] += (0.7071067811865475*f[2]*zeta[2]+0.7071067811865475*f[0]*zeta[0])*volFact; 
  outM2[1] += (0.7071067811865475*zeta[2]*f[3]+0.7071067811865475*zeta[0]*f[1])*volFact; 
} 
void MomentCalc1x1vSer_FiveMoments_P2(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double tempM0[3], tempM1i[3]; 

  tempM0[0] = 1.414213562373095*f[0]*volFact; 
  tempM0[1] = 1.414213562373095*f[1]*volFact; 
  tempM0[2] = 1.414213562373095*f[4]*volFact; 

  tempM1i[0] = 0.408248290463863*f[2]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.408248290463863*f[3]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.408248290463863*f[6]*dv1*volFact+tempM0[2]*wx1; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM0[2] += tempM0[2]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM1i[2] += tempM1i[2]; 
  outM2[0] += (0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact-1.0*tempM0[0]*wx1_sq+2.0*tempM1i[0]*wx1; 
  outM2[1] += (0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact-1.0*tempM0[1]*wx1_sq+2.0*tempM1i[1]*wx1; 
  outM2[2] += 0.1178511301977579*f[4]*dv1_sq*volFact-1.0*tempM0[2]*wx1_sq+2.0*tempM1i[2]*wx1; 
} 
void MomentCalc1x1vSer_FiveMoments_P3(const double *w, const double *dxv, const double *f, double *outM0, double *outM1i, double *outM2) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  double tempM0[4], tempM1i[4]; 

  tempM0[0] = 1.414213562373095*f[0]*volFact; 
  tempM0[1] = 1.414213562373095*f[1]*volFact; 
  tempM0[2] = 1.414213562373095*f[4]*volFact; 
  tempM0[3] = 1.414213562373095*f[8]*volFact; 

  tempM1i[0] = 0.408248290463863*f[2]*dv1*volFact+tempM0[0]*wx1; 
  tempM1i[1] = 0.408248290463863*f[3]*dv1*volFact+tempM0[1]*wx1; 
  tempM1i[2] = 0.408248290463863*f[6]*dv1*volFact+tempM0[2]*wx1; 
  tempM1i[3] = 0.4082482904638629*f[10]*dv1*volFact+tempM0[3]*wx1; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM0[2] += tempM0[2]; 
  outM0[3] += tempM0[3]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM1i[2] += tempM1i[2]; 
  outM1i[3] += tempM1i[3]; 
  outM2[0] += (0.105409255338946*f[5]*dv1_sq+0.1178511301977579*f[0]*dv1_sq)*volFact-1.0*tempM0[0]*wx1_sq+2.0*tempM1i[0]*wx1; 
  outM2[1] += (0.105409255338946*f[7]*dv1_sq+0.1178511301977579*f[1]*dv1_sq)*volFact-1.0*tempM0[1]*wx1_sq+2.0*tempM1i[1]*wx1; 
  outM2[2] += 0.1178511301977579*f[4]*dv1_sq*volFact-1.0*tempM0[2]*wx1_sq+2.0*tempM1i[2]*wx1; 
  outM2[3] += 0.1178511301977579*f[8]*dv1_sq*volFact-1.0*tempM0[3]*wx1_sq+2.0*tempM1i[3]*wx1; 
} 
