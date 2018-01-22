#include <DistFuncMomentCalcModDecl.h> 
void MomentCalc1x1vMax_M0_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
} 
void MomentCalc1x1vMax_M0_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact; 
} 
void MomentCalc1x1vMax_M1i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += 1.414213562373095*f[0]*volFact*wx1+0.408248290463863*f[2]*dv1*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact*wx1; 
} 
void MomentCalc1x1vMax_M1i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  out[0] += 1.414213562373095*f[0]*volFact*wx1+0.408248290463863*f[2]*dv1*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact*wx1+0.408248290463863*f[3]*dv1*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact*wx1; 
} 
void MomentCalc1x1vMax_M2ij_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += 1.414213562373095*f[0]*volFact*wx1_sq+0.8164965809277261*f[2]*dv1*volFact*wx1+0.1178511301977579*f[0]*dv1_sq*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact*wx1_sq+0.1178511301977579*f[1]*dv1_sq*volFact; 
} 
void MomentCalc1x1vMax_M2ij_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += 1.414213562373095*f[0]*volFact*wx1_sq+0.8164965809277261*f[2]*dv1*volFact*wx1+0.105409255338946*f[5]*dv1_sq*volFact+0.1178511301977579*f[0]*dv1_sq*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact*wx1_sq+0.8164965809277261*f[3]*dv1*volFact*wx1+0.1178511301977579*f[1]*dv1_sq*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact*wx1_sq+0.1178511301977579*f[4]*dv1_sq*volFact; 
} 
void MomentCalc1x1vMax_M2_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += 1.414213562373095*f[0]*volFact*wx1_sq+0.8164965809277261*f[2]*dv1*volFact*wx1+0.1178511301977579*f[0]*dv1_sq*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact*wx1_sq+0.1178511301977579*f[1]*dv1_sq*volFact; 
} 
void MomentCalc1x1vMax_M2_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  out[0] += 1.414213562373095*f[0]*volFact*wx1_sq+0.8164965809277261*f[2]*dv1*volFact*wx1+0.105409255338946*f[5]*dv1_sq*volFact+0.1178511301977579*f[0]*dv1_sq*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact*wx1_sq+0.8164965809277261*f[3]*dv1*volFact*wx1+0.1178511301977579*f[1]*dv1_sq*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact*wx1_sq+0.1178511301977579*f[4]*dv1_sq*volFact; 
} 
void MomentCalc1x1vMax_M3i_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += 1.414213562373095*f[0]*volFact*wx1*wx1_sq+1.224744871391589*f[2]*dv1*volFact*wx1_sq+0.3535533905932737*f[0]*dv1_sq*volFact*wx1+0.06123724356957942*f[2]*dv1*dv1_sq*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact*wx1*wx1_sq+0.3535533905932737*f[1]*dv1_sq*volFact*wx1; 
} 
void MomentCalc1x1vMax_M3i_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[1]/2; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  out[0] += 1.414213562373095*f[0]*volFact*wx1*wx1_sq+1.224744871391589*f[2]*dv1*volFact*wx1_sq+0.3162277660168379*f[5]*dv1_sq*volFact*wx1+0.3535533905932737*f[0]*dv1_sq*volFact*wx1+0.06123724356957942*f[2]*dv1*dv1_sq*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact*wx1*wx1_sq+1.224744871391589*f[3]*dv1*volFact*wx1_sq+0.3535533905932737*f[1]*dv1_sq*volFact*wx1+0.06123724356957942*f[3]*dv1*dv1_sq*volFact; 
  out[2] += 1.414213562373095*f[4]*volFact*wx1*wx1_sq+0.3535533905932737*f[4]*dv1_sq*volFact*wx1; 
} 
