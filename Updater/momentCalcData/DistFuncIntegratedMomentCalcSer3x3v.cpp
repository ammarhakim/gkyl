#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc3x3vSer_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*dxv[5]/64; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 8.0*f[0]*volFact; 
  out[4] += 8.0*f[0]*wx3_sq*volFact+4.618802153517007*f[6]*dv3*wx3*volFact+8.0*f[0]*wx2_sq*volFact+4.618802153517007*f[5]*dv2*wx2*volFact+8.0*f[0]*wx1_sq*volFact+4.618802153517007*f[4]*dv1*wx1*volFact+0.6666666666666666*f[0]*dv3_sq*volFact+0.6666666666666666*f[0]*dv2_sq*volFact+0.6666666666666666*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc3x3vSer_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*dxv[5]/64; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 8.0*f[0]*volFact; 
  out[4] += 8.0*f[0]*wx3_sq*volFact+4.618802153517007*f[6]*dv3*wx3*volFact+8.0*f[0]*wx2_sq*volFact+4.618802153517007*f[5]*dv2*wx2*volFact+8.0*f[0]*wx1_sq*volFact+4.618802153517007*f[4]*dv1*wx1*volFact+0.5962847939999438*f[27]*dv3_sq*volFact+0.6666666666666666*f[0]*dv3_sq*volFact+0.5962847939999438*f[26]*dv2_sq*volFact+0.6666666666666666*f[0]*dv2_sq*volFact+0.5962847939999438*f[25]*dv1_sq*volFact+0.6666666666666666*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc3x3vSer_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*dxv[5]/64; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 8.0*f[0]*volFact; 
  out[4] += 8.0*f[0]*wx3_sq*volFact+4.618802153517007*f[6]*dv3*wx3*volFact+8.0*f[0]*wx2_sq*volFact+4.618802153517007*f[5]*dv2*wx2*volFact+8.0*f[0]*wx1_sq*volFact+4.618802153517007*f[4]*dv1*wx1*volFact+0.5962847939999438*f[27]*dv3_sq*volFact+0.6666666666666666*f[0]*dv3_sq*volFact+0.5962847939999438*f[26]*dv2_sq*volFact+0.6666666666666666*f[0]*dv2_sq*volFact+0.5962847939999438*f[25]*dv1_sq*volFact+0.6666666666666666*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc3x3vSer_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*dxv[5]/64; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[5], dv3 = dxv[5]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 8.0*f[0]*volFact; 
  out[4] += 8.0*f[0]*wx3_sq*volFact+4.618802153517007*f[6]*dv3*wx3*volFact+8.0*f[0]*wx2_sq*volFact+4.618802153517007*f[5]*dv2*wx2*volFact+8.0*f[0]*wx1_sq*volFact+4.618802153517007*f[4]*dv1*wx1*volFact+0.5962847939999438*f[27]*dv3_sq*volFact+0.6666666666666666*f[0]*dv3_sq*volFact+0.5962847939999438*f[26]*dv2_sq*volFact+0.6666666666666666*f[0]*dv2_sq*volFact+0.5962847939999438*f[25]*dv1_sq*volFact+0.6666666666666666*f[0]*dv1_sq*volFact; 
} 
