#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc1x3vSer_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]/16; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 4.0*f[0]*volFact; 
  out[4] += 4.0*f[0]*wx3_sq*volFact+2.309401076758503*f[4]*dv3*wx3*volFact+4.0*f[0]*wx2_sq*volFact+2.309401076758503*f[3]*dv2*wx2*volFact+4.0*f[0]*wx1_sq*volFact+2.309401076758503*f[2]*dv1*wx1*volFact+0.3333333333333333*f[0]*dv3_sq*volFact+0.3333333333333333*f[0]*dv2_sq*volFact+0.3333333333333333*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x3vSer_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]/16; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 4.0*f[0]*volFact; 
  out[4] += 4.0*f[0]*wx3_sq*volFact+2.309401076758503*f[4]*dv3*wx3*volFact+4.0*f[0]*wx2_sq*volFact+2.309401076758503*f[3]*dv2*wx2*volFact+4.0*f[0]*wx1_sq*volFact+2.309401076758503*f[2]*dv1*wx1*volFact+0.2981423969999719*f[14]*dv3_sq*volFact+0.3333333333333333*f[0]*dv3_sq*volFact+0.2981423969999719*f[13]*dv2_sq*volFact+0.3333333333333333*f[0]*dv2_sq*volFact+0.2981423969999719*f[12]*dv1_sq*volFact+0.3333333333333333*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x3vSer_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]/16; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 4.0*f[0]*volFact; 
  out[4] += 4.0*f[0]*wx3_sq*volFact+2.309401076758503*f[4]*dv3*wx3*volFact+4.0*f[0]*wx2_sq*volFact+2.309401076758503*f[3]*dv2*wx2*volFact+4.0*f[0]*wx1_sq*volFact+2.309401076758503*f[2]*dv1*wx1*volFact+0.2981423969999719*f[14]*dv3_sq*volFact+0.3333333333333333*f[0]*dv3_sq*volFact+0.2981423969999719*f[13]*dv2_sq*volFact+0.3333333333333333*f[0]*dv2_sq*volFact+0.2981423969999719*f[12]*dv1_sq*volFact+0.3333333333333333*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x3vSer_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]/16; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 4.0*f[0]*volFact; 
  out[4] += 4.0*f[0]*wx3_sq*volFact+2.309401076758503*f[4]*dv3*wx3*volFact+4.0*f[0]*wx2_sq*volFact+2.309401076758503*f[3]*dv2*wx2*volFact+4.0*f[0]*wx1_sq*volFact+2.309401076758503*f[2]*dv1*wx1*volFact+0.2981423969999719*f[14]*dv3_sq*volFact+0.3333333333333333*f[0]*dv3_sq*volFact+0.2981423969999719*f[13]*dv2_sq*volFact+0.3333333333333333*f[0]*dv2_sq*volFact+0.2981423969999719*f[12]*dv1_sq*volFact+0.3333333333333333*f[0]*dv1_sq*volFact; 
} 
