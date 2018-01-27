#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc2x3vSer_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/32; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[4] += 5.656854249492382*f[0]*wx3_sq*volFact+3.265986323710906*f[5]*dv3*wx3*volFact+5.656854249492382*f[0]*wx2_sq*volFact+3.265986323710906*f[4]*dv2*wx2*volFact+5.656854249492382*f[0]*wx1_sq*volFact+3.265986323710906*f[3]*dv1*wx1*volFact+0.4714045207910317*f[0]*dv3_sq*volFact+0.4714045207910317*f[0]*dv2_sq*volFact+0.4714045207910317*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc2x3vSer_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/32; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[4] += 5.656854249492382*f[0]*wx3_sq*volFact+3.265986323710906*f[5]*dv3*wx3*volFact+5.656854249492382*f[0]*wx2_sq*volFact+3.265986323710906*f[4]*dv2*wx2*volFact+5.656854249492382*f[0]*wx1_sq*volFact+3.265986323710906*f[3]*dv1*wx1*volFact+0.421637021355784*f[20]*dv3_sq*volFact+0.4714045207910317*f[0]*dv3_sq*volFact+0.421637021355784*f[19]*dv2_sq*volFact+0.4714045207910317*f[0]*dv2_sq*volFact+0.421637021355784*f[18]*dv1_sq*volFact+0.4714045207910317*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc2x3vSer_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/32; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[4] += 5.656854249492382*f[0]*wx3_sq*volFact+3.265986323710906*f[5]*dv3*wx3*volFact+5.656854249492382*f[0]*wx2_sq*volFact+3.265986323710906*f[4]*dv2*wx2*volFact+5.656854249492382*f[0]*wx1_sq*volFact+3.265986323710906*f[3]*dv1*wx1*volFact+0.421637021355784*f[20]*dv3_sq*volFact+0.4714045207910317*f[0]*dv3_sq*volFact+0.421637021355784*f[19]*dv2_sq*volFact+0.4714045207910317*f[0]*dv2_sq*volFact+0.421637021355784*f[18]*dv1_sq*volFact+0.4714045207910317*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc2x3vSer_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]/32; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[4] += 5.656854249492382*f[0]*wx3_sq*volFact+3.265986323710906*f[5]*dv3*wx3*volFact+5.656854249492382*f[0]*wx2_sq*volFact+3.265986323710906*f[4]*dv2*wx2*volFact+5.656854249492382*f[0]*wx1_sq*volFact+3.265986323710906*f[3]*dv1*wx1*volFact+0.421637021355784*f[20]*dv3_sq*volFact+0.4714045207910317*f[0]*dv3_sq*volFact+0.421637021355784*f[19]*dv2_sq*volFact+0.4714045207910317*f[0]*dv2_sq*volFact+0.421637021355784*f[18]*dv1_sq*volFact+0.4714045207910317*f[0]*dv1_sq*volFact; 
} 
