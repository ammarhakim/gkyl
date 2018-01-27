#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc1x1vMax_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
 
  out[0] += 2.0*f[0]*volFact; 
  out[4] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x1vMax_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
 
  out[0] += 2.0*f[0]*volFact; 
  out[4] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.149071198499986*f[5]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x1vMax_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
 
  out[0] += 2.0*f[0]*volFact; 
  out[4] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.149071198499986*f[5]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x1vMax_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
 
  out[0] += 2.0*f[0]*volFact; 
  out[4] += 2.0*f[0]*wx1_sq*volFact+1.154700538379252*f[2]*dv1*wx1*volFact+0.149071198499986*f[5]*dv1_sq*volFact+0.1666666666666667*f[0]*dv1_sq*volFact; 
} 
