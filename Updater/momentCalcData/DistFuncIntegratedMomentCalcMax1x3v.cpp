#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc1x3vMax_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]/16; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 4.0*(%e^((-vDrift^2/(2*vt^2))+(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact; 
  out[4] += 4.0*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx3_sq+2.309401076758503*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[4]*dv3*volFact*wx3+4.0*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx2_sq+2.309401076758503*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[3]*dv2*volFact*wx2+4.0*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx1_sq+2.309401076758503*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[2]*dv1*volFact*wx1+0.3333333333333333*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*dv3_sq*volFact+0.3333333333333333*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*dv2_sq*volFact+0.3333333333333333*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x3vMax_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]/16; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx3 = w[3], dv3 = dxv[3]; 
  const double wx3_sq = wx3*wx3, dv3_sq = dv3*dv3; 
 
  out[0] += 4.0*(%e^((-vDrift^2/(2*vt^2))+(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact; 
  out[4] += 4.0*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx3_sq+2.309401076758503*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[4]*dv3*volFact*wx3+4.0*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx2_sq+2.309401076758503*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[3]*dv2*volFact*wx2+4.0*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx1_sq+2.309401076758503*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[2]*dv1*volFact*wx1+0.2981423969999719*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[14]*dv3_sq*volFact+0.3333333333333333*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*dv3_sq*volFact+0.2981423969999719*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[13]*dv2_sq*volFact+0.3333333333333333*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*dv2_sq*volFact+0.2981423969999719*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[12]*dv1_sq*volFact+0.3333333333333333*(%e^((-vDrift^2/(2*vt^2))+(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(3*vDrift)/vt^2-9/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*dv1_sq*volFact; 
} 
