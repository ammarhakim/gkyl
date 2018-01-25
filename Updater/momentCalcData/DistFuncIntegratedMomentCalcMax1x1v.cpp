#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc1x1vMax_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
 
  out[0] += 2.0*(%e^((-vDrift^2/(2*vt^2))+(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact; 
  out[4] += 2.0*(%e^((-vDrift^2/(2*vt^2))+vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx1_sq+1.154700538379252*(%e^((-vDrift^2/(2*vt^2))+vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[2]*dv1*volFact*wx1+0.1666666666666667*(%e^((-vDrift^2/(2*vt^2))+vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x1vMax_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
 
  out[0] += 2.0*(%e^((-vDrift^2/(2*vt^2))+(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact; 
  out[4] += 2.0*(%e^((-vDrift^2/(2*vt^2))+vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx1_sq+1.154700538379252*(%e^((-vDrift^2/(2*vt^2))+vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[2]*dv1*volFact*wx1+0.149071198499986*(%e^((-vDrift^2/(2*vt^2))+vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[5]*dv1_sq*volFact+0.1666666666666667*(%e^((-vDrift^2/(2*vt^2))+vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-vDrift/vt^2-1/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*dv1_sq*volFact; 
} 
