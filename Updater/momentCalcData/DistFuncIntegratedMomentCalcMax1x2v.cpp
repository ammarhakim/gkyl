#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc1x2vMax_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  out[0] += 2.828427124746191*(%e^((-vDrift^2/(2*vt^2))+(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact; 
  out[4] += 2.828427124746191*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx2_sq+1.632993161855453*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[3]*dv2*volFact*wx2+2.828427124746191*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx1_sq+1.632993161855453*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[2]*dv1*volFact*wx1+0.2357022603955158*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[0]*dv2_sq*volFact+0.2357022603955158*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x2vMax_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  out[0] += 2.828427124746191*(%e^((-vDrift^2/(2*vt^2))+(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(v*vDrift)/vt^2-v^2/(2*vt^2))/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact; 
  out[4] += 2.828427124746191*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx2_sq+1.632993161855453*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[3]*dv2*volFact*wx2+2.828427124746191*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[0]*volFact*wx1_sq+1.632993161855453*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[2]*dv1*volFact*wx1+0.210818510677892*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[9]*dv2_sq*volFact+0.2357022603955158*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[0]*dv2_sq*volFact+0.210818510677892*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[8]*dv1_sq*volFact+0.2357022603955158*(%e^((-vDrift^2/(2*vt^2))+(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt)+%e^((-vDrift^2/(2*vt^2))-(2*vDrift)/vt^2-2/vt^2)/(2^(3/2)*sqrt(%pi)*vt))[0]*dv1_sq*volFact; 
} 
