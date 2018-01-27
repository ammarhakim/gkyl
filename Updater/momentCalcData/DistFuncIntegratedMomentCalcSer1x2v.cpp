#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc1x2vSer_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[4] += 2.828427124746191*f[0]*wx2_sq*volFact+1.632993161855453*f[3]*dv2*wx2*volFact+2.828427124746191*f[0]*wx1_sq*volFact+1.632993161855453*f[2]*dv1*wx1*volFact+0.2357022603955158*f[0]*dv2_sq*volFact+0.2357022603955158*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x2vSer_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[4] += 2.828427124746191*f[0]*wx2_sq*volFact+1.632993161855453*f[3]*dv2*wx2*volFact+2.828427124746191*f[0]*wx1_sq*volFact+1.632993161855453*f[2]*dv1*wx1*volFact+0.210818510677892*f[9]*dv2_sq*volFact+0.2357022603955158*f[0]*dv2_sq*volFact+0.210818510677892*f[8]*dv1_sq*volFact+0.2357022603955158*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x2vSer_P3(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[4] += 2.828427124746191*f[0]*wx2_sq*volFact+1.632993161855453*f[3]*dv2*wx2*volFact+2.828427124746191*f[0]*wx1_sq*volFact+1.632993161855453*f[2]*dv1*wx1*volFact+0.210818510677892*f[9]*dv2_sq*volFact+0.2357022603955158*f[0]*dv2_sq*volFact+0.210818510677892*f[8]*dv1_sq*volFact+0.2357022603955158*f[0]*dv1_sq*volFact; 
} 
void IntMomentCalc1x2vSer_P4(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]/8; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[4] += 2.828427124746191*f[0]*wx2_sq*volFact+1.632993161855453*f[3]*dv2*wx2*volFact+2.828427124746191*f[0]*wx1_sq*volFact+1.632993161855453*f[2]*dv1*wx1*volFact+0.210818510677892*f[9]*dv2_sq*volFact+0.2357022603955158*f[0]*dv2_sq*volFact+0.210818510677892*f[8]*dv1_sq*volFact+0.2357022603955158*f[0]*dv1_sq*volFact; 
} 
