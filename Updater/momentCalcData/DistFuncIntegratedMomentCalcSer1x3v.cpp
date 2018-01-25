#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc1x3vSer_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 4.0*f[0]*volFact; 
} 
void IntMomentCalc1x3vSer_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]/8; 
  out[0] += 4.0*f[0]*volFact; 
} 
