#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc1x1vSer_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]/2; 
  out[0] += 2.0*f[0]*volFact; 
} 
void IntMomentCalc1x1vSer_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]/2; 
  out[0] += 2.0*f[0]*volFact; 
} 
