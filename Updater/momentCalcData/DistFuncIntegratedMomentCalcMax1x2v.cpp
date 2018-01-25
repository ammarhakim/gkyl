#include <DistFuncMomentCalcModDecl.h> 
void IntMomentCalc1x2vMax_P1(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]/4; 
  out[0] += 2.828427124746191*f[0]*volFact; 
} 
void IntMomentCalc1x2vMax_P2(const double *w, const double *dxv, const double *f, double *out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]/4; 
  out[0] += 2.828427124746191*f[0]*volFact; 
} 
