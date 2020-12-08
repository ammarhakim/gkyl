#ifndef TWIST_SHIFT_INTERP_MOD_DECL_H 
#define TWIST_SHIFT_INTERP_MOD_DECL_H 
 
#include <cmath>
#include <algorithm>
 
extern "C" { 
 
void TwistShiftInterp2xSer_P1(const double *xLimLo, const double *xLimUp, const double *fldSrc, double *fldDest);
void TwistShiftInterp_limXvarYfixed2xSer_P1(const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double *fldSrc, double *fldDest);
void TwistShiftInterp_xLimDG2xSer_P1(const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyLim, const double ycLim, const double *fldSrc, double *fldDest);


} 
#endif 
