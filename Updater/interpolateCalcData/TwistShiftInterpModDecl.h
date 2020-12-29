#ifndef TWIST_SHIFT_INTERP_MOD_DECL_H 
#define TWIST_SHIFT_INTERP_MOD_DECL_H 
 
#include <cmath>
#include <algorithm>
 
extern "C" { 
 
void TwistShiftInterp_xLimDG2xSer_P1(const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyLim, const double ycLim, const double dyDo, const double yOff, const double *ySh, const double *fldSrc, double *fldDest);
void TwistShiftInterp_yLimDG2xSer_P1(const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dxLim, const double xcLim, const double dyDo, const double yOff, const double *ySh, const double *fldSrc, double *fldDest);
void TwistShiftInterp_mTwoCorners2xSer_P1(const double *xLimLoL, const double *xLimUpL, const double yLimLoL, const double yLimUpL, const double dyLimL, const double ycLimL, const double *xLimLoR, const double *xLimUpR, const double yLimLoR, const double yLimUpR, const double dyLimR, const double ycLimR, const double dyDo, const double yOff, const double *ySh, const double *fldSrc, double *fldDest);


} 
#endif 
