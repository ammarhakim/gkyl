#ifndef GKSHEATH_MOD_DECL_H 
#define GKSHEATH_MOD_DECL_H 
#include <cmath> 
#include <algorithm>
#include <CartFieldBinOpModDecl.h>

extern "C" { 
double calcSheathDeltaPhi1xSer_P1(const double *phi, const double *phiWall, const double zVal);
void calcSheathPartialReflectionWeakEquiv1x1vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
void calcSheathPartialReflectionScaled1x1vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
void calcSheathReflection1x1vSer_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
double calcSheathDeltaPhi1xSer_P2(const double *phi, const double *phiWall, const double zVal);
void calcSheathPartialReflectionWeakEquiv1x2vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
void calcSheathPartialReflectionScaled1x2vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
void calcSheathReflection1x2vSer_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
double calcSheathDeltaPhi3xSer_P2(const double *phi, const double *phiWall, const double zVal);
void calcSheathPartialReflectionWeakEquiv3x2vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
void calcSheathPartialReflectionScaled3x2vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
void calcSheathReflection3x2vSer_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
} 
#endif 
