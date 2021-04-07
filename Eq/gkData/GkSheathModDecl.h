#ifndef GKSHEATH_MOD_DECL_H 
#define GKSHEATH_MOD_DECL_H 

#include <cmath>
#include <algorithm>
#include <CartFieldBinOpModDecl.h>

// approximation for inverse Langevin function 
template <typename T> T invL(T x) {
  // from Kroger 
  return (3.*x-x*x*x*(6. + x*x - 2.*x*x*x*x)/5.)/(1.-x*x); 
}

extern "C" { 

  double calcSheathDeltaPhi1xSer_P1(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv1x1vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled1x1vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection1x1vSer_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi1xSer_P1(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv1x2vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled1x2vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection1x2vSer_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi3xSer_P1(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv3x2vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled3x2vSer_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection3x2vSer_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi1xSer_P2(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv1x1vSer_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled1x1vSer_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection1x1vSer_P2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi1xSer_P2(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv1x2vSer_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled1x2vSer_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection1x2vSer_P2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi3xSer_P2(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv3x2vSer_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled3x2vSer_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection3x2vSer_P2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi1xTensor_P1(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv1x1vTensor_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled1x1vTensor_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection1x1vTensor_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi1xTensor_P1(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv1x2vTensor_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled1x2vTensor_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection1x2vTensor_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi3xTensor_P1(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv3x2vTensor_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled3x2vTensor_P1(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection3x2vTensor_P1(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi1xTensor_P2(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv1x1vTensor_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled1x1vTensor_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection1x1vTensor_P2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi1xTensor_P2(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv1x2vTensor_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled1x2vTensor_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection1x2vTensor_P2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 
  double calcSheathDeltaPhi3xTensor_P2(const double *phi, const double *phiWall, const double zVal);
  void calcSheathPartialReflectionWeakEquiv3x2vTensor_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathPartialReflectionScaled3x2vTensor_P2(binOpData_t* data, const double wv, const double dv, const double zVal, const double vcut, const double *f, double *fhat);
  void calcSheathReflection3x2vTensor_P2(const double wv, const double dv, const double vlowerSq, const double vupperSq, const double zVal, const double q_, const double m_, const double *phi, const double *phiWall, const double *f, double *fRefl); 

} 
#endif 
