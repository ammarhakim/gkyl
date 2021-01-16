// -- Gkyl ---------------------------------------------------------------------
//
// C++ header for TwistShift kernels.
//
//    _______     ___
// + 6 @ |||| # P ||| +
// -----------------------------------------------------------------------------
#ifndef TWIST_SHIFT_MOD_DECL_H 
#define TWIST_SHIFT_MOD_DECL_H 

#include <cmath>
#include <algorithm>

#include <Eigen/Dense>
#include <vector>

struct tsStruct;

extern "C" { 

  void* twistShift_alloc(const int numMats, const int matN);
  void twistShift_allocCellMat(tsStruct *tsData, const int cellIdx, const int numMats);
  void twistShift_delete(tsStruct *tsData);

  void twistShift_matVecMult2xSerP1(tsStruct *tsData, const int xIdx, const int matIdx, const double *fldDo, double *fldTar);
  void twistShift_xLimDG2xSerP1_yShP1(const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const int matIdx);
  void twistShift_yLimDG2xSerP1_yShP1(const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const int matIdx);
  void twistShift_M2Corners2xSerP1_yShP1(const double *xLimLoL, const double *xLimUpL, const double yLimLoL, const double yLimUpL, const double *xLimLoR, const double *xLimUpR, const double yLimLoR, const double yLimUpR, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const int matIdx);


  void twistShift_matVecMult2xSerP2(tsStruct *tsData, const int xIdx, const int matIdx, const double *fldDo, double *fldTar);
  void twistShift_xLimDG2xSerP2_yShP2(const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const int matIdx);
  void twistShift_yLimDG2xSerP2_yShP2(const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const int matIdx);
  void twistShift_M2Corners2xSerP2_yShP2(const double *xLimLoL, const double *xLimUpL, const double yLimLoL, const double yLimUpL, const double *xLimLoR, const double *xLimUpR, const double yLimLoR, const double yLimUpR, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const int matIdx);


}

struct tsStruct
{
  tsStruct(const int numMats, const int matN);
  std::vector<std::vector<Eigen::MatrixXd>> cellMat;
  Eigen::MatrixXd mat;
  Eigen::VectorXd vecDo;
  Eigen::VectorXd vecTar;
};

#endif 
