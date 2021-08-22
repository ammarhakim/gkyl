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

struct tsStruct
{
  tsStruct(const int numMats, const int matN);
  std::vector<std::vector<Eigen::MatrixXd>> cellMat;
  Eigen::MatrixXd mat;
  Eigen::VectorXd vecDo;
  Eigen::VectorXd vecTar;
};

extern "C" { 

  void* twistShift_alloc(const int numMats, const int matN);
  void twistShift_allocCellMat(tsStruct *tsData, const int cellIdx, const int numMats);
  void twistShift_delete(tsStruct *tsData);

  void twistShift_matVecMult_2xSerP1(tsStruct *tsData, const int xIdx, const int matIdx, const double *fldDo, double *fldTar);
  void twistShift_xLimDG_yShP1_2xSerP1(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP1_2xSerP1(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP1_2xSerP1(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_xLimDG_yShP2_2xSerP1(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP2_2xSerP1(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP2_2xSerP1(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_matVecMult_2xSerP2(tsStruct *tsData, const int xIdx, const int matIdx, const double *fldDo, double *fldTar);
  void twistShift_xLimDG_yShP1_2xSerP2(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP1_2xSerP2(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP1_2xSerP2(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_xLimDG_yShP2_2xSerP2(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP2_2xSerP2(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP2_2xSerP2(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_xLimDG_yShP3_2xSerP2(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP3_2xSerP2(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP3_2xSerP2(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);


  void twistShift_matVecMult_3xSerP1(tsStruct *tsData, const int xIdx, const int matIdx, const double *fldDo, double *fldTar);
  void twistShift_xLimDG_yShP1_3xSerP1(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP1_3xSerP1(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP1_3xSerP1(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_xLimDG_yShP2_3xSerP1(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP2_3xSerP1(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP2_3xSerP1(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_matVecMult_3xSerP2(tsStruct *tsData, const int xIdx, const int matIdx, const double *fldDo, double *fldTar);
  void twistShift_xLimDG_yShP1_3xSerP2(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP1_3xSerP2(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP1_3xSerP2(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_xLimDG_yShP2_3xSerP2(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP2_3xSerP2(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP2_3xSerP2(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_xLimDG_yShP3_3xSerP2(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP3_3xSerP2(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP3_3xSerP2(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);


  void twistShift_matVecMult_3x2vSerP1(tsStruct *tsData, const int xIdx, const int matIdx, const double *fldDo, double *fldTar);
  void twistShift_xLimDG_yShP1_3x2vSerP1(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP1_3x2vSerP1(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP1_3x2vSerP1(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_xLimDG_yShP2_3x2vSerP1(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP2_3x2vSerP1(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP2_3x2vSerP1(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_matVecMult_3x2vSerP2(tsStruct *tsData, const int xIdx, const int matIdx, const double *fldDo, double *fldTar);
  void twistShift_xLimDG_yShP1_3x2vSerP2(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP1_3x2vSerP2(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP1_3x2vSerP2(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_xLimDG_yShP2_3x2vSerP2(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP2_3x2vSerP2(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP2_3x2vSerP2(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);

  void twistShift_xLimDG_yShP3_3x2vSerP2(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_yLimDG_yShP3_3x2vSerP2(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);
  void twistShift_fullCell_yShP3_3x2vSerP2(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew);


}

#endif 
