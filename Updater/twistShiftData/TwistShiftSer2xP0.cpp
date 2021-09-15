#include <TwistShiftModDecl.h> 
 
void twistShift_xLimDG_yShP1_2xSerP0(const double sFac, const double *xLimLo, const double *xLimUp, const double yLimLo, const double yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew) 
{ 
  // sFac:   scale factor, typically +/- 1.0.
  // xLimUp: 1D DG expansion of the function xLimLo(y) giving the lower limit of the x integral.
  // xLimLo: 1D DG expansion of the function xLimUp(y) giving the upper limit of the x integral.
  // yLimLo: lower limit of the y integral in [-1,1] logical space.
  // yLimUp: upper limit of the y integral in [-1,1] logical space.
  // dyDo:   y cell length of the donor cell.
  // yOff:   yOffset between the donor and target cells (in the direction of the shift).
  // ySh:    yShift translated to [-1,1] logical space.
  // tsData: twist-shift matrices and vectors.
  // xIdx:   x-index of the cell.
  // pushNew: push_back a new matrix (true) or add contribution to the last matrix (false).

  // Length of the subregion in which the DG expansion of the yLimLo/yLimUp functions are defined.
  double dyLim = yLimUp-yLimLo;

  // Logical center of the subregion in which the DG expansion of the yLimLo/yLimUp functions are defined.
  double ycLim = 0.5*(yLimUp+yLimLo);

  tsData->mat.setZero();
  double yLimLoR2 = std::pow(yLimLo,2);
  double yLimUpR2 = std::pow(yLimUp,2);

  tsData->mat(0,0) += -(0.125*(((4.898979485566357*xLimUp[1]-4.898979485566357*xLimLo[1])*sFac*yLimUp+(4.898979485566357*xLimLo[1]-4.898979485566357*xLimUp[1])*sFac*yLimLo)*ycLim+(2.449489742783178*xLimLo[1]-2.449489742783178*xLimUp[1])*sFac*yLimUpR2+(1.414213562373095*xLimLo[0]-1.414213562373095*xLimUp[0])*dyLim*sFac*yLimUp+(2.449489742783178*xLimUp[1]-2.449489742783178*xLimLo[1])*sFac*yLimLoR2+(1.414213562373095*xLimUp[0]-1.414213562373095*xLimLo[0])*dyLim*sFac*yLimLo))/dyLim; 

  if (pushNew) {
    tsData->cellMat[xIdx-1].push_back(tsData->mat);
  } else {
    tsData->cellMat[xIdx-1].back() += tsData->mat;
  };

}

void twistShift_yLimDG_yShP1_2xSerP0(const double sFac, const double xLimLo, const double xLimUp, const double *yLimLo, const double *yLimUp, const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew) 
{ 
  // sFac:    scale factor, typically +/- 1.0.
  // xLimUp:  lower limit of the x integral in [-1,1] logical space.
  // xLimLo:  upper limit of the x integral in [-1,1] logical space.
  // yLimUp:  1D DG expansion of the function yLimLo(x) giving the lower limit of the y integral.
  // yLimLo:  1D DG expansion of the function yLimUp(x) giving the upper limit of the y integral.
  // dyDo:    y cell length of the donor cell.
  // yOff:    yOffset between the donor and target cells (in the direction of the shift).
  // ySh:     yShift translated to [-1,1] logical space.
  // tsData:  twist-shift matrices and vectors.
  // xIdx:    x-index of the cell.
  // pushNew: push_back a new matrix (true) or add contribution to the last matrix (false).

  // Length of the subregion in which the DG expansion of the yLimLo/yLimUp functions are defined.
  double dxLim = xLimUp-xLimLo;

  // Logical center of the subregion in which the DG expansion of the yLimLo/yLimUp functions are defined.
  double xcLim = 0.5*(xLimUp+xLimLo);

  tsData->mat.setZero();
  double xLimLoR2 = std::pow(xLimLo,2);
  double xLimUpR2 = std::pow(xLimUp,2);

  tsData->mat(0,0) += -(0.125*(((4.898979485566357*yLimUp[1]-4.898979485566357*yLimLo[1])*sFac*xLimUp+(4.898979485566357*yLimLo[1]-4.898979485566357*yLimUp[1])*sFac*xLimLo)*xcLim+(2.449489742783178*yLimLo[1]-2.449489742783178*yLimUp[1])*sFac*xLimUpR2+(1.414213562373095*yLimLo[0]-1.414213562373095*yLimUp[0])*dxLim*sFac*xLimUp+(2.449489742783178*yLimUp[1]-2.449489742783178*yLimLo[1])*sFac*xLimLoR2+(1.414213562373095*yLimUp[0]-1.414213562373095*yLimLo[0])*dxLim*sFac*xLimLo))/dxLim; 

  if (pushNew) {
    tsData->cellMat[xIdx-1].push_back(tsData->mat);
  } else {
    tsData->cellMat[xIdx-1].back() += tsData->mat;
  };

}

void twistShift_fullCell_yShP1_2xSerP0(const double dyDo, const double yOff, const double *ySh, tsStruct *tsData, const int xIdx, const bool pushNew) 
{ 
  // dyDo:   y cell length of the donor cell.
  // yOff:   yOffset between the donor and target cells (in the direction of the shift).
  // ySh:    yShift translated to [-1,1] logical space.
  // tsData: twist-shift matrices and vectors.
  // xIdx:   x-index of the cell.
  // pushNew: push_back a new matrix (true) or add contribution to the last matrix (false).

  tsData->mat.setZero();

  tsData->mat(0,0) += 1.0; 

  if (pushNew) {
    tsData->cellMat[xIdx-1].push_back(tsData->mat);
  } else {
    tsData->cellMat[xIdx-1].back() += tsData->mat;
  };

}

void twistShift_matVecMult_2xSerP0(tsStruct *tsData, const int xIdx, const int matIdx, const double *fldDo, double *fldTar) 
{ 
  // tsData: twist-shift matrices and vectors.
  // xIdx:   x-index of the cell.
  // matIdx: index of the matrix to be assigned.
  // fldDo:  donor field.
  // fldTar: target field.

  tsData->vecDo << fldDo[0];

  tsData->vecTar = tsData->cellMat[xIdx-1][matIdx-1] * tsData->vecDo;

  Eigen::Map<Eigen::VectorXd>(fldTar,1,1) += tsData->vecTar;

}

