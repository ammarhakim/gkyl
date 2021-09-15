// -- Gkyl ---------------------------------------------------------------------
//
// C++ back-end for allocating matrices and vectors in twist-shift struct.
//
//    _______     ___
// + 6 @ |||| # P ||| +
// -----------------------------------------------------------------------------

#include <TwistShiftModDecl.h>

using namespace Eigen;

tsStruct::tsStruct(const int numCells, const int matN) {
  mat = Eigen::MatrixXd::Zero(matN,matN);
  vecDo = Eigen::VectorXd::Zero(matN);
  vecTar = Eigen::VectorXd::Zero(matN);
  cellMat.resize(numCells);
}
extern "C" void* twistShift_alloc(const int numCells, const int matN)
{
  tsStruct* m = new tsStruct(numCells, matN);
  return reinterpret_cast<void*>(m);
}
extern "C" void twistShift_allocCellMat(tsStruct *tsData, const int cellIdx, const int numMats)
{
  tsData->cellMat[cellIdx-1].reserve(numMats);
}
extern "C" void twistShift_delete(tsStruct *tsData)
{
  delete tsData;
}
