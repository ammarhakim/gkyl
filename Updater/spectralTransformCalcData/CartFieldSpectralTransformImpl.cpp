#include <CartFieldSpectralTransformImpl.h>
 
spectralTransform::spectralTransform(const int nModes, const int nCells, const int pOrder, const int nSurfB) 
{ 
  // nModes: number of spectral modes represented. 
  // nSurfB: number of surface (un-transformed dimensions) basis elements. 
 
  // Declare stiffness matrix A as an Eigen Matrix (emA) containing the 
  // inner product of the DG bases with the global spectral bases. 
  // We will also store the inverse of this matrix in-place. 
  emA = Eigen::MatrixXd::Zero(nCells*(pOrder+1),nModes); 
  // Store the decomposition of matrix emA in-place. Need the object emAlu.
  // NOTE: various methods exist for doing the decomposition and
  //       more testing is need to understand their impact and best choice.
  //Eigen::FullPivHouseholderQR<Eigen::Ref<Eigen::MatrixXd> > emAlu(emA);
  emAlu = Eigen::FullPivHouseholderQR<Eigen::MatrixXd>(); // emAlu(emA);
  // Declare right-hand side vector B as an Eigen Vector (evB) 
  // which contains the coefficients of the DG expansion. 
  // Create one vector for each surface basis element. 
  emB = Eigen::MatrixXd::Zero(nCells*(pOrder+1), nSurfB); 
  // Declare Eigen Vector 'U' (evU) with solution to system of equations. 
  emU = Eigen::MatrixXd::Zero(nModes, nSurfB); 
}
 
void spectralTransform::assignMassMatrixSer(const int pOrder, const int cellIdx, const int spectralIdx, const double *spectralBasisIn) 
{ 
  // pOrder:          polynomial order. 
  // cellIdx:         index of current cell (1-indexed because 0 is a ghost cell). 
  // spectralIdx:     index of current spectral basis elements (0-indexed). 
  // spectralBasisIn: DG coefficients of spectral basis elements. 
 
  // Element in row (p+1)i+k and column m corresponds to the k-th coefficient in the i-th cell of the projection of the m-th spectral basis. 
  if (pOrder == 1) { 
    emA(2*(cellIdx-1),spectralIdx) = spectralBasisIn[0]; 
    emA(2*(cellIdx-1)+1,spectralIdx) = spectralBasisIn[1]; 
  } else if (pOrder == 2) { 
    emA(3*(cellIdx-1),spectralIdx) = spectralBasisIn[0]; 
    emA(3*(cellIdx-1)+1,spectralIdx) = spectralBasisIn[1]; 
    emA(3*(cellIdx-1)+2,spectralIdx) = spectralBasisIn[2]; 
  } 
 
} 
 
void spectralTransform::getMassMatrixInverse() 
{ 
 
  //emA = emA.inverse();
  //emA = emA.ldlt();
  emAlu.compute(emA);
  //Eigen::FullPivHouseholderQR<Eigen::Ref<Eigen::MatrixXd> > emAlu(emA);
 
} 
 
void spectralTransform::assignSourceMatrix1x1vSer_P1OpDir1(const int cellIdx, const double *fDG) 
{ 
  // cellIdx: index of current cell (1-indexed because 0 is a ghost cell).
  // fDG:     DG coefficients of function we wish to transform.
 
  // Element in row (p+1)i+k and column n corresponds to the k-th coefficient in the i-th cell of the projection of the DG function onto the n-th surface basis. 
  emB(2*(cellIdx-1),0) = fDG[0]; 
  emB(2*(cellIdx-1)+1,0) = fDG[1]; 
  emB(2*(cellIdx-1),1) = fDG[2]; 
  emB(2*(cellIdx-1)+1,1) = fDG[3]; 
} 
 
void spectralTransform::assignSourceMatrix1x1vSer_P1OpDir2(const int cellIdx, const double *fDG) 
{ 
  // cellIdx: index of current cell (1-indexed because 0 is a ghost cell).
  // fDG:     DG coefficients of function we wish to transform.
 
  // Element in row (p+1)i+k and column n corresponds to the k-th coefficient in the i-th cell of the projection of the DG function onto the n-th surface basis. 
  emB(2*(cellIdx-1),0) = fDG[0]; 
  emB(2*(cellIdx-1)+1,0) = fDG[2]; 
  emB(2*(cellIdx-1),1) = fDG[1]; 
  emB(2*(cellIdx-1)+1,1) = fDG[3]; 
} 
 
void spectralTransform::assignSourceMatrix1x1vSer_P2OpDir1(const int cellIdx, const double *fDG) 
{ 
  // cellIdx: index of current cell (1-indexed because 0 is a ghost cell).
  // fDG:     DG coefficients of function we wish to transform.
 
  // Element in row (p+1)i+k and column n corresponds to the k-th coefficient in the i-th cell of the projection of the DG function onto the n-th surface basis. 
  emB(3*(cellIdx-1),0) = fDG[0]; 
  emB(3*(cellIdx-1)+1,0) = fDG[1]; 
  emB(3*(cellIdx-1)+2,0) = fDG[4]; 
  emB(3*(cellIdx-1),1) = fDG[2]; 
  emB(3*(cellIdx-1)+1,1) = fDG[3]; 
  emB(3*(cellIdx-1)+2,1) = 1.0*fDG[6]; 
  emB(3*(cellIdx-1),2) = fDG[5]; 
  emB(3*(cellIdx-1)+1,2) = 1.0*fDG[7]; 
} 
 
void spectralTransform::assignSourceMatrix1x1vSer_P2OpDir2(const int cellIdx, const double *fDG) 
{ 
  // cellIdx: index of current cell (1-indexed because 0 is a ghost cell).
  // fDG:     DG coefficients of function we wish to transform.
 
  // Element in row (p+1)i+k and column n corresponds to the k-th coefficient in the i-th cell of the projection of the DG function onto the n-th surface basis. 
  emB(3*(cellIdx-1),0) = fDG[0]; 
  emB(3*(cellIdx-1)+1,0) = fDG[2]; 
  emB(3*(cellIdx-1)+2,0) = fDG[5]; 
  emB(3*(cellIdx-1),1) = fDG[1]; 
  emB(3*(cellIdx-1)+1,1) = fDG[3]; 
  emB(3*(cellIdx-1)+2,1) = 1.0*fDG[7]; 
  emB(3*(cellIdx-1),2) = fDG[4]; 
  emB(3*(cellIdx-1)+1,2) = 1.0*fDG[6]; 
} 
 
 
void spectralTransform::solveLinearProblem() 
{ 
 
  //emU = emA.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(emB);
  //emU = emA.fullPivHouseholderQr().solve(emB);
  emU = emAlu.solve(emB);
 
} 
 
void spectralTransform::redistributeSolution1x1vSer_P1(const int cellIdx, double *fSpectral) 
{ 
  // cellIdx:   index of current cell (1-indexed because 0 is a ghost cell).
  // fSpectral: spectral coefficients of transformed function.
 
  Eigen::Map<Eigen::MatrixXd>(fSpectral,2,2) = emU.block(2*(cellIdx-1),0,2,2); 
 
} 
void spectralTransform::redistributeSolution1x1vSer_P2(const int cellIdx, double *fSpectral) 
{ 
  // cellIdx:   index of current cell (1-indexed because 0 is a ghost cell).
  // fSpectral: spectral coefficients of transformed function.
 
  Eigen::Map<Eigen::MatrixXd>(fSpectral,3,3) = emU.block(3*(cellIdx-1),0,3,3); 
 
} 
 
// C wrappers for interfacing with spectralTransform class. 
extern "C" void* new_spectralTransform(int nModes, const int nCells, const int pOrder, int nSurfB) 
{ 
  spectralTransform* b = new spectralTransform(nModes, nCells, pOrder, nSurfB); 
  return reinterpret_cast<void*>(b); 
} 
 
extern "C" void assignLHSMatrixSer(spectralTransform *sTransObj, const int pOrder, const int cellIdx, const int spectralIdx, const double *spectralBasisIn) 
{ 
  sTransObj->assignMassMatrixSer(pOrder, cellIdx, spectralIdx, spectralBasisIn); 
} 
 
extern "C" void getLHSMatrixInverse(spectralTransform *sTransObj) 
{ 
  sTransObj->getMassMatrixInverse(); 
} 
 
extern "C" void assignRHSMatrix1x1vSer_P1OpDir1(spectralTransform *sTransObj, const int cellIdx, const double *fDG) 
{ 
  sTransObj->assignSourceMatrix1x1vSer_P1OpDir1(cellIdx, fDG); 
} 
 
extern "C" void assignRHSMatrix1x1vSer_P1OpDir2(spectralTransform *sTransObj, const int cellIdx, const double *fDG) 
{ 
  sTransObj->assignSourceMatrix1x1vSer_P1OpDir2(cellIdx, fDG); 
} 
 
extern "C" void assignRHSMatrix1x1vSer_P2OpDir1(spectralTransform *sTransObj, const int cellIdx, const double *fDG) 
{ 
  sTransObj->assignSourceMatrix1x1vSer_P2OpDir1(cellIdx, fDG); 
} 
 
extern "C" void assignRHSMatrix1x1vSer_P2OpDir2(spectralTransform *sTransObj, const int cellIdx, const double *fDG) 
{ 
  sTransObj->assignSourceMatrix1x1vSer_P2OpDir2(cellIdx, fDG); 
} 
 
 
extern "C" void solveTransform(spectralTransform *sTransObj) 
{ 
  sTransObj->solveLinearProblem(); 
} 
 
extern "C" void getSolution1x1vSer_P1(spectralTransform *sTransObj, const int cellIdx, double *fSpectral) 
{ 
  sTransObj->redistributeSolution1x1vSer_P1(cellIdx, fSpectral); 
} 
 
extern "C" void getSolution1x1vSer_P2(spectralTransform *sTransObj, const int cellIdx, double *fSpectral) 
{ 
  sTransObj->redistributeSolution1x1vSer_P2(cellIdx, fSpectral); 
} 
 
 
