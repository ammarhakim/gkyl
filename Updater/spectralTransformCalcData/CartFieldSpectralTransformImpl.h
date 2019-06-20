#ifndef CART_FIELD_SPECTRAL_TRANSFORM_H
#define CART_FIELD_SPECTRAL_TRANSFORM_H
 
// Eigen include statements.
#include <Eigen/Dense>
 
class spectralTransform;
 
extern "C" {
 
  // C wrappers that interface with spectralTransform class.
  void* new_spectralTransform(const int nModes, const int nSurfB);
  // Assign elements (coefficients of projection of spectral basis)
  // of the left-side mass matrix.
  void assignLHSMatrixSer(spectralTransform *sTransObj, const int pOrder, const int cellIdx, const int spectralIdx, const double *spectralBasisIn);
  // Obtain the inverse of the mass matrix.
  void getLHSMatrixInverse(spectralTransform *sTransObj);
  // The right-side (source) vector is actually a matrix: one vector
  // for each element of the surface (un-transformed) basis.
  void assignRHSMatrix1x1vSer_P1OpDir1(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
  void assignRHSMatrix1x1vSer_P1OpDir2(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
  void assignRHSMatrix1x1vSer_P2OpDir1(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
  void assignRHSMatrix1x1vSer_P2OpDir2(spectralTransform *sTransObj, const int cellIdx, const double *fDG);
  // Solve the spectral transform linear problem.
  void solveTransform(spectralTransform *sTransObj);
  // Take solution from matrix and put it in a Gkeyll CartField.
  void getSolution1x1vSer_P1(spectralTransform *sTransObj, const int cellIdx, double *fSpectral);
  void getSolution1x1vSer_P2(spectralTransform *sTransObj, const int cellIdx, double *fSpectral);
 
}
class spectralTransform
{
  public:
  spectralTransform(const int nModes, const int nSurfB);
  ~spectralTransform();
  
  void assignMassMatrixSer(const int pOrder, const int cellIdx, const int spectralIdx, const double *spectralBasisIn);
  void getMassMatrixInverse();
  void assignSourceMatrix1x1vSer_P1OpDir1(const int cellIdx, const double *fDG);
  void assignSourceMatrix1x1vSer_P1OpDir2(const int cellIdx, const double *fDG);
  void assignSourceMatrix1x1vSer_P2OpDir1(const int cellIdx, const double *fDG);
  void assignSourceMatrix1x1vSer_P2OpDir2(const int cellIdx, const double *fDG);
  void solveLinearProblem();
  void redistributeSolution1x1vSer_P1(const int cellIdx, double *fSpectral);
  void redistributeSolution1x1vSer_P2(const int cellIdx, double *fSpectral);

  private:
  Eigen::MatrixXd emA;
  Eigen::MatrixXd emB;
  Eigen::MatrixXd emU;
};

#endif
