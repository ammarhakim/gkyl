// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for 2D FEM Poisson solver
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef FEM_POISSON_H
#define FEM_POISSON_H

#include <stdint.h>

// eigen includes
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/SparseExtra>

#include <mpi.h>

class FemGyroaverage;

extern "C" {
/** Structure to store BC data. */
  typedef struct 
  {
/** Flag to indicate if Bc was set */
    bool isSet;
/** Boundary condition type: one of 0 (for Dirichlet), 1 (for Neumann) */
    int type;
/** Value to apply */
    double value;

    int istart[3];
    int iend[3];
    int cornerstart[3];
    int cornerend[3];
  } bcdata_t;

// C wrappers for interfacing with FemGyroaverage class
  void* new_FemGyroaverage(int nx, int ny, int ndim, int polyOrder, bool periodicFlgs[2], bcdata_t bc[2][2], bool writeMatrix);
  void delete_FemGyroaverage(FemGyroaverage* f);
  void makeGlobalStiffGy(FemGyroaverage* f, double *modifierWeight, int idx, int idy);
  void makeGyavgMatrix(FemGyroaverage* f, double *rho1, double *rho2, double *rho3, int idx);
  void finishGlobalStiffGy(FemGyroaverage* f);
  void createGlobalSrcGy(FemGyroaverage* f, double* localSrcPtr, int idx, int idy);
  void zeroGlobalSrcGy(FemGyroaverage* f);
  void allreduceGlobalSrcGy(FemGyroaverage* f, MPI_Comm comm);
  void allgatherGlobalStiffGy(FemGyroaverage* f, MPI_Comm comm);
  void getSolutionGy(FemGyroaverage* f, double* localSolPtr, int idx, int idy);
  void getNodalSolutionGy(FemGyroaverage* f, double* localSolPtr, int idx, int idy);
  void solveGy(FemGyroaverage* f);
}

class FemGyroaverage
{
 public:
  FemGyroaverage(int nx, int ny, int ndim, int polyOrder, 
             bool periodicFlgs[2], bcdata_t bc[2][2], bool writeMatrix);
  ~FemGyroaverage();
  void createGlobalSrcGy(double* ptr, int idx, int idy);
  void zeroGlobalSrcGy();
  void allreduceGlobalSrcGy(MPI_Comm comm);
  void allgatherGlobalStiffGy(MPI_Comm comm);
  void makeGlobalPerpStiffnessMatrix(double *modifierWeight, int idx, int idy);
  void makeGyavgMatrix(double *rho1, double *rho2, double *rho3, int idx);
  void finishGlobalPerpStiffnessMatrix();
  void solveGy();
  void getSolutionGy(double* ptr, int idx, int idy);
  void getNodalSolutionGy(double* ptr, int idx, int idy);

 private:
  const int nx, ny, ndim, polyOrder;
  const bool writeMatrix;
  bcdata_t bc[2][2];
  bool periodicFlgs[2];
  bool allPeriodic;
  double cornerval;
  MPI_Datatype MPI_triplet_t;
  MPI_Op MPI_vectorSumGy_op;
  std::vector<Eigen::Triplet<double> > stiffTripletList;
  /** Eigen sparse matrix to store stiffness matrix */
  Eigen::SparseMatrix<double,Eigen::ColMajor> stiffMat;
  /** Eigen vectors for source and dirichlet modifications to source*/
  Eigen::VectorXd globalSrc, sourceModVec;
  /** Eigen vector for solution */
  Eigen::VectorXd x;
  /** Eigen solver method */
  //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  //Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
  Eigen::MatrixXd *localGyavgModToNod;
  Eigen::MatrixXd localNodToMod, localModToNod;
  bool analyzed_; // flag so that stiffness matrix only analyzed once
  
  int getNumPerpGlobalNodes(int nx, int ny, int ndim, int p, bool periodicFlgs[2]);
  int getNumLocalNodes(int ndim, int p);
  
  void setupBoundaryIndices(bcdata_t bc[2][2], int ndim, int polyOrder);
  void getPerpLocalToGlobalInteriorBLRT(std::vector<int>& lgMap, int idx, int idy, int nx, int ny, int ndim, int p, bool periodicFlgs[2]);
  
  /**
    * For local-to-global interior-boundary mapping, p=1
    */
  int F1_func(int nx, int ny, int ninterior, int ix, int iy, bool periodicFlgs[2]);
  /**
    * For local-to-global interior-boundary mapping, p=2 
    */
  int F2_func(int nx, int ny, int ninterior, int ix, int iy, bool periodicFlgs[2]);
  /**
    * For local-to-global interior-boundary mapping, p=2 
    */
  int G2_func(int nx, int ny, int ninterior, int ix, int iy, bool periodicFlgs[2]);
};


#endif // FEM_POISSON_H
