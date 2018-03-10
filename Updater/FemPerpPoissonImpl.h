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
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/SparseExtra>

#include <mpi.h>

class FemPerpPoisson;

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
  } bcdata_t;

// C wrappers for interfacing with FemPerpPoisson class
  void* new_FemPerpPoisson(int nx, int ny, int ndim, int polyOrder, double dx, double dy, bool periodicFlgs[2], bcdata_t bc[2][2], bool writeMatrix, double laplacianWeight, double modifierConstant);
  void delete_FemPerpPoisson(FemPerpPoisson* f);
  void createGlobalSrc(FemPerpPoisson* f, double* localSrcPtr, int idx, int idy, double intSrcVol);
  void zeroGlobalSrc(FemPerpPoisson* f);
  void allreduceGlobalSrc(FemPerpPoisson* f, MPI_Comm comm);
  void getSolution(FemPerpPoisson* f, double* localSolPtr, int idx, int idy);
  void getNodalSolution(FemPerpPoisson* f, double* localSolPtr, int idx, int idy);
}

class FemPerpPoisson
{
 public:
  FemPerpPoisson(int nx, int ny, int ndim, int polyOrder, 
             double dx, double dy, bool periodicFlgs[2],
             bcdata_t bc[2][2], bool writeMatrix,
             double laplacianWeight, double modifierConstant) ;
  ~FemPerpPoisson();
  void createGlobalSrc(double* ptr, int idx, int idy, double intSrcVol);
  void zeroGlobalSrc();
  void allreduceGlobalSrc(MPI_Comm comm);
  void solve();
  void getSolution(double* ptr, int idx, int idy);
  void getNodalSolution(double* ptr, int idx, int idy);

 private:
  const int nx, ny, ndim, polyOrder;
  const double dx, dy;
  const bool writeMatrix;
  const double laplacianWeight, modifierConstant;
  bcdata_t bc[2][2], bc2d[2][2], bc2d_z0[2][2];
  bool periodicFlgs[2];
  bool allPeriodic;
  MPI_Datatype MPI_vector_t;
  MPI_Op MPI_vectorSum_op;
  /** Eigen sparse matrix to store stiffness matrix */
  Eigen::SparseMatrix<double,Eigen::ColMajor> stiffMat, stiffMat_z0;
  /** Eigen vectors for source and dirichlet modifications to source*/
  Eigen::VectorXd globalSrc, sourceModVec, sourceModVec_z0;
  /** Eigen vector for solution */
  Eigen::VectorXd x;
  /** Eigen solver method */
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver, solver_z0;
  Eigen::MatrixXd localMassModToNod, localNodToMod, localModToNod;
  
  int getNumPerpGlobalNodes(int nx, int ny, int ndim, int p, bool periodicFlgs[2]);
  int getNumLocalNodes(int ndim, int p);
  
  void setupBoundaryIndices(bcdata_t bc[2][2], int ndim, int polyOrder);
  void makeGlobalPerpStiffnessMatrix(
     Eigen::SparseMatrix<double,Eigen::ColMajor>& stiffMat,
     Eigen::VectorXd& sourceModVec,
     int ndim, int polyOrder, bcdata_t bc[2][2]);
  void getPerpStiffnessMatrix(Eigen::MatrixXd& localStiff, int ndim, int p, double dx, double dy);
  void getNodToModMatrix(Eigen::MatrixXd& localNodToMod, int ndim, int p);
  void getMassMatrix(Eigen::MatrixXd& localMass, int ndim, int p);
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
