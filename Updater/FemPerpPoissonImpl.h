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
    int cornerstart[3];
    int cornerend[3];
  } bcdata_t;

// C wrappers for interfacing with FemPerpPoisson class
  void* new_FemPerpPoisson(int nx, int ny, int ndim, int polyOrder, double dx, double dy, bool periodicFlgs[2], bcdata_t bc[2][2], bool writeMatrix, bool adjustSource);
  void delete_FemPerpPoisson(FemPerpPoisson* f);
  void makeGlobalStiff(FemPerpPoisson* f, double *laplacianWeight, double *modifierWeight, double *gxx, double *gxy, double *gyy, int idx, int idy);
  void finishGlobalStiff(FemPerpPoisson* f);
  void createGlobalSrc(FemPerpPoisson* f, double* localSrcPtr, int idx, int idy, double intSrcVol);
  void zeroGlobalSrc(FemPerpPoisson* f);
  void allreduceGlobalSrc(FemPerpPoisson* f, MPI_Comm comm);
  void allgatherGlobalStiff(FemPerpPoisson* f, MPI_Comm comm);
  void getSolution(FemPerpPoisson* f, double* localSolPtr, int idx, int idy);
  void getNodalSolution(FemPerpPoisson* f, double* localSolPtr, int idx, int idy);
}

class FemPerpPoisson
{
 public:
  FemPerpPoisson(int nx, int ny, int ndim, int polyOrder, 
             double dx, double dy, bool periodicFlgs[2],
             bcdata_t bc[2][2], bool writeMatrix,
             bool adjustSource);
  ~FemPerpPoisson();
  void createGlobalSrc(double* ptr, int idx, int idy, double intSrcVol);
  void zeroGlobalSrc();
  void allreduceGlobalSrc(MPI_Comm comm);
  void allgatherGlobalStiff(MPI_Comm comm);
  void makeGlobalPerpStiffnessMatrix(double *laplacianWeight, double *modifierWeight, double *gxx, double *gxy, double *gyy, int idx, int idy);
  void finishGlobalPerpStiffnessMatrix();
  void solve();
  void getSolution(double* ptr, int idx, int idy);
  void getNodalSolution(double* ptr, int idx, int idy);

 private:
  const int nx, ny, ndim, polyOrder;
  const double dx, dy;
  const bool writeMatrix;
  bcdata_t bc[2][2], bc2d[2][2], bc2d_z0[2][2];
  bool periodicFlgs[2];
  bool allPeriodic;
  const bool adjustSource;
  double cornerval;
  MPI_Datatype MPI_triplet_t;
  MPI_Op MPI_vectorSum_op;
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
  Eigen::MatrixXd localMassModToNod, localNodToMod, localModToNod;
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
