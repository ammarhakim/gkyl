// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for 2D FEM Poisson solver
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef DISCONT_POISSON_H
#define DISCONT_POISSON_H

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

class DiscontPoisson;

extern "C" {
/** Structure to store BC data. */
  typedef struct {
    /** Flag to indicate if Bc was set */
    bool isSet;
    /** Boundary condition type: one of 0 (for Dirichlet), 1 (for Neumann) */
    int type;
    /** Value to apply */
    double value;

    int istart[8];
    int iend[8];
    int cornerstart[8];
    int cornerend[8];
  } bcdata_t;

  // C wrappers for interfacing with DiscontPoisson class
  void* new_DiscontPoisson(int nx, int ny, int ndim, int polyOrder, double dx, double dy, bool periodicFlgs[2], bcdata_t bc[2][2]);
  void delete_DiscontPoisson(DiscontPoisson* f);

  void discontPoisson_pushTripletSet(DiscontPoisson* f, int idxX, int idxY);
  void discontPoisson_pushTriplet(DiscontPoisson* f, int i, int j, double val);
  void discontPoisson_constructStiffMatrix(DiscontPoisson* f);
  void discontPoisson_pushSource(DiscontPoisson* f, int idx, double* src);
  void discontPoisson_getSolution(DiscontPoisson* f, int idx, double* sol);
  void discontPoisson_solve(DiscontPoisson* f);
}

class DiscontPoisson
{
 public:
  DiscontPoisson(int nx, int ny, int ndim, int polyOrder, 
                 double dx, double dy, bool periodicFlgs[2],
                 bcdata_t bc[2][2]);
  ~DiscontPoisson();

  int getNumLocalNodes(int ndim, int p);
  unsigned indexer2Dto1D(int i, int j, int basisIdx);
  void pushTripletSet(int idxX, int idxY);
  void pushTriplet(int i, int j, double val);
  void constructStiffMatrix();
  void pushSource(int idx, double* src);
  void getSolution(int idx, double* sol);
  void solve();
  
 private:
  const int nx, ny, ndim, polyOrder;
  const double dx, dy;
  int nb, N;
  
  bcdata_t bc[2][2], bc2d[2][2], bc2d_z0[2][2];
  bool periodicFlgs[2];
  bool allPeriodic;
  
  bool _first = true;
  std::vector<int> sizes, displs;
  MPI_Datatype MPI_triplet_t;
  MPI_Op MPI_vectorSum_op;

  // Eigen triplets
  std::vector<Eigen::Triplet<double> > stiffTripletList;
  std::vector<Eigen::Triplet<double> > stiffTripletListGathered;
  
  // Eigen sparse matrix to store stiffness matrix 
  Eigen::SparseMatrix<double,Eigen::ColMajor> stiffMat;
  /** row major copy of stiffness matrix */
  Eigen::SparseMatrix<double,Eigen::RowMajor> stiffMatRowMajor;
  
  // Matrix for modifying source when there are Dirichlet BCs
  Eigen::SparseMatrix<double,Eigen::ColMajor> sourceModMat;
  Eigen::SparseMatrix<double,Eigen::ColMajor> dirichletIdentity;
  
  // Eigen vectors for source
  Eigen::VectorXd globalSrc;
  
  /** Eigen vector for solution */
  Eigen::VectorXd x;
  
  /** Eigen solver method */
  //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;

  bool analyzed_; // flag so that stiffness matrix only analyzed once
};


#endif // DISCONT_POISSON_H
