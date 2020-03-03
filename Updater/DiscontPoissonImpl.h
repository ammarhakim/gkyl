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
  // C wrappers for interfacing with DiscontPoisson class
  void* new_DiscontPoisson(int ncell[3], int ndim, int nbasis,
                           int nnonzero, int polyOrder, bool writeMatrix);
  void delete_DiscontPoisson(DiscontPoisson* f);
  void discontPoisson_pushTriplet(DiscontPoisson* f, int i, int j, double val);
  void discontPoisson_constructStiffMatrix(DiscontPoisson* f);
  void discontPoisson_pushSource(DiscontPoisson* f, int idx, double* src, double* srcMod);
  void discontPoisson_getSolution(DiscontPoisson* f, int idx, double* sol);
  void discontPoisson_solve(DiscontPoisson* f);
}

class DiscontPoisson
{
 public:
  DiscontPoisson(int ncell[3], int ndim, int nbasis,
                 int nnonzero, int polyOrder, bool writeMatrix);
  ~DiscontPoisson();

  int getNumLocalNodes(int ndim, int p);
  void pushTriplet(int i, int j, double val);
  void constructStiffMatrix();
  void pushSource(int idx, double* src, double* srcMod);
  void getSolution(int idx, double* sol);
  void solve();
  
 private:
  const int ndim, polyOrder, nbasis, nnonzero;
  const bool writeMatrix;
  int ncell[3];
  int N;
  
  // Eigen triplets
  std::vector<Eigen::Triplet<double> > stiffTripletList;
  //std::vector<Eigen::Triplet<double> > stiffTripletListGathered;
  
  // Eigen sparse matrix to store stiffness matrix 
  Eigen::SparseMatrix<double,Eigen::ColMajor> stiffMat;
  /** row major copy of stiffness matrix */
  Eigen::SparseMatrix<double,Eigen::RowMajor> stiffMatRowMajor;
  
  // Eigen vectors for source
  Eigen::VectorXd globalSrc;
  
  /** Eigen vector for solution */
  Eigen::VectorXd x;
  
  /** Eigen solver method */
  //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver;
  //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
};


#endif // DISCONT_POISSON_H
