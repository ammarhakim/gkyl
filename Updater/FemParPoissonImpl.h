// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for 2D FEM Poisson solver
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef FEM_PAR_POISSON_H
#define FEM_PAR_POISSON_H

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

class FemParPoisson;

extern "C" {
/** Structure to store BC data. */
  typedef struct 
  {
/** Flag to indicate if Bc was set */
    bool isSet;
/** Boundary condition type: one of 0 (for Dirichlet), 1 (for Neumann) */
    unsigned type;
/** Value to apply */
    double value;

    int istart;
    int iend;
  } bcdataPar_t;

// C wrappers for interfacing with FemParPoisson class
  void* new_FemParPoisson(int nz, int ndim, int polyOrder, double dz, bool periodicFlg, bcdataPar_t bc[2], bool writeMatrix);
  void delete_FemParPoisson(FemParPoisson* f);
  void makeParGlobalStiff(FemParPoisson* f, double *laplacianWeight, double *modifierWeight, int idz);
  void finishParGlobalStiff(FemParPoisson* f);
  void createParGlobalSrc(FemParPoisson* f, double* localSrcPtr, int idz, double intSrcVol);
  void zeroParGlobalSrc(FemParPoisson* f);
  void allreduceParGlobalSrc(FemParPoisson* f, MPI_Comm comm);
  void allgatherParGlobalStiff(FemParPoisson* f, MPI_Comm comm);
  void getSolutionPar(FemParPoisson* f, double* localSolPtr, int idz);
  void getNodalSolutionPar(FemParPoisson* f, double* localSolPtr, int idz);
}

class FemParPoisson
{
 public:
  FemParPoisson(int nz, int ndim, int polyOrder, 
             double dz, bool periodicFlgs,
             bcdataPar_t bc[2], bool writeMatrix);
  ~FemParPoisson();
  void createGlobalSrc(double* ptr, int idz, double intSrcVol);
  void zeroGlobalSrc();
  void allreduceGlobalSrc(MPI_Comm comm);
  void allgatherGlobalStiff(MPI_Comm comm);
  void makeGlobalParStiffnessMatrix(double *laplacianWeight, double *modifierWeight, int idz);
  void finishGlobalParStiffnessMatrix();
  void solve();
  void getSolution(double* ptr, int idz);
  void getNodalSolution(double* ptr, int idz);

 private:
  const int nz, ndim, polyOrder;
  const double dz;
  const bool z_periodic;
  const bool writeMatrix;
  bcdataPar_t bc[2], bc1d[2];
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
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
  Eigen::MatrixXd localMassModToNod, localNodToMod, localModToNod;
  bool analyzed_; // flag so that stiffness matrix only analyzed once
  
  int getNumParGlobalNodes(int nz, int ndim, int p, bool periodicFlg);
  int getNumLocalNodes(int ndim, int p);
  
  void setupBoundaryIndices(bcdataPar_t bc[2], int ndim, int polyOrder);
  void getParLocalToGlobalInteriorBoundary(std::vector<int>& lgMap, int idz, int nz, int ndim, int p, bool periodicFlg);
};


#endif // FEM_PAR_POISSON_H
