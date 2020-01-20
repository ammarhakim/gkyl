// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for 2D discontinuous Poisson solver
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <DiscontPoissonImpl.h>

// std includes
#include <string>
#include <vector>

using namespace Eigen;
static const int DIRICHLET_BC = 0;
static const int DIRICHLET_VARIABLE_BC = 2;
static const int DX = 0;
static const int DY = 1;
static const int LO = 0;
static const int HI = 1;


DiscontPoisson::DiscontPoisson(int ncell_[3], int ndim_, int nbasis_,
                               int nnonzero_, int polyOrder_, double dx_[3])
  : ndim(ndim_), nbasis(nbasis_),
    nnonzero(nnonzero_), polyOrder(polyOrder_)
{

  int vol = 1;
  for (int d = 0; d < ndim; ++d) {
    ncell[d] = ncell_[d];
    dx[d] = dx_[d];
    vol = vol*ncell[d];
  }
  N = vol*nbasis;

  stiffMatRowMajor = SparseMatrix<double,RowMajor>(N, N);
  stiffTripletList.reserve(vol*nnonzero); // estimate number of nonzero elements
  globalSrc = VectorXd::Zero(N);
}

DiscontPoisson::~DiscontPoisson() 
{
  stiffMat.resize(0,0);
  globalSrc.resize(0);
  x.resize(0);
}

void DiscontPoisson::pushTriplet(int i, int j, double val) {
  stiffTripletList.push_back(Triplet<double>(i, j, val));
}

void DiscontPoisson::constructStiffMatrix() {
  stiffMatRowMajor.setFromTriplets(stiffTripletList.begin(), stiffTripletList.end());
  stiffTripletList.resize(0);

  // create column major copy of stiffMat so that we can zero columns
  stiffMat = SparseMatrix<double,ColMajor>(stiffMatRowMajor);
  // de-allocate row major copy
  stiffMatRowMajor.resize(0,0);
  solver.analyzePattern(stiffMat);
  saveMarket(stiffMat, "matrix");
  solver.factorize(stiffMat);
}

void DiscontPoisson::pushSource(int idx, double* src) {
  for (int k = 0; k < nbasis; ++k) {
    globalSrc.coeffRef(idx+k) = src[k];
  }
}

void DiscontPoisson::getSolution(int idx, double* sol) {
  for (int k = 0; k < nbasis; ++k) {
    sol[k] = x.coeffRef(idx+k);
  }
}

void DiscontPoisson::solve() {
  //saveMarket(globalSrc, "source");
  x = VectorXd::Zero(N);
  x = solver.solve(globalSrc);
  //saveMarket(x, "solution");
}


// C wrappers for interfacing with DiscontPoisson class
extern "C" void* new_DiscontPoisson(int ncell[3], int ndim, int nbasis, int nnonzero, int polyOrder, double dx[3])
{
  DiscontPoisson *f = new DiscontPoisson(ncell, ndim, nbasis, nnonzero, polyOrder, dx);
  return reinterpret_cast<void*>(f);
}

extern "C" void delete_DiscontPoisson(DiscontPoisson* f)
{
  delete f;
}

extern "C" void discontPoisson_solve(DiscontPoisson* f)
{
  f->solve();
}

extern "C" void discontPoisson_pushTriplet(DiscontPoisson* f, int i, int j, double val)
{
  f->pushTriplet(i, j, val);
}

extern "C" void discontPoisson_constructStiffMatrix(DiscontPoisson* f)
{
  f->constructStiffMatrix();
}

extern "C" void discontPoisson_pushSource(DiscontPoisson* f, int idx, double* src)
{
  f->pushSource(idx, src);
}

extern "C" void discontPoisson_getSolution(DiscontPoisson* f, int idx, double* sol)
{
  f->getSolution(idx, sol);
}


/*void DiscontPoisson::pushTripletSet(int idxX, int idxY) {
  stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 0), indexer2Dto1D(idxX, idxY, 0), -3.0));
  stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 1), indexer2Dto1D(idxX, idxY, 1), -23.0/2.0));
  stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 2), indexer2Dto1D(idxX, idxY, 2), -23.0/2.0));
  stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 3), indexer2Dto1D(idxX, idxY, 3), -20.0));

  if (idxX > 0) {
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 0), indexer2Dto1D(idxX-1, idxY, 0), 3.0/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 0), indexer2Dto1D(idxX-1, idxY, 1), sqrtf(3.0)/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 1), indexer2Dto1D(idxX-1, idxY, 0), -5.0*sqrtf(3.0)/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 1), indexer2Dto1D(idxX-1, idxY, 1), -5.0/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 2), indexer2Dto1D(idxX-1, idxY, 2), 3.0/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 2), indexer2Dto1D(idxX-1, idxY, 3), sqrtf(3.0)/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 3), indexer2Dto1D(idxX-1, idxY, 2), -5.0*sqrtf(3.0)/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 3), indexer2Dto1D(idxX-1, idxY, 3), -5.0/2.0));
  }

  if (idxX < nx-1) {
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 0), indexer2Dto1D(idxX+1, idxY, 0), 3.0/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 0), indexer2Dto1D(idxX+1, idxY, 1), -sqrtf(3.0)/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 1), indexer2Dto1D(idxX+1, idxY, 0), 5.0*sqrtf(3.0)/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 1), indexer2Dto1D(idxX+1, idxY, 1), -5.0/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 2), indexer2Dto1D(idxX+1, idxY, 2), 3.0/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 2), indexer2Dto1D(idxX+1, idxY, 3), -sqrtf(3.0)/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 3), indexer2Dto1D(idxX+1, idxY, 2), 5.0*sqrt(3.0)/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 3), indexer2Dto1D(idxX+1, idxY, 3), -5.0/2.0));
  }

  if (idxY < ny-1) {
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 0), indexer2Dto1D(idxX, idxY+1, 0), 3.0/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 0), indexer2Dto1D(idxX, idxY+1, 2), -sqrtf(3.0)/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 1), indexer2Dto1D(idxX, idxY+1, 1), 3.0/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 1), indexer2Dto1D(idxX, idxY+1, 3), -sqrt(3.0)/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 2), indexer2Dto1D(idxX, idxY+1, 0), 5.0*sqrtf(3.0)/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 2), indexer2Dto1D(idxX, idxY+1, 2), -5.0/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 3), indexer2Dto1D(idxX, idxY+1, 1), 5.0*sqrt(3.0)/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 3), indexer2Dto1D(idxX, idxY+1, 3), -5.0/2.0));
  }

  if (idxY > 0) {
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 0), indexer2Dto1D(idxX, idxY-1, 0), 3.0/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 0), indexer2Dto1D(idxX, idxY-1, 2), sqrtf(3.0)/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 1), indexer2Dto1D(idxX, idxY-1, 1), 3.0/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 1), indexer2Dto1D(idxX, idxY-1, 3), sqrt(3.0)/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 2), indexer2Dto1D(idxX, idxY-1, 0), -5.0*sqrt(3.0)/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 2), indexer2Dto1D(idxX, idxY-1, 2), -5.0/2.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 3), indexer2Dto1D(idxX, idxY-1, 1), -5.0*sqrt(3.0)/4.0));
    stiffTripletList.push_back(Triplet<double>(indexer2Dto1D(idxX, idxY, 3), indexer2Dto1D(idxX, idxY-1, 3), -5.0/2.0));
  }
}*/

