// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for the discontinuous Poisson solver.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <DiscontPoissonImpl.h>

// std includes
#include <string>
#include <vector>

using namespace Eigen;

DiscontPoisson::DiscontPoisson(int _ncell[3], int _ndim, int _nbasis,
                               int _nnonzero, int _polyOrder, bool _writeMatrix)
  : ndim(_ndim), nbasis(_nbasis),
    nnonzero(_nnonzero), polyOrder(_polyOrder),
    writeMatrix(_writeMatrix)
{

  int vol = 1;
  for (int d = 0; d < ndim; ++d) {
    ncell[d] = _ncell[d];
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
  if (writeMatrix) {
    saveMarket(stiffMat, "stiffMat.mm");
  }
  // de-allocate row major copy
  stiffMatRowMajor.resize(0,0);
  solver.analyzePattern(stiffMat);
  //saveMarket(stiffMat, "matrix");
  solver.factorize(stiffMat);
}

void DiscontPoisson::pushSource(int idx, double* src, double* srcMod) {
  for (int k = 0; k < nbasis; ++k) {
    globalSrc.coeffRef(idx+k) = -src[k] + srcMod[k];
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
extern "C" void* new_DiscontPoisson(int ncell[3], int ndim, int nbasis,
                                    int nnonzero, int polyOrder, bool writeMatrix)
{
  DiscontPoisson *f = new DiscontPoisson(ncell, ndim, nbasis,
                                         nnonzero, polyOrder, writeMatrix);
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

extern "C" void discontPoisson_pushSource(DiscontPoisson* f, int idx, double* src, double* srcMod)
{
  f->pushSource(idx, src, srcMod);
}

extern "C" void discontPoisson_getSolution(DiscontPoisson* f, int idx, double* sol)
{
  f->getSolution(idx, sol);
}
