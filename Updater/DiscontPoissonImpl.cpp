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

static double take_last(const double &in,const double &b) { return b; }
static void vectorSum(double *in, double *inout, int *len, MPI_Datatype *dptr)
{
  int i;
  for(i=0; i< *len; ++i) {
    inout[i] += in[i];
  }
}

DiscontPoisson::DiscontPoisson(int nx_, int ny_, int ndim_, int polyOrder_, 
                               double dx_, double dy_, bool periodicFlgs_[2], 
                               bcdata_t bc_[2][2])
  : nx(nx_), ny(ny_), ndim(ndim_), 
    polyOrder(polyOrder_), dx(dx_), dy(dy_)
{
  nb = getNumLocalNodes(ndim, polyOrder);
  N = nx*ny*nb;

  stiffMatRowMajor = SparseMatrix<double,RowMajor>(N, N);
  stiffTripletList.reserve(nx*ny*36); // estimate number of nonzero elements
  globalSrc = VectorXd::Zero(N);
  
  // copy to input to class structures
  for(int i=0; i<2; i++) {
    periodicFlgs[i] = periodicFlgs_[i];
    for(int j=0; j<2; j++) {
      bc[i][j] = bc_[i][j];
      if(!bc[i][j].isSet) {
        bc[i][j].type = -1;
      }
    }
  }

  allPeriodic = false;
  if(periodicFlgs[0] && periodicFlgs[1]) allPeriodic = true;
  
  analyzed_ = false; // flag so that stiffness matrix only analyzed once
}

DiscontPoisson::~DiscontPoisson() 
{
  stiffMat.resize(0,0);
  globalSrc.resize(0);
  x.resize(0);
}

int DiscontPoisson::getNumLocalNodes(int ndim, int p) 
{
  int numNodes = -1;
  if(ndim==2) {
    if(p==1) numNodes = 4;
    else if(p==2) numNodes = 8;
  }
  else if(ndim==3) {
    if(p==1) numNodes = 8;
    else if(p==2) numNodes = 20;
  }
  return numNodes;
}

unsigned DiscontPoisson::indexer2Dto1D(int i, int j, int basisIdx) {
  return i*ny*nb + j*nb + basisIdx;
}

void DiscontPoisson::pushTripletSet(int idxX, int idxY) {
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
}

void DiscontPoisson::constructStiffMatrix() {
  stiffMatRowMajor.setFromTriplets(stiffTripletList.begin(), stiffTripletList.end());
  stiffTripletList.resize(0);

  // create column major copy of stiffMat so that we can zero columns
  stiffMat = SparseMatrix<double,ColMajor>(stiffMatRowMajor);
  // de-allocate row major copy
  stiffMatRowMajor.resize(0,0);
  
  if(!analyzed_) {
    // only do analyzePattern once, assuming structure of matrix doesn't change
    solver.analyzePattern(stiffMat);
    analyzed_ = true;
  }
  saveMarket(stiffMat, "matrix");
  solver.factorize(stiffMat);
}

void DiscontPoisson::pushSource(int idxX, int idxY, double* src) {
  for (int k = 0; k < nb; ++k) {
    globalSrc.coeffRef(indexer2Dto1D(idxX, idxY, k)) = src[k];
  }
}

void DiscontPoisson::getSolution(int idxX, int idxY, double* sol) {
  for (int k = 0; k < nb; ++k) {
    sol[k] = x.coeffRef(indexer2Dto1D(idxX, idxY, k));
  }
}

void DiscontPoisson::solve() {
  saveMarket(globalSrc, "source");
  x = VectorXd::Zero(N);
  x = solver.solve(globalSrc);
  saveMarket(x, "solution");
}


// C wrappers for interfacing with DiscontPoisson class
extern "C" void* new_DiscontPoisson(int nx, int ny, int ndim, int polyOrder, double dx, double dy, bool periodicFlgs[2], bcdata_t bc[2][2])
{
  DiscontPoisson *f = new DiscontPoisson(nx, ny, ndim, polyOrder, dx, dy, periodicFlgs, bc);
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

extern "C" void discontPoisson_pushTripletSet(DiscontPoisson* f, int idxX, int idxY)
{
  f->pushTripletSet(idxX, idxY);
}

extern "C" void discontPoisson_constructStiffMatrix(DiscontPoisson* f)
{
  f->constructStiffMatrix();
}

extern "C" void discontPoisson_pushSource(DiscontPoisson* f, int idxX, int idxY, double* src)
{
  f->pushSource(idxX, idxY, src);
}

extern "C" void discontPoisson_getSolution(DiscontPoisson* f, int idxX, int idxY, double* sol)
{
  f->getSolution(idxX, idxY, sol);
}
