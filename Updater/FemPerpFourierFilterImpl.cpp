// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for 2D FEM FourierFilter solver
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <FemPerpFourierFilterImpl.h>
#include <FemPerpPoissonImpl.h>
#include <FemMatrices.h>
#include <fftw3.h>

// std includes
#include <string>
#include <vector>
using namespace Eigen;

FemPerpFourierFilter::FemPerpFourierFilter(int nx_, int ny_, int ndim_, int *ikyFilter, int numFilter)
  : nx(nx_), ny(ny_), ndim(ndim_)
{
  data_r = (double*) fftw_malloc(sizeof(double)*nx*ny);
  data_c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*(ny/2+1));
  mask = (int*) malloc(sizeof(int)*(ny/2+1));

  if(ndim==2) {
    plan_r2c = fftw_plan_dft_r2c_2d(nx, ny, data_r, data_c, FFTW_ESTIMATE);
    plan_c2r = fftw_plan_dft_c2r_2d(nx, ny, data_c, data_r, FFTW_ESTIMATE);
    nlocal = 4;
  }

  localNodToMod = MatrixXd::Zero(nlocal, nlocal);
  localModToNod = MatrixXd::Zero(nlocal, nlocal);
  getPerpNodToModMatrix(localNodToMod, ndim, 1);
  localModToNod = localNodToMod.transpose().inverse().transpose();

  for (unsigned i=0; i<ny/2+1; i++) {
    mask[i] = 0;
  }
  for (unsigned n=0; n<numFilter; n++) {
    mask[ikyFilter[n]] = 1;
  }
}

FemPerpFourierFilter::~FemPerpFourierFilter() 
{
}

void FemPerpFourierFilter::assembleGlobalSrc(double *modalSrc, int idx, int idy) {
  // convert modal data to nodal data. only store first node.
  int index = idy + ny*idx;
  double norm = 1.0/(nx*ny);
  data_r[index] = 0.;
  for (unsigned m=0; m<nlocal; ++m) {
    data_r[index] += localModToNod(0,m)*modalSrc[m]*norm;
  }
}

void FemPerpFourierFilter::getFilteredSolution(double *modalSol, int idx, int idy) {
  // convert nodal data to modal data. 
  // since only bottom left corner node of each cell is stored in data_r,
  // need to accumulate data from multiple cells to compute modal expansion
  int index = idy + ny*idx;
  for (unsigned k=0; k<nlocal; ++k) {
    modalSol[k] = localNodToMod(k,0)*data_r[index];
    modalSol[k] += localNodToMod(k,1)*data_r[(idy%ny)+ny*((idx+1)%nx)];
    modalSol[k] += localNodToMod(k,2)*data_r[(idy+1)%ny+ny*(idx%nx)];
    modalSol[k] += localNodToMod(k,3)*data_r[(idy+1)%ny+ny*((idx+1)%nx)];
  }
}

void FemPerpFourierFilter::fft_r2c() {
  fftw_execute(plan_r2c);
}

void FemPerpFourierFilter::fft_c2r() {
  fftw_execute(plan_c2r);
}

void FemPerpFourierFilter::filter() {
  for (unsigned idy=0; idy<(ny/2+1); idy++) {
    for (unsigned idx=0; idx<nx; idx++) {
      unsigned i = idy + (ny/2+1)*idx;
      data_c[i][0] *= mask[idy];
      data_c[i][1] *= mask[idy];
    }
  }
}

// C wrappers for interfacing with FemPerpFourierFilter class
extern "C" void* new_FemPerpFourierFilter(int nx, int ny, int ndim, int *ikyFilter, int numFilter)
{
  FemPerpFourierFilter *f = new FemPerpFourierFilter(nx, ny, ndim, ikyFilter, numFilter);
  return reinterpret_cast<void*>(f);
}

extern "C" void delete_FemPerpFourierFilter(FemPerpFourierFilter* f)
{
  delete f;
}

extern "C" void assembleGlobalSrc(FemPerpFourierFilter* f, double *data, int idx, int idy) {
  f->assembleGlobalSrc(data, idx, idy);
}

extern "C" void getFilteredSolution(FemPerpFourierFilter* f, double *data, int idx, int idy) {
  f->getFilteredSolution(data, idx, idy);
}

extern "C" void fft_r2c(FemPerpFourierFilter* f) {
  f->fft_r2c();
}

extern "C" void fft_c2r(FemPerpFourierFilter* f) {
  f->fft_c2r();
}

extern "C" void filter(FemPerpFourierFilter* f) {
  f->filter();
}
