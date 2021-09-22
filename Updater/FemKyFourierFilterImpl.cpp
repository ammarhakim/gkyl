// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for 2D FEM FourierFilter solver
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <FemKyFourierFilterImpl.h>
#include <FemPerpPoissonImpl.h>
#include <FemMatrices.h>
#include <fftw3.h>

// std includes
#include <string>
#include <vector>
using namespace Eigen;

FemKyFourierFilter::FemKyFourierFilter(int nx_, int ny_, int nz_, int ndim_, int *ikyFilter, int numFilter)
  : nx(nx_), ny(ny_), nz(nz_), ndim(ndim_)
{
  if(ndim==2) {
    nlocal = 4;
  } else if(ndim==3) {
    nlocal = 8;
  }

  data_r = (double*) fftw_malloc(sizeof(double)*nx*ny*nz);
  data_c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*nx*(ny/2+1)*nz);
  mask = (int*) malloc(sizeof(int)*(ny/2+1));

  int rank = 1; // only compute 1d transform (in y)
  int n[] = {ny}; // 1d transform in y of length ny
  int howmany = nx*nz;
  int idist = ny;
  int odist = ny/2+1;
  int istride = 1;
  int ostride = 1;
  int *inembed = n;
  int *onembed = n;
  plan_r2c = fftw_plan_many_dft_r2c(rank, n, howmany, data_r, NULL, istride, idist, data_c, NULL, ostride, odist, FFTW_ESTIMATE);
  plan_c2r = fftw_plan_many_dft_c2r(rank, n, howmany, data_c, NULL, ostride, odist, data_r, NULL, istride, idist, FFTW_ESTIMATE);

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

FemKyFourierFilter::~FemKyFourierFilter() 
{
}

void FemKyFourierFilter::assembleGlobalSrc(double *modalSrc, int idx, int idy, int idz) {
  // convert modal data to nodal data. only store first node.
  if(ndim==2) idz = 0;
  int index = idy + ny*idx + ny*nx*idz; // use ordering with y the fastest index
  double norm = 1.0/(ny);
  data_r[index] = 0.;
  for (unsigned m=0; m<nlocal; ++m) {
    data_r[index] += localModToNod(0,m)*modalSrc[m]*norm;
  }

  // handle x and z upper ghosts
  if(idx==nx-2) { // if on last non-ghost x cell row
    index = idy + ny*(idx+1) + ny*nx*idz;
    data_r[index] = 0.;
    for (unsigned m=0; m<nlocal; ++m) {
      data_r[index] += localModToNod(1,m)*modalSrc[m]*norm;
    }
  }
  if(idz==nz-2) { // if on last non-ghost z cell row
    index = idy + ny*idx + ny*nx*(idz+1);
    data_r[index] = 0.;
    for (unsigned m=0; m<nlocal; ++m) {
      data_r[index] += localModToNod(4,m)*modalSrc[m]*norm;
    }
  }
  if(idx==nx-2 && idz==nz-2) { // if on upper x-z corner cell
    index = idy + ny*(idx+1) + ny*nx*(idz+1);
    data_r[index] = 0.;
    for (unsigned m=0; m<nlocal; ++m) {
      data_r[index] += localModToNod(5,m)*modalSrc[m]*norm;
    }
  }
}

void FemKyFourierFilter::getFilteredSolution(double *modalSol, int idx, int idy, int idz) {
  // convert nodal data to modal data. 
  // since only bottom left corner node of each cell is stored in data_r,
  // need to accumulate data from multiple cells to compute modal expansion
  // assume y is periodic, so use modulo ny when necessary
  if(ndim==2) idz = 0;
  int index = idy + ny*idx + ny*nx*idz;
  for (unsigned k=0; k<nlocal; ++k) {
    modalSol[k] = localNodToMod(k,0)*data_r[index];
    modalSol[k] += localNodToMod(k,1)*data_r[(idy%ny)+ny*(idx+1)+ny*nx*idz];
    modalSol[k] += localNodToMod(k,2)*data_r[(idy+1)%ny+ny*idx+ny*nx*idz];
    modalSol[k] += localNodToMod(k,3)*data_r[(idy+1)%ny+ny*(idx+1)+ny*nx*idz];
    if(ndim==3) {
      modalSol[k] += localNodToMod(k,4)*data_r[(idy%ny)+ny*(idx)+ny*nx*(idz+1)];
      modalSol[k] += localNodToMod(k,5)*data_r[(idy%ny)+ny*(idx+1)+ny*nx*(idz+1)];
      modalSol[k] += localNodToMod(k,6)*data_r[(idy+1)%ny+ny*idx+ny*nx*(idz+1)];
      modalSol[k] += localNodToMod(k,7)*data_r[(idy+1)%ny+ny*(idx+1)+ny*nx*(idz+1)];
    }
  }
}

void FemKyFourierFilter::fft_r2c() {
  fftw_execute(plan_r2c);
}

void FemKyFourierFilter::fft_c2r() {
  fftw_execute(plan_c2r);
}

void FemKyFourierFilter::filter() {
  for (unsigned idy=0; idy<(ny/2+1); idy++) {
    for (unsigned idxz=0; idxz<nx*nz; idxz++) {
      unsigned i = idy + (ny/2+1)*idxz;
      //printf("%d  %d  %f\n", idx, idy, data_c[i][0]);
      data_c[i][0] *= mask[idy];
      data_c[i][1] *= mask[idy];
    }
  }
}

// C wrappers for interfacing with FemKyFourierFilter class
extern "C" void* new_FemKyFourierFilter(int nx, int ny, int nz, int ndim, int *ikyFilter, int numFilter)
{
  FemKyFourierFilter *f = new FemKyFourierFilter(nx, ny, nz, ndim, ikyFilter, numFilter);
  return reinterpret_cast<void*>(f);
}

extern "C" void delete_FemKyFourierFilter(FemKyFourierFilter* f)
{
  delete f;
}

extern "C" void assembleGlobalSrc(FemKyFourierFilter* f, double *data, int idx, int idy, int idz) {
  f->assembleGlobalSrc(data, idx, idy, idz);
}

extern "C" void getFilteredSolution(FemKyFourierFilter* f, double *data, int idx, int idy, int idz) {
  f->getFilteredSolution(data, idx, idy, idz);
}

extern "C" void fft_r2c(FemKyFourierFilter* f) {
  f->fft_r2c();
}

extern "C" void fft_c2r(FemKyFourierFilter* f) {
  f->fft_c2r();
}

extern "C" void filter(FemKyFourierFilter* f) {
  f->filter();
}
