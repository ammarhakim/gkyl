// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for FemKyFourierFilter solver
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef FEM_FOURIER_FILTER_H
#define FEM_FOURIER_FILTER_H

#include <stdint.h>
#include <Eigen/Core>

// fftw includes
#include <fftw3.h>
#include <mpi.h>

class FemKyFourierFilter;

extern "C" {
// C wrappers for interfacing with FemKyFourierFilter class
  void* new_FemKyFourierFilter(int nx, int ny, int nz, int ndim, int *ikyFilter, int numFilter);
  void delete_FemKyFourierFilter(FemKyFourierFilter* f);
  void assembleGlobalSrc(FemKyFourierFilter* f, double *modalSrc, int idx, int idy, int idz);
  void getFilteredSolution(FemKyFourierFilter* f, double *modalSol, int idx, int idy, int idz);
  void fft_r2c(FemKyFourierFilter* f);
  void fft_c2r(FemKyFourierFilter* f);
  void filter(FemKyFourierFilter* f);
}

class FemKyFourierFilter
{
  public:
    FemKyFourierFilter(int nx, int ny, int nz, int ndim, int *ikyFilter, int numFilter);
    ~FemKyFourierFilter();
    void assembleGlobalSrc(double *modalSrc, int idx, int idy, int idz);
    void getFilteredSolution(double *modalSol, int idx, int idy, int idz);
    void fft_r2c();
    void fft_c2r();
    void filter();
  
  private:
    const int nx, ny, nz, ndim;
    int nlocal;
    fftw_plan plan_r2c, plan_c2r;
    fftw_complex *data_c;
    double *data_r;
    int *mask;
    Eigen::MatrixXd localNodToMod, localModToNod;
};


#endif // FEM_FOURIER_FILTER_H
