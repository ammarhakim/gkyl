// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for 2D FEM Fourier filter 
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

class FemPerpFourierFilter;

extern "C" {
// C wrappers for interfacing with FemPerpFourierFilter class
  void* new_FemPerpFourierFilter(int nx, int ny, int nz, int ndim, int *ikyFilter, int numFilter);
  void delete_FemPerpFourierFilter(FemPerpFourierFilter* f);
  void assembleGlobalSrc(FemPerpFourierFilter* f, double *modalSrc, int idx, int idy, int idz);
  void getFilteredSolution(FemPerpFourierFilter* f, double *modalSol, int idx, int idy, int idz);
  void fft_r2c(FemPerpFourierFilter* f);
  void fft_c2r(FemPerpFourierFilter* f);
  void filter(FemPerpFourierFilter* f);
}

class FemPerpFourierFilter
{
  public:
    FemPerpFourierFilter(int nx, int ny, int nz, int ndim, int *ikyFilter, int numFilter);
    ~FemPerpFourierFilter();
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
