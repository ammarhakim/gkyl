// Gkyl ------------------------------------------------------------------------
//
// Basic equation object
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

#include <GkylCudaConfig.h>

extern "C"  {

    // Function pointer types for volume and surface terms
    typedef double (*volTermFunc_t)(const void* __restrict__ self, 
      const double* __restrict__ xc, const double* __restrict__ dx, const int* __restrict__ idx, const double* __restrict__ qIn, double *qRhsOut);
    
    typedef double (*surfTermFunc_t)(const void* __restrict__ self, int dir,
      const double* __restrict__ xcL, const double* __restrict__ xcR, const double* __restrict__ dxL, const double* __restrict__ dxR,
      double maxsOld, const int* __restrict__ idxL, const int* __restrict__ idxR,
      const double* __restrict__ qInL, const double* __restrict__ qInR, double *qRhsOutL, double *qRhsOutR);

    typedef double (*boundarySurfTermFunc_t)(const void* __restrict__ self, int dir,
      const double* __restrict__ xcL, const double* __restrict__ xcR, const double* __restrict__ dxL, const double* __restrict__ dxR,
      double maxsOld, const int* __restrict__ idxL, const int* __restrict__ idxR,
      const double* __restrict__ qInL, const double* __restrict__ qInR, double *qRhsOutL, double *qRhsOutR);

    typedef struct {
        // Pointer to specific equation object
        void *equation;
        // pointer to equation volTerm function
        volTermFunc_t equationVolTerm;
        // pointer to equation surfTerm function
        surfTermFunc_t equationSurfTerm;
         // pointer to equation boundarySurfTerm function
        boundarySurfTermFunc_t equationBoundarySurfTerm;
        
        __host__ __device__ double volTerm(const double* __restrict__ xc, const double* __restrict__ dx, const int* __restrict__ idx, const double* __restrict__ qIn, double *qRhsOut) {
          return equationVolTerm(equation, xc, dx, idx, qIn, qRhsOut);
        }

        __host__ __device__ double surfTerm(int dir,
          const double* __restrict__ xcL, const double* __restrict__ xcR, const double* __restrict__ dxL, const double* __restrict__ dxR,
          double maxsOld, const int* __restrict__ idxL, const int* __restrict__ idxR,
          const double* __restrict__ qInL, const double* __restrict__ qInR, double *qRhsOutL, double *qRhsOutR) {
          return equationSurfTerm(equation, dir,
            xcL, xcR, dxL, dxR,
            maxsOld, idxL, idxR,
            qInL, qInR, qRhsOutL, qRhsOutR);
        }

        __host__ __device__ double boundarySurfTerm(int dir,
          const double* __restrict__ xcL, const double* __restrict__ xcR, const double* __restrict__ dxL, const double* __restrict__ dxR,
          double maxsOld, const int* __restrict__ idxL, const int* __restrict__ idxR,
          const double* __restrict__ qInL, const double* __restrict__ qInR, double *qRhsOutL, double *qRhsOutR) {
          return equationBoundarySurfTerm(equation, dir,
            xcL, xcR, dxL, dxR,
            maxsOld, idxL, idxR,
            qInL, qInR, qRhsOutL, qRhsOutR);
        }
    } GkylEquation_t;
}
