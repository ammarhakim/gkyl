#ifndef GKYL_EULER_H
#define GKYL_EULER_H

// Gkyl includes
#include <GkylCudaConfig.h>
#include <GkylCartField.h>

// std includes
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace Gkyl {

  class Euler;

  /* C wrappers to member functions, so that they can be called from Lua */
  extern "C" {
    void* new_Euler();
    void* new_Euler_onDevice(Euler *v);
  }

  // FIXME correct specifier?
  __device__ __constant__ const int dirShuffle[3][4] = {{0,1,2,3}, {0,2,3,1}, {0,3,1,2}};

  class Euler {
    public:
      __host__ __device__ Euler() {}
      ~Euler() = default;

      __host__ __device__ int getNumWaves();

      __host__ __device__ int getNumEquations();

      __host__ __device__ double gasGamma();

      __host__ __device__ double pressure(const double *q);

      __host__ __device__ void rp(
          const unsigned dir, const double *delta, const double *ql,
          const double *qr, double *waves, double *s);

      __host__ __device__ void qFluctuations(
          const unsigned dir, const double *ql, const double *qr,
          const double *waves, const double *s, double *amdq,
          double *apdq);

    private:
      /* dimension and basis parameters */
      const double _gasGamma = 5./3.;
      unsigned _numWaves = 3;
      const unsigned _numEquations = 5;


      double _fl[5], _fr[5]; /* Storage for left/right fluxes */

      __host__ __device__ void rpLax(
          const unsigned dir, const double *delta, const double *ql,
          const double *qr, double *waves, double *s);

      __host__ __device__ void qFluctuationsLax(
          const unsigned dir, const double *ql, const double *qr,
          const double *waves, const double *s, double *amdq,
          double *apdq);
  };
}
#endif
