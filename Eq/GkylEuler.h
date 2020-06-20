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
    void* new_Euler(const double gasGamma);
    void* new_Euler_onDevice(Euler *v);

    int numEquations_Euler(Euler *v);
    void rp_Euler(
        Euler *v, const int dir, const double *delta, const double *ql,
        const double *qr, double *waves, double *s);
    void qFluctuations_Euler(
        Euler *eq, const int dir, const double *ql, const double *qr,
        const double *waves, const double *s, double *amdq, double *apdq);
    void flux_Euler(Euler *eq, const int dir, const double *qIn, double *fOut);
  }

  // FIXME correct specifier?
  __device__ __constant__ const int dirShuffle[3][4] = {{0,1,2,3}, {0,2,3,1}, {0,3,1,2}};

  class Euler {
    public:
      __host__ __device__ Euler(const double gasGamma);
      ~Euler() = default;

      __host__ __device__ int numWaves();

      __host__ __device__ int numEquations();

      __host__ __device__ double gasGamma();

      __host__ __device__ double pressure(const double *q);

      __host__ __device__ void rp(
          const int dir, const double *delta, const double *ql,
          const double *qr, double *waves, double *s);

      __host__ __device__ void qFluctuations(
          const int dir, const double *ql, const double *qr,
          const double *waves, const double *s, double *amdq,
          double *apdq);

      __host__ __device__ void flux(
          const int dir, const double *qIn, double *fOut);

    private:
      int _numWaves = 1;
      const int _numEquations = 5;
      const double _gasGamma = 5./3.;

      __host__ __device__ void rpLax(
          const int dir, const double *delta, const double *ql,
          const double *qr, double *waves, double *s);

      __host__ __device__ void qFluctuationsLax(
          const int dir, const double *ql, const double *qr,
          const double *waves, const double *s, double *amdq,
          double *apdq);
  };
}
#endif
