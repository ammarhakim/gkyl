// Gkyl ------------------------------------------------------------------------
//
// Basic equation object for wave-propagation finite-volume method
//  _______   ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

#include <GkylCudaConfig.h>
#include <GkylCudaFuncs.h>

extern "C"  {
  __device__ __constant__ const int dirShuffle[3][4] = {
    {0,1,2,3}, {0,2,3,1}, {0,3,1,2}};

  typedef void (*rp_t)(
      const void * __restrict__ eqn,
      const int dir,
      const double * __restrict__ ql,
      const double * __restrict__ qr,
      double * __restrict__ waves,
      double * __restrict__ s);

  typedef void (*qFluctuations_t)(
      const void * __restrict__ eqn,
      const int dir,
      const double * __restrict__ ql,
      const double * __restrict__ qr,
      const double * __restrict__ waves,
      const double * __restrict__ s,
      double * __restrict__ amdq,
      double * __restrict__ apdq);
    
  typedef void (*flux_t)(
      const void * __restrict__ eqn,
      const int dir,
      const double * __restrict__ qIn,
      double * __restrict__ fOut);

  typedef struct {
    int numWaves;
    int numEquations;
    void *equation;
    rp_t equationRp;
    qFluctuations_t equationQFluctuations;
    flux_t equationFlux;

    __host__ __device__ void rp(
        const int dir,
        const double * __restrict__ ql,
        const double * __restrict__ qr,
        double * __restrict__ waves,
        double * __restrict__ s)
    {
      equationRp(equation, dir, ql, qr, waves, s);
    }

    __host__ __device__ void qFluctuations(
        const int dir,
        const double * __restrict__ ql,
        const double * __restrict__ qr,
        const double * __restrict__ waves,
        const double * __restrict__ s,
        double * __restrict__ amdq,
        double * __restrict__ apdq)
    {
      equationQFluctuations(equation, dir, ql, qr, waves, s, amdq, apdq);
    }

    __host__ __device__ void flux(
        const int dir,
        const double * __restrict__ qIn,
        double * __restrict__ fOut)
    {
      equationFlux(equation, dir, qIn, fOut);
    }
  } GkylEquationFv_t;
}
