// Gkyl ------------------------------------------------------------------------
//
// Euler equation object for wave-propagation finite-volume method
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

#include <GkylEquationFv.h>

extern "C" {
  typedef struct {
    int numWaves;
    int numEquations;
    double gasGamma;
  } GkylEquationFvEuler_t;

  __host__ __device__ inline double Euler_pressure(
      const double gasGamma,
      const double * __restrict__ q);

  __host__ __device__ void Euler_rp(
      const void * __restrict__ eqn,
      const int dir,
      const double * __restrict__ ql,
      const double * __restrict__ qr,
      double * __restrict__ waves,
      double * __restrict__ speeds);

  __host__ __device__ void Euler_qFluctuations(
      const void * __restrict__ eqn,
      const int dir,
      const double * __restrict__ ql,
      const double * __restrict__ qr,
      const double * __restrict__ waves,
      const double * __restrict__ speeds,
      double * __restrict__ amdq,
      double * __restrict__ apdq);

  __host__ __device__ void Euler_flux(
      const void * __restrict__ eqn,
      const int dir,
      const double * __restrict__ qIn,
      double * __restrict__ fOut);

  GkylEquationFvEuler_t *new_EquationFvEulerOnHost(const double gasGamma);

  GkylEquationFv_t *new_EquationFvEulerOnDevice(const double gasGamma);
}
