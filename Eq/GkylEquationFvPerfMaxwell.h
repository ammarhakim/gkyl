// Gkyl ------------------------------------------------------------------------
//
// Perfectly-Hyperbolic Maxwell equations object for wave-propagation
// finite-volume method
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

#include <GkylEquationFv.h>

extern "C" {
  typedef struct {
    int numWaves;
    int numEquations;
    double lightSpeed;
    double elcErrorSpeedFactor;
    double mgnErrorSpeedFactor;
  } GkylEquationFvPerfMaxwell_t;

  __host__ __device__ void PerfMaxwell_rp(
      const void * __restrict__ eqn,
      const int dir,
      const double * __restrict__ ql,
      const double * __restrict__ qr,
      double * __restrict__ waves,
      double * __restrict__ speeds);

  __host__ __device__ void PerfMaxwell_qFluctuations(
      const void * __restrict__ eqn,
      const int dir,
      const double * __restrict__ ql,
      const double * __restrict__ qr,
      const double * __restrict__ waves,
      const double * __restrict__ speeds,
      double * __restrict__ amdq,
      double * __restrict__ apdq);

  __host__ __device__ void PerfMaxwell_flux(
      const void * __restrict__ eqn,
      const int dir,
      const double * __restrict__ qIn,
      double * __restrict__ fOut);

  GkylEquationFvPerfMaxwell_t *new_EquationFvPerfMaxwellOnHost(
      const double lightSpeed, const double elcErrorSpeedFactor,
      const double mgnErrorSpeedFactor);

  GkylEquationFv_t *new_EquationFvPerfMaxwellOnDevice(
      const double lightSpeed, const double elcErrorSpeedFactor,
      const double mgnErrorSpeedFactor);
}
