// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for geometric source for axisymmetric five-moment (Euler) model
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_AXISYMMETRIC_FIVE_MOMENT_SRC_H
#define GK_AXISYMMETRIC_FIVE_MOMENT_SRC_H

#include <stdint.h>

extern "C" {
  typedef struct {
    bool evolve;
  } AxisymmetricFluidData_t;

  typedef struct {
    int nFluids;
    double gasGamma;
    bool hasPressure;
  } AxisymmetricFiveMomentSrcData_t;

  void gkylAxisymmetricFiveMomentSrcForwardEuler(
    const AxisymmetricFiveMomentSrcData_t *sd,
    const AxisymmetricFluidData_t *fd,
    const double dt,
    const double *xc,
    double **fPtrs);

  void
  gkylAxisymmetricFiveMomentSrcRk3(const AxisymmetricFiveMomentSrcData_t *sd,
                                   const AxisymmetricFluidData_t *fd,
                                   const double dt,
                                   const double *xc,
                                   double **fPtrs);
}

#endif // GK_AXISYMMETRIC_FIVE_MOMENT_SRC_H
