// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for geometric source for axisymmetric perfectly-hyperbolic
// Maxwell's equations
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_AXISYMMETRIC_FIVE_MOMENT_SRC_H
#define GK_AXISYMMETRIC_FIVE_MOMENT_SRC_H

#include <stdint.h>

extern "C" {
  typedef struct {
    double epsilon0;
    double mu0;
    double chi_e;
    double chi_m;
  } AxisymmetricPhMaxwellSrcData_t;

  void gkylAxisymmetricPhMaxwellSrcForwardEuler(
    const AxisymmetricPhMaxwellSrcData_t *sd,
    const double dt,
    const double *xc,
    double *fPtr);

  void
  gkylAxisymmetricPhMaxwellSrcRk3(const AxisymmetricPhMaxwellSrcData_t *sd,
                                  const double dt,
                                  const double *xc,
                                  double *fPtr);
}

#endif // GK_AXISYMMETRIC_FIVE_MOMENT_SRC_H
