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
    double gasGamma;
  } GkylEquationFvEuler_t;

  GkylEquationFv_t *new_EquationFvEulerOnDevice(const double gasGamma);
}
