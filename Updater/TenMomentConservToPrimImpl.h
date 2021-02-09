// Gkyl ----------------------------------------------------------------------------
//
// C++ back-end for ten-moment conservative to primitive (and vice versa) conversion
//    _______     ___
// + 6 @ |||| # P ||| +
//----------------------------------------------------------------------------------

#ifndef GK_TEN_MOMENT_CONSERV_TO_PRIM_H
#define GK_TEN_MOMENT_CONSERV_TO_PRIM_H

#include <stdint.h>

extern "C" {
  void gkylTenMomentConservativeToPrimitive(double *f);
  void gkylTenMomentPrimitiveToConservative(double *f);
}

#endif // GK_TEN_MOMENT_CONSERV_TO_PRIM_H
