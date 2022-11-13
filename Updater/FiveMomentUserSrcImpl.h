//------------------------------------------------------------------------------
// C++ back-end for inter-species friction for the five-moment (Euler) model.
//------------------------------------------------------------------------------

#ifndef GK_FIVE_MOMENT_FRICTION_SRC_H
#define GK_FIVE_MOMENT_FRICTION_SRC_H

#include <stdint.h>

extern "C" {
typedef struct {
  bool evolve;
  double mass;
} FiveMomentUserSrcFluidData_t;

typedef struct {
  int nFluids;
  double gasGamma;
  double kBoltzmann;
} FiveMomentUserSrcData_t;

void
gkylFiveMomentUserSrcForwardEuler(
  const FiveMomentUserSrcData_t *sd,
  const FiveMomentUserSrcFluidData_t *fd,
  const double dt,
  double *qPtrs[],
  const double *sourcePtrs[]);
}

#endif // GK_FIVE_MOMENT_FRICTION_SRC_H
