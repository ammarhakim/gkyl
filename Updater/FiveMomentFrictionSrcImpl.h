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
  } FiveMomentFrictionFluidData_t;

  typedef struct {
    int nFluids;
    double gasGamma;
    bool updatePressure;
    double nuBase[5*4/2];
  } FiveMomentFrictionSrcData_t;

  void gkylFiveMomentFrictionSrcForwardEuler(
    const FiveMomentFrictionSrcData_t *sd,
    const FiveMomentFrictionFluidData_t *fd,
    const double dt,
    double **fPtrs);

  void gkylFiveMomentFrictionSrcTimeCentered(
    const FiveMomentFrictionSrcData_t *sd,
    const FiveMomentFrictionFluidData_t *fd,
    const double dt,
    double **fPtrs);

  void gkylFiveMomentFrictionSrcExact(
    const FiveMomentFrictionSrcData_t *sd,
    const FiveMomentFrictionFluidData_t *fd,
    const double dt,
    double **fPtrs);
}

#endif // GK_FIVE_MOMENT_FRICTION_SRC_H
