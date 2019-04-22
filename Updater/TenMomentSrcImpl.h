// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for ten-moment source terms
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_TEN_MOMENT_SRC_H
#define GK_TEN_MOMENT_SRC_H

#include <stdint.h>
#include <MomentSrcCommon.h>

extern "C" {
    /* Method to update fluids and flield using explicit RK3 method */
    void gkylTenMomentSrcRk3(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em);
    /* Method to update fluids and flield using time-centered implicit method */
    void gkylTenMomentSrcTimeCentered(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm, double *sigma);
    void gkylTenMomentSrcTimeCenteredDirect2(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma);
    void gkylTenMomentSrcTimeCenteredDirect(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma);
    void gkylTenMomentSrcExact(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm, double *sigma);
}

#endif // GK_TEN_MOMENT_SRC_H
