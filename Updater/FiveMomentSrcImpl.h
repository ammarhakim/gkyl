// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for five-moment source terms
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_FIVE_MOMENT_SRC_H
#define GK_FIVE_MOMENT_SRC_H

#include <stdint.h>
#include <MomentSrcCommon.h>

extern "C" {
    /* Method to update fluids and flield using explicit RK3 method */
    void gkylFiveMomentSrcRk3(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em);
    /* Method to update fluids and flield using time-centered implicit method */
    void gkylFiveMomentSrcTimeCentered(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm);
    /* Method to update fluids and flield using time-centered implicit method in the exact form */
    void gkylFiveMomentSrcTimeCenteredDirect(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm);
    /* Method to update fluids and flield using time-centered implicit method in the exact form; the older version*/
    void gkylFiveMomentSrcTimeCenteredDirect2(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm);
    /* Method to update fluids and flield using the exact ODE solutions */
    void gkylFiveMomentSrcExact(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm);
}

#endif // GK_FIVE_MOMENT_SRC_H
