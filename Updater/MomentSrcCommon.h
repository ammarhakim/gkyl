// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for five-moment source terms
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_MOMENT_SRC_H
#define GK_MOMENT_SRC_H

#include <stdint.h>

extern "C" {
    typedef struct {
        double charge, mass; /* Charge and mass */
        bool evolve;

        double qbym;
    } FluidData_t;

    typedef struct {
        int8_t nFluids; /* Number of fluids */
        double epsilon0; /* Permittivity of free space */
        double chi_e, chi_m; /* Propagation speed factor for electric field error potential */
        int8_t gravityDir; /* Direction of gravity force */
        double gravity; /* Gravitational acceleration */
        bool hasStatic, hasPressure; /* Flag to indicate if there is: static EB field, pressure */
        int8_t linSolType; /* Flag to indicate linear solver type for implicit method */
    } MomentSrcData_t;

    /* Method to update fluids and flield using explicit RK3 method */
    void gkylMomentSrcRk3(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em);
    /* Method to update fluids and flield using time-centered implicit method */
    void gkylMomentSrcTimeCentered(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm);
    /* Method to update fluids and flield using time-centered implicit method in the exact form */
    void gkylMomentSrcTimeCenteredDirect(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm);
    /* Method to update fluids and flield using time-centered implicit method in the exact form; the older version*/
    void gkylMomentSrcTimeCenteredDirect2(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **f, double *em, double *staticEm);
    /* Method to update fluids and flield using the exact ODE solutions */
    void gkylMomentSrcExact(MomentSrcData_t *sd, FluidData_t *fd, double dt, double **ff, double *em, double *staticEm);
}

#endif // GK_MOMENT_SRC_H
