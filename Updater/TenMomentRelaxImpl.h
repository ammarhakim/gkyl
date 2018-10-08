// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for ten-moment relaxation
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_TEN_MOMENT_RELAX_H
#define GK_TEN_MOMENT_RELAX_H

#include <stdint.h>

extern "C" {
    typedef struct {
        double charge, mass; /* Charge and mass */
    } FluidData_t;

    typedef struct {
        double k;
        bool hasEm; /* Flag to indicate if there is: EB field */
        bool hasStatic; /* Flag to indicate if there is: static EB field */
    } TenMomentRelaxData_t;

    /* Use explicit RK3 method */
    void gkylTenMomentRelaxRk3(TenMomentRelaxData_t *sd, FluidData_t *fd, double dt, double *f, double *em);
    /* Use time-centered implicit method */
    void gkylTenMomentRelaxExplicit(TenMomentRelaxData_t *sd, FluidData_t *fd, double dt, double *f, double *em, double *staticEm);
}

#endif // GK_TEN_MOMENT_RELAX_H
