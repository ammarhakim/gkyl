// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for ten-moment gradient-based closure for heat flux
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_TEN_MOMENT_GRAD_H
#define GK_TEN_MOMENT_GRAD_H

#include <stdint.h>

/* Helper indices for indexing the temperature tensor and symmetrized heat flux tensor */
static const unsigned T11 = 0;
static const unsigned T12 = 1;
static const unsigned T13 = 2;
static const unsigned T22 = 3;
static const unsigned T23 = 4;
static const unsigned T33 = 5;

static const unsigned Q111 = 0;
static const unsigned Q112 = 1;
static const unsigned Q113 = 2;
static const unsigned Q122 = 3;
static const unsigned Q123 = 4;
static const unsigned Q133 = 5;
static const unsigned Q222 = 6;
static const unsigned Q223 = 7;
static const unsigned Q233 = 8;
static const unsigned Q333 = 9;

extern "C" {
   /* Helper function to obtain symmetrized heat flux */  
   void gkylTenMomentHeatFlux(const double alpha, const double* dT1, const double* dT2, const double* dT3, const double* f, double* q);
   /* Helper function to accumulate divergence of the heat flux tensor onto corresponding components of the stress tensor */
   void gkylTenMomentAccumulateGradClosure(const double dt, const double* divQ1, const double* divQ2, const double* divQ3, double* f);
   /* Compute the gradient of the temperature tensor in direction dir */
   void gkylTenMomentGradT(const int dir, const double* dxv, const double* fL, const double* fR, double* dT);
   /* Compute the components of the divergence of the heat flux tensor in X, Y, and Z */
   void gkylTenMomentDivQX(const double* dxv, const double* qL, const double* qR, double* divQ);
   void gkylTenMomentDivQY(const double* dxv, const double* qL, const double* qR, double* divQ);
   void gkylTenMomentDivQZ(const double* dxv, const double* qL, const double* qR, double* divQ);
}

#endif // GK_TEN_MOMENT_GRAD_H
