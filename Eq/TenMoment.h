// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for ten-moment core functions
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef HAVE_TENMOM_RP_H
#define HAVE_TENMOM_RP_H

extern "C" {
    void gkylTenMomentRp(int dir, double delta[10], double ql[10], double qr[10], double *waves, double s[5]);
    void gkylTenMomentQFluct(double *waves, double *s, double *amdq, double *apdq);
}

#endif // HAVE_TENMOM_RP_H
