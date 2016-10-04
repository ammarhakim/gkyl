// Gkyl ------------------------------------------------------------------------
//
// Euler RP in C
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_EULER_EQN_H
#define GK_EULER_EQN_H

extern "C" {
    /* Euler data */
    typedef struct { double gasGamma; } EulerEqn_t;

    /**
     * Reimann solver: `delta` is the vector we wish to split,
     * `ql`/`qr` the left/right states. On output, `waves` and `s`
     * contain the waves and speeds. waves is a mwave X meqn
     * matrix. See LeVeque's book for explanations.
     */
    void Euler_rp(EulerEqn_t *e, int dir, double *delta, double *ql, double *qr, double *waves, double *s);
}

#endif // GK_EULER_EQN_H
