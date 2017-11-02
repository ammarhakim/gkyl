// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for CartField
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_CART_FIELD_H
#define GK_CART_FIELD_H

extern "C" {
    // s: start index. sv: number of values to copy
    void gkylCartFieldAccumulate(unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldAssign(unsigned s, unsigned nv, double fact, const double *inp, double *out);
}

#endif // GK_CART_FIELD_H
