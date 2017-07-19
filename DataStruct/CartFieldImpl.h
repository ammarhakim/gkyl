// Gkyl ------------------------------------------------------------------------
//
// C++ back-end for CartField
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_CART_FIELD_H
#define GK_CART_FIELD_H

extern "C" {
    void gkylCartFieldAccumulate(unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldAssign(unsigned nv, double fact, const double *inp, double *out);
}

#endif // GK_CART_FIELD_H
