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
    void gkylCartFieldScale(unsigned s, unsigned nv, double fact, double *out);
    void gkylCartFieldScaleByCell(unsigned s, unsigned nv, unsigned ncomp, double *fact, double *out);
    void gkylCartFieldAbs(unsigned s, unsigned nv, double *out);

    // ncopy: number of components to copy (size of cInp and cOut arrays)
    // ncInp: number of components in input field. cInp: list of components to copy from
    // ncOut: number of components in output field. cOut: list of components to copy to
    void gkylCartFieldGenAccumulate(unsigned s, unsigned nv, double fact, int ncopy,
      int ncInp, const int *cInp, const double *inp,
      int ncOut, const int *cOut, double *out);

    // copy component data from/to field
    void gkylCopyFromField(double *data, double *f, unsigned numComponents, unsigned c);
    void gkylCopyToField(double *f, double *data, unsigned numComponents, unsigned c);

    // assign all elements to specified value
    void gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out);
}

#endif // GK_CART_FIELD_H
