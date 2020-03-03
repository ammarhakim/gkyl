// Gkyl ------------------------------------------------------------------------
//
// CUDA back-end for CartField
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_CART_FIELD_DEVICE_H
#define GK_CART_FIELD_DEVICE_H

extern "C" {
    // s: start index. sv: number of values to copy
    void gkylCartFieldDeviceAccumulate(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceAssign(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceScale(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, double *out);
    void gkylCartFieldDeviceAbs(int numBlocks, int numThreads, unsigned s, unsigned nv, double *out);

    // copy component data from/to field
    void gkylCopyFromFieldDevice(int numBlocks, int numThreads, double *data, double *f, unsigned numComponents, unsigned c);
    void gkylCopyToFieldDevice(int numBlocks, int numThreads, double *f, double *data, unsigned numComponents, unsigned c);

    // assign all elements to specified value
    void gkylCartFieldDeviceAssignAll(int numBlocks, int numThreads, unsigned s, unsigned nv, double val, double *out);
}

#endif // GK_CART_FIELD_DEVICE_H
