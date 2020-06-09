// Gkyl ------------------------------------------------------------------------
//
// CUDA back-end for CartField
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <GkylCudaFuncs.h>
#include <GkylCartField.h>
#include <GkylCudaReduce.h>

#ifndef GKYL_CART_FIELD_DEVICE_H
#define GKYL_CART_FIELD_DEVICE_H

extern "C" {
    // s: start index. sv: number of values to copy
    void gkylCartFieldDeviceAccumulate(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceAssign(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceScale(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, double *out);
    void gkylCartFieldDeviceAbs(int numBlocks, int numThreads, unsigned s, unsigned nv, double *out);

    // copy component data from/to field
    void gkylCopyFromFieldDevice(int numBlocks, int numThreads, double *data, double *f, unsigned numComponents, unsigned c);
    void gkylCopyToFieldDevice(int numBlocks, int numThreads, double *f, double *data, unsigned numComponents, unsigned c);

    // Assign all elements to specified value.
    void gkylCartFieldDeviceAssignAll(int numBlocks, int numThreads, unsigned s, unsigned nv, double val, double *out);

    // Reduction to a single value.
    void gkylCartFieldDeviceReduce(baseReduceOp *redOp, int numCellsTot, int numBlocks, int numThreads, int maxBlocks, int maxThreads,
       GkDeviceProp *prop, GkylCartField_t *fIn, double *blockOut, double *intermediate, double *out);
}

#endif // GKYL_CART_FIELD_DEVICE_H
