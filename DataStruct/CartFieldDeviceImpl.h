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
    void gkylCartFieldDeviceAccumulateOffset(int numBlocks, int numThreads, unsigned sInp, unsigned sOut, unsigned nCells, unsigned compStart, unsigned nCompInp, unsigned nCompOut, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceAssign(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, const double *inp, double *out);
    void gkylCartFieldDeviceScale(int numBlocks, int numThreads, unsigned s, unsigned nv, double fact, double *out);
    void gkylCartFieldDeviceAbs(int numBlocks, int numThreads, unsigned s, unsigned nv, double *out);

    // copy periodic boundary conditions when using only one GPU
    void gkylPeriodicCopy(int numBlocks, int numThreads, GkylRange_t *rangeSkin, GkylRange_t *rangeGhost, GkylCartField_t *f, unsigned numComponents);

    // Assign all elements to specified value.
    void gkylCartFieldDeviceAssignAll(int numBlocks, int numThreads, unsigned s, unsigned nv, double val, double *out);

    // Reduction to a single value.
    void gkylCartFieldDeviceReduce(baseReduceOp *redOp, int numCellsTot, int numComponents, int numBlocks, int numThreads, int maxBlocks, int maxThreads, GkDeviceProp *prop, GkylCartField_t *fIn, double *blockOut, double *intermediate, double *out);
}

#endif // GKYL_CART_FIELD_DEVICE_H
