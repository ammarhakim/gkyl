// Gkyl ------------------------------------------------------------------------
//
// Header to C functions calling CUDA kernels that compute moments of the
// distribution function. 
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------
//
#include <GkylRectCart.h>
#include <GkylRange.h>
#include <GkylCudaFuncs.h>

#ifndef DIST_FUNC_MOMENT_CALC_DEVICE_WRAPPERS_H 
#define DIST_FUNC_MOMENT_CALC_DEVICE_WRAPPERS_H 

extern "C" { 

void cuda_MomentCalc1x1vSer_M0_P1(GkylRectCart_t *grid, GkylRange_t *pRange, GkylRange_t *cRange, GkDeviceProp *prop, int numBlocks, int numThreads, const double *fIn, double *out); 

} 
#endif 


