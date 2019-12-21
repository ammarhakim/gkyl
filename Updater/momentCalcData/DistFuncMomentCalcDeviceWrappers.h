// Gkyl ------------------------------------------------------------------------
//
// Header to C functions calling CUDA kernels that compute moments of the
// distribution function. 
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------
//
#include <RectCartDeviceImpl.h>
#include <RangeDeviceImpl.h>
#include <GkCudaFuncs.h>

#ifndef DIST_FUNC_MOMENT_CALC_DEVICE_WRAPPERS_H 
#define DIST_FUNC_MOMENT_CALC_DEVICE_WRAPPERS_H 

extern "C" { 

void calcMom1x1vSer_M0_P1(RectCart_t *grid, Range_t *range, GkDeviceProp *prop, int numBlocks, int numThreads, const double *fIn, double *out); 

} 
#endif 


