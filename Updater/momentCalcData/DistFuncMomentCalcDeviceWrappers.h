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
#include <GkylCartField.h>

#ifndef DIST_FUNC_MOMENT_CALC_DEVICE_WRAPPERS_H 
#define DIST_FUNC_MOMENT_CALC_DEVICE_WRAPPERS_H 

extern "C" { 

void cuda_MomentCalc1x1vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 
void cuda_MomentCalc1x1vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out); 

} 
#endif 


