/* -*- c++ -*- */
// Gkyl ------------------------------------------------------------------------
//
// CUDA wrappers for moment-calculating kernels.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include "DistFuncMomentCalcDeviceWrappers.h"
#include "DistFuncMomentCalcDeviceCommon.cu"

void cuda_MomentCalc1x1vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 1, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}

void cuda_MomentCalc1x1vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 1, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}

void cuda_MomentCalc1x1vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 1, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}

void cuda_MomentCalc1x1vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 1, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}

