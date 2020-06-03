/* -*- c++ -*- */
// Gkyl ------------------------------------------------------------------------
//
// C wrappers to moment-calculating CUDA kernels.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include "DistFuncMomentCalcDeviceModDecl.h"
#include "DistFuncMomentCalcDeviceCommon.h"

 void cuda_MomentCalc1x1vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 1, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 1, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 1, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 1, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 1, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 1, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 1, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 1, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 1, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 1, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 1, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 1, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 1, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 1, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vSer_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 1, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}

// void cuda_MomentCalc1x2vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
//  int warpSize = prop->warpSize;
//
//  d_calcM0<1, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
//}
// void cuda_MomentCalc1x2vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
//  int warpSize = prop->warpSize;
//
//  d_calcM1i<1, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
//}
// void cuda_MomentCalc1x2vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
//  int warpSize = prop->warpSize;
//
//  d_calcM2ij<1, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
//}
// void cuda_MomentCalc1x2vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
//  int warpSize = prop->warpSize;
//
//  d_calcM2<1, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
//}
// void cuda_MomentCalc1x2vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
//  int warpSize = prop->warpSize;
//
//  d_calcM3i<1, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
//}



