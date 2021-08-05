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


 void cuda_MomentCalc1x2vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 2, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 2, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 2, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 2, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 2, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 2, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 2, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 2, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 2, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vSer_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 2, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 3, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 9*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 3, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 12*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 3, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 9*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 3, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 12*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 3, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 9*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 3, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 12*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 3, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 9*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 3, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 12*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 3, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 9*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x3vSer_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 3, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 12*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<2, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<2, 2, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 16*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M0_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<2, 2, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<2, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<2, 2, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 16*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M1i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<2, 2, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<2, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<2, 2, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 16*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M2ij_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<2, 2, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<2, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<2, 2, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 16*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M2_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<2, 2, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<2, 2, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<2, 2, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 16*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vSer_M3i_P3(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<2, 2, 3, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x3vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<2, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 12*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x3vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<2, 3, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x3vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<2, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 12*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x3vSer_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<2, 3, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x3vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<2, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 12*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x3vSer_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<2, 3, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x3vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<2, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 12*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x3vSer_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<2, 3, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x3vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<2, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 12*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x3vSer_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<2, 3, 2, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc3x3vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<3, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc3x3vSer_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<3, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc3x3vSer_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<3, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc3x3vSer_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<3, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc3x3vSer_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<3, 3, 1, Gkyl::G_SERENDIPITY_C><<<numBlocks, numThreads, 24*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vTensor_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 1, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vTensor_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 1, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vTensor_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 1, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vTensor_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 1, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vTensor_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 1, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vTensor_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 1, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vTensor_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 1, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vTensor_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 1, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vTensor_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 1, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x1vTensor_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 1, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vTensor_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 2, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vTensor_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<1, 2, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vTensor_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 2, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vTensor_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<1, 2, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vTensor_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 2, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vTensor_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<1, 2, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vTensor_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 2, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vTensor_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<1, 2, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vTensor_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 2, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 4*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc1x2vTensor_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<1, 2, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 6*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vTensor_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<2, 2, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vTensor_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM0<2, 2, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 18*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vTensor_M1i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<2, 2, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vTensor_M1i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM1i<2, 2, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 18*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vTensor_M2ij_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<2, 2, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vTensor_M2ij_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2ij<2, 2, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 18*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vTensor_M2_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<2, 2, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vTensor_M2_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM2<2, 2, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 18*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vTensor_M3i_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<2, 2, 1, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 8*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 void cuda_MomentCalc2x2vTensor_M3i_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcM3i<2, 2, 2, Gkyl::G_TENSOR_C><<<numBlocks, numThreads, 18*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}

