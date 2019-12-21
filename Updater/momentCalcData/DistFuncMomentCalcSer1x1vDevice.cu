/* -*- c++ -*- */
// Gkyl ------------------------------------------------------------------------
//
// CUDA wrappers for moment-calculating kernels.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include "DistFuncMomentCalcDeviceWrappers.h"
#include "DistFuncMomentCalcModDecl.h"
#include "DistFuncMomentCalcDeviceCommon.cu"


 __global__ void d_calcMom1x1vSer_M0_P1(RectCart_t *grid, Range_t *range, double *fIn, double *out) {
  // In computing moments we will first assign whole configuration-space cells to a single block.
  // Then one must perform a reduction across a block for each conf-space basis coefficient.

  // Index of the first phase-space memory address to access.
  unsigned int phaseLinIdx = blockIdx.x*blockDim.x + threadIdx.x;
  GkylCuda::Indexer<2> phaseIndxr(range);
  int phaseMultiDimIdx[2];
  phaseIndxr.invIndex(phaseLinIdx, phaseMultiDimIdx);

  double localSum[2];
  for (unsigned int k = 0; k < 2; k++) {
    localSum[k]=0.0;
  }

  // Pointers to quantities expected by the moment kernel.
  double *distF       = &fIn[phaseLinIdx*4];
  double cellw[2];
  grid->cellCenter(phaseMultiDimIdx,cellw);
  double *cellCenter  = &cellw[0];
  double *cellSize    = &grid->dx[0];
  double *localSumPtr = &localSum[0];

  MomentCalc1x1vSer_M0_P1(cellCenter, cellSize, distF, localSumPtr);
  
  blockReduceComponentsSum(localSumPtr, 2);

  // Configuration space indexes.
  GkylCuda::Indexer<1> confIndxr(range);
  int confMultiDimIdx[1];
  confMultiDimIdx[0]=phaseMultiDimIdx[0];
  int confLinMemIdx = 2*(confIndxr.index(confMultiDimIdx[0]));

  if (threadIdx.x==0) {
    out[confLinMemIdx]=localSumPtr[0];
    out[confLinMemIdx+1]=localSumPtr[1];
  }

}

 void calcMom1x1vSer_M0_P1(RectCart_t *grid, Range_t *range, GkDeviceProp *prop, int numBlocks, int numThreads, double *fIn, double *out) {

  int warpSize = prop->warpSize;

  d_calcMom1x1vSer_M0_P1<<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(grid, range, fIn, out);
}

