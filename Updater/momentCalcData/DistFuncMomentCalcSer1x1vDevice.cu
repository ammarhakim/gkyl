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


 __global__ void d_calcMom1x1vSer_M0_P1(GkylCartField_t *fIn, GkylCartField_t *out) {
  // In computing moments we will first assign whole configuration-space cells to a single block.
  // Then one must perform a reduction across a block for each conf-space basis coefficient.

  // Index of the first phase-space memory address to access.
  unsigned int linearIdx       = blockIdx.x*blockDim.x + threadIdx.x;
  GkylRange_t *localPhaseRange = fIn->localRange;
  Gkyl::GenIndexer localPhaseIdxr(localPhaseRange);
  Gkyl::GenIndexer fIdxr = fIn->genIndexer();

  int phaseIdx[2];
  localPhaseIdxr.invIndex(linearIdx, phaseIdx);
  int phaseLinIdx = fIdxr.index(phaseIdx);

  double localSum[2];
  for (unsigned int k = 0; k < 2; k++) {
    localSum[k]=0.0;
  }

  // Pointers to quantities expected by the moment kernel.
  const double *distF  = fIn->getDataPtrAt(phaseLinIdx);
  GkylRectCart_t *grid = fIn->grid;
  double cellxc[2];
  grid->cellCenter(phaseIdx,cellxc);
  double *cellCenter  = &cellxc[0];
  double *cellSize    = &grid->dx[0];
  double *localSumPtr = &localSum[0];

  MomentCalc1x1vSer_M0_P1(cellCenter, cellSize, distF, localSumPtr);

  blockReduceComponentsSum(localSumPtr, 2);

  // Configuration space indexes.
  GkylRange_t *localConfRange = out->localRange;
  Gkyl::GenIndexer localConfIdxr(localConfRange);
  Gkyl::GenIndexer outIdxr = out->genIndexer();
  int confIdx[1];
  for (unsigned int k = 0; k < 1; k++) {
    confIdx[k] = phaseIdx[k];
  }
  int confLinIdx = outIdxr.index(phaseIdx);
  double *mom    = out->getDataPtrAt(confLinIdx);
  unsigned int numComponents = out->numComponents;

  if (threadIdx.x==0) {
    for (unsigned int k = 0; k < numComponents; k++) {
      mom[confLinIdx+k] = localSumPtr[k];
    }
  }

}


 void cuda_MomentCalc1x1vSer_M0_P1(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcMom1x1vSer_M0_P1<<<numBlocks, numThreads, 2*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}


 __global__ void d_calcMom1x1vSer_M0_P2(GkylCartField_t *fIn, GkylCartField_t *out) {
  // In computing moments we will first assign whole configuration-space cells to a single block.
  // Then one must perform a reduction across a block for each conf-space basis coefficient.

  // Index of the first phase-space memory address to access.
  unsigned int linearIdx       = blockIdx.x*blockDim.x + threadIdx.x;
  GkylRange_t *localPhaseRange = fIn->localRange;
  Gkyl::GenIndexer localPhaseIdxr(localPhaseRange);
  Gkyl::GenIndexer fIdxr = fIn->genIndexer();

  int phaseIdx[2];
  localPhaseIdxr.invIndex(linearIdx, phaseIdx);
  int phaseLinIdx = fIdxr.index(phaseIdx);

  double localSum[3];
  for (unsigned int k = 0; k < 3; k++) {
    localSum[k]=0.0;
  }

  // Pointers to quantities expected by the moment kernel.
  const double *distF  = fIn->getDataPtrAt(phaseLinIdx);
  GkylRectCart_t *grid = fIn->grid;
  double cellxc[2];
  grid->cellCenter(phaseIdx,cellxc);
  double *cellCenter  = &cellxc[0];
  double *cellSize    = &grid->dx[0];
  double *localSumPtr = &localSum[0];

  MomentCalc1x1vSer_M0_P2(cellCenter, cellSize, distF, localSumPtr);

  blockReduceComponentsSum(localSumPtr, 3);

  // Configuration space indexes.
  GkylRange_t *localConfRange = out->localRange;
  Gkyl::GenIndexer localConfIdxr(localConfRange);
  Gkyl::GenIndexer outIdxr = out->genIndexer();
  int confIdx[1];
  for (unsigned int k = 0; k < 1; k++) {
    confIdx[k] = phaseIdx[k];
  }
  int confLinIdx = outIdxr.index(phaseIdx);
  double *mom    = out->getDataPtrAt(confLinIdx);
  unsigned int numComponents = out->numComponents;

  if (threadIdx.x==0) {
    for (unsigned int k = 0; k < numComponents; k++) {
      mom[confLinIdx+k] = localSumPtr[k];
    }
  }

}


 void cuda_MomentCalc1x1vSer_M0_P2(GkDeviceProp *prop, int numBlocks, int numThreads, GkylCartField_t *fIn, GkylCartField_t *out) {
  int warpSize = prop->warpSize;

  d_calcMom1x1vSer_M0_P2<<<numBlocks, numThreads, 3*(numThreads/warpSize)*sizeof(double)>>>(fIn, out);
}
