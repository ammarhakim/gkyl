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
#include "MomentTmplModDecl.h"

template <int CDIM, int VDIM, int POLYORDER, int BASISTYPE>
__global__ void d_calcM0(GkylCartField_t *fIn, GkylCartField_t *out) {
  // In computing moments we will first assign whole configuration-space cells to a single block.
  // Then one must perform a reduction across a block for each conf-space basis coefficient.

  // Index of the first phase-space memory address to access.
  unsigned int linearIdx       = blockIdx.x*blockDim.x + threadIdx.x;
  GkylRange_t *localPhaseRange = fIn->localRange;
  Gkyl::GenIndexer localPhaseIdxr(localPhaseRange);
  Gkyl::GenIndexer fIdxr = fIn->genIndexer();

  int phaseIdx[CDIM+VDIM];
  localPhaseIdxr.invIndex(linearIdx, phaseIdx);
  int phaseLinIdx = fIdxr.index(phaseIdx);

  double localSum[Gkyl::BasisCount<CDIM+VDIM, POLYORDER, BASISTYPE>::numBasis()];
  unsigned int numComponents = out->numComponents;
  for (unsigned int k = 0; k < numComponents; k++) {
    localSum[k] = 0.0;
  }

  // Pointers to quantities expected by the moment kernel.
  const double *distF  = fIn->getDataPtrAt(phaseLinIdx);
  GkylRectCart_t *grid = fIn->grid;
  double cellxc[CDIM+VDIM];
  grid->cellCenter(phaseIdx,cellxc);
  double *cellCenter  = &cellxc[0];
  double *cellSize    = &grid->dx[0];
  double *localSumPtr = &localSum[0];

  MomentModDecl<CDIM, VDIM, POLYORDER, BASISTYPE>::calcM0(cellCenter, cellSize, distF, localSumPtr);
  
  blockReduceComponentsSum(localSumPtr, numComponents);

  // Configuration space indexes.
  GkylRange_t *localConfRange = out->localRange;
  Gkyl::GenIndexer localConfIdxr(localConfRange);
  Gkyl::GenIndexer outIdxr = out->genIndexer();
  int confIdx[CDIM];
  for (unsigned int k = 0; k < CDIM; k++) {
    confIdx[k] = phaseIdx[k];
  }
  int confLinIdx = outIdxr.index(confIdx);
  double *mom    = out->getDataPtrAt(confLinIdx);

  if (threadIdx.x==0) {
    for (unsigned int k = 0; k < numComponents; k++) {
      mom[k] = localSumPtr[k];
    }
  }

}

template <int CDIM, int VDIM, int POLYORDER, int BASISTYPE>
__global__ void d_calcM1i(GkylCartField_t *fIn, GkylCartField_t *out) {
  // In computing moments we will first assign whole configuration-space cells to a single block.
  // Then one must perform a reduction across a block for each conf-space basis coefficient.

  // Index of the first phase-space memory address to access.
  unsigned int linearIdx       = blockIdx.x*blockDim.x + threadIdx.x;
  GkylRange_t *localPhaseRange = fIn->localRange;
  Gkyl::GenIndexer localPhaseIdxr(localPhaseRange);
  Gkyl::GenIndexer fIdxr = fIn->genIndexer();

  int phaseIdx[CDIM+VDIM];
  localPhaseIdxr.invIndex(linearIdx, phaseIdx);
  int phaseLinIdx = fIdxr.index(phaseIdx);

  double localSum[Gkyl::BasisCount<CDIM+VDIM, POLYORDER, BASISTYPE>::numBasis()];
  unsigned int numComponents = out->numComponents;
  for (unsigned int k = 0; k < numComponents; k++) {
    localSum[k]=0.0;
  }

  // Pointers to quantities expected by the moment kernel.
  const double *distF  = fIn->getDataPtrAt(phaseLinIdx);
  GkylRectCart_t *grid = fIn->grid;
  double cellxc[CDIM+VDIM];
  grid->cellCenter(phaseIdx,cellxc);
  double *cellCenter  = &cellxc[0];
  double *cellSize    = &grid->dx[0];
  double *localSumPtr = &localSum[0];

  MomentModDecl<CDIM, VDIM, POLYORDER, BASISTYPE>::calcM1i(cellCenter, cellSize, distF, localSumPtr);
  
  blockReduceComponentsSum(localSumPtr, numComponents);

  // Configuration space indexes.
  GkylRange_t *localConfRange = out->localRange;
  Gkyl::GenIndexer localConfIdxr(localConfRange);
  Gkyl::GenIndexer outIdxr = out->genIndexer();
  int confIdx[CDIM];
  for (unsigned int k = 0; k < CDIM; k++) {
    confIdx[k] = phaseIdx[k];
  }
  int confLinIdx = outIdxr.index(confIdx);
  double *mom    = out->getDataPtrAt(confLinIdx);

  if (threadIdx.x==0) {
    for (unsigned int k = 0; k < numComponents; k++) {
      mom[k] = localSumPtr[k];
    }
  }

}

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

