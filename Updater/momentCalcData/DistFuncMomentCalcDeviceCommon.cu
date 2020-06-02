/* -*- c++ -*- */
// Gkyl ------------------------------------------------------------------------
//
// Device kernels common to moment-calculating CUDA kernels.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------
#include "DistFuncMomentCalcModDecl.h"
#include "DistFuncMomentCalcTmpl.h"
#include "DistFuncMomentCalcDeviceWrappers.h"

__inline__ __device__ void warpReduceComponentsSum(double *vals, int nComps) {
  // Perform 'nComps' independent (sum) reductions across a warp,
  // one for each component in 'vals'.

  for (unsigned int k = 0; k < nComps; k++) {
    for (int offset = warpSize/2; offset > 0; offset /= 2) {
      vals[0+k] += __shfl_down_sync(0xffffffff, vals[0+k], offset, warpSize);
    }
  }
}

__inline__ __device__ void blockReduceComponentsSum(double *vals, int nComps) {
  // Perform 'nComps' independent (sum) reductions across a block,
  // one for each component in 'nComps'.

  extern __shared__ double warpSum[]; // Stores partial sums.
  int lane   = threadIdx.x % warpSize;
  int warpID = threadIdx.x / warpSize;

  warpReduceComponentsSum(vals, nComps);            // Each warp performs partial reduction.

  if (blockDim.x > warpSize) {    // If statement needed in case number of velocity cells<warpSize.
    if (lane==0) {
      // Write reduced value to shared memory.
      for (unsigned int k = 0; k < nComps; k++) {
        warpSum[warpID*nComps+k] = vals[k];
      }
    }
  
    __syncthreads();                     // Wait for all partial reductions.
  
    // Read from shared memory (only for by the first warp).
    for (unsigned int k = 0; k < nComps; k++) {
      vals[k] = (threadIdx.x < blockDim.x / warpSize) ? warpSum[lane*nComps+k] : 0;
    }
  
    if (warpID==0) warpReduceComponentsSum(vals, nComps); // Final reduce within first warp.
  }

}

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
  double *cellSize    = &grid->dx[0];
  double *localSumPtr = &localSum[0];

  MomentModDecl<CDIM, VDIM, POLYORDER, BASISTYPE>::calcM0(cellxc, cellSize, distF, localSumPtr);

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
  double *cellSize    = &grid->dx[0];
  double *localSumPtr = &localSum[0];

  MomentModDecl<CDIM, VDIM, POLYORDER, BASISTYPE>::calcM1i(cellxc, cellSize, distF, localSumPtr);

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
