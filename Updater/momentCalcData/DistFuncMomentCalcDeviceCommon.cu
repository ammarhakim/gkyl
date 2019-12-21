/* -*- c++ -*- */
// Gkyl ------------------------------------------------------------------------
//
// Device kernels common to moment-calculating CUDA kernels.
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

// a#include "DistFuncMomentCalcDeviceCommon.h"

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
