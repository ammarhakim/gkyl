/* -*- c++ -*- */
// Gkyl ------------------------------------------------------------------------
//
// CUDA wrappers for moment-calculating kernels.
// 
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include "DistFuncMomentCalcDevice.h"
#include "DistFuncMomentCalcModDeclDevice.cu"

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
      warpSum[warpID*2+k] = vals[k];
    }
  }

  __syncthreads();                     // Wait for all partial reductions.

  // Read from shared memory (only for by the first warp).
  for (unsigned int k = 0; k < nComps; k++) {
    vals[k] = (threadIdx.x < blockDim.x / warpSize) ? warpSum[lane*nComps+k] : 0;
  }

  if (warpID==0) warpReduceComponentsSum(vals, nComps); // Final reduce within first warp.

}

__host__ __device__ void MomentCalc1x1vSer_M0_P1(const double *w, const double *dxv, const double *f, double *out)
{
  const double volFact = dxv[1]/2;
  out[0] += 1.414213562373095*f[0]*volFact;
  out[1] += 1.414213562373095*f[1]*volFact;
}

__global__ void d_calcMom1x1vSer_M0_P1(int *nCells, double *w, double *dxv, double *fIn, double *out) {
  // Calculate the zeroth moment of the distribution function. We will first assign
  // whole configuration-space cells to a single block. Then one must perform a reduction
  // across a block for each conf-space basis coefficient.
  // Index of the current phase-space cell.
  unsigned int phaseGridIdx = blockIdx.x*blockDim.x + threadIdx.x;

  // Configuration and velocity space indexes.
  unsigned int confIdx = phaseGridIdx/nCells[1];
  unsigned int velIdx  = phaseGridIdx-confIdx*nCells[1];

  // Index of the first phase-space memory address to access.
  unsigned int phaseFldIdx = phaseGridIdx*4;

  double localSum[2];
  for (unsigned int k = 0; k < 2; k++) {
    localSum[k] = 0.0;  // Need to zero this out because kernel below increments.
  }

  // Pointers to quantities expected by the moment kernel.
  double *distF       = &fIn[phaseFldIdx];
  double *cellCenter  = &w[phaseGridIdx*2];
  double *cellSize    = &dxv[phaseGridIdx*2];
  double *localSumPtr = &localSum[0];

  MomentCalc1x1vSer_M0_P1(cellCenter, cellSize, distF, localSumPtr);

  blockReduceComponentsSum(localSumPtr, 2);
  if (threadIdx.x==0) {
    out[confIdx*2]   = localSumPtr[0];
    out[confIdx*2+1] = localSumPtr[1];
  }
}

void calcMom1x1vSer_M0_P1(int numBlocks, int numThreads, int *nCells, double *w, double *dxv, double *fIn, double *out)
{
  d_calcMom1x1vSer_M0_P1<<<numBlocks, numThreads, nCells[0]*2*sizeof(double)>>>(nCells, w, dxv, fIn, out);
}
