// Gkyl ------------------------------------------------------------------------
//
// Functions to compute reductions in GPU (Cuda).
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <cuda_runtime.h>
#include <cooperative_groups.h>

#include <GkylCudaFuncs.h>
#include <GkylCartField.h>
#include <GkylRange.h>

namespace cg = cooperative_groups;

extern "C"
{
  void getNumBlocksAndThreads(GkDeviceProp *prop, int numElements, int maxBlocks,
                              int maxThreads, int &blocks, int &threads);
  void cartFieldMax(int numCellsTot, int numBlocks, int numThreads, int maxBlocks, int maxThreads, GkDeviceProp *prop,
                    GkylCartField_t *fIn, double *blockOut, double *intermediate, double *out);
}

bool isPow2(unsigned int x) { return ((x & (x - 1)) == 0); }

unsigned int nextPow2(unsigned int x) {
  --x;
  x |= x >> 1;
  x |= x >> 2;
  x |= x >> 4;
  x |= x >> 8;
  x |= x >> 16;
  return ++x;
}

#ifndef MIN
#define MIN(x, y) ((x < y) ? x : y)
#endif

#ifndef MAX
#define MAX(x, y) ((x > y) ? x : y)
#endif

// Compute the number of threads and blocks to use for the given reduction
// kerne. We set threads/block to the minimum of maxThreads and n/2.
// We observe the maximum specified number of blocks, because
// each thread in the kernel can process a variable number of elements.
void getNumBlocksAndThreads(GkDeviceProp *prop, int numElements, int maxBlocks,
                            int maxThreads, int &blocks, int &threads) {

  threads = (numElements < maxThreads * 2) ? nextPow2((numElements + 1) / 2) : maxThreads;
  blocks  = (numElements + (threads * 2 - 1)) / (threads * 2);

  if ((float)threads * blocks >
      (float)(prop->maxGridSize)[0] * prop->maxThreadsPerBlock) {
    printf("n is too large, please choose a smaller number!\n");
  }

  if (blocks > (prop->maxGridSize)[0]) {
    printf("Grid size <%d> exceeds the device capability <%d>, set block size as %d (original %d)\n",
        blocks, (prop->maxGridSize)[0], threads * 2, threads);

    blocks  /= 2;
    threads *= 2;
  }

  blocks = MIN(maxBlocks, blocks);
}

// This version adds multiple elements per thread sequentially. This reduces
// the overall cost of the algorithm while keeping the work complexity O(n) and
// the step complexity O(log n). (Brent's Theorem optimization)
// Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
// In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
// If blockSize > 32, allocate blockSize*sizeof(T) bytes.
template <unsigned int BLOCKSIZE, bool nIsPow2>
__global__ void d_reduceCartField(GkylCartField_t *fIn, double *redPerBlock) {
  // Handle to thread block group.
  cg::thread_block cgThreadBlock = cg::this_thread_block();
  extern __shared__ double sdata[];  // Stores partial reductions.

  // Perform first level of reduction, reading from global memory, writing to shared memory.
  unsigned int tID       = threadIdx.x;
  unsigned int linearIdx = blockIdx.x * BLOCKSIZE * 2 + threadIdx.x;
  unsigned int gridSize  = BLOCKSIZE * 2 * gridDim.x;

  double myMax = 0;

  GkylRange_t *localRange  = fIn->localRange;
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr   = fIn->genIndexer();
  unsigned int numCellsTot = localRange->volume();

  // We reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridDim). More blocks will result
  // in a larger gridSize and therefore fewer elements per thread.
  while (linearIdx < numCellsTot) {
    int idx[2];  // Should be CDIM+VDIM, but onerous to template just for this.
    localIdxr.invIndex(linearIdx, idx);
    int linIdx        = fIdxr.index(idx);
    const double *fld = fIn->getDataPtrAt(linIdx);

    myMax = MAX(myMax, fld[0]);

    // Ensure we don't read out of bounds (optimized away for powerOf2 sized arrays).
    unsigned int newLinearIdx = linearIdx+BLOCKSIZE;
    if (nIsPow2 || newLinearIdx<numCellsTot) {
      localIdxr.invIndex(newLinearIdx, idx);
      linIdx = fIdxr.index(idx);
      fld    = fIn->getDataPtrAt(linIdx);

      myMax = MAX(myMax, fld[0]);
    }

    linearIdx += gridSize;
  }

  // Each thread puts its local max into shared memory.
  sdata[tID] = myMax;
  cg::sync(cgThreadBlock);

  // Do reduction in shared mem.
  if ((BLOCKSIZE >= 512) && (tID < 256)) {
    sdata[tID] = myMax = MAX(myMax, sdata[tID + 256]);
  }

  cg::sync(cgThreadBlock);

  if ((BLOCKSIZE >= 256) && (tID < 128)) {
    sdata[tID] = myMax = MAX(myMax, sdata[tID + 128]);
  }

  cg::sync(cgThreadBlock);

  if ((BLOCKSIZE >= 128) && (tID < 64)) {
    sdata[tID] = myMax = MAX(myMax, sdata[tID + 64]);
  }

  cg::sync(cgThreadBlock);

  cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cgThreadBlock);

  if (cgThreadBlock.thread_rank() < 32) {
    // Fetch final intermediate max from 2nd warp.
    if (BLOCKSIZE >= 64) myMax = MAX(myMax, sdata[tID + 32]);
    // Reduce final warp using shuffle.
    for (int offset = tile32.size()/2; offset > 0; offset /= 2) {
      double shflMax = tile32.shfl_down(myMax, offset);
      myMax = MAX(myMax, shflMax);
    }
  }

  // Write result for this block to global mem.
  if (cgThreadBlock.thread_rank() == 0) {
//    printf("max found = %14.12f\n",myMax);
    redPerBlock[blockIdx.x] = myMax;
  }
}

// This algorithm reduces multiple elements per thread sequentially. This reduces
// the overall cost of the algorithm while keeping the work complexity O(n) and
// the step complexity O(log n). (Brent's Theorem optimization)
// Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
// In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
// If blockSize > 32, allocate blockSize*sizeof(T) bytes.
template <unsigned int BLOCKSIZE, bool nIsPow2>
__global__ void d_reduceMax(double *dataIn, double *out, unsigned int nElements) {
  // Handle to thread block group.
  cg::thread_block cgThreadBlock = cg::this_thread_block();
  extern __shared__ double sdata[];  // Stores partial reductions.

  // Perform first level of reduction, reading from global memory, writing to shared memory.
  unsigned int tID       = threadIdx.x;
  unsigned int linearIdx = blockIdx.x * BLOCKSIZE * 2 + threadIdx.x;
  unsigned int gridSize  = BLOCKSIZE * 2 * gridDim.x;

  double myMax = 0;

  // We reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridDim).  More blocks will result
  // in a larger gridSize and therefore fewer elements per thread
  while (linearIdx < nElements) {
    myMax = MAX(myMax, dataIn[linearIdx]);

    // Ensure we don't read out of bounds (optimized away for powerOf2 sized arrays)/
    if (nIsPow2 || linearIdx+BLOCKSIZE<nElements) myMax = MAX(myMax, dataIn[linearIdx+BLOCKSIZE]);

    linearIdx += gridSize;
  }

  // Each thread puts its local reduction into shared memory.
  sdata[tID] = myMax;
  cg::sync(cgThreadBlock);

  // Do reduction in shared mem.
  if ((BLOCKSIZE >= 512) && (tID < 256)) {
    sdata[tID] = myMax = MAX(myMax, sdata[tID + 256]);
  }

  cg::sync(cgThreadBlock);

  if ((BLOCKSIZE >= 256) && (tID < 128)) {
    sdata[tID] = myMax = MAX(myMax, sdata[tID + 128]);
  }

  cg::sync(cgThreadBlock);

  if ((BLOCKSIZE >= 128) && (tID < 64)) {
    sdata[tID] = myMax = MAX(myMax, sdata[tID + 64]);
  }

  cg::sync(cgThreadBlock);

  cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cgThreadBlock);

  if (cgThreadBlock.thread_rank() < 32) {
    // Fetch final intermediate reduction from 2nd warp.
    if (BLOCKSIZE >= 64) myMax = MAX(myMax, sdata[tID + 32]);
    // Reduce final warp using shuffle.
    for (int offset = tile32.size() / 2; offset > 0; offset /= 2) {
      double shflMax = tile32.shfl_down(myMax, offset);
      myMax = MAX(myMax, shflMax);
    }
  }

  // Write result for this block to global mem.
  if (cgThreadBlock.thread_rank() == 0) {
    out[blockIdx.x] = myMax;
  }
}

void reduceCartField(int numCellsTot, int blocks, int threads, GkylCartField_t *fIn, double *blockRed) {
  // Launch the device kernel that reduces a CartField to an array,
  // each element of the array being the reduction performed by a block.

  // When there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = (threads <= 32) ? 2 * threads * sizeof(double) : threads * sizeof(double);

  if (isPow2(numCellsTot)) {
    switch (threads) {
      case 512:
        d_reduceCartField<512, true><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 256:
        d_reduceCartField<256, true><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 128:
        d_reduceCartField<128, true><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 64:
        d_reduceCartField<64, true><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 32:
        d_reduceCartField<32, true><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 16:
        d_reduceCartField<16, true><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 8:
        d_reduceCartField<8, true><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 4:
        d_reduceCartField<4, true><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 2:
        d_reduceCartField<2, true><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 1:
        d_reduceCartField<1, true><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
    }
  } else {
    switch (threads) {
      case 512:
        d_reduceCartField<512, false><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 256:
        d_reduceCartField<256, false><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 128:
        d_reduceCartField<128, false><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 64:
        d_reduceCartField<64, false><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 32:
        d_reduceCartField<32, false><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 16:
        d_reduceCartField<16, false><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 8:
        d_reduceCartField<8, false><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 4:
        d_reduceCartField<4, false><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 2:
        d_reduceCartField<2, false><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
      case 1:
        d_reduceCartField<1, false><<<blocks, threads, smemSize>>>(fIn, blockRed);
        break;
    }
  }
}

void reduceMax(int numElements, int blocks, int threads, double *d_dataIn, double *d_dataOut) {
  // Launch the device kernel that reduces a device array 'd_dataIn'
  // containing 'numElements' elements.

  // When there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = (threads <= 32) ? 2 * threads * sizeof(double) : threads * sizeof(double);

  if (isPow2(numElements)) {
    switch (threads) {
      case 512:
        d_reduceMax<512, true><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 256:
        d_reduceMax<256, true><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 128:
        d_reduceMax<128, true><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 64:
        d_reduceMax<64, true><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 32:
        d_reduceMax<32, true><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 16:
        d_reduceMax<16, true><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 8:
        d_reduceMax<8, true><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 4:
        d_reduceMax<4, true><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 2:
        d_reduceMax<2, true><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 1:
        d_reduceMax<1, true><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
    }
  } else {
    switch (threads) {
      case 512:
        d_reduceMax<512, false><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 256:
        d_reduceMax<256, false><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 128:
        d_reduceMax<128, false><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 64:
        d_reduceMax<64, false><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 32:
        d_reduceMax<32, false><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 16:
        d_reduceMax<16, false><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 8:
        d_reduceMax<8, false><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 4:
        d_reduceMax<4, false><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 2:
        d_reduceMax<2, false><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
      case 1:
        d_reduceMax<1, false><<<blocks, threads, smemSize>>>(d_dataIn, d_dataOut, numElements);
        break;
    }
  }
}

void cartFieldMax(int numCellsTot, int numBlocks, int numThreads, int maxBlocks, int maxThreads, GkDeviceProp *prop,
                  GkylCartField_t *fIn, double *blockOut, double *intermediate, double *out) {
  // Find the maximum in the CartField 'fIn' (type double) and place
  // it in the device-memory variable 'out'.
  // This function follows 'reduce6' (using Cooperative Groups) in cuda-samples:
  //   https://github.com/NVIDIA/cuda-samples/tree/master/Samples/reduction
  // The algorithm uses two other temporary variables: 'blockOut' and 
  // 'intermediate' have size=numBlocks and were allocated already.


  // Call the kernel that reduces a CartField (fIn) into a device array (blockOut)
  // which contains the reduction performed by each block.
  reduceCartField(numCellsTot, numBlocks, numThreads, fIn, blockOut);

  // Reduce partial block reductions on GPU.
  int newNum = numBlocks;
  while (newNum > 1) {
    int threads = 0, blocks = 0;

    getNumBlocksAndThreads(prop, newNum, maxBlocks, maxThreads, blocks, threads);

    checkCudaErrors(cudaMemcpy(intermediate, blockOut, newNum * sizeof(double), cudaMemcpyDeviceToDevice));

    reduceMax(newNum, blocks, threads, intermediate, blockOut);

    newNum = (newNum + (threads*2-1))/(threads*2);
  }

  cudaDeviceSynchronize();
  // Copy final reduction to output variable.
  checkCudaErrors(cudaMemcpy(out, blockOut, sizeof(double), cudaMemcpyDeviceToDevice));
}

