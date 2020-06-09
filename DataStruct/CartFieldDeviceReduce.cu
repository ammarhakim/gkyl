// Gkyl ------------------------------------------------------------------------
//
// Functions to compute reductions of a CartField in GPU (Cuda).
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <CartFieldDeviceImpl.h>
#include <GkylRange.h>

extern "C" bool isPow2(unsigned int x);

// This version adds multiple elements per thread sequentially. This reduces
// the overall cost of the algorithm while keeping the work complexity O(n) and
// the step complexity O(log n). (Brent's Theorem optimization)
// Note, this kernel needs a minimum of 64*sizeof(T) bytes of shared memory.
// In other words if blockSize <= 32, allocate 64*sizeof(T) bytes.
// If blockSize > 32, allocate blockSize*sizeof(T) bytes.
template <unsigned int BLOCKSIZE, bool nIsPow2>
__global__ void d_reduceCartField(baseReduceOp *redOpIn, GkylCartField_t *fIn, double *redPerBlock) {
  // Handle to thread block group.
  cg::thread_block cgThreadBlock = cg::this_thread_block();
  extern __shared__ double sdata[];  // Stores partial reductions.

  // Perform first level of reduction, reading from global memory, writing to shared memory.
  unsigned int tID       = threadIdx.x;
  unsigned int linearIdx = blockIdx.x * BLOCKSIZE * 2 + threadIdx.x;
  unsigned int gridSize  = BLOCKSIZE * 2 * gridDim.x;

  double myReduc = redOpIn->initValue;

  GkylRange_t *localRange  = fIn->localRange;
  Gkyl::GenIndexer localIdxr(localRange);
  Gkyl::GenIndexer fIdxr   = fIn->genIndexer();
  unsigned int numCellsTot = localRange->volume();

  // We reduce multiple elements per thread.  The number is determined by the
  // number of active thread blocks (via gridDim). More blocks will result
  // in a larger gridSize and therefore fewer elements per thread.
  while (linearIdx < numCellsTot) {
    int idx[6];  // Should be CDIM+VDIM, but onerous to template just for this.
    localIdxr.invIndex(linearIdx, idx);
    int linIdx        = fIdxr.index(idx);
    const double *fld = fIn->getDataPtrAt(linIdx);

    myReduc = redOpIn->reduce(myReduc, fld[0]);

    // Ensure we don't read out of bounds (optimized away for powerOf2 sized arrays).
    unsigned int newLinearIdx = linearIdx+BLOCKSIZE;
    if (nIsPow2 || newLinearIdx<numCellsTot) {
      localIdxr.invIndex(newLinearIdx, idx);
      linIdx = fIdxr.index(idx);
      fld    = fIn->getDataPtrAt(linIdx);

      myReduc = redOpIn->reduce(myReduc, fld[0]);
    }

    linearIdx += gridSize;
  }

  // Each thread puts its local reduction into shared memory.
  sdata[tID] = myReduc;
  cg::sync(cgThreadBlock);

  // Do reduction in shared mem.
  if ((BLOCKSIZE >= 512) && (tID < 256)) {
    sdata[tID] = myReduc = redOpIn->reduce(myReduc, sdata[tID + 256]);
  }

  cg::sync(cgThreadBlock);

  if ((BLOCKSIZE >= 256) && (tID < 128)) {
    sdata[tID] = myReduc = redOpIn->reduce(myReduc, sdata[tID + 128]);
  }

  cg::sync(cgThreadBlock);

  if ((BLOCKSIZE >= 128) && (tID < 64)) {
    sdata[tID] = myReduc = redOpIn->reduce(myReduc, sdata[tID + 64]);
  }

  cg::sync(cgThreadBlock);

  cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cgThreadBlock);

  if (cgThreadBlock.thread_rank() < 32) {
    // Fetch final intermediate reduction from 2nd warp.
    if (BLOCKSIZE >= 64) myReduc = redOpIn->reduce(myReduc, sdata[tID + 32]);
    // Reduce final warp using shuffle.
    for (int offset = tile32.size()/2; offset > 0; offset /= 2) {
      double shflMax = tile32.shfl_down(myReduc, offset);
      myReduc = redOpIn->reduce(myReduc, shflMax);
    }
  }

  // Write result for this block to global mem.
  if (cgThreadBlock.thread_rank() == 0) { redPerBlock[blockIdx.x] = myReduc; }
}

void reduceCartField(baseReduceOp *opIn, int numCellsTot, int blocks, int threads, GkylCartField_t *fIn, double *blockRed) {
  // Launch the device kernel that reduces a CartField to an array,
  // each element of the array being the reduction performed by a block.

  // When there is only one warp per block, we need to allocate two warps
  // worth of shared memory so that we don't index shared memory out of bounds
  int smemSize = (threads <= 32) ? 2 * threads * sizeof(double) : threads * sizeof(double);

  if (isPow2(numCellsTot)) {
    switch (threads) {
      case 512:
        d_reduceCartField<512,true><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 256:
        d_reduceCartField<256,true><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 128:
        d_reduceCartField<128,true><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 64:
        d_reduceCartField<64,true><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 32:
        d_reduceCartField<32,true><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 16:
        d_reduceCartField<16,true><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 8:
        d_reduceCartField<8,true><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 4:
        d_reduceCartField<4,true><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 2:
        d_reduceCartField<2,true><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 1:
        d_reduceCartField<1,true><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
    }
  } else {
    switch (threads) {
      case 512:
        d_reduceCartField<512,false><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 256:
        d_reduceCartField<256,false><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 128:
        d_reduceCartField<128,false><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 64:
        d_reduceCartField<64,false><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 32:
        d_reduceCartField<32,false><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 16:
        d_reduceCartField<16,false><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 8:
        d_reduceCartField<8,false><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 4:
        d_reduceCartField<4,false><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 2:
        d_reduceCartField<2,false><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
      case 1:
        d_reduceCartField<1,false><<<blocks, threads, smemSize>>>(opIn, fIn, blockRed);
        break;
    }
  }
}

void gkylCartFieldDeviceReduce(baseReduceOp *redOp, int numCellsTot, int numBlocks, int numThreads, int maxBlocks, int maxThreads,
                  GkDeviceProp *prop, GkylCartField_t *fIn, double *blockOut, double *intermediate, double *out) {
  // Reduce the CartField 'fIn' (type double) according to the operation 'redOp'
  // and place it in the device-memory variable 'out'.
  // This function follows 'reduce6' (using Cooperative Groups) in cuda-samples:
  //   https://github.com/NVIDIA/cuda-samples/tree/master/Samples/reduction
  // The algorithm uses two other temporary variables: 'blockOut' and 
  // 'intermediate' have size=numBlocks and were allocated already.


  // Call the kernel that reduces a CartField (fIn) into a device array (blockOut)
  // which contains the reduction performed by each block.
  reduceCartField(redOp, numCellsTot, numBlocks, numThreads, fIn, blockOut);

  // Reduce partial block reductions on GPU.
  int newNum = numBlocks;
  while (newNum > 1) {
    int threads = 0, blocks = 0;

    reductionBlocksAndThreads(prop, newNum, maxBlocks, maxThreads, blocks, threads);

    checkCudaErrors(cudaMemcpy(intermediate, blockOut, newNum * sizeof(double), cudaMemcpyDeviceToDevice));

    reduceDeviceArray(redOp, newNum, blocks, threads, intermediate, blockOut);

    newNum = (newNum + (threads*2-1))/(threads*2);
  }

  cudaDeviceSynchronize();
  // Copy final reduction to output variable.
  checkCudaErrors(cudaMemcpy(out, blockOut, sizeof(double), cudaMemcpyDeviceToDevice));
}

