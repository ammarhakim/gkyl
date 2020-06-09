// Gkyl ------------------------------------------------------------------------
//
// Functions to compute reductions in GPU (Cuda).
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#pragma once

#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <float.h>

#include <GkylCudaFuncs.h>

#ifndef MIN
#define MIN(x, y) ((x < y) ? x : y)
#endif

#ifndef MAX
#define MAX(x, y) ((x > y) ? x : y)
#endif

#ifndef SUM
#define SUM(x,y) (x+y)
#endif

extern "C" {
    void reductionBlocksAndThreads(GkDeviceProp *prop, int numElements, int maxBlocks,
      int maxThreads, int &blocks, int &threads);
    bool isPow2(unsigned int x);
    unsigned int nextPow2(unsigned int x);
    void reduceDeviceArray(int opIn, int numElements, int blocks, int threads, double *d_dataIn, double *d_dataOut);
}

namespace Gkyl {

  enum class BinOp {
    binOpMin = 1,
    binOpMax = 2,
    binOpSum = 3
  };

  template <Gkyl::BinOp binOpType>
  __inline__ __device__ double binOp(double a, double b) {
    if (binOpType == BinOp::binOpMin) {
      return ((a < b) ? a : b);
    } else if (binOpType == BinOp::binOpMax) {
      return ((a > b) ? a : b);
    } else if (binOpType == BinOp::binOpSum) {
      return a+b;
    }
    return 0;
  }
}
