// Gkyl ------------------------------------------------------------------------
//
// Functions to compute reductions in GPU (Cuda).
//
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <cuda_runtime.h>
#include <cooperative_groups.h>
#include <float.h>

#include <GkylCudaFuncs.h>

namespace cg = cooperative_groups;

#ifndef binOpMin
#define binOpMin 1
#endif

#ifndef binOpMax
#define binOpMax 2
#endif

#ifndef binOpSum
#define binOpSum 3
#endif

#ifndef MIN
#define MIN(x, y) ((x < y) ? x : y)
#endif

#ifndef MAX
#define MAX(x, y) ((x > y) ? x : y)
#endif

#ifndef SUM
#define SUM(x,y) (x+y)
#endif

#ifndef GKYL_DEVICE_REDUCE_H
#define GKYL_DEVICE_REDUCE_H

extern "C" {
  void reductionBlocksAndThreads(GkDeviceProp *prop, int numElements, int maxBlocks,
                                 int maxThreads, int &blocks, int &threads);
  bool isPow2(unsigned int x);
  unsigned int nextPow2(unsigned int x);
  void reduceDeviceArray(int opIn, int numElements, int blocks, int threads, double *d_dataIn, double *d_dataOut);
}

template <unsigned int binOpType>
__inline__ __device__ double binOp(double a, double b) {
  double result = 0.0;
  if (binOpType==binOpMin) {
    result = MIN(a,b);
  } else if (binOpType==binOpMax) {
    result = MAX(a,b);
  } else if (binOpType==binOpSum) {
    result = SUM(a,b);
  }
  return result;
}

#endif  // GKYL_DEVICE_REDUCE_H
