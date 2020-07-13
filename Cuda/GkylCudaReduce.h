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

extern "C" {
  void reductionBlocksAndThreads(GkDeviceProp *prop, int numElements, int maxBlocks,
                                 int maxThreads, int &blocks, int &threads);
  bool isPow2(unsigned int x);
  unsigned int nextPow2(unsigned int x);

  // Type of function that operates on two numbers.
  typedef double (*redBinOpFunc_t)(double a, double b);
  // Base reduction operator object.
  struct baseReduceOp {
    double initValue;               // Initial value to start with.
    redBinOpFunc_t reduceFunc { };  // Provided by "children".
  
    __host__ __device__ double reduce(double a, double b) {
      return reduceFunc(a, b);
    }
  };
  
  redBinOpFunc_t getRedMinFuncFromDevice();
  redBinOpFunc_t getRedMaxFuncFromDevice();
  redBinOpFunc_t getRedSumFuncFromDevice();
  
  void reduceDeviceArray(baseReduceOp *opIn, int numElements, int blocks,
                         int threads, double *d_dataIn, double *d_dataOut);
}

namespace Gkyl {
  __inline__ __device__ double MIN(double a, double b) {
    return ((a < b) ? a : b);
  }
  __inline__ __device__ double MAX(double a, double b) {
    return ((a > b) ? a : b);
  }
  __inline__ __device__ double SUM(double a, double b) {
    return a+b;
  }
}
