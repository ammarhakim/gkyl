// Gkyl ------------------------------------------------------------------------
//
// Functions for use in Cuda LuaJIT binding
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <iostream>

#include <cuda_runtime.h>
#include <GkCudaFuncs.h>

// Error codes (enum cudaError)
GET_CUDA_OBJECT(int, cudaSuccess);
GET_CUDA_OBJECT(int, cudaErrorInvalidValue);
GET_CUDA_OBJECT(int, cudaErrorMemoryAllocation);
