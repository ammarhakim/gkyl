// Gkyl ------------------------------------------------------------------------
//
// Functions for use in Cuda LuaJIT binding
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#ifndef GK_CUDA_FUNCS_H
#define GK_CUDA_FUNCS_H

#include <GkCudaMacros.h>

extern "C" {
    // Pre-defined objects and constants

    // Error codes (enum cudaError)
    DECL_GET_CUDA_OBJECT(int, cudaSuccess);
    DECL_GET_CUDA_OBJECT(int, cudaErrorInvalidValue);
    DECL_GET_CUDA_OBJECT(int, cudaErrorMemoryAllocation);
    DECL_GET_CUDA_OBJECT(int, cudaErrorInvalidMemcpyDirection);

    // Codes for CudaMemcpy (enum cudaMemcpyKind)
    DECL_GET_CUDA_OBJECT(int, cudaMemcpyHostToHost);
    DECL_GET_CUDA_OBJECT(int, cudaMemcpyHostToDevice);
    DECL_GET_CUDA_OBJECT(int, cudaMemcpyDeviceToHost);
    DECL_GET_CUDA_OBJECT(int, cudaMemcpyDeviceToDevice);
    DECL_GET_CUDA_OBJECT(int, cudaMemcpyDefault);
}

#endif // GK_CUDA_FUNCS_H
