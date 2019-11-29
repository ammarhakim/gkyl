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
    DECL_GET_CUDA_OBJECT(int, cudaSuccess);
    DECL_GET_CUDA_OBJECT(int, cudaErrorInvalidValue);
    DECL_GET_CUDA_OBJECT(int, cudaErrorMemoryAllocation);
}

#endif // GK_CUDA_FUNCS_H
