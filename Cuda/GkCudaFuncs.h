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

typedef struct {
    char name[256];
    size_t totalGlobalMem;
    size_t sharedMemPerBlock;
    int regsPerBlock;
    int warpSize;
    size_t memPitch;
    int maxThreadsPerBlock;
    int maxThreadsDim[3];
    int maxGridSize[3];
    int clockRate;
    size_t totalConstMem;
    int major;
    int minor;
    size_t textureAlignment;
    size_t texturePitchAlignment;
    int deviceOverlap;
    int multiProcessorCount;
    int kernelExecTimeoutEnabled;
    int integrated;
    int canMapHostMemory;
    int computeMode;
    int concurrentKernels;
    int ECCEnabled;
    int asyncEngineCount;
    int unifiedAddressing;
    int memoryClockRate;
    int memoryBusWidth;
    int l2CacheSize;
    int maxThreadsPerMultiProcessor;
    int streamPrioritiesSupported;
    int globalL1CacheSupported;
    int localL1CacheSupported;
    size_t sharedMemPerMultiprocessor;
    int regsPerMultiprocessor;
    int managedMemory;
    int isMultiGpuBoard;
    int multiGpuBoardGroupID;
    int singleToDoublePrecisionPerfRatio;
    int pageableMemoryAccess;
    int concurrentManagedAccess;
    int computePreemptionSupported;
    int canUseHostPointerForRegisteredMem;
    int cooperativeLaunch;
    int cooperativeMultiDeviceLaunch;
    int pageableMemoryAccessUsesHostPageTables;
    int directManagedMemAccessFromHost;
} GkDeviceProp;

    int GkCuda_GetDeviceProp(GkDeviceProp *prop, int dev);
}

#endif // GK_CUDA_FUNCS_H
