// Gkyl ------------------------------------------------------------------------
//
// Functions for use in Cuda LuaJIT binding
//    _______     ___
// + 6 @ |||| # P ||| +
//------------------------------------------------------------------------------

#include <iostream>
#include <cstring>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_types.h>

#include <GkCudaFuncs.h>

// Error codes (enum cudaError)
GET_CUDA_OBJECT(int, cudaSuccess);
GET_CUDA_OBJECT(int, cudaErrorInvalidValue);
GET_CUDA_OBJECT(int, cudaErrorMemoryAllocation);
GET_CUDA_OBJECT(int, cudaErrorInvalidMemcpyDirection);

// Codes for CudaMemcpy (enum cudaMemcpyKind)
GET_CUDA_OBJECT(int, cudaMemcpyHostToHost);
GET_CUDA_OBJECT(int, cudaMemcpyHostToDevice);
GET_CUDA_OBJECT(int, cudaMemcpyDeviceToHost);
GET_CUDA_OBJECT(int, cudaMemcpyDeviceToDevice);
GET_CUDA_OBJECT(int, cudaMemcpyDefault);

GET_CUDA_OBJECT(unsigned, cudaMemAttachGlobal);
GET_CUDA_OBJECT(unsigned, cudaMemAttachHost);

int
GkCuda_GetDeviceProp(GkDeviceProp *prop, int dev) {
  cudaDeviceProp cuProp;
  int err = cudaGetDeviceProperties(&cuProp, dev);

  strncpy(prop->name, cuProp.name, 256);
  prop->totalGlobalMem = cuProp.totalGlobalMem;
  prop->sharedMemPerBlock = cuProp.sharedMemPerBlock ;
  prop->regsPerBlock = cuProp.regsPerBlock ;
  prop->warpSize = cuProp.warpSize ;
  prop->memPitch = cuProp.memPitch ;
  prop->maxThreadsPerBlock = cuProp.maxThreadsPerBlock ;
  for (auto i=0; i<3; ++i) {
    prop->maxThreadsDim[i] = cuProp.maxThreadsDim[i];
    prop->maxGridSize[i] = cuProp.maxGridSize[i];
  }
  prop->clockRate = cuProp.clockRate ;
  prop->totalConstMem = cuProp.totalConstMem ;
  prop->major = cuProp.major ;
  prop->minor = cuProp.minor ;
  prop->textureAlignment = cuProp.textureAlignment ;
  prop->texturePitchAlignment = cuProp.texturePitchAlignment ;
  prop->deviceOverlap = cuProp.deviceOverlap ;
  prop->multiProcessorCount = cuProp.multiProcessorCount ;
  prop->kernelExecTimeoutEnabled = cuProp.kernelExecTimeoutEnabled ;
  prop->integrated = cuProp.integrated ;
  prop->canMapHostMemory = cuProp.canMapHostMemory ;
  prop->computeMode = cuProp.computeMode ;
  prop->concurrentKernels = cuProp.concurrentKernels ;
  prop->ECCEnabled = cuProp.ECCEnabled ;
  prop->asyncEngineCount = cuProp.asyncEngineCount ;
  prop->unifiedAddressing = cuProp.unifiedAddressing ;
  prop->memoryClockRate = cuProp.memoryClockRate ;
  prop->memoryBusWidth = cuProp.memoryBusWidth ;
  prop->l2CacheSize = cuProp.l2CacheSize ;
  prop->maxThreadsPerMultiProcessor = cuProp.maxThreadsPerMultiProcessor ;
  prop->streamPrioritiesSupported = cuProp.streamPrioritiesSupported ;
  prop->globalL1CacheSupported = cuProp.globalL1CacheSupported ;
  prop->localL1CacheSupported = cuProp.localL1CacheSupported ;
  prop->sharedMemPerMultiprocessor = cuProp.sharedMemPerMultiprocessor ;
  prop->regsPerMultiprocessor = cuProp.regsPerMultiprocessor ;
  prop->managedMemory = cuProp.managedMemory ;
  prop->isMultiGpuBoard = cuProp.isMultiGpuBoard ;
  prop->multiGpuBoardGroupID = cuProp.multiGpuBoardGroupID ;
  prop->singleToDoublePrecisionPerfRatio = cuProp.singleToDoublePrecisionPerfRatio ;
  prop->pageableMemoryAccess = cuProp.pageableMemoryAccess ;
  prop->concurrentManagedAccess = cuProp.concurrentManagedAccess ;
  prop->computePreemptionSupported = cuProp.computePreemptionSupported ;
  prop->canUseHostPointerForRegisteredMem = cuProp.canUseHostPointerForRegisteredMem ;
  prop->cooperativeLaunch = cuProp.cooperativeLaunch ;
  prop->cooperativeMultiDeviceLaunch = cuProp.cooperativeMultiDeviceLaunch ;
  prop->pageableMemoryAccessUsesHostPageTables = cuProp.pageableMemoryAccessUsesHostPageTables ;
  prop->directManagedMemAccessFromHost = cuProp.directManagedMemAccessFromHost ;

  return err;
}
