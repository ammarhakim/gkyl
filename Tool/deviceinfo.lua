-- Gkyl ------------------------------------------------------------------------
--
-- Displays device (GPU) info
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

if not GKYL_HAVE_CUDA then
  print("Gkeyll not build with CUDA!")
end

local cuda = require "Cuda.RunTime"

-- fetch device and print information
local devNum, _ = cuda.GetDevice() 
local prop, err = cuda.GetDeviceProperties(devNum)

print("totalGlobalMem", prop.totalGlobalMem)
print("sharedMemPerBlock", prop.sharedMemPerBlock)
print("regsPerBlock", prop.regsPerBlock)
print("warpSize", prop.warpSize)
print("memPitch", prop.memPitch)
print("maxThreadsPerBlock", prop.maxThreadsPerBlock)
for d = 1, 2 do
  print(string.format(" maxThreadsDim[%d] = %d", d, prop.maxThreadsDim[d-1]))
end
for d = 1, 3 do
  print(string.format(" maxGridSize[%d] = %d", d, prop.maxThreadsDim[d-1]))
end
print("clockRate", prop.clockRate)
print("totalConstMem", prop.totalConstMem)
print(string.format("Major = %d, Minor = %d", prop.major, prop.minor))
print("deviceOverlap", prop.deviceOverlap)
print("multiProcessorCount", prop.multiProcessorCount)
print("kernelExecTimeoutEnabled", prop.kernelExecTimeoutEnabled)
print("integrated", prop.integrated)
print("canMapHostMemory", prop.canMapHostMemory)
print("computeMode", prop.computeMode)
print("concurrentKernels", prop.concurrentKernels)
print("ECCEnabled", prop.ECCEnabled)
print("asyncEngineCount", prop.asyncEngineCount)
print("unifiedAddressing", prop.unifiedAddressing)
print("memoryClockRate", prop.memoryClockRate)
print("memoryBusWidth", prop.memoryBusWidth)
print("l2CacheSize", prop.l2CacheSize)
print("maxThreadsPerMultiProcessor", prop.maxThreadsPerMultiProcessor)
print("streamPrioritiesSupported", prop.streamPrioritiesSupported)
print("globalL1CacheSupported", prop.globalL1CacheSupported)
print("localL1CacheSupported", prop.localL1CacheSupported)
print("sharedMemPerMultiprocessor", prop.sharedMemPerMultiprocessor)
print("regsPerMultiprocessor", prop.regsPerMultiprocessor)
print("managedMemory", prop.managedMemory)
print("isMultiGpuBoard", prop.isMultiGpuBoard)
print("multiGpuBoardGroupID", prop.multiGpuBoardGroupID)
print("singleToDoublePrecisionPerfRatio", prop.singleToDoublePrecisionPerfRatio)
print("pageableMemoryAccess", prop.pageableMemoryAccess)
print("concurrentManagedAccess", prop.concurrentManagedAccess)
print("computePreemptionSupported", prop.computePreemptionSupported)
print("canUseHostPointerForRegisteredMem", prop.canUseHostPointerForRegisteredMem)
print("cooperativeLaunch", prop.cooperativeLaunch)
print("cooperativeMultiDeviceLaunch", prop.cooperativeMultiDeviceLaunch)
print("pageableMemoryAccessUsesHostPageTables", prop.pageableMemoryAccessUsesHostPageTables)
print("directManagedMemAccessFromHost", prop.directManagedMemAccessFromHost)