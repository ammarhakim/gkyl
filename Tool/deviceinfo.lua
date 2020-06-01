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

local fmt = "%-38s: %s"

print(string.format(fmt, "totalGlobalMem", prop.totalGlobalMem))
print(string.format(fmt, "sharedMemPerBlock", prop.sharedMemPerBlock))
print(string.format(fmt, "regsPerBlock", prop.regsPerBlock))
print(string.format(fmt, "warpSize", prop.warpSize))
print(string.format(fmt, "memPitch", prop.memPitch))
print(string.format(fmt, "maxThreadsPerBlock", prop.maxThreadsPerBlock))
for d = 1, 3 do
  print(string.format(fmt, string.format(" maxThreadsDim[%d]", d), prop.maxThreadsDim[d-1]))
end
for d = 1, 3 do
  print(string.format(fmt, string.format(" maxGridSize[%d]", d), prop.maxGridSize[d-1]))
end
print(string.format(fmt, "clockRate", prop.clockRate))
print(string.format(fmt, "totalConstMem", prop.totalConstMem))
print(string.format(fmt, "Major:Minor", string.format("%s:%s", prop.major, prop.minor)))
print(string.format(fmt, "deviceOverlap", prop.deviceOverlap))
print(string.format(fmt, "multiProcessorCount", prop.multiProcessorCount))
print(string.format(fmt, "kernelExecTimeoutEnabled", prop.kernelExecTimeoutEnabled))
print(string.format(fmt, "integrated", prop.integrated))
print(string.format(fmt, "canMapHostMemory", prop.canMapHostMemory))
print(string.format(fmt, "computeMode", prop.computeMode))
print(string.format(fmt, "concurrentKernels", prop.concurrentKernels))
print(string.format(fmt, "ECCEnabled", prop.ECCEnabled))
print(string.format(fmt, "asyncEngineCount", prop.asyncEngineCount))
print(string.format(fmt, "unifiedAddressing", prop.unifiedAddressing))
print(string.format(fmt, "memoryClockRate", prop.memoryClockRate))
print(string.format(fmt, "memoryBusWidth", prop.memoryBusWidth))
print(string.format(fmt, "l2CacheSize", prop.l2CacheSize))
print(string.format(fmt, "maxThreadsPerMultiProcessor", prop.maxThreadsPerMultiProcessor))
print(string.format(fmt, "streamPrioritiesSupported", prop.streamPrioritiesSupported))
print(string.format(fmt, "globalL1CacheSupported", prop.globalL1CacheSupported))
print(string.format(fmt, "localL1CacheSupported", prop.localL1CacheSupported))
print(string.format(fmt, "sharedMemPerMultiprocessor", prop.sharedMemPerMultiprocessor))
print(string.format(fmt, "regsPerMultiprocessor", prop.regsPerMultiprocessor))
print(string.format(fmt, "managedMemory", prop.managedMemory))
print(string.format(fmt, "isMultiGpuBoard", prop.isMultiGpuBoard))
print(string.format(fmt, "multiGpuBoardGroupID", prop.multiGpuBoardGroupID))
print(string.format(fmt, "singleToDoublePrecisionPerfRatio", prop.singleToDoublePrecisionPerfRatio))
print(string.format(fmt, "pageableMemoryAccess", prop.pageableMemoryAccess))
print(string.format(fmt, "concurrentManagedAccess", prop.concurrentManagedAccess))
print(string.format(fmt, "computePreemptionSupported", prop.computePreemptionSupported))
print(string.format(fmt, "canUseHostPointerForRegisteredMem", prop.canUseHostPointerForRegisteredMem))
print(string.format(fmt, "cooperativeLaunch", prop.cooperativeLaunch))
print(string.format(fmt, "cooperativeMultiDeviceLaunch", prop.cooperativeMultiDeviceLaunch))
print(string.format(fmt, "pageableMemoryAccessUsesHostPageTables", prop.pageableMemoryAccessUsesHostPageTables))
print(string.format(fmt, "directManagedMemAccessFromHost", prop.directManagedMemAccessFromHost))