-- Gkyl ------------------------------------------------------------------------
-- CUDA runtime functions used in Gkeyll
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- don't bother if we don't have CUDA
assert(GKYL_HAVE_CUDA, "Gkyl was not built with CUDA!")

local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, typeof = xsys.from(ffi,
     "new, typeof")

local _M = {}

-- CUDA types and functions
ffi.cdef [[
  // typedefs for use in API below
  typedef int cudaMemcpyKind;

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

  // Run-time information
  int cudaDriverGetVersion(int *driverVersion);
  int cudaGetDevice(int* device);
  int cudaGetDeviceCount ( int* count );
  int GkCuda_GetDeviceProp(GkDeviceProp *prop, int dev);

  // Device management
  int cudaDeviceSynchronize();

  // Memory management
  int cudaMalloc(void **devPtr, size_t size);
  int cudaFree(void *ptr);
  int cudaMemcpy(void* dst, const void* src, size_t count, cudaMemcpyKind kind);
  int cudaMallocManaged(void** devPtr, size_t size, unsigned int flags);

  // error codes
  int get_cudaSuccess();
  int get_cudaErrorInvalidValue();
  int get_cudaErrorMemoryAllocation();
  int get_cudaErrorInvalidMemcpyDirection();

  // MemcpyKind flags
  int get_cudaMemcpyHostToHost();
  int get_cudaMemcpyHostToDevice();
  int get_cudaMemcpyDeviceToHost();
  int get_cudaMemcpyDeviceToDevice();
  int get_cudaMemcpyDefault();

  // Flags for cudaMallocManaged
  unsigned get_cudaMemAttachGlobal();
  unsigned get_cudaMemAttachHost();
]]

-- CUDA runtime error codes
_M.Success = ffi.C.get_cudaSuccess()
_M.ErrorInvalidValue = ffi.C.get_cudaErrorInvalidValue()
_M.ErrorMemoryAllocation = ffi.C.get_cudaErrorMemoryAllocation()
_M.ErrorInvalidMemcpyDirection = ffi.C.get_cudaErrorInvalidMemcpyDirection()

-- CUDA MemcpyKind flags
_M.MemcpyHostToHost = ffi.C.get_cudaMemcpyHostToHost()
_M.MemcpyHostToDevice = ffi.C.get_cudaMemcpyHostToDevice()
_M.MemcpyDeviceToHost = ffi.C.get_cudaMemcpyDeviceToHost()
_M.MemcpyDeviceToDevice = ffi.C.get_cudaMemcpyDeviceToDevice()
_M.MemcpyDefault = ffi.C.get_cudaMemcpyDefault()

-- Flags for cudaMallocManaged
_M.MemAttachGlobal = ffi.C.get_cudaMemAttachGlobal()
_M.MemAttachHost = ffi.C.get_cudaMemAttachHost()

-- some types for use in CUDA functions
local int_1 = typeof("int[1]")
local uint_1 = typeof("unsigned[1]")
local voidp = typeof("void *[1]")

-- Note: Most functions below return error code as last return value

-- cudaDriverGetVersion
function _M.DriverGetVersion()
    local r = int_1()
    local err = ffiC.cudaDriverGetVersion(r)
    return r[0], err
end
-- cudaGetDevice
function _M.GetDevice()
   local r = int_1()
   local err = ffiC.cudaGetDevice(r)
   return r[0], err
end
-- cudaGetDeviceCount
function _M.GetDeviceCount()
   local r = int_1()
   local err = ffiC.cudaGetDeviceCount(r)
   return r[0], err
end
-- cudaGetDeviceProperties
function _M.GetDeviceProperties(device)
   local prop = ffi.new("GkDeviceProp[1]")
   local err = ffiC.GkCuda_GetDeviceProp(prop, device)
   return prop[0], err
end

-- cudaDeviceSynchronize
function _M.DeviceSynchronize()
   local err = ffiC.cudaDeviceSynchronize()
   return err
end

-- cudaMalloc
function _M.Malloc(sz)
   local devPtr = voidp()
   local err = ffiC.cudaMalloc(devPtr, sz)
   return devPtr[0], err
end
-- cudaFree
function _M.Free(d)
   if d then
      ffiC.cudaFree(d); d = nil
   end
end
-- cudaMallocManaged
function _M.MallocManaged(sz, kind)
   local devPtr = voidp()
   local myKind = kind
   -- CUDA specs say that default is cudaMemAttacedGlobal
   if kind == nil then myKind = _M.MemAttachGlobal end
   local err = ffiC.cudaMallocManaged(devPtr, sz, myKind)
   return devPtr[0], err
end

-- cudaMemcpy
function _M.Memcpy(dst, src, count, kind)
    return ffiC.cudaMemcpy(dst, src, count, kind)
end

return _M