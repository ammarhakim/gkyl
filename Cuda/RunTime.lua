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

ffi.cdef [[
  // typedefs for use in API below
  typedef int cudaMemcpyKind;

  // Run-time information
  int cudaDriverGetVersion(int *driverVersion);

  // Memory management
  int cudaMalloc(void **devPtr, size_t size);
  int cudaFree(void *ptr);
  int cudaMemcpy(void* dst, const void* src, size_t count, cudaMemcpyKind kind);

  // error codes
  int get_cudaSuccess();
  int get_cudaErrorInvalidValue();
  int get_cudaErrorMemoryAllocation();

  // MemcpyKind flags
  int get_cudaMemcpyHostToHost();
  int get_cudaMemcpyHostToDevice();
  int get_cudaMemcpyDeviceToHost();
  int get_cudaMemcpyDeviceToDevice();
  int get_cudaMemcpyDefault();
]]


-- CUDA runtime error codes
_M.Success = ffi.C.get_cudaSuccess()
_M.ErrorInvalidValue = ffi.C.get_cudaErrorInvalidValue()
_M.ErrorMemoryAllocation = ffi.C.get_cudaErrorMemoryAllocation()

-- CUDA MemcpyKind flags
_M.MemcpyHostToHost = ffi.C.get_cudaMemcpyHostToHost()
_M.MemcpyHostToDevice = ffi.C.get_cudaMemcpyHostToDevice()
_M.MemcpyDeviceToHost = ffi.C.get_cudaMemcpyDeviceToHost()
_M.MemcpyDeviceToDevice = ffi.C.get_cudaMemcpyDeviceToDevice()
_M.MemcpyDefault = ffi.C.get_cudaMemcpyDefault()

-- some types for use in CUDA functions
local int_1 = typeof("int[1]")
local uint_1 = typeof("unsigned[1]")
local voidp = typeof("void *[1]")

-- cudaMalloc
function _M.Malloc(sz)
   local devPtr = voidp()
   local err = ffiC.cudaMalloc(devPtr, sz)
   assert(err == _M.Success, "Unable to allocate device memory")
   return devPtr[0]
end
-- cudaFree
function _M.Free(d)
   if d then
      ffiC.cudaFree(d); d = nil
   end
end

-- cudaDriverGetVersion
function _M.DriverGetVersion()
    local r = int_1()
    local _ = ffiC.cudaDriverGetVersion(r)
    return r[0]
end

-- cudaMemcpy
function _M.Memcpy(dst, src, count, kind)
    return ffiC.cudaMemcpy(dst, src, count, kind)
end

return _M