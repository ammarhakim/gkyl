-- Gkyl ------------------------------------------------------------------------
-- CUDA runtime functions used in Gkeyll
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- don't bother if we don't have CUDA
assert(GKYL_HAVE_CUDA, "Gkyl was not built with CUDA!")

local ffi  = require "ffi"
local ffiC = ffi.C

local _M = {}

ffi.cdef [[
  // Run-time information
  int cudaDriverGetVersion(int *driverVersion);

  // Memory management
  int cudaMalloc(void **devPtr, size_t size);
  int cudaFree(void *ptr);
]]

-- some types for use in CUDA functions
local int_1 = typeof("int[1]")
local uint_1 = typeof("unsigned[1]")
local voidp = typeof("void *[1]")

-- cudaMalloc
function _M.Malloc(sz)
   local devPtr = voidp()
   local err = ffiC.cudaMalloc(devPtr, sz)
   
   assert(d, "cudaMalloc: Unable to allocate memory!")
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

return _M