-- Gkyl ------------------------------------------------------------------------
--
-- Utility functions to make dealing with C-structs easier
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
if GKYL_HAVE_CUDA then
   cuda = require "Cuda.RunTime"
end

-- This function returns a table that implements some methods used in
-- struct wrappers. This is mainly to allow adding cloneOnDevice
-- function to C-structs to which a metatype may be already
-- attached. The table returned by this function should be appended to
-- the metatype of a struct to provide the cloneOnDevice method
local function getCloneOnDeviceFunctionTable(elct)
   return {
      cloneOnDevice = GKYL_HAVE_CUDA and
	 function (self)
	    local cuObj, err = cuda.Malloc(ffi.sizeof(elct))
	    assert(cuda.Success == err, "StructUtils.Struct: unable to allocate device memory!")
	    cuda.Memcpy(cuObj, self, ffi.sizeof(elct), cuda.MemcpyHostToDevice)
	    return cuObj
	 end or
	 function (self)
	    -- when not on device, return clone on host
	    local v = ffi.new(elct)
	    ffi.copy(v, self, ffi.sizeof(elct))
	    return v
	 end
   }
end

-- Struct ----------------------------------------------------------------------
--
-- A wrapper around a C-struct to make some functionality cleaner
--------------------------------------------------------------------------------

local function Struct(elct)
   local struct_mf = {
      size = function(self)
	 return ffi.sizeof(elct)
      end,      
      clone = function(self)
	 -- clone does a shallow copy: i.e. pointers in original and
	 -- cloned struct point to the same address
	 local v = ffi.new(self)
	 ffi.copy(v, self, ffi.sizeof(elct))
	 return v
      end,
      copy = function(self, src)
	 -- copy does a shallow copy: i.e. pointers in original and
	 -- copied struct point to the same address
	 ffi.copy(self, src, ffi.sizeof(elct))
      end,
      cloneOnDevice = GKYL_HAVE_CUDA and
	 function (self)
	    local cuObj, err = cuda.Malloc(ffi.sizeof(elct))
	    assert(cuda.Success == err, "StructUtils.Struct: unable to allocate device memory!")
	    cuda.Memcpy(cuObj, self, ffi.sizeof(elct), cuda.MemcpyHostToDevice)
	    return cuObj
	 end or
	 function (self)
	    -- when not on device, return host clone
	    return self.clone()
	 end
   }
   local struct_mt = {
      __new = function(ct, init)
	 -- 'init' is a table of values values to initialize
	 -- struct. Values must be in same order as entries in
	 -- C-struct
	 if init then
	    return ffi.new(ct, init)
	 end
	 return ffi.new(ct)
      end,
      __index = function (self, k)
	 return struct_mf[k]
      end,
   }
   return ffi.metatype(elct, struct_mt)
end

return {
   Struct = Struct,
   getCloneOnDeviceFunctionTable = getCloneOnDeviceFunctionTable,
}
