-- Gkyl ------------------------------------------------------------------------
--
-- Allocator for memory that is shared between processes. The
-- processes must live on a communicator created using MPI-SHM calls.
-- 
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

local Mpi = require "Mpi"

-- AllocShared -----------------------------------------------------------------
--
-- A 1D shared array type
--------------------------------------------------------------------------------

local function AllocShared_meta_ctor(elct)
   local elmSz = sizeof(elct) -- element size in bytes
   -- copy function for non-numeric types: this is used in methods
   -- that set array values if the element type stored is not numeric
   local isNumberType = false
   local copyElemFunc = nil
   if ffi.istype(new(elct), new("double")) then
      isNumberType = true
   elseif ffi.istype(new(elct), new("float")) then
      isNumberType = true
   elseif ffi.istype(new(elct), new("int")) then
      isNumberType = true
   else
      isNumberType = false
      copyElemFunc = function (dst, src) ffi.copy(dst, src, elmSz) end
   end   

   local alloc_funcs = {
      elemType = function (self)
	 return elct
      end,
      elemSize = function (self)
	 return sizeof(elct)
      end,
   }
   local alloc_mt = {
      __new = function (ct, num)
	 if num then
	    return alloc(ct, num)
	 else
	    return alloc(ct, 0)
	 end
      end,
      __index = function (self, k)
	 -- for non-numeric types, a reference will be returned: this
	 -- is a little dangerous, but it allows for more "natural"
	 -- style of programming when working with C-structs
	 if type(k) == "number" then
	    return self._data[k-1]
	 else
	    return alloc_funcs[k]
	 end
      end,
      __newindex = copyElemFunc and
	 function (self, k, v)
	    copyElemFunc(self._data[k-1], v)
	 end or
	 function (self, k, v)
	    self._data[k-1] = v
	 end,
      __gc = function (self)
	 self:delete()
      end,
   }
   return metatype(typeof("struct { int32_t _size; $* _data; }", elct), alloc_mt)
end

-- function to create an allocator for custom type
local function createAllocator(typeStr)
   return AllocShared_meta_ctor(ffi.typeof(typeStr))
end

return {
   AllocShared_meta_ctor = AllocShared_meta_ctor,
   Double = createAllocator("double"),
   Float = createAllocator("float"),
   Int = createAllocator("int"),
}
