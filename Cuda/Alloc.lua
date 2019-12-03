-- Gkyl ------------------------------------------------------------------------
--
-- Allocators for memory on GPU device. The data returned by this
-- can't be diretly indexed as it is device memory
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local cuda = require "Cuda.RunTime"
local ffi  = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Alloc ----------------------------------------------------------------------
--
-- A GPU device memory object
--------------------------------------------------------------------------------

local function Alloc_meta_ctor(elct)
   local elmSz = sizeof(elct) -- element size in bytes
   -- block size in bytes
   local blockSz = 2048*elmSz
   -- function to adjust allocation bytes to block-size: elmSz is the
   -- size of each element, num is number of such elements
   local function calcAdjustedSize(num)
      local adjBytes = (math.floor(elmSz*num/blockSz)+1)*blockSz
      local adjNum = math.floor(adjBytes/elmSz)
      return adjBytes, adjNum
   end

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
   elseif ffi.istype(new(elct), new("complex")) then
      isNumberType = true
   else
      isNumberType = false
   end   

   -- Use cuda memory allocator to allocate the memory on device
   --
   local function alloc(ct, num)
      local adjBytes, adjNum = calcAdjustedSize(num)
      
      local v = new(ct)
      v._capacity = 0
      local err
      v._data, err = cuda.Malloc(adjNum*elmSz)
      -- not exactly sure what to do when cudaMalloc fails: aborting
      -- perhaps is the best solution but may lead to issues with the
      -- device
      assert(cuda.Success == err, "Cuda.Alloc: unable to allocate device memory!")
      v._capacity = adjNum
      v._size = num
      return v
   end
   
   local alloc_funcs = {
      elemType = function (self)
	 return elct
      end,
      elemSize = function (self)
	 return sizeof(elct)
      end,
      data = function (self)
	 return self._data
      end,
      capacity = function (self)
	 return self._capacity
      end,
      size = function(self)
	 return self._size
      end,
      copyToHost = function(self, hd)
	 return cuda.Memcpy(hd:data(), self:data(), self:elemSize()*self:size(), cuda.MemcpyDeviceToHost)
      end,
      copyFromHost = function(self, hd)
	 return cuda.Memcpy(self:data(), hd:data(), self:elemSize()*self:size(), cuda.MemcpyHostToDevice)
      end,
      delete = function (self)
	 if self._capacity == 0 then
	    return false
	 else
	    cuda.Free(self._data)
	    self._capacity = 0
	 end
	 return true
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
	 return alloc_funcs[k]
      end,
      __gc = function (self)
	 self:delete()
      end,
   }
   return metatype(typeof("struct { int32_t _size; int32_t _capacity; $* _data; }", elct), alloc_mt)
end

-- function to create an allocator for custom type
local function createAllocator(typeStr)
   return Alloc_meta_ctor(ffi.typeof(typeStr))
end

return {
   Alloc_meta_ctor = Alloc_meta_ctor,
   Double = createAllocator("double"),
   Float = createAllocator("float"),
   Int = createAllocator("int"),
   createAllocator = createAllocator,
}
