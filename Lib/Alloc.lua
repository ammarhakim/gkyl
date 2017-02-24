-- Gkyl ------------------------------------------------------------------------
--
-- Allocators for memory not managed by Lua garbage-collector: allows
-- allocation of blocks not bound by Lua GC limitations.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Declare malloc/free functions from C library
ffi.cdef [[
  void* malloc(size_t size);
  void* calloc(size_t nitems, size_t size);
  void *realloc(void *ptr, size_t new_size);
  void free(void *ptr);
]]

-- Wrapper around core memory functions
local function malloc(sz)
   local d = ffi.C.malloc(sz)
   assert(d, "Alloc.malloc: Unable to allocate memory!")
   return d
end
local function calloc(nitems, sz)
   local d = ffi.C.calloc(nitems, sz)
   assert(d, "Alloc.calloc: Unable to allocate memory!")
   return d
end
local function realloc(dold, sz)
   local d = ffi.C.realloc(dold, sz)
   assert(d, "Alloc.realloc: Unable to allocate memory!")
   return d
end
local function free(d)
   if d then
      ffi.C.free(d); d = nil
   end
end

 -- block size in bytes
local blockSz = 16*1024
-- function to adjust allocation bytes to block-size: elmSz is the
-- size of each element, num is number of such elements
local function calcAdjustedSize(elmSz, num)
   local adjBytes = (math.floor(elmSz*num/blockSz)+1)*blockSz
   local adjNum = math.floor(adjBytes/elmSz)
   return adjBytes, adjNum
end

-- Use calloc to allocate 'num' elements of memory for type
-- 'elct'. The type 'ct' is a struct holding size and pointer to
-- allocated memory. This allocator rounds out the memory to "blockSz"
-- bytes to prevent fragmentation of the C memory allocator.
--
-- NOTE: We use "calloc" and not malloc as calloc zeros out the
-- memory. Uninitialized memory can cause random crashes in LuaJIT.
local function alloc(ct, elct, num)
   local adjBytes, adjNum = calcAdjustedSize(sizeof(elct), num)

   local v = new(ct)
   v._capacity = 0
   v._data = calloc(adjNum, sizeof(elct))
   v._capacity = adjNum
   v._size = num
   return v
end

local function Alloc_meta_ctor(elct)
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
      expand = function (self, nSize)
	 -- don't do anything if requested size is less than capacity
	 if nSize < self._capacity then
	    self._size = nSize
	    return 
	 end

	 -- reallocate memory
	 local adjBytes, adjNum = calcAdjustedSize(self:elemSize(), nSize)
	 self._data = realloc(self._data, adjBytes)
	 self._capacity = adjNum
	 self._size = nSize
      end,
      clear = function(self)
	 local adjBytes, adjNum = calcAdjustedSize(self:elemSize(), 0) -- single block
	 self._data = realloc(self._data, adjBytes)
	 self._capacity = adjNum
	 self._size = 0
      end,
      fill = function (self, v)
	 for i = 0, self._size-1 do
	    self._data[i] = v
	 end
      end,
      push = function (self, v)
	 self:expand(self._size+1)
	 self._data[self._size-1] = v
      end,
      last = function (self)
	 return self._data[self._size-1]
      end,
      delete = function (self)
	 if self._capacity == 0 then
	    return false
	 else
	    free(self._data)
	    self._capacity = 0
	 end
	 return true
      end,
   }
   local alloc_mt = {
      __new = function (ct, num)
	 if num then
	    return alloc(ct, elct, num)
	 else
	    return alloc(ct, elct, 0)
	 end
      end,
      __index = function (self, k)
	 if type(k) == "number" then
	    return self._data[k-1]
	 else
	    return alloc_funcs[k]
	 end
      end,
      __newindex = function (self, k, v)
	 self._data[k-1] = v
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
   malloc = malloc,
   realloc = realloc,
   calloc = calloc,
   free = free,
   createAllocator = createAllocator,
}
