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
  void free(void *ptr);
  void *memset(void *s, int c, size_t n);
]]

-- Wrapper around core memory functions
local function malloc(sz)
   local d = ffi.C.malloc(sz)
   assert(d, "malloc: Unable to allocate memory!")
   return d
end
local function free(d)
   if d then
      ffi.C.free(d); d = nil
   end
end

 -- block size in bytes
local blockSz = 16*1024
-- Use malloc to allocate 'num' elements of memory for type
-- 'elct'. The type 'ct' is a struct holding size and pointer to
-- allocated memory. This allocator rounds out the memory to "blockSz"
-- bytes to prevent fragmentation of the C memory allocator
local function alloc(ct, elct, num)
   local sz = sizeof(elct)
   local adjBytes = (math.floor(sz*num/blockSz)+1)*blockSz
   local adjNum = math.floor(adjBytes/sz)

   local v = new(ct)
   v._capacity = 0
   v._data = ffi.C.malloc(adjBytes)
   assert(v._data, "Alloc.alloc: Unable to allocate memory!")
   v._capacity = adjNum
   return v
end

local function Alloc_meta_ctor(elct)
   local alloc_funcs = {
      elemType = function (self)
	 return elct
      end,
      data = function (self)
	 return self._data
      end,
      capacity = function (self)
	 return self._capacity
      end,
      fill = function (self, v)
	 for i = 0, self._capacity-1 do
	    self._data[i] = v
	 end
      end,
      delete = function (self)
	 if self._capacity == 0 then
	    return false
	 else
	    ffi.C.free(self._data); self._data = nil
	    self._capacity = 0
	 end
	 return true
      end,
   }
   local alloc_mt = {
      __new = function (ct, num)
	 return alloc(ct, elct, num)
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
   return metatype(typeof("struct { uint32_t _capacity; $* _data; }", elct), alloc_mt)
end

return {
   Alloc_meta_ctor = Alloc_meta_ctor,
   Double = Alloc_meta_ctor(typeof("double")),
   Float = Alloc_meta_ctor(typeof("float")),
   Int = Alloc_meta_ctor(typeof("int")),
   malloc = malloc,
   free = free,
   memset = memset,
}
