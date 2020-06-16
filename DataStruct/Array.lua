-- Gkyl ------------------------------------------------------------------------
--
-- General multidimensional array object. All Gkeyll datastructures
-- are based on array.
---
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

local Alloc = require "Lib.Alloc"

-- Array object is a low-level data containter for storing data. The
-- interface to this object is delibrately kept low-level and via the
-- C-struct definition of the array directly. This is done so as to
-- make communication with C-code (including CUDA) easy and via a
-- simple interface.
--
-- Arrays are reference counted. However, for this to work one must be
-- disciplined. If you need to store a pointer to the C-struct in
-- another object or C-struct, please obtain it with the aquire()
-- method. When you are done with the pointer release it using the
-- release() method. Release can be automated by proper use of __gc
-- method in the struct metamethods. If you do not do this memory can
-- leak and, very seriously, can disappear under you.
-- 
-- In general, no ordering is imposed on the data. Ordering needs to
-- be provided elsewhere. For example, there is no mechanism to loop
-- over the array in row/col major order. If one needs to do this one
-- can use the Range object to impose an order and do looping.

-- Array structure
ffi.cdef [[
  typedef struct {
    // rank, shape, total elements, reference count, element size
    long r, *s, N, c, sz;
    void* d; // pointer to data
  } GkylArray_t;
]]

-- (Meta) function to create Array constructor for given type
local function Array_meta_ctor(elct)
   local sz = sizeof(elct)
   local arrayCt = typeof("GkylArray_t")
   local arrayPtrCt = typeof("struct { $*p; }", elct)

   local array_fn = {
      data = function (self)
	 local sptr = arrayPtrCt(self.d)
	 return sptr.p
      end,
      clone = function (self)
	 local a = new(self)
	 a.s = Alloc.malloc(a.r*sizeof("long"))
	 a.d = Alloc.calloc(a.N, sz)
	 a.r, a.N, a.sz = self.r, self.N, self.sz
	 a.c = 1 -- newly created
	 ffi.copy(a.s, self.s, a.r*sizeof("long"))
	 ffi.copy(a.d, self.d, a.N*sz)
	 return a
      end,
      aquire = function (self)
	 self.c = self.c + 1
	 return self
      end,
      release = function (self)
	 self.c = self.c - 1
	 if self.c == 0 then
	    Alloc.free(self.s); Alloc.free(self.d)
	 end
      end,
   }
   
   local array_mt = {
      __new = function (self, shape)
	 local a = new(self)
	 a.r = #shape
	 a.s = Alloc.malloc(a.r*sizeof("long"))
	 a.N = 1
	 for i = 1, #shape do
	    a.s[i-1] = shape[i]
	    a.N = a.N*shape[i]
	 end
	 a.c = 1
	 a.sz = sz
	 a.d = Alloc.calloc(a.N, sz)
	 return a
      end,
      __gc = function(self)
	 self:release()
      end,
      __index = array_fn
   }
   return metatype(arrayCt, array_mt)
end

return {
   Array_meta_ctor,
   Array = Array_meta_ctor(typeof("double")),
}
