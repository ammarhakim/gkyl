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

_M = {}

_M.__help = [=[

Array object is a data containter for storing multi-dimensional
data. The interface to this object is delibrately kept low-level and
is accessible via the C-struct definition of the array directly. This
is done so as to make communication with C-code (including CUDA) easy
and via a simple interface.

Arrays are reference counted. However, for this to work one must be
disciplined. If you need to store a pointer to the array in another
object or C-struct, please obtain it with the aquire() method. When
you are done with the pointer release it using the release()
method. Release can be automated by proper use of __gc method in the
struct metamethods. If you do not do this memory can leak and, very
seriously, can disappear under you.

In general, no ordering is imposed on the data. Ordering needs to be
provided elsewhere. For example, there is no mechanism to loop over
the array in row/col major order. If one needs to do this one can use
the Range object to impose an order and do looping.

Array structure. The "d" pointer points data stored in the array. This
is stored as void* and so needs to be cast to the proper C-pointer
before use. To do this, if 'arr' is the array and it stores doubles
then do:

local d = ffi.cast("double*", arr.d)

At this point, d[0] ... d[n-1] give you access to the values stored in
the array.

Recall that the size and shape are stored as long int. This is needed
as int are too small to store information for very large array we may
need in many simulations. Hence, you need to use tonumber() Lua method
to convert the size and shape to something Lua can use (for example in
loops).

 ]=]

ffi.cdef [[
  typedef struct {
    int r, c, sz; // rank, use-count, element-size
    long n, *s; // total number of elements and shape
    void *d; // pointer to data
  } GkylArray_t;
]]

local longSz = sizeof("long")

-- Various types for arrays of basic C-types
_M.float = typeof('float')
_M.double = typeof('double')
_M.int = typeof('int')
_M.long = typeof('long')
_M.char = typeof('char')

-- Array ctype
local ArrayCt = typeof("GkylArray_t")

local array_fn = {
   rank = function (self)
      return self.r
   end,
   shape = function (self)
      local s = {}
      for i = 1, self.r do
	 s[i] = tonumber(self.s[i-1])
      end
      return s
   end,
   size = function (self)
      return tonumber(self.n)
   end,
   elemSize = function (self)
      return self.sz
   end,
   clone = function (self)
      local a = new(self)
      a.r, a.n, a.sz = self.r, self.n, self.sz
      
      a.d = Alloc.calloc(self.n, self.sz)
      ffi.copy(a.d, self.d, self.n*self.sz)

      a.s = Alloc.malloc(self.r*longSz)
      ffi.copy(a.s, self.s, self.r*longSz)

      a.c = 1 -- newly created

      return a
   end,
   aquire = function (self)
      self.c = self.c+1
      return self
   end,
   release = function (self)
      self.c = self.c-1
      if self.c == 0 then
	 Alloc.free(self.s)
	 Alloc.free(self.d)
      end
   end,
}

local array_mt = {
   __new = function (self, shape, atype)
      local a = new(self)
      a.r = #shape
      -- compute total number of elements
      a.n = 1
      for i = 1, #shape do a.n = a.n*shape[i] end
      a.c, a.sz = 1, sizeof(atype)

      -- allocate memory
      a.d = Alloc.calloc(a.n, a.sz)

      -- set shape of array
      a.s = Alloc.malloc(a.r*longSz)
      for i = 1, #shape do a.s[i-1] = shape[i] end
      return a
   end,
   __gc = function(self)
      self:release()
   end,
   __index = array_fn
}
local ArrayCtor = metatype(ArrayCt, array_mt)

-- Construct array of given shape and type
_M.Array = function (shape, atype)
   return ArrayCtor(shape, atype)
end

return _M
