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

-- Array object is a data containter for storing multi-dimensional
-- data. The interface to this object is delibrately kept low-level
-- and is accessible via the C-struct definition of the array
-- directly. This is done so as to make communication with C-code
-- (including CUDA) easy and via a simple interface.
--
-- Arrays are reference counted. However, for this to work one must be
-- disciplined. If you need to store a reference in another object or
-- C-struct, please obtain it with the aquire() method. When you are
-- done with the pointer release it using the release()
-- method. Release can be automated by proper use of __gc method in
-- the struct metamethods. If you do not do this memory can leak and,
-- very seriously, can disappear under you.
--
-- In general, no ordering is imposed on the data. Ordering needs to
-- be provided elsewhere. For example, there is no mechanism to loop
-- over the array in row/col major order. If one needs to do this one
-- can use the Range object to impose an order and do looping.

-- Array structure. The "d" pointer points to the shape information
-- and also the data stored in the array. d[0] ... d[r-1] store the
-- shape of the array and the rest of the elements of 'd' store the
-- actual stored data. For accessing this data the pointer (d+r) needs
-- to be cast to the proper C-pointer.

ffi.cdef [[
  typedef struct {
    long r, n, c, s, *d; // rank, size, use-count, sizeof element, data
  } GkylArray_t;
]]

-- Various types for arrays of basic C-types
_M.float = { typeof('float'), typeof('float*') }
_M.double = { typeof('double'), typeof('double*') }
_M.int = { typeof('int'), typeof('int*') }
_M.long = { typeof('long'), typeof('long*') }
_M.char = { typeof('char'), typeof('char*') }

-- Array ctype
local ArrayCt = typeof("GkylArray_t")

-- compute size in bytes to allocate
local calcSize = function (r, N, sz)
   return r*sizeof("long") + N*sz
end

local array_fn = {
   clone = function (self)
      local a = new(self)
      local totalSz = calcSize(self.r, self.n, self.s)
      a.d = Alloc.malloc(totalSz)
      a.r, a.n, a.s = self.r, self.n, self.s
      a.c = 1 -- newly created
      ffi.copy(a.d, self.d, totalSz)
      return a
   end,
   aquire = function (self)
      self.c = self.c+1
      return self
   end,
   release = function (self)
      self.c = self.c-1
      if self.c == 0 then
	 Alloc.free(self.d);
	 self.d = 0; self.p = 0
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
      a.c, a.s = 1, sizeof(atype)
      -- allocate memory
      local totalSz = calcSize(a.r, a.n, a.s)
      a.d = Alloc.malloc(totalSz); ffi.fill(a.d, totalSz)
      -- set shape of array
      for i = 1, #shape do a.d[i-1] = shape[i] end
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
   if type(atype) == "table" then
      return ArrayCtor(shape, atype[1])
   else
      return ArrayCtor(shape, atype)
   end
end

-- Get data-pointer to data in array, case to proper C-type
_M.dataPtr = function (arr, atype)
   return ffi.cast(atype[2], arr.d+arr.r)
end

return _M