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

-- Array structure
local ArrayStruct = [[
  typedef struct {
    int r;   // rank
    long *s; // shape
    long N;  // total size
    int c;   // number of references
    int sz;  // size of individual elements
    $* d;    // pointer to data
  } GkylArray_t;
]]

-- (Meta) function to create Array constructor for given type
local function Array_meta_ctor(elct)
   local sz = sizeof(elct)

   local array_fn = {
      aquire = function (self)
	 self.c = self.c + 1
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
      __new = function (self, shape)
	 local a = ffi.new(self)
	 a.r = #shape
	 a.s = Alloc.malloc(a.r*size("long"))
	 a.N = 1
	 for i = 1, #shape do
	    a.s[i-1] = shape[i]
	    a.N = a.N*shape[i]
	 end
	 a.c = 1
	 a.sz = sz
	 a.d = Alloc.calloc(a.N, sz)
      end,
      __gc = function(self)
	 self:release()
      end,
      __index = array_fn
   }
   return metatype(typeof("GkylArray_t", elct), array_mt) 
end

return {
   Array_meta_ctor,
   Array = Array_meta_ctor(ffi.typeof("double")),
}
