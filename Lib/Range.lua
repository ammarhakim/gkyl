-- Gkyl ------------------------------------------------------------------------
--
-- Range object
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi  = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries
local Lin = require "Lib.Linalg"

_M = {}

-- Range ----------------------------------------------------------------------
--
-- A range object, representing a ndim integer index set. 
--------------------------------------------------------------------------------

ffi.cdef [[ typedef struct { uint8_t _ndim; int _lower[6]; int _upper[6]; } Range_t; ]]

-- generic iterator function creator: only difference between row- and
-- col-major order is the order in which the indices are incremented
local function make_range_iter(si, ei, incr)
   return function (iterState)
      if iterState.isEmpty then return nil end -- nothing to do for empty range
      
      local idx = iterState.currIdx
      if iterState.isFirst then
	 -- first time just return start index
	 iterState.isFirst = false
	 return idx
      end
      
      local range = iterState.range
      local ndim = range:ndim()

      for dir = si, ei, incr do
	 -- bump index till we can bump no more
	 idx[dir] = idx[dir]+1
	 if idx[dir] > range:upper(dir) then
	    -- run out of indices to bump, reset
	    idx[dir] = range:lower(dir)
	 else
	    return idx
	 end
      end
   end
end

-- Range object meta-type
local range_mt = {
   __new = function (self, lower, upper)
      -- lower and upper are tables of integers. Indices are inclusive
      local r = new(self)
      r._ndim = #lower
      for d = 1, #lower do
	 r._lower[d-1] = lower[d]
	 r._upper[d-1] = upper[d]
      end
      for d = 1, #lower do
	 -- adjust to give zero volume range if upper is less than lower
	 if r._upper[d-1] < r._lower[d-1] then
	    r._upper[d-1] = r._lower[d-1]-1
	 end
      end
      return r
   end,
   __index = {
      ndim = function (self)
	 return self._ndim
      end,
      lower = function (self, dir)
	 return self._lower[dir-1]
      end,
      upper = function (self, dir)
	 return self._upper[dir-1]
      end,
      shape = function (self, dir)
	 return self._upper[dir-1]-self._lower[dir-1]+1
      end,
      volume = function (self)
	 local v = 1
	 for dir = 1, self._ndim do
	    v = v*self:shape(dir)
	 end
	 return v
      end,
      extend = function (self, lExt, uExt)
	 local r = new(typeof("Range_t"))
	 r._ndim = self:ndim()
	 for dir = 1, self:ndim() do
	    r._lower[dir-1], r._upper[dir-1] = self:lower(dir)-lExt, self:upper(dir)+uExt
	 end
	 return r
      end,
      _iter = function (self, iter_func)
	 local iterState = { isFirst = true, isEmpty = self:volume() == 0 and true or false }
	 iterState.currIdx = Lin.IntVec(self:ndim())
	 for dir = 1, self:ndim() do
	    iterState.currIdx[dir] = self:lower(dir)
	 end
	 iterState.range = self
	 return iter_func, iterState
      end,
      colMajorIter = function (self)
	 return self:_iter( make_range_iter(1, self:ndim(), 1) )
      end,
      rowMajorIter = function (self)
	 return self:_iter( make_range_iter(self:ndim(), 1, -1) )
      end,      
   }
}
-- construct Range object, attaching meta-type to it
_M.Range = metatype(typeof("Range_t"), range_mt)

-- Indexers --------------------------------------------------------------------
--
-- Linear indexers, mapping n-dimensional index to a linear index
--------------------------------------------------------------------------------

local function getIndex1(ac, i1)
   return ac[0] + i1*ac[1]
end
local function getIndex2(ac, i1, i2)
   return ac[0]+i1*ac[1]+i2*ac[2]
end
local function getIndex3(ac, i1, i2, i3)
   return ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]
end
local function getIndex4(ac, i1, i2, i3, i4)
   return ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]
end
local function getIndex5(ac, i1, i2, i3, i4, i5)
   return ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]+i5*ac[5]
end
local function getIndex6(ac, i1, i2, i3, i4, i5, i6)
   return ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]+i5*ac[5]+i6*ac[6]
end
local function getIndex7(ac, i1, i2, i3, i4, i5, i6, i7)
   return ac[0]+i1*ac[1]+i2*ac[2]+i3*ac[3]+i4*ac[4]+i5*ac[5]+i6*ac[6]+i7*ac[7]
end
-- package these up into a table
_M.indexerFunctions = {getIndex1, getIndex2, getIndex3, getIndex4, getIndex5, getIndex6, getIndex7}

-- create a row-major indexer  given "range" object
function _M.makeRowMajorIndexer(range)
   local ac = new("double[7]")
   local ndim = range:ndim()
   ac[ndim] = 1
   for i = ndim-1, 1, -1 do
      ac[i] = ac[i+1]*range:shape(i+1)
   end
   local start = 0
   for i = 1, ndim do
      start = start + ac[i]*range:lower(i)
   end
   ac[0] = -start
   return ac
end

-- create a column-major indexer  given "range" object
function _M.makeColMajorIndexer(range)
   local ac = new("double[7]")
   local ndim = range:ndim()
   ac[1] = 1
   for i = 2, ndim do
      ac[i] = ac[i-1]*range:shape(i-1)
   end
   local start = 0
   for i = 1, ndim do
      start = start + ac[i]*range:lower(i)
   end
   ac[0] = -start
   return ac
end

return _M
