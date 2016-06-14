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

-- Range ----------------------------------------------------------------------
--
-- A range object, representing a ndim integer index set. 
--------------------------------------------------------------------------------

-- generic iterator function creator: only difference between row- and
-- col-major order is the order in which the indices are incremented
local function make_range_iter(si, ei, incr)
   return function (iterState)
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
      _iter = function (self, iter_func)
	 local iterState = { isFirst = true }
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
local rct = typeof("struct { uint8_t _ndim; int _lower[6]; int _upper[6]; }")
local Range = metatype(rct, range_mt)

return {
   Range = Range
}
