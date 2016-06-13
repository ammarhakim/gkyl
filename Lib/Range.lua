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

-- Range ----------------------------------------------------------------------
--
-- A range object, representing a ndim integer index set. 
--------------------------------------------------------------------------------

-- Range object
local function new_Range_ct()
   local rct = typeof("struct { uint8_t _ndim; int _lower[6]; int _upper[6]; }")
   local range_mf = {
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
   }
   local range_mt = {
      __new = function (self, lower, upper)
	 -- lower and upper are tables of integers. Indices are inclusive
	 local r = new(rct)
	 r._ndim = #lower
	 for d = 1, #lower do
	    r._lower[d-1] = lower[d]
	    r._upper[d-1] = upper[d]
	 end
	 return r
      end,
      __index = range_mf
   }
   return metatype(rct, range_mt)
end

return {
   Range = new_Range_ct()
}
