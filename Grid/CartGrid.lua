--------------------------------------------------------------------------------
-- Cartesian grids, uniform and non-uniform
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- define C interfaces
ffi.cdef [[
      typedef struct { 
        uint8_t _ndim; 
        double _lower[6], _upper[6]; 
        int32_t _numCells[6]; 
      } UniformCartGrid_t;
]]

-- Uniform cartesian mesh meta-object creator
local function new_UniformCartGrid_ct()
   local uc_type = typeof("UniformCartGrid_t")
   local uni_cart_mf = {
      ndim = function (self, i)
	 return self._ndim
      end,
      lower = function (self, i)
	 return self._lower[i]
      end,
      upper = function (self, i)
	 return self._upper[i]
      end,
      numCells = function (self, i)
	 return self._numCells[i]
      end,
      dx = function (self, i)
	 return (self._upper[i]-self._lower[i])/self._numCells[i]
      end,
   }
   local uni_cart_mt = {
      __new = function (self, tbl)
	 local lo, up, cells = tbl["lower"], tbl["upper"], tbl["cells"]
	 local g = new(uc_type)
	 g._ndim = #lo
	 for d = 1, #lo do
	    g._lower[d-1] = lo[d]
	    g._upper[d-1] = up[d]
	    g._numCells[d-1] = cells[d]
	 end
	 return g
      end,
      __index = uni_cart_mf,
   }
   return metatype(uc_type, uni_cart_mt)
end

return {
   CartGrid = new_UniformCartGrid_ct(),
}
