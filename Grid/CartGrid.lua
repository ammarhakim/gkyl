-- Gkyl ------------------------------------------------------------------------
--
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

/* Uniform cartesian grid */
typedef struct {
    uint8_t _ndim; 
    double _lower[6], _upper[6]; 
    int32_t _numCells[6]; 
} UniformCartGrid_t;      

/* Node coordinates for use in non-uniform cartesian grid */
typedef struct { uint32_t _n; double _x[?]; } NodeCoord_t;
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
	 local lo = tbl.lower or {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
	 local up = tbl.upper or {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
	 local cells = tbl.cells
	 local g = new(uc_type)
	 g._ndim = #cells
	 for d = 1, #cells do
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

-- Non-uniform cartesian mesh meta-object creator: this is basically
-- an object which stores the nodal coordinates of each node in 1D
-- arrays. The methods are similar to the uniform cartesian grid
-- object to allow transparent use of uniform meshes in updaters which
-- work for general non-uniform meshes.
local function new_NonUniformCartGrid_ct()
end

return {
   CartGrid = new_UniformCartGrid_ct(),
   NonUniformCartGrid = new_NonUniformCartGrid_ct(),
}
