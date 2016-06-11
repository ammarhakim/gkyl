-- Gkyl ------------------------------------------------------------------------
--
-- Cartesian grids, uniform and non-uniform
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries
local Lin = require "Lib.Linalg"

-- define C interfaces
ffi.cdef [[

/* Uniform cartesian grid */
typedef struct {
    uint8_t _ndim; 
    double _lower[6], _upper[6]; 
    int32_t _numCells[6];
    int32_t _currIdx[6];
} RectCartGrid_t;  

]]

-- RectCart ----------------------------------------------------------------------
--
-- A uniform Cartesian grid
--------------------------------------------------------------------------------

local function new_RectCartGrid_ct()
   local uc_type = typeof("RectCartGrid_t")
   local uni_cart_mf = {
      ndim = function (self)
	 return self._ndim
      end,
      lower = function (self, dir)
	 return self._lower[dir]
      end,
      upper = function (self, dir)
	 return self._upper[dir]
      end,
      numCells = function (self, dir)
	 return self._numCells[dir]
      end,
      setIndex = function(self, idx)
	 for d = 1, self._ndim do
	    self._currIdx[d-1] = idx[d]
	 end
      end,
      dx = function (self, dir)
	 return (self._upper[dir]-self._lower[dir])/self._numCells[dir]
      end,
      cellVolume = function (self)
	 local v = 1.0
	 for i = 1, self._ndim do
	    v = v*self:dx(i-1)
	 end
	 return v
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
-- create object factory for uniform cartesian grids
local RectCart = new_RectCartGrid_ct()

-- NonUniformRectCartGrid ----------------------------------------------------------------------
--
-- Stores the nodal coordinates of each node in 1D vectors. The
-- methods names are the same as the uniform cartesian grid object to
-- allow transparent use of uniform meshes in updaters which work for
-- general non-uniform meshes.
--------------------------------------------------------------------------------

-- Fill array ncoords such that cell spacing is uniform
local function fillWithUniformNodeCoords(lower, upper, ncell, vcoords)
   local dx = (upper-lower)/ncell
   for n = 1, ncell+1 do
      vcoords[n] = lower+dx*(n-1)
   end
end

-- Compute nodal coordinates using mapping function
local function fillWithMappedNodeCoords(lower, upper, ncell, mapFunc, vcoords)
   local dx = (upper-lower)/ncell
   for n = 1, ncell+1 do
      vcoords[n] = mapFunc(lower+dx*(n-1))
   end
end

local function new_NonUniformRectCartGrid_ct()
   local nonUni_type = {}
   -- methods 
   local nonUni_mf = {
      ndim = function (self)
	 return self._compGrid._ndim
      end,
      lower = function (self, dir)
	 return self._nodeCoords[dir+1][1]
      end,
      upper = function (self, dir)
	 return self._nodeCoords[dir+1][self:numCells(dir)+1]
      end,
      numCells = function (self, dir)
	 return self._compGrid:numCells(dir)
      end,
      nodeCoords = function (self, dir)
	 return self._nodeCoords[dir+1]
      end,      
      setIndex = function(self, idx)
	 for d = 1, self:ndim() do
	    self._compGrid._currIdx[d-1] = idx[d]
	 end
      end,
      dx = function (self, dir)
	 local nodeCoords, idx = self:nodeCoords(dir), self._compGrid._currIdx
	 return nodeCoords[idx[dir]+1]-nodeCoords[idx[dir]]
      end,
      cellVolume = function (self)
	 local v = 1.0
	 for i = 1, self:ndim() do
	    v = v*self:dx(i-1)
	 end
	 return v
      end,
   }
   nonUni_type.__index = nonUni_mf

   -- constructor to make a new non-uniform grid
   function nonUni_type:new(tbl)
      local self = setmetatable({}, nonUni_type)
      self._compGrid = RectCart(tbl) -- computational space grid

      local ndim = self._compGrid._ndim
      -- set grid index to first cell in domain
      for d = 0, ndim-1 do
	 self._compGrid._currIdx[d] = 1;
      end

      self._nodeCoords = {} -- nodal coordinates in each direction      
      -- initialize nodes to be uniform (will be over-written if mappings are provided)
      for d = 0, ndim-1 do
	 -- allocate space: one node extra than cells
	 local v =  Lin.vec(self._compGrid:numCells(d)+1)
	 fillWithUniformNodeCoords(
	    self._compGrid._lower[d], self._compGrid._upper[d], self._compGrid:numCells(d), v)
	 self._nodeCoords[d+1] = v
      end

      -- compute nodal coordinates
      if tbl.mappings then
	 -- loop over mapping functions, using them to set nodal coordinates
	 for i, mapFunc in next, tbl.mappings, nil do
	    if i > ndim then break end -- break out if too many functions provided
	    local d = i-1 -- directions start from 0
	    fillWithMappedNodeCoords(
	       self._compGrid._lower[d], self._compGrid._upper[d], self._compGrid:numCells(d),
	       mapFunc, self._nodeCoords[i])
	 end
      end
      
      return self
   end

   -- make object callable, and redirect call to the :new method
   setmetatable (nonUni_type, {
		    __call = function (self, o)
		       return self.new(self, o)
		    end,
   })

   return nonUni_type
end
-- create object factory for non-uniform cartesian grids
local NonUniformRectCartGrid = new_NonUniformRectCartGrid_ct()

return {
   RectCart = RectCart,
   NonUniformRectCartGrid = NonUniformRectCartGrid,
}
