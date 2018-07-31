-- Gkyl ------------------------------------------------------------------------
--
-- Non-uniform Cartesian grids
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"

-- Gkyl libraries
local DecompRegionCalc = require "Lib.CartDecomp"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local RectCart = require "Grid.RectCart"

-- NonUniformRectCartGrid ------------------------------------------------------
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

local NonUniformRectCart = Proto(RectCart) -- extends RectCart

function NonUniformRectCart:init(tbl)
   NonUniformRectCart.super.init(self, tbl)

   local ndim = self._ndim
   -- set grid index to first cell in domain
   for d = 1, ndim do
      self._currIdx[d] = 1;
   end
   
   self._nodeCoords = {} -- nodal coordinates in each direction      
   -- initialize nodes to be uniform (will be over-written if mappings are provided)
   for d = 1, ndim do
      -- allocate space: one node extra than cells
      local v =  Lin.Vec(self._numCells[d]+1)
      fillWithUniformNodeCoords(
	 self._lower[d], self._upper[d], self._numCells[d], v)
      self._nodeCoords[d] = v
   end
   
   -- compute nodal coordinates
   if tbl.mappings then
      self.mappings = tbl.mappings
      -- loop over mapping functions, using them to set nodal coordinates
      for d, mapFunc in next, tbl.mappings, nil do
	 if d > ndim then break end -- break out if too many functions provided
	 fillWithMappedNodeCoords(
	    self._lower[d], self._upper[d], self._numCells[d],
	    mapFunc, self._nodeCoords[d])
      end
   end
end

-- member functions
function NonUniformRectCart:id() return "nonuniform" end
function NonUniformRectCart:lower(dir) return self._nodeCoords[dir][1] end
function NonUniformRectCart:upper(dir) return self._nodeCoords[dir][self:numCells(dir)+1] end
function NonUniformRectCart:nodeCoords(dir) return self._nodeCoords[dir] end
function NonUniformRectCart:dx(dir)
   local nodeCoords, idx = self:nodeCoords(dir), self._currIdx
   return nodeCoords[idx[dir]+1]-nodeCoords[idx[dir]]
end
function NonUniformRectCart:cellCenterInDir(d)
   local nodeCoords = self:nodeCoords(d)
   return 0.5*(nodeCoords[idx[d]+1]+nodeCoords[idx[d]])
end
function NonUniformRectCart:cellCenter(xc)
   local idx = self._currIdx
   for d = 1, self._ndim do
      local nodeCoords = self:nodeCoords(d)
      xc[d] = 0.5*(nodeCoords[idx[d]+1]+nodeCoords[idx[d]])
   end
end
function NonUniformRectCart:cellVolume()
   local v = 1.0
   for d = 1, self:ndim() do
      v = v*self:dx(d)
   end
   return v
end
function NonUniformRectCart:getMappings(dir)
   if dir == nil then
      return self.mappings
   else 
      return self.mappings[dir]
   end
end

return NonUniformRectCart
