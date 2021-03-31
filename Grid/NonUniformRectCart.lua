-- Gkyl ------------------------------------------------------------------------
--
-- Non-uniform Cartesian grids.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- System libraries.
local ffi = require "ffi"

-- Gkyl libraries.
local DataStruct       = require "DataStruct"
local DecompRegionCalc = require "Lib.CartDecomp"
local Lin              = require "Lib.Linalg"
local Mpi              = require "Comm.Mpi"
local Proto            = require "Lib.Proto"
local Range            = require "Lib.Range"
local RectCart         = require "Grid.RectCart"

-- NonUniformRectCartGrid ------------------------------------------------------
--
-- Stores the nodal coordinates of each node in 1D vectors. The
-- methods names are the same as the uniform cartesian grid object to
-- allow transparent use of uniform meshes in updaters which work for
-- general non-uniform meshes.
--------------------------------------------------------------------------------

-- Fill array ncoords such that cell spacing is uniform.
local function fillWithUniformNodeCoords(lower, upper, ncell, vcoords)
   local dx = (upper-lower)/ncell
   for n = 1, ncell+1 do vcoords[n] = lower+dx*(n-1) end
end

-- Compute nodal coordinates using mapping function.
local function fillWithMappedNodeCoords(lower, upper, ncell, mapFunc, vcoords)
   local dx = (upper-lower)/ncell
   for n = 1, ncell+1 do vcoords[n] = mapFunc(lower+dx*(n-1)) end
end

local NonUniformRectCart = Proto(RectCart)   -- Extends RectCart.

function NonUniformRectCart:init(tbl)
   NonUniformRectCart.super.init(self, tbl)

   local ndim = self._ndim
   -- Set grid index to first cell in domain.
   for d = 1, ndim do self._currIdx[d] = 1 end

   self._nodeCoords = {}   -- Nodal coordinates in each direction.
   -- Initialize nodes to be uniform (will be over-written if mappings are provided).
   for d = 1, ndim do
      -- Allocate space: one node extra than cells.
      local v = Lin.Vec(self._numCells[d]+1)
      fillWithUniformNodeCoords(
         self._lower[d], self._upper[d], self._numCells[d], v)
      self._nodeCoords[d] = v
   end
   
   -- Compute nodal coordinates.
   if tbl.mappings then
      self.mappings = tbl.mappings
      -- Loop over mapping functions, using them to set nodal coordinates.
      for d, mapFunc in next, tbl.mappings, nil do
         if d > ndim then break end   -- Break out if too many functions provided.
         fillWithMappedNodeCoords(
            self._lower[d], self._upper[d], self._numCells[d],
            mapFunc, self._nodeCoords[d])
      end
   end

   self._gridVol = 1.0
   for d = 1, self._ndim do
      self._gridVol = self._gridVol*(self._nodeCoords[d][self._numCells[d]+1]-self._nodeCoords[d][1])
   end
   
   -- Precompute cell length, center and boundary coords.
   self._cellDx, self._cellCenter, self._cellLower, self._cellUpper = {}, {}, {}, {}
   for d = 1, ndim do
      self._cellDx[d]     = Lin.Vec(self._numCells[d])
      self._cellCenter[d] = Lin.Vec(self._numCells[d])
      self._cellLower[d]  = Lin.Vec(self._numCells[d])
      self._cellUpper[d]  = Lin.Vec(self._numCells[d])
      local nodeCoords = self._nodeCoords[d]
      for idx = 1, self._numCells[d] do
         self._cellDx[d][idx]     = nodeCoords[idx+1]-nodeCoords[idx]
         self._cellCenter[d][idx] = 0.5*(nodeCoords[idx+1]+nodeCoords[idx])
         self._cellLower[d][idx]  = self._cellCenter[d][idx]-0.5*self._cellDx[d][idx]
         self._cellUpper[d][idx]  = self._cellCenter[d][idx]+0.5*self._cellDx[d][idx]
      end
   end
end

-- Member functions.
function NonUniformRectCart:id() return "mapped" end
function NonUniformRectCart:lower(dir) return self._nodeCoords[dir][1] end
function NonUniformRectCart:mid(dir) return self:lower(dir) + (self:upper(dir)-self:lower(dir))/2 end
function NonUniformRectCart:upper(dir) return self._nodeCoords[dir][self:numCells(dir)+1] end
function NonUniformRectCart:logicalLower(dir) return self._lower[dir] end
function NonUniformRectCart:logicalMid(dir) return self._lower[dir] + (self._upper[dir]-self._lower[dir])/2 end
function NonUniformRectCart:logicalUpper(dir) return self._upper[dir] end
function NonUniformRectCart:nodeCoords(dir) return self._nodeCoords[dir] end
function NonUniformRectCart:dx(dir)
   local idxInDir = self._currIdx[dir]
   if idxInDir == 0 then
      idxInDir = 1
   elseif idxInDir > self:localRange():upper(dir) then
      idxInDir = self:localRange():upper(dir)
   end
   return self._cellDx[dir][idxInDir]
end
function NonUniformRectCart:getDx(dxOut) 
   for d = 1, self:ndim() do dxOut[d] = self:dx(d) end
end
function NonUniformRectCart:cellCenterInDir(dir)
   local idxInDir = self._currIdx[dir]
   if idxInDir == 0 then
      idxInDir = 1
   elseif idxInDir > self:localRange():upper(dir) then
      idxInDir = self:localRange():upper(dir)
   end
   return self._cellCenter[dir][idxInDir]
end
function NonUniformRectCart:cellLowerInDir(dir)
   local idxInDir = self._currIdx[dir]
   return self._cellLower[dir][idxInDir]
end
function NonUniformRectCart:cellUpperInDir(dir)
   local idxInDir = self._currIdx[dir]
   return self._cellUpper[dir][idxInDir]
end
function NonUniformRectCart:cellCenter(xc)
   for d = 1, self._ndim do xc[d] = self:cellCenterInDir(d) end
end
function NonUniformRectCart:cellVolume()
   local v = 1.0
   for d = 1, self:ndim() do v = v*self:dx(d) end
   return v
end
function NonUniformRectCart:getMappings(dir)
   if dir == nil then
      return self.mappings
   else 
      return self.mappings[dir]
   end
end

function NonUniformRectCart:write(fName, tmStamp, metaData)
   -- Write a file containing the grid node coordinates. 

   -- Create a grid over nodes and a field to store nodal coordinates.
   local cells, lower, upper, dx = {}, {}, {}, {}
   for d = 1, self:ndim() do
      cells[d] = self:numCells(d)+1   -- One more layer of nodes than cells.
      -- This ensures cell-center of nodal grid lie at nodes of original grid.
      lower[d] = self:lower(d) - 0.5*self:dx(d)
      upper[d] = self:upper(d) + 0.5*self:dx(d)
      dx[d]    = self:dx(d)
   end
   -- WILL NEED TO MAKE THIS WORK IN PARALLEL .. EVENTUALLY
   local grid = RectCart {
      lower = lower,
      upper = upper,
      cells = cells,
   }
   local nodalCoords = DataStruct.Field {
      onGrid = grid,
      numComponents = self:ndim(),
      metaData = metaData
   }

   local xnc, xnp = Lin.Vec(self:ndim()), Lin.Vec(self:ndim())

   local localRange = nodalCoords:localRange()
   local indexer    = nodalCoords:genIndexer()
   for idx in localRange:rowMajorIter() do
      grid:setIndex(idx)
      local nitr = nodalCoords:get(indexer(idx))
      for d = 1, self:ndim() do
         nodeCoords = self:nodeCoords(d)
         nitr[d] = nodeCoords[idx[d]]
      end
   end

   nodalCoords:write(fName, tmStamp)
end

return NonUniformRectCart
