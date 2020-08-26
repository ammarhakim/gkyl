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
local DataStruct = require "DataStruct"
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

-- Fill array ncoords such that cell spacing is uniform.
local function fillWithUniformNodeCoords(lower, upper, ncell, vcoords)
   local dx = (upper-lower)/ncell
   for n = 1, ncell+1 do
      vcoords[n] = lower+dx*(n-1)
   end
end

-- Compute nodal coordinates using mapping function.
local function fillWithMappedNodeCoords(lower, upper, ncell, mapFunc, vcoords)
   local dx = (upper-lower)/ncell
   for n = 1, ncell+1 do
      vcoords[n] = mapFunc(lower+dx*(n-1))
   end
end

local NonUniformRectCart = Proto(RectCart)   -- Extends RectCart.

function NonUniformRectCart:init(tbl)
   NonUniformRectCart.super.init(self, tbl)

   local ndim = self._ndim
   -- Set grid index to first cell in domain.
   for d = 1, ndim do self._currIdx[d] = 1 end

   self._mappings = tbl.mappings
   
   self._nodeCoords = {}   -- Nodal coordinates in each direction.
   -- Initialize nodes to be uniform (will be over-written if mappings are provided).
   for d = 1, ndim do
      -- Allocate space: one node extra than cells.
      local v =  Lin.Vec(self._numCells[d]+1)
      fillWithUniformNodeCoords(
         self._lower[d], self._upper[d], self._numCells[d], v)
      self._nodeCoords[d] = v
   end
   
   -- Compute nodal coordinates.
   if tbl.mappings then
      -- Loop over mapping functions, using them to set nodal coordinates.
      for d, mapFunc in next, tbl.mappings, nil do
         if d > ndim then break end   -- Break out if too many functions provided.
         fillWithMappedNodeCoords(
            self._lower[d], self._upper[d], self._numCells[d],
            mapFunc, self._nodeCoords[d])
      end
   end
end

-- Determine local domain index. This is complicated by the fact that
-- when using MPI-SHM the processes that live in shmComm all have the
-- same domain index. Hence, we need to use rank from nodeComm and
-- broadcast it to everyone in shmComm.
local function getSubDomIndex(nodeComm, shmComm)
   local idx = ffi.new("int[1]")
   if Mpi.Is_comm_valid(nodeComm) then
      idx[0] = Mpi.Comm_rank(nodeComm)+1 -- sub-domains are indexed from 1
   end
   -- send from rank 0 of shmComm to everyone else: this works as rank
   -- 0 of shmComm is contained in nodeComm
   Mpi.Bcast(idx, 1, Mpi.INT, 0, shmComm)
   return idx[0]
end

-- Member functions.
function NonUniformRectCart:id() return "mapped" end
function NonUniformRectCart:lower(dir) return self._nodeCoords[dir][1] end
function NonUniformRectCart:upper(dir) return self._nodeCoords[dir][self:numCells(dir)+1] end
function NonUniformRectCart:nodeCoords(dir) return self._nodeCoords[dir] end
function NonUniformRectCart:dx(dir)
   local nodeCoords, idx = self:nodeCoords(dir), self._currIdx
   return nodeCoords[idx[dir]+1]-nodeCoords[idx[dir]]
end
function NonUniformRectCart:getDx(dxOut) 
   for d = 1, self:ndim() do
      dxOut[d] = self:dx(d)
   end
end
function NonUniformRectCart:cellCenterInDir(dir)
   local nodeCoords, idx = self:nodeCoords(dir), self._currIdx
   return 0.5*(nodeCoords[idx[dir]+1]+nodeCoords[idx[dir]])
end
function NonUniformRectCart:cellLowerInDir(dir)
   return self:cellCenterInDir(dir)-0.5*self:dx(dir)
end
function NonUniformRectCart:cellUpperInDir(dir)
   return self:cellCenterInDir(dir)+0.5*self:dx(dir)
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

function NonUniformRectCart:write(fName, tmStamp, metaData)
   -- Write a file containing the grid node coordinates. 

   -- Create a grid over nodes and a field to store nodal coordinates.
   local cells, lower, upper, dx = {}, {}, {}, {}
   for d = 1, self:ndim() do
      cells[d] = self:numCells(d)+1   -- One more layer of nodes than cells.
      -- This ensures cell-center of nodal grid lie at nodes of original grid.
      lower[d] = self:lower(d) - 0.5*self:dx(d)
      upper[d] = self:upper(d) + 0.5*self:dx(d)
      dx[d] = self:dx(d)
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
   local indexer = nodalCoords:genIndexer()
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
