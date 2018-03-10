-- Gkyl ------------------------------------------------------------------------
--
-- Cartesian grids, uniform and non-uniform
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

-- RectCart --------------------------------------------------------------------
--
-- A uniform Cartesian grid
--------------------------------------------------------------------------------

-- Uniform Cartesian grid
local RectCart = Proto()

function RectCart:init(tbl)
   local cells = tbl.cells
   self._ndim = #cells -- dimensions
   -- default grid extents
   local lo = tbl.lower and tbl.lower or {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
   local up = tbl.upper and tbl.upper or {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}

   -- periodic directions
   self._periodicDirs = tbl.periodicDirs and tbl.periodicDirs or {}
   -- construct table to indicate if directions are periodic  (default no directions are periodic)
   self._isDirPeriodic = {false, false, false, false, false, false}
   for _, dir in ipairs(self._periodicDirs) do
      self._isDirPeriodic[dir] = true
   end

   -- set stuff up
   self._lower = Lin.Vec(self._ndim)
   self._upper = Lin.Vec(self._ndim)
   self._dx = Lin.Vec(self._ndim)
   self._numCells = Lin.IntVec(self._ndim)
   self._currIdx = Lin.IntVec(self._ndim)

   self._vol = 1.0
   for d = 1, #cells do
      self._lower[d], self._upper[d] = lo[d], up[d]
      self._numCells[d] = cells[d]
      self._dx[d] = (up[d]-lo[d])/cells[d]
      self._vol = self._vol*self._dx[d]
   end

   -- compute global range
   local l, u = {}, {}
   for d = 1, #cells do
      l[d], u[d] = 1, cells[d]
   end
   self._globalRange = Range.Range(l, u)   
   self._localRange = Range.Range(l, u)
   self._block = 1 -- block number for use in parallel communications
   self._isShared = false -- is grid shared
   
   local decomp = tbl.decomposition and tbl.decomposition or nil  -- decomposition
   if decomp then
      assert(decomp:ndim() == self._ndim, "Decomposition dimensions must be same as grid dimensions!")

      self._isShared = decomp:isShared()
      -- in parallel, we need to adjust local range      
      self._commSet = decomp:commSet()
      self._decomposedRange = decomp:decompose(self._globalRange)
      local subDomIdx = getSubDomIndex(self._commSet.nodeComm, self._commSet.sharedComm)
      self._block = subDomIdx
      local localRange = self._decomposedRange:subDomain(subDomIdx)
      self._localRange:copy(localRange)
      self._cuts = {}
      for i = 1, self._ndim do self._cuts[i] = decomp:cuts(i) end
   else
      -- create a dummy decomp and use it to set the grid
      local cuts = {}
      for i = 1, self._ndim do cuts[i] = 1 end
      local dec1 = DecompRegionCalc.CartProd { cuts = cuts, useShared = true }
      self._commSet = dec1:commSet()
      self._decomposedRange = dec1:decompose(self._globalRange)
      self._block = 1
      self._cuts = cuts
   end
end

-- member functions
function RectCart:commSet() return self._commSet end 
function RectCart:isShared() return self._isShared end
function RectCart:subGridId() return self._block end
function RectCart:numSharedProcs() return Mpi.Comm_size(self._commSet.sharedComm) end
function RectCart:decomposedRange() return self._decomposedRange end
function RectCart:ndim() return self._ndim end
function RectCart:lower(dir) return self._lower[dir] end
function RectCart:upper(dir) return self._upper[dir] end
function RectCart:numCells(dir) return self._numCells[dir] end
function RectCart:localRange() return self._localRange end
function RectCart:globalRange() return self._globalRange end
function RectCart:isDirPeriodic(dir) return self._isDirPeriodic[dir] end
function RectCart:cuts(dir) return self._cuts[dir] end
function RectCart:setIndex(idx)
   for d = 1, self._ndim do
      self._currIdx[d] = idx[d]
   end
end
function RectCart:dx(dir) return self._dx[dir] end
function RectCart:cellCenterInDir(d)
   return self:lower(d) + (self._currIdx[d]-0.5)*self:dx(d)
end
function RectCart:cellCenter(xc)
   for d = 1, self._ndim do
      xc[d] = self:lower(d) + (self._currIdx[d]-0.5)*self:dx(d)
   end
end
function RectCart:cellVolume() return self._vol end

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

local NonUniformRectCart = Proto(RectCart) -- extends RectCart
-- ctor
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

return {
   RectCart = RectCart,
   NonUniformRectCart = NonUniformRectCart,
}
