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
local DecompRegionCalc = require "Lib.CartDecomp"
local Range = require "Lib.Range"
local Mpi = require "Comm.Mpi"

-- Determine local domain index. This is complicated by the fact that
-- when using MPI-SHM the processes that live in shmComm all have the
-- same domain index. Hence, we need to use rank from nodeComm and
-- broadcast it to everyone in shmComm.
local function getSubDomIndex(nodeComm, shmComm)
   local idx = new("int[1]")
   if Mpi.Is_comm_valid(nodeComm) then
      idx[0] = Mpi.Comm_rank(nodeComm)+1 -- sub-domains are indexed from 1
   end
   -- send from rank 0 of shmComm to everyone else: this works as rank
   -- 0 of shmComm is contained in nodeComm
   Mpi.Bcast(idx, 1, Mpi.INT, 0, shmComm)
   return idx[0]
end

-- Method to initialize basic elements in a cartesian grid. Used by
-- uniform and non-uniform grid objects
local function initCartGrid(self, tbl)
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
   else
      -- create a dummy decomp and use it to set the grid
      local cuts = {}
      for i = 1, self._ndim do cuts[i] = 1 end
      local dec1 = DecompRegionCalc.CartProd { cuts = cuts, useShared = true }
      self._commSet = dec1:commSet()
      self._decomposedRange = dec1:decompose(self._globalRange)
      self._block = 1
   end
end

-- RectCart --------------------------------------------------------------------
--
-- A uniform Cartesian grid
--------------------------------------------------------------------------------

-- Uniform Cartesian grid
local RectCart = {}
-- constructor to make a new uniform grid
function RectCart:new(tbl)
   local self = setmetatable({}, RectCart)
   initCartGrid(self, tbl)
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(RectCart, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
RectCart.__index = {
   commSet = function (self)
      return self._commSet
   end,
   isShared = function (self)
      return self._isShared
   end,
   subGridId = function (self)
      return self._block
   end,
   numSharedProcs = function (self)
      return Mpi.Comm_size(self._commSet.sharedComm)
   end,
   decomposedRange = function (self)
      return self._decomposedRange
   end,
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
   localRange = function (self)
      return self._localRange
   end,
   globalRange = function (self)
      return self._globalRange
   end,
   isDirPeriodic = function (self, dir)
      return self._isDirPeriodic[dir]
   end,
   setIndex = function(self, idx)
      for d = 1, self._ndim do
      	 self._currIdx[d] = idx[d]
      end
   end,
   dx = function (self, dir)
      return self._dx[dir]
   end,
   cellCenter = function (self, xc)
      for d = 1, self._ndim do
	 xc[d] = self:lower(d) + (self._currIdx[d]-0.5)*self:dx(d)
      end
   end,
   cellVolume = function (self)
      return self._vol
   end,
}

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

local NonUniformRectCart = {}
-- constructor to make a new non-uniform grid
function NonUniformRectCart:new(tbl)
   local self = setmetatable({}, NonUniformRectCart)

   initCartGrid(self, tbl) -- initialize basic stuff

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
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable (NonUniformRectCart, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
NonUniformRectCart.__index = {
   commSet = function (self)
      return self._commSet
   end,
   isShared = function (self)
      return self._isShared
   end,   
   subGridId = function (self)
      return self._block
   end,
   numSharedProcs = function (self)
      return Mpi.Comm_size(self._commSet.sharedComm)
   end,
   decomposedRange = function (self)
      return self._decomposedRange
   end,
   ndim = function (self)
      return self._ndim
   end,
   lower = function (self, dir)
      return self._nodeCoords[dir][1]
   end,
   upper = function (self, dir)
      return self._nodeCoords[dir][self:numCells(dir)+1]
   end,
   numCells = function (self, dir)
      return self._numCells[dir]
   end,
   nodeCoords = function (self, dir)
      return self._nodeCoords[dir]
   end,
   localRange = function (self)
      return self._localRange
   end,
   globalRange = function (self)
      return self._globalRange
   end,
   isDirPeriodic = function (self, dir)
      return self._isDirPeriodic[dir]
   end,
   setIndex = function(self, idx)
      for d = 1, self:ndim() do
	 self._currIdx[d] = idx[d]
      end
   end,
   dx = function (self, dir)
      local nodeCoords, idx = self:nodeCoords(dir), self._currIdx
      return nodeCoords[idx[dir]+1]-nodeCoords[idx[dir]]
   end,
   cellCenter = function (self, xc)
      local idx = self._currIdx
      for d = 1, self._ndim do
	 local nodeCoords = self:nodeCoords(d)
	 xc[d] = 0.5*(nodeCoords[idx[d]+1]+nodeCoords[idx[d]])
      end
   end,
   cellVolume = function (self)
      local v = 1.0
      for d = 1, self:ndim() do
	 v = v*self:dx(d)
      end
      return v
   end,
}

return {
   RectCart = RectCart,
   NonUniformRectCart = NonUniformRectCart,
}
