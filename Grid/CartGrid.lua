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
local Range = require "Lib.Range"
local Mpi = require "Comm.Mpi"

-- RectCart --------------------------------------------------------------------
--
-- A uniform Cartesian grid
--------------------------------------------------------------------------------

-- Meta-type of base-cartesian grid: both uniform and non-uniform
-- grids use this object
local uni_cart_mt = {
   __new = function (self, tbl)
      local lo = tbl.lower or {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
      local up = tbl.upper or {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}
      local cells = tbl.cells
      local g = new(self)
      g._ndim = #cells
      for d = 1, #cells do
	 g._lower[d-1] = lo[d]
	 g._upper[d-1] = up[d]
	 g._numCells[d-1] = cells[d]
      end

      local l, u = {}, {}
      for d = 1, #cells do
	 l[d], u[d] = 1, cells[d]
      end
      g._globalRange = Range.Range(l, u)
      g._localRange = Range.Range(l, u)
      if tbl.decomposition then
	 local comm = tbl.decomposition:comm()
	 local localRange = tbl.decomposition:subDomain(Mpi.Comm_rank(self._comm))
	 g._localRange:copy(localRange)
      end

      return g
   end,
   __index =  {
      ndim = function (self)
	 return self._ndim
      end,
      lower = function (self, dir)
	 return self._lower[dir-1]
      end,
      upper = function (self, dir)
	 return self._upper[dir-1]
      end,
      numCells = function (self, dir)
	 return self._numCells[dir-1]
      end,
      localRange = function (self)
	 return self._localRange
      end,
      globalRange = function (self)
	 return self._globalRange
      end,      
      setIndex = function(self, idx)
	 for d = 1, self._ndim do
	    self._currIdx[d-1] = idx[d]
	 end
      end,
      dx = function (self, dir)
	 return (self:upper(dir)-self:lower(dir))/self:numCells(dir)
      end,
      cellVolume = function (self)
	 local v = 1.0
	 for i = 1, self._ndim do
	    v = v*self:dx(i)
	 end
	 return v
      end,
   }
}

-- Uniform Cartesian grid
local RectCart = {}
-- constructor to make a new uniform grid
function RectCart:new(tbl)
   local self = setmetatable({}, RectCart)
   local cells = tbl.cells
   self._ndim = #cells -- dimensions
   -- default grid extents
   local lo = tbl.lower or {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
   local up = tbl.upper or {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}

   self._lower = new("int[?]", self._ndim)
   self._upper = new("int[?]", self._ndim)
   self._numCells = new("int[?]", self._ndim)
   self._currIdx = new("int[?]", self._ndim)
   
   for d = 1, #cells do
      self._lower[d-1], self._upper[d-1] = lo[d], up[d]
      self._numCells[d-1] = cells[d]
   end
   -- compute global range
   local l, u = {}, {}
   for d = 1, #cells do
      l[d], u[d] = 1, cells[d]
   end
   self._globalRange = Range.Range(l, u)   
   self._localRange = Range.Range(l, u)

   local decomp = tbl.decomposition and tbl.decomposition or nil  -- decomposition
   if decomp then
      -- in parallel, we need to adjust local range      
      self._comm = decomp:comm()
      self._decomposedRange = decomp:decompose(self._globalRange)
      local subDomIdx = Mpi.Comm_rank(self._comm)+1 -- sub-domains are indexed from 1
      local localRange = self._decomposedRange:subDomain(subDomIdx)
      self._localRange:copy(localRange)
   end
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(RectCart, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
RectCart.__index = {
   comm = function (self)
      return self._comm
   end,
   ndim = function (self)
      return self._ndim
   end,
   lower = function (self, dir)
      return self._lower[dir-1]
   end,
   upper = function (self, dir)
      return self._upper[dir-1]
   end,
   numCells = function (self, dir)
      return self._numCells[dir-1]
   end,
   localRange = function (self)
      return self._localRange
   end,
   globalRange = function (self)
      return self._globalRange
   end,
   setIndex = function(self, idx)
      for d = 1, self._ndim do
	 self._currIdx[d-1] = idx[d]
      end
   end,
   dx = function (self, dir)
      return (self:upper(dir)-self:lower(dir))/self:numCells(dir)
   end,
   cellVolume = function (self)
      local v = 1.0
      for i = 1, self._ndim do
	 v = v*self:dx(i)
      end
      return v
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

   local cells = tbl.cells
   self._ndim = #cells -- dimensions
   -- default grid extents
   local lo = tbl.lower or {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
   local up = tbl.upper or {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}

   self._lower = new("int[?]", self._ndim)
   self._upper = new("int[?]", self._ndim)
   self._numCells = new("int[?]", self._ndim)
   self._currIdx = new("int[?]", self._ndim)
   
   for d = 1, #cells do
      self._lower[d-1], self._upper[d-1] = lo[d], up[d]
      self._numCells[d-1] = cells[d]
   end
   -- compute global range
   local l, u = {}, {}
   for d = 1, #cells do
      l[d], u[d] = 1, cells[d]
   end
   self._globalRange = Range.Range(l, u)   
   self._localRange = Range.Range(l, u)

   local decomp = tbl.decomposition and tbl.decomposition or nil  -- decomposition
   if decomp then
      -- in parallel, we need to adjust local range
      self._decomposedRange = decomp:decompose(self._globalRange)
      local subDomIdx = Mpi.Comm_rank(self._comm)+1 -- sub-domains are indexed from 1
      local localRange = self._decomposedRange:subDomain(subDomIdx)
      self._localRange:copy(localRange)
   end   

   local ndim = self._ndim
   -- set grid index to first cell in domain
   for d = 1, ndim do
      self._currIdx[d-1] = 1;
   end
   
   self._nodeCoords = {} -- nodal coordinates in each direction      
   -- initialize nodes to be uniform (will be over-written if mappings are provided)
   for d = 1, ndim do
      -- allocate space: one node extra than cells
      local v =  Lin.Vec(self._numCells[d-1]+1)
      fillWithUniformNodeCoords(
	 self._lower[d-1], self._upper[d-1], self._numCells[d-1], v)
      self._nodeCoords[d] = v
   end
   
   -- compute nodal coordinates
   if tbl.mappings then
      -- loop over mapping functions, using them to set nodal coordinates
      for d, mapFunc in next, tbl.mappings, nil do
	 if d > ndim then break end -- break out if too many functions provided
	 fillWithMappedNodeCoords(
	    self._lower[d-1], self._upper[d-1], self._numCells[d-1],
	    mapFunc, self._nodeCoords[d])
      end
   end
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable (NonUniformRectCart, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
NonUniformRectCart.__index = {
   comm = function (self)
      return self._comm
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
      return self._numCells[dir-1]
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
   setIndex = function(self, idx)
      for d = 1, self:ndim() do
	 self._currIdx[d-1] = idx[d]
      end
   end,
   dx = function (self, dir)
      local nodeCoords, idx = self:nodeCoords(dir), self._currIdx
      return nodeCoords[idx[dir-1]+1]-nodeCoords[idx[dir-1]]
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
