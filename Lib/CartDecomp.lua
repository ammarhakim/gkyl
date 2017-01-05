-- Gkyl ------------------------------------------------------------------------
--
-- Decomposition algorithms for rectangular grids.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"
local xsys = require "xsys"
local new, copy, fill, sizeof, typeof, metatype = xsys.from(ffi,
     "new, copy, fill, sizeof, typeof, metatype")

-- Gkyl libraries
local Range = require "Lib.Range"
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"

-- create constructor to store vector of Range objects
local RangeVec = Lin.new_vec_ct(ffi.typeof("Range_t"))
-- create constructor to store integer pairs
local PairsVec = Lin.new_vec_ct(ffi.typeof("struct { int32_t lower, upper; } "))

-- DecomposedRange --------------------------------------------------------------
--
-- The decomposition object's decompose() method creates a
-- DecomposedRange object
--------------------------------------------------------------------------------

local DecomposedRange = {}
-- constructor to make object that stores decomposed ranges
function DecomposedRange:new(decomp)
   local self = setmetatable({}, DecomposedRange)

   local ones, upper = {}, {}
   for d = 1, decomp:ndim() do
      ones[d] = 1; upper[d] = decomp:cuts(d)
   end

   self._cutsRange = Range.Range(ones, upper)
   self._comm = decomp._comm   
   self._domains = RangeVec(self._cutsRange:volume())

   -- for periodic directions store "skeleton" sub-domains, i.e. those
   -- that touch the domain boundary
   self._periodicDomPairs = {} -- table of size 'ndim'

   for d = 1, decomp:ndim() do
      local sz = self._cutsRange:shorten(d):volume() -- number of skeleton cells
      self._periodicDomPairs[d] = PairsVec(sz) -- store as integer pair
   end
   
   return self   
end
-- make object callable, and redirect call to the :new method
setmetatable(DecomposedRange, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
DecomposedRange.__index = {
   comm = function (self)
      return self._comm
   end,
   ndim = function (self)
      return self._cutsRange:ndim()
   end,
   cuts = function (self, dir)
      return self._cutsRange:upper(dir)
   end,
   numSubDomains = function (self)
      return self._cutsRange:volume()
   end,
   subDomain = function (self, k)
      return self._domains[k]
   end
}

-- CartProdDecomp --------------------------------------------------------------
--
-- Decomposition of grid into smaller boxes, specified by number of
-- cuts in each direction
--------------------------------------------------------------------------------

local CartProdDecomp = {}
-- constructor to make new cart-prod decomp
function CartProdDecomp:new(tbl)
   local self = setmetatable({}, CartProdDecomp)
   
   local ones = {}
   for d = 1, #tbl.cuts do ones[d] = 1 end

   self._cutsRange = Range.Range(ones, tbl.cuts)

   if not tbl.__serTesting then
      self._comm = Mpi.COMM_WORLD -- default is to use MPI_COMM_WORLD   
      if Mpi.Comm_size(self._comm) ~= self._cutsRange:volume() then
	 assert(false,
		string.format(
		   "CartProdDecomp: Number of sub-domains (%d) and processors (%d) must match",
		   self._cutsRange:volume(), Mpi.Comm_size(self._comm)))
	 
      end
   end
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(CartProdDecomp, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
CartProdDecomp.__index = {
   comm = function (self)
      return self._comm
   end,
   ndim = function (self)
      return self._cutsRange:ndim()
   end,
   cuts = function (self, dir)
      return self._cutsRange:upper(dir)
   end,
   decompose = function (self, range) -- decompose range
      local decompRgn = DecomposedRange(self)
      
      local shapes, starts = {}, {} -- to store shapes, start in each direction
      -- compute shape and start index in each direction
      for dir = 1, range:ndim() do
	 shapes[dir] = Lin.IntVec(self:cuts(dir))
	 local baseShape = math.floor(range:shape(dir)/self:cuts(dir))
	 local remCells = range:shape(dir) % self:cuts(dir) -- extra cells left over
	 for c = 1, self:cuts(dir) do
	    shapes[dir][c] = baseShape + (remCells>0 and 1 or 0) -- add extra cell if any remain
	    remCells = remCells-1
	 end

	 starts[dir] = Lin.IntVec(self:cuts(dir))
	 starts[dir][1] = range:lower(dir)
	 for c = 2, self:cuts(dir) do
	    starts[dir][c] = starts[dir][c-1] + shapes[dir][c-1]
	 end
      end

      local c = 1
      -- add regions
      for idx in self._cutsRange:colMajorIter() do
	 local l, u = {}, {}
	 for dir = 1, range:ndim() do
	    l[dir] = starts[dir][idx[dir]]
	    u[dir] = l[dir]+shapes[dir][idx[dir]]-1
	 end
	 decompRgn._domains[c] = Range.Range(l, u)
	 c = c+1
      end

      local cutIdxr = Range.makeColMajorGenIndexer(self._cutsRange)
      -- loop over each direction, adding boundary sub-regions to
      -- pair-list to allow communication in periodic directions
      for dir = 1, range:ndim() do
	 -- loop over "shortened" range which are basically the
	 -- sub-domains that lie on the lower domain boundary
	 local shortRange = self._cutsRange:shorten(dir)	 
	 local c = 1
	 for idx in shortRange:colMajorIter() do
	    decompRgn._periodicDomPairs[dir][c].lower = cutIdxr(idx) -- lower sub-domain on boundary
	    idx[dir] = idx[dir]+self:cuts(dir) 
	    decompRgn._periodicDomPairs[dir][c].upper = cutIdxr(idx) -- upper sub-domain on boundary
	    c = c+1
	 end
      end

      return decompRgn
   end,
}

return {
   CartProd = CartProdDecomp,
}
