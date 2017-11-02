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
   self._commSet = decomp:commSet() -- set of communicators
   self._domains = RangeVec(self._cutsRange:volume())

   -- for periodic directions store "skeleton" sub-domains, i.e. those
   -- that touch the domain boundary
   self._periodicDomPairs = {} -- table of size 'ndim'
   for d = 1, decomp:ndim() do
      local sz = self._cutsRange:shorten(d):volume() -- number of skeleton cells
      self._periodicDomPairs[d] = PairsVec(sz) -- store as integer pair
   end

   -- (we can do this here as skeleton sub-domain IDs only depends on
   -- cuts and not on domain to be decomposed)
   local cutsIdxr = Range.makeColMajorGenIndexer(self._cutsRange)
   -- loop over each direction, adding boundary sub-regions to
   -- pair-list to allow communication in periodic directions
   for dir = 1, decomp:ndim() do
      -- loop over "shortened" range which are basically
      -- sub-domains that lie on the lower domain boundary
      local shortRange = self._cutsRange:shorten(dir)
      local c = 1
      for idx in shortRange:colMajorIter() do
	 local idxp = idx:copy()
	 idxp[dir] = idxp[dir]+decomp:cuts(dir)-1 -- upper index

	 self._periodicDomPairs[dir][c].lower = cutsIdxr(idx) -- lower sub-domain on boundary
	 self._periodicDomPairs[dir][c].upper = cutsIdxr(idxp) -- upper sub-domain on boundary
	 c = c+1
      end
   end   
   
   return self   
end
-- make object callable, and redirect call to the :new method
setmetatable(DecomposedRange, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
DecomposedRange.__index = {
   commSet = function (self)
      return self._commSet
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
   end,
   boundarySubDomainIds = function (self, dir)
      return self._periodicDomPairs[dir]
   end,
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
   self._useShared = tbl.useShared and tbl.useShared or false

   local comm, shmComm = Mpi.COMM_WORLD, nil
   -- create various communicators
   if self._useShared then
      shmComm = Mpi.Comm_split_type(comm, Mpi.COMM_TYPE_SHARED, 0, Mpi.INFO_NULL)
   else
      -- when not using MPI-SHM, make each processor its own "SHM"
      -- communicator (i.e each SHM comm has 1 proc)
      local ranks = Lin.IntVec(1)
      ranks[1] = Mpi.Comm_rank(comm) -- only local rank is "shared"
      shmComm = Mpi.Split_comm(comm, ranks)
   end

   local worldSz = Mpi.Comm_size(comm)
   local shmSz = Mpi.Comm_size(shmComm)
   
   -- get mapping from shared -> global communicator ranks
   local worldGrp = Mpi.Comm_group(comm)
   local shmGrp = Mpi.Comm_group(shmComm)
   local shmRanks = Lin.IntVec(shmSz)
   for d = 1, shmSz do shmRanks[d] = d-1 end -- ranks are indexed from 0
   local worldRanks = Mpi.Group_translate_ranks(shmGrp, shmRanks, worldGrp)

   -- store world rank corresponding to rank zero of shmComm
   local nodeSz = worldSz/shmSz
   local localZeroRanks = Lin.IntVec(nodeSz)
   for d = 1, #localZeroRanks do localZeroRanks[d] = 0 end

   -- collect global ranks corresponding to rank 0 of shmComm
   local worldRank = Mpi.Comm_rank(comm) -- global rank of process
   if worldRank % shmSz == 0 then
      localZeroRanks[worldRank/shmSz+1] = worldRanks[1]
   end
   local zeroRanks = Lin.IntVec(nodeSz)
   Mpi.Allreduce(localZeroRanks:data(), zeroRanks:data(), nodeSz, Mpi.INT, Mpi.SUM, comm)
   -- (the above Allreduce ensures that the zeroRanks vector is
   -- identical on processes in comm. It works as the localZeroRanks
   -- only has a single non-zero element that at a unique index
   -- location).

   -- now create nodeComm from the collected rank zeros
   local nodeComm = Mpi.Split_comm(comm, zeroRanks)

   -- check if total number of domains specified by 'cuts' matches
   -- number of MPI ranks in nodeComm
   if not tbl.__serTesting then -- skip check if just testing
      if Mpi.Is_comm_valid(nodeComm) then
	 if Mpi.Comm_size(nodeComm) ~= self._cutsRange:volume() then
	    assert(false,
		   string.format(
		      "CartProdDecomp: Number of sub-domains (%d) and processors (%d) must match",
		      self._cutsRange:volume(), Mpi.Comm_size(nodeComm)))
	 end
      end
   end

   -- store various communicators in a table
   self._commSet = {
      comm = comm, sharedComm = shmComm, nodeComm = nodeComm
   }
   
   return self
end
-- make object callable, and redirect call to the :new method
setmetatable(CartProdDecomp, { __call = function (self, o) return self.new(self, o) end })

-- set callable methods
CartProdDecomp.__index = {
   commSet = function (self)
      return self._commSet
   end,
   ndim = function (self)
      return self._cutsRange:ndim()
   end,
   cuts = function (self, dir)
      return self._cutsRange:upper(dir)
   end,
   isShared = function (self)
      return self._useShared
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

      return decompRgn
   end,
}

return {
   CartProd = CartProdDecomp,
}
