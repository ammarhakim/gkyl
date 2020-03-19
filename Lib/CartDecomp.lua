-- Gkyl ------------------------------------------------------------------------
--
-- Decomposition algorithms for rectangular grids.
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- system libraries
local ffi = require "ffi"

-- Gkyl libraries
local Lin = require "Lib.Linalg"
local Mpi = require "Comm.Mpi"
local PrimeFactor = require "Lib.PrimeFactor"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local xsys = require "xsys"

-- create constructor to store vector of Range objects
local RangeVec = Lin.new_vec_ct(ffi.typeof("Range_t"))
-- create constructor to store integer pairs
local PairsVec = Lin.new_vec_ct(ffi.typeof("struct { int32_t lower, upper; } "))

-- DecomposedRange --------------------------------------------------------------
--
-- The decomposition object's decompose() method creates a
-- DecomposedRange object
--------------------------------------------------------------------------------

local DecomposedRange = Proto()
-- constructor to make object that stores decomposed ranges
function DecomposedRange:init(decomp)
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
end

-- set callable methods
function DecomposedRange:commSet() return self._commSet end
function DecomposedRange:ndim() return self._cutsRange:ndim() end
function DecomposedRange:cuts(dir) return self._cutsRange:upper(dir) end
function DecomposedRange:numSubDomains() return self._cutsRange:volume() end
function DecomposedRange:subDomain(k) return self._domains[k] end
function DecomposedRange:boundarySubDomainIds(dir) return self._periodicDomPairs[dir] end
function DecomposedRange:cutsIndexer() return Range.makeColMajorGenIndexer(self._cutsRange) end
function DecomposedRange:cutsInvIndexer() return Range.makeColMajorInvIndexer(self._cutsRange) end

-- CartProdDecomp --------------------------------------------------------------
--
-- Decomposition of grid into smaller boxes, specified by number of
-- cuts in each direction
--------------------------------------------------------------------------------

local CartProdDecomp = Proto()
-- constructor to make new cart-prod decomp
function CartProdDecomp:init(tbl)
   local ones = {}
   for d = 1, #tbl.cuts do ones[d] = 1 end

   self._cutsRange = Range.Range(ones, tbl.cuts)
   self._useShared = xsys.pickBool(tbl.useShared, false)

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
   -- only has a single non-zero element at a unique index location).

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
end

-- set callable methods
function CartProdDecomp:commSet() return self._commSet end
function CartProdDecomp:ndim() return self._cutsRange:ndim() end
function CartProdDecomp:cuts(dir) return self._cutsRange:upper(dir) end
function CartProdDecomp:isShared() return self._useShared end

function CartProdDecomp:decompose(range) -- decompose range
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
end

--------------------------------------------------------------------------------
--
-- Set utility functions to construct a reasonable decomposition when
-- user has not specified one explicitly
--
--------------------------------------------------------------------------------

-- check if numSubDoms^(1/ndim) is an integer
local function isEvenDecomp(ndim, numSubDoms)

   -- This code is tricky as round off errors can mess the simple (and
   -- obvious) test for even decomposition, i.e. math.ceil(v) =
   -- math.floor(v) where v = numSubDoms^(1/ndim)
   
   local c = math.pow(numSubDoms, 1/ndim)
   local cdown = math.floor(c)
   local cup = math.ceil(c)

   local cuts = {}
   if cdown^ndim == numSubDoms then
      for i = 1, ndim do cuts[i] = cdown end
      return true, cuts
   end
   if cup^ndim == numSubDoms then
      for i = 1, ndim do cuts[i] = cup end
      return true, cuts      
   end
   return false, cuts
end

local function calcCuts(ndim, numSubDoms)

   -- first check if we can create an "even" decomposition
   local evFlag, evCuts = isEvenDecomp(ndim, numSubDoms)
   if evFlag then return evCuts end
   
   -- Compute prime factorization of numSubDoms and then use this to
   -- compute "cuts". This algorithm may not produce the best possible
   -- decomposition, but it is hard to tell what "best" means. For
   -- most reasonable cases (say numSubDoms is not a prime or small
   -- multiple of a prime) the decompositions are probably okay.
   
   local factors = PrimeFactor.all(numSubDoms)
   local div, rem = math.floor(#factors/ndim), #factors%ndim
   local cuts, idx = {}, 1
   for d = 1, ndim do
      cuts[d] = 1
      for i = 1, div do
	 cuts[d] = cuts[d]*factors[idx]
	 idx = idx+1
      end
      if rem > 0 then
	 cuts[d] = cuts[d]*factors[idx]
	 idx = idx+1
	 rem = rem-1
      end
   end
   return cuts
end

-- Construct a reasonable decomposition (cuts) given dimensions,
-- number of sub-domains and cells along each dimension. The returned
-- cuts table can be used in CartProdDecomp to decompose a given range
-- object.
local function makeCuts(ndim, numSubDoms, cells)
   -- This initial sort is required as calcCuts does not return cuts
   -- in any particular order. The sort ensures that the rearrangement
   -- of cuts (done in the second for-loop) then matches order of
   -- cells table
   local rawCuts = calcCuts(ndim, numSubDoms); table.sort(rawCuts)

   -- The idea here is that we need cuts arranged in same order as
   -- cell sizes. This is needed as "rawCuts" is not returned in any
   -- particular order, and needs to be rearranged to give a
   -- reasonable decomposition
   
   local cellsWithIdx = {}
   for i = 1, #cells do
      cellsWithIdx[i] = { cells[i], i }
   end
   table.sort(cellsWithIdx, function(l,r) return l[1]<r[1] end)

   local cuts = {}
   for i = 1, #cells do
      cuts[i] = rawCuts[cellsWithIdx[i][2]]
   end
   
   return cuts
end

return {
   CartProd = CartProdDecomp,
   makeCuts = makeCuts,
}
