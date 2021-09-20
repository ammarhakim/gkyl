-- Gkyl ------------------------------------------------------------------------
--
-- Updater to solve FourierFilter (or rather Helmholtz) equation
--
--   div{ epsilon nabla_perp{ phi } } + beta * phi = sigma
--
-- in perpendicular directions with FEM scheme. Perpendicular directions assumed
-- to be first two configuration-space directions
--
-- We use the following terminology:
--   epsilon: Laplacian weight.
--   beta: modifier weight.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc       = require "Lib.Alloc"
local UpdaterBase = require "Updater.Base"
local Lin         = require "Lib.Linalg"
local Proto       = require "Lib.Proto"
local Range       = require "Lib.Range"
local ffi         = require "ffi"
local ffiC        = ffi.C
local xsys        = require "xsys"
local ProjectOnBasis = require "Updater.ProjectOnBasis"
local DataStruct     = require "DataStruct"
local Time           = require "Lib.Time"
local Logger         = require "Lib.Logger"
local Mpi
if GKYL_HAVE_MPI then Mpi = require "Comm.Mpi" end

ffi.cdef[[
  typedef struct FemPerpFourierFilter FemPerpFourierFilter;
  FemPerpFourierFilter* new_FemPerpFourierFilter(int nx, int ny, int ndim, int *ikyFilter, int numFilter);
  void delete_FemPerpFourierFilter(FemPerpFourierFilter* f);
  void assembleGlobalSrc(FemPerpFourierFilter* f, double *data, int idx, int idy);
  void getFilteredSolution(FemPerpFourierFilter* f, double *data, int idx, int idy);
  void fft_r2c(FemPerpFourierFilter* f);
  void fft_c2r(FemPerpFourierFilter* f);
  void filter(FemPerpFourierFilter* f);
]]

-- FEM FourierFilter solver updater object.
local FemPerpFourierFilter = Proto(UpdaterBase)

function FemPerpFourierFilter:init(tbl)
   FemPerpFourierFilter.super.init(self, tbl)

   self._grid  = assert(tbl.onGrid, "Updater.FemPerpFourierFilter: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.FemPerpFourierFilter: Must specify basis functions to use using 'basis'")
   self._ikyTbl = tbl.kyfilter

   self._ndim = self._grid:ndim()
   self._nx = self._grid:numCells(1)
   self._ny = self._grid:numCells(2)
   self._p  = self._basis:polyOrder()

   assert(self._basis:id()=="serendipity", "Updater.FemPerpFourierFilter: only implemented for modal serendipity basis")
   assert(self._ndim == self._basis:ndim(), "Updater.FemPerpFourierFilter: dimensions of basis and grid must match")
   assert(self._p==1, "Updater.FemPerpFourierFilter: only implemented for polyOrder = 1")
   assert(self._ndim==2 or self._ndim==3, "Updater.FemPerpFourierFilter: only implemented for 2D or 3D (with no filter in 3rd dimension)")

   self._onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._first = true

   --if GKYL_HAVE_MPI then
   --   -- Split communicators in z.
   --   local commSet   = self._grid:commSet()
   --   local worldComm = commSet.comm
   --   local nodeComm  = commSet.nodeComm
   --   local nodeRank  = Mpi.Comm_rank(nodeComm)
   --   local zrank     = 0
   --   if self._ndim==3 then zrank = math.floor(nodeRank/self._grid:cuts(1)/self._grid:cuts(2)) end
   --   self._zcomm = Mpi.Comm_split(worldComm, zrank, nodeRank)
   --end

   --local localRange = self._grid:localRange()
   ---- Create region that is effectively 2d and global in x-y directions.
   --self.local_z_lower = 1
   --self.local_z_upper = 1
   --if (self._ndim == 3) then
   --   self.local_z_lower = localRange:lower(3)
   --   self.local_z_upper = localRange:upper(3)
   --end

   self._numFilter = #self._ikyTbl
   self._iky = ffi.new("int[?]", self._numFilter)
   for i = 0, self._numFilter-1 do
      self._iky[i] = self._ikyTbl[i+1]
   end

   self._filter = ffiC.new_FemPerpFourierFilter(self._nx, self._ny, self._ndim, self._iky, self._numFilter)

   return self
end

function FemPerpFourierFilter:_advance(tCurr, inFld, outFld) 
   -- Advance method. Assembles the linear problem and solves it.
   local ndim = self._ndim
   local grid = self._grid

   local src = assert(inFld[1], "FemPerpFourierFilter.advance: Must specify an input field")
   local sol = assert(outFld[1], "FemPerpFourierFilter.advance: Must specify an output field")
   local localRange = sol:localRange()
   local localExtRange = sol:localExtRange()

   local srcPtr = src:get(1)
   local solPtr = sol:get(1)
   local indexer = src:genIndexer()

   for idx in localRange:rowMajorIter() do
      src:fill(indexer(idx), srcPtr)
      ffiC.assembleGlobalSrc(self._filter, srcPtr:data(), idx[1]-1, idx[2]-1)
   end

   ffiC.fft_r2c(self._filter)

   ffiC.filter(self._filter)

   ffiC.fft_c2r(self._filter)

   for idx in localRange:rowMajorIter() do
      sol:fill(indexer(idx), solPtr)
      ffiC.getFilteredSolution(self._filter, solPtr:data(), idx[1]-1, idx[2]-1)
   end
end

function FemPerpFourierFilter:delete()
   for idz=self.local_z_lower, self.local_z_upper do
      ffiC.delete_FemPerpFourierFilter(self._poisson[idz])
   end
end

return FemPerpFourierFilter
