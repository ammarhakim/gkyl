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
  FemPerpFourierFilter* new_FemPerpFourierFilter(int nx, int ny, int nz, int ndim, int *ikyFilter, int numFilter);
  void delete_FemPerpFourierFilter(FemPerpFourierFilter* f);
  void assembleGlobalSrc(FemPerpFourierFilter* f, double *data, int idx, int idy, int idz);
  void getFilteredSolution(FemPerpFourierFilter* f, double *data, int idx, int idy, int idz);
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
   self._ndim = self._grid:ndim()
   self._p  = self._basis:polyOrder()

   assert(self._basis:id()=="serendipity", "Updater.FemPerpFourierFilter: only implemented for modal serendipity basis")
   assert(self._ndim == self._basis:ndim(), "Updater.FemPerpFourierFilter: dimensions of basis and grid must match")
   assert(self._p==1, "Updater.FemPerpFourierFilter: only implemented for polyOrder = 1")
   assert(self._ndim==2 or self._ndim==3, "Updater.FemPerpFourierFilter: only implemented for 2D or 3D (with no filter in 3rd dimension)")

   self._nx = self._grid:localNumCells(1)
   self._ny = self._grid:numCells(2)
   if self._ndim > 2 then self._nz = self._grid:localNumCells(3) else self._nz = 0 end

   self._onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._first = true

   self._kyTbl = tbl.kyfilter
   self._numFilter = #self._kyTbl
   self._iky = ffi.new("int[?]", self._numFilter)
   local Ly = self._grid:upper(2) - self._grid:lower(2)
   for i = 0, self._numFilter-1 do
      self._iky[i] = self._kyTbl[i+1]*Ly/(2*math.pi)
   end

   -- use nx+1 and nz+1 to include node at end of x and z dimensions
   -- for y, we will assume periodicity, so we do not include this last node
   self._filter = ffiC.new_FemPerpFourierFilter(self._nx+1, self._ny, self._nz+1, self._ndim, self._iky, self._numFilter)

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

   for idx in localExtRange:rowMajorIter() do
      if (idx[1]>0 and idx[2]>0 and idx[2]<=self._ny and (ndim==2 or idx[3]>0)) then
        src:fill(indexer(idx), srcPtr)
        ffiC.assembleGlobalSrc(self._filter, srcPtr:data(), idx[1]-1, idx[2]-1, idx[3]-1)
      end
   end

   ffiC.fft_r2c(self._filter)

   ffiC.filter(self._filter)

   ffiC.fft_c2r(self._filter)

   for idx in localRange:rowMajorIter() do
      sol:fill(indexer(idx), solPtr)
      ffiC.getFilteredSolution(self._filter, solPtr:data(), idx[1]-1, idx[2]-1, idx[3]-1)
   end
end

function FemPerpFourierFilter:delete()
   for idz=self.local_z_lower, self.local_z_upper do
      ffiC.delete_FemPerpFourierFilter(self._poisson[idz])
   end
end

return FemPerpFourierFilter
