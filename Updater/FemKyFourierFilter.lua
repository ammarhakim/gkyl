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
  typedef struct FemKyFourierFilter FemKyFourierFilter;
  FemKyFourierFilter* new_FemKyFourierFilter(int nx, int ny, int nz, int ndim, int *ikyFilter, int numFilter);
  void delete_FemKyFourierFilter(FemKyFourierFilter* f);
  void assembleGlobalSrc(FemKyFourierFilter* f, double *data, int idx, int idy, int idz);
  void getFilteredSolution(FemKyFourierFilter* f, double *data, int idx, int idy, int idz);
  void fft_r2c(FemKyFourierFilter* f);
  void fft_c2r(FemKyFourierFilter* f);
  void filter(FemKyFourierFilter* f);
]]

-- FEM FourierFilter solver updater object.
local FemKyFourierFilter = Proto(UpdaterBase)

function FemKyFourierFilter:init(tbl)
   FemKyFourierFilter.super.init(self, tbl)

   self._grid  = assert(tbl.onGrid, "Updater.FemKyFourierFilter: Must provide grid object using 'onGrid'")
   self._basis = assert(tbl.basis, "Updater.FemKyFourierFilter: Must specify basis functions to use using 'basis'")
   self._ndim = self._grid:ndim()
   self._p  = self._basis:polyOrder()

   assert(self._basis:id()=="serendipity", "Updater.FemKyFourierFilter: only implemented for modal serendipity basis")
   assert(self._ndim == self._basis:ndim(), "Updater.FemKyFourierFilter: dimensions of basis and grid must match")
   assert(self._p==1, "Updater.FemKyFourierFilter: only implemented for polyOrder = 1")
   assert(self._ndim==2 or self._ndim==3, "Updater.FemKyFourierFilter: only implemented for 2D or 3D (with no filter in 3rd dimension)")

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
      self._iky[i] = math.floor(self._kyTbl[i+1]*Ly/(2*math.pi)+0.5)
   end

   -- use nx+1 and nz+1 to include node at end of x and z dimensions
   -- for y, we will assume periodicity, so we do not include this last node
   self._filter = ffiC.new_FemKyFourierFilter(self._nx+1, self._ny, self._nz+1, self._ndim, self._iky, self._numFilter)

   return self
end

function FemKyFourierFilter:_advance(tCurr, inFld, outFld) 
   -- Advance method. Assembles the linear problem and solves it.
   local ndim = self._ndim
   local grid = self._grid

   local src = assert(inFld[1], "FemKyFourierFilter.advance: Must specify an input field")
   local sol = assert(outFld[1], "FemKyFourierFilter.advance: Must specify an output field")
   local localRange = sol:localRange()
   local localExtRange = sol:localExtRange()

   local srcPtr = src:get(1)
   local solPtr = sol:get(1)
   local indexer = src:genIndexer()

   for idx in localExtRange:rowMajorIter() do
      -- get local x,y,z indices for use in fft data layout (these start from 0 on each proc)
      local ix = idx[1] - localExtRange:lower(1) - 1
      local iy = idx[2] - localExtRange:lower(2) - 1
      local iz = 0
      if ndim==3 then
        iz = idx[3] - localExtRange:lower(3) - 1
      end
      if (ix>=0 and iy>=0 and iy<self._ny and (ndim==2 or iz>=0)) then
        src:fill(indexer(idx), srcPtr)
        ffiC.assembleGlobalSrc(self._filter, srcPtr:data(), ix, iy, iz)
      end
   end

   ffiC.fft_r2c(self._filter)

   ffiC.filter(self._filter)

   ffiC.fft_c2r(self._filter)

   for idx in localRange:rowMajorIter() do
      sol:fill(indexer(idx), solPtr)
      -- get local x,y,z indices for use in fft data layout (these start from 0 on each proc)
      local ix = idx[1] - localExtRange:lower(1) - 1
      local iy = idx[2] - localExtRange:lower(2) - 1
      local iz = 0
      if ndim==3 then
        iz = idx[3] - localExtRange:lower(3) - 1
      end
      ffiC.getFilteredSolution(self._filter, solPtr:data(), ix, iy, iz)
   end
end

function FemKyFourierFilter:delete()
   for idz=self.local_z_lower, self.local_z_upper do
      ffiC.delete_FemKyFourierFilter(self._poisson[idz])
   end
end

return FemKyFourierFilter
