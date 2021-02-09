-- Gkyl --------------------------------------------------------------------------
--
-- Updater to convert conservative variables to primitive variables and vice versa
-- For 10 moments system converts conservative variables to primitive variables:
-- rho, ux, uy, uz, Pxx, Pxy, Pxz, Pyy, Pyz, Pzz
-- 
-- And for the transformation back, the result is:
-- rho, rho*ux, rho*uy, rho*uz,  Pxx + rho*ux*ux, Pxy + rho*ux*uy, 
-- Pxz + rho*ux*uz, Pyy + rho*uy*uy, Pyz + rho*uy*uz, Pzz + rho*uz*uz
--    _______     ___
-- + 6 @ |||| # P ||| +
----------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

-- system libraries
local ffi = require "ffi"

-- Define C types for storing private data for use in updater
ffi.cdef [[
  void gkylTenMomentConservativeToPrimitive(double *f);
  void gkylTenMomentPrimitiveToConservative(double *f);
]]

-- Explicit, SSP RK3 scheme
local function updateConservativeToPrimitive(fPtr)
   ffi.C.gkylTenMomentConservativeToPrimitive(fPtr)
end

-- Use an implicit scheme to update momentum and electric field
local function updatePrimitiveToConservative(fPtr)
   ffi.C.gkylTenMomentPrimitiveToConservative(fPtr)
end

-- Ten-moment source updater object
local TenMomentConservToPrim = Proto(UpdaterBase)

function TenMomentConservToPrim:init(tbl)
   TenMomentConservToPrim.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid, "Updater.TenMomentConservToPrim: Must provide grid object using 'onGrid'")
   self._ndim = self._onGrid:ndim()
   local mode = tbl.mode and tbl.mode or "ToPrimitive"
   self._update = nil
   if mode == "ToPrimitive" then
      self._update = updateConservativeToPrimitive
   elseif mode == "ToConservative" then
      self._update = updatePrimitiveToConservative
   end
end

-- advance method
function TenMomentConservToPrim:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local ndim = self._ndim

   local f = assert(outFld[1], "TenMomentGrad.advance: Must specify an output field")
   local fIdxr = f:genIndexer() -- indexer functions into fields
   local idxc = Lin.IntVec(ndim) -- index on left/center/right

   -- pointers for (re)use in update
   local fC = f:get(1)
   local localRange = f:localRange()  
   for idx in localRange:rowMajorIter() do
      idx:copyInto(idxc)
      -- get grid information and set pointers for current cell
      grid:setIndex(idxc)
      f:fill(fIdxr(idxc), fC)
      self._update(fC:data())
   end

   return true, GKYL_MAX_DOUBLE
end

return TenMomentConservToPrim