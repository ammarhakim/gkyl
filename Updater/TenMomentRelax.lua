-- Gkyl ------------------------------------------------------------------------
--
-- Updater to update ten-moment relaxation. This updater allows
-- both explicit and implicit updates.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"

local COL_PIV_HOUSEHOLDER_QR = 0;
local PARTIAL_PIV_LU = 1;

-- system libraries
local ffi = require "ffi"

-- Define C types for storing private data for use in updater
ffi.cdef [[
typedef struct {
  double charge, mass; /* Charge and mass */
} FluidData_t;

typedef struct {
  double k;
  bool hasEm; /* Flag to indicate if there is: EB field */
  bool hasStatic; /* Flag to indicate if there is: static EB field */
} TenMomentRelaxData_t;

  void gkylTenMomentRelaxRk3(TenMomentRelaxData_t *sd, FluidData_t *fd, double dt, double *f, double *em);
  void gkylTenMomentRelaxExplicit(TenMomentRelaxData_t *sd, FluidData_t *fd, double dt, double *f, double *em, double *staticEm);
]]

-- Explicit, SSP RK3 scheme
local function updateRelaxRk3(self, dt, fPtr, emPtr)
   ffi.C.gkylTenMomentRelaxRk3(self._sd, self._fd, dt, fPtr, emPtr)
end

-- Use an implicit scheme to update momentum and electric field
local function updateRelaxExplicit(self, dt, fPtr, emPtr, staticEmPtr)
   ffi.C.gkylTenMomentRelaxExplicit(self._sd, self._fd, dt, fPtr, emPtr, staticEmPtr)
end

-- Ten-moment source updater object
local TenMomentRelax = Proto(UpdaterBase)

function TenMomentRelax:init(tbl)
   TenMomentRelax.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid, "Updater.TenMomentRelax: Must provide grid object using 'onGrid'")

   self._sd = ffi.new(ffi.typeof("TenMomentRelaxData_t"))   

   self._sd.k = tbl.k ~= nil and tbl.k or 0.
   self._sd.hasEm = tbl.hasEmField ~= nil and tbl.hasEmField or false
   self._sd.hasStatic = tbl.hasStaticField ~= nil and tbl.hasStaticField or false
   
   self._fd = ffi.new("FluidData_t[?]", 1)
    self._fd[0].charge = tbl.charge
    self._fd[0].mass = tbl.mass

   local scheme = tbl.scheme and tbl.scheme or "ssp-rk3"
   self._updateRelax = nil
   if scheme == "ssp-rk3" then
      self._updateRelax = updateRelaxRk3
   elseif scheme == "explicit" then
      self._updateRelax = updateRelaxExplicit
   end
end

-- advance method
function TenMomentRelax:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid
 
   local fFld = outFld[1]
   local emFld
   local staticEmFld
   if (self._sd.hasEm) then
      emFld= inFld[1] -- EM field
   end
   if (self._sd.hasStatic) then
      staticEmFld = inFld[2] -- static EM field
   end

   local fIdxr = outFld[1]:genIndexer()
   local emIdxr
   local staticEmIdxr
   if (self._sd.hasEm) then
      emIdxr = emFld:genIndexer()
   end
   if (self._sd.hasStatic) then
      staticEmIdxr = staticEmFld:genIndexer()
   end

   local fDp = ffi.new("double*")
   local emDp = ffi.new("double*")
   local staticEmDp = ffi.new("double*")

   local localRange = fFld:localRange()   
   for idx in localRange:colMajorIter() do
	    fDp = fFld:getDataPtrAt(fIdxr(idx))
      if (self._sd.hasEm) then
         emDp = emFld:getDataPtrAt(emIdxr(idx))
      end
      if (self._sd.hasStatic) then
         staticEmDp = emFld:getDataPtrAt(staticEmIdxr(idx))
      end

      self._updateRelax(self, dt, fDp, emDp, staticEmDp)
   end

   return true, GKYL_MAX_DOUBLE
end

return TenMomentRelax
