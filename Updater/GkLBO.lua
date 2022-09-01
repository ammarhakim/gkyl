--------------------------------------------------------------------------------
--
-- This updater wraps g0's dg_updater_lbo_gyrokinetic to advance the LBO collisions for Gk
-- equations with Discontinuous Galerkin scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc        = require "Lib.Alloc"
local Basis        = require "Basis.BasisCdef"
local DataStruct   = require "DataStruct"
local EqBase       = require "Eq.EqBase"
local Grid         = require "Grid.RectCart"
local CartField    = require "DataStruct.CartField"
local Lin          = require "Lib.Linalg"
local Mpi          = require "Comm.Mpi"
local Proto        = require "Lib.Proto"
local Range        = require "Lib.Range"
local UpdaterBase  = require "Updater.Base"
local ffi          = require "ffi"
local xsys         = require "xsys"

local ffiC = ffi.C
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

ffi.cdef [[ 

typedef struct gkyl_dg_updater_lbo_gyrokinetic gkyl_dg_updater_lbo_gyrokinetic;

gkyl_dg_updater_lbo_gyrokinetic*
gkyl_dg_updater_lbo_gyrokinetic_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_range, double mass, bool use_gpu);

void
gkyl_dg_updater_lbo_gyrokinetic_advance(gkyl_dg_updater_lbo_gyrokinetic *lbo,
  const struct gkyl_range *update_rng, const struct gkyl_array *bmagInv,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq,
  const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

void
gkyl_dg_updater_lbo_gyrokinetic_advance_cu(gkyl_dg_updater_lbo_gyrokinetic *lbo,
  const struct gkyl_range *update_rng, const struct gkyl_array *bmagInv,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq,
  const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

void gkyl_dg_updater_lbo_gyrokinetic_release(gkyl_dg_updater_lbo_gyrokinetic *lbo);
]]

-- GkLBO DG solver updater object
local GkLBO = Proto(UpdaterBase)

function GkLBO:init(tbl)
   GkLBO.super.init(self, tbl)

   -- Read data from input file.
   self._onGrid     = assert(tbl.onGrid, "Updater.GkLBO: Must provide grid object using 'onGrid'")
   self._phaseBasis = assert(tbl.phaseBasis, "Updater.GkLBO: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis  = assert(tbl.confBasis, "Updater.GkLBO: Must specify conf-space basis functions to use using 'confBasis'")
   self._confRange  = assert(tbl.confRange, "Updater.GkLBO: Must specify conf-space range using 'confRange'")
   self._mass       = assert(tbl.mass, "Updater.GkLBO: Must specify species mass using 'mass'")

   self._zero = ffi.gc(ffiC.gkyl_dg_updater_lbo_gyrokinetic_new(self._onGrid._zero, self._confBasis._zero, self._phaseBasis._zero,
                                                                self._confRange, self._mass, GKYL_USE_GPU or 0),
                       ffiC.gkyl_dg_updater_lbo_gyrokinetic_release)

   return self
end

-- advance method
function GkLBO:_advance(tCurr, inFld, outFld)

   local fIn     = assert(inFld[1], "GkLBO.advance: Must pass input distf")
   local bmagInv = assert(inFld[2], "GkLBO.advance: Must pass input bmagInv")
   local nu_u    = assert(inFld[3], "GkLBO.advance: Must pass nu_u")
   local nu_vtsq = assert(inFld[4], "GkLBO.advance: Must pass nu_vtsq")
   local nu_sum  = assert(inFld[5], "GkLBO.advance: Must pass nu_sum")
 
   local fRhsOut       = assert(outFld[1], "GkLBO.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "GkLBO.advance: Must pass cflRate field in output table")

   local localRange = fRhsOut:localRange()
   ffiC.gkyl_dg_updater_lbo_gyrokinetic_advance(self._zero, localRange, bmagInv._zero, nu_sum._zero, nu_u._zero, nu_vtsq._zero, fIn._zero, cflRateByCell._zero, fRhsOut._zero)

end

-- advance method
function GkLBO:_advanceOnDevice(tCurr, inFld, outFld)

   local fIn     = assert(inFld[1], "GkLBO.advance: Must pass input distf")
   local bmagInv = assert(inFld[2], "GkLBO.advance: Must pass input bmagInv")
   local nu_u    = assert(inFld[3], "GkLBO.advance: Must pass nu_u")
   local nu_vtsq = assert(inFld[4], "GkLBO.advance: Must pass nu_vtsq")
   local nu_sum  = assert(inFld[5], "GkLBO.advance: Must pass nu_sum")
 
   local fRhsOut       = assert(outFld[1], "GkLBO.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "GkLBO.advance: Must pass cflRate field in output table")

   local localRange = fRhsOut:localRange()
   ffiC.gkyl_dg_updater_lbo_gyrokinetic_advance_cu(self._zero, localRange, bmagInv._zeroDevice, nu_sum._zeroDevice, nu_u._zeroDevice, nu_vtsq._zeroDevice, fIn._zeroDevice, cflRateByCell._zeroDevice, fRhsOut._zeroDevice)

end

return GkLBO
