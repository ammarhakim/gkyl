--
-- This updater wraps g0's dg_updater_lbo_vlasov to advance the LBO collisions for Vlasov
-- equations with Discontinuous Galerkin scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local Basis = require "Basis.BasisCdef"
local DataStruct = require "DataStruct"
local EqBase = require "Eq.EqBase"
local Grid = require "Grid.RectCart"
local CartField = require "DataStruct.CartField"
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

ffi.cdef [[ 

typedef struct gkyl_dg_updater_lbo_vlasov gkyl_dg_updater_lbo_vlasov;

gkyl_dg_updater_lbo_vlasov*
gkyl_dg_updater_lbo_vlasov_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, const struct gkyl_range *conf_range, bool use_gpu);

void
gkyl_dg_updater_lbo_vlasov_advance(gkyl_dg_updater_lbo_vlasov *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq,
  const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

void
gkyl_dg_updater_lbo_vlasov_advance_cu(gkyl_dg_updater_lbo_vlasov *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_u, const struct gkyl_array *nu_vthsq,
  const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

void gkyl_dg_updater_lbo_vlasov_release(gkyl_dg_updater_lbo_vlasov *lbo);
]]

-- VlasovLBO DG solver updater object
local VlasovLBO = Proto(UpdaterBase)

function VlasovLBO:init(tbl)
   VlasovLBO.super.init(self, tbl)

   -- Read data from input file.
   self._onGrid = assert(tbl.onGrid, "Updater.VlasovLBO: Must provide grid object using 'onGrid'")
   self._phaseBasis = assert(tbl.phaseBasis, "Updater.VlasovLBO: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(tbl.confBasis, "Updater.VlasovLBO: Must specify conf-space basis functions to use using 'confBasis'")

   self._confRange = assert(tbl.confRange, "Updater.VlasovLBO: Must specify conf-space range using 'confRange'")

   self._zero = ffi.gc(ffiC.gkyl_dg_updater_lbo_vlasov_new(self._onGrid._zero, self._confBasis._zero, self._phaseBasis._zero, self._confRange, GKYL_USE_GPU or 0),
                       ffiC.gkyl_dg_updater_lbo_vlasov_release)

   return self
end

-- advance method
function VlasovLBO:_advance(tCurr, inFld, outFld)

   local fIn = assert(inFld[1], "VlasovLBO.advance: Must pass input distf")
   local nu_u = assert(inFld[2], "VlasovLBO.advance: Must pass nu_u")
   local nu_vthsq = assert(inFld[3], "VlasovLBO.advance: Must pass nu_vthsq")
   local nu_sum = assert(inFld[4], "VlasovLBO.advance: Must pass nu_sum")
 
   local fRhsOut = assert(outFld[1], "VlasovLBO.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "VlasovLBO.advance: Must pass cflRate field in output table")

   local localRange = fRhsOut:localRange()
   ffiC.gkyl_dg_updater_lbo_vlasov_advance(self._zero, localRange, nu_sum._zero, nu_u._zero, nu_vthsq._zero, fIn._zero, cflRateByCell._zero, fRhsOut._zero)

end

-- advance method
function VlasovLBO:_advanceOnDevice(tCurr, inFld, outFld)

   local fIn = assert(inFld[1], "VlasovLBO.advance: Must pass input distf")
   local nu_u = assert(inFld[2], "VlasovLBO.advance: Must pass nu_u")
   local nu_vthsq = assert(inFld[3], "VlasovLBO.advance: Must pass nu_vthsq")
   local nu_sum = assert(inFld[4], "VlasovLBO.advance: Must pass nu_sum")
 
   local fRhsOut = assert(outFld[1], "VlasovLBO.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "VlasovLBO.advance: Must pass cflRate field in output table")

   local localRange = fRhsOut:localRange()
   ffiC.gkyl_dg_updater_lbo_vlasov_advance_cu(self._zero, localRange, nu_sum._zeroDevice, nu_u._zeroDevice, nu_vthsq._zeroDevice, fIn._zeroDevice, cflRateByCell._zeroDevice, fRhsOut._zeroDevice)

end

return VlasovLBO