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
// Object type
typedef struct gkyl_dg_updater_collisions gkyl_dg_updater_collisions;

/**
 * Create new updater to update lbo equations using hyper dg.
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Config space range
 * @param mass Species mass
 * @return New LBO updater object
 */
gkyl_dg_updater_lbo_gyrokinetic* gkyl_dg_updater_lbo_gyrokinetic_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, double mass, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param lbo LBO updater object
 * @param update_rng Range on which to compute.
 * @param bmag Magnitude of magnetic field
 * @param nu_sum Sum of coll freq
 * @param nu_prim_moms Sum of coll freq*u and freq*vtsq
 * @param m2self 2nd velocity moment of this species.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_lbo_gyrokinetic_advance(gkyl_dg_updater_lbo_gyrokinetic *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *bmag_inv,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_prim_moms,
  const struct gkyl_array *m2self,
  const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

void gkyl_dg_updater_lbo_gyrokinetic_advance_cu(gkyl_dg_updater_lbo_gyrokinetic *lbo,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *bmag_inv,
  const struct gkyl_array *nu_sum, const struct gkyl_array *nu_prim_moms,
  const struct gkyl_array *m2self,
  const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

/**
 * Delete updater.
 *
 * @param lbo Updater to delete.
 */
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
   local nu_prim_moms = assert(inFld[3], "GkLBO.advance: Must pass nu_prim_moms")
   local nu_sum  = assert(inFld[4], "GkLBO.advance: Must pass nu_sum")
   local m2self  = assert(inFld[5], "GkLBO.advance: Must pass m2self") 
 
   local fRhsOut       = assert(outFld[1], "GkLBO.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "GkLBO.advance: Must pass cflRate field in output table")

   local localRange = fRhsOut:localRange()
   ffiC.gkyl_dg_updater_lbo_gyrokinetic_advance(self._zero, localRange, bmagInv._zero, nu_sum._zero, nu_prim_moms._zero, m2self._zero, fIn._zero, cflRateByCell._zero, fRhsOut._zero)

end

-- advance method
function GkLBO:_advanceOnDevice(tCurr, inFld, outFld)

   local fIn     = assert(inFld[1], "GkLBO.advance: Must pass input distf")
   local bmagInv = assert(inFld[2], "GkLBO.advance: Must pass input bmagInv")
   local nu_prim_moms = assert(inFld[3], "GkLBO.advance: Must pass nu_prim_moms")
   local nu_sum  = assert(inFld[4], "GkLBO.advance: Must pass nu_sum")
   local m2self  = assert(inFld[5], "GkLBO.advance: Must pass m2self") 
 
   local fRhsOut       = assert(outFld[1], "GkLBO.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "GkLBO.advance: Must pass cflRate field in output table")

   local localRange = fRhsOut:localRange()
   ffiC.gkyl_dg_updater_lbo_gyrokinetic_advance_cu(self._zero, localRange, bmagInv._zeroDevice, nu_sum._zeroDevice, nu_prim_moms._zeroDevice, m2self._zeroDevice, fIn._zeroDevice, cflRateByCell._zeroDevice, fRhsOut._zeroDevice)

end

return GkLBO
