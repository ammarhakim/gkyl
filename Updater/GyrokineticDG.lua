--
-- This updater wraps g0's dg_updater_gyrokinetic to advance the Gyrokinetic
-- equations with Discontinuous Galerkin scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto       = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local ffi         = require "ffi"
local xsys        = require "xsys"

local ffiC = ffi.C
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

ffi.cdef [[ 
// Object type
typedef struct gkyl_dg_updater_gyrokinetic gkyl_dg_updater_gyrokinetic;

// Boundary condition types.
enum gkyl_gkeqn_id {
  GKYL_GK_DEFAULT=0,    // Default GK equation (kernels).
//  GKYL_GK_SOFTMIRROR,
};

/**
 * Create new updater to update gyrokinetic equations using hyper dg.
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Config space range
 * @param vel_range Velocity space range
 * @param field_id Enum identifier for field type (see gkyl_eqn_type.h)
 *
 * @return New gyrokinetic updater object
 */
gkyl_dg_updater_gyrokinetic* gkyl_dg_updater_gyrokinetic_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range,
  enum gkyl_gkeqn_id eqn_id, double charge, double mass, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up gyrokinetic updater object
 * @param field_id Enum identifier for field type (see gkyl_eqn_type.h)
 * @param update_rng Range on which to compute.
 * @param bmag Magnetic field magnitude.
 * @param jacobtot_inv Reciprocal of the total (conf * guiding center) Jacobian.
 * @param cmag Multiplicative factor in Clebsch-like form of the magnetic field.
 * @param b_i Covariant components of hat{b}.
 * @param phi Electrostatic potential.
 * @param apar Parallel component of magnetic vector potential.
 * @param apardot Time rate of change of apar.
 * @param fIn Input to updater.
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_gyrokinetic_advance(gkyl_dg_updater_gyrokinetic *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array *bmag, const struct gkyl_array *jacobtot_inv,
  const struct gkyl_array *cmag, const struct gkyl_array *b_i,
  const struct gkyl_array *phi, const struct gkyl_array *apar,
  const struct gkyl_array *apardot, const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_updater_gyrokinetic_release(gkyl_dg_updater_gyrokinetic* up);
]]

-- Gyrokinetic DG solver updater object.
local GyrokineticDG = Proto(UpdaterBase)

function GyrokineticDG:init(tbl)
   GyrokineticDG.super.init(self, tbl)

   -- Read data from input file.
   self._onGrid = assert(tbl.onGrid, "Updater.GyrokineticDG: Must provide grid object using 'onGrid'")
   self._phaseBasis = assert(tbl.phaseBasis, "Updater.GyrokineticDG: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(tbl.confBasis, "Updater.GyrokineticDG: Must specify conf-space basis functions to use using 'confBasis'")

   self._confRange = assert(tbl.confRange, "Updater.GyrokineticDG: Must specify conf-space range using 'confRange'")
   assert(self._confRange:isSubRange()==1, "Eq.Gyrokinetic: confRange must be a sub-range") 

   local charge = assert(tbl.charge, "GyrokineticDG: must specify charge using 'charge' ")
   local mass   = assert(tbl.mass, "GyrokineticDG: must specify mass using 'mass' ")

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   local equationId = "GKYL_GK_DEFAULT"

   self._zero = ffi.gc(ffiC.gkyl_dg_updater_gyrokinetic_new(self._onGrid._zero, self._confBasis._zero,
                       self._phaseBasis._zero, self._confRange, nil, equationId, charge, mass, self._useGPU),
                       ffiC.gkyl_dg_updater_gyrokinetic_release)

   return self
end

-- Advance method.
function GyrokineticDG:_advance(tCurr, inFld, outFld)
   local qIn        = assert(inFld[1], "GyrokineticDG.advance: Must specify an input field.")
   local potentials = assert(inFld[2], "GyrokineticDG.advance: Must specify electromagnetic potentials.")
   local geo        = assert(inFld[3], "GyrokineticDG.advance: Must specify geometry fields.")

   local qRhsOut       = assert(outFld[1], "GyrokineticDG.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "GyrokineticDG.advance: Must pass cflRate field in output table")

   local localRange = qRhsOut:localRange()
   ffiC.gkyl_dg_updater_gyrokinetic_advance(self._zero, localRange, geo.bmag._zero, geo.jacobTotInv._zero,
      geo.cmag._zero, geo.b_i._zero, potentials.phi._zero, potentials.apar._zero, potentials.dApardt._zero,
      qIn._zero, cflRateByCell._zero, qRhsOut._zero)
end

function GyrokineticDG:_advanceOnDevice(tCurr, inFld, outFld)
   local qIn        = assert(inFld[1], "GyrokineticDG.advance: Must specify an input field.")
   local potentials = assert(inFld[2], "GyrokineticDG.advance: Must specify electromagnetic potentials.")
   local geo        = assert(inFld[3], "GyrokineticDG.advance: Must specify geometry fields.")

   local qRhsOut       = assert(outFld[1], "GyrokineticDG.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "GyrokineticDG.advance: Must pass cflRate field in output table")

   local localRange = qRhsOut:localRange()
   ffiC.gkyl_dg_updater_gyrokinetic_advance(self._zero, localRange, geo.bmag._zeroDevice, geo.jacobTotInv._zeroDevice,
      geo.cmag._zeroDevice, geo.b_i._zeroDevice, potentials.phi._zeroDevice, potentials.apar._zeroDevice, potentials.dApardt._zeroDevice,
      qIn._zeroDevice, cflRateByCell._zeroDevice, qRhsOut._zeroDevice)
end

-- Set up pointers to dt and cflRateByCell.
function GyrokineticDG:setDtAndCflRate(dt, cflRateByCell)
   GyrokineticDG.super.setDtAndCflRate(self, dt, cflRateByCell)

   if self._onDevice then
      ffiC.setDtAndCflRate(self._onDevice, dt, cflRateByCell._onDevice)
   end
end

-- Fetch equation updater.
function GyrokineticDG:getEquation() return self._zero end

return GyrokineticDG
