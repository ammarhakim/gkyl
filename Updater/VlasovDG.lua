--
-- This updater wraps g0's dg_updater_vlasov to advance the Vlasov
-- equations with Discontinuous Galerkin scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc       = require "Lib.Alloc"
local Basis       = require "Basis.BasisCdef"
local DataStruct  = require "DataStruct"
local EqBase      = require "Eq.EqBase"
local Grid        = require "Grid.RectCart"
local CartField   = require "DataStruct.CartField"
local Lin         = require "Lib.Linalg"
local Mpi         = require "Comm.Mpi"
local Proto       = require "Lib.Proto"
local Range       = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local ffi         = require "ffi"
local xsys        = require "xsys"

local ffiC = ffi.C
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

ffi.cdef [[ 
// Identifiers for specific field object types
enum gkyl_field_id {
  GKYL_FIELD_E_B = 0, // Maxwell (E, B). This is default
  GKYL_FIELD_PHI = 1, // Poisson (only phi)
  GKYL_FIELD_PHI_A = 2, // Poisson with static B = curl(A) (phi, A)
  GKYL_FIELD_NULL = 3, // no field is present
};

// Identifiers for subsidary models
// These are used to distinguish things like special relativistic from non-relativistic
// or the parallel-kinetic-perpendicular-moment model
enum gkyl_model_id {
  GKYL_MODEL_DEFAULT = 0, // No subsidiary model specified.
  GKYL_MODEL_SR = 1,
  GKYL_MODEL_GEN_GEO = 2,
  GKYL_MODEL_PKPM = 3,
  GKYL_MODEL_SR_PKPM = 4,
};

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_vlasov_auxfields { 
  const struct gkyl_array *field; // q/m*(E,B) for Maxwell's, q/m*phi for Poisson's (gradient calculated in kernel)
  const struct gkyl_array *ext_field; // constant q/m*A for Poisson's (curl calculated in kernel)
  const struct gkyl_array *cot_vec; // cotangent vectors (e^i) used in volume term if general geometry enabled
  const struct gkyl_array *alpha_geo; // alpha^i (e^i . alpha) used in surface term if general geometry enabled
};

// Object type
typedef struct gkyl_dg_updater_vlasov gkyl_dg_updater_vlasov;

/**
 * Create new updater to update Vlasov equations using hyper dg.
 * Supports Vlasov-Maxwell, Vlasov-Poisson (with and without vector potential A)
 * and special relativistic Vlasov-Maxwell
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Config space range
 * @param vel_range Velocity space range
 * @param phase_range Phase space range
 * @param model_id Enum identifier for model type (e.g., SR, General Geometry, see gkyl_eqn_type.h)
 * @param field_id Enum identifier for field type (e.g., Maxwell's, Poisson, see gkyl_eqn_type.h)
 * @param aux_inp Void pointer to auxiliary fields. Void to be flexible to different auxfields structs
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 * 
 * @return New vlasov updater object
 */
gkyl_dg_updater_vlasov* gkyl_dg_updater_vlasov_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range, const struct gkyl_range *phase_range,
  enum gkyl_model_id model_id, enum gkyl_field_id field_id, void *aux_inp, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param vlasov vlasov updater object
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_vlasov_advance(gkyl_dg_updater_vlasov *vlasov,
  const struct gkyl_range *update_rng, const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

/**
 * Delete updater.
 *
 * @param vlasov Updater to delete.
 */
void gkyl_dg_updater_vlasov_release(gkyl_dg_updater_vlasov* vlasov);

]]

-- Vlasov DG solver updater object
local VlasovDG = Proto(UpdaterBase)

function VlasovDG:init(tbl)
   VlasovDG.super.init(self, tbl)

   -- Read data from input file.
   self._onGrid = assert(tbl.onGrid, "Updater.VlasovDG: Must provide grid object using 'onGrid'")
   self._phaseBasis = assert(tbl.phaseBasis,
     "Updater.VlasovDG: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(tbl.confBasis,
     "Updater.VlasovDG: Must specify conf-space basis functions to use using 'confBasis'")

   self._confRange = assert(tbl.confRange,
     "Updater.VlasovDG: Must specify conf-space range using 'confRange'")
   assert(self._confRange:isSubRange()==1, "Eq.Vlasov: confRange must be a sub-range") 

   self._phaseRange = assert(tbl.phaseRange,
     "Updater.VlasovDG: Must specify phase-space range using 'phaseRange'")
   assert(self._phaseRange:isSubRange()==1, "Eq.Vlasov: phaseRange must be a sub-range") 

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)
   
   -- Check if we have an electric and magnetic field.
   local hasElcField    = xsys.pickBool(tbl.hasElectricField, false)
   local hasMagField    = xsys.pickBool(tbl.hasMagneticField, false)
   local hasExtForce    = xsys.pickBool(tbl.hasExtForce, false)
   self._plasmaMagField = xsys.pickBool(tbl.plasmaMagField, false)

   self._modelId = "GKYL_MODEL_DEFAULT"
   
   if hasElcField and self._plasmaMagField then 
      self._fieldId = "GKYL_FIELD_E_B"
   elseif hasElcField then
      self._fieldId = "GKYL_FIELD_PHI"
   else
      self._fieldId = "GKYL_FIELD_NULL"
   end

   -- Create aux fields struct for Vlasov and set pointers based on input fields
   -- Only supports Vlasov-Maxwell and Vlasov-Poisson for now
   self._auxfieldsC = ffi.new("struct gkyl_dg_vlasov_auxfields")
   if self._useGPU then
      self._auxfieldsC.field = tbl.fldPtrs[1]._zeroDevice
   else
      self._auxfieldsC.field = tbl.fldPtrs[1]._zero
   end

   self._zero = ffi.gc(
      ffiC.gkyl_dg_updater_vlasov_new(self._onGrid._zero, 
        self._confBasis._zero, self._phaseBasis._zero, 
        self._confRange, nil, self._phaseRange,
        self._modelId, self._fieldId, self._auxfieldsC, self._useGPU),
      ffiC.gkyl_dg_updater_vlasov_release
   )

   return self
end

-- advance method
function VlasovDG:_advance(tCurr, inFld, outFld)

   local qIn = assert(inFld[1], "VlasovDG.advance: Must specify an input field")
   
   local qRhsOut = assert(outFld[1], "VlasovDG.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2],
     "VlasovDG.advance: Must pass cflRate field in output table")

   local localRange = qRhsOut:localRange()
   ffiC.gkyl_dg_updater_vlasov_advance(self._zero, localRange,
     qIn._zero, cflRateByCell._zero, qRhsOut._zero)

end

function VlasovDG:_advanceOnDevice(tCurr, inFld, outFld)

   local qIn = assert(inFld[1], "VlasovDG.advance: Must specify an input field")
   
   local qRhsOut = assert(outFld[1], "VlasovDG.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2],
     "VlasovDG.advance: Must pass cflRate field in output table")

   local localRange = qRhsOut:localRange()
   ffiC.gkyl_dg_updater_vlasov_advance(self._zero, localRange,
     qIn._zeroDevice, cflRateByCell._zeroDevice, qRhsOut._zeroDevice)

end

-- Fetch equation updater.
function VlasovDG:getEquation() return self._zero end

-- set up pointers to dt and cflRateByCell
function VlasovDG:setDtAndCflRate(dt, cflRateByCell)
   VlasovDG.super.setDtAndCflRate(self, dt, cflRateByCell)

   if self._onDevice then
      ffiC.setDtAndCflRate(self._onDevice, dt, cflRateByCell._onDevice)
   end
end

return VlasovDG
