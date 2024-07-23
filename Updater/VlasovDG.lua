--
-- This updater wraps g0's dg_updater_vlasov to advance the Vlasov
-- equations with Discontinuous Galerkin scheme.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto       = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local Lin         = require "Lib.Linalg"
local ffi         = require "ffi"
local xsys        = require "xsys"

local ffiC = ffi.C
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

ffi.cdef [[ 
// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_vlasov_auxfields { 
  const struct gkyl_array *field; // q/m*(E,B) for Maxwell's, q/m*phi for Poisson's (gradient calculated in kernel)
  const struct gkyl_array *cot_vec; // cotangent vectors (e^i) used in volume term if general geometry enabled
  const struct gkyl_array *alpha_geo; // alpha^i (e^i . alpha) used in surface term if general geometry enabled
};

// Struct containing the pointers to auxiliary fields.
struct gkyl_dg_vlasov_sr_auxfields { 
  const struct gkyl_array *qmem; // q/m * EM
  const struct gkyl_array *gamma; // gamma = sqrt(1 + p^2), particle Lorentz factor
};

// Object type
typedef struct gkyl_dg_updater_vlasov gkyl_dg_updater_vlasov;

/**
 * Create new updater to update vlasov equations using hyper dg.
 * Supports Vlasov-Maxwell, Vlasov-Poisson (with and without vector potential A)
 * and special relativistic Vlasov-Maxwell
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase space basis function
 * @param conf_range Configuration space range
 * @param vel_range Velocity space range
 * @param phase_range Phase space range
 * @param is_zero_flux_dir True in directions with (lower and upper) zero flux BCs.
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
  const bool *is_zero_flux_dir, enum gkyl_model_id model_id, enum gkyl_field_id field_id, void *aux_inp, bool use_gpu);

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
   self._phaseBasis = assert(
     tbl.phaseBasis, "Updater.VlasovDG: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(
     tbl.confBasis, "Updater.VlasovDG: Must specify conf-space basis functions to use using 'confBasis'")

   self._confRange = assert(
     tbl.confRange, "Updater.VlasovDG: Must specify conf-space range using 'confRange'")
   assert(self._confRange:isSubRange()==1, "Updater.VlasovDG: confRange must be a sub-range") 

   self._velRange = assert(
     tbl.velRange, "Updater.VlasovDG: Must specify velocity-space range using 'velRange'")
   assert(self._velRange:isSubRange()==1, "Updater.VlasovDG: velRange must be a sub-range") 

   self._phaseRange = assert(
     tbl.phaseRange, "Updater.VlasovDG: Must specify phase-space range using 'phaseRange'")
   assert(self._phaseRange:isSubRange()==1, "Updater.VlasovDG: phaseRange must be a sub-range") 

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   self._modelId = assert(tbl.model_id, "Updater.VlasovDG: Must provide model ID using 'model_id'")
   self._fieldId = assert(tbl.field_id, "Updater.VlasovDG: Must provide field ID using 'field_id'")

   local model_id
   if self._modelId == "GKYL_MODEL_DEFAULT" then model_id = 0
   elseif self._modelId == "GKYL_MODEL_SR"  then model_id = 1
   elseif self._modelId == "GKYL_MODEL_GEN_GEO"  then model_id = 2 
   end   

   local field_id
   if self._fieldId == "GKYL_FIELD_E_B" then field_id = 0
   elseif self._fieldId == "GKYL_FIELD_PHI"  then field_id = 1
   elseif self._fieldId == "GKYL_FIELD_PHI_A"  then field_id = 2
   elseif self._fieldId == "GKYL_FIELD_NULL"  then field_id = 3
   end 

   -- Create aux fields struct for Vlasov and set pointers based on input fields
   -- Only supports Vlasov-Maxwell, Vlasov-Poisson (no external B), and SR Vlasov for now
   if model_id == 1 then 
      self._auxfields = ffi.new("struct gkyl_dg_vlasov_sr_auxfields")
      if self._useGPU then
         self._auxfields.qmem = tbl.fldPtrs[1]._zeroDevice
         self._auxfields.gamma = tbl.fldPtrs[2]._zeroDevice
      else
         self._auxfields.qmem = tbl.fldPtrs[1]._zero
         self._auxfields.gamma = tbl.fldPtrs[2]._zero
      end
   else
      self._auxfields = ffi.new("struct gkyl_dg_vlasov_auxfields")
      if self._useGPU then
         self._auxfields.field = tbl.fldPtrs[1]._zeroDevice
      else
         self._auxfields.field = tbl.fldPtrs[1]._zero
      end
   end

   local cdim, pdim = self._confBasis:ndim(), self._phaseBasis:ndim()
   local is_zfd = Lin.BoolVec(pdim)
   for d = 1, pdim do is_zfd[d] = d>cdim and true or false end
   local zfd = tbl.zeroFluxDirs -- Directions in which to specify zero flux BCs.
   if zfd then
      for i = 1, #zfd do is_zfd[zfd[i]] = true end
   end

   self._zero = ffi.gc(
      ffiC.gkyl_dg_updater_vlasov_new(self._onGrid._zero, 
        self._confBasis._zero, self._phaseBasis._zero, 
        self._confRange, self._velRange, self._phaseRange,
        is_zfd:data(), model_id, field_id, self._auxfields, self._useGPU),
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
