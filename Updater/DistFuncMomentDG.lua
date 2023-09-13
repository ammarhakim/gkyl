-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments or integrated moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Lin          = require "Lib.Linalg"
local Proto        = require "Lib.Proto"
local UpdaterBase  = require "Updater.Base"
local lume         = require "Lib.lume"
local ffi         = require "ffi"
local xsys         = require "xsys"

local ffiC = ffi.C
require "Lib.ZeroUtil"

ffi.cdef [[
// Struct containing the pointers to auxiliary fields.
struct gkyl_mom_vlasov_sr_auxfields { 
  const struct gkyl_array *p_over_gamma; // p/gamma (velocity)
  const struct gkyl_array *gamma; // gamma = sqrt(1 + p^2)
  const struct gkyl_array *gamma_inv; // gamma_inv = 1/gamma = 1/sqrt(1 + p^2)
  const struct gkyl_array *V_drift; // bulk fluid velocity (computed from M0*V_drift = M1i with weak division)
  const struct gkyl_array *GammaV2; // Gamma^2 = 1/(1 - V_drift^2/c^2), Lorentz boost factor squared from bulk fluid velocity
  const struct gkyl_array *GammaV_inv; // Gamma_inv = sqrt(1 - V_drift^2/c^2), inverse Lorentz boost factor from bulk fluid velocity
};

// Object type
typedef struct gkyl_dg_updater_moment gkyl_dg_updater_moment;  

/**
 * Create new updater to compute moments of distribution function.
 * Supports Vlasov-Maxwell, special relativistic Vlasov-Maxwell,
 * and parallel-kinetic-perpendicular-moment (pkpm) Vlasov
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase-space basis function
 * @param conf_range Config space range
 * @param vel_range Velocity space range
 * @param model_id Enum identifier for model type (e.g., SR, PKPM, see gkyl_eqn_type.h)
 * @param aux_inp Void pointer to auxiliary fields. Void to be flexible to different auxfields structs
 * @param is_integrated Boolean for if the moment is an integrated moment
 * @param mom Name of moment
 * @param mass Mass of species 
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 * 
 * @return New moment updater object
 */
struct gkyl_dg_updater_moment*
gkyl_dg_updater_moment_new(const struct gkyl_rect_grid *grid, 
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, 
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range,
  enum gkyl_model_id model_id, void *aux_inp, 
  const char *mom, bool is_integrated, double mass, bool use_gpu);

/**
 * Compute moment. The update_phase_rng and update_conf_rng MUST be a sub-range of the
 * be a sub-range of the range on which the array is defined. 
 * That is, it must be either the same range as the array range, or one created using the
 * or one created using the gkyl_sub_range_init method.
 *
 * @param moment moment updater object
 * @param update_phase_rng Phase space range on which to compute.
 * @param update_conf_rng Configuration space range on which to compute.
 * @param fIn Input to updater
 * @param mout Output moment
 */
void
gkyl_dg_updater_moment_advance(struct gkyl_dg_updater_moment *moment,
  const struct gkyl_range *update_phase_rng, const struct gkyl_range *update_conf_rng,
  const struct gkyl_array* fIn, struct gkyl_array* mout);

/**
 * Delete updater.
 *
 * @param moment Updater to delete.
 */
void gkyl_dg_updater_moment_release(struct gkyl_dg_updater_moment* moment);

]]

-- Moment updater object
local DistFuncMomentDG = Proto(UpdaterBase)

function DistFuncMomentDG:init(tbl)
   DistFuncMomentDG.super.init(self, tbl)    -- Setup base object.

   self._onGrid = assert(tbl.onGrid, "Updater.DistFuncMomentDG: Must specify grid to use with 'onGrid'.")
   self._confBasis = assert(tbl.confBasis, "Updater.DistFuncMomentDG: Must specify configuration space basis with 'confBasis'.")
   self._phaseBasis = assert(tbl.phaseBasis, "Updater.DistFuncMomentDG: Must specify phase space basis with 'confBasis'.")

   -- Get the ranges needed for certain moment routines if they were passed as inputs
   self._velRange = nil
   if tbl.velRange then 
      self._velRange = tbl.velRange
      assert(self._velRange:isSubRange()==1, "Updater.DistFuncMomentDG: velRange must be a sub-range") 
   end
   self._confRange = nil
   if tbl.confRange then 
      self._confRange = tbl.confRange
      assert(self._confRange:isSubRange()==1, "Updater.DistFuncMomentDG: confRange must be a sub-range") 
   end

   self._moment = assert(tbl.moment, "Updater.DistFuncMomentDG: Must specify moment type with 'moment'.")
   self._isIntegrated = tbl.isIntegrated or false

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU)

   self._modelId = assert(tbl.model_id, "Updater.DistFuncMomentDG: Must provide model ID using 'model_id'")
   
   local model_id
   if self._modelId == "GKYL_MODEL_DEFAULT" then model_id = 0
   elseif self._modelId == "GKYL_MODEL_SR"  then model_id = 1
   elseif self._modelId == "GKYL_MODEL_GEN_GEO"  then model_id = 2 
   elseif self._modelId == "GKYL_MODEL_PKPM"  then model_id = 3
   elseif self._modelId == "GKYL_MODEL_SR_PKPM"  then model_id = 4
   end 

   if model_id == 1 then 
      -- Create aux fields struct for special relativistic moments
      self._auxfields = ffi.new("struct gkyl_mom_vlasov_sr_auxfields")
      if useGPU then
         self._auxfields.p_over_gamma = tbl.momPtrs[1]._zeroDevice
         self._auxfields.gamma = tbl.momPtrs[2]._zeroDevice
         self._auxfields.gamma_inv = tbl.momPtrs[3]._zeroDevice
      else
         self._auxfields.p_over_gamma = tbl.momPtrs[1]._zero
         self._auxfields.gamma = tbl.momPtrs[2]._zero
         self._auxfields.gamma_inv = tbl.momPtrs[3]._zero
      end
   end

   self._zero = ffi.gc(ffiC.gkyl_dg_updater_moment_new(self._onGrid._zero, self._confBasis._zero, self._phaseBasis._zero, 
                       self._confRange, self._velRange, model_id, self._auxfields, 
                       self._moment, self._isIntegrated, 0.0, useGPU or 0),
                       ffiC.gkyl_dg_updater_moment_release)
end

-- Advance method.
function DistFuncMomentDG:_advance(tCurr, inFld, outFld)
   local qIn = assert(inFld[1], "DistFuncMomentDG.advance: Must specify an input field")
   local mOut = assert(outFld[1], "DistFuncMomentDG.advance: Must specify an output field")

   local confRange = inFld[2] or mOut:localExtRange()
   local phaseRange = inFld[3] or qIn:localExtRange()

   ffiC.gkyl_dg_updater_moment_advance(self._zero, phaseRange, confRange, qIn._zero, mOut._zero)
end

function DistFuncMomentDG:_advanceOnDevice(tCurr, inFld, outFld)
   local qIn = assert(inFld[1], "DistFuncMomentDG.advance: Must specify an input field")
   local mOut = assert(outFld[1], "DistFuncMomentDG.advance: Must specify an output field")

   local confRange = inFld[2] or mOut:localExtRange()
   local phaseRange = inFld[3] or qIn:localExtRange()
   
   ffiC.gkyl_dg_updater_moment_advance(self._zero, phaseRange, confRange, qIn._zeroDevice, mOut._zeroDevice)
end


return DistFuncMomentDG
