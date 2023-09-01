--------------------------------------------------------------------------------
--
-- Wrap g0's dg_updater_diffusion to apply diffusion using DG.
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
// Identifiers for specific diffusion object types
enum gkyl_diffusion_id {
  GKYL_DIFFUSION_NONE = 0,  // No diffusion. This is default.
  GKYL_DIFFUSION_GEN,       // General diffusion tensor.
  GKYL_DIFFUSION_DIAGONAL_CONST,  // Diagonal const (hyper)diffusion for a scalar equation.
  GKYL_DIFFUSION_DIAGONAL_CONST_EULER,       // Euler.
  GKYL_DIFFUSION_DIAGONAL_CONST_EULER_ISO,   // Isothermal Euler.
  GKYL_DIFFUSION_DIAGONAL_CONST_VLASOV,      // Vlasov.
  GKYL_DIFFUSION_DIAGONAL_CONST_GYROKINETIC, // Gyrokinetic.
  GKYL_DIFFUSION_DIAGONAL_CONST_PKPM,        // PKPM.
  GKYL_DIFFUSION_DIAGONAL_VAR,  // Diagonal spatially varying (hyper)diffusion for a scalar equation.
  GKYL_DIFFUSION_DIAGONAL_VAR_EULER,       // Euler.
  GKYL_DIFFUSION_DIAGONAL_VAR_EULER_ISO,   // Isothermal Euler.
  GKYL_DIFFUSION_DIAGONAL_VAR_VLASOV,      // Vlasov.
  GKYL_DIFFUSION_DIAGONAL_VAR_GYROKINETIC, // Gyrokinetic.
  GKYL_DIFFUSION_DIAGONAL_VAR_PKPM,        // PKPM.
};

// Object type
typedef struct gkyl_dg_updater_diffusion gkyl_dg_updater_diffusion;

// return type for drag and diffusion timers
struct gkyl_dg_updater_diffusion_tm {
  double diffusion_tm; // time for diffusion updates
};

/**
 * Create new updater to update diffusion equations using hyper dg or hyper dg gen stencil.
 *
 * @param grid Grid object.
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param diffusion_id Diffusion type (constant/varying & fluid/vlasov/etc).
 * @param diff_in_dir Whether to apply diffusion in each direction.
 * @param diff_order Diffusion order.
 * @param conf_range Configuration space range (to index diff coefficient).
 * @param use_gpu Whether to run on host or device.
 * @return New diff updater object
 */
struct gkyl_dg_updater_diffusion* gkyl_dg_updater_diffusion_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis,
  enum gkyl_diffusion_id diffusion_id, bool *diff_in_dir, int diff_order,
  const struct gkyl_range *conf_range, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up diffusion updater object.
 * @param update_rng Range on which to compute.
 * @param coeff Diffusion coefficient/tensor.
 * @param fIn Input to updater.
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param rhs RHS output
 */
void gkyl_dg_updater_diffusion_advance(struct gkyl_dg_updater_diffusion *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *coeff,
  const struct gkyl_array* fIn, struct gkyl_array* cflrate,
  struct gkyl_array* rhs);

/**
 * Return total time spent in diffusion terms
 *
 * @param diffusion Updater object
 * @return timers
 */
struct gkyl_dg_updater_diffusion_tm gkyl_dg_updater_diffusion_get_tm(const struct gkyl_dg_updater_diffusion *up);

/**
 * Delete updater.
 *
 * @param diffusion Updater to delete.
 */
void gkyl_dg_updater_diffusion_release(struct gkyl_dg_updater_diffusion *up);
]]

-- Diffusion DG solver updater object.
local DiffusionDG = Proto(UpdaterBase)

function DiffusionDG:init(tbl)
   DiffusionDG.super.init(self, tbl)

   -- Read data from input file.
   local onGrid = assert(tbl.onGrid, "Updater.DiffusionDG: Must provide grid object using 'onGrid'")
   local onBasis = assert(tbl.onBasis, "Updater.DiffusionDG: Must specify basis functions to use using 'onBasis'")
   local confBasis = tbl.confBasis or onBasis

   local diffType = assert(tbl.modelType, "Updater.DiffusionDG: Must indicate the diffusion model type using 'modelType'")
   local diffDirs = tbl.directions  -- Directions to apply diffusion in.
   local diffOrder = tbl.order and tbl.order or 2 -- Order of diffusion to apply.

   local confRange = tbl.confRange  -- Used to index a spatially varying diffusion coefficient.

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   local diffusionId
   if     diffType == 'general'                       then diffusionId = 1
   elseif diffType == 'diagonal_constant'             then diffusionId = 2
   elseif diffType == 'diagonal_constant_euler'       then diffusionId = 3
   elseif diffType == 'diagonal_constant_euler_iso'   then diffusionId = 4
   elseif diffType == 'diagonal_constant_vlasov'      then diffusionId = 5
   elseif diffType == 'diagonal_constant_gyrokinetic' then diffusionId = 6
   elseif diffType == 'diagonal_constant_pkpm'        then diffusionId = 7
   elseif diffType == 'diagonal_varying'              then diffusionId = 8
   elseif diffType == 'diagonal_varying_euler'        then diffusionId = 9
   elseif diffType == 'diagonal_varying_euler_iso'    then diffusionId = 10
   elseif diffType == 'diagonal_varying_vlasov'       then diffusionId = 11
   elseif diffType == 'diagonal_varying_gyrokinetic'  then diffusionId = 12
   elseif diffType == 'diagonal_varying_pkpm'         then diffusionId = 13
   end

   if diffusionId>7 then assert(confRange:isSubRange()==1, "Eq.Diffusion: confRange must be a sub-range") end

   local cdim = confBasis:ndim()
   local indirs, indirsPtr
   if diffDirs then
      indirs = Lin.BoolVec(cdim)
      for i=1,cdim do indirs[i] = false end
      for i=1,#diffDirs do
         assert((diffDirs[i]>0) and (diffDirs[i]<cdim+1), "Updater.DiffusionDG: diffusive directions must be >1 and <cdim+1.")
         indirs[diffDirs[i]] = true
      end
      indirsPtr = indirs:data()
   end

   self._zero = ffi.gc(ffiC.gkyl_dg_updater_diffusion_new(onGrid._zero, onBasis._zero,
                       confBasis._zero, diffusionId, indirsPtr, diffOrder, confRange, useGPU),
                       ffiC.gkyl_dg_updater_diffusion_release)

   return self
end

-- Advance method.
function DiffusionDG:_advance(tCurr, inFld, outFld)
   local fIn   = assert(inFld[1], "DiffusionDG.advance: Must specify an input field.")
   local diffD = assert(inFld[2], "DiffusionDG.advance: Must specify the diffusion coefficient (a CartField).")

   local fRhsOut       = assert(outFld[1], "DiffusionDG.advance: Must specify an output field.")
   local cflRateByCell = assert(outFld[2], "DiffusionDG.advance: Must pass cflRate field in output table.")

   local localRange = fRhsOut:localRange()

   ffiC.gkyl_dg_updater_diffusion_advance(self._zero, localRange, diffD._zero, fIn._zero, cflRateByCell._zero, fRhsOut._zero) 
end

function DiffusionDG:_advanceOnDevice(tCurr, inFld, outFld)
   local fIn   = assert(inFld[1], "DiffusionDG.advance: Must specify an input field.")
   local diffD = assert(inFld[2], "DiffusionDG.advance: Must specify the diffusion coefficient (a CartField).")

   local fRhsOut       = assert(outFld[1], "DiffusionDG.advance: Must specify an output field.")
   local cflRateByCell = assert(outFld[2], "DiffusionDG.advance: Must pass cflRate field in output table.")

   local localRange = fRhsOut:localRange()

   ffiC.gkyl_dg_updater_diffusion_advance(self._zero, localRange, diffD._zeroDevice, fIn._zeroDevice, cflRateByCell._zeroDevice, fRhsOut._zeroDevice) 
end

return DiffusionDG
