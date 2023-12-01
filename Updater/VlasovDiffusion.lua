--------------------------------------------------------------------------------
--
-- Wrap g0's dg_updater_diffusion_vlasov to apply diffusion using DG.
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
// Object type
typedef struct gkyl_dg_updater_diffusion_vlasov gkyl_dg_updater_diffusion_vlasov;

// return type for drag and diffusion timers
struct gkyl_dg_updater_diffusion_vlasov_tm {
  double diffusion_tm; // time for diffusion updates
};

/**
 * Create new updater to update diffusion equations using hyper dg or hyper dg gen stencil.
 *
 * @param grid Grid object.
 * @param basis Basis functions of the equation system.
 * @param cbasis Configuration space basis.
 * @param is_diff_constant If diffusion coefficient spatially constant.
 * @param diff_in_dir Whether to apply diffusion in each direction.
 * @param diff_order Diffusion order.
 * @param diff_range Range object to index the diffusion coefficient.
 * @param is_zero_flux_dir True in directions with (lower and upper) zero flux BCs.
 * @param use_gpu Whether to run on host or device.
 * @return New diff updater object
 */
struct gkyl_dg_updater_diffusion_vlasov* gkyl_dg_updater_diffusion_vlasov_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *basis, const struct gkyl_basis *cbasis, bool is_diff_const, const bool *diff_in_dir,
  int diff_order, const struct gkyl_range *diff_range, const bool *is_zero_flux_dir, bool use_gpu);

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
void gkyl_dg_updater_diffusion_vlasov_advance(struct gkyl_dg_updater_diffusion_vlasov *up,
  const struct gkyl_range *update_rng, const struct gkyl_array *coeff,
  const struct gkyl_array* fIn, struct gkyl_array* cflrate,
  struct gkyl_array* rhs);

/**
 * Return total time spent in diffusion terms
 *
 * @param diffusion Updater object
 * @return timers
 */
struct gkyl_dg_updater_diffusion_vlasov_tm gkyl_dg_updater_diffusion_vlasov_get_tm(const struct gkyl_dg_updater_diffusion_vlasov *up);

/**
 * Delete updater.
 *
 * @param diffusion Updater to delete.
 */
void gkyl_dg_updater_diffusion_vlasov_release(struct gkyl_dg_updater_diffusion_vlasov *up);
]]

-- Diffusion DG solver updater object.
local VlasovDiffusion = Proto(UpdaterBase)

function VlasovDiffusion:init(tbl)
   VlasovDiffusion.super.init(self, tbl)

   -- Read data from input file.
   local onGrid = assert(tbl.onGrid, "Updater.VlasovDiffusion: Must provide grid object using 'onGrid'")
   local onBasis = assert(tbl.onBasis, "Updater.VlasovDiffusion: Must specify basis functions to use using 'onBasis'")
   local confBasis = tbl.confBasis or onBasis

   local constDiff = tbl.constDiff
   assert(tbl.constDiff~=nil, "Updater.VlasovDiffusion: Must indicate if the coefficient is constant in 'constDiff'")

   local diffDirs = tbl.directions  -- Directions to apply diffusion in.
   local diffOrder = tbl.order and tbl.order or 2 -- Order of diffusion to apply.

   local confRange = tbl.confRange  -- Used to index a spatially varying diffusion coefficient.

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   local cdim = confBasis:ndim()
   local indirs, indirsPtr
   if diffDirs then
      indirs = Lin.BoolVec(cdim)
      for i=1,cdim do indirs[i] = false end
      for i=1,#diffDirs do
         assert((diffDirs[i]>0) and (diffDirs[i]<cdim+1), "Updater.VlasovDiffusion: diffusive directions must be >1 and <cdim+1.")
         indirs[diffDirs[i]] = true
      end
      indirsPtr = indirs:data()
   end

   local pdim = onBasis:ndim()
   local is_zfd = Lin.BoolVec(pdim)
   for d = 1, pdim do is_zfd[d] = false end 
   local zfd = tbl.zeroFluxDirs -- Directions in which to specify zero flux BCs.
   if zfd then
      for i = 1, #zfd do is_zfd[zfd[i]] = true end
   end

   self._zero = ffi.gc(ffiC.gkyl_dg_updater_diffusion_vlasov_new(onGrid._zero, onBasis._zero,
                       confBasis._zero, constDiff, indirsPtr, diffOrder, confRange, is_zfd:data(), useGPU),
                       ffiC.gkyl_dg_updater_diffusion_vlasov_release)

   return self
end

-- Advance method.
function VlasovDiffusion:_advance(tCurr, inFld, outFld)
   local fIn   = assert(inFld[1], "VlasovDiffusion.advance: Must specify an input field.")
   local diffD = assert(inFld[2], "VlasovDiffusion.advance: Must specify the diffusion coefficient (a CartField).")

   local fRhsOut       = assert(outFld[1], "VlasovDiffusion.advance: Must specify an output field.")
   local cflRateByCell = assert(outFld[2], "VlasovDiffusion.advance: Must pass cflRate field in output table.")

   local localRange = fRhsOut:localRange()

   ffiC.gkyl_dg_updater_diffusion_vlasov_advance(self._zero, localRange, diffD._zero, fIn._zero, cflRateByCell._zero, fRhsOut._zero)
end

function VlasovDiffusion:_advanceOnDevice(tCurr, inFld, outFld)
   local fIn   = assert(inFld[1], "VlasovDiffusion.advance: Must specify an input field.")
   local diffD = assert(inFld[2], "VlasovDiffusion.advance: Must specify the diffusion coefficient (a CartField).")

   local fRhsOut       = assert(outFld[1], "VlasovDiffusion.advance: Must specify an output field.")
   local cflRateByCell = assert(outFld[2], "VlasovDiffusion.advance: Must pass cflRate field in output table.")

   local localRange = fRhsOut:localRange()

   ffiC.gkyl_dg_updater_diffusion_vlasov_advance(self._zero, localRange, diffD._zeroDevice, fIn._zeroDevice, cflRateByCell._zeroDevice, fRhsOut._zeroDevice) 
end

return VlasovDiffusion
