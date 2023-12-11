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
struct gkyl_dg_vlasov_poisson_auxfields {
  const struct gkyl_array *field; // q/m*(phi_tot,A_ext)
};

// Object type
typedef struct gkyl_dg_updater_vlasov_poisson gkyl_dg_updater_vlasov_poisson;

/**
 * Create new updater to update Vlasov Poisson equations using hyper dg.
 *
 * @param grid Grid object
 * @param cbasis Configuration space basis functions
 * @param pbasis Phase space basis function
 * @param conf_range Configuration space range
 * @param vel_range Velocity space range
 * @param vmap Velocity mapping projected onto DG basis.
 * @param is_zero_flux_dir True in directions with (lower and upper) zero flux BCs.
 * @param field_id Enum identifier for field type (e.g., Maxwell's, Poisson, see gkyl_eqn_type.h)
 * @param aux_inp Void pointer to auxiliary fields. Void to be flexible to different auxfields structs
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 *
 * @return New vlasov updater object
 */
gkyl_dg_updater_vlasov_poisson* gkyl_dg_updater_vlasov_poisson_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range, const struct gkyl_array *vmap,
  const bool *is_zero_flux_dir, enum gkyl_field_id field_id, void *aux_inp, bool use_gpu);

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
void gkyl_dg_updater_vlasov_poisson_advance(gkyl_dg_updater_vlasov_poisson *vlasov,
  const struct gkyl_range *update_rng, const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

/**
 * Delete updater.
 *
 * @param vlasov Updater to delete.
 */
void gkyl_dg_updater_vlasov_poisson_release(gkyl_dg_updater_vlasov_poisson *vlasov);
]]

-- Vlasov DG solver updater object
local VlasovPoissonDG = Proto(UpdaterBase)

function VlasovPoissonDG:init(tbl)
   VlasovPoissonDG.super.init(self, tbl)

   -- Read data from input file.
   self._onGrid = assert(tbl.onGrid, "Updater.VlasovPoissonDG: Must provide grid object using 'onGrid'")
   self._phaseBasis = assert(
     tbl.phaseBasis, "Updater.VlasovPoissonDG: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(
     tbl.confBasis, "Updater.VlasovPoissonDG: Must specify conf-space basis functions to use using 'confBasis'")

   self._confRange = assert(
     tbl.confRange, "Updater.VlasovPoissonDG: Must specify conf-space range using 'confRange'")
   assert(self._confRange:isSubRange()==1, "Updater.VlasovPoissonDG: confRange must be a sub-range") 

   self._velRange = assert(
     tbl.velRange, "Updater.VlasovPoissonDG: Must specify velocity-space range using 'velRange'")
   assert(self._velRange:isSubRange()==1, "Updater.VlasovPoissonDG: velRange must be a sub-range") 

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   self._fieldId = assert(tbl.field_id, "Updater.VlasovPoissonDG: Must provide field ID using 'field_id'")

   -- Velocity coordinate mapping used for nonuniform velocity grids.
   local vMap = nil
   if tbl.velocityMap then
      vMap = self._useGPU and tbl.velocityMap._zeroDevice or tbl.velocityMap._zero
   end

   local field_id
   if     self._fieldId == "GKYL_FIELD_PHI"    then field_id = 1
   elseif self._fieldId == "GKYL_FIELD_PHI_A"  then field_id = 2
   end 

   -- Create aux fields struct for Vlasov and set pointers based on input fields
   -- Only supports Vlasov-Maxwell, Vlasov-Poisson (no external B), and SR Vlasov for now
   self._auxfields = ffi.new("struct gkyl_dg_vlasov_poisson_auxfields")
   if self._useGPU then
      self._auxfields.field = tbl.fldPtrs[1]._zeroDevice
   else
      self._auxfields.field = tbl.fldPtrs[1]._zero
   end

   local cdim, pdim = self._confBasis:ndim(), self._phaseBasis:ndim()
   local is_zfd = Lin.BoolVec(pdim)
   for d = 1, pdim do is_zfd[d] = d>cdim and true or false end
   local zfd = tbl.zeroFluxDirs -- Directions in which to specify zero flux BCs.
   if zfd then
      for i = 1, #zfd do is_zfd[zfd[i]] = true end
   end

   self._zero = ffi.gc(
      ffiC.gkyl_dg_updater_vlasov_poisson_new(self._onGrid._zero, 
        self._confBasis._zero, self._phaseBasis._zero, 
        self._confRange, self._velRange, vMap,
        is_zfd:data(), field_id, self._auxfields, self._useGPU),
      ffiC.gkyl_dg_updater_vlasov_poisson_release
   )

   return self
end

-- advance method
function VlasovPoissonDG:_advance(tCurr, inFld, outFld)

   local qIn = assert(inFld[1], "VlasovPoissonDG.advance: Must specify an input field")
   
   local qRhsOut = assert(outFld[1], "VlasovPoissonDG.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2],
     "VlasovPoissonDG.advance: Must pass cflRate field in output table")

   local localRange = qRhsOut:localRange()
   ffiC.gkyl_dg_updater_vlasov_poisson_advance(self._zero, localRange,
     qIn._zero, cflRateByCell._zero, qRhsOut._zero)

end

function VlasovPoissonDG:_advanceOnDevice(tCurr, inFld, outFld)

   local qIn = assert(inFld[1], "VlasovPoissonDG.advance: Must specify an input field")
   
   local qRhsOut = assert(outFld[1], "VlasovPoissonDG.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2],
     "VlasovPoissonDG.advance: Must pass cflRate field in output table")

   local localRange = qRhsOut:localRange()
   ffiC.gkyl_dg_updater_vlasov_poisson_advance(self._zero, localRange,
     qIn._zeroDevice, cflRateByCell._zeroDevice, qRhsOut._zeroDevice)

end

-- Fetch equation updater.
function VlasovPoissonDG:getEquation() return self._zero end

-- set up pointers to dt and cflRateByCell
function VlasovPoissonDG:setDtAndCflRate(dt, cflRateByCell)
   VlasovPoissonDG.super.setDtAndCflRate(self, dt, cflRateByCell)

   if self._onDevice then
      ffiC.setDtAndCflRate(self._onDevice, dt, cflRateByCell._onDevice)
   end
end

return VlasovPoissonDG
