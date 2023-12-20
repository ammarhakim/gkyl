-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute boundary fluxes in the ghost cells.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Lin          = require "Lib.Linalg"
local Proto        = require "Lib.Proto"
local UpdaterBase  = require "Updater.Base"
local lume         = require "Lib.lume"
local ffi          = require "ffi"
local xsys         = require "xsys"

local ffiC = ffi.C
require "Lib.ZeroUtil"

ffi.cdef [[
// Dependent object type
typedef struct gkyl_dg_updater_vlasov gkyl_dg_updater_vlasov;

// Object type
typedef struct gkyl_dg_updater_bflux_vlasov gkyl_dg_updater_bflux_vlasov;

/**
 * Create new updater to compute the boundary flux contributions from
 * Vlasov equations. These are DG surface terms using hyper dg.
 * Supports Vlasov-Maxwell, special relativistic Vlasov-Maxwell and
 * neutral Vlasov.
 *
 * @param grid Grid object.
 * @param cdim Number of configuration space dimensions.
 * @param vlasov objects which performs the Vlasov DG update.
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 *
 * @return New boundary flux updater object.
 */
gkyl_dg_updater_bflux_vlasov* gkyl_dg_updater_bflux_vlasov_new(const struct gkyl_rect_grid *grid,
  int cdim, const gkyl_dg_updater_vlasov *vlasov, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Boundary flux updater.
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater.
 * @param rhs RHS output.
 */
void gkyl_dg_updater_bflux_vlasov_advance(gkyl_dg_updater_bflux_vlasov *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* fIn, struct gkyl_array* rhs);

void gkyl_dg_updater_bflux_vlasov_advance_cu(gkyl_dg_updater_bflux_vlasov *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* fIn, struct gkyl_array* rhs);

/**
 * Delete updater.
 *
 * @param up Boundary flux updater.
 */
void gkyl_dg_updater_bflux_vlasov_release(gkyl_dg_updater_bflux_vlasov* up);

// Dependent object type
typedef struct gkyl_dg_updater_vlasov_poisson gkyl_dg_updater_vlasov_poisson;

// Object type
typedef struct gkyl_dg_updater_bflux_vlasov_poisson gkyl_dg_updater_bflux_vlasov_poisson;

/**
 * Create new updater to compute the boundary flux contributions from
 * Vlasov Poisson equations. These are DG surface terms using hyper dg.
 *
 * @param grid Grid object.
 * @param cdim Number of configuration space dimensions.
 * @param vlasov objects which performs the Vlasov-Poisson DG update.
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 *
 * @return New boundary flux updater object.
 */
struct gkyl_dg_updater_bflux_vlasov_poisson* gkyl_dg_updater_bflux_vlasov_poisson_new(const struct gkyl_rect_grid *grid,
  int cdim, const gkyl_dg_updater_vlasov_poisson *vlasov, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Boundary flux updater.
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater.
 * @param rhs RHS output.
 */
void gkyl_dg_updater_bflux_vlasov_poisson_advance(struct gkyl_dg_updater_bflux_vlasov_poisson *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* fIn, struct gkyl_array* rhs);

void gkyl_dg_updater_bflux_vlasov_poisson_advance_cu(struct gkyl_dg_updater_bflux_vlasov_poisson *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* fIn, struct gkyl_array* rhs);

/**
 * Delete updater.
 *
 * @param up Boundary flux updater.
 */
void gkyl_dg_updater_bflux_vlasov_poisson_release(struct gkyl_dg_updater_bflux_vlasov_poisson *up);

// Dependent object type
typedef struct gkyl_dg_updater_gyrokinetic gkyl_dg_updater_gyrokinetic;

// Object type
typedef struct gkyl_dg_updater_bflux_gyrokinetic gkyl_dg_updater_bflux_gyrokinetic;

/**
 * Create new updater to compute the boundary flux contributions from
 * gyrokinetic equations. These are DG surface terms using hyper dg.
 *
 * @param grid Grid object.
 * @param cdim Number of configuration space dimensions.
 * @param gyrokinetic objects which performs the gyrokinetic DG update.
 * @param use_gpu Boolean to determine whether struct objects are on host or device
 *
 * @return New boundary flux updater object.
 */
gkyl_dg_updater_bflux_gyrokinetic* gkyl_dg_updater_bflux_gyrokinetic_new(const struct gkyl_rect_grid *grid,
  int cdim, const gkyl_dg_updater_gyrokinetic *gyrokinetic, bool use_gpu);

/**
 * Compute RHS of DG update. The update_rng MUST be a sub-range of the
 * range on which the array is defined. That is, it must be either the
 * same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param up Boundary flux updater.
 * @param update_rng Range on which to compute.
 * @param fIn Input to updater.
 * @param rhs RHS output.
 */
void gkyl_dg_updater_bflux_gyrokinetic_advance(gkyl_dg_updater_bflux_gyrokinetic *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* fIn, struct gkyl_array* rhs);

void gkyl_dg_updater_bflux_gyrokinetic_advance_cu(gkyl_dg_updater_bflux_gyrokinetic *up,
  const struct gkyl_range *update_rng,
  const struct gkyl_array* fIn, struct gkyl_array* rhs);

/**
 * Delete updater.
 *
 * @param up Boundary flux updater.
 */
void gkyl_dg_updater_bflux_gyrokinetic_release(gkyl_dg_updater_bflux_gyrokinetic* up);
]]

-- Boundary fluxes updater object.
local BoundaryFluxCalc = Proto(UpdaterBase)

function BoundaryFluxCalc:init(tbl)
   BoundaryFluxCalc.super.init(self, tbl)    -- Setup base object.

   local grid = assert(tbl.onGrid, "Updater.BoundaryFluxCalc: Must specify grid to use with 'onGrid'.")
   local cdim = assert(tbl.cdim, "Updater.BoundaryFluxCalc: Must specify number of configuration space dimensions to use with 'cdim'.")
   local equation_id = assert(tbl.equation_id, "Updater.BoundaryFluxCalc: Must specify equation type to use with 'equation_id'.")

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU)

   local equation = assert(tbl.equation, "Updater.BoundaryFluxCalc: Must specify equation to use with 'equation'.")

   if equation_id == "vlasov" then
      self._zero = ffi.gc(ffiC.gkyl_dg_updater_bflux_vlasov_new(grid._zero, cdim, equation, useGPU or 0),
                          ffiC.gkyl_dg_updater_bflux_vlasov_release)
      self._bfluxAdvance = function(updateRange, fIn, fOut)
         ffiC.gkyl_dg_updater_bflux_vlasov_advance(self._zero, updateRange, fIn._zero, fOut._zero)
      end
      self._bfluxAdvanceOnDevice = function(updateRange, fIn, fOut)
         ffiC.gkyl_dg_updater_bflux_vlasov_advance_cu(self._zero, updateRange, fIn._zeroDevice, fOut._zeroDevice)
      end
   elseif equation_id == "vlasov_poisson" then
      self._zero = ffi.gc(ffiC.gkyl_dg_updater_bflux_vlasov_poisson_new(grid._zero, cdim, equation, useGPU or 0),
                          ffiC.gkyl_dg_updater_bflux_vlasov_poisson_release)
      self._bfluxAdvance = function(updateRange, fIn, fOut)
         ffiC.gkyl_dg_updater_bflux_vlasov_poisson_advance(self._zero, updateRange, fIn._zero, fOut._zero)
      end
      self._bfluxAdvanceOnDevice = function(updateRange, fIn, fOut)
         ffiC.gkyl_dg_updater_bflux_vlasov_poisson_advance_cu(self._zero, updateRange, fIn._zeroDevice, fOut._zeroDevice)
      end
   elseif equation_id == "gyrokinetic" then
      self._zero = ffi.gc(ffiC.gkyl_dg_updater_bflux_gyrokinetic_new(grid._zero, cdim, equation, useGPU or 0),
                          ffiC.gkyl_dg_updater_bflux_gyrokinetic_release)
      self._bfluxAdvance = function(updateRange, fIn, fOut)
         ffiC.gkyl_dg_updater_bflux_gyrokinetic_advance(self._zero, updateRange, fIn._zero, fOut._zero)
      end
      self._bfluxAdvanceOnDevice = function(updateRange, fIn, fOut)
         ffiC.gkyl_dg_updater_bflux_gyrokinetic_advance_cu(self._zero, updateRange, fIn._zeroDevice, fOut._zeroDevice)
      end
   end
end

-- Advance method.
function BoundaryFluxCalc:_advance(tCurr, inFld, outFld)
   local qIn = assert(inFld[1], "BoundaryFluxCalc.advance: Must specify an input field")
   local qRhsOut = assert(outFld[1], "BoundaryFluxCalc.advance: Must specify an output field")
   local phaseRange = qRhsOut:localExtRange()

   self._bfluxAdvance(phaseRange, qIn, qRhsOut)
end

function BoundaryFluxCalc:_advanceOnDevice(tCurr, inFld, outFld)
   local qIn = assert(inFld[1], "BoundaryFluxCalc.advance: Must specify an input field")
   local qRhsOut = assert(outFld[1], "BoundaryFluxCalc.advance: Must specify an output field")
   local phaseRange = qRhsOut:localExtRange()

   self._bfluxAdvanceOnDevice(phaseRange, qIn, qRhsOut)
end


return BoundaryFluxCalc
