-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments of distribution function.
--
-- For LBO collisions this updater also computes the boundary corrections and,
-- if using a piecewise polynomial basis, the star moments.
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

typedef struct gkyl_dg_updater_vlasov gkyl_dg_updater_vlasov;

// Object type
typedef struct gkyl_dg_updater_bflux_vlasov gkyl_dg_updater_bflux_vlasov;

gkyl_dg_updater_bflux_vlasov* gkyl_dg_updater_bflux_vlasov_new(const struct gkyl_rect_grid *grid, 
  int cdim, const gkyl_dg_updater_vlasov *vlasov, bool use_gpu);

void gkyl_dg_updater_bflux_vlasov_advance(gkyl_dg_updater_bflux_vlasov *bflux, const struct gkyl_range *update_rng, const struct gkyl_array* fIn, struct gkyl_array* rhs);

void gkyl_dg_updater_bflux_vlasov_advance_cu(gkyl_dg_updater_bflux_vlasov *bflux, const struct gkyl_range *update_rng, const struct gkyl_array* fIn, struct gkyl_array* rhs);

void gkyl_dg_updater_bflux_vlasov_release(gkyl_dg_updater_bflux_vlasov* bflux);
]]

-- Boundary fluxes updater object.
local BoundaryFluxCalc = Proto(UpdaterBase)

function BoundaryFluxCalc:init(tbl)
   BoundaryFluxCalc.super.init(self, tbl)    -- Setup base object.

   self._onGrid = assert(tbl.onGrid, "Updater.BoundaryFluxCalc: Must specify grid to use with 'onGrid'.")
   self._cdim = assert(tbl.cdim, "Updater.BoundaryFluxCalc: Must specify number of configuration space dimensions to use with 'cdim'.")
   self._equation_id = assert(tbl.equation_id, "Updater.BoundaryFluxCalc: Must specify equation type to use with 'equation_id'.")

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU)

   self._equation = assert(tbl.equation, "Updater.BoundaryFluxCalc: Must specify equation to use with 'equation'.")

   if self._equation_id == "vlasov" then
      self._zero = ffi.gc(ffiC.gkyl_dg_updater_bflux_vlasov_new(self._onGrid._zero, self._cdim, self._equation, useGPU or 0),
                          ffiC.gkyl_dg_updater_bflux_vlasov_release)
   end
end

-- Advance method.
function BoundaryFluxCalc:_advance(tCurr, inFld, outFld)
   local qIn = assert(inFld[1], "BoundaryFluxCalc.advance: Must specify an input field")
   local qRhsOut = assert(outFld[1], "BoundaryFluxCalc.advance: Must specify an output field")
   local phaseRange = qRhsOut:localExtRange()

   ffiC.gkyl_dg_updater_bflux_vlasov_advance(self._zero, phaseRange, qIn._zero, qRhsOut._zero)
end

function BoundaryFluxCalc:_advanceOnDevice(tCurr, inFld, outFld)
   local qIn = assert(inFld[1], "BoundaryFluxCalc.advance: Must specify an input field")
   local qRhsOut = assert(outFld[1], "BoundaryFluxCalc.advance: Must specify an output field")
   local phaseRange = qRhsOut:localExtRange()

   ffiC.gkyl_dg_updater_bflux_vlasov_advance_cu(self._zero, phaseRange, qIn._zeroDevice, qRhsOut._zeroDevice)
end


return BoundaryFluxCalc
