--
-- This updater wraps g0's dg_updater_fluid to advance the fluid
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

 local cuda = nil
 if GKYL_HAVE_CUDA then
    cuda    = require "Cuda.RunTime"
    cuAlloc = require "Cuda.Alloc"
 end

ffi.cdef [[

// Identifiers for various equation systems TODO: we should move this to some equation type lua file in eq folder as it will be used multiple times
enum gkyl_eqn_type {
  GKYL_EQN_EULER, // Euler equations
  GKYL_EQN_SR_EULER, // SR Euler equations
  GKYL_EQN_ISO_EULER, // Isothermal Euler equations
  GKYL_EQN_TEN_MOMENT, // Ten-moment (with pressure tensor)
  GKYL_EQN_MAXWELL, // Maxwell equations
  GKYL_EQN_MHD,  // Ideal MHD equations
  GKYL_EQN_BURGERS, // Burgers equations
  GKYL_EQN_ADVECTION, // Scalar advection equation
  GKYL_EQN_EULER_PKPM, // Euler equations with parallel-kinetic-perpendicular-moment (pkpm) model
};

typedef struct gkyl_dg_updater_fluid gkyl_dg_updater_fluid;

gkyl_dg_updater_fluid*
gkyl_dg_updater_fluid_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_range *conf_range,
  enum gkyl_eqn_type eqn_id, double param, bool use_gpu);

void
gkyl_dg_updater_fluid_advance(gkyl_dg_updater_fluid *fluid,
  enum gkyl_eqn_type eqn_id, const struct gkyl_range *update_rng,
  const struct gkyl_array *u_i, struct gkyl_array *p_ij,
  const struct gkyl_array *aux3,
  const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

void
gkyl_dg_updater_fluid_advance_cu(gkyl_dg_updater_fluid *fluid,
  enum gkyl_eqn_type eqn_id, const struct gkyl_range *update_rng,
  const struct gkyl_array *u_i, struct gkyl_array *p_ij,
  const struct gkyl_array *aux3,
  const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

void gkyl_dg_updater_fluid_release(gkyl_dg_updater_fluid *fluid);

]]

-- Fluid DG solver updater object
local FluidDG = Proto(UpdaterBase)

function FluidDG:init(tbl)
   FluidDG.super.init(self, tbl)

   -- Read data from input file.
   self._onGrid = assert(tbl.onGrid, "Updater.FluidDG: Must provide grid object using 'onGrid'")
   self._confBasis = assert(tbl.confBasis, "Updater.FluidDG: Must specify conf-space basis functions to use using 'confBasis'")
   self._confRange = assert(tbl.confRange, "Updater.FluidDG: Must specify conf-space range using 'confRange'")
   assert(self._confRange:isSubRange()==1, "Eq.Fluid: confRange must be a sub-range")

   self._eqnId = assert(tbl.eqnId, "Updater.FluidDG: Must specify equation ID")

   -- By default, clear output field before incrementing with vol/surf updates.
   self._clearOut = xsys.pickBool(tbl.clearOut, true)

   assert(self._onGrid:ndim() == self._confBasis:ndim(), "Dimensions of basis and grid must match")
   self._ndim = self._onGrid:ndim()

   -- TODO: fix 'param' being passed for GKYL_EQN_EULER and GKYL_EQN_ISO_EULER (for now we just pass 1)
   self._zero = ffi.gc(
                  ffiC.gkyl_dg_updater_fluid_new(self._onGrid._zero, self._confBasis._zero, self._confRange, self._eqnId, 1., GKYL_USE_GPU or 0),
                  ffiC.gkyl_dg_updater_fluid_release
                )

   return self
end

-- advance method
function FluidDG:_advance(tCurr, inFld, outFld)

   local qIn = assert(inFld[1], "FluidDG.advance: Must specify an input field")

   assert(self._eqnId, "FluidDG.advance: Must specify eqn ID")
   if self._eqnId == "GKYL_EQN_EULER_PKPM" or self._eqnId == "GKYL_EQN_EULER" or self._eqnId == "GKYL_EQN_ISO_EULER" then
      local aux_uvar = inFld[2]._zeroDevice
   end
   if self._eqnId == "GKYL_EQN_EULER_PKPM" then
     aux_p_ij = inFld[3]._zeroDevice
   end

   local qRhsOut = assert(outFld[1], "FluidDG.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "FluidDG.advance: Must pass cflRate field in output table")

   local localRange = qRhsOut[1]:localRange()
   print("TODO: FIX THIS CALL (DONT PASS aux_uvar 3 times)")
   ffiC.gkyl_dg_updater_fluid_advance(self._zero, self._eqnId, localRange, aux_uvar, aux_uvar, aux_uvar, qIn._zero, cflRateByCell._zero, qRhsOut._zero)

end

function FluidDG:_advanceOnDevice(tCurr, inFld, outFld)

  local qIn = assert(inFld[1], "FluidDG.advance: Must specify an input field")
  local aux_uvar = inFld[2]._zeroDevice
  if self._eqnId == "GKYL_EQN_EULER_PKPM" then
    aux_p_ij = inFld[3]._zeroDevice
  end

  local qRhsOut = assert(outFld[1], "FluidDG.advance: Must specify an output field")
  local cflRateByCell = assert(outFld[2], "FluidDG.advance: Must pass cflRate field in output table")

  local localRange = qRhsOut:localRange()
  ffiC.gkyl_dg_updater_fluid_advance_cu(self._zero, self._eqnId, localRange, aux_uvar, aux_p_ij, qIn._zero, cflRateByCell._zero, qRhsOut._zero)

end

-- set up pointers to dt and cflRateByCell
function FluidDG:setDtAndCflRate(dt, cflRateByCell)
   FluidDG.super.setDtAndCflRate(self, dt, cflRateByCell)

   if self._onDevice then
      ffiC.setDtAndCflRate(self._onDevice, dt, cflRateByCell._onDevice)
   end
end

return FluidDG
