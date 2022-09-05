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
  GKYL_FIELD_SR_E_B, // Maxwell (E, B) with special relativity
  GKYL_FIELD_PHI, // Poisson (only phi)
  GKYL_FIELD_PHI_A, // Poisson with static B = curl(A) (phi, A)
  GKYL_FIELD_NULL, // no field is present
  GKYL_FIELD_SR_NULL // no field is present, special relativistic Vlasov
};

typedef struct gkyl_dg_updater_vlasov gkyl_dg_updater_vlasov;

gkyl_dg_updater_vlasov*
gkyl_dg_updater_vlasov_new(const struct gkyl_rect_grid *grid,
  const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis,
  const struct gkyl_range *conf_range, const struct gkyl_range *vel_range,
  enum gkyl_field_id field_id, bool use_gpu);

void
gkyl_dg_updater_vlasov_advance(gkyl_dg_updater_vlasov *vlasov,
  enum gkyl_field_id field_id, const struct gkyl_range *update_rng,
  const struct gkyl_array *aux1, const struct gkyl_array *aux2,
  const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

void
gkyl_dg_updater_vlasov_advance_cu(gkyl_dg_updater_vlasov *vlasov,
  enum gkyl_field_id field_id, const struct gkyl_range *update_rng,
  const struct gkyl_array *aux1, const struct gkyl_array *aux2,
  const struct gkyl_array* fIn,
  struct gkyl_array* cflrate, struct gkyl_array* rhs);

void gkyl_dg_updater_vlasov_release(gkyl_dg_updater_vlasov *vlasov);

]]

-- Vlasov DG solver updater object
local VlasovDG = Proto(UpdaterBase)

function VlasovDG:init(tbl)
   VlasovDG.super.init(self, tbl)

   -- Read data from input file.
   self._onGrid = assert(tbl.onGrid, "Updater.VlasovDG: Must provide grid object using 'onGrid'")
   self._phaseBasis = assert(tbl.phaseBasis, "Updater.VlasovDG: Must specify phase-space basis functions to use using 'phaseBasis'")
   self._confBasis = assert(tbl.confBasis, "Updater.VlasovDG: Must specify conf-space basis functions to use using 'confBasis'")

   self._confRange = assert(tbl.confRange, "Updater.VlasovDG: Must specify conf-space range using 'confRange'")
   assert(self._confRange:isSubRange()==1, "Eq.Vlasov: confRange must be a sub-range") 

   -- Check if we have an electric and magnetic field.
   local hasElcField    = xsys.pickBool(tbl.hasElectricField, false)
   local hasMagField    = xsys.pickBool(tbl.hasMagneticField, false)
   local hasExtForce    = xsys.pickBool(tbl.hasExtForce, false)
   self._plasmaMagField = xsys.pickBool(tbl.plasmaMagField, false)

   if hasElcField and self._plasmaMagField then 
      self._fieldId = "GKYL_FIELD_E_B"
   elseif hasElcField then
      self._fieldId = "GKYL_FIELD_PHI"
   else
      self._fieldId = "GKYL_FIELD_NULL"
   end

   self._zero = ffi.gc(
                  ffiC.gkyl_dg_updater_vlasov_new(self._onGrid._zero, self._confBasis._zero, self._phaseBasis._zero, self._confRange, nil, self._fieldId, GKYL_USE_GPU or 0),
                  ffiC.gkyl_dg_updater_vlasov_release
                )

   return self
end

-- advance method
function VlasovDG:_advance(tCurr, inFld, outFld)

   local qIn = assert(inFld[1], "VlasovDG.advance: Must specify an input field")
   local aux1, aux2
   if self._fieldId == "GKYL_FIELD_PHI" then
      aux1 = inFld[2]._zero
      aux2 = nil
   elseif self._fieldId == "GKYL_FIELD_E_B" then
      aux1 = inFld[2]._zero
      aux2 = inFld[3]._zero
   end
   local qRhsOut = assert(outFld[1], "VlasovDG.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "VlasovDG.advance: Must pass cflRate field in output table")

   local localRange = qRhsOut:localRange()
   ffiC.gkyl_dg_updater_vlasov_advance(self._zero, self._fieldId, localRange, aux1, aux2, qIn._zero, cflRateByCell._zero, qRhsOut._zero)

end

function VlasovDG:_advanceOnDevice(tCurr, inFld, outFld)

   local qIn = assert(inFld[1], "VlasovDG.advance: Must specify an input field")
   local aux1, aux2
   if self._fieldId == "GKYL_FIELD_PHI" then
      aux1 = inFld[2]._zeroDevice
      aux2 = inFld[2]._zeroDevice
   elseif self._fieldId == "GKYL_FIELD_E_B" then
      aux1 = inFld[2]._zeroDevice
      aux2 = inFld[3]._zeroDevice
   end
   local qRhsOut = assert(outFld[1], "VlasovDG.advance: Must specify an output field")
   local cflRateByCell = assert(outFld[2], "VlasovDG.advance: Must pass cflRate field in output table")

   local localRange = qRhsOut:localRange()
   ffiC.gkyl_dg_updater_vlasov_advance_cu(self._zero, self._fieldId, localRange, aux1, aux2, qIn._zeroDevice, cflRateByCell._zeroDevice, qRhsOut._zeroDevice)

end

-- set up pointers to dt and cflRateByCell
function VlasovDG:setDtAndCflRate(dt, cflRateByCell)
   VlasovDG.super.setDtAndCflRate(self, dt, cflRateByCell)

   if self._onDevice then
      ffiC.setDtAndCflRate(self._onDevice, dt, cflRateByCell._onDevice)
   end
end

return VlasovDG
