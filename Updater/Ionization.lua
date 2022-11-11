-- Gkyl ------------------------------------------------------------------------
--
-- Updater for Voronov ionization calculations:
--     temperature (vtSqIz) for electrons or Voronov reaction rate (coefIz)
--
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
local Time        = require "Lib.Time"
local UpdaterBase = require "Updater.Base"
local ffi         = require "ffi"
local xsys        = require "xsys"

local ffiC = ffi.C
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

ffi.cdef [[ 
// Identifiers for different ionization types
enum gkyl_dg_iz_type
{
  GKYL_H, // Hydrogen plasma
  GKYL_AR, // Argon plasma
};

// Object type
typedef struct gkyl_dg_iz gkyl_dg_iz;

/**
 * Create new updater to calculate ionization temperature or reaction rate
 * @param cbasis Configuration-space basis-functions
 * @param elem_charge Elementary charge value
 * @param mass_elc Mass of the electron value
 * @param type_ion Enum for type of ion for ionization (support H^+ and Ar^+)
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_iz* gkyl_dg_iz_new(const struct gkyl_basis* cbasis, 
  double elem_charge, double mass_elc, enum gkyl_dg_iz_type type_ion, 
  bool use_gpu);

/**
 * Compute ionization temperature for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param iz Ionization object.
 * @param update_rng Update range (Configuration space)
 * @param vth_sq_elc Input electron vth^2 = T/m
 * @param vth_sq_iz Output ionization vth^2 = T_iz/m
 */

void gkyl_dg_iz_temp(const struct gkyl_dg_iz *iz,
  const struct gkyl_range *update_rng, const struct gkyl_array *vth_sq_elc,
  struct gkyl_array *vth_sq_iz);

/**
 * Compute reaction rate coefficient for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param iz Ionization object.
 * @param update_rng Update range (Configuration space)
 * @param phase_rng Phase-space range (for indexing cflrate for stable timestep)
 * @param n_neut Input neutral density 
 * @param vth_sq_neut Input neutral vth^2 = T/m
 * @param vth_sq_elc Input electron vth^2 = T/m
 * @param cflrate CFL scalar rate (frequency) array (units of 1/[T])
 * @param coef_iz Output reaction rate coefficient
 */

void gkyl_dg_iz_react_rate(const struct gkyl_dg_iz *viz,
  const struct gkyl_range *update_rng, const struct gkyl_range *phase_rng, 
  const struct gkyl_array *n_neut, const struct gkyl_array *vth_sq_neut, const struct gkyl_array *vth_sq_elc,
  struct gkyl_array *cflrate, struct gkyl_array *coef_iz);

/**
 * Delete updater.
 *
 * @param iz Updater to delete.
 */
void gkyl_dg_iz_release(gkyl_dg_iz *iz);
]]
-- Voronov Collisions updater object.
local Ionization = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function Ionization:init(tbl)
   Ionization.super.init(self, tbl) -- setup base object
   
   self._confBasis = assert(tbl.confBasis, "Updater.Ionization: Must provide configuration space basis object using 'confBasis'")
   self._elemCharge = assert(tbl.elemCharge, "Updater.Ionization: Must provide elementary charge using 'elemCharge'")
   self._elcMass = assert(tbl.elcMass, "Updater.Ionization: Must provide electron mass using 'elcMass'")
   self._plasma = assert(tbl.plasma, "Updater.Ionization: Must provide ion type using 'plasma'")

   if self._plasma == "H" then
      self._ionType = "GKYL_H"
   else if self._plasma == "Ar" then
      self._ionType = "GKYL_AR"
   end

   self._zero = ffi.gc(
                  ffiC.gkyl_dg_iz_new(self._confBasis._zero, self._elemCharge, self._elcMass, self._ionType, GKYL_USE_GPU or 0),
                  ffiC.gkyl_dg_iz_release
                )

   self._tmEvalMom = 0.0

   return self
end

function Ionization:ionizationTemp(tCurr, inFld, outFld)
   local tmEvalMomStart = Time.clock()

   local elcVtSq = inFld[1]
   local vtSqIz  = outFld[1]

   local localRange = vtSqIz:localRange()
   ffiC.gkyl_dg_iz_temp(self._zero, localRange, elcVtSq._zero, vtSqIz._zero)

   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function Ionization:reactRateCoef(tCurr, inFld, outFld) 
   local tmEvalMomStart = Time.clock()

   local neutM0   = inFld[1]
   local neutVtSq = inFld[2]
   local elcVtSq  = inFld[3]
   local coefIz   = outFld[1]
   local cflRateByCell = outFld[2]

   local localRange = coefIz:localRange()
   -- Reaction rate sets a cflrate for a stable timestep
   -- Note: This is a phase space cflrate, localRange for cflRateByCell is phaseRange
   local cflRange = cflRateByCell:localRange()

   ffiC.gkyl_dg_iz_react_rate(self._zero, localRange, cflRange, 
      neutM0._zero, neutVtSq._zero, elcVtSq._zero, cflRateByCell._zero, coefIz._zero)

   self._tmEvalMom = self._tmEvalMom + Time.clock() - tmEvalMomStart
end

function Ionization:evalMomTime() return self._tmEvalMom end

return Ionization
