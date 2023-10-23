-- Gkyl ------------------------------------------------------------------------
--
-- Updater for ADAS ionization calculations.
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

local g0_share_prefix = os.getenv("HOME") .. "/gkylsoft/gkylzero/share"

ffi.cdef [[ 
// Identifiers for different ionization types
enum gkyl_dg_iz_type
{
  GKYL_IZ_H,  // Hydrogen ions
  GKYL_IZ_HE, // Helium ions
  GKYL_IZ_LI, // Lithium ions
  GKYL_IZ_BE, // Beryllium ions
  GKYL_IZ_B,  // Boron ions
  GKYL_IZ_C,  // Carbon ions
  GKYL_IZ_N,  // Nitrogen ions
  GKYL_IZ_O,  // Oxygen ions
};

// Identifiers for self species to determine form of collision operator
enum gkyl_dg_iz_self
{
  GKYL_IZ_ELC, // Electron species
  GKYL_IZ_ION, // Resulting ion species (increases charge state)
  GKYL_IZ_DONOR, // Reacting species (donates electron)
};

// Object type
typedef struct gkyl_dg_iz gkyl_dg_iz;

/**
 * Create new updater to calculate ionization temperature or reaction rate
 * @param grid Grid object needed for fmax
 * @param cbasis Configuration-space basis-functions
 * @param pbasis Phase-space basis-functions (GK)
 * @param conf_rng Configuration range
 * @param phase_rng Phase range
 * @param elem_charge Elementary charge value
 * @param mass_elc Mass of the electron 
 * @param mass_elc Mass of the ion
 * @param type_ion Enum for type of ion for ionization (support H, He, Li)
 * @param charge_state Int for ion charge state
 * @param type_self Enum for species type (electron, ion or neutral)
 * @param all_gk Boolean to indicate if all 3 species are gyrokinetic
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_iz* gkyl_dg_iz_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, double elem_charge,
  double mass_elc, double mass_ion, enum gkyl_dg_iz_type type_ion, int charge_state, enum gkyl_dg_iz_self type_self,
  bool all_gk, const char *base, bool use_gpu); 

/**
 * Create new ionization updater type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_iz* gkyl_dg_iz_cu_dev_new(struct gkyl_rect_grid* grid, struct gkyl_basis* cbasis, struct gkyl_basis* pbasis,
  const struct gkyl_range *conf_rng, const struct gkyl_range *phase_rng, double elem_charge,
  double mass_elc, double mass_ion, enum gkyl_dg_iz_type type_ion, int charge_state, enum gkyl_dg_iz_self type_self, 
  bool all_gk, const char *base); 

/**
 * Compute ionization collision term for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param iz Ionization object.
 * @param moms_elc Input electron moments
 * @param moms_donor Input neutral moments
 * @param bmag Magnetic field used for GK fmax 
 * @param jacob_tot Total Jacobian used for GK fmax
 * @param bhat_vec Unit bmag vector in Cartesian (X,Y,Z) components
 * @param distf_self Species self distribution function
 * @param coll_iz Output reaction rate coefficient
 */

void gkyl_dg_iz_coll(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_donor,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *distf_self, struct gkyl_array *coll_iz, struct gkyl_array *cflrate);

void gkyl_dg_iz_coll_cu(const struct gkyl_dg_iz *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_donor,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *distf_self, struct gkyl_array *coll_iz, struct gkyl_array *cflrate);

/**
 * Delete updater.
 *
 * @param up Updater to delete.
 */
void gkyl_dg_iz_release(gkyl_dg_iz up);
]]
-- Voronov Collisions updater object.
local Ionization = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function Ionization:init(tbl)
   Ionization.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid, "Updater.Ionization: Must provide grid object using 'onGrid'")
   self._confBasis = assert(tbl.confBasis, "Updater.Ionization: Must provide configuration space basis object using 'confBasis'")
   self._phaseBasis = assert(tbl.phaseBasis, "Updater.Ionization: Must provide phase space basis object using 'phaseBasis'")
   self._confRange = assert(tbl.confRange, "Updater.Ionization: Must provide configuration space range object using 'confRange'")
   self._phaseRange = assert(tbl.phaseRange, "Updater.Ionization: Must provide phase space range object using 'phaseRange'")
   self._elemCharge = assert(tbl.elemCharge, "Updater.Ionization: Must provide elementary charge using 'elemCharge'")
   self._elcMass = assert(tbl.elcMass, "Updater.Ionization: Must provide electron mass using 'elcMass'")
   self._ionMass = assert(tbl.ionMass, "Updater.Ionization: Must provide ion mass using 'ionMass'")
   self._plasma = assert(tbl.plasma, "Updater.Ionization: Must provide ion element type using 'plasma'")
   self._chargeState = assert(tbl.chargeState, "Updater.Ionization: Must provide charge state of donor species using 'chargeState'")
   self._selfSpecies = assert(tbl.selfSpecies, "Updater.Ionization: Must provide self species type using 'selfSpecies'")
   
   if self._plasma == "H" then
      self._ionType = "GKYL_IZ_H"
   elseif self._plasma == "He" then
      self._ionType = "GKYL_IZ_HE"
   elseif self._plasma == "Li" then
      self._ionType = "GKYL_IZ_LI"
   elseif self._plasma == "Be" then
      self._ionType = "GKYL_IZ_BE"
   elseif self._plasma == "B" then
      self._ionType = "GKYL_IZ_B"
   elseif self._plasma == "C" then
      self._ionType = "GKYL_IZ_C"
   elseif self._plasma == "N" then
      self._ionType = "GKYL_IZ_N"
   elseif self._plasma == "O" then
      self._ionType = "GKYL_IZ_O"
   else error("Updater.Ionization: 'ionType' must be one of 'H','He','Li','Be','B','C','N','O'. Was " .. self._plasma .. " instead") end

   if self._selfSpecies == 'elc' then
      self._selfType = "GKYL_IZ_ELC"
   elseif self._selfSpecies == 'ion' then
      self._selfType = "GKYL_IZ_ION"
   elseif self._selfSpecies == 'donor' then
      self._selfType = "GKYL_IZ_DONOR"
   else error("Updater.Ionization: 'selfSpecies' must be one of 'elc', 'ion', or 'donor'. Was " .. self._plasma .. " instead") end

   local allGK = true --fix this
      
   self._zero = ffi.gc(
                  ffiC.gkyl_dg_iz_new(self._onGrid._zero, self._confBasis._zero, self._phaseBasis._zero, self._confRange, self._phaseRange, self._elemCharge, self._elcMass, self._ionMass, self._ionType, self._chargeState, self._selfType, allGK, g0_share_prefix, GKYL_USE_GPU or 0),
                  ffiC.gkyl_dg_iz_release
                )
   return self
end

function Ionization:advance(tCurr, inFld, outFld)

   local momsElc = assert(inFld[1], "Ionization.advance: Must pass input momsElc")
   local momsDonor = assert(inFld[2], "Ionization.advance: Must pass input momsDonor")
   local bmag = assert(inFld[3], "Ionization.advance: Must pass input bmag")
   local jacobTot = assert(inFld[4], "Ionization.advance: Must pass input jacobTot")
   local b_i = assert(inFld[5], "Ionization.advance: Must pass input b_i")
   local distfSelf = assert(inFld[6], "Ionization.advance: Must pass input distfSelf")
   
   local collIz  = assert(outFld[1], "Ionization.advance: Must specifiy output field collIz")
   local cflRateByCell = assert(outFld[2], "Ionization.advance: Must pass cflRate field in output table")
   
   ffiC.gkyl_dg_iz_coll(self._zero, momsElc._zero, momsDonor._zero, bmag._zero, jacobTot._zero, b_i._zero, distfSelf._zero, collIz._zero, cflRateByCell._zero)

end

function Ionization:evalMomTime() return self._tmEvalMom end

return Ionization
