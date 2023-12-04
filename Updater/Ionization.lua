-- Gkyl ------------------------------------------------------------------------
--
-- Updater for ADAS ionization calculations.
--
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto       = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local ffi         = require "ffi"
local xsys        = require "xsys"
local DataStruct  = require "DataStruct"

local ffiC = ffi.C
local new, sizeof, typeof, metatype = xsys.from(ffi,
     "new, sizeof, typeof, metatype")

local g0_share_prefix = os.getenv("HOME") .. "/gkylsoft/gkylzero/share"

ffi.cdef [[ 
// Identifiers for different ionization types
enum gkyl_dg_iz_type
{
  GKYL_IZ_H = 0,  // Hydrogen ions
  GKYL_IZ_HE = 1, // Helium ions
  GKYL_IZ_LI = 2, // Lithium ions
  GKYL_IZ_BE = 3, // Beryllium ions
  GKYL_IZ_B = 4,  // Boron ions
  GKYL_IZ_C = 5,  // Carbon ions
  GKYL_IZ_N = 6,  // Nitrogen ions
  GKYL_IZ_O = 7,  // Oxygen ions
  GKYL_IZ_AR = 8, // Argon ions
};

// Identifiers for self species to determine form of collision operator
enum gkyl_dg_iz_self
{
  GKYL_IZ_ELC = 0, // Electron species
  GKYL_IZ_ION = 1, // Resulting ion species (increases charge state)
  GKYL_IZ_DONOR = 2, // Reacting species (donates electron)
};

struct gkyl_dg_iz_inp {
  const struct gkyl_rect_grid* grid; // Grid object needed for fmax
  struct gkyl_basis* cbasis; // Configuration-space basis-functions
  struct gkyl_basis* pbasis; // Phase-space basis-functions
  const struct gkyl_range *conf_rng; // Configuration range
  const struct gkyl_range *phase_rng; // Phase range
  double mass_ion; // Mass of the ion 
  enum gkyl_dg_iz_type type_ion; // Enum for type of ion for ionization (H thru 0)
  int charge_state; // Ion charge state
  enum gkyl_dg_iz_self type_self; // Species type (ion, electron or donor)
  bool all_gk; // To indicate if all 3 interacting species are GK or not
  const char* base; // File path to locate adas data
};

// Object type
typedef struct gkyl_dg_iz gkyl_dg_iz;

/**
 * Create new updater to calculate ionization temperature or reaction rate
 * @param gkyl_dg_iz_inp
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_iz* gkyl_dg_iz_new(struct gkyl_dg_iz_inp *inp, bool use_gpu);

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
void gkyl_dg_iz_release(struct gkyl_dg_iz *up);
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
   self._ionMass = assert(tbl.ionMass, "Updater.Ionization: Must provide ion mass using 'ionMass'")
   self._plasma = assert(tbl.plasma, "Updater.Ionization: Must provide ion element type using 'plasma'")
   self._chargeState = assert(tbl.chargeState, "Updater.Ionization: Must provide charge state of donor species using 'chargeState'")
   self._selfSpecies = assert(tbl.selfSpecies, "Updater.Ionization: Must provide self species type using 'selfSpecies'")
   self._allGk = tbl.donorGk

   local ion_type, self_type
   if self._plasma == "H" then ion_type = 0
   elseif self._plasma == "He" then ion_type = 1
   elseif self._plasma == "Li" then ion_type = 2
   elseif self._plasma == "Be" then ion_type = 3
   elseif self._plasma == "B" then ion_type = 4
   elseif self._plasma == "C" then ion_type = 5
   elseif self._plasma == "N" then ion_type = 6
   elseif self._plasma == "O" then ion_type = 7
   elseif self._plasma == "Ar" then ion_type = 8
   else error("Updater.Ionization: 'ionType' must be one of 'H','He','Li','Be','B','C','N','O','Ar'. Was " .. self._plasma .. " instead") end

   if self._selfSpecies == 'elc' then self_type = 0
   elseif self._selfSpecies == 'ion' then self_type = 1
   elseif self._selfSpecies == 'donor' then self_type = 2
   else error("Updater.Ionization: 'selfSpecies' must be one of 'elc', 'ion', or 'donor'. Was " .. self._plasma .. " instead") end

   local izInp  = ffi.new("struct gkyl_dg_iz_inp")
      izInp.grid = self._onGrid._zero
      izInp.cbasis = self._confBasis._zero
      izInp.pbasis = self._phaseBasis._zero
      izInp.conf_rng = self._confRange
      izInp.phase_rng = self._phaseRange
      izInp.mass_ion = self._ionMass
      izInp.type_ion = ion_type
      izInp.charge_state = self._chargeState
      izInp.type_self = self_type
      izInp.all_gk = self._allGk
      izInp.base = g0_share_prefix
   
   self._zero = ffi.gc(
      ffiC.gkyl_dg_iz_new(izInp, GKYL_USE_GPU or 0),
      ffiC.gkyl_dg_iz_release
   )
   return self
end

function Ionization:_advance(tCurr, inFld, outFld)

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

function Ionization:_advanceOnDevice(tCurr, inFld, outFld)

   local momsElc = assert(inFld[1], "Ionization.advance: Must pass input momsElc")
   local momsDonor = assert(inFld[2], "Ionization.advance: Must pass input momsDonor")
   local bmag = assert(inFld[3], "Ionization.advance: Must pass input bmag")
   local jacobTot = assert(inFld[4], "Ionization.advance: Must pass input jacobTot")
   local b_i = assert(inFld[5], "Ionization.advance: Must pass input b_i")
   local distfSelf = assert(inFld[6], "Ionization.advance: Must pass input distfSelf")
   
   local collIz  = assert(outFld[1], "Ionization.advance: Must specifiy output field collIz")
   local cflRateByCell = assert(outFld[2], "Ionization.advance: Must pass cflRate field in output table")
   
   ffiC.gkyl_dg_iz_coll(self._zero, momsElc._zeroDevice, momsDonor._zeroDevice, bmag._zeroDevice, jacobTot._zeroDevice, b_i._zeroDevice, distfSelf._zeroDevice, collIz._zeroDevice, cflRateByCell._zeroDevice)

end

function Ionization:evalMomTime() return self._tmEvalMom end

return Ionization
