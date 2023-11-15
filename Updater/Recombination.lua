-- Gkyl ------------------------------------------------------------------------
--
-- Updater for ADAS recombination calculations.
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
// Identifiers for different ionrecombation types
enum gkyl_dg_recomb_type
{
  GKYL_RECOMB_H = 0,  // Hydrogen ions
  GKYL_RECOMB_HE = 1, // Helium ions
  GKYL_RECOMB_LI = 2, // Lithium ions
  GKYL_RECOMB_BE = 3, // Beryllium ions
  GKYL_RECOMB_B = 4,  // Boron ions
  GKYL_RECOMB_C = 5,  // Carbon ions
  GKYL_RECOMB_N = 6,  // Nitrogen ions
  GKYL_RECOMB_O = 7,  // Oxygen ions
};

// Identifiers for self species to determine form of collision operator
enum gkyl_dg_recomb_self
{
  GKYL_RECOMB_ELC = 0, // Electron species
  GKYL_RECOMB_ION = 1, // Reacting ion species (increases charge state)
  GKYL_RECOMB_RECVR = 2, // Resulting species (receives electron)
};

struct gkyl_dg_recomb_inp {
  const struct gkyl_rect_grid* grid; // Grid object needed for fmax
  struct gkyl_basis* cbasis; // Configuration-space basis-functions
  struct gkyl_basis* pbasis; // Phase-space basis-functions
  const struct gkyl_range *conf_rng; // Configuration range
  const struct gkyl_range *phase_rng; // Phase range
  double mass_self; // Mass of the species
  enum gkyl_dg_recomb_type type_ion; // Enum for type of ion for ionization (H thru 0)
  int charge_state; // Ion charge state
  enum gkyl_dg_recomb_self type_self; // Species type (ion, electron or receiver)
  bool all_gk; // To indicate if all 3 interacting species are GK or not
  const char* base; // File path to locate adas data
};


// Object type
typedef struct gkyl_dg_recomb gkyl_dg_recomb;

/**
 * Create new updater to calculate recomb temperature or reaction rate
 * @param gkyl_dg_recomb_inp
 * @param use_gpu Boolean for whether struct is on host or device
 */
struct gkyl_dg_recomb* gkyl_dg_recomb_new(struct gkyl_dg_recomb_inp inp, bool use_gpu); 

/**
 * Create new ionrecombation updater type object on NV-GPU: 
 * see new() method above for documentation.
 */
struct gkyl_dg_recomb* gkyl_dg_recomb_cu_dev_new(struct gkyl_dg_recomb_inp inp); 

/**
 * Compute recombation collision term for use in neutral reactions. 
 * The update_rng MUST be a sub-range of the
 * range on which the array is defined.  That is, it must be either
 * the same range as the array range, or one created using the
 * gkyl_sub_range_init method.
 *
 * @param recomb Recombation object.
 * @param moms_elc Input electron moments
 * @param moms_ion Input ion moments
 * @param bmag Magnetic field used for GK fmax 
 * @param jacob_tot Total Jacobian used for GK fmax
 * @param b_i Unit bmag vector in Cartesian (X,Y,Z) components
 * @param distf_self Species self distribution function
 * @param coll_recomb Output reaction rate coefficient
 */
void gkyl_dg_recomb_coll(const struct gkyl_dg_recomb *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_ion,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i, 
  const struct gkyl_array *f_self, struct gkyl_array *coll_recomb, struct gkyl_array *cflrate);

void gkyl_dg_recomb_coll_cu(const struct gkyl_dg_recomb *up,
  const struct gkyl_array *moms_elc, const struct gkyl_array *moms_ion,
  const struct gkyl_array *bmag, const struct gkyl_array *jacob_tot, const struct gkyl_array *b_i,
  const struct gkyl_array *f_self, struct gkyl_array *coll_recomb, struct gkyl_array *cflrate);
/**
 * Delete updater.
 *
 * @param recomb Updater to delete.
 */
void gkyl_dg_recomb_release(gkyl_dg_recomb *recomb);
]]
-- Recomination collisions updater object.
local Recombination = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function Recombination:init(tbl)
   Recombination.super.init(self, tbl) -- setup base object

   self._onGrid = assert(tbl.onGrid, "Updater.Recombination: Must provide grid object using 'onGrid'")
   self._confBasis = assert(tbl.confBasis, "Updater.Recombination: Must provide configuration space basis object using 'confBasis'")
   self._phaseBasis = assert(tbl.phaseBasis, "Updater.Recombination: Must provide phase space basis object using 'phaseBasis'")
   self._confRange = assert(tbl.confRange, "Updater.Recombination: Must provide configuration space range object using 'confRange'")
   self._phaseRange = assert(tbl.phaseRange, "Updater.Recombination: Must provide phase space range object using 'phaseRange'")
   self._selfMass = assert(tbl.selfMass, "Updater.Recombination: Must provide ion mass using 'selfMass'")
   self._plasma = assert(tbl.plasma, "Updater.Recombination: Must provide ion element type using 'plasma'")
   self._chargeState = assert(tbl.chargeState, "Updater.Recombination: Must provide charge state of donor species using 'chargeState'")
   self._selfSpecies = assert(tbl.selfSpecies, "Updater.Recombination: Must provide self species type using 'selfSpecies'")
   self._allGk = tbl.recvrGk

   local ion_type, self_type
   if self._plasma == "H" then ion_type = 0
   elseif self._plasma == "He" then ion_type = 1
   elseif self._plasma == "Li" then ion_type = 2
   elseif self._plasma == "Be" then ion_type = 3
   elseif self._plasma == "B" then ion_type = 4
   elseif self._plasma == "C" then ion_type = 5
   elseif self._plasma == "N" then ion_type = 6
   elseif self._plasma == "O" then ion_type = 7
   else error("Updater.Recombination: 'ionType' must be one of 'H','He','Li','Be','B','C','N','O'. Was " .. self._plasma .. " instead") end

   if self._selfSpecies == 'elc' then self_type = 0
   elseif self._selfSpecies == 'ion' then self_type = 1
   elseif self._selfSpecies == 'recvr' then self_type = 2
   else error("Updater.Recombination: 'selfSpecies' must be one of 'elc', 'ion', or 'recvr'. Was " .. self._plasma .. " instead") end

   local recInp  = ffi.new("struct gkyl_dg_recomb_inp")
      recInp.grid = self._onGrid._zero
      recInp.cbasis = self._confBasis._zero
      recInp.pbasis = self._phaseBasis._zero
      recInp.conf_rng = self._confRange
      recInp.phase_rng = self._phaseRange
      recInp.mass_self = self._selfMass
      recInp.type_ion = ion_type
      recInp.charge_state = self._chargeState
      recInp.type_self = self_type
      recInp.all_gk = self._allGk
      recInp.base = g0_share_prefix
   
   self._zero = ffi.gc(
      ffiC.gkyl_dg_recomb_new(recInp, GKYL_USE_GPU or 0),
      ffiC.gkyl_dg_recomb_release
   )
   return self
end

function Recombination:advance(tCurr, inFld, outFld)

   local momsElc = assert(inFld[1], "Recombination.advance: Must pass input momsElc")
   local momsIon = assert(inFld[2], "Recombination.advance: Must pass input momsIon")
   local bmag = assert(inFld[3], "Recombination.advance: Must pass input bmag")
   local jacobTot = assert(inFld[4], "Recombination.advance: Must pass input jacobTot")
   local b_i = assert(inFld[5], "Recombination.advance: Must pass input b_i") -- check geo fields...
   local distfSelf = assert(inFld[6], "Recombination.advance: Must pass input distfSelf")
   
   local collRecomb  = assert(outFld[1], "Recombination.advance: Must specifiy output field collRecom")
   local cflRateByCell = assert(outFld[2], "Recombination.advance: Must pass cflRate field in output table")
   
   ffiC.gkyl_dg_recomb_coll(self._zero, momsElc._zero, momsIon._zero, bmag._zero, jacobTot._zero, b_i._zero, distfSelf._zero, collRecomb._zero, cflRateByCell._zero)

end

function Recombination:evalMomTime() return self._tmEvalMom end

return Recombination
