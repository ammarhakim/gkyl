-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the increment due to BGK collisions.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Proto       = require "Proto"
local UpdaterBase = require "Updater.Base"
local xsys        = require "xsys"
local ffi         = require "ffi"

local ffiC = ffi.C
ffi.cdef [[
// Object type
typedef struct gkyl_bgk_collisions gkyl_bgk_collisions;

/**
 * Create new updater to compute the increment due to a BGK
 * collision operator
 *   nu*f_M - nu*f)
 * where nu*f_M=sum_r nu_sr*f_Msr and nu=sum_r nu_sr, in
 * order to support multispecies collisions. The quantities
 * nu*f_M and nu must be computed elsewhere.
 *
 * @param cbasis Basis object (configuration space).
 * @param pbasis Basis object (phase space).
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_bgk_collisions* gkyl_bgk_collisions_new(const struct gkyl_basis *cbasis,
  const struct gkyl_basis *pbasis, bool use_gpu);

/**
 * Advance BGK operator (compute the BGK contribution to df/dt).
 *
 * @param up Spizer collision frequency updater object.
 * @param crange Config-space range.
 * @param prange Phase-space range.
 * @param nu Sum of collision frequencies.
 * @param nufM Sum of collision frequencies times their respective Maxwellian.
 * @param fin Input distribution function.
 * @param out BGK contribution to df/dt.
 * @param cflfreq Output CFL frequency.
 */
void gkyl_bgk_collisions_advance(const gkyl_bgk_collisions *up,
  const struct gkyl_range *crange, const struct gkyl_range *prange,
  const struct gkyl_array *nu, const struct gkyl_array *nufM, const struct gkyl_array *fin,
  struct gkyl_array *out, struct gkyl_array *cflfreq);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_bgk_collisions_release(gkyl_bgk_collisions* up);
]]

-- BGK Collisions updater object.
local BGKcollisions = Proto(UpdaterBase)

function BGKcollisions:init(tbl)
   BGKcollisions.super.init(self, tbl) -- Setup base object.

   self._phaseBasis = assert(tbl.phaseBasis,
      "Updater.BGKcollisions: Must provide phase space basis object using 'phaseBasis'")
   self._confBasis  = assert(tbl.confBasis,
      "Updater.BGKcollisions: Must provide configuration space basis object using 'confBasis'")

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   self._zero = ffi.gc(ffiC.gkyl_bgk_collisions_new(self._confBasis._zero, self._phaseBasis._zero, self._useGPU),
                       ffiC.gkyl_bgk_collisions_release)
end

function BGKcollisions:_advance(tCurr, inFld, outFld)
   local sumNu         = assert(inFld[1],
      "BGKcollisions.advance: Must specify the sum of collisionalities as input[1].")
   local sumNufMaxwell = assert(inFld[2],
      "BGKcollisions.advance: Must specify sum(nu*Maxwellian) field as input[2].")
   local fIn           = assert(inFld[3],
      "BGKcollisions.advance: Must specify an input distribution function field as input[3].")

   local fRhsOut = assert(outFld[1], "BGKcollisions.advance: Must specify an output field in outFld[1].")
   local cflFreq = assert(outFld[2], "BGKcollisions.advance: Must specify the CFL rate in outFld[2].")

   local confRange  = self.onGhosts and sumNu:localExtRange() or sumNu:localRange()
   local phaseRange = self.onGhosts and fRhsOut:localExtRange() or fRhsOut:localRange()

   ffiC.gkyl_bgk_collisions_advance(self._zero, confRange, phaseRange, sumNu._zero,
      sumNufMaxwell._zero, fIn._zero, fRhsOut._zero, cflFreq._zero)
end

function BGKcollisions:_advanceOnDevice(tCurr, inFld, outFld)
   local sumNu         = assert(inFld[1],
      "BGKcollisions.advance: Must specify the sum of collisionalities as input[1].")
   local sumNufMaxwell = assert(inFld[2],
      "BGKcollisions.advance: Must specify sum(nu*Maxwellian) field as input[2].")
   local fIn           = assert(inFld[3],
      "BGKcollisions.advance: Must specify an input distribution function field as input[3].")

   local fRhsOut = assert(outFld[1], "BGKcollisions.advance: Must specify an output field in outFld[1].")
   local cflFreq = assert(outFld[2], "BGKcollisions.advance: Must specify the CFL rate in outFld[2].")

   local confRange  = self.onGhosts and sumNu:localExtRange() or sumNu:localRange()
   local phaseRange = self.onGhosts and fRhsOut:localExtRange() or fRhsOut:localRange()

   ffiC.gkyl_bgk_collisions_advance(self._zero, confRange, phaseRange, sumNu._zeroDevice,
      sumNufMaxwell._zeroDevice, fIn._zeroDevice, fRhsOut._zeroDevice, cflFreq._zeroDevice)
end

return BGKcollisions
