-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the Spitzer collision frequency.
-- There are two options.
-- a) Scale a normalized nu by n_r/(v_ts^2+v_tr^2)^(3/2).
-- c) Build nu in SI units from scrath.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase = require "Updater.Base"
local Proto       = require "Lib.Proto"
local xsys        = require "xsys"
local ffi         = require "ffi"

local ffiC = ffi.C
ffi.cdef [[
// Object type
typedef struct gkyl_spitzer_coll_freq gkyl_spitzer_coll_freq;

/**
 * Create new updater to either compute the Spitzer collision frequency from
 * scratch based on local parameters, or scale a normalized collision frequency
 * by the local n_r/(v_ts^2+v_tr^2)^(3/2).
 *
 * @param basis Basis object (configuration space).
 * @param num_quad Number of quadrature nodes.
 * @param use_gpu Boolean indicating whether to use the GPU.
 * @param nufrac Fraction to multiply the collision frequency by.
 * @param eps0 Permittivity of vacuum.
 * @param hbar Planck's constant divided by 2*pi.
 * @return New updater pointer.
 */
gkyl_spitzer_coll_freq* gkyl_spitzer_coll_freq_new(const struct gkyl_basis *basis,
  int num_quad, double nufrac, double eps0, double hbar, bool use_gpu);

/**
 * Scale the normalized collision frequency, normNu, by
 * n_r/(v_ts^2+v_tr^2)^(3/2) and project it on to the basis.
 *
 * @param up Spizer collision frequency updater object.
 * @param range Config-space rang.e
 * @param vtSqSelf Thermal speed squared of this species.
 * @param m0Other Thermal speed squared of the other species.
 * @param vtSqOther Thermal speed squared of the other species.
 * @param normNu Normalized collision frequency to scale.
 * @param nuOut Output collision frequency.
 */
void gkyl_spitzer_coll_freq_advance_normnu(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *vtSqSelf,
  const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther,
  double normNu, struct gkyl_array *nuOut);

/**
 * Compute the Spitzer collision frequency from scratch. Coulomb Logarithm
 * is computed using cell averaged values.
 *
 * @param up Spizer collision frequency updater object.
 * @param range Config-space range.
 * @param bmag Magnetic field magnitude.
 * @param qSelf Charge of this species.
 * @param mSelf Mass of this species.
 * @param m0Self Thermal speed squared of the other species.
 * @param vtSqSelf Thermal speed squared of this species.
 * @param qOther Charge of this species.
 * @param mOther Mass of this species.
 * @param m0Other Thermal speed squared of the other species.
 * @param vtSqOther Thermal speed squared of the other species.
 * @param nuOut Output collision frequency.
 */
void gkyl_spitzer_coll_freq_advance(const gkyl_spitzer_coll_freq *up,
  const struct gkyl_range *range, const struct gkyl_array *bmag,
  double qSelf, double mSelf, const struct gkyl_array *m0Self, const struct gkyl_array *vtSqSelf,
  double qOther, double mOther, const struct gkyl_array *m0Other, const struct gkyl_array *vtSqOther,
  struct gkyl_array *nuOut);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_spitzer_coll_freq_release(gkyl_spitzer_coll_freq* up);
]]

-- Updater object.
local SpitzerCollisionality = Proto(UpdaterBase)

function SpitzerCollisionality:init(tbl)
   SpitzerCollisionality.super.init(self, tbl) -- Setup base object.

   local confBasis = assert(tbl.confBasis,
      "Updater.SpitzerCollisionality: Must provide the configuration basis object using 'confBasis'.")

   local isNormNu = tbl.willInputNormNu

   local epsilon0 = assert(tbl.epsilon0,
      "Updater.SpitzerCollisionality: Must specify vacuum permittivity ('epsilon0') to build Spitzer collisionality.")
   local hBar = assert(tbl.hBar,
      "Updater.SpitzerCollisionality: Must specify Planck's constant h divided by 2pi ('hBar') to build Spitzer collisionality.")
   self._nuFrac = tbl.nuFrac and tbl.nuFrac or 1.0

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   local numQuad = confBasis:polyOrder()+1
   self._zero = ffi.gc(ffiC.gkyl_spitzer_coll_freq_new(confBasis._zero, numQuad, self._nuFrac, epsilon0, hBar, self._useGPU),
                       ffiC.gkyl_spitzer_coll_freq_release)

   if isNormNu then
      self.advanceFunc = function(tCurr, inFld, outFld)
         SpitzerCollisionality["_advance_normNu"](self, tCurr, inFld, outFld)
      end
      self.advanceOnDeviceFunc = function(tCurr, inFld, outFld)
         SpitzerCollisionality["_advanceOnDevice_normNu"](self, tCurr, inFld, outFld)
      end
   else
      self.advanceFunc = function(tCurr, inFld, outFld)
         SpitzerCollisionality["_advance_build"](self, tCurr, inFld, outFld)
      end
      self.advanceOnDeviceFunc = function(tCurr, inFld, outFld)
         SpitzerCollisionality["_advanceOnDevice_build"](self, tCurr, inFld, outFld)
      end
   end

end

-- Advance methods.
function SpitzerCollisionality:_advance_normNu(tCurr, inFld, outFld)
   local chargeSelf, massSelf   = inFld[1], inFld[2]
   local m0Self, vtSqSelf       = inFld[3], inFld[4]
   local chargeOther, massOther = inFld[5], inFld[6]
   local m0Other, vtSqOther     = inFld[7], inFld[8]
   local normNu = inFld[9]

   local nuOut = outFld[1]

   local range = self.onGhosts and nuOut:localExtRange() or nuOut:localRange()

   ffiC.gkyl_spitzer_coll_freq_advance_normnu(self._zero, range,
      vtSqSelf._zero, m0Other._zero, vtSqOther._zero, normNu*self._nuFrac, nuOut._zero)
end

function SpitzerCollisionality:_advanceOnDevice_normNu(tCurr, inFld, outFld)
   local chargeSelf, massSelf   = inFld[1], inFld[2]
   local m0Self, vtSqSelf       = inFld[3], inFld[4]
   local chargeOther, massOther = inFld[5], inFld[6]
   local m0Other, vtSqOther     = inFld[7], inFld[8]
   local normNu = inFld[9]

   local nuOut = outFld[1]

   local range = self.onGhosts and nuOut:localExtRange() or nuOut:localRange()

   ffiC.gkyl_spitzer_coll_freq_advance_normnu(self._zero, range,
      vtSqSelf._zeroDevice, m0Other._zeroDevice, vtSqOther._zeroDevice, normNu*self._nuFrac, nuOut._zeroDevice)
end

function SpitzerCollisionality:_advance_build(tCurr, inFld, outFld)
   local chargeSelf, massSelf   = inFld[1], inFld[2]
   local m0Self, vtSqSelf       = inFld[3], inFld[4]
   local chargeOther, massOther = inFld[5], inFld[6]
   local m0Other, vtSqOther     = inFld[7], inFld[8]
   local Bmag = inFld[10]

   local nuOut = outFld[1]

   local range = self.onGhosts and nuOut:localExtRange() or nuOut:localRange()

   ffiC.gkyl_spitzer_coll_freq_advance(self._zero, range, Bmag._zero, 
      chargeSelf, massSelf, m0Self._zero, vtSqSelf._zero,
      chargeOther, massOther, m0Other._zero, vtSqOther._zero, nuOut._zero)
end

function SpitzerCollisionality:_advanceOnDevice_build(tCurr, inFld, outFld)
   local chargeSelf, massSelf   = inFld[1], inFld[2]
   local m0Self, vtSqSelf       = inFld[3], inFld[4]
   local chargeOther, massOther = inFld[5], inFld[6]
   local m0Other, vtSqOther     = inFld[7], inFld[8]
   local Bmag = inFld[10]

   local nuOut = outFld[1]

   local range = self.onGhosts and nuOut:localExtRange() or nuOut:localRange()

   ffiC.gkyl_spitzer_coll_freq_advance(self._zero, range, Bmag._zeroDevice, 
      chargeSelf, massSelf, m0Self._zeroDevice, vtSqSelf._zeroDevice,
      chargeOther, massOther, m0Other._zeroDevice, vtSqOther._zeroDevice, nuOut._zeroDevice)
end

function SpitzerCollisionality:_advance(tCurr, inFld, outFld)
   self.advanceFunc(tCurr, inFld, outFld)
end

function SpitzerCollisionality:_advanceOnDevice(tCurr, inFld, outFld)
   self.advanceOnDeviceFunc(tCurr, inFld, outFld)
end

return SpitzerCollisionality
