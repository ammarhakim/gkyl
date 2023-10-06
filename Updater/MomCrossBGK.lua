-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute the cross moments for the BGK collision operator.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Proto = require "Proto"
local UpdaterBase = require "Updater.Base"
local xsys = require "xsys"
local ffi = require "ffi"

local ffiC = ffi.C
ffi.cdef[[
// Type for storing preallocated memory
typedef struct gkyl_mom_cross_bgk_gyrokinetic gkyl_mom_cross_bgk_gyrokinetic;

/**
 * Allocate memory for use in cross moments calculation.
 *
 * @ param phase_basis Phase space basis functions
 * @ param conf_basis Configuration space basis functions
 * @ param use_gpu Boolian to determine if on GPU
 */
gkyl_mom_cross_bgk_gyrokinetic* gkyl_mom_cross_bgk_gyrokinetic_new(
  const struct gkyl_basis *phase_basis, const struct gkyl_basis *conf_basis, bool use_gpu);

/**
 * Compute the cross moments with moments of each species.
 *
 * @param up Cross moments updater
 * @param conf_rng Configuration space range
 * @param beta Arbitrary parameter
 * @param m_self Mass of the self species
 * @param moms_self Moments of the self species
 * @param m_other Mass of the other species
 * @param moms_other Moments of the other species
 * @param nu_sr Cross collisionality, self with other
 * @param nu_rs Cross collisionality, other with self
 * @param moms_cross Six output moments
 */
void gkyl_mom_cross_bgk_gyrokinetic_advance(
  gkyl_mom_cross_bgk_gyrokinetic *up,
  const struct gkyl_range *conf_rng, const double beta,
  const double m_self, const struct gkyl_array *moms_self,
  const double m_other, const struct gkyl_array *moms_other,
  const double nu_sr, const double nu_rs, struct gkyl_array *moms_cross);
/*
void gkyl_mom_cross_bgk_gyrokinetic_advance_cu(
  gkyl_mom_cross_bgk_gyrokinetic *up,
  const struct gkyl_range *conf_rng, const double beta,
  const double m_self, const struct gkyl_array *moms_self,
  const double m_other, const struct gkyl_array *moms_other,
  const double nu_sr, const double nu_rs, struct gkyl_array *moms_cross);
*/
/**
 * Release memory needed in the cross moments calculation.
 *
 * @param up Memory to release.
 */
void gkyl_mom_cross_bgk_gyrokinetic_release(gkyl_mom_cross_bgk_gyrokinetic *up);
]]

-- Inherit the base Updater from UpdaterBase updater object.
local MomCrossBGK = Proto(UpdaterBase)

function MomCrossBGK:init(tbl)
   MomCrossBGK.super.init(self, tbl)  -- set up base object
   local phaseBasis = assert(tbl.phaseBasis, "Updater.MomCrossBGK: Must provide phase space basis object 'phaseBasis'")
   local confBasis = assert(tbl.confBasis, "Updater.MomCrossBGK: Must provide configuration space basis object 'confBasis'")
   
   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)
  
   self.beta   = tbl.beta
   self.mass_s = tbl.mass_s
   self.mass_r = tbl.mass_r
   self.nu_sr  = tbl.nu_sr
   self.nu_rs  = tbl.nu_rs
 
   self._zero = ffi.gc(ffiC.gkyl_mom_cross_bgk_gyrokinetic_new(phaseBasis._zero, confBasis._zero, self._useGPU), ffiC.gkyl_mom_cross_bgk_gyrokinetic_release)
end

function MomCrossBGK:_advance(tCurr, inFld, outFld)
   local moms_s = assert(inFld[1], "MomCrossBGK.advance: Must specify moments of the self species in 'inFld[1]'")
   local moms_r = assert(inFld[2], "MomCrossBGK.advance: Must specify moments of the other species in 'inFld[2]'")
   local moms_cross =  assert(outFld[1], "MomCrossBGK.advance: Must specify an output field in 'outFld[1]'")

   local confRange = moms_cross:localRange() 

   ffiC.gkyl_mom_cross_bgk_gyrokinetic_advance(self._zero, confRange, self.beta, self.mass_s, moms_s._zero, self.mass_r, moms_r._zero, self.nu_sr, self.nu_rs, moms_cross._zero)
end
--[[
function MomCrossBGK:_advanceOnDevice(tCurr, inFld, outFld)
   local moms_s = assert(inFld[1], "MomCrossBGK.advance: Must specify moments of the self species in 'inFld[1]'")
   local moms_r = assert(inFld[2], "MomCrossBGK.advance: Must specify moments of the other species in 'inFld[2]'")
   local moms_cross =  assert(outFld[1], "MomCrossBGK.advance: Must specify an output field in 'outFld[1]'")

   local confRange = moms_cross:localRange() 

   ffiC.gkyl_mom_cross_bgk_gyrokinetic_advance(self._zero, confRange, self.beta, self.mass_s, moms_s._zeroDevice, self.mass_r, moms_r._zeroDevice, self.nu_sr, self.nu_rs, moms_cross._zeroDevice)
end
--]]
return MomCrossBGK

