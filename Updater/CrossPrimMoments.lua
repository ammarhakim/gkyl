-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the cross-primitive moments for the cross species 
-- collisions given the masses, collisionalities, coupling moments, primitive
-- moments and boundary corrections (and maybe Greene's beta) of two species.
--
-- For electron-ion plasmas and LBO, these are u_ei, u_ie, vtSq_ei=T_ei/m_e, 
-- and vtSq_ie=T_ie/m_i. Here the subscript ei means that it corresponds to
-- the effect on electrons due to collisions with ions, and vice-versa for ie.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local UpdaterBase     = require "Updater.Base"
local Lin             = require "Lib.Linalg"
local Proto           = require "Lib.Proto"
local PrimMomentsDecl = require "Updater.primMomentsCalcData.PrimMomentsModDecl"
local xsys            = require "xsys"
local ffi             = require "ffi"

local ffiC = ffi.C
ffi.cdef [[
// Object type
typedef struct gkyl_prim_lbo_cross_calc gkyl_prim_lbo_cross_calc;

/**
 * Compute cross-primitive moments of distribution function. The conf_rng
 * MUST be a sub-range of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Primitive moment calculator updater to run
 * @param cbasis_rng Config-space basis functions
 * @param conf_rng Config-space range
 * @param greene Greene's factor
 * @param self_m Mass of the species
 * @param self_moms Moments of distribution function (Zeroth, First, and Second)
 * @param self_prim_moms Drift velocity & thermal speed squared of this species
 * @param other_m Mass of the colliding species
 * @param other_moms Moments of distribution function (Zeroth, First, and Second)
 * @param other_prim_moms Drift velocity & thermal speed squared of the colliding species
 * @param boundary_corrections Momentum and Energy boundary corrections
 * @param prim_moms_out Output drift velocity and thermal speed squared
 */
void gkyl_prim_lbo_cross_calc_advance(gkyl_prim_lbo_cross_calc* calc,
  struct gkyl_basis cbasis, const struct gkyl_range *conf_rng,
  const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_moms, const struct gkyl_array *self_prim_moms,
  double other_m, const struct gkyl_array *other_moms, const struct gkyl_array *other_prim_moms,
  const struct gkyl_array *boundary_corrections,
  struct gkyl_array *prim_moms_out);

void gkyl_prim_lbo_cross_calc_advance_cu(gkyl_prim_lbo_cross_calc* calc,
  struct gkyl_basis cbasis, const struct gkyl_range *conf_rng,
  const struct gkyl_array *greene,
  double self_m, const struct gkyl_array *self_moms, const struct gkyl_array *self_prim_moms,
  double other_m, const struct gkyl_array *other_moms, const struct gkyl_array *other_prim_moms,
  const struct gkyl_array *boundary_corrections,
  struct gkyl_array *prim_moms_out);

/**
 * Delete pointer to primitive moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_prim_lbo_cross_calc_release(gkyl_prim_lbo_cross_calc* calc);

// "derived" class constructors
gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_vlasov_cross_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_gyrokinetic_cross_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_vlasov_cross_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_cross_calc*
gkyl_prim_lbo_gyrokinetic_cross_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);


// Object type
typedef struct gkyl_prim_cross_m0deltas gkyl_prim_cross_m0deltas;

/**
 * Create a new updater that computes the m0_s*delta_s*(beta+1) prefactor
 * in the calculation of cross-primitive moments for LBO and BGK
 * collisions. That is:
 *   m0_s*delta_s*(beta+1) = m0_s*2*m_r*m0_r*nu_rs*(beta+1)/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs)
 *
 * @param basis Basis object (configuration space).
 * @param range Range in which we'll compute m0_s*delta_s.
 * @param betap1 beta+1 parameter in Greene's formulism.
 * @param use_gpu boolean indicating whether to use the GPU.
 * @return New updater pointer.
 */
gkyl_prim_cross_m0deltas* gkyl_prim_cross_m0deltas_new(
  const struct gkyl_basis *basis, const struct gkyl_range *range,
  double betap1, bool use_gpu);

/**
 * Compute m0_s*delta_s.
 *
 * @param up Struct defining this updater..
 * @param basis Basis object (configuration space).
 * @param massself Mass of this species, m_s.
 * @param m0self Number density of this species, m0_s.
 * @param nuself Cross collision frequency of this species, nu_sr.
 * @param massother Mass of the other species, m_r.
 * @param m0other Number density of the other species, m0_r.
 * @param nuother Cross collision frequency of the other species, nu_rs.
 * @param prem0s Prefactor, M_0s for LBO, 1 for BGK.
 * @param range Range in which we'll compute m0_s*delta_s.
 * @param out Output array.
 * @return New updater pointer.
 */
void gkyl_prim_cross_m0deltas_advance(gkyl_prim_cross_m0deltas *up, struct gkyl_basis basis,
  double massself, const struct gkyl_array* m0self, const struct gkyl_array* nuself,
  double massother, const struct gkyl_array* m0other, const struct gkyl_array* nuother,
  const struct gkyl_array* prem0s, const struct gkyl_range *range, struct gkyl_array* out);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_prim_cross_m0deltas_release(gkyl_prim_cross_m0deltas* up);

/**
 * Calculate the cross primitive moments, u_{sr,i} and v_{tsr}^2,
 * for cross-species collisions with the BGK operator.
 *
 *   u_sr = u_si-0.5*delta_s*(beta+1)*(u_si - u_ri)
 *
 *   v_tsr^2 = v_ts^2-(u_sri-u_ri).(u_sri-u_si)/vdim_phys
 *     -(delta_s*(beta+1)/(m_s+m_r))*
 *      (m_s*v_ts^2-m_r*v_tr^2+((m_s-m_r)/(2*vdim_phys))*(u_si-u_ri)^2)
 *
 * @param basis Basis object (configuration space).
 * @param vdim_phys Physical number of velocity dimensions represented.
 * @param m0sdeltas Prefactor m_0s*delta_s*(beta+1) (and m_0s should be 1 here).
 * @param massself Mass of this species.
 * @param primsself Primitive moments of this species, u_s and v_ts^2.
 * @param massother Mass of the other species.
 * @param primsother Primitive moments of the other species, u_r and v_tr^2.
 * @param range Range in which we'll compute m0_s*delta_s.
 * @return crossprims Cross primitive moments, u_sri and v_tsr^2.
 */
void gkyl_prim_bgk_cross_calc_advance(struct gkyl_basis basis,
  int vdim_phys, const struct gkyl_array* m0sdeltas,
  double massself, const struct gkyl_array* primsself,
  double massother, const struct gkyl_array* primsother,
  const struct gkyl_range *range, struct gkyl_array* crossprims);

void gkyl_prim_bgk_cross_calc_advance_cu(struct gkyl_basis basis,
  int vdim_phys, const struct gkyl_array* m0sdeltas,
  double massself, const struct gkyl_array* primsself,
  double massother, const struct gkyl_array* primsother,
  const struct gkyl_range *range, struct gkyl_array* crossprims);
]]

-- Function to check if operator option is correct.
local function isOperatorGood(nm)
   if nm == "VmLBO" or nm == "GkLBO" or nm == "VmBGK" or nm == "GkBGK" then
      return true
   end
   return false
end

-- Primitive moments updater object.
local CrossPrimMoments = Proto(UpdaterBase)

function CrossPrimMoments:init(tbl)
   CrossPrimMoments.super.init(self, tbl) -- Setup base object.

   local onGrid = assert(
      tbl.onGrid, "Updater.CrossPrimMoments: Must provide grid object using 'onGrid'.")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.CrossPrimMoments: Must provide the phase basis object using 'phaseBasis'.")

   self.confBasis = assert(
      tbl.confBasis, "Updater.CrossPrimMoments: Must provide the configuration basis object using 'confBasis'.")

   local operator = assert(
      tbl.operator, "Updater.CrossPrimMoments: Must specify the collision operator (VmLBO, GkLBO, VmBGK or GkBGK) using 'operator'.")
   assert(isOperatorGood(operator), string.format("CrossPrimMoments: Operator option must be 'VmLBO', 'GkLBO', 'VmBGK' or 'GkBGK'. Requested %s instead.", operator))

   self.m0deltas = assert(tbl.m0s_deltas, "Updater.CrossPrimMoments: Must provide the film to hold m0_s*delta_s in 'm0s_deltas'.")

   -- Free-parameter in John Greene's equations.
   local beta   = tbl.betaGreene
   self._betaP1 = beta+1

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   -- Ensure sanity.
   assert(phaseBasis:polyOrder() == self.confBasis:polyOrder(),
          "Updater.CrossPrimMoments: Polynomial orders of phase and conf basis must match.")
   assert((phaseBasis:id() == self.confBasis:id()) or
         ((phaseBasis:id()=="hybrid" or phaseBasis:id()=="gkhybrid") and self.confBasis:id()=="serendipity"),
          "Updater.CrossPrimMoments: Type of phase and conf basis must match.")

   local isLBO, isGK = true, false
   if operator == "VmBGK" or operator=="GkBGK" then isLBO = false end
   if operator == "GkLBO" or operator=="GkBGK"  then isGK = true end

   -- Determine configuration and velocity space dims.
   local cDim = self.confBasis:ndim()
   local vDim = phaseBasis:ndim() - cDim
   self._vDimPhys = isGK and 2*vDim-1 or vDim

   if self._useGPU then
      if operator=="GkLBO" then
         self._zero = ffi.gc(ffiC.gkyl_prim_lbo_gyrokinetic_cross_calc_cu_dev_new(onGrid._zero, self.confBasis._zero, phaseBasis._zero),
                             ffiC.gkyl_prim_lbo_cross_calc_release)
      elseif operator=="VmLBO" then
         self._zero = ffi.gc(ffiC.gkyl_prim_lbo_vlasov_cross_calc_cu_dev_new(onGrid._zero, self.confBasis._zero, phaseBasis._zero),
                             ffiC.gkyl_prim_lbo_cross_calc_release)
      end
   else
      if operator=="GkLBO" then
         self._zero = ffi.gc(ffiC.gkyl_prim_lbo_gyrokinetic_cross_calc_new(onGrid._zero, self.confBasis._zero, phaseBasis._zero),
                             ffiC.gkyl_prim_lbo_cross_calc_release)
      elseif operator=="VmLBO" then
         self._zero = ffi.gc(ffiC.gkyl_prim_lbo_vlasov_cross_calc_new(onGrid._zero, self.confBasis._zero, phaseBasis._zero),
                             ffiC.gkyl_prim_lbo_cross_calc_release)
      end
   end
   
   local onRange = self.onGhosts and self.m0deltas:localExtRange() or self.m0deltas:localRange()
   self._zero_m0deltas = ffi.gc(ffiC.gkyl_prim_cross_m0deltas_new(self.confBasis._zero, onRange, self._betaP1, self._useGPU),
                                ffiC.gkyl_prim_cross_m0deltas_release)

   if isLBO then
      self.advanceFunc = function(tCurr, inFld, outFld)
         CrossPrimMoments["_advanceLBO"](self, tCurr, inFld, outFld) end
      self.advanceOnDeviceFunc = function(tCurr, inFld, outFld)
         CrossPrimMoments["_advanceOnDeviceLBO"](self, tCurr, inFld, outFld) end
   else
      self.advanceFunc = function(tCurr, inFld, outFld)
         CrossPrimMoments["_advanceBGK"](self, tCurr, inFld, outFld) end
      self.advanceOnDeviceFunc = function(tCurr, inFld, outFld)
         CrossPrimMoments["_advanceOnDeviceBGK"](self, tCurr, inFld, outFld) end
   end
end

-- Advance method.
function CrossPrimMoments:_advanceLBO(tCurr, inFld, outFld)
   local mSelf, nuSelf = inFld[1], inFld[2]
   local momsSelf      = inFld[3]
   local primMomsSelf  = inFld[4]
   local bCorrsSelf    = inFld[5]

   local mOther, nuOther = inFld[6], inFld[7]
   local momsOther       = inFld[8]
   local primMomsOther   = inFld[9]

   local m0sPreFac = inFld[10]

   local primMomsCrossSelf = outFld[1]


   -- Compose the pre-factor:
   --   m0_s*delta_s*(1+beta)
   --     = m0_s*(2*m_r*m0_r*nu_rs/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs))*(1+beta)
   ffiC.gkyl_prim_cross_m0deltas_advance(self._zero_m0deltas, self.confBasis._zero,
      mSelf, momsSelf._zero, nuSelf._zero, mOther, momsOther._zero, nuOther._zero,
      m0sPreFac._zero, momsSelf:localRange(), self.m0deltas._zero)

   -- Compute u and vtsq.
   ffiC.gkyl_prim_lbo_cross_calc_advance(self._zero, self.confBasis._zero, primMomsSelf:localRange(), self.m0deltas._zero,
      mSelf, momsSelf._zero, primMomsSelf._zero, mOther, momsOther._zero, primMomsOther._zero, 
      bCorrsSelf._zero, primMomsCrossSelf._zero)
end

function CrossPrimMoments:_advanceBGK(tCurr, inFld, outFld)
   local mSelf, nuSelf = inFld[1], inFld[2]
   local momsSelf      = inFld[3]
   local primMomsSelf  = inFld[4]

   local mOther, nuOther = inFld[5], inFld[6]
   local momsOther       = inFld[7]
   local primMomsOther   = inFld[8]

   local m0sPreFac = inFld[9]

   local primMomsCrossSelf = outFld[1]

   -- Compose the pre-factor:
   --   m0_s*delta_s*(1+beta)
   --     = m0_s*(2*m_r*m0_r*nu_rs/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs))*(1+beta)
   ffiC.gkyl_prim_cross_m0deltas_advance(self._zero_m0deltas, self.confBasis._zero,
      mSelf, momsSelf._zero, nuSelf._zero, mOther, momsOther._zero, nuOther._zero,
      m0sPreFac._zero, momsSelf:localRange(), self.m0deltas._zero)

   -- Compute u and vtsq.
   ffiC.gkyl_prim_bgk_cross_calc_advance(self.confBasis._zero, self._vDimPhys,
      self.m0deltas._zero, mSelf, primMomsSelf._zero, mOther, primMomsOther._zero,
      primMomsSelf:localRange(), primMomsCrossSelf._zero)
end

function CrossPrimMoments:_advance(tCurr, inFld, outFld)
   self.advanceFunc(tCurr, inFld, outFld)
end

function CrossPrimMoments:_advanceOnDeviceLBO(tCurr, inFld, outFld)
   local mSelf, nuSelf = inFld[1], inFld[2]
   local momsSelf      = inFld[3]
   local primMomsSelf  = inFld[4]
   local bCorrsSelf    = inFld[5]

   local mOther, nuOther = inFld[6], inFld[7]
   local momsOther       = inFld[8]
   local primMomsOther   = inFld[9]

   local m0sPreFac = inFld[10]

   local primMomsCrossSelf = outFld[1]

   -- Compose the pre-factor:
   --   m0_s*delta_s*(1+beta)
   --     = m0_s*(2*m_r*m0_r*nu_rs/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs))*(1+beta)
   ffiC.gkyl_prim_cross_m0deltas_advance(self._zero_m0deltas, self.confBasis._zero,
      mSelf, momsSelf._zeroDevice, nuSelf._zeroDevice, mOther, momsOther._zeroDevice, nuOther._zeroDevice,
      m0sPreFac._zeroDevice, momsSelf:localRange(), self.m0deltas._zeroDevice)

   -- Compute u and vtsq.
   ffiC.gkyl_prim_lbo_cross_calc_advance_cu(self._zero, self.confBasis._zero, primMomsSelf:localRange(), self.m0deltas._zeroDevice,
      mSelf, momsSelf._zeroDevice, primMomsSelf._zeroDevice,
      mOther, momsOther._zeroDevice, primMomsOther._zeroDevice, 
      bCorrsSelf._zeroDevice, primMomsCrossSelf._zeroDevice)
end

function CrossPrimMoments:_advanceOnDeviceBGK(tCurr, inFld, outFld)
   local mSelf, nuSelf = inFld[1], inFld[2]
   local momsSelf      = inFld[3]
   local primMomsSelf  = inFld[4]

   local mOther, nuOther = inFld[5], inFld[6]
   local momsOther       = inFld[7]
   local primMomsOther   = inFld[8]

   local m0sPreFac = inFld[9]

   local primMomsCrossSelf = outFld[1]

   -- Compose the pre-factor:
   --   m0_s*delta_s*(1+beta)
   --     = m0_s*(2*m_r*m0_r*nu_rs/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs))*(1+beta)
   ffiC.gkyl_prim_cross_m0deltas_advance(self._zero_m0deltas, self.confBasis._zero,
      mSelf, momsSelf._zeroDevice, nuSelf._zeroDevice, mOther, momsOther._zeroDevice, nuOther._zeroDevice,
      m0sPreFac._zeroDevice, momsSelf:localRange(), self.m0deltas._zeroDevice)

   -- Compute u and vtsq.
   ffiC.gkyl_prim_bgk_cross_calc_advance_cu(self.confBasis._zero, self._vDimPhys,
      self.m0deltas._zeroDevice, mSelf, primMomsSelf._zeroDevice, mOther, primMomsOther._zeroDevice,
      primMomsSelf:localRange(), primMomsCrossSelf._zeroDevice)
end

function CrossPrimMoments:_advanceOnDevice(tCurr, inFld, outFld)
   self.advanceOnDeviceFunc(tCurr, inFld, outFld)
end

return CrossPrimMoments
