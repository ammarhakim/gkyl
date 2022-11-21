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
 * Create a new updater that computes the m0_s*delta_s prefactor
 * in the calculation of cross-primitive moments for LBO and BGK
 * collisions. That is:
 *   m0_s*delta_s = m0_s*2*m_r*m0_r*nu_rs/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs)
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
 * @param range Range in which we'll compute m0_s*delta_s.
 * @param out Output array.
 * @return New updater pointer.
 */
void gkyl_prim_cross_m0deltas_advance(gkyl_prim_cross_m0deltas *up, struct gkyl_basis basis,
  double massself, struct gkyl_array* m0self, struct gkyl_array* nuself,
  double massother, struct gkyl_array* m0other, struct gkyl_array* nuother,
  const struct gkyl_range *range, struct gkyl_array* out);

/**
 * Delete updater.
 *
 * @param pob Updater to delete.
 */
void gkyl_prim_cross_m0deltas_release(gkyl_prim_cross_m0deltas* up);
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

   self._onGrid = assert(
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
   self._beta   = tbl.betaGreene
   self._betaP1 = self._beta+1

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   -- Ensure sanity.
   assert(phaseBasis:polyOrder() == self.confBasis:polyOrder(),
          "Updater.CrossPrimMoments: Polynomial orders of phase and conf basis must match.")
   assert((phaseBasis:id() == self.confBasis:id()) or
         ((phaseBasis:id()=="hybrid" or phaseBasis:id()=="gkhybrid") and self.confBasis:id()=="serendipity"),
          "Updater.CrossPrimMoments: Type of phase and conf basis must match.")

   if operator=="VmBGK"or operator=="GkBGK" then
      -- Determine configuration and velocity space dims.
      local cDim = self.confBasis:ndim()
      local vDim = phaseBasis:ndim() - cDim
      self._uDim = vDim
      if operator=="GkLBO" or operator=="GkBGK" then self._uDim = 1 end

      -- Need two Eigen matrices: one to divide by (ms*nusr*m0s+mr*nurs*m0r)
      -- and one to compute cross-primitive moments.
      local numBasisC = self.confBasis:numBasis()
      self._binOpData    = ffiC.new_binOpData_t(numBasisC*2*(self._uDim+1), 0)
      self._binOpDataDiv = ffiC.new_binOpData_t(numBasisC, 0)

      -- Select C kernel to be used.
      local basisID, polyOrder = self.confBasis:id(), self.confBasis:polyOrder()
      self._crossPrimMomentsCalc = PrimMomentsDecl.selectCrossPrimMomentsCalc(operator, basisID, cDim, vDim, polyOrder)

      self._uCrossOtherBuf    = Lin.Vec(self._uDim*self.confBasis:numBasis())
      self._vtSqCrossOtherBuf = Lin.Vec(self.confBasis:numBasis())

      -- To obtain the cell average, multiply the zeroth coefficient by this factor.
      self._cellAvFac = 1.0/math.sqrt(2.0^cDim)
   end

   if self._useGPU then
      if operator=="GkLBO" then
         self._zero = ffi.gc(ffiC.gkyl_prim_lbo_gyrokinetic_cross_calc_cu_dev_new(self._onGrid._zero, self.confBasis._zero, phaseBasis._zero),
                             ffiC.gkyl_prim_lbo_cross_calc_release)
      elseif operator=="VmLBO" then
         self._zero = ffi.gc(ffiC.gkyl_prim_lbo_vlasov_cross_calc_cu_dev_new(self._onGrid._zero, self.confBasis._zero, phaseBasis._zero),
                             ffiC.gkyl_prim_lbo_cross_calc_release)
      else -- MF 2022/11/19: BGK not implemented on GPU.
         self._advanceFunc = self._advanceNoDeviceImpl
      end
   else
      if operator=="GkLBO" then
         self._zero = ffi.gc(ffiC.gkyl_prim_lbo_gyrokinetic_cross_calc_new(self._onGrid._zero, self.confBasis._zero, phaseBasis._zero),
                             ffiC.gkyl_prim_lbo_cross_calc_release)
      elseif operator=="VmLBO" then
         self._zero = ffi.gc(ffiC.gkyl_prim_lbo_vlasov_cross_calc_new(self._onGrid._zero, self.confBasis._zero, phaseBasis._zero),
                             ffiC.gkyl_prim_lbo_cross_calc_release)
      end
   end
   
   if operator=="VmLBO" or operator=="GkLBO" then
      local onRange = self.onGhosts and self.m0deltas:localExtRange() or self.m0deltas:localRange()
      self._zero_m0deltas = ffi.gc(ffiC.gkyl_prim_cross_m0deltas_new(self.confBasis._zero, onRange, self._betaP1, self._useGPU),
                                   ffiC.gkyl_prim_cross_m0deltas_release)
   end
end

-- Advance method.
function CrossPrimMoments:_advance(tCurr, inFld, outFld)

   if self._zero then

      local mSelf, nuSelf = inFld[1], inFld[2]
      local momsSelf      = inFld[3]
      local primMomsSelf  = inFld[4]
      local bCorrsSelf    = inFld[5]

      local mOther, nuOther = inFld[6], inFld[7]
      local momsOther       = inFld[8]
      local primMomsOther   = inFld[9]

      local primMomsCrossSelf = outFld[1]

      -- Compose the pre-factor:
      --   m0_s*delta_s*(1+beta)
      --     = m0_s*(2*m_r*m0_r*nu_rs/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs))*(1+beta)
      ffiC.gkyl_prim_cross_m0deltas_advance(self._zero_m0deltas, self.confBasis._zero,
         mSelf, momsSelf._zero, nuSelf._zero, mOther, momsOther._zero, nuOther._zero,
         momsSelf:localRange(), self.m0deltas._zero)

      -- Compute u and vtsq.
      ffiC.gkyl_prim_lbo_cross_calc_advance(self._zero, self.confBasis._zero, primMomsSelf:localRange(), self.m0deltas._zero,
         mSelf, momsSelf._zero, primMomsSelf._zero, mOther, momsOther._zero, primMomsOther._zero, 
         bCorrsSelf._zero, primMomsCrossSelf._zero)

      return
   end

   local mSelf, nuSelf = inFld[1], inFld[2]
   local momsSelf      = inFld[3]
   local primMomsSelf  = inFld[4]

   local mOther, nuOther = inFld[5], inFld[6]
   local momsOther       = inFld[7]
   local primMomsOther   = inFld[8]

   local m0sdeltas = inFld[9]

   local primMomsCrossSelf = outFld[1]

   local confRange = self.onGhosts and nuSelf:localExtRange() or nuSelf:localRange()

   local confIndexer = nuSelf:genIndexer()

   local momsSelfItr     = momsSelf:get(1)
   local primMomsSelfItr = primMomsSelf:get(1)

   local momsOtherItr     = momsOther:get(1)
   local primMomsOtherItr = primMomsOther:get(1)

   local primMomsCrossSelfItr = primMomsCrossSelf:get(1)

   local nuSelf0, nuSelfItr
   local nuOther0, nuOtherItr
   nuSelfItr  = nuSelf:get(1) 
   nuOtherItr = nuOther:get(1)

   for cIdx in confRange:rowMajorIter() do
      momsSelf:fill(confIndexer(cIdx), momsSelfItr)
      primMomsSelf:fill(confIndexer(cIdx), primMomsSelfItr)
      local uSelfItr    = primMomsSelfItr:data()
      local vtSqSelfItr = primMomsSelfItr:data()+self._uDim*self.confBasis:numBasis()
   
      momsOther:fill(confIndexer(cIdx), momsOtherItr)
      primMomsOther:fill(confIndexer(cIdx), primMomsOtherItr)
      local uOtherItr    = primMomsOtherItr:data()
      local vtSqOtherItr = primMomsOtherItr:data()+self._uDim*self.confBasis:numBasis()
   
      primMomsCrossSelf:fill(confIndexer(cIdx), primMomsCrossSelfItr)
      local uCrossSelfItr     = primMomsCrossSelfItr:data()
      local vtSqCrossSelfItr  = primMomsCrossSelfItr:data()+self._uDim*self.confBasis:numBasis()
      local uCrossOtherItr    = self._uCrossOtherBuf:data()
      local vtSqCrossOtherItr = self._vtSqCrossOtherBuf:data()
   
      nuSelf:fill(confIndexer(cIdx), nuSelfItr)
      nuOther:fill(confIndexer(cIdx), nuOtherItr)
      local nuSelf0  = nuSelfItr[1]*self._cellAvFac 
      local nuOther0 = nuOtherItr[1]*self._cellAvFac 
   
      self._crossPrimMomentsCalc(self._binOpDataDiv, self._betaP1, mSelf, nuSelf0, momsSelfItr:data(), uSelfItr, vtSqSelfItr, mOther, nuOther0, momsOtherItr:data(), uOtherItr, vtSqOtherItr, uCrossSelfItr, vtSqCrossSelfItr, uCrossOtherItr, vtSqCrossOtherItr)
   end

end

function CrossPrimMoments:_advanceOnDevice(tCurr, inFld, outFld)

   local mSelf, nuSelf = inFld[1], inFld[2]
   local momsSelf      = inFld[3]
   local primMomsSelf  = inFld[4]
   local bCorrsSelf    = inFld[5]

   local mOther, nuOther = inFld[6], inFld[7]
   local momsOther       = inFld[8]
   local primMomsOther   = inFld[9]

   local primMomsCrossSelf = outFld[1]

   -- Compose the pre-factor:
   --   m0_s*delta_s*(1+beta)
   --     = m0_s*(2*m_r*m0_r*nu_rs/(m_s*m0_s*nu_sr+m_r*m0_r*nu_rs))*(1+beta)
   ffiC.gkyl_prim_cross_m0deltas_advance(self._zero_m0deltas, self.confBasis._zero,
      mSelf, momsSelf._zeroDevice, nuSelf._zeroDevice, mOther, momsOther._zeroDevice, nuOther._zeroDevice,
      momsSelf:localRange(), self.m0deltas._zeroDevice)

   -- Compute u and vtsq.
   ffiC.gkyl_prim_lbo_cross_calc_advance_cu(self._zero, self.confBasis._zero, primMomsSelf:localRange(), self.m0deltas._zeroDevice,
      mSelf, momsSelf._zeroDevice, primMomsSelf._zeroDevice,
      mOther, momsOther._zeroDevice, primMomsOther._zeroDevice, 
      bCorrsSelf._zeroDevice, primMomsCrossSelf._zeroDevice)

end

return CrossPrimMoments
