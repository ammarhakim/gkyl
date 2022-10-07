-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the primitive moments, u and vtSq=sqrt(T/m), given
-- the distribution function and its moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase     = require "Updater.Base"
local Lin             = require "Lib.Linalg"
local Proto           = require "Lib.Proto"
local PrimMomentsDecl = require "Updater.primMomentsCalcData.PrimMomentsModDecl"
local xsys            = require "xsys"
local ffi             = require "ffi"

local ffiC = ffi.C
require "Lib.ZeroUtil"

ffi.cdef [[
// Object type
typedef struct gkyl_prim_lbo_calc gkyl_prim_lbo_calc;

/**
 * Compute primitive moments of distribution function. The phase_rng and conf_rng
 * MUST be a sub-ranges of the range on which the distribution
 * function and the moments are defined. These ranges must be
 * on_dev-consistently constructed.
 *
 * @param calc Primitive moment calculator updater to run
 * @param cbasis_rng Config-space basis functions
 * @param conf_rng Config-space range
 * @param moms Moments of distribution function (Zeroth, First, and Second)
 * @param boundary_corrections Momentum and Energy boundary corrections
 * @param uout Output drift velocity primitive moment array
 * @param vtSqout Output thermal velocity primitive moment array
 */
void gkyl_prim_lbo_calc_advance(gkyl_prim_lbo_calc* calc, struct gkyl_basis cbasis,
  struct gkyl_range *conf_rng,
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections,
  struct gkyl_array *uout, struct gkyl_array *vtSqout);

void gkyl_prim_lbo_calc_advance_cu(gkyl_prim_lbo_calc* calc, struct gkyl_basis cbasis,
  struct gkyl_range *conf_rng,
  const struct gkyl_array *moms, const struct gkyl_array *boundary_corrections,
  struct gkyl_array* uout, struct gkyl_array* vtSqout);

/**
 * Delete pointer to primitive moment calculator updater.
 *
 * @param calc Updater to delete.
 */
void gkyl_prim_lbo_calc_release(gkyl_prim_lbo_calc* calc);

// "derived" class constructors
gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_calc*
gkyl_prim_lbo_gyrokinetic_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_with_fluid_calc_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, const struct gkyl_range *conf_rng);

gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_calc*
gkyl_prim_lbo_gyrokinetic_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis);

gkyl_prim_lbo_calc*
gkyl_prim_lbo_vlasov_with_fluid_calc_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis *cbasis, const struct gkyl_basis *pbasis, const struct gkyl_range *conf_rng);


// Object type
typedef struct gkyl_mom_calc_bcorr gkyl_mom_calc_bcorr;

/**
 * Compute boundary correction moments.
 *
 * @param bcorr Boundary correction updater object
 * @param phase_rng Phase space range on which to compute.
 * @param conf_rng Configuration space range on which to compute.
 * @param fIn Input to updater
 * @param out Output
 */
void gkyl_mom_calc_bcorr_advance(gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fIn, struct gkyl_array *out);

void gkyl_mom_calc_bcorr_advance_cu(gkyl_mom_calc_bcorr *bcorr,
  const struct gkyl_range *phase_rng, const struct gkyl_range *conf_rng,
  const struct gkyl_array *fIn, struct gkyl_array *out);

/**
 * Delete updater.
 *
 * @param bcorr Updater to delete.
 */
void gkyl_mom_calc_bcorr_release(gkyl_mom_calc_bcorr* up);

// "derived" class constructors
gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_vlasov_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double* vBoundary);

gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_gyrokinetic_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double* vBoundary, double mass);

gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_vlasov_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double* vBoundary);

gkyl_mom_calc_bcorr*
gkyl_mom_calc_bcorr_lbo_gyrokinetic_cu_dev_new(const struct gkyl_rect_grid *grid, const struct gkyl_basis* cbasis, const struct gkyl_basis* pbasis, const double* vBoundary, double mass);
]]

-- function to check if operator option is correct
local function isOperatorGood(nm)
   if nm == "VmLBO" or nm == "GkLBO" then
      return true
   end
   return false
end

-- Moments updater object.
local SelfPrimMoments = Proto(UpdaterBase)

function SelfPrimMoments:init(tbl)
   SelfPrimMoments.super.init(self, tbl) -- setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.SelfPrimMoments: Must provide grid object using 'onGrid'.")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.SelfPrimMoments: Must provide the phase basis object using 'phaseBasis'.")

   self.confBasis = assert(
      tbl.confBasis, "Updater.SelfPrimMoments: Must provide the configuration basis object using 'confBasis'.")
   local confBasis = self.confBasis

   local operator = assert(
      tbl.operator, "Updater.SelfPrimMoments: Must specify the collision operator (VmLBO or GkLBO for now) using 'operator'.")

   -- Dimension of phase space.
   self._pDim = phaseBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID   = confBasis:id()
   self._polyOrder = confBasis:polyOrder()

   -- Ensure sanity.
   assert(phaseBasis:polyOrder() == self._polyOrder,
          "Polynomial orders of phase and conf basis must match.")
   assert((phaseBasis:id() == self._basisID) or
          ((phaseBasis:id()=="hybrid" or phaseBasis:id()=="gkhybrid") and self._basisID=="serendipity"),
          "Type of phase and conf basis must match.")
   -- Determine configuration and velocity space dims.
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Number of basis functions. Used to compute number of vector components.
   self._numBasisP = phaseBasis:numBasis()
   self._numBasisC = confBasis:numBasis()

   local uDim       = self._vDim  -- Dimensionality of flow velocity vector.
   self._kinSpecies = "Vm"        -- Vlasov-Maxwell species.
   self._isGkLBO    = false
   if isOperatorGood(operator) then
     if operator == "GkLBO" then
        self._isGkLBO    = true
        self._kinSpecies = "Gk"    -- Gyrokinetic species.
        uDim             = 1       -- A (vpar,mu) simulation has 3 physical velocity dimensions.
     end
   else
      assert(false, string.format(
                "SelfPrimMoments: Operator option must be 'VmLBO' or 'GkLBO'. Requested %s instead.", operator))
   end

   self._SelfPrimMomentsCalc = PrimMomentsDecl.selectSelfPrimMomentsCalc(self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)
   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._binOpData = ffiC.new_binOpData_t(self._numBasisC*(uDim+1), 0) 

   local vbounds = assert(tbl.vbounds, "SelfPrimMoments: must pass vbounds")
   local mass
   if self._isGkLBO then mass = assert(tbl.mass, "SelfPrimMoments: must pass species mass for GkLBO") end
   if GKYL_USE_GPU then
      if self._isGkLBO then
         self._zero_bcorr_calc = ffi.gc(ffiC.gkyl_mom_calc_bcorr_lbo_gyrokinetic_cu_dev_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero, vbounds, mass),
                                        ffiC.gkyl_mom_calc_bcorr_release)
         self._zero_prim_calc  = ffi.gc(ffiC.gkyl_prim_lbo_gyrokinetic_calc_cu_dev_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero),
                                        ffiC.gkyl_prim_lbo_calc_release)
      else
         self._zero_bcorr_calc = ffi.gc(ffiC.gkyl_mom_calc_bcorr_lbo_vlasov_cu_dev_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero, vbounds),
                                        ffiC.gkyl_mom_calc_bcorr_release)
         self._zero_prim_calc  = ffi.gc(ffiC.gkyl_prim_lbo_vlasov_calc_cu_dev_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero),
                                        ffiC.gkyl_prim_lbo_calc_release)
      end
   else
      if self._isGkLBO then
         self._zero_bcorr_calc = ffi.gc(ffiC.gkyl_mom_calc_bcorr_lbo_gyrokinetic_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero, vbounds, mass),
                                        ffiC.gkyl_mom_calc_bcorr_release)
         self._zero_prim_calc  = ffi.gc(ffiC.gkyl_prim_lbo_gyrokinetic_calc_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero),
                                        ffiC.gkyl_prim_lbo_calc_release)
      else
         self._zero_bcorr_calc = ffi.gc(ffiC.gkyl_mom_calc_bcorr_lbo_vlasov_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero, vbounds),
                                        ffiC.gkyl_mom_calc_bcorr_release)
         self._zero_prim_calc  = ffi.gc(ffiC.gkyl_prim_lbo_vlasov_calc_new(self._onGrid._zero, confBasis._zero, phaseBasis._zero),
                                        ffiC.gkyl_prim_lbo_calc_release)
      end
   end
end

-- Advance method.
function SelfPrimMoments:_advance(tCurr, inFld, outFld)

   if self._zero_prim_calc then

      local moments, fIn = inFld[1], inFld[2]
      local boundaryCorrections, u, vtsq = outFld[1], outFld[2], outFld[3]

      -- Compute boundary corrections.
      ffiC.gkyl_mom_calc_bcorr_advance(self._zero_bcorr_calc, fIn:localRange(), u:localRange(), fIn._zero, boundaryCorrections._zero)

      -- Compute u and vtsq.
      ffiC.gkyl_prim_lbo_calc_advance(self._zero_prim_calc, self.confBasis._zero, u:localRange(), moments._zero, boundaryCorrections._zero, u._zero, vtsq._zero)

      return
   end

   local grid = self._onGrid

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local uOut           = outFld[1]
   local vtSqOut        = outFld[2]

   local uOutItr        = uOut:get(1)
   local vtSqOutItr     = vtSqOut:get(1)

   -- Moments used for all polyOrders.
   local m0, m1 = inFld[1], inFld[2]
   local confIndexer = m0:genIndexer()
   local m0Itr       = m0:get(1)
   local m1Itr       = m1:get(1)

   -- Boundary corrections.
   local cMomB, cEnergyB = inFld[4], inFld[5]
   local cMomBItr        = cMomB:get(1)
   local cEnergyBItr     = cEnergyB:get(1)

   local m2, m2Itr
   local m0Star, m0StarItr
   local m1Star, m1StarItr
   local m2Star, m2StarItr
   if self._polyOrder > 1 then
      m2    = inFld[3]
      m2Itr = m2:get(1)
   else
      m0Star, m1Star, m2Star = inFld[6], inFld[7], inFld[8]
      m0StarItr = m0Star:get(1)
      m1StarItr = m1Star:get(1)
      m2StarItr = m2Star:get(1)
   end

   local confRange = self.onGhosts and m0:localExtRange() or m0:localRange()

   -- polyOrder=1 and >1 each use separate velocity grid loops to
   -- avoid evaluating (if polyOrder==1) at each cell.
   if self._polyOrder > 1 then

      for cIdx in confRange:rowMajorIter() do
         grid:setIndex(cIdx)

         m0:fill(confIndexer(cIdx), m0Itr)
         m1:fill(confIndexer(cIdx), m1Itr)
         m2:fill(confIndexer(cIdx), m2Itr)
         cMomB:fill(confIndexer(cIdx), cMomBItr)
         cEnergyB:fill(confIndexer(cIdx), cEnergyBItr)

         uOut:fill(confIndexer(cIdx), uOutItr)
         vtSqOut:fill(confIndexer(cIdx), vtSqOutItr)

         self._SelfPrimMomentsCalc(self._binOpData, m0Itr:data(), m1Itr:data(), m2Itr:data(), cMomBItr:data(), cEnergyBItr:data(), uOutItr:data(), vtSqOutItr:data())
      end

   else

      for cIdx in confRange:rowMajorIter() do
         grid:setIndex(cIdx)

         m0:fill(confIndexer(cIdx), m0Itr)
         m1:fill(confIndexer(cIdx), m1Itr)
         cMomB:fill(confIndexer(cIdx), cMomBItr)
         cEnergyB:fill(confIndexer(cIdx), cEnergyBItr)
         m0Star:fill(confIndexer(cIdx), m0StarItr)
         m1Star:fill(confIndexer(cIdx), m1StarItr)
         m2Star:fill(confIndexer(cIdx), m2StarItr)

         uOut:fill(confIndexer(cIdx), uOutItr)
         vtSqOut:fill(confIndexer(cIdx), vtSqOutItr)

         self._SelfPrimMomentsCalc(self._binOpData, m0Itr:data(), m1Itr:data(), m0StarItr:data(), m1StarItr:data(), m2StarItr:data(), cMomBItr:data(), cEnergyBItr:data(), uOutItr:data(), vtSqOutItr:data())
      end

   end
end

function SelfPrimMoments:_advanceOnDevice(tCurr, inFld, outFld)

   local moments, fIn = inFld[1], inFld[2]
   local boundaryCorrections, u, vtsq = outFld[1], outFld[2], outFld[3]

   -- Compute boundary corrections.
   ffiC.gkyl_mom_calc_bcorr_advance_cu(self._zero_bcorr_calc, fIn:localRange(), u:localRange(), fIn._zeroDevice, boundaryCorrections._zeroDevice)

   -- Compute u and vtsq.
   ffiC.gkyl_prim_lbo_calc_advance_cu(self._zero_prim_calc, self.confBasis._zero, u:localRange(), moments._zeroDevice, boundaryCorrections._zeroDevice, u._zeroDevice, vtsq._zeroDevice)

end
   

return SelfPrimMoments
