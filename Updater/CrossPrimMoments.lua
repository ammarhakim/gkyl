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
local LinearDecomp    = require "Lib.LinearDecomp"
local Proto           = require "Lib.Proto"
local PrimMomentsDecl = require "Updater.primMomentsCalcData.PrimMomentsModDecl"
local xsys            = require "xsys"
local ffi             = require "ffi"

local ffiC = ffi.C
ffi.cdef [[
  void gkylCopyToField(double *f, double *data, unsigned numComponents, unsigned c);
  void gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out);
]]

-- Function to check if operator option is correct.
local function isOperatorGood(nm)
   if nm == "VmLBO" or nm == "GkLBO" then
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

   local confBasis = assert(
      tbl.confBasis, "Updater.CrossPrimMoments: Must provide the configuration basis object using 'confBasis'.")

   local operator = assert(
      tbl.operator, "Updater.CrossPrimMoments: Must specify the collision operator (VmLBO, or GkLBO) using 'operator'.")
   assert(isOperatorGood(operator), string.format("CrossPrimMoments: Operator option must be 'VmLBO' or 'GkLBO'. Requested %s instead.", operator))
   if operator=="VmLBO" then
      self._kinSpecies = "Vm"
   elseif operator=="GkLBO" then
      self._kinSpecies = "Gk"
   end

   -- Free-parameter in John Greene's equations.
   self._beta   = tbl.betaGreene
   self._betaP1 = self._beta+1

   -- Dimension of spaces.
   self._pDim = phaseBasis:ndim()
   -- Ensure sanity.
   assert(phaseBasis:polyOrder() == confBasis:polyOrder(),
          "Updater.CrossPrimMoments: Polynomial orders of phase and conf basis must match.")
   assert(phaseBasis:id() == confBasis:id(),
          "Updater.CrossPrimMoments: Type of phase and conf basis must match.")
   -- Determine configuration and velocity space dims.
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   self._uDim = self._vDim
   if self._kinSpecies == "Gk" then
      self._uDim = 1
   end

   self._numBasisC = confBasis:numBasis()

   self._basisID, self._polyOrder = confBasis:id(), confBasis:polyOrder()

   -- Need two Eigen matrices: one to divide by (ms*nusr*m0s+mr*nurs*m0r)
   -- and one to compute cross-primitive moments.
   self._binOpData    = ffiC.new_binOpData_t(self._numBasisC*2*(self._uDim+1), 0)
   self._binOpDataDiv = ffiC.new_binOpData_t(self._numBasisC, 0)

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)
end

-- Advance method.
function CrossPrimMoments:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   local crossPrimMomentsCalc = PrimMomentsDecl.selectCrossPrimMomentsCalc(self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)

   local mSelf, nuSelf, m0Self, m1Self, m2Self, uSelf, vtSqSelf, m1CorrectionSelf, m2CorrectionSelf
   local mOther, nuOther, m0Other, m1Other, m2Other, uOther, vtSqOther, m1CorrectionOther, m2CorrectionOther
   mSelf, nuSelf                        = inFld[1], inFld[2]
   m0Self, m1Self, m2Self               = inFld[3][1], inFld[3][2], inFld[3][3]
   uSelf, vtSqSelf                      = inFld[4][1], inFld[4][2]
   m1CorrectionSelf, m2CorrectionSelf   = inFld[5][1], inFld[5][2]

   mOther, nuOther                      = inFld[7], inFld[8]
   m0Other, m1Other, m2Other            = inFld[9][1], inFld[9][2], inFld[9][3]
   uOther, vtSqOther                    = inFld[10][1], inFld[10][2]
   m1CorrectionOther, m2CorrectionOther = inFld[11][1], inFld[11][2]

   local m0StarSelf, m1StarSelf, m2StarSelf
   local m0StarOther, m1StarOther, m2StarOther
   local m0StarSelfItr, m1StarSelfItr, m2StarSelfItr
   local m0StarOtherItr, m1StarOtherItr, m2StarOtherItr
   if self._polyOrder == 1 then
      m0StarSelf, m1StarSelf, m2StarSelf = inFld[6][1], inFld[6][2], inFld[6][3]
      m0StarSelfItr     = m0StarSelf:get(1)
      m1StarSelfItr     = m1StarSelf:get(1)
      m2StarSelfItr     = m2StarSelf:get(1)

      m0StarOther, m1StarOther, m2StarOther = inFld[12][1], inFld[12][2], inFld[12][3]
      m0StarOtherItr     = m0StarOther:get(1)
      m1StarOtherItr     = m1StarOther:get(1)
      m2StarOtherItr     = m2StarOther:get(1)
   end

   local uCrossSelf     = outFld[1]
   local vtSqCrossSelf  = outFld[2]
   local uCrossOther    = outFld[3]
   local vtSqCrossOther = outFld[4]

   local confRange = uSelf:localRange()
   if self.onGhosts then confRange = uSelf:localExtRange() end

   -- Construct ranges for nested loops.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = confRange:selectFirst(self._cDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId() -- Local thread ID.

   local confIndexer   = m0Self:genIndexer()

   local m0SelfItr   = m0Self:get(1)
   local m1SelfItr   = m1Self:get(1)
   local m2SelfItr   = m2Self:get(1)
   local uSelfItr    = uSelf:get(1)
   local vtSqSelfItr = vtSqSelf:get(1)
   local m1CorrectionSelfItr = m1CorrectionSelf:get(1)
   local m2CorrectionSelfItr = m2CorrectionSelf:get(1)

   local m0OtherItr   = m0Other:get(1)
   local m1OtherItr   = m1Other:get(1)
   local m2OtherItr   = m2Other:get(1)
   local uOtherItr    = uOther:get(1)
   local vtSqOtherItr = vtSqOther:get(1)
   local m1CorrectionOtherItr = m1CorrectionOther:get(1)
   local m2CorrectionOtherItr = m2CorrectionOther:get(1)

   local uCrossSelfItr     = uCrossSelf:get(1)
   local vtSqCrossSelfItr  = vtSqCrossSelf:get(1)
   local uCrossOtherItr    = uCrossOther:get(1)
   local vtSqCrossOtherItr = vtSqCrossOther:get(1)

   -- polyOrder=1 and >1 each use separate velocity grid loops to
   -- avoid evaluating (if polyOrder==1) at each cell.
   if self._polyOrder > 1 then

      for cIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(cIdx)
   
         m0Self:fill(confIndexer(cIdx), m0SelfItr)
         m1Self:fill(confIndexer(cIdx), m1SelfItr)
         m2Self:fill(confIndexer(cIdx), m2SelfItr)
         uSelf:fill(confIndexer(cIdx), uSelfItr)
         vtSqSelf:fill(confIndexer(cIdx), vtSqSelfItr)
         m1CorrectionSelf:fill(confIndexer(cIdx), m1CorrectionSelfItr)
         m2CorrectionSelf:fill(confIndexer(cIdx), m2CorrectionSelfItr)
   
         m0Other:fill(confIndexer(cIdx), m0OtherItr)
         m1Other:fill(confIndexer(cIdx), m1OtherItr)
         m2Other:fill(confIndexer(cIdx), m2OtherItr)
         uOther:fill(confIndexer(cIdx), uOtherItr)
         vtSqOther:fill(confIndexer(cIdx), vtSqOtherItr)
         m1CorrectionOther:fill(confIndexer(cIdx), m1CorrectionOtherItr)
         m2CorrectionOther:fill(confIndexer(cIdx), m2CorrectionOtherItr)
   
         uCrossSelf:fill(confIndexer(cIdx), uCrossSelfItr)
         vtSqCrossSelf:fill(confIndexer(cIdx), vtSqCrossSelfItr)
         uCrossOther:fill(confIndexer(cIdx), uCrossOtherItr)
         vtSqCrossOther:fill(confIndexer(cIdx), vtSqCrossOtherItr)
   
         crossPrimMomentsCalc(self._binOpData, self._binOpDataDiv, self._betaP1, mSelf, nuSelf, m0SelfItr:data(), m1SelfItr:data(), m2SelfItr:data(), uSelfItr:data(), vtSqSelfItr:data(), m1CorrectionSelfItr:data(), m2CorrectionSelfItr:data(), mOther, nuOther, m0OtherItr:data(), m1OtherItr:data(), m2OtherItr:data(), uOtherItr:data(), vtSqOtherItr:data(), m1CorrectionOtherItr:data(), m2CorrectionOtherItr:data(), uCrossSelfItr:data(), vtSqCrossSelfItr:data(), uCrossOtherItr:data(), vtSqCrossOtherItr:data())
      end

   else

      for cIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(cIdx)
   
         m0Self:fill(confIndexer(cIdx), m0SelfItr)
         m1Self:fill(confIndexer(cIdx), m1SelfItr)
         m2Self:fill(confIndexer(cIdx), m2SelfItr)
         uSelf:fill(confIndexer(cIdx), uSelfItr)
         vtSqSelf:fill(confIndexer(cIdx), vtSqSelfItr)
         m1CorrectionSelf:fill(confIndexer(cIdx), m1CorrectionSelfItr)
         m2CorrectionSelf:fill(confIndexer(cIdx), m2CorrectionSelfItr)
         m0StarSelf:fill(confIndexer(cIdx), m0StarSelfItr)
         m1StarSelf:fill(confIndexer(cIdx), m1StarSelfItr)
         m2StarSelf:fill(confIndexer(cIdx), m2StarSelfItr)
   
         m0Other:fill(confIndexer(cIdx), m0OtherItr)
         m1Other:fill(confIndexer(cIdx), m1OtherItr)
         m2Other:fill(confIndexer(cIdx), m2OtherItr)
         uOther:fill(confIndexer(cIdx), uOtherItr)
         vtSqOther:fill(confIndexer(cIdx), vtSqOtherItr)
         m1CorrectionOther:fill(confIndexer(cIdx), m1CorrectionOtherItr)
         m2CorrectionOther:fill(confIndexer(cIdx), m2CorrectionOtherItr)
         m0StarOther:fill(confIndexer(cIdx), m0StarOtherItr)
         m1StarOther:fill(confIndexer(cIdx), m1StarOtherItr)
         m2StarOther:fill(confIndexer(cIdx), m2StarOtherItr)
   
         uCrossSelf:fill(confIndexer(cIdx), uCrossSelfItr)
         vtSqCrossSelf:fill(confIndexer(cIdx), vtSqCrossSelfItr)
         uCrossOther:fill(confIndexer(cIdx), uCrossOtherItr)
         vtSqCrossOther:fill(confIndexer(cIdx), vtSqCrossOtherItr)
   
         crossPrimMomentsCalc(self._binOpData, self._binOpDataDiv, self._betaP1, mSelf, nuSelf, m0SelfItr:data(), m1SelfItr:data(), m2SelfItr:data(), uSelfItr:data(), vtSqSelfItr:data(), m1CorrectionSelfItr:data(), m2CorrectionSelfItr:data(), m0StarSelfItr:data(), m1StarSelfItr:data(), m2StarSelfItr:data(), mOther, nuOther, m0OtherItr:data(), m1OtherItr:data(), m2OtherItr:data(), uOtherItr:data(), vtSqOtherItr:data(), m1CorrectionOtherItr:data(), m2CorrectionOtherItr:data(), m0StarOtherItr:data(), m1StarOtherItr:data(), m2StarOtherItr:data(), uCrossSelfItr:data(), vtSqCrossSelfItr:data(), uCrossOtherItr:data(), vtSqCrossOtherItr:data())
      end

   end

end

return CrossPrimMoments
