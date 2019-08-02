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

   local confBasis = assert(
      tbl.confBasis, "Updater.CrossPrimMoments: Must provide the configuration basis object using 'confBasis'.")

   self._operator = assert(
      tbl.operator, "Updater.CrossPrimMoments: Must specify the collision operator (VmLBO, GkLBO, VmBGK or GkBGK) using 'operator'.")
   assert(isOperatorGood(self._operator), string.format("CrossPrimMoments: Operator option must be 'VmLBO', 'GkLBO', 'VmBGK' or 'GkBGK'. Requested %s instead.", self._operator))

   -- Indicate if collisionality is spatially varying and if it is cell-wise constant.
   self._varNu       = tbl.varyingNu
   self._cellConstNu = tbl.useCellAverageNu 

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

   self._numBasisC = confBasis:numBasis()

   self._basisID, self._polyOrder = confBasis:id(), confBasis:polyOrder()

   local uDim = self._vDim
   if self._operator=="GkLBO" or self._operator=="GkBGK" then
      uDim = 1
   end
   self._isLBO = false
   if self._operator=="VmLBO" or self._operator=="GkLBO" then
      self._isLBO = true
   end

   -- Need two Eigen matrices: one to divide by (ms*nusr*m0s+mr*nurs*m0r)
   -- and one to compute cross-primitive moments.
   self._binOpData    = ffiC.new_binOpData_t(self._numBasisC*2*(uDim+1), 0)
   self._binOpDataDiv = ffiC.new_binOpData_t(self._numBasisC, 0)

   -- To obtain the cell average, multiply the zeroth coefficient by this factor.
   self._cellAvFac = 1.0/math.sqrt(2.0^self._cDim)

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)
end

-- Advance method.
function CrossPrimMoments:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   local crossPrimMomentsCalc = PrimMomentsDecl.selectCrossPrimMomentsCalc(self._operator, self._basisID, self._cDim, self._vDim, self._polyOrder)

   local mSelf, nuSelfIn, m0Self, uSelf, vtSqSelf
   local mOther, nuOtherIn, m0Other, uOther, vtSqOther
   mSelf, nuSelfIn   = inFld[1], inFld[2]
   m0Self            = inFld[3][1]
   uSelf, vtSqSelf   = inFld[4][1], inFld[4][2]

   mOther, nuOtherIn = inFld[5], inFld[6]
   m0Other           = inFld[7][1]
   uOther, vtSqOther = inFld[8][1], inFld[8][2]

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
   local uSelfItr    = uSelf:get(1)
   local vtSqSelfItr = vtSqSelf:get(1)

   local m0OtherItr   = m0Other:get(1)
   local uOtherItr    = uOther:get(1)
   local vtSqOtherItr = vtSqOther:get(1)

   local uCrossSelfItr     = uCrossSelf:get(1)
   local vtSqCrossSelfItr  = vtSqCrossSelf:get(1)
   local uCrossOtherItr    = uCrossOther:get(1)
   local vtSqCrossOtherItr = vtSqCrossOther:get(1)

   -- For LBO need a few more inputs. Also, in order to avoid evaluation
   -- of if-statements in the spatial loop, we have separate loops for each
   -- operator, and for each of cellConstNu. Makes the code longer but slightly more efficient.
   local m1Self, m2Self, m1CorrectionSelf, m2CorrectionSelf
   local m1Other, m2Other, m1CorrectionOther, m2CorrectionOther
   local m0StarSelf, m1StarSelf, m2StarSelf
   local m0StarOther, m1StarOther, m2StarOther
   local m0StarSelfItr, m1StarSelfItr, m2StarSelfItr
   local m0StarOtherItr, m1StarOtherItr, m2StarOtherItr
   local m1SelfItr, m2SelfItr, m1CorrectionSelfItr, m2CorrectionSelfItr
   local m1OtherItr, m2OtherItr, m1CorrectionOtherItr, m2CorrectionOtherItr
   local nuSelf, nuSelfItr
   local nuOther, nuOtherItr
   if self._cellConstNu then

      if self._varNu then
         nuSelfItr  = nuSelfIn:get(1) 
         nuOtherItr = nuOtherIn:get(1)
      end

      if self._isLBO then
         m1Self, m2Self                       = inFld[3][2], inFld[3][3]
         m1Other, m2Other                     = inFld[7][2], inFld[7][3]
         m1CorrectionSelf, m2CorrectionSelf   = inFld[9][1], inFld[9][2]
         m1CorrectionOther, m2CorrectionOther = inFld[11][1], inFld[11][2]

         if self._polyOrder == 1 then
            m0StarSelf, m1StarSelf, m2StarSelf = inFld[10][1], inFld[10][2], inFld[10][3]
            m0StarSelfItr = m0StarSelf:get(1)
            m1StarSelfItr = m1StarSelf:get(1)
            m2StarSelfItr = m2StarSelf:get(1)

            m0StarOther, m1StarOther, m2StarOther = inFld[12][1], inFld[12][2], inFld[12][3]
            m0StarOtherItr = m0StarOther:get(1)
            m1StarOtherItr = m1StarOther:get(1)
            m2StarOtherItr = m2StarOther:get(1)
         end

         m1SelfItr            = m1Self:get(1)
         m2SelfItr            = m2Self:get(1)
         m1CorrectionSelfItr  = m1CorrectionSelf:get(1)
         m2CorrectionSelfItr  = m2CorrectionSelf:get(1)
         m1OtherItr           = m1Other:get(1)
         m2OtherItr           = m2Other:get(1)
         m1CorrectionOtherItr = m1CorrectionOther:get(1)
         m2CorrectionOtherItr = m2CorrectionOther:get(1)

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

               if self._varNu then
                  nuSelfIn:fill(confIndexer(cIdx), nuSelfItr)
                  nuOtherIn:fill(confIndexer(cIdx), nuOtherItr)

                  nuSelf  = nuSelfItr[0]*self._cellAvFac 
                  nuOther = nuOtherItr[0]*self._cellAvFac 
               else
                  nuSelf  = nuSelfIn
                  nuOther = nuOtherIn
               end
         
               crossPrimMomentsCalc(self._binOpData, self._binOpDataDiv, self._betaP1, mSelf, nuSelf, m0SelfItr:data(), m1SelfItr:data(), m2SelfItr:data(), uSelfItr:data(), vtSqSelfItr:data(), m1CorrectionSelfItr:data(), m2CorrectionSelfItr:data(), mOther, nuOther, m0OtherItr:data(), m1OtherItr:data(), m2OtherItr:data(), uOtherItr:data(), vtSqOtherItr:data(), m1CorrectionOtherItr:data(), m2CorrectionOtherItr:data(), uCrossSelfItr:data(), vtSqCrossSelfItr:data(), uCrossOtherItr:data(), vtSqCrossOtherItr:data())
            end

         else    -- Piecewise linear polynomial basis below.

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
         
               if self._varNu then
                  nuSelfIn:fill(confIndexer(cIdx), nuSelfItr)
                  nuOtherIn:fill(confIndexer(cIdx), nuOtherItr)

                  nuSelf  = nuSelfItr[0]*self._cellAvFac 
                  nuOther = nuOtherItr[0]*self._cellAvFac 
               else
                  nuSelf  = nuSelfIn
                  nuOther = nuOtherIn
               end
         
               crossPrimMomentsCalc(self._binOpData, self._binOpDataDiv, self._betaP1, mSelf, nuSelf, m0SelfItr:data(), m1SelfItr:data(), m2SelfItr:data(), uSelfItr:data(), vtSqSelfItr:data(), m1CorrectionSelfItr:data(), m2CorrectionSelfItr:data(), m0StarSelfItr:data(), m1StarSelfItr:data(), m2StarSelfItr:data(), mOther, nuOther, m0OtherItr:data(), m1OtherItr:data(), m2OtherItr:data(), uOtherItr:data(), vtSqOtherItr:data(), m1CorrectionOtherItr:data(), m2CorrectionOtherItr:data(), m0StarOtherItr:data(), m1StarOtherItr:data(), m2StarOtherItr:data(), uCrossSelfItr:data(), vtSqCrossSelfItr:data(), uCrossOtherItr:data(), vtSqCrossOtherItr:data())
            end

         end    -- end if polyOrder>1.
      else    -- BGK operator below (needs less inputs).

         for cIdx in confRangeDecomp:rowMajorIter(tId) do
            grid:setIndex(cIdx)
         
            m0Self:fill(confIndexer(cIdx), m0SelfItr)
            uSelf:fill(confIndexer(cIdx), uSelfItr)
            vtSqSelf:fill(confIndexer(cIdx), vtSqSelfItr)
         
            m0Other:fill(confIndexer(cIdx), m0OtherItr)
            uOther:fill(confIndexer(cIdx), uOtherItr)
            vtSqOther:fill(confIndexer(cIdx), vtSqOtherItr)
         
            uCrossSelf:fill(confIndexer(cIdx), uCrossSelfItr)
            vtSqCrossSelf:fill(confIndexer(cIdx), vtSqCrossSelfItr)
            uCrossOther:fill(confIndexer(cIdx), uCrossOtherItr)
            vtSqCrossOther:fill(confIndexer(cIdx), vtSqCrossOtherItr)
         
            if self._varNu then
               nuSelfIn:fill(confIndexer(cIdx), nuSelfItr)
               nuOtherIn:fill(confIndexer(cIdx), nuOtherItr)

               nuSelf  = nuSelfItr[0]*self._cellAvFac 
               nuOther = nuOtherItr[0]*self._cellAvFac 
            else
               nuSelf  = nuSelfIn
               nuOther = nuOtherIn
            end
         
            crossPrimMomentsCalc(self._binOpDataDiv, self._betaP1, mSelf, nuSelf, m0SelfItr:data(), uSelfItr:data(), vtSqSelfItr:data(), mOther, nuOther, m0OtherItr:data(), uOtherItr:data(), vtSqOtherItr:data(), uCrossSelfItr:data(), vtSqCrossSelfItr:data(), uCrossOtherItr:data(), vtSqCrossOtherItr:data())
         end

      end    -- end if self._isLBO.
   else    -- Below: collisionality is not cell-wise constant.
   end

end

return CrossPrimMoments
