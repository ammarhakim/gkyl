-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the cross-primitive moments, u_ei, u_ie, vtSq_ei=T_ei/m_e,
-- and vtSq_ie=T_ie/m_i for the cross species collisions given the primitive
-- moments of the species.
-- Here the subscript ei means that it corresponds to the contribution to
-- electrons due to collisions with ions, and vice-versa for ie.
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

-- function to check if collide option is correct
local function isCollideGood(nm)
   if nm == "Vmei" or nm == "Vmie" or nm == "Vmall" or
      nm == "Gkei" or nm == "Gkie" or nm == "Gkall" then
      return true
   end
   return false
end

-- function to check if operator option is correct
local function isOperatorGood(nm)
   if nm == "BGK" or nm == "LBO" then
      return true
   end
   return false
end

-- Primitive moments updater object.
local CrossPrimMoments = Proto(UpdaterBase)

function CrossPrimMoments:init(tbl)
   CrossPrimMoments.super.init(self, tbl) -- setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.CrossPrimMoments: Must provide grid object using 'onGrid'.")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.CrossPrimMoments: Must provide the phase basis object using 'phaseBasis'.")

   local confBasis = assert(
      tbl.confBasis, "Updater.CrossPrimMoments: Must provide the configuration basis object using 'confBasis'.")

   self._massRat = assert(
      tbl.massRatio, "Updater.CrossPrimMoments: Must provide the mass ratio (mi/me) using 'massRatio'.")

   local operator = assert(
      tbl.operator, "Updater.CrossPrimMoments: Must specify the collision operator (BGK or LBO) using 'operator'.")

   local collide = assert(
      tbl.collide, "Updater.CrossPrimMoments: Must specify which species to collide (Vmei, Vmie, Vmall, Gkei, Gkie, Gkall) using 'collide'.")

   -- dimension of spaces.
   self._pDim = phaseBasis:ndim()
   -- ensure sanity.
   assert(phaseBasis:polyOrder() == confBasis:polyOrder(),
          "Updater.CrossPrimMoments: Polynomial orders of phase and conf basis must match.")
   assert(phaseBasis:id() == confBasis:id(),
          "Updater.CrossPrimMoments: Type of phase and conf basis must match.")
   -- determine configuration and velocity space dims.
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   self._numBasisP = phaseBasis:numBasis()
   self._numBasisC = confBasis:numBasis()

   local id, polyOrder = confBasis:id(), confBasis:polyOrder()

   if isOperatorGood(operator) then
      if isCollideGood(collide) then
         self._eiColl = false
         self._ieColl = false
         if ((collide == "Vmall") or (collide == "Vmei")) or 
            ((collide == "Gkall") or (collide == "Gkei")) then self._eiColl = true end
         if ((collide == "Vmall") or (collide == "Vmie")) or
            ((collide == "Gkall") or (collide == "Gkie")) then self._ieColl = true end
         self._allColl = self._eiColl and self._ieColl

         if self._eiColl then
           self._eiCrossPrimMomentsCalc = PrimMomentsDecl.selectCrossPrimMomentsCalc(operator, collide, id, self._cDim, self._vDim, polyOrder)
         end
         if self._ieColl then
           self._ieCrossPrimMomentsCalc = PrimMomentsDecl.selectCrossPrimMomentsCalc(operator, collide, id, self._cDim, self._vDim, polyOrder)
         end
      else
         assert(false, string.format(
                   "CrossPrimMoments: Collide option must be 'Vmei', 'Vmie', 'Vmall', 'Gkei', 'Gkie' or 'Gkall'. Requested %s instead.", collide))
      end
   else
      assert(false, string.format(
                "CrossPrimMoments: Operator option must be 'BGK' or 'LBO'. Requested %s instead.", operator))
   end

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)

end

-- advance method
function CrossPrimMoments:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   -- Need to allow the cases of ei and ie collisions,
   -- ei collisions alone, and ie collisions alone.
   if self._allColl then
      local UeiOut    = outFld[1]
      local vtSqeiOut = outFld[2]
      local UieOut    = outFld[3]
      local vtSqieOut = outFld[4]
   elseif self._eiColl then
      local UeiOut    = outFld[1]
      local vtSqeiOut = outFld[2]
   elseif self._ieColl then
      local UieOut    = outFld[1]
      local vtSqieOut = outFld[2]
   end

   local neFld, UeFld, vtSqeFld, niFld, UiFld, vtSqiFld
   neFld, UeFld, vtSqeFld = inFld[1], inFld[2], inFld[3]
   niFld, UiFld, vtSqiFld = inFld[4], inFld[5], inFld[6]

   local confRange = UeFld:localRange()
   if self.onGhosts then confRange = UeFld:localExtRange() end

   local neFldIndexer    = neFld:genIndexer()
   local UeFldIndexer    = UeFld:genIndexer()
   local vtSqeFldIndexer = vtSqeFld:genIndexer()
   local niFldIndexer    = niFld:genIndexer()
   local UiFldIndexer    = UiFld:genIndexer()
   local vtSqiFldIndexer = vtSqiFld:genIndexer()

   local neFldItr     = neFld:get(1)
   local UeFldItr     = UeFld:get(1)
   local vtSqeFldItr  = vtSqeFld:get(1)
   local niFldItr     = niFld:get(1)
   local UiFldItr     = UiFld:get(1)
   local vtSqiFldItr  = vtSqiFld:get(1)

   if self._allColl or self._eiColl then
      local UeiOutIndexer    = UeiOut:genIndexer()
      local vtSqeiOutIndexer = vtSqeiOut:genIndexer()

      local UeiOutItr    = UeiOut:get(1)
      local vtSqeiOutItr = vtSqeiOut:get(1)
   end
   if self._allColl or self._ieColl then
      local UieOutIndexer    = UieOut:genIndexer()
      local vtSqieOutIndexer = vtSqieOut:genIndexer()

      local UieOutItr    = UieOut:get(1)
      local vtSqieOutItr = vtSqieOut:get(1)
   end

   -- configuration space loop, computing primitive moments in each cell
   for confIdx in confRange:colMajorIter() do
      grid:setIndex(confIdx)

      neFld:fill(neFldIndexer(confIdx), neFldItr)
      UeFld:fill(UeFldIndexer(confIdx), UeFldItr)
      vtSqeFld:fill(vtSqeFldIndexer(confIdx), vtSqeFldItr)
      niFld:fill(niFldIndexer(confIdx), niFldItr)
      UiFld:fill(UiFldIndexer(confIdx), UiFldItr)
      vtSqiFld:fill(vtSqiFldIndexer(confIdx), vtSqiFldItr)
      
      if self._eiColl then
         UeiOut:fill(UeiOutIndexer(confIdx), UeiOutItr)
         vtSqeiOut:fill(vtSqeiOutIndexer(confIdx), vtSqeiOutItr)

         self._eiCrossPrimMomentsCalc(self._massRat, neFldItr:data(), UeFldItr:data(), vtSqeFldItr:data(), niFldItr:data(), UiFldItr:data(), vtSqiFldItr:data(), UeiOutItr:data(), vtSqeiOutItr:data())
      end

      if self._ieColl then
         UieOut:fill(UieOutIndexer(confIdx), UieOutItr)
         vtSqieOut:fill(vtSqieOutIndexer(confIdx), vtSqieOutItr)

         self._ieCrossPrimMomentsCalc(self._massRat, neFldItr:data(), UeFldItr:data(), vtSqeFldItr:data(), niFldItr:data(), UiFldItr:data(), vtSqiFldItr:data(), UieOutItr:data(), vtSqieOutItr:data())
      end
   end
end

return CrossPrimMoments
