-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the cross-primitive moments for the cross species 
-- collisions given the masses and primitive moments (and maybe Greene's beta)
-- of two species.
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

-- Function to check if operator option is correct.
local function isOperatorGood(nm)
   if nm == "VmLBO" or nm == "GkLBO" then
      return true
   end
   return false
end

-- Function to check if formulas option is correct.
local function isFormulasGood(nm)
   if nm=="HeavyIons" or nm=="GreenSmallAngle" or
      nm=="GreeneSmallAngleLimit" or nm=="Greene" then
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

   self._formulas = assert(
      tbl.formulas, "Updater.CrossPrimMoments: Must specify which formulas to use (HeavyIons, Greene, GreeneSmallAngle) using 'formulas'.")
   self._formulasKernel     = self._formulas

   self._isGreeneSmallAngle = false    -- MF: Better to check a boolean than compare a string in advance method.
   if self._formulas == "GreeneSmallAngle" then
      self._formulasKernel     = "Greene"
      self._isGreeneSmallAngle = true
   end
   self._isHeavyIons = false    -- MF: Better to check a boolean than compare a string in advance method.
   if self._formulas == "HeavyIons" then
      self._isHeavyIons = true
   end

   -- Free-parameter in John Greene's equations.
   self._beta = tbl.betaGreene

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

   self._numBasisP = phaseBasis:numBasis()
   self._numBasisC = confBasis:numBasis()

   self._id, self._polyOrder = confBasis:id(), confBasis:polyOrder()

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)
end

-- Advance method.
function CrossPrimMoments:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   -- Subscripts 1 and 2 refer to first and second species.
   -- Species 1 is the negative-charge species, and 2 the positive-charge species.
   -- For electron-ion, for example, one should think of 1 as the electrons
   -- and 2 as the ions.
   local collTermSub = inFld[1]    -- Collision term subscripts, 12 or 21.
   local crossPrimMomentsCalc = PrimMomentsDecl.selectCrossPrimMomentsCalc(self._kinSpecies, self._formulasKernel, self._id, self._cDim, self._vDim, self._polyOrder)

   -- Ratio of mass of other species to self species mass.
   -- or if "HeavyIons", ratio of lighter to heavier species' mass.
   local mRat = inFld[2]
   if self._isHeavyIons and collTermSub=="12" then
      mRat = 1.0
   end

   -- n1, u1, vtSq1 correspond to the self-species and n2, u2, vtSq2 
   -- to the other species (or lighter and heavier, respectively, if
   -- using "HeavyIons").
   local n1Fld, n2Fld, u1Fld, u2Fld, vtSq1Fld, vtSq2Fld
   n1Fld, u1Fld, vtSq1Fld = inFld[3], inFld[4], inFld[5]
   n2Fld, u2Fld, vtSq2Fld = inFld[6], inFld[7], inFld[8]

   local uCross    = outFld[1]
   local vtSqCross = outFld[2]

   local confRange = u1Fld:localRange()
   if self.onGhosts then confRange = u1Fld:localExtRange() end

   -- construct ranges for nested loops
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = confRange:selectFirst(self._cDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId() -- local thread ID

   local u1FldIndexer    = u1Fld:genIndexer()
   local vtSq1FldIndexer = vtSq1Fld:genIndexer()
   local u2FldIndexer    = u2Fld:genIndexer()
   local vtSq2FldIndexer = vtSq2Fld:genIndexer()

   local u1FldItr     = u1Fld:get(1)
   local vtSq1FldItr  = vtSq1Fld:get(1)
   local u2FldItr     = u2Fld:get(1)
   local vtSq2FldItr  = vtSq2Fld:get(1)

   local uCrossIndexer    = uCross:genIndexer()
   local vtSqCrossIndexer = vtSqCross:genIndexer()

   local uCrossItr    = uCross:get(1)
   local vtSqCrossItr = vtSqCross:get(1)

   if self._isGreeneSmallAngle then

      -- In this case number densities are used to compute beta.
      local n1FldIndexer = n1Fld:genIndexer()
      local n2FldIndexer = n2Fld:genIndexer()
      local n1FldItr     = n1Fld:get(1)
      local n2FldItr     = n2Fld:get(1)

      -- To obtain the cell average, multiply the zeroth coefficient by this factor.
      local massRatFac = 2.0*(1.0+mRat) 

      for confIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(confIdx)

         n1Fld:fill(n1FldIndexer(confIdx), n1FldItr)
         u1Fld:fill(u1FldIndexer(confIdx), u1FldItr)
         vtSq1Fld:fill(vtSq1FldIndexer(confIdx), vtSq1FldItr)
         n2Fld:fill(n2FldIndexer(confIdx), n2FldItr)
         u2Fld:fill(u2FldIndexer(confIdx), u2FldItr)
         vtSq2Fld:fill(vtSq2FldIndexer(confIdx), vtSq2FldItr)
         
         uCross:fill(uCrossIndexer(confIdx), uCrossItr)
         vtSqCross:fill(vtSqCrossIndexer(confIdx), vtSqCrossItr)

         self._beta   = (n1FldItr[1]/n2FldItr[1])*massRatFac/math.sqrt((1.0+vtSq2FldItr[1]/vtSq1FldItr[1])^3)-1.0

         crossPrimMomentsCalc(mRat, self._beta, u1FldItr:data(), vtSq1FldItr:data(), u2FldItr:data(), vtSq2FldItr:data(), uCrossItr:data(), vtSqCrossItr:data())
      end

   else

      -- Configuration space loop, computing cross-primitive moments in each cell.
      for confIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(confIdx)

         u1Fld:fill(u1FldIndexer(confIdx), u1FldItr)
         vtSq1Fld:fill(vtSq1FldIndexer(confIdx), vtSq1FldItr)
         u2Fld:fill(u2FldIndexer(confIdx), u2FldItr)
         vtSq2Fld:fill(vtSq2FldIndexer(confIdx), vtSq2FldItr)
         
         uCross:fill(uCrossIndexer(confIdx), uCrossItr)
         vtSqCross:fill(vtSqCrossIndexer(confIdx), vtSqCrossItr)

         crossPrimMomentsCalc(mRat, self._beta, u1FldItr:data(), vtSq1FldItr:data(), u2FldItr:data(), vtSq2FldItr:data(), uCrossItr:data(), vtSqCrossItr:data())
      end

   end
end

return CrossPrimMoments
