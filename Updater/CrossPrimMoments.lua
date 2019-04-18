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

-- Function to check if formulas option is correct.
local function isFormulasGood(nm)
   if nm=="GreenSmallAngle" or
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

   self._phaseGrid = assert(
      tbl.phaseGrid, "Updater.SelfPrimMoments: Must provide phase grid object using 'phaseGrid'.")

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
      tbl.formulas, "Updater.CrossPrimMoments: Must specify which formulas to use (Greene, GreeneSmallAngle) using 'formulas'.")
   self._formulasKernel     = self._formulas

   self._isGreeneSmallAngle = false    -- MF: Better to check a boolean than compare a string in advance method.
   if self._formulas == "GreeneSmallAngle" then
      self._formulasKernel     = "Greene"
      self._isGreeneSmallAngle = true
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

   self._numBasisP = phaseBasis:numBasis()
   self._numBasisC = confBasis:numBasis()

   self._basisID, self._polyOrder = confBasis:id(), confBasis:polyOrder()

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)

   self._isFirst   = true
   self._perpRangeSelf  = {} -- perp ranges in velocity directions.
   self._perpRangeOther = {} -- perp ranges in velocity directions.

   self._uDim       = self._vDim  -- Dimensionality of flow velocity vector.
   self._isGkLBO = false
   self._kinSpecies = "Vm"        -- Vlasov-Maxwell species.
   self._phaseGridSelf = self._phaseGrid
--   self._phaseGridOther = self._phaseGrid
   self._SelfPrimMomentsCalc = PrimMomentsDecl.selectSelfPrimMomentsCalc(self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)

   self._binOpDataSelf = ffiC.new_binOpData_t(self._numBasisC*(self._uDim+1), 0)
   self._binOpData = ffiC.new_binOpData_t(self._numBasisC*2*(self._uDim+1), 0)
end

-- Advance method.
function CrossPrimMoments:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   -- Subscripts 1 and 2 refer to first and second species.
   -- Species 1 is the negative-charge species, and 2 the positive-charge species.
   -- For electron-ion, for example, one should think of 1 as the electrons
   -- and 2 as the ions.
   local crossPrimMomentsCalc = PrimMomentsDecl.selectCrossPrimMomentsCalc(self._kinSpecies, self._formulasKernel, self._basisID, self._cDim, self._vDim, self._polyOrder)

   -- n1, u1, vtSq1 correspond to the self-species and n2, u2, vtSq2 
   -- to the other species.
   local mSelf, nuSelf, m0Self, m1Self, m2Self, uSelf, vtSqSelf, fInSelf
   local mOther, nuOther, m0Other, m1Other, m2Other, uOther, vtSqOther, fInOther
   mSelf, nuSelf          = inFld[1], inFld[2]
   m0Self, m1Self, m2Self = inFld[3], inFld[4], inFld[5]
   uSelf, vtSqSelf        = inFld[6], inFld[7]
   fInSelf                = inFld[8]
   mOther, nuOther           = inFld[9], inFld[10]
   m0Other, m1Other, m2Other = inFld[11], inFld[12], inFld[13]
   uOther, vtSqOther         = inFld[14], inFld[15]
   fInOther                  = inFld[16]
   self._phaseGridOther = inFld[17]

   local mRat = mOther/mSelf

   local uCross    = outFld[1]
   local vtSqCross = outFld[2]

   local phaseRangeSelf  = fInSelf:localRange()
   local phaseRangeOther = fInOther:localRange()

   local confRange = uSelf:localRange()
   if self.onGhosts then confRange = uSelf:localExtRange() end

   -- construct ranges for nested loops
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = confRange:selectFirst(self._cDim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId() -- local thread ID

   local m0SelfIndexer    = m0Self:genIndexer()
   local m1SelfIndexer    = m1Self:genIndexer()
   local m2SelfIndexer    = m2Self:genIndexer()
   local uSelfIndexer     = uSelf:genIndexer()
   local vtSqSelfIndexer  = vtSqSelf:genIndexer()
   local phaseIndexerSelf = fInSelf:genIndexer()

   local m0OtherIndexer    = m0Other:genIndexer()
   local m1OtherIndexer    = m1Other:genIndexer()
   local m2OtherIndexer    = m2Other:genIndexer()
   local uOtherIndexer     = uOther:genIndexer()
   local vtSqOtherIndexer  = vtSqOther:genIndexer()
   local phaseIndexerOther = fInOther:genIndexer()

   local m0SelfItr   = m0Self:get(1)
   local m1SelfItr   = m1Self:get(1)
   local m2SelfItr   = m2Self:get(1)
   local uSelfItr    = uSelf:get(1)
   local vtSqSelfItr = vtSqSelf:get(1)
   local fInSelfItrP = fInSelf:get(1)
   local fInSelfItrM = fInSelf:get(1)

   local m0OtherItr   = m0Other:get(1)
   local m1OtherItr   = m1Other:get(1)
   local m2OtherItr   = m2Other:get(1)
   local uOtherItr    = uOther:get(1)
   local vtSqOtherItr = vtSqOther:get(1)
   local fInOtherItrP = fInOther:get(1)
   local fInOtherItrM = fInOther:get(1)

   local uCrossIndexer    = uCross:genIndexer()
   local vtSqCrossIndexer = vtSqCross:genIndexer()

   local uCrossItr    = uCross:get(1)
   local vtSqCrossItr = vtSqCross:get(1)

   -- These store the corrections to momentum and energy
   -- from the (velocity) boundary surface integrals.
   local cMomBSelf     = Lin.Vec(self._numBasisC*self._uDim)
   local cEnergyBSelf  = Lin.Vec(self._numBasisC)
   local cMomBOther    = Lin.Vec(self._numBasisC*self._uDim)
   local cEnergyBOther = Lin.Vec(self._numBasisC)

   -- Cell index, center and length left and right of a cell-boundary.
   local idxM, idxP = Lin.IntVec(pDim), Lin.IntVec(pDim)
   local xcM, xcP   = Lin.Vec(pDim), Lin.Vec(pDim)
   local dxM, dxP   = Lin.Vec(pDim), Lin.Vec(pDim)

   if self._isGreeneSmallAngle then

      -- In this case number densities are used to compute beta.
      local m0SelfIndexer = m0Self:genIndexer()
      local m0OtherIndexer = m0Other:genIndexer()
      local m0SelfItr     = m0Self:get(1)
      local m0OtherItr     = m0Other:get(1)

      -- To obtain the cell average, multiply the zeroth coefficient by this factor.
      local massRatFac = 2.0*(1.0+mRat) 

      for confIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(confIdx)

         m0Self:fill(m0SelfIndexer(confIdx), m0SelfItr)
         uSelf:fill(uSelfIndexer(confIdx), uSelfItr)
         vtSqSelf:fill(vtSqSelfIndexer(confIdx), vtSqSelfItr)
         m0Other:fill(m0OtherIndexer(confIdx), m0OtherItr)
         uOther:fill(uOtherIndexer(confIdx), uOtherItr)
         vtSqOther:fill(vtSqOtherIndexer(confIdx), vtSqOtherItr)
         
         uCross:fill(uCrossIndexer(confIdx), uCrossItr)
         vtSqCross:fill(vtSqCrossIndexer(confIdx), vtSqCrossItr)

         self._beta   = (m0SelfItr[1]/m0OtherItr[1])*massRatFac/math.sqrt((1.0+vtSqOtherItr[1]/vtSqSelfItr[1])^3)-1.0

         crossPrimMomentsCalc(mRat, self._beta, uSelfItr:data(), vtSqSelfItr:data(), uOtherItr:data(), vtSqOtherItr:data(), uCrossItr:data(), vtSqCrossItr:data())
      end

   else

      -- Configuration space loop, computing cross-primitive moments in each cell.
      for confIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(confIdx)

         m0Self:fill(m0SelfIndexer(confIdx), m0SelfItr)
         m1Self:fill(m1SelfIndexer(confIdx), m1SelfItr)
         m2Self:fill(m2SelfIndexer(confIdx), m2SelfItr)
         uSelf:fill(uSelfIndexer(confIdx), uSelfItr)
         vtSqSelf:fill(vtSqSelfIndexer(confIdx), vtSqSelfItr)

         m0Other:fill(m0OtherIndexer(confIdx), m0OtherItr)
         m1Other:fill(m1OtherIndexer(confIdx), m1OtherItr)
         m2Other:fill(m2OtherIndexer(confIdx), m2OtherItr)
         uOther:fill(uOtherIndexer(confIdx), uOtherItr)
         vtSqOther:fill(vtSqOtherIndexer(confIdx), vtSqOtherItr)
         
         uCross:fill(uCrossIndexer(confIdx), uCrossItr)
         vtSqCross:fill(vtSqCrossIndexer(confIdx), vtSqCrossItr)

         -- Compute the corrections to u and vtSq due to
         -- finite velocity grid to ensure conservation.
         ffiC.gkylCartFieldAssignAll(0, self._numBasisC, 0.0, cEnergyBSelf:data())
         ffiC.gkylCartFieldAssignAll(0, self._uDim*self._numBasisC, 0.0, cMomBSelf:data())
         ffiC.gkylCartFieldAssignAll(0, self._numBasisC, 0.0, cEnergyBOther:data())
         ffiC.gkylCartFieldAssignAll(0, self._uDim*self._numBasisC, 0.0, cMomBOther:data())

         -- Only when the contributions to m0Star from the first direction
         -- are collected, do we collect contributions to m1Star and m2Star.
         -- Also, since Gk velocities are organized as (vpar,mu) the velocity
         -- correction is only computed for the first velocity direction.
         local firstDir = true

         -- isLo=true current cell is the lower boundary cell.
         -- isLo=false current cell is the upper boundary cell.
         local isLo = true

         -- Boundary corrections from this species.
         for vDir = 1, vDim do
            if (not self._isGkLBO) or (self._isGkLBO and firstDir) then
               self._uCorrection    = PrimMomentsDecl.selectBoundaryFintegral(vDir, self._kinSpecies, self._basisID, cDim, vDim, self._polyOrder)
            end
            self._vtSqCorrection = PrimMomentsDecl.selectBoundaryVFintegral(vDir, self._kinSpecies, self._basisID, cDim, vDim, self._polyOrder)

            -- Lower/upper bounds in direction 'vDir': cell indices.
            local dirLoIdx, dirUpIdx = phaseRangeSelf:lower(cDim+vDir), phaseRangeSelf:upper(cDim+vDir)

            if self._isFirst then
               self._perpRangeSelf[vDir] = phaseRangeSelf
               for cd = 1, cDim do
                  self._perpRangeSelf[vDir] = self._perpRangeSelf[vDir]:shorten(cd) -- shorten configuration range.
               end
               self._perpRangeSelf[vDir] = self._perpRangeSelf[vDir]:shorten(cDim+vDir) -- velocity range orthogonal to 'vDir'.
            end
            local perpRange = self._perpRangeSelf[vDir]

            for vPerpIdx in perpRange:rowMajorIter() do
               vPerpIdx:copyInto(idxP)
               for d = 1, cDim do idxP[d] = confIdx[d] end

               for _, i in ipairs({dirLoIdx, dirUpIdx}) do     -- This loop is over edges.
                  idxP[cDim+vDir] = i

                  self._phaseGridSelf:setIndex(idxP)
                  self._phaseGridSelf:getDx(dxP)
                  self._phaseGridSelf:cellCenter(xcP)

                  fInSelf:fill(phaseIndexerSelf(idxP), fInSelfItrP)

                  local vBound = 0.0
                  if isLo then
                     vBound = self._phaseGridSelf:cellLowerInDir(cDim + vDir)
                  else
                     vBound = self._phaseGridSelf:cellUpperInDir(cDim + vDir)
                  end

                  if (self._isGkLBO) then
                     if (firstDir) then
                       self._uCorrection(isLo, self._intFac[1], vBound, dxP:data(), fInSelfItrP:data(), cMomBSelf:data())
                     end
                     self._vtSqCorrection(isLo, self._intFac[vDir], vBound, dxP:data(), fInSelfItrP:data(), cEnergyBSelf:data())
                  else
                     self._uCorrection(isLo, vBound, dxP:data(), fInSelfItrP:data(), cMomBSelf:data())
                     self._vtSqCorrection(isLo, vBound, dxP:data(), fInSelfItrP:data(), cEnergyBSelf:data())
                  end

                  isLo = not isLo
               end
            end
            firstDir = false
         end

         -- Boundary corrections from the other species.
         for vDir = 1, vDim do
            if (not self._isGkLBO) or (self._isGkLBO and firstDir) then
               self._uCorrection    = PrimMomentsDecl.selectBoundaryFintegral(vDir, self._kinSpecies, self._basisID, cDim, vDim, self._polyOrder)
            end
            self._vtSqCorrection = PrimMomentsDecl.selectBoundaryVFintegral(vDir, self._kinSpecies, self._basisID, cDim, vDim, self._polyOrder)

            -- Lower/upper bounds in direction 'vDir': cell indices.
            local dirLoIdx, dirUpIdx = phaseRangeOther:lower(cDim+vDir), phaseRangeOther:upper(cDim+vDir)

            if self._isFirst then
               self._perpRangeOther[vDir] = phaseRangeOther
               for cd = 1, cDim do
                  self._perpRangeOther[vDir] = self._perpRangeOther[vDir]:shorten(cd) -- shorten configuration range.
               end
               self._perpRangeOther[vDir] = self._perpRangeOther[vDir]:shorten(cDim+vDir) -- velocity range orthogonal to 'vDir'.
            end
            local perpRange = self._perpRangeOther[vDir]

            for vPerpIdx in perpRange:rowMajorIter() do
               vPerpIdx:copyInto(idxP)
               for d = 1, cDim do idxP[d] = confIdx[d] end

               for _, i in ipairs({dirLoIdx, dirUpIdx}) do     -- This loop is over edges.
                  idxP[cDim+vDir] = i

                  self._phaseGridOther:setIndex(idxP)
                  self._phaseGridOther:getDx(dxP)
                  self._phaseGridOther:cellCenter(xcP)

                  fInOther:fill(phaseIndexerOther(idxP), fInOtherItrP)

                  local vBound = 0.0
                  if isLo then
                     vBound = self._phaseGridOther:cellLowerInDir(cDim + vDir)
                  else
                     vBound = self._phaseGridOther:cellUpperInDir(cDim + vDir)
                  end

                  if (self._isGkLBO) then
                     if (firstDir) then
                       self._uCorrection(isLo, self._intFac[1], vBound, dxP:data(), fInOtherItrP:data(), cMomBOther:data())
                     end
                     self._vtSqCorrection(isLo, self._intFac[vDir], vBound, dxP:data(), fInOtherItrP:data(), cEnergyBOther:data())
                  else
                     self._uCorrection(isLo, vBound, dxP:data(), fInOtherItrP:data(), cMomBOther:data())
                     self._vtSqCorrection(isLo, vBound, dxP:data(), fInOtherItrP:data(), cEnergyBOther:data())
                  end

                  isLo = not isLo
               end
            end
            firstDir = false
         end

         self._SelfPrimMomentsCalc(self._binOpDataSelf, m0OtherItr:data(), m1OtherItr:data(), m2OtherItr:data(), cMomBOther:data(), cEnergyBOther:data(), uOtherItr:data(), vtSqOtherItr:data())

--         crossPrimMomentsCalc(mRat, self._beta, uSelfItr:data(), vtSqSelfItr:data(), uOtherItr:data(), vtSqOtherItr:data(), uCrossItr:data(), vtSqCrossItr:data())
         crossPrimMomentsCalc(self._binOpData, self._betaP1, mSelf, nuSelf, m0SelfItr:data(), m1SelfItr:data(), m2SelfItr:data(), uSelfItr:data(), vtSqSelfItr:data(), cMomBSelf:data(), cEnergyBSelf:data(), mOther, nuOther, m0OtherItr:data(), m1OtherItr:data(), m2OtherItr:data(), uOtherItr:data(), vtSqOtherItr:data(), cMomBOther:data(), cEnergyBOther:data(), uCrossItr:data(), vtSqCrossItr:data())
      end

   end
   self._isFirst = false
end

return CrossPrimMoments
