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

   self._phaseGrid = assert(
      tbl.phaseGrid, "Updater.SelfPrimMoments: Must provide phase grid object using 'phaseGrid'.")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.SelfPrimMoments: Must provide the phase basis object using 'phaseBasis'.")

   local confBasis = assert(
      tbl.confBasis, "Updater.SelfPrimMoments: Must provide the configuration basis object using 'confBasis'.")

   local operator = assert(
      tbl.operator, "Updater.SelfPrimMoments: Must specify the collision operator (VmLBO or GkLBO for now) using 'operator'.")

   -- dimension of phase space.
   self._pDim = phaseBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID   = confBasis:id()
   self._polyOrder = confBasis:polyOrder()

   -- ensure sanity.
   assert(phaseBasis:polyOrder() == self._polyOrder,
          "Polynomial orders of phase and conf basis must match.")
   assert(phaseBasis:id() == self._basisID,
          "Type of phase and conf basis must match.")
   -- determine configuration and velocity space dims.
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Number of basis functions. Used to compute number of vector components.
   self._numBasisP = phaseBasis:numBasis()
   self._numBasisC = confBasis:numBasis()

   -- Extract vLower and vUpper for corrections to u and vtSq
   -- that conserve momentum and energy.
   self._vLower = Lin.Vec(self._vDim)
   self._vUpper = Lin.Vec(self._vDim)
   for d = 1, self._vDim do
      -- Pass absolute value of vLower because kernel already
      -- evaluates v*f at -1 (see Maxima scripts).
      self._vLower[d] = self._phaseGrid:lower(self._cDim + d)
      self._vUpper[d] = self._phaseGrid:upper(self._cDim + d) 
   end

   -- If vDim>1, intFac=2*pi/m or 4*pi/m for Gk. This prevents 
   -- the need for separate kernels for Gk.
   self._intFac       = Lin.Vec(self._vDim)
   for d = 1,self._vDim do 
     self._intFac[d]  = 1.0
   end
   self._physVdim     = self._vDim
   local vCorrections = self._vDim
   if isOperatorGood(operator) then
     if operator == "GkLBO" then
        self._isGkLBO   = true
        local mass = assert(
           tbl.mass, "Updater.SelfPrimMoments: Must provide species mass using 'mass'.")
        if self._vDim > 1 then -- A (vpar,mu) simulation has 3 physical velocity dimensions.
           self._intFac[1] = 2.0*math.pi/mass
           self._intFac[2] = 4.0*math.pi/mass
           self._physVdim  = 3
           vCorrections    = 1
        end
     else -- VmLBO
        self._isGkLBO = false
     end
   else
      assert(false, string.format(
                "SelfPrimMoments: Operator option must be 'VmLBO' or 'GkLBO'. Requested %s instead.", operator))
   end

   self._SelfPrimMomentsCalc = PrimMomentsDecl.selectSelfPrimMomentsCalc(self._basisID, self._cDim, vCorrections, self._polyOrder)

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)
   
   self._isFirst = true
   self._perpRange = {} -- perp ranges in velocity directions.
   if self._polyOrder == 1 then
     self._perpRangeR = {} -- restricted perp ranges in velocity directions.
   end

end

-- advance method
function SelfPrimMoments:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid

   local uOut    = outFld[1]
   local vtSqOut = outFld[2]

   local m0fld, m1fld, m2fld, fIn, m1Sfld, m2Sfld = inFld[1], inFld[2], inFld[3], inFld[4], inFld[5], inFld[6]

   local confRange  = m0fld:localRange()
   if self.onGhosts then confRange = m0fld:localExtRange() end

   local m0fldIndexer   = m0fld:genIndexer()
   local m1fldIndexer   = m1fld:genIndexer()
   local m2fldIndexer   = m2fld:genIndexer()
   local uOutIndexer    = uOut:genIndexer()
   local vtSqOutIndexer = vtSqOut:genIndexer()
   local phaseIndexer   = fIn:genIndexer()
   local m1SfldIndexer  = m1Sfld:genIndexer()
   local m2SfldIndexer  = m2Sfld:genIndexer()

   local m0fldItr   = m0fld:get(1)
   local m1fldItr   = m1fld:get(1)
   local m2fldItr   = m2fld:get(1)
   local uOutItr    = uOut:get(1)
   local vtSqOutItr = vtSqOut:get(1)
   local fInItr     = fIn:get(1)
   local m1SfldItr  = m1Sfld:get(1)

   -- obtain f at edges of velocity domain.
   local phaseRange           = fIn:localRange()
   local vLowerIdx, vUpperIdx = {}, {}
   for d = 1, self._vDim do
      vLowerIdx[d] = phaseRange:lower(self._cDim + d)
      vUpperIdx[d] = phaseRange:upper(self._cDim + d)
   end

   -- Distribution functions left and right of a cell-boundary.
   local fInP = Lin.Vec(self._numBasisP)
   local fInM = Lin.Vec(self._numBasisP)
   -- Distribution functions at boundaries of velocity space.
   local fVmin = Lin.Vec(self._numBasisP)
   local fVmax = Lin.Vec(self._numBasisP)

   -- These store the corrections to momentum and energy
   -- from the (valocity) boundary surface integrals.
   local cMomB    = Lin.Vec(self._numBasisC*self._vDim)
   local cEnergyB = Lin.Vec(self._numBasisC)
   local m0Star   = Lin.Vec(self._numBasisC)

   local phaseIdx = Lin.IntVec(self._pDim)
   local dxv      = Lin.Vec(self._pDim) -- cell shape.
   -- Cell index, center and length left and right of a cell-boundary.
   local idxM, idxP = Lin.IntVec(self._pDim), Lin.IntVec(self._pDim)
   local xcM, xcP   = Lin.Vec(self._pDim), Lin.Vec(self._pDim)
   local dxM, dxP   = Lin.Vec(self._pDim), Lin.Vec(self._pDim)

   for confIdx in confRange:colMajorIter() do
      grid:setIndex(confIdx)

      m0fld:fill(m0fldIndexer(confIdx), m0fldItr)
      m1fld:fill(m1fldIndexer(confIdx), m1fldItr)
      uOut:fill(uOutIndexer(confIdx), uOutItr)
      vtSqOut:fill(vtSqOutIndexer(confIdx), vtSqOutItr)
      if self._polyOrder == 1 then
         -- Piecewise linear uses star moments.
         m1Sfld:fill(m1SfldIndexer(confIdx), m1SfldItr)
         m2Sfld:fill(m2SfldIndexer(confIdx), m2fldItr)
      else
         -- p>1 uses regular moments only.
         m1fld:fill(m1fldIndexer(confIdx), m1SfldItr)
         m2fld:fill(m2fldIndexer(confIdx), m2fldItr)
      end

      -- Compute the corrections to u and vtSq due to 
      -- finite velocity grid to ensure conservation.
      for k = 1, self._numBasisC do
        cEnergyB[k] = 0 
        for vd = 1, self._vDim do
          cMomB[(vd-1)*self._numBasisC + k] = 0
        end
         -- To have energy conservation with piece-wise linear, we must use
         -- star moments in the second equation of the weak system solved
         -- in SelfPrimMoments.
         m0Star[k] = 0 
      end

      if self._polyOrder == 1 then
         -- To have energy conservation with piece-wise linear, we must use
         -- star moments in the second equation of the weak system solved
         -- in SelfPrimMoments.
         for vDir = 1, self._vDim do
            self._starMom = PrimMomentsDecl.selectStarMoms(vDir, self._basisID, self._cDim, self._vDim)

            -- Lower/upper bounds in direction 'vDir': inner edge indices.
            local dirLoIdx, dirUpIdx = phaseRange:lower(self._cDim+vDir)+1, phaseRange:upper(self._cDim+vDir)

            if self._isFirst then
               -- Restricted velocity range.
               -- Velocity integral in m0Star does not include last cell.
               self._perpRangeR[vDir] = phaseRange
               for cd = 1, self._cDim do
                  self._perpRangeR[vDir] = self._perpRangeR[vDir]:shorten(cd) -- shorten configuration range.
               end
               self._perpRangeR[vDir] = self._perpRangeR[vDir]:shorten(self._cDim+vDir) -- velocity range orthogonal to 'vDir'.
            end
            local perpRangeR = self._perpRangeR[vDir]

            -- Outer loop is over directions orthogonal to 'vDir' and
            -- inner loop is over 1D slice in `vDir`.
            for idx in perpRangeR:colMajorIter() do
               idx:copyInto(idxM); idx:copyInto(idxP)
               for d = 1, self._cDim do idxM[d] = confIdx[d] end
               for d = 1, self._cDim do idxP[d] = confIdx[d] end

               for i = dirLoIdx, dirUpIdx do     -- This loop is over edges.
                  idxM[self._cDim+vDir], idxP[self._cDim+vDir]  = i-1, i -- Cell left/right of edge 'i'.

	          self._phaseGrid:setIndex(idxM)
	          for d = 1, self._pDim do dxM[d] = self._phaseGrid:dx(d) end
	          self._phaseGrid:cellCenter(xcM)

	          self._phaseGrid:setIndex(idxP)
	          for d = 1, self._pDim do dxP[d] = self._phaseGrid:dx(d) end
	          self._phaseGrid:cellCenter(xcP)

                  fIn:fill(phaseIndexer(idxM), fInItr)
                  for k = 1, self._numBasisP do fInM[k] = fInItr[k] end

                  fIn:fill(phaseIndexer(idxP), fInItr)
                  for k = 1, self._numBasisP do fInP[k] = fInItr[k] end

                  self._starMom(self._intFac[vDir], xcM:data(), xcP:data(), dxM:data(), dxP:data(), fInM:data(), fInP:data(), m0Star:data())
               end
            end
         end
      end

      for vDir = 1, self._vDim do
         self._uCorrection = PrimMomentsDecl.selectBoundaryFintegral(vDir, self._basisID, self._cDim, self._vDim, self._polyOrder)
         self._vtSqCorrection = PrimMomentsDecl.selectBoundaryVFintegral(vDir, self._basisID, self._cDim, self._vDim, self._polyOrder)

         if self._isFirst then
            self._perpRange[vDir] = phaseRange
            for cd = 1, self._cDim do
               self._perpRange[vDir] = self._perpRange[vDir]:shorten(cd) -- shorten configuration range.
            end
            self._perpRange[vDir] = self._perpRange[vDir]:shorten(self._cDim+vDir) -- velocity range orthogonal to 'vDir'.
         end
         local perpRange = self._perpRange[vDir]

         for vPerpIdx in perpRange:colMajorIter() do
            vPerpIdx:copyInto(phaseIdx)
            for d = 1, self._cDim do phaseIdx[d] = confIdx[d] end

            phaseIdx[self._cDim + vDir] = vLowerIdx[vDir]
            fIn:fill(phaseIndexer(phaseIdx), fInItr)
            for k = 1, self._numBasisP do fVmin[k] = fInItr[k] end

            phaseIdx[self._cDim + vDir] = vUpperIdx[vDir]
            fIn:fill(phaseIndexer(phaseIdx), fInItr)
            for k = 1, self._numBasisP do fVmax[k] = fInItr[k] end

            -- BEWARE: This assumes that dxv is the same at vmin and vmax. So do correction kernels.
            self._phaseGrid:setIndex(phaseIdx)
            for d = 1, self._pDim do dxv[d] = self._phaseGrid:dx(d) end

            if (not self._isGkLBO) or (self._isGkLBO and vDir==1) then
               self._uCorrection(self._intFac[vDir], self._vLower[vDir], self._vUpper[vDir], dxv:data(), fVmin:data(), fVmax:data(), cMomB:data())
            end

            self._vtSqCorrection(self._intFac[vDir], self._vLower[vDir], self._vUpper[vDir], dxv:data(), fVmin:data(), fVmax:data(), cEnergyB:data())
         end
      end

      self._SelfPrimMomentsCalc(self._physVdim, m0fldItr:data(), m1fldItr:data(), m2fldItr:data(), m0Star:data(), m1SfldItr:data(), cMomB:data(), cEnergyB:data(), uOutItr:data(), vtSqOutItr:data())

   end
   self._isFirst = false

   return true, GKYL_MAX_DOUBLE
end

return SelfPrimMoments
