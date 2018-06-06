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

-- Moments updater object.
local SelfPrimMoments = Proto(UpdaterBase)

function SelfPrimMoments:init(tbl)
   SelfPrimMoments.super.init(self, tbl) -- setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.SelfPrimMoments: Must provide grid object using 'onGrid'.")

   local phaseGrid = assert(
      tbl.phaseGrid, "Updater.SelfPrimMoments: Must provide phase grid object using 'phaseGrid'.")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.SelfPrimMoments: Must provide the phase basis object using 'phaseBasis'.")

   local confBasis = assert(
      tbl.confBasis, "Updater.SelfPrimMoments: Must provide the configuration basis object using 'confBasis'.")

   -- dimension of spaces.
   self._pDim = phaseBasis:ndim()
   -- ensure sanity.
   assert(phaseBasis:polyOrder() == confBasis:polyOrder(),
          "Polynomial orders of phase and conf basis must match.")
   assert(phaseBasis:id() == confBasis:id(),
          "Type of phase and conf basis must match.")
   -- determine configuration and velocity space dims.
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Number of basis functions. Used to compute number of vector components.
   self._numBasisP = phaseBasis:numBasis()
   self._numBasisC = confBasis:numBasis()

   local id, polyOrder = confBasis:id(), confBasis:polyOrder()

   -- Extract vLower and vUpper for corrections to u and vtSq
   -- that preserve momentum and energy conservation.
   self._vLower = Lin.Vec(self._vDim)
   self._vUpper = Lin.Vec(self._vDim)
   for d = 1, self._vDim do
      -- Pass absolute value of vLower because kernel already
      -- evaluates v*f at -1 (see Maxima scripts).
      self._vLower[d] = math.abs(phaseGrid:lower(self._cDim + d))
      self._vUpper[d] = phaseGrid:upper(self._cDim + d) 
   end

   self._SelfPrimMomentsCalc = PrimMomentsDecl.selectSelfPrimMomentsCalc(id, self._cDim, self._vDim, polyOrder)

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)

end

-- advance method
function SelfPrimMoments:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid

   local uOut    = outFld[1]
   local vtSqOut = outFld[2]

   local m0fld, m1fld, m2fld, fIn
   m0fld, m1fld, m2fld, fIn = inFld[1], inFld[2], inFld[3], inFld[4]

   local confRange  = m0fld:localRange()
   if self.onGhosts then confRange = m0fld:localExtRange() end

   local m0fldIndexer   = m0fld:genIndexer()
   local m1fldIndexer   = m1fld:genIndexer()
   local m2fldIndexer   = m2fld:genIndexer()
   local uOutIndexer    = uOut:genIndexer()
   local vtSqOutIndexer = vtSqOut:genIndexer()
   local phaseIndexer   = fIn:genIndexer()

   local m0fldItr   = m0fld:get(1)
   local m1fldItr   = m1fld:get(1)
   local m2fldItr   = m2fld:get(1)
   local uOutItr    = uOut:get(1)
   local vtSqOutItr = vtSqOut:get(1)
   local fInItr     = fIn:get(1)

   -- obtain f at edges of velocity domain.
   local phaseIdx     = {}
   local phaseRange   = fIn:localRange()
   local vLowerIdx, vUpperIdx = {}, {}
   for d = 1, self._vDim do
      vLowerIdx[d]  = phaseRange:lower(self._cDim + d)
      vUpperIdxV[d] = phaseRange:upper(self._cDim + d)
   end

   local fvmin = Lin.Vec(self._numBasisP*self._vDim)
   local fvmax = Lin.Vec(self._numBasisP*self._vDim)

   -- loop, computing moments in each cell
   for confIdx in confRange:colMajorIter() do
      grid:setIndex(confIdx)

      m0fld:fill(m0fldIndexer(confIdx), m0fldItr)
      m1fld:fill(m1fldIndexer(confIdx), m1fldItr)
      m2fld:fill(m2fldIndexer(confIdx), m2fldItr)
      uOut:fill(uOutIndexer(confIdx), uOutItr)
      vtSqOut:fill(vtSqOutIndexer(confIdx), vtSqOutItr)

      -- Need last cells (along velocity dimensions) of distribution
      -- function for LBO to be momentum and energy conserving.
      for d = 1, self._cDim do phaseIdx[d] = confIdx[d] end

      for d = 1, self._vDim do phaseIdx[d + self._cDim] = vLowerIdx[d] end
      fIn:fill(phaseIndexer(phaseIdx), fInItr)
      for k = 1, self._numBasisP do fvmin[k] = fInItr[k] end

      for d = 1, self._vDim do phaseIdx[d + self._cDim] = vUpperIdx[d] end
      fIn:fill(phaseIndexer(phaseIdx), fInItr)
      for k = 1, self._numBasisP do fvmax[k] = fInItr[k] end

      self._SelfPrimMomentsCalc(m0fldItr:data(), m1fldItr:data(), m2fldItr:data(), fvmin:data(), fvmax:data(), self._vLower:data(), self._vUpper:data(), uOutItr:data(), vtSqOutItr:data())
   end

   return true, GKYL_MAX_DOUBLE
end

return SelfPrimMoments
