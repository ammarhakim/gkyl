-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments of distribution function on a
-- rectangular (but potentially non-uniform) grid.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"
local Proto = require "Lib.Proto"
local MomDecl = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"

-- function to check if moment name is correct
local function isMomentNameGood(nm)
   if nm == "M0" or nm == "M1i" or nm == "M2ij" or nm == "M2" or nm == "M3i" then
      return true
   end
   return false
end

-- Moments updater object
local DistFuncMomentCalc = Proto(UpdaterBase)

function DistFuncMomentCalc:init(tbl)
   DistFuncMomentCalc.super.init(self, tbl) -- setup base object

   self._onGrid = assert(
      tbl.onGrid, "Updater.DistFuncMomentCalc: Must provide grid object using 'onGrid'")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.DistFuncMomentCalc: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis, "Updater.DistFuncMomentCalc: Must provide configuration-space basis object using 'confBasis'")

   -- dimension of spaces
   self._pDim = phaseBasis:ndim() 
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- ensure sanity
   assert(phaseBasis:polyOrder() == confBasis:polyOrder(),
	  "Polynomial orders of phase-space and config-space basis must match")
   assert(phaseBasis:id() == confBasis:id(),
	  "Type of phase-space and config-space basis must match")

   local id, polyOrder = phaseBasis:id(), phaseBasis:polyOrder()
   
   local mom = assert(
      tbl.moment, "Updater.DistFuncMomentCalc: Must provide moment to compute using 'moment'.")

   -- function to compute specified moment
   if isMomentNameGood(mom) then
      self._momCalcFun = MomDecl.selectMomCalc(mom, id, self._cDim, self._vDim, polyOrder)
   else
      assert(false, "DistFuncMomentCalc: Moments must be one of M0, M1i, M2ij, M2, M3i")
   end

   -- for use in _advance() method
   self.dxv = Lin.Vec(self._pDim) -- cell shape
   self.w = Lin.Vec(self._pDim) -- phase-space cell center
end

-- advance method
function DistFuncMomentCalc:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local distf = assert(inFld[1], "DistFuncMomentCalc:advance: Must supply input distribution function")
   local mom = assert(outFld[1], "DistFuncMomentCalc:advance: Must supply output moment field")

   -- sanity checks
   assert(distf:ndim() == self._pDim,
   	  "DistFuncMomentCalc:advance: Phase-space dimensions do not match input distribution function ")
   assert(mom:ndim() == self._cDim,
   	  "DistFuncMomentCalc:advance: Config-space dimensions do not match output moment field")

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local localRange = distf:localRange()
   local distfIndexer = distf:genIndexer()
   local momIndexer = mom:genIndexer()

   local distfItr, momItr = distf:get(1), mom:get(1)

   mom:scale(0.0) -- zero out moments
   -- loop, computing moments in each cell
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)
      grid:cellCenter(self.w)
      for d = 1, pDim do self.dxv[d] = grid:dx(d) end
      
      distf:fill(distfIndexer(idx), distfItr)
      mom:fill(momIndexer(idx), momItr)
      self._momCalcFun(self.w:data(), self.dxv:data(), distfItr:data(), momItr:data())
   end
   
   return true, GKYL_MAX_DOUBLE
end

return DistFuncMomentCalc
