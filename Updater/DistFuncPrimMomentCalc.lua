-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute primitive moments given the standard moments computed
-- by kinetic solvers. For example, compute the flow speed 'u' given the number
-- density 'numDen' and the first moment 'mom1=numDen*u'.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local PrimMomentDecl = require "Updater.primMomentCalcData.DistFuncPrimMomentCalcModDecl"

-- function to check if moment name is correct
local function isMomentNameGood(nm)
   if nm == "Ui" or nm == "VtSq" then
      return true
   end
   return false
end
local function isGkMomentNameGood(nm)
   if nm == "GkUpar" or nm == "GkVtSq" then
      return true
   end
   return false
end

-- Moments updater object
local DistFuncPrimMomentCalc = Proto(UpdaterBase)

function DistFuncPrimMomentCalc:init(tbl)
   DistFuncPrimMomentCalc.super.init(self, tbl) -- setup base object

   self._onGrid = assert(
      tbl.onGrid, "Updater.DistFuncPrimMomentCalc: Must provide grid object using 'onGrid'")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.DistFuncPrimMomentCalc: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis, "Updater.DistFuncPrimMomentCalc: Must provide configuration-space basis object using 'confBasis'")

   local mom = assert(
      tbl.moment, "Updater.DistFuncPrimMomentCalc: Must provide moment to compute using 'moment'.")

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

   -- function to compute specified moment
   if isMomentNameGood(mom) then
      self._PmomCalc = PrimMomentDecl.selectMomCalc(mom, id, self._cDim, self._vDim, polyOrder)
   elseif isGkMomentNameGood(mom) then
      self._PmomCalc = PrimMomentDecl.selectGkMomCalc(mom, id, self._cDim, self._vDim, polyOrder)
   else
      print("DistFuncPrimMomentCalc: Requested moment is", mom)
      assert(false, "DistFuncPrimMomentCalc: Moments must be one of Ui, VtSq")
   end

end

-- advance method
function DistFuncPrimMomentCalc:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local mom0, mom1, mom2, Pmom = inFld[1], inFld[2], inFld[3], outFld[1]

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   -- assuming the localRange is the same for mom0, mom1, mom2.
   local localRange  = mom0:localRange()

   local mom0Indexer = mom0:genIndexer()
   local mom1Indexer = mom1:genIndexer()
   local mom2Indexer = mom2:genIndexer()
   local PmomIndexer = Pmom:genIndexer()

   local mom0Itr     = mom0:get(1)
   local mom1Itr     = mom1:get(1)
   local mom2Itr     = mom2:get(1)
   local PmomItr     = Pmom:get(1)

   Pmom:scale(0.0) -- zero out moments
   -- loop, computing moments in each cell
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)

      mom0:fill(mom0Indexer(idx), mom0Itr)
      mom1:fill(mom1Indexer(idx), mom1Itr)
      mom2:fill(mom2Indexer(idx), mom2Itr)

      Pmom:fill(PmomIndexer(idx), PmomItr)

      self._PmomCalc(mom0Itr:data(), mom1Itr:data(), mom2Itr:data(), PmomItr:data())
   end

   return true, GKYL_MAX_DOUBLE
end

return DistFuncPrimMomentCalc
