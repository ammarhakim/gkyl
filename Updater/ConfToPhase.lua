-- Gkyl ------------------------------------------------------------------------
--
-- Updater to project a configuration space field
-- onto a phase space basis
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local ConfToPhaseDecl = require "Updater.confToPhaseData.ConfToPhaseModDecl"

-- updater object
local ConfToPhase = Proto(UpdaterBase)

-- function to check if operation name is correct
local function isOperationNameGood(nm)
   if nm == "accumulate" or nm == "assign" then
      return true
   end
   return false
end

function ConfToPhase:init(tbl)
   ConfToPhase.super.init(self, tbl) -- setup base object

   self._onGrid = assert(
      tbl.onGrid, "Updater.ConfToPhase: Must provide grid object using 'onGrid'")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.ConfToPhase: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis, "Updater.ConfToPhase: Must provide configuration-space basis object using 'confBasis'")

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
   
   local op = assert(
      tbl.operation, "Updater.ConfToPhase: Must provide operation (accumulate or assign) to compute using 'operation'.")

   -- function to compute specified moment
   if isOperationNameGood(op) then
      self._confToPhaseFun = ConfToPhaseDecl.selectConfToPhase(op, id, self._cDim, self._vDim, polyOrder)
   else
      print("ConfToPhase: Requested operation is", op)
      assert(false, "ConfToPhase: Operation must be one of accumulate, assign")
   end
end

-- advance method
function ConfToPhase:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local fact, fconf = inFld[1], inFld[2]
   local fphase = outFld[1]

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local localRange = fphase:localRange()
   local fphaseIndexer = fphase:genIndexer()
   local fconfIndexer = fconf:genIndexer()

   local fphaseItr, fconfItr = fphase:get(1), fconf:get(1)

   -- loop, computing projection in each cell
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)
      
      fphase:fill(fphaseIndexer(idx), fphaseItr)
      fconf:fill(fconfIndexer(idx), fconfItr)
      self._confToPhaseFun(fact, fconfItr:data(), fphaseItr:data())
   end
end

return ConfToPhase
