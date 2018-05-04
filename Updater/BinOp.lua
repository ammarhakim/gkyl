-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the binary operaation BinOp(A,B) of two fields A, B.
-- Currently only supporting the following:
--   1) Weak division (A divided by B. B has to be a scalar function).
--   2) Weak multiplication.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local BinOpDecl = require "Updater.binOpCalcData.BinOpModDecl"

-- function to check if moment name is correct
local function isOpNameGood(nm)
   if nm == "Multiply" or nm == "Divide" then
      return true
   end
   return false
end

-- Moments updater object
local BinOp = Proto(UpdaterBase)

function BinOp:init(tbl)
   BinOp.super.init(self, tbl) -- setup base object

   self._onGrid = assert(
      tbl.onGrid, "Updater.BinOp: Must provide grid object using 'onGrid'")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.BinOp: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis, "Updater.BinOp: Must provide configuration-space basis object using 'confBasis'")

   local op = assert(
      tbl.operation, "Updater.BinOp: Must provide an operation using 'operation'.")

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
   if isOpNameGood(op) then
      self._BinOpCalc = BinOpDecl.selectBinOpCalc(op, id, self._cDim, polyOrder)
   else
      print("BinOp: Requested operation is", op)
      assert(false, "BinOp: Operation must be one of Multiply, Divide")
   end

end

-- advance method
function BinOp:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   -- deno: "denominator field.
   -- nume: "numerator" field.
   -- uOut: ouput field.
   local deno, nume, uOut = inFld[1], inFld[2], outFld[1]

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   -- assuming the localRange is the same for nume and deno.
   local localRange  = deno:localRange()

   local denoIndexer = deno:genIndexer()
   local numeIndexer = nume:genIndexer()
   local uOutIndexer = uOut:genIndexer()

   local denoItr     = deno:get(1)
   local numeItr     = nume:get(1)
   local uOutItr     = uOut:get(1)

   uOut:scale(0.0) -- zero out moments
   -- loop, computing moments in each cell
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)

      deno:fill(denoIndexer(idx), denoItr)
      nume:fill(numeIndexer(idx), numeItr)
      uOut:fill(uOutIndexer(idx), uOutItr)

      self._BinOpCalc(denoItr:data(), numeItr:data(), vDim, uOutItr:data())
   end

   return true, GKYL_MAX_DOUBLE
end

return BinOp
