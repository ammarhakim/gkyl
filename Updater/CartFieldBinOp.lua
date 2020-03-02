-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the binary operaation BinOp(A,B) of two fields A, B.
-- Currently only supporting the following:
--   1) Weak division (B divided by A. A must be a scalar function).
--   2) Weak multiplication.
--   3) Dot product of two vector fields.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local BinOpDecl    = require "Updater.binOpCalcData.CartFieldBinOpModDecl"
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto        = require "Lib.Proto"
local UpdaterBase  = require "Updater.Base"
local xsys         = require "xsys"
local ffi          = require "ffi"
local ffiC         = ffi.C

-- Function to check if moment name is correct.
local function isOpNameGood(nm)
   if nm == "Multiply" or nm == "Divide" or nm=="DotProduct" then
      return true
   end
   return false
end

-- Moments updater object.
local CartFieldBinOp = Proto(UpdaterBase)

function CartFieldBinOp:init(tbl)
   CartFieldBinOp.super.init(self, tbl) -- Setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.CartFieldBinOp: Must provide grid object using 'onGrid'.")

   local weakBasis = assert(
      tbl.weakBasis, "Updater.CartFieldBinOp: Must provide the weak basis object using 'weakBasis'.")
   local fieldBasis = tbl.fieldBasis

   local op = assert(
      tbl.operation, "Updater.CartFieldBinOp: Must provide an operation using 'operation'.")
 
   local inclVx2 = xsys.pickBool(tbl.inclVx2, false)

   local applyPositivity = xsys.pickBool(tbl.positivity,false)   -- Positivity preserving option.

   -- Dimension of spaces.
   self._wDim = weakBasis:ndim()
   if fieldBasis then
     -- Dealing with phase space simulation.
     -- Ensure sanity.
     assert(weakBasis:polyOrder() == fieldBasis:polyOrder(),
            "Polynomial orders of weak and field basis must match.")
     assert(weakBasis:id() == fieldBasis:id(),
            "Type of weak and field basis must match.")
     -- Determine configuration and velocity space dims.
     self._cDim = fieldBasis:ndim()
     self._vDim = self._wDim - self._cDim
   else
     self._cDim = self._wDim
   end

   -- Number of basis functions. Used to compute number of vector components.
   self._numBasis = weakBasis:numBasis()

   local id, polyOrder = weakBasis:id(), weakBasis:polyOrder()

   -- Function to compute specified operation.
   if isOpNameGood(op) then
      self._BinOpCalcS = BinOpDecl.selectBinOpCalcS(op, id, self._cDim, polyOrder, applyPositivity)
      if fieldBasis then self._BinOpCalcD = BinOpDecl.selectBinOpCalcD(op, id, self._cDim, self._vDim, polyOrder, inclVx2) end
   else
      assert(false, string.format(
		"CartFieldBinOp: Operation must be one of Multiply, Divide, DotProduct. Requested %s instead.", op))
   end

   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)

   -- Create struct containing allocated binOp arrays.
   if fieldBasis then 
      self._binOpData = ffiC.new_binOpData_t(fieldBasis:numBasis(), self._numBasis) 
   else 
      self._binOpData = ffiC.new_binOpData_t(self._numBasis, 0) 
   end
end

-- Advance method.
function CartFieldBinOp:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   -- Multiplication: Afld * Bfld (can be scalar*scalar, vector*scalar or scalar*vector,
   --                              but in the latter Afld must be the scalar).
   -- Division:       Bfld/Afld (Afld must be a scalar function).
   -- DotProduct:     Afld . Bfld (both vector fields).

   local uOut = outFld[1]
   -- Remove SOME burden from the user in ordering the inputs. Order them here so that
   -- in scalar-vector and conf-phase operations they enter the kernels as
   -- BinOp(scalar,vector) and BinOp(conf field,phase field).
   local Afld, Bfld
   if inFld[1]:numComponents() <= inFld[2]:numComponents() then
      Afld, Bfld = inFld[1], inFld[2]
   elseif inFld[1]:numComponents() > inFld[2]:numComponents() then
      Bfld, Afld = inFld[1], inFld[2]
   end

   -- Either the localRange is the same for Bfld and Afld,
   -- or just use the range of the phase space field,
   local localBRangeDecomp
   if self.onGhosts then
      localBRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = Bfld:localExtRange(), numSplit = grid:numSharedProcs() }
   else
      localBRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = Bfld:localRange(), numSplit = grid:numSharedProcs() }
   end

   local AfldIndexer = Afld:genIndexer()
   local BfldIndexer = Bfld:genIndexer()
   local uOutIndexer = uOut:genIndexer()

   local AfldItr = Afld:get(1)
   local BfldItr = Bfld:get(1)
   local uOutItr = uOut:get(1)

   -- Number of vectorial components.
   local nComp = Bfld:numComponents()/self._numBasis

   local nCompEq
   -- This factor is used in the kernel to differentiate the case of
   -- scalar-vector multiplication from the vector-vector multiplication.
   if Afld:numComponents() == Bfld:numComponents() then
     nCompEq = 1
   else
     nCompEq = 0
   end
   if (Afld:ndim() == Bfld:ndim()) then
     -- If the dimensions are the same assume we are doing vector-vector or
     -- scalar-vector operations in configuration space.
     self._BinOpCalc = self._BinOpCalcS
   else
     -- If dimensions are different assume we are doing multiplication of
     -- a conf space function (scalar?) with a phase space function.
     self._BinOpCalc = self._BinOpCalcD
   end

   local tId = grid:subGridSharedId() -- Local thread ID.
   -- Loop, computing binOp in each cell.
   for idx in localBRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(idx)

      Afld:fill(AfldIndexer(idx), AfldItr)
      Bfld:fill(BfldIndexer(idx), BfldItr)
      uOut:fill(uOutIndexer(idx), uOutItr)

      self._BinOpCalc(self._binOpData, AfldItr:data(), BfldItr:data(), nComp, nCompEq, uOutItr:data())
   end
end

return CartFieldBinOp
