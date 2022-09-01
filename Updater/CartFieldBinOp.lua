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
local ffi = require "ffi"
local ffiC = ffi.C
require "Lib.ZeroUtil"

ffi.cdef [[
// Type for storing preallocating memory needed in various batch
// operations
typedef struct gkyl_dg_bin_op_mem gkyl_dg_bin_op_mem;

/**
 * Allocate memory for use in bin op (division operator). Free using
 * release method.
 *
 * @param nbatch Batch size
 * @param neqn Number of equations in each batch
 */
gkyl_dg_bin_op_mem* gkyl_dg_bin_op_mem_new(size_t nbatch, size_t neqn);
// Same as above, except for GPUs
gkyl_dg_bin_op_mem *gkyl_dg_bin_op_mem_cu_dev_new(size_t nbatch, size_t neqn);

/**
 * Release memory needed in the bin ops.
 *
 * @param mem Memory to release
 */
void gkyl_dg_bin_op_mem_release(gkyl_dg_bin_op_mem *mem);

/**
 * Same as gkyl_dg_mul_op, except operator is applied only on
 * specified range (sub-range of range containing the DG fields).
 *
 * @param basis Basis functions used in expansions
 * @param c_oop Component of output field in which to store product
 * @param out Output DG field
 * @param c_lop Component of left operand to use in product
 * @param lop Left operand DG field
 * @param c_rop Component of right operand to use in product
 * @param rop Right operand DG field
 * @param range Range to apply multiplication operator
 */
void gkyl_dg_mul_op_range(struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range *range);

/**
 * Compute pout = cop*pop on specified range (sub-range of range
 * containing the DG fields), where pout and pop are phase-space
 * operands, and cop is a conf-space operand.
 *
 * @param cbasis Configuration space basis functions used in expansions.
 * @param pbasis Phase space basis functions used in expansions.
 * @param pout Output phase-space DG field.
 * @param cop Conf-space operand DG field.
 * @param pop Phase-space operand DG field.
 * @param crange Conf-space range to apply multiplication operator.
 * @param prange Phase-space range to apply multiplication operator.
 */
void gkyl_dg_mul_conf_phase_op_range(struct gkyl_basis *cbasis,
  struct gkyl_basis *pbasis, struct gkyl_array* pout,
  const struct gkyl_array* cop, const struct gkyl_array* pop,
  struct gkyl_range *crange, struct gkyl_range *prange);

/**
 * Same as gkyl_dg_div_op, except operator is applied only on
 * specified range (sub-range of range containing the DG fields).
 *
 * @param mem Pre-allocated space for use in the division
 * @param basis Basis functions used in expansions
 * @param c_oop Component of output field in which to store product
 * @param out Output DG field
 * @param c_lop Component of left operand to use in product
 * @param lop Left operand DG field
 * @param c_rop Component of right operand to use in product
 * @param rop Right operand DG field
 * @param range Range to apply multiplication operator
 */
void gkyl_dg_div_op_range(gkyl_dg_bin_op_mem *mem, struct gkyl_basis basis,
  int c_oop, struct gkyl_array* out,
  int c_lop, const struct gkyl_array* lop,
  int c_rop, const struct gkyl_array* rop, struct gkyl_range *range);
]]

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

   self._weakBasis = assert(
      tbl.weakBasis, "Updater.CartFieldBinOp: Must provide the weak basis object using 'weakBasis'.")
   local op = assert(
      tbl.operation, "Updater.CartFieldBinOp: Must provide an operation using 'operation'.")

   self._fieldBasis = tbl.fieldBasis
   self.onGhosts    = xsys.pickBool(tbl.onGhosts, false)
   self._useGPU     = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU)

   local weakBasis, fieldBasis = self._weakBasis, self._fieldBasis

   -- Positivity option disabled for now. Need to investigate bugs in kernels.
   local applyPositivity = false -- xsys.pickBool(tbl.positivity,false)   -- Positivity preserving option.

   -- Dimension of spaces.
   self._wDim = weakBasis:ndim()
   if fieldBasis then
     -- Dealing with phase space simulation.
     -- Ensure sanity.
     assert(weakBasis:polyOrder() == fieldBasis:polyOrder(),
            "Polynomial orders of weak and field basis must match.")
     assert((weakBasis:id() == fieldBasis:id()) or
            ((weakBasis:id()=="hybrid" or weakBasis:id()=="gkhybrid") and fieldBasis:id()=="serendipity"),
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
   if (id=="hybrid" or id=="gkhybrid") then id="serendipity" end

   -- Function to compute specified operation.
   if isOpNameGood(op) then
      self._BinOpCalcS = BinOpDecl.selectBinOpCalcS(op, id, self._cDim, self._vDim, polyOrder, applyPositivity)
      if fieldBasis then self._BinOpCalcD = BinOpDecl.selectBinOpCalcD(op, id, self._cDim, self._vDim, polyOrder) end
   else
      assert(false, string.format(
		"CartFieldBinOp: Operation must be one of Multiply, Divide, DotProduct. Requested %s instead.", op))
   end

   -- Create struct containing allocated binOp arrays.
   if fieldBasis then 
      self._binOpData = ffiC.new_binOpData_t(fieldBasis:numBasis(), self._numBasis) 
   else 
      self._binOpData = ffiC.new_binOpData_t(self._numBasis, 0) 
   end

   if op == "Multiply" or op == "Divide" then self._zero_op = op end

   if op == "Divide" then
      local onRange = assert(tbl.onRange, "Updater.CartFieldBinOp: Must provide the range to perform division in 'onRange'.")
      if self._useGPU then
         self._mem = ffi.gc(ffiC.gkyl_dg_bin_op_mem_cu_dev_new(onRange:volume(), self._numBasis),
                            ffiC.gkyl_dg_bin_op_mem_release)
      else
         self._mem = ffi.gc(ffiC.gkyl_dg_bin_op_mem_new(onRange:volume(), self._numBasis),
                            ffiC.gkyl_dg_bin_op_mem_release)
      end
   end
end

-- Advance method.
function CartFieldBinOp:_advance(tCurr, inFld, outFld)
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

   local localuRange = self.onGhosts and uOut:localExtRange() or uOut:localRange()
   if self._zero_op and self._zero_op == "Multiply" then
      if self._fieldBasis then
         -- Conf-space * phase-space multiplication.
         local localARange = self.onGhosts and Afld:localExtRange() or Afld:localRange()
         ffiC.gkyl_dg_mul_conf_phase_op_range(self._fieldBasis._zero, self._weakBasis._zero,
                                              uOut._zero, Afld._zero, Bfld._zero, localARange, localuRange)
      else
         -- Conf-space scalar * scalar or scalar * vector multiplication.
         local nVecComp = Bfld:numComponents()/self._numBasis
         for d = 0, nVecComp-1 do
            ffiC.gkyl_dg_mul_op_range(self._weakBasis._zero, d, uOut._zero, 0, Afld._zero, d, Bfld._zero, localuRange)
         end
      end
   
      return
   elseif self._zero_op and self._zero_op == "Divide" then
      -- Conf-space scalar / scalar or vector / scalar division.
      local nVecComp = Bfld:numComponents()/self._numBasis
      for d = 0, nVecComp-1 do
         ffiC.gkyl_dg_div_op_range(self._mem, self._weakBasis._zero, d, uOut._zero, d, Bfld._zero, 0, Afld._zero, localuRange)
      end
      return
   end

   -- g2 implementation below. To be deleted eventually.

   local grid = Afld:grid()

   -- Either the localRange is the same for Bfld and Afld,
   -- or just use the range of the phase space field,
   local localBRange = self.onGhosts and Bfld:localExtRange() or Bfld:localRange()

   local AfldIndexer = Afld:genIndexer()
   local BfldIndexer = Bfld:genIndexer()
   local uOutIndexer = uOut:genIndexer()

   local AfldItr = Afld:get(1)
   local BfldItr = Bfld:get(1)
   local uOutItr = uOut:get(1)

   -- Number of vectorial components.
   local nComp = Bfld:numComponents()/self._numBasis

   -- This factor is used in the kernel to differentiate the case of
   -- scalar-vector multiplication from the vector-vector multiplication.
   local nCompEq = Afld:numComponents() == Bfld:numComponents() and 1 or 0
   -- If the dimensions are the same assume we are doing vector-vector or
   -- scalar-vector operations in configuration space.
   -- If dimensions are different assume we are doing multiplication of
   -- a conf space function (scalar?) with a phase space function.
   self._BinOpCalc = (Afld:ndim() == Bfld:ndim()) and self._BinOpCalcS or self._BinOpCalcD

   -- Loop, computing binOp in each cell.
   for idx in localBRange:rowMajorIter() do
      grid:setIndex(idx)

      Afld:fill(AfldIndexer(idx), AfldItr)
      Bfld:fill(BfldIndexer(idx), BfldItr)
      uOut:fill(uOutIndexer(idx), uOutItr)

      self._BinOpCalc(self._binOpData, AfldItr:data(), BfldItr:data(), nComp, nCompEq, uOutItr:data())
   end
end

function CartFieldBinOp:_advanceOnDevice(tCurr, inFld, outFld)
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

   local localuRange = self.onGhosts and uOut:localExtRange() or uOut:localRange()
   if self._zero_op and self._zero_op == "Multiply" then
      if self._fieldBasis then
         -- Conf-space * phase-space multiplication.
         local localARange = self.onGhosts and Afld:localExtRange() or Afld:localRange()
         ffiC.gkyl_dg_mul_conf_phase_op_range(self._fieldBasis._zero, self._weakBasis._zero,
                                              uOut._zeroDevice, Afld._zeroDevice, Bfld._zeroDevice, localARange, localuRange)
      else
         -- Conf-space scalar * scalar or scalar * vector multiplication.
         local nVecComp = Bfld:numComponents()/self._numBasis
         for d = 0, nVecComp-1 do
            ffiC.gkyl_dg_mul_op_range(self._weakBasis._zero, d, uOut._zeroDevice,
	                              0, Afld._zeroDevice, d, Bfld._zeroDevice, localuRange)
         end
      end
   
      return
   elseif self._zero_op and self._zero_op == "Divide" then
      -- Conf-space scalar / scalar or vector / scalar division.
      local nVecComp = Bfld:numComponents()/self._numBasis
      for d = 0, nVecComp-1 do
         ffiC.gkyl_dg_div_op_range(self._mem, self._weakBasis._zero, d, uOut._zeroDevice,
	                           d, Bfld._zeroDevice, 0, Afld._zeroDevice, localuRange)
      end
      return
   else 
      -- NYI
      assert(false, "GPU bin op NYI")

      return
   end

end

return CartFieldBinOp
