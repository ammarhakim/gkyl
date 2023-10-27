-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute integrated quantities on a Cartesian grid.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Proto       = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local ffi         = require "ffi"
local xsys        = require "xsys"
local Lin         = require "Lib.Linalg"
local cuda
if GKYL_HAVE_CUDA then
   cuda = require "Cuda.RunTime"
end

local ffiC = ffi.C
ffi.cdef [[
// Object type
typedef struct gkyl_array_integrate gkyl_array_integrate;

enum gkyl_array_integrate_op {
  GKYL_ARRAY_INTEGRATE_OP_NONE = 0,
  GKYL_ARRAY_INTEGRATE_OP_ABS,
  GKYL_ARRAY_INTEGRATE_OP_SQ,
};

/**
 * Create a new updater that integrates a gkyl_array, with an option to perform
 * additional operations during the integration (e.g. square, abs, etc). *
 *
 * @param grid Grid array is defined on.
 * @param basis Basis array is defined on.
 * @param num_comp Number of (vector) components in the array).
 * @param op Additional operator to apply in very cell.
 * @param use_gpu Indicate whether to perform integral on the device.
 */
struct gkyl_array_integrate*
gkyl_array_integrate_new(const struct gkyl_rect_grid* grid, const struct gkyl_basis* basis,
  int num_comp, enum gkyl_array_integrate_op op, bool use_gpu);

/**
 * Compute the array integral.
 *
 * @param up array_integrate updater.
 * @param arr Input gkyl_array.
 * @param weight Factor to multiply by.
 * @param range Range we'll integrate over.
 * @return out Output integral result(s). On device memory if use_gpu=true.
 */
void gkyl_array_integrate_advance(gkyl_array_integrate *up, const struct gkyl_array *arr,
  double weight, const struct gkyl_range *range, double *out);

/**
 * Release memory associated with this updater.
 *
 * @param up array_integrate updater.
 */
void gkyl_array_integrate_release(gkyl_array_integrate *up);
]]

-- Integrated quantities calculator.
local CartFieldIntegrate = Proto(UpdaterBase)

function CartFieldIntegrate:init(tbl)
   CartFieldIntegrate.super.init(self, tbl)    -- Setup base object.

   -- Grid and basis.
   self.onGrid = assert(tbl.onGrid, 
      "Updater.CartFieldIntegrate: Must provide grid object using 'onGrid'.")
   self.basis = assert(tbl.basis, 
      "Updater.CartFieldIntegrate: Must provide phase-space basis object using 'basis'.")

   -- Number of components to set.
   self.numComponents = tbl.numComponents and tbl.numComponents or 1

   local useGPU = xsys.pickBool(tbl.useDevice, GKYL_USE_GPU or false)

   local opIn = tbl.operator or "none"
   local integ_op = nil
   -- These must match the definitions in the gkyl_array_integrate header file.
   if     opIn == "none" then integ_op = 0
   elseif opIn == "abs"  then integ_op = 1
   elseif opIn == "sq"   then integ_op = 2
   end
   assert(integ_op, "CartFieldIntegrate: operator not implemented. See available options.")

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   self._zero = ffi.gc(ffiC.gkyl_array_integrate_new(self.onGrid._zero, self.basis._zero, self.numComponents, integ_op, useGPU),
                       ffiC.gkyl_array_integrate_release)

   if useGPU then
      self.elmSz = ffi.sizeof("double")
      self.localVals, self.globalVals = cuda.Malloc(self.elmSz*self.numComponents), cudaMalloc(self.elmSz*self.numComponents)
      self.globalVals_ho = Lin.Vec(self.numComponents)
   else
      self.localVals, self.globalVals = Lin.Vec(self.numComponents), Lin.Vec(self.numComponents)
   end
end   

-- Advance method.
function CartFieldIntegrate:_advance(tCurr, inFld, outFld)
   local field, multfac = inFld[1], inFld[2]
   local dynVecOut      = outFld[1]

   local onRange = self.onGhosts and field:localExtRange() or field:localRange()

   ffiC.gkyl_array_integrate_advance(self._zero, field._zero, multfac or 1., onRange, self.localVals:data())

   local nvals = self.numComponents
   for i = 1, nvals do self.globalVals[i] = 0. end

   -- All-reduce across processors.
   local messenger = self.onGrid:getMessenger()
   if messenger then
      messenger:Allreduce(self.localVals:data(), self.globalVals:data(), nvals, "double", "sum", self.onGrid:commSet().comm)
   else
      for i = 1, nvals do self.globalVals[i] = self.localVals[i] end
   end

   -- Push result into DynVector.
   dynVecOut:appendData(tCurr, self.globalVals)
end

function CartFieldIntegrate:_advanceOnDevice(tCurr, inFld, outFld)
   local field, multfac = inFld[1], inFld[2]
   local dynVecOut      = outFld[1]

   local onRange = self.onGhosts and field:localExtRange() or field:localRange()

   ffiC.gkyl_array_integrate_advance(self._zero, field._zeroDevice, multfac or 1., onRange, self.localVals)

   local _ = cuda.Memset(self.globalVals, 0, self.numComponents*self.elmSz)

   -- All-reduce across processors.
   local messenger = self.onGrid:getMessenger()
   if messenger then
      messenger:Allreduce(self.localVals, self.globalVals, nvals, "double", "sum", self.onGrid:commSet().comm)
   else
      local _ = cuda.Memcpy(self.globalVals, self.localVals, nvals*self.elmSz, cuda.MemcpyDeviceToDevice)
   end

   -- Copy result to host.
   local err = cuda.Memcpy(self.globalVals_ho:data(), self.globalVals, nvals*self.elmSz, cuda.MemcpyDeviceToHost)
   assert_equal(cuda.Success, err, "CartFieldIntegrate: Memcpy failed.")

   -- Push result into DynVector.
   dynVecOut:appendData(tCurr, self.globalVals_ho)
end

return CartFieldIntegrate
