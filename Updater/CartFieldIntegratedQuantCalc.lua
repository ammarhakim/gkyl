-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute integrated quantities on a Cartesian grid.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Alloc        = require "Lib.Alloc"
local LinearDecomp = require "Lib.LinearDecomp"
local Lin          = require "Lib.Linalg"
local Mpi          = require "Comm.Mpi"
local Proto        = require "Lib.Proto"
local Range        = require "Lib.Range"
local UpdaterBase  = require "Updater.Base"
local ffi          = require "ffi"

local ffiC = ffi.C
ffi.cdef [[
    void gkylCartFieldIntQuantV(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
    void gkylCartFieldIntQuantV2(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
    void gkylCartFieldIntQuantAbsV(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
    void gkylCartFieldIntQuantGradPerpV2(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
]]

-- Integrated quantities calculator.
local CartFieldIntegratedQuantCalc = Proto(UpdaterBase)

function CartFieldIntegratedQuantCalc:init(tbl)
   CartFieldIntegratedQuantCalc.super.init(self, tbl)    -- Setup base object.

   -- Grid and basis.
   self.onGrid = assert(
      tbl.onGrid, "Updater.CartFieldIntegratedQuantCalc: Must provide grid object using 'onGrid'")
   self.basis = assert(
      tbl.basis,
      "Updater.CartFieldIntegratedQuantCalc: Must provide phase-space basis object using 'basis'")

   -- Number of components to set.
   self.numComponents = tbl.numComponents and tbl.numComponents or 1

   assert(tbl.quantity == "V" or tbl.quantity == "V2" or tbl.quantity == "AbsV" or tbl.quantity == "GradPerpV2",
	  "CartFieldIntegratedQuantCalc: quantity must be one of V, V2, AbsV, GradPerpV2")
   self.updateFunc = ffiC["gkylCartFieldIntQuant"..tbl.quantity]

   if tbl.quantity == "GradPerpV2" then assert(self.numComponents==1 and self.basis:polyOrder()==1, 
          "CartFieldIntegratedQuantCalc: GradPerpV2 currently only implemented for p=1 and numComponents=1")
   end

   -- For use in advance method.
   self.dxv        = Lin.Vec(self.basis:ndim())    -- Cell shape.
   self.localVals  = Lin.Vec(self.numComponents)
   self.globalVals = Lin.Vec(self.numComponents)
end   

-- Advance method.
function CartFieldIntegratedQuantCalc:_advance(tCurr, inFld, outFld)
   local grid        = self.onGrid
   local field, vals = inFld[1], outFld[1]
   local multfac     = inFld[2]

   local ndim  = self.basis:ndim()
   local nvals = self.numComponents

   local fieldIndexer = field:genIndexer()
   local fieldItr     = field:get(1)

   -- Clear local values.
   for i = 1, nvals do
      self.localVals[i]  = 0.0
      self.globalVals[i] = 0.0
   end

   -- Construct range for shared memory.
   local fieldRange = field:localRange()
   local fieldRangeDecomp = LinearDecomp.LinearDecompRange {
      range = fieldRange:selectFirst(ndim), numSplit = grid:numSharedProcs() }
   local tId = grid:subGridSharedId()    -- Local thread ID.

   -- Loop, computing integrated moments in each cell.
   for idx in fieldRangeDecomp:rowMajorIter(tId) do
      grid:setIndex(idx)
      grid:getDx(self.dxv)

      field:fill(fieldIndexer(idx), fieldItr)
      -- Compute integrated quantities.
      self.updateFunc(
	 ndim, nvals, self.basis:numBasis(), self.dxv:data(), fieldItr:data(), self.localVals:data())
   end

   -- All-reduce across processors and push result into dyn-vector.
   Mpi.Allreduce(
      self.localVals:data(), self.globalVals:data(), nvals, Mpi.DOUBLE, Mpi.SUM, self:getComm())

   if multfac then 
      for i = 1, nvals do
         self.globalVals[i] = self.globalVals[i]*multfac
      end
   end

   vals:appendData(tCurr, self.globalVals)
end

return CartFieldIntegratedQuantCalc
