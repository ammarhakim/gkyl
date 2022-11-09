-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute integrated quantities on a Cartesian grid.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Alloc       = require "Lib.Alloc"
local Lin         = require "Lib.Linalg"
local Mpi         = require "Comm.Mpi"
local Proto       = require "Lib.Proto"
local Range       = require "Lib.Range"
local UpdaterBase = require "Updater.Base"
local ffi         = require "ffi"
local xsys        = require "xsys"

local ffiC = ffi.C
ffi.cdef [[
    void gkylCartFieldIntQuantV(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
    void gkylCartFieldIntQuantV2(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
    void gkylCartFieldIntQuantAbsV(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *fIn, double *out);
    void gkylCartFieldIntQuantGradPerpV2_2x_p1(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *inw, const double *fIn, double *out);
    void gkylCartFieldIntQuantGradPerpV2_2x_p2(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *inw, const double *fIn, double *out);
    void gkylCartFieldIntQuantGradPerpV2_3x_p1(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *inw, const double *fIn, double *out);
    void gkylCartFieldIntQuantGradPerpV2_3x_p2(
      int ndim, unsigned nc, unsigned nb, const double *dxv, const double *inw, const double *fIn, double *out);
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

   self.sqrt = false

   assert(tbl.quantity == "V" or tbl.quantity == "V2" or tbl.quantity == "AbsV" or tbl.quantity == "GradPerpV2" or tbl.quantity == "RmsV",
	  "CartFieldIntegratedQuantCalc: quantity must be one of V, V2, AbsV, GradPerpV2")

   if tbl.quantity == "RmsV" then self.sqrt = true; tbl.quantity = "V2" end

   if tbl.quantity == "GradPerpV2" then
      assert(self.numComponents==1 and self.basis:polyOrder()<3, 
             "CartFieldIntegratedQuantCalc: GradPerpV2 currently only implemented for p<3 and numComponents=1")
      tbl.quantity = tbl.quantity .. "_" .. tostring(self.onGrid:ndim()) .. "x_p" .. tostring(self.basis:polyOrder())

      self.advFunc = function(tCurr, inFld, outFld)
         self:advanceGradPerpV2(tCurr, inFld, outFld)
      end
   else
      self.advFunc = function(tCurr, inFld, outFld)
         self:advanceBasic(tCurr, inFld, outFld)
      end
   end

   self.updateFunc = ffiC["gkylCartFieldIntQuant"..tbl.quantity]

   self.timeIntegrate = xsys.pickBool(tbl.timeIntegrate, false)

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   -- For use in advance method.
   self.dxv        = Lin.Vec(self.basis:ndim())    -- Cell shape.
   self.localVals  = Lin.Vec(self.numComponents)
   self.globalVals = Lin.Vec(self.numComponents)
end   

function CartFieldIntegratedQuantCalc:advanceGradPerpV2(tCurr, inFld, outFld)
   local grid        = self.onGrid
   local field, vals = inFld[1], outFld[1]
   local multfac     = assert(inFld[2], "Updater.CartFieldIntegratedQuantCalc: GradPerpV2 requires that you pass the position-space jacobian.")

   local ndim  = self.basis:ndim()
   local nvals = self.numComponents

   local fieldIndexer = field:genIndexer()
   local fieldItr     = field:get(1)

   -- Clear local values.
   for i = 1, nvals do
      self.localVals[i]  = 0.0
      self.globalVals[i] = 0.0
   end

   -- Assume multfac is a field with the same number of cells and components as field.
   local multfacItr = multfac:get(1)
   -- Loop, computing integrated moments in each cell.
   local fieldRange = self.onGhosts and field:localExtRange() or field:localRange()
   for idx in fieldRange:rowMajorIter() do
      grid:setIndex(idx)
      grid:getDx(self.dxv)

      field:fill(fieldIndexer(idx), fieldItr)
      multfac:fill(fieldIndexer(idx), multfacItr)

      -- Compute integrated quantities.
      self.updateFunc(
         ndim, nvals, self.basis:numBasis(), self.dxv:data(), multfacItr:data(), fieldItr:data(), self.localVals:data())
   end

   -- All-reduce across processors and push result into dyn-vector.
   Mpi.Allreduce(
      self.localVals:data(), self.globalVals:data(), nvals, Mpi.DOUBLE, Mpi.SUM, self.onGrid:commSet().comm)

   vals:appendData(tCurr, self.globalVals)
end

function CartFieldIntegratedQuantCalc:advanceBasic(tCurr, inFld, outFld)
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

   -- Loop, computing integrated moments in each cell.
   local fieldRange = self.onGhosts and field:localExtRange() or field:localRange()
   for idx in fieldRange:rowMajorIter() do
      grid:setIndex(idx)
      grid:getDx(self.dxv)

      field:fill(fieldIndexer(idx), fieldItr)

      -- Compute integrated quantities.
      self.updateFunc(
         ndim, nvals, self.basis:numBasis(), self.dxv:data(), fieldItr:data(), self.localVals:data())
   end

   -- All-reduce across processors and push result into dyn-vector.
   Mpi.Allreduce(
      self.localVals:data(), self.globalVals:data(), nvals, Mpi.DOUBLE, Mpi.SUM, self.onGrid:commSet().comm)

   if multfac or self.sqrt then 
      for i = 1, nvals do
         if self.sqrt then self.globalVals[i] = math.sqrt(self.globalVals[i]) end
         if multfac then self.globalVals[i] = self.globalVals[i]*multfac end
      end
   end

   if self.timeIntegrate then
      local tLast, lastVals = vals:lastData()
      for i = 1, nvals do
         self.globalVals[i] = lastVals[i] + (tCurr - tLast)*self.globalVals[i]
      end
   end

   vals:appendData(tCurr, self.globalVals)
end

-- Advance method.
function CartFieldIntegratedQuantCalc:_advance(tCurr, inFld, outFld)
   self.advFunc(tCurr, inFld, outFld)
end

return CartFieldIntegratedQuantCalc
