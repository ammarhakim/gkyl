-- Gkyl ------------------------------------------------------------------------
--
-- Updater to rescale field at nodes to ensure positivity
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local DataStruct = require "DataStruct"
local Proto = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"
local ffiC = ffi.C

ffi.cdef[[
  double findMinNodalValue(double *fIn, int ndim); 
  double rescale(const double *fIn, double *fOut, int ndim, int numBasis, int *idx, double tCurr);
]]

local PositivityRescale = Proto(UpdaterBase)

function PositivityRescale:init(tbl)
   PositivityRescale.super.init(self, tbl) -- setup base object

   -- grid and basis
   self.onGrid = assert(
      tbl.onGrid, "Updater.PositivityRescale: Must provide grid object using 'onGrid'")
   self.basis = assert(
      tbl.basis,
      "Updater.PositivityRescale: Must provide basis object using 'basis'")
   assert(self.basis:polyOrder()==1, "Updater.PositivityRescale only implemented for p=1")

   -- number of components to set
   self.numComponents = tbl.numComponents and tbl.numComponents or 1
   assert(self.numComponents == 1, "Updater.PositivityRescale only implemented for fields with numComponents = 1")

   self.delChangeL = DataStruct.DynVector {
      numComponents = 1,
   }
   self.delChangeG = DataStruct.DynVector {
      numComponents = 1,
   }
   self.rescaledCellsL = DataStruct.DynVector {
      numComponents = 1,
   }
   self.delChange = 0.0
   self.rescaledCells = 0.0
   self.tCurrOld = 0.0
end   

-- advance method
function PositivityRescale:advance(tCurr, inFld, outFld, computeDiagnostics, zeroOut)
   local grid = self.onGrid
   local fIn, fOut = inFld[1], outFld[1]

   local ndim = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   local fInIndexer = fIn:genIndexer()
   local fInPtr = fIn:get(1)
   local fOutIndexer = fOut:genIndexer()
   local fOutPtr = fOut:get(1)

   if computeDiagnostics == nil then computeDiagnostics = true end

   if zeroOut then 
      if computeDiagnostics and tCurr > 0 then 
         self.delChangeL:appendData(self.tCurrOld, {self.delChange}) 
         self.rescaledCellsL:appendData(self.tCurrOld, {self.rescaledCells}) 
      end
      self.delChange = 0.
      self.rescaledCells = 0.0
      self.tCurrOld = tCurr
   end
 
   local localRange = fIn:localRange()   
   local del2Change = 0.
   for idx in localRange:rowMajorIter() do
      grid:setIndex(idx)

      fIn:fill(fInIndexer(idx), fInPtr)
      fOut:fill(fOutIndexer(idx), fOutPtr)

      local del2ChangeCell = ffiC.rescale(fInPtr:data(), fOutPtr:data(), ndim, numBasis, idx:data(), tCurr)
      if computeDiagnostics then 
         del2Change = del2Change + del2ChangeCell*grid:cellVolume()
         if del2ChangeCell ~= 0. then self.rescaledCells = self.rescaledCells + 1 end
      end
   end

   self.delChange = self.delChange + math.sqrt(del2Change)

   --self.del2ChangeG:appendData(tCurr, {0.0})
end

function PositivityRescale:write(tm, frame, nm)
   --Mpi.Allreduce(self.del2ChangeL:data():data(), self.del2ChangeG:data():data(), self.del2ChangeG:size(),
   --              Mpi.DOUBLE, Mpi.SUM, self.confGrid:commSet().comm)
   self.delChangeL:write(string.format("%s_%s_%d.bp", nm, "delChange", frame), tm, frame, true)
   self.rescaledCellsL:write(string.format("%s_%s_%d.bp", nm, "rescaledCells", frame), tm, frame, true)
   --self.del2ChangeL:clear(0.0)
   --self.del2ChangeG:clear(0.0)
end

return PositivityRescale
