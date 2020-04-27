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
local Mpi = require "Comm.Mpi"
local ffi = require "ffi"
local ffiC = ffi.C

ffi.cdef[[
  double findMinNodalValue(double *fIn, int ndim); 
  double rescale(const double *fIn, double *fOut, int ndim, int numBasis, int *idx, double tCurr);
  bool check(const double *fIn, int ndim, int numBasis, int *idx, double tCurr, int rkIdx);
]]

local PositivityCheck = Proto(UpdaterBase)

function PositivityCheck:init(tbl)
   PositivityCheck.super.init(self, tbl) -- setup base object

   -- grid and basis
   self.onGrid = assert(
      tbl.onGrid, "Updater.PositivityCheck: Must provide grid object using 'onGrid'")
   self.basis = assert(
      tbl.basis,
      "Updater.PositivityCheck: Must provide basis object using 'basis'")
   assert(self.basis:polyOrder()==1, "Updater.PositivityCheck only implemented for p=1")

   -- number of components to set
   self.numComponents = tbl.numComponents and tbl.numComponents or 1
   assert(self.numComponents == 1, "Updater.PositivityCheck only implemented for fields with numComponents = 1")

   self.delChangeL = DataStruct.DynVector {
      numComponents = 1,
   }
   self.delChangeG = DataStruct.DynVector {
      numComponents = 1,
   }
   self.del2Change = 0.0
end   

-- check positivity of cell averages and control nodes. prints messages if violations. (only prints first cell avg violation)
-- status returns true if cell averages are positive
function PositivityCheck:_advance(tCurr, inFld, outFld)
   local fIn = inFld[1]
   local rkIdx = outFld[1] or 0

   local ndim = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   local fInIndexer = fIn:genIndexer()
   local fInPtr = fIn:get(1)

   local status = true

   local localRange = fIn:localRange()   
   for idx in localRange:rowMajorIter() do
      fIn:fill(fInIndexer(idx), fInPtr)
      myStatus = ffiC.check(fInPtr:data(), ndim, numBasis, idx:data(), tCurr, rkIdx)
      status = status and myStatus 
   end

   return status
end

-- rescale to ensure positivity of control nodes. potentially diffusive.
function PositivityCheck:rescale(tCurr, inFld, outFld, write)
   local fIn, fOut = inFld[1], outFld[1]

   local ndim = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   local fInIndexer = fIn:genIndexer()
   local fInPtr = fIn:get(1)
   local fOutIndexer = fOut:genIndexer()
   local fOutPtr = fOut:get(1)

   if write == nil then write = true end

   self.del2Change = 0.0
 
   local localRange = fIn:localRange()   
   for idx in localRange:rowMajorIter() do
      fIn:fill(fInIndexer(idx), fInPtr)
      fOut:fill(fOutIndexer(idx), fOutPtr)

      local del2ChangeCell = ffiC.rescale(fInPtr:data(), fOutPtr:data(), ndim, numBasis, idx:data(), tCurr)
      self.del2Change = self.del2Change + del2ChangeCell
   end

   if write then self.delChangeL:appendData(tCurr, {math.sqrt(self.del2Change)}) end
   --self.del2ChangeG:appendData(tCurr, {0.0})
end


function PositivityCheck:write(tm, frame, nm)
   --Mpi.Allreduce(self.del2ChangeL:data():data(), self.del2ChangeG:data():data(), self.del2ChangeL:size(),
                 --Mpi.DOUBLE, Mpi.SUM, self.onGrid:commSet().comm)
   self.delChangeL:write(string.format("%s_%s_%d.bp", nm, "delChange", frame), tm, frame, true)
   --self.del2ChangeL:clear(0.0)
   --self.del2ChangeG:clear(0.0)
end

return PositivityCheck
