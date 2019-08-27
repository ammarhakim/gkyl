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
  bool check(const double *fIn, int ndim, int numBasis, int *idx, double tCurr);
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

   self.del2ChangeL = DataStruct.DynVector {
      numComponents = 1,
   }
   self.del2ChangeG = DataStruct.DynVector {
      numComponents = 1,
   }
   self.del2Change = 0.0
end   

-- advance method
-- check positivity of cell averages and report status = false if there are negative cell avgs
function PositivityCheck:_advance(tCurr, inFld, outFld)
   local fIn = inFld[1]

   local ndim = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   local fInIndexer = fIn:genIndexer()
   local fInPtr = fIn:get(1)

   local status = true

   local localRange = fIn:localRange()   
   for idx in localRange:rowMajorIter() do
      fIn:fill(fInIndexer(idx), fInPtr)
      status = status and fInPtr:data()[0] >= 0.
   end

   return status
end

-- check positivity of cell averages and control nodes. prints messages if violations.
function PositivityCheck:checkControlNodes(tCurr, inFld, outFld)
   local fIn = inFld[1]

   local ndim = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   local fInIndexer = fIn:genIndexer()
   local fInPtr = fIn:get(1)

   local status = true

   local localRange = fIn:localRange()   
   for idx in localRange:rowMajorIter() do
      fIn:fill(fInIndexer(idx), fInPtr)
      status = status and ffiC.check(fInPtr:data(), ndim, numBasis, idx:data(), tCurr)
   end
end

-- rescale to ensure positivity of control nodes. potentially diffusive.
function PositivityCheck:rescale(tCurr, inFld, outFld, pr)
   local fIn, fOut = inFld[1], outFld[1]

   local ndim = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   local fInIndexer = fIn:genIndexer()
   local fInPtr = fIn:get(1)
   local fOutIndexer = fOut:genIndexer()
   local fOutPtr = fOut:get(1)

   self.del2Change = 0.0
 
   -- this should be ext range since rescaling might be done after applyBc
   local localRange = fIn:localExtRange()   
   for idx in localRange:rowMajorIter() do
      fIn:fill(fInIndexer(idx), fInPtr)
      fOut:fill(fOutIndexer(idx), fOutPtr)

      local del2ChangeCell = ffiC.rescale(fInPtr:data(), fOutPtr:data(), ndim, numBasis, idx:data(), tCurr)
      self.del2Change = self.del2Change + del2ChangeCell
   end

   if pr then print(self.del2Change) end
   --self.del2ChangeL:appendData(tCurr, {self.del2Change})
   --self.del2ChangeG:appendData(tCurr, {0.0})
end


function PositivityCheck:write(tm, frame, nm)
   Mpi.Allreduce(self.del2ChangeL:data():data(), self.del2ChangeG:data():data(), self.del2ChangeG:size(),
                 Mpi.DOUBLE, Mpi.SUM, self.confGrid:commSet().comm)
   self.del2ChangeG:write(string.format("%s_%s_%d.bp", nm, "del2Change", frame), tm, frame, true)
   self.del2ChangeL:clear(0.0)
   self.del2ChangeG:clear(0.0)
end

return PositivityCheck
