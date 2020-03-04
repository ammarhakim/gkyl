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

return PositivityCheck
