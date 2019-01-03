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

ffi.cdef[[
  double limiter(const double *fIn, double *fHat, double *fOut, int ndim, int numBasis);
]]

local PositivityVolLimiter = Proto(UpdaterBase)

function PositivityVolLimiter:init(tbl)
   PositivityVolLimiter.super.init(self, tbl) -- setup base object

   -- grid and basis
   self.onGrid = assert(
      tbl.onGrid, "Updater.PositivityVolLimiter: Must provide grid object using 'onGrid'")
   self.basis = assert(
      tbl.basis,
      "Updater.PositivityVolLimiter: Must provide basis object using 'basis'")
   assert(self.basis:polyOrder()==1, "Updater.PositivityVolLimiter only implemented for p=1")

   -- number of components to set
   self.numComponents = tbl.numComponents and tbl.numComponents or 1
   assert(self.numComponents == 1, "Updater.PositivityVolLimiter only implemented for fields with numComponents = 1")
end   

-- advance method
function PositivityVolLimiter:_advance(tCurr, inFld, outFld)
   local grid = self.onGrid
   local fIn, fHat, fOut = inFld[1], inFld[2], outFld[1]

   local ndim = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   local fInIndexer = fIn:genIndexer()
   local fInPtr = fIn:get(1)
   local fHatIndexer = fHat:genIndexer()
   local fHatPtr = fHat:get(1)
   local fOutIndexer = fOut:genIndexer()
   local fOutPtr = fOut:get(1)

   local localRange = fIn:localRange()   
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)

      fIn:fill(fInIndexer(idx), fInPtr)
      fHat:fill(fHatIndexer(idx), fHatPtr)
      fOut:fill(fOutIndexer(idx), fOutPtr)

      ffi.C.limiter(fInPtr:data(), fHatPtr:data(), fOutPtr:data(), ndim, numBasis)
   end
end

return PositivityVolLimiter
