local LinearDecomp     = require "Lib.LinearDecomp"
local Proto            = require "Proto"
local UpdaterBase      = require "Updater.Base"

local SeparateVectorComponents = Proto(UpdaterBase)

function SeparateVectorComponents:init(tbl)
   SeparateVectorComponents.super.init(self, tbl) -- Setup base object.

   self._onGrid  = assert(tbl.onGrid, "Updater.SeparateVectorComponents: Must provide grid object using 'onGrid'")
   self._basis   = assert(tbl.basis, "Updater.SeparateVectorComponents: Must specify basis functions to use using 'basis'")

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")

   local ndim      = self._basis:ndim()
   local polyOrder = self._basis:polyOrder()
end

-- Advance method.
function SeparateVectorComponents:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local ndim      = grid:ndim()
   local numBasis  = self._basis:numBasis()
   local polyOrder = self._basis:polyOrder()
   local fIn = assert(inFld[1], "SeparateVectorComponents.advance: Must specify a vector input field")
   assert(inFld[1]:numComponents()/numBasis == #outFld, "SeparateVectorComponents.advance: Number of vector components of input field must match number of output fields")
   local numVals = #outFld

   local tId = grid:subGridSharedId() -- Local thread ID.
   -- Object to iterate over only region owned by local SHM thread.
   local localRangeDecomp = LinearDecomp.LinearDecompRange {
	 range = fIn:localRange(), numSplit = grid:numSharedProcs() }

   local indexer = fIn:genIndexer()
   local outIndexer = outFld[1]:genIndexer()
   local fItr = fIn:get(1)
   for idx in localRangeDecomp:colMajorIter(tId) do
      fIn:fill(indexer(idx), fItr)
      for j = 1, numVals do
         local outItr = outFld[j]:get(1)
         outFld[j]:fill(outIndexer(idx), outItr)
         for i = 1, numBasis do
            outItr[i] = fItr[(j-1)*numBasis+i]
         end
      end
   end
end

return SeparateVectorComponents
