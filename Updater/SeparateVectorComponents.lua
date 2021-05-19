-- Gkyl ------------------------------------------------------------------------
--
-- Given a CartField with several vector components (each a single scalar DG
-- field), separate them into individual CartFields.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local LinearDecomp = require "Lib.LinearDecomp"
local Proto        = require "Proto"
local UpdaterBase  = require "Updater.Base"

local SeparateVectorComponents = Proto(UpdaterBase)

function SeparateVectorComponents:init(tbl)
   SeparateVectorComponents.super.init(self, tbl) -- Setup base object.

   self._onGrid = assert(tbl.onGrid, "Updater.SeparateVectorComponents: Must provide grid object using 'onGrid'")
   self._basis  = assert(tbl.basis, "Updater.SeparateVectorComponents: Must specify basis functions to use using 'basis'")

   assert(self._onGrid:ndim() == self._basis:ndim(), "Dimensions of basis and grid must match")

   local ndim      = self._basis:ndim()
   local polyOrder = self._basis:polyOrder()
end

-- Advance method.
function SeparateVectorComponents:_advance(tCurr, inFld, outFld)
   local numBasis  = self._basis:numBasis()
   local fIn       = assert(inFld[1], "SeparateVectorComponents.advance: Must specify a vector input field")

   assert(inFld[1]:numComponents()/numBasis == #outFld, "SeparateVectorComponents.advance: Number of vector components of input field must match number of output fields")

   local numVals = #outFld

   for fI = 1, numVals do
      outFld[fI]:combineOffset(1., fIn, (fI-1)*numBasis)
   end

end

return SeparateVectorComponents
