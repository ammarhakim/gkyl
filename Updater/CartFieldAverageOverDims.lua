-- Gkyl ------------------------------------------------------------------------
--
-- Average a CartField over dimensions other than all the dimensions in the grid.
-- For example, average over y in a n-dimensional grid, producing an
-- (n-1)-dimensional quantity.
--
-- The user must supply the lower dimensional grid, basis and CartField.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto        = require "Lib.Proto"
local UpdaterBase  = require "Updater.Base"
local lume         = require "Lib.lume"
local Range        = require "Lib.Range"
local ModDecl      = require "Updater.averageDimsCalcData.CartFieldAverageOverDimsModDecl"

-- Moments updater object.
local CartFieldAverageOverDirs = Proto(UpdaterBase)

function CartFieldAverageOverDirs:init(tbl)
   CartFieldAverageOverDirs.super.init(self, tbl)    -- Setup base object.

   self.onGrid = assert(
      tbl.onGrid, "Updater.CartFieldAverageOverDirs: Must provide grid object using 'onGrid'.")
   local onBasis = assert(
      tbl.onBasis, "Updater.CartFieldAverageOverDirs: Must provide basis object using 'onBasis'.")
   self.avgDirs = assert(
      tbl.averageDirs, "Updater.CartFieldAverageOverDirs: Must provide directions to average in 'averageDirs'.")
   self.childGrid = assert(
      tbl.childGrid, "Updater.CartFieldAverageOverDirs: Must provide lower dim grid object using 'childGrid'.")
   local childBasis = assert(
      tbl.childBasis, "Updater.CartFieldAverageOverDirs: Must provide lower dim basis object using 'childBasis'.")

   local basisID, polyOrder = onBasis:id(), onBasis:polyOrder()
   assert(polyOrder == childBasis:polyOrder(),
          "Updater.CartFieldAverageOverDirs: Polynomial orders of phase-space and config-space basis must match")
   assert(basisID == childBasis:id(),
          "Updater.CartFieldAverageOverDirs: Type of phase-space and config-space basis must match")

   self.dim      = onBasis:ndim()
   self.childDim = childBasis:ndim()
   self.avgDim   = #self.avgDirs
   assert(self.dim > self.childDim and self.avgDim+self.childDim==self.dim, "Updater.CartFieldAverageOverDirs: can only integrate over some and not all of the dimensions. Use CartFieldIntegratedQuantCalc in order to integrate over all dimensions.")

   self.childDirs = {} -- Directions kept after integration.
   local c, k = 0, 0
   for d = 1, self.dim do
      if not lume.any(self.avgDirs, function(e) return e==d end) then
         c = c+1
         self.childDirs[c] = d
      end
   end
   assert(#self.childDirs == self.childDim, "Updater.CartFieldAverageOverDirs: self.avgDirs and childBasis are inconsistent.")

   -- Number of basis functions.
   self.numBasis = onBasis:numBasis()
   self.childNumBasis = childBasis:numBasis()

   self.avgOverDirsKernel = ModDecl.selectAverageOverDirsKernel(self.dim, basisID, polyOrder, self.avgDirs)

   self.idx = Lin.IntVec(self.dim)

   -- Factor multiplying contributions from each cell.
   self.rNumCells = 1
   for _, d in ipairs(self.avgDirs) do self.rNumCells = self.rNumCells/self.onGrid:numCells(d) end

end

-- Advance method.
function CartFieldAverageOverDirs:_advance(tCurr, inFld, outFld)

   local parentFld, childFld = inFld[1], outFld[1]

   local grid, childGrid = self.onGrid, self.childGrid

   local localRange, childRange = grid:localRange(), childGrid:localRange()

   local avgRangeLo, avgRangeUp = {}, {}
   local c = 0
   for _, d in ipairs(self.avgDirs) do
      c = c+1
      avgRangeLo[c], avgRangeUp[c] = localRange:lower(d), localRange:upper(d)
   end
   local avgRange = Range.Range(avgRangeLo, avgRangeUp)

   local childRangeDecomp = LinearDecomp.LinearDecompRange {
      range = childRange, numSplit = grid:numSharedProcs(), threadComm = grid:commSet().sharedComm }
   local tId = grid:subGridSharedId()   -- Local thread ID.

   local childIdxr, childPtr = childFld:genIndexer() , childFld:get(1)
   local parentIdxr, parentPtr = parentFld:genIndexer() , parentFld:get(1)

   childFld:clear(0.)

   -- Loop over the child dimensions (dimensions kept after averaging).
   for cidx in childRangeDecomp:rowMajorIter(tId) do

      for d = 1, self.childDim do self.idx[self.childDirs[d]] = cidx[d] end

      childFld:fill(childIdxr(cidx), childPtr)

      -- Loop over the dimensions we wish to average.
      for aidx in avgRange:rowMajorIter() do
         for d = 1, self.avgDim do self.idx[self.avgDirs[d]] = aidx[d] end

         parentFld:fill(parentIdxr(self.idx), parentPtr)

         self.avgOverDirsKernel(self.rNumCells, parentPtr:data(), childPtr:data())
      end
   end

end

return CartFieldAverageOverDirs
