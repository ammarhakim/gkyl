-- Gkyl ------------------------------------------------------------------------
--
-- Updater to rescale field at nodes to ensure positivity
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Proto = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"

ffi.cdef[[
  double findMinNodalValue(double *fIn, int ndim); 
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
end   

-- advance method
function PositivityRescale:_advance(tCurr, inFld, outFld)
   local grid = self.onGrid
   local fIn, fOut = inFld[1], outFld[1]

   local ndim = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   local fInIndexer = fIn:genIndexer()
   local fInPtr = fIn:get(1)
   local fOutIndexer = fOut:genIndexer()
   local fOutPtr = fOut:get(1)
 
   local localRange = fIn:localExtRange()   
   -- loop, computing integrated moments in each cell
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)

      fIn:fill(fInIndexer(idx), fInPtr)
      fOut:fill(fOutIndexer(idx), fOutPtr)
      
      local f0 = 1.0/math.sqrt(2.0)^ndim*fInPtr[1] -- cell average
      if f0 < 0 then
        if ndim == 1 then
          print(string.format("WARNING: negative cell avg %e in cell %d, tCurr = %e", f0, idx[1], tCurr))
        elseif ndim == 2 then
          print(string.format("WARNING: negative cell avg %e in cell %d %d, tCurr = %e", f0, idx[1], idx[2], tCurr))
        elseif ndim == 3 then
          print(string.format("WARNING: negative cell avg %e in cell %d %d %d, tCurr = %e", f0, idx[1], idx[2], idx[3], tCurr))
        elseif ndim == 4 then
          print(string.format("WARNING: negative cell avg %e in cell %d %d %d %d, tCurr = %e", f0, idx[1], idx[2], idx[3], idx[4], tCurr))
        elseif ndim == 5 then
          print(string.format("WARNING: negative cell avg %e in cell %d %d %d %d %d, tCurr = %e", f0, idx[1], idx[2], idx[3], idx[4], idx[5], tCurr))
        end
      end
      
      local fmin = ffi.C.findMinNodalValue(fInPtr:data(), ndim)
    
      -- compute scaling factor
      local theta = math.min(1, f0/(f0 - fmin + GKYL_EPSILON))

      -- modify moments. note no change to cell average
      fOutPtr[1] = fInPtr[1]
      for i = 2, numBasis do
         fOutPtr[i] = fInPtr[i]*theta
      end
   end
end

return PositivityRescale
