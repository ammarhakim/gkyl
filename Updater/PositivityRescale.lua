-- Gkyl ------------------------------------------------------------------------
--
-- Updater to rescale field at nodes to ensure positivity.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local DataStruct  = require "DataStruct"
local Proto       = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local Mpi         = require "Comm.Mpi"
local ffi         = require "ffi"
local ffiC        = ffi.C

ffi.cdef[[
  double findMinNodalValue(double *fIn, int ndim, int polyOrder); 
  double rescale(const double *fIn, double *fOut, int ndim, int polyOrder, int numBasis, int *idx, double tCurr);
  double calcVolTermRescale(const double tCurr, const double dt, const double *fIn, const double weight, const double *fRhsSurf, double *fRhsVol, int ndim, int polyOrder, int numBasis, int *idx, int printWarnings);
  double rescaleVolTerm(const double tCurr, const double dt, const double *fIn, const double weight, const double *fRhsSurf, double *fRhsVol, int ndim, int polyOrder, int numBasis, int *idx, int printWarnings);
]]

local PositivityRescale = Proto(UpdaterBase)

function PositivityRescale:init(tbl)
   PositivityRescale.super.init(self, tbl) -- Setup base object.

   -- Grid and basis
   self.grid = assert(
      tbl.onGrid, "Updater.PositivityRescale: Must provide grid object using 'onGrid'")
   self.basis = assert(
      tbl.basis,
      "Updater.PositivityRescale: Must provide basis object using 'basis'")
   assert(self.basis:polyOrder()<=2, "Updater.PositivityRescale only implemented for p=1 and p=2")

   -- Number of components to set.
   self.numComponents = tbl.numComponents and tbl.numComponents or 1
   assert(self.numComponents == 1, "Updater.PositivityRescale only implemented for fields with numComponents = 1")

   self.del2ChangeByCell =  DataStruct.Field {
      onGrid        = self.grid,
      numComponents = 4,
      ghost         = {1, 1},
   }

   self.del2ChangeL = {
      DataStruct.DynVector {
         numComponents = 1,
      },
      DataStruct.DynVector {
         numComponents = 1,
      },
      DataStruct.DynVector {
         numComponents = 1,
      },
      DataStruct.DynVector {
         numComponents = 1,
      },
   }
   self.del2ChangeG = {
      DataStruct.DynVector {
         numComponents = 1,
      },
      DataStruct.DynVector {
         numComponents = 1,
      },
      DataStruct.DynVector {
         numComponents = 1,
      },
      DataStruct.DynVector {
         numComponents = 1,
      },
   }
   self.delChangeG = DataStruct.DynVector {
      numComponents = 1,
   }
   self.rescaledCellsL = DataStruct.DynVector {
      numComponents = 1,
   }
   self.rescaledCellsG = DataStruct.DynVector {
      numComponents = 1,
   }
   self.del2Change    = {0., 0., 0., 0.}
   self.delChange     = 0.0
   self.rescaledCells = 0.0
   self.tCurrOld      = 0.0
   self.rkIdx         = 1

   -- print positivity warnings when surface terms make control points negative?
   self.printWarnings = false

   self._first = true
end   

-- advance method: for diffusive (non-conservative) rescaling of f at end of timestep
function PositivityRescale:advance(tCurr, inFld, outFld)
   local grid      = self.grid
   local fIn, fOut = inFld[1], outFld[1]

   local ndim     = self.basis:ndim()
   local polyOrder = self.basis:polyOrder()
   local numBasis = self.basis:numBasis()

   if self._first then 
      self.fInIndexer        = fIn:genIndexer()
      self.fInPtr            = fIn:get(1)
      self.fOutIndexer       = fOut:genIndexer()
      self.fOutPtr           = fOut:get(1)
      self._first            = false
   end

   local localRange = fIn:localRange()   

   for idx in localRange:rowMajorIter() do
      grid:setIndex(idx)

      fIn:fill(self.fInIndexer(idx), self.fInPtr)
      fOut:fill(self.fOutIndexer(idx), self.fOutPtr)

      local del2ChangeCell = ffiC.rescale(self.fInPtr:data(), self.fOutPtr:data(), ndim, polyOrder, numBasis, idx:data(), tCurr)
   end
end

function PositivityRescale:calcVolTermRescale(tCurr, dtApprox, fIn, weights, weightDirs, fRhsSurf, fRhsVol, volRescale)
   -- Compute the rescale factor of the volume term needed to
   -- preserve positivity and store it in a CartField (volRescale).
   local ndim     = self.basis:ndim()
   local polyOrder = self.basis:polyOrder()
   local numBasis = self.basis:numBasis()

   local localRange     = fRhsSurf:localRange()   

   local fIn_ptr        = fIn:get(1)
   local fRhsSurf_ptr   = fRhsSurf:get(1)
   local fRhsVol_ptr    = fRhsVol:get(1)
   local volRescale_ptr = volRescale:get(1)

   local phaseIndexer   = fIn:genIndexer()

   local weights_ptr, weightsIndexer
   if weights then
      weights_ptr    = weights:get(1)
      weightsIndexer = weights:genIndexer()
   end

   for idx in localRange:rowMajorIter() do
      fIn:fill(phaseIndexer(idx), fIn_ptr)
      fRhsSurf:fill(phaseIndexer(idx), fRhsSurf_ptr)
      fRhsVol:fill(phaseIndexer(idx), fRhsVol_ptr)
      volRescale:fill(phaseIndexer(idx), volRescale_ptr)

      local weightfac = 1.0
      if weights then
         weightfac = 0.0
         weights:fill(weightsIndexer(idx), weights_ptr)
         for i, d in ipairs(weightDirs) do
            weightfac = weightfac + weights_ptr:data()[d]
         end
         weightfac = weightfac/weights_ptr:data()[0]
         if weightfac > 1 then weightfac = 1 end
      end
      
      volRescale_ptr[1] = ffiC.calcVolTermRescale(tCurr, dtApprox, fIn_ptr:data(), weightfac, fRhsSurf_ptr:data(),
                                                  fRhsVol_ptr:data(), ndim, polyOrder, numBasis, idx:data(), self.printWarnings)
   end
end

function PositivityRescale:rescaleVolTerm(tCurr, dtApprox, fIn, weights, weightDirs, fRhsSurf, fRhsVol)
   local ndim     = self.basis:ndim()
   local polyOrder = self.basis:polyOrder()
   local numBasis = self.basis:numBasis()

   local localRange      = fRhsSurf:localRange()   
   local fIn_ptr         = fIn:get(1)
   local fRhsSurf_ptr    = fRhsSurf:get(1)
   local fRhsVol_ptr     = fRhsVol:get(1)
   local fInIndexer      = fIn:genIndexer()
   local fRhsSurfIndexer = fRhsSurf:genIndexer()
   local fRhsVolIndexer  = fRhsVol:genIndexer()

   local weights_ptr, weightsIndexer
   if weights then
      weights_ptr    = weights:get(1)
      weightsIndexer = weights:genIndexer()
   end

   for idx in localRange:rowMajorIter() do
      fIn:fill(fInIndexer(idx), fIn_ptr)
      fRhsSurf:fill(fRhsSurfIndexer(idx), fRhsSurf_ptr)
      fRhsVol:fill(fRhsVolIndexer(idx), fRhsVol_ptr)

      local weightfac = 1.0
      if weights then
         weightfac = 0.0
         weights:fill(weightsIndexer(idx), weights_ptr)
         for i, d in ipairs(weightDirs) do
            weightfac = weightfac + weights_ptr:data()[d]
         end
         weightfac = weightfac/weights_ptr:data()[0]
         if weightfac > 1 then weightfac = 1 end
      end
      
      local scaler = ffiC.rescaleVolTerm(tCurr, dtApprox, fIn_ptr:data(), weightfac, fRhsSurf_ptr:data(), fRhsVol_ptr:data(), ndim, polyOrder, numBasis, idx:data(), self.printWarnings)
   end
end

return PositivityRescale
