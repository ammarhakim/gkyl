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
  double findMinNodalValue(double *fIn, int ndim); 
  double rescale(const double *fIn, double *fOut, int ndim, int numBasis, int *idx, double tCurr);
  double calcVolTermRescale(const double tCurr, const double dtApprox, const double *fIn, const double weight, const double *fRhsSurf, double *fRhsVol, int ndim, int numBasis, int *idx);
  double rescaleVolTerm(const double tCurr, const double dtApprox, const double *fIn, const double weight, const double *fRhsSurf, double *fRhsVol, int ndim, int numBasis, int *idx);
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
   assert(self.basis:polyOrder()==1, "Updater.PositivityRescale only implemented for p=1")

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

   self._first = true
end   

-- advance method
function PositivityRescale:advance(tCurr, inFld, outFld, computeDiagnostics, zeroOut)
   local grid      = self.grid
   local fIn, fOut = inFld[1], outFld[1]

   local ndim     = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   if self._first then 
      self.fInIndexer        = fIn:genIndexer()
      self.fInPtr            = fIn:get(1)
      self.fOutIndexer       = fOut:genIndexer()
      self.fOutPtr           = fOut:get(1)
      self.del2ChangeIndexer = self.del2ChangeByCell:genIndexer()
      self.del2ChangePtr     = self.del2ChangeByCell:get(1)
      self._first            = false
   end

   if computeDiagnostics == nil then computeDiagnostics = true end

   if zeroOut then 
      if computeDiagnostics and tCurr > 0 then 
         self.del2ChangeL[1]:appendData(self.tCurrOld, {self.del2Change[1]}) 
         self.del2ChangeG[1]:appendData(self.tCurrOld, {0})
         self.del2ChangeL[2]:appendData(self.tCurrOld, {self.del2Change[2]}) 
         self.del2ChangeG[2]:appendData(self.tCurrOld, {0})
         self.del2ChangeL[3]:appendData(self.tCurrOld, {self.del2Change[3]}) 
         self.del2ChangeG[3]:appendData(self.tCurrOld, {0})
         self.del2ChangeL[4]:appendData(self.tCurrOld, {self.del2Change[4]}) 
         self.del2ChangeG[4]:appendData(self.tCurrOld, {0})
         self.rescaledCellsL:appendData(self.tCurrOld, {self.rescaledCells}) 
         self.rescaledCellsG:appendData(self.tCurrOld, {0}) 
      end
      self.delChange     = 0.
      self.rescaledCells = 0.0
      self.tCurrOld      = tCurr
      self.rkIdx         = 1
   end
 
   local localRange = fIn:localRange()   
   self.del2Change[self.rkIdx] = 0.

   for idx in localRange:rowMajorIter() do
      grid:setIndex(idx)

      fIn:fill(self.fInIndexer(idx), self.fInPtr)
      fOut:fill(self.fOutIndexer(idx), self.fOutPtr)

      local del2ChangeCell = ffiC.rescale(self.fInPtr:data(), self.fOutPtr:data(), ndim, numBasis, idx:data(), tCurr)
      if computeDiagnostics then 
         self.del2ChangeByCell:fill(self.del2ChangeIndexer(idx), self.del2ChangePtr)
         self.del2ChangePtr:data()[self.rkIdx-1] = del2ChangeCell
         self.del2Change[self.rkIdx] = self.del2Change[self.rkIdx] + del2ChangeCell*grid:cellVolume()
         if del2ChangeCell ~= 0. then self.rescaledCells = self.rescaledCells + 1 end
      end
   end
   self.rkIdx = self.rkIdx+1
end

function PositivityRescale:calcVolTermRescale(tCurr, dtApprox, fIn, weights, weightDirs, fRhsSurf, fRhsVol, volRescale)
   -- Compute the rescale factor of the volume term needed to
   -- preserve positivity and store it in a CartField (volRescale).
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
         for i, d in pairs(weightDirs) do
            weightfac = weightfac + weights_ptr:data()[d]
         end
         weightfac = weightfac/weights_ptr:data()[0]
      end
      
      volRescale_ptr[1] = ffiC.calcVolTermRescale(tCurr, dtApprox, fIn_ptr:data(), weightfac, fRhsSurf_ptr:data(),
                                                  fRhsVol_ptr:data(), self.basis:ndim(), self.basis:numBasis(), idx:data())
   end
end

function PositivityRescale:scaleByCell(volRescale, fIn, fScaledOut)
   -- Multiply the (1-component) volume term rescale factor by a phase space field (e.g. distribution function).
   local localRange     = fIn:localRange()

   local volRescale_ptr = volRescale:get(1)
   local fIn_ptr        = fIn:get(1)
   local fScaledOut_ptr = fScaledOut:get(1)

   local phaseIndexer   = fIn:genIndexer()

   for idx in localRange:rowMajorIter() do
      volRescale:fill(phaseIndexer(idx), volRescale_ptr)
      fIn:fill(phaseIndexer(idx), fIn_ptr)
      fScaledOut:fill(phaseIndexer(idx), fScaledOut_ptr)

      for k = 1, self.basis:numBasis() do
         fScaledOut_ptr[k] = volRescale_ptr[1]*fIn_ptr[k]
      end
   end
end

function PositivityRescale:rescaleVolTerm(tCurr, dtApprox, fIn, weights, weightDirs, fRhsSurf, fRhsVol)
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
         for i, d in pairs(weightDirs) do
            weightfac = weightfac + weights_ptr:data()[d]
         end
         weightfac = weightfac/weights_ptr:data()[0]
      end
      
      local scaler = ffiC.rescaleVolTerm(tCurr, dtApprox, fIn_ptr:data(), weightfac, fRhsSurf_ptr:data(), fRhsVol_ptr:data(), self.basis:ndim(), self.basis:numBasis(), idx:data())
   end
end

function PositivityRescale:write(tm, frame, nm)
   self.del2ChangeByCell:write(string.format("%s_%s_%d.bp", nm, "del2ChangeByCell", frame), tm, frame, false)

   Mpi.Allreduce(self.rescaledCellsL:data():data(), self.rescaledCellsG:data():data(), self.rescaledCellsG:size()*2,
                 Mpi.DOUBLE, Mpi.SUM, self.grid:commSet().comm)
   for i=1, 4 do
      Mpi.Allreduce(self.del2ChangeL[i]:data():data(), self.del2ChangeG[i]:data():data(), self.del2ChangeG[i]:size()*2,
                    Mpi.DOUBLE, Mpi.SUM, self.grid:commSet().comm)
   end
   for j=1, self.del2ChangeG[1]:size() do
      local delChange = 0.
      for i=1, 4 do
         delChange = delChange + math.sqrt(self.del2ChangeG[i]:data()[j][1])
      end
      self.delChangeG:appendData(self.del2ChangeG[1]:timeMesh():data()[j-1], {delChange})
   end
   
   self.delChangeG:write(string.format("%s_%s_%d.bp", nm, "delChange", frame), tm, frame, true)
   self.rescaledCellsG:write(string.format("%s_%s_%d.bp", nm, "rescaledCells", frame), tm, frame, true)
   for i=1, 4 do
      self.del2ChangeL[i]:clear()
      self.del2ChangeG[i]:clear()
   end
   self.rescaledCellsL:clear()
   self.rescaledCellsG:clear()
   self.delChangeG:clear()
end

return PositivityRescale
