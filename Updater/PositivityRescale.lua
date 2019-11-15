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
  double findMinNodalRatio(const double *fNum, const double *fDenom, double fac, int ndim);
  double rescale(const double *fIn, double *fOut, int ndim, int numBasis, int *idx, double tCurr);
  double rescaleVolTerm(const double *fOutSurf, const double dt, double *fVol, int ndim, int numBasis, int *idx);
]]

local PositivityRescale = Proto(UpdaterBase)

function PositivityRescale:init(tbl)
   PositivityRescale.super.init(self, tbl) -- setup base object

   -- grid and basis
   self.grid = assert(
      tbl.onGrid, "Updater.PositivityRescale: Must provide grid object using 'onGrid'")
   self.basis = assert(
      tbl.basis,
      "Updater.PositivityRescale: Must provide basis object using 'basis'")
   assert(self.basis:polyOrder()==1, "Updater.PositivityRescale only implemented for p=1")

   -- number of components to set
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
   self.del2Change = {0., 0., 0., 0.}
   self.delChange = 0.0
   self.rescaledCells = 0.0
   self.tCurrOld = 0.0
   self.rkIdx = 1

   self._first = true
end   

-- advance method
function PositivityRescale:advance(tCurr, inFld, outFld, computeDiagnostics, zeroOut)
   local grid = self.grid
   local fIn, fOut = inFld[1], outFld[1]

   local ndim = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   if self._first then 
      self.fInIndexer = fIn:genIndexer()
      self.fInPtr = fIn:get(1)
      self.fOutIndexer = fOut:genIndexer()
      self.fOutPtr = fOut:get(1)
      self.del2ChangeIndexer = self.del2ChangeByCell:genIndexer()
      self.del2ChangePtr = self.del2ChangeByCell:get(1)
      self._first = false
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
      self.delChange = 0.
      self.rescaledCells = 0.0
      self.tCurrOld = tCurr
      self.rkIdx=1
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

function PositivityRescale:rescaleVolTerm(tCurr, fOutSurf, dt, fVol)
   local localRange = fOutSurf:localRange()   
   local fOutSurfPtr = fOutSurf:get(1)
   local fVolPtr = fVol:get(1)
   local fOutSurfIndexer = fOutSurf:genIndexer()
   local fVolIndexer = fVol:genIndexer()

   for idx in localRange:rowMajorIter() do
      fOutSurf:fill(fOutSurfIndexer(idx), fOutSurfPtr)
      fVol:fill(fVolIndexer(idx), fVolPtr)

      local minOutSurf = ffiC.findMinNodalValue(fOutSurfPtr:data(), self.basis:ndim())
      --local f0 = fOutSurfPtr:data()[0]/math.sqrt(2)^self.basis:ndim()
      --if minOutSurf < -GKYL_EPSILON*math.abs(minOutSurf-2*f0)*4 then print("warning: at time ", tCurr, ", surface terms making control node negative, with value = ", minOutSurf, 
      --                                                                     " at idx = ", idx[1], idx[2],  idx[3], idx[4], idx[5], "eps = ", GKYL_EPSILON*math.abs(minOutSurf-2*f0)*4) end
      local scaler = ffiC.rescaleVolTerm(fOutSurfPtr:data(), dt, fVolPtr:data(), self.basis:ndim(), self.basis:numBasis(), idx:data())
   end
end

function PositivityRescale:findMinNodalRatio(fInPtr, fRhsSurfPtr, fac, ndim, idxPtr)
   return ffiC.findMinNodalRatio(fInPtr, fRhsSurfPtr, fac, ndim, idxPtr)
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
