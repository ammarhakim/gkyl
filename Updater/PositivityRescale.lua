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
   self.i = 1
end   

-- advance method
function PositivityRescale:advance(tCurr, inFld, outFld, computeDiagnostics, zeroOut)
   local grid = self.onGrid
   local fIn, fOut = inFld[1], outFld[1]

   local ndim = self.basis:ndim()
   local numBasis = self.basis:numBasis()

   local fInIndexer = fIn:genIndexer()
   local fInPtr = fIn:get(1)
   local fOutIndexer = fOut:genIndexer()
   local fOutPtr = fOut:get(1)

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
      end
      self.delChange = 0.
      self.rescaledCells = 0.0
      self.tCurrOld = tCurr
      self.i=1
   end
 
   local localRange = fIn:localRange()   
   self.del2Change[self.i] = 0.
   for idx in localRange:rowMajorIter() do
      grid:setIndex(idx)

      fIn:fill(fInIndexer(idx), fInPtr)
      fOut:fill(fOutIndexer(idx), fOutPtr)

      local del2ChangeCell = ffiC.rescale(fInPtr:data(), fOutPtr:data(), ndim, numBasis, idx:data(), tCurr)
      if computeDiagnostics then 
         self.del2Change[self.i] = self.del2Change[self.i] + del2ChangeCell*grid:cellVolume()
         if del2ChangeCell ~= 0. then self.rescaledCells = self.rescaledCells + 1 end
      end
   end

   self.i = self.i+1
end

function PositivityRescale:write(tm, frame, nm)
   Mpi.Allreduce(self.rescaledCellsL:data():data(), self.rescaledCellsG:data():data(), self.rescaledCellsG:size()*2,
                 Mpi.DOUBLE, Mpi.SUM, self.onGrid:commSet().comm)
   for i=1, 4 do
      Mpi.Allreduce(self.del2ChangeL[i]:data():data(), self.del2ChangeG[i]:data():data(), self.del2ChangeG[i]:size()*2,
                    Mpi.DOUBLE, Mpi.SUM, self.onGrid:commSet().comm)
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
end

return PositivityRescale
