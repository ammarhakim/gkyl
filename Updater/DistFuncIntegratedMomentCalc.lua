-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute integrated moments (total particle, total
-- kinetic energy and total particle energy is computed)
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Alloc = require "Lib.Alloc"
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Range = require "Lib.Range"
local Proto = require "Lib.Proto"
local MomDecl = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"

-- Moments updater object
local DistFuncIntegratedMomentCalc = Proto(UpdaterBase)

function DistFuncIntegratedMomentCalc:init(tbl)
   DistFuncIntegratedMomentCalc.super.init(self, tbl) -- setup base object

   self._onGrid = assert(
      tbl.onGrid, "Updater.DistFuncIntegratedMomentCalc: Must provide grid object using 'onGrid'")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.DistFuncIntegratedMomentCalc: Must provide phase-space basis object using 'phaseBasis'")

   -- dimension of spaces
   self._pDim = phaseBasis:ndim() 

   local id, polyOrder = phaseBasis:id(), phaseBasis:polyOrder()
   
   -- function to compute integrated moments
   self._intMomCalcFun = MomDecl.selectIntMomCalc(id, self._vDim, polyOrder)

   -- for use in _advance() method
   self.dxv = Lin.Vec(self._pDim) -- cell shape
   self.w = Lin.Vec(self._pDim) -- phase-space cell center
   self.localMom = Lin.Vec(5)
   self.globalMom = Lin.Vec(5)
end

-- advance method
function DistFuncIntegratedMomentCalc:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local distf, mom = inFld[1], outFld[1]

   local pDim = self._pDim
   local distfIndexer = distf:genIndexer()
   local distfItr = distf:get(1)

   local localRange = distf:localRange()   
   -- loop, computing integrated moments in each cell
   for idx in localRange:colMajorIter() do
      grid:setIndex(idx)
      grid:cellCenter(self.w)
      for d = 1, pDim do self.dxv[d] = grid:dx(d) end
      distf:fill(distfIndexer(idx), distfItr)
      --self._intMomCalcFun(self.w:data(), self.dxv:data(), distfItr:data(), self.localMom)
   end

   -- all-reduce across processors
   local nodeComm = self:getNodeComm()
   Mpi.Allreduce(self.localMom:data(), self:globalMom(), 5, Mpi.DOUBLE, Mpi.SUM, nodeComm)
   
   return true, GKYL_MAX_DOUBLE
end

return DistFuncIntegratedMomentCalc
