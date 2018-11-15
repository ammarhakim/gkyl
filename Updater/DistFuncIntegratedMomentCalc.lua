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
local Lin = require "Lib.Linalg"
local MomDecl = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"
local Mpi = require "Comm.Mpi"
local Proto = require "Lib.Proto"
local Range = require "Lib.Range"
local UpdaterBase = require "Updater.Base"

-- Integrated moments updater object
local DistFuncIntegratedMomentCalc = Proto(UpdaterBase)

function DistFuncIntegratedMomentCalc:init(tbl)
   DistFuncIntegratedMomentCalc.super.init(self, tbl) -- setup base object

   self._onGrid = assert(
      tbl.onGrid, "Updater.DistFuncIntegratedMomentCalc: Must provide grid object using 'onGrid'")

   local phaseBasis = assert(
      tbl.phaseBasis,
      "Updater.DistFuncIntegratedMomentCalc: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis,
      "Updater.DistFuncIntegratedMomentCalc: Must provide configuration-space basis object using 'confBasis'")

   -- ensure sanity
   assert(phaseBasis:polyOrder() == confBasis:polyOrder(),
	  "Polynomial orders of phase-space and config-space basis must match")
   assert(phaseBasis:id() == confBasis:id(),
	  "Type of phase-space and config-space basis must match")

   -- dimension of spaces
   self._pDim = phaseBasis:ndim()
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   local id, polyOrder = phaseBasis:id(), phaseBasis:polyOrder()
   
   -- function to compute integrated moments
   self._intMomCalcFun = MomDecl.selectIntMomCalc(id, self._cDim, self._vDim, polyOrder)

   -- for use in _advance() method
   self.dxv = Lin.Vec(self._pDim) -- cell shape
   self.w = Lin.Vec(self._pDim) -- phase-space cell center
   self.localMom = Lin.Vec(5)
   self.globalMom = Lin.Vec(5)
end

-- advance method
function DistFuncIntegratedMomentCalc:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local distf, mom = inFld[1], outFld[1]

   local pDim = self._pDim
   local distfIndexer = distf:genIndexer()
   local distfItr = distf:get(1)

   -- clear local values
   self.localMom[1] = 0.0
   self.localMom[2] = 0.0
   self.localMom[3] = 0.0
   self.localMom[4] = 0.0
   self.localMom[5] = 0.0

   -- loop, computing integrated moments in each cell
   for idx in distf:localRangeIter() do
      grid:setIndex(idx)
      grid:cellCenter(self.w)
      for d = 1, pDim do self.dxv[d] = grid:dx(d) end
      distf:fill(distfIndexer(idx), distfItr)
      self._intMomCalcFun(self.w:data(), self.dxv:data(), distfItr:data(), self.localMom:data())
   end

   -- all-reduce across processors and push result into dyn-vector
   Mpi.Allreduce(
      self.localMom:data(), self.globalMom:data(), 5, Mpi.DOUBLE, Mpi.SUM, self:getComm())
   mom:appendData(tCurr, self.globalMom)
end

return DistFuncIntegratedMomentCalc
