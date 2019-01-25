-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments of distribution function on a
-- rectangular (but potentially non-uniform) grid.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local MomDecl = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"
local Proto = require "Lib.Proto"
local UpdaterBase = require "Updater.Base"
local lume = require "Lib.lume"
local xsys = require "xsys"

-- Moments updater object
local DistFuncMomentCalc = Proto(UpdaterBase)

-- Valid moment names for Vlasov and GK equations
local goodMomNames = {
   "M0", "M1i", "M2ij", "M2", "M3i", "FiveMoments"
}
local goodGkMomNames = {
   "GkM0", "GkM1", "GkM1proj", "GkM2par", "GkM2perp", "GkM2", "GkM3par", "GkM3perp",
   "GkThreeMoments"
}

function DistFuncMomentCalc:isMomentNameGood(nm)
   if lume.find(goodMomNames, nm) then
      return true
   end
   return false
end

function DistFuncMomentCalc:isGkMomentNameGood(nm)
   if lume.find(goodGkMomNames, nm) then
      return true
   end
   return false
end

function DistFuncMomentCalc:init(tbl)
   DistFuncMomentCalc.super.init(self, tbl) -- setup base object

   self._onGrid = assert(
      tbl.onGrid, "Updater.DistFuncMomentCalc: Must provide grid object using 'onGrid'")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.DistFuncMomentCalc: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis, "Updater.DistFuncMomentCalc: Must provide configuration-space basis object using 'confBasis'")

   -- dimension of spaces
   self._pDim = phaseBasis:ndim() 
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- ensure sanity
   assert(phaseBasis:polyOrder() == confBasis:polyOrder(),
	  "Polynomial orders of phase-space and config-space basis must match")
   assert(phaseBasis:id() == confBasis:id(),
	  "Type of phase-space and config-space basis must match")

   local id, polyOrder = phaseBasis:id(), phaseBasis:polyOrder()
   
   local mom = assert(
      tbl.moment, "Updater.DistFuncMomentCalc: Must provide moment to compute using 'moment'.")

   if mom == "FiveMoments" or mom == "GkThreeMoments" then
      self._fivemoments = true
   end

   -- function to compute specified moment
   self._isGK = false
   if self:isMomentNameGood(mom) then
      self._momCalcFun = MomDecl.selectMomCalc(mom, id, self._cDim, self._vDim, polyOrder)
   elseif self:isGkMomentNameGood(mom) then
      self._momCalcFun = MomDecl.selectGkMomCalc(mom, id, self._cDim, self._vDim, polyOrder)
      self._isGK = true
      assert(tbl.gkfacs, [[DistFuncMomentCalc: must provide a gkfacs table 
                        containing the species mass and the background magnetic field
                        to calculate a Gk moment]])
   else
      assert(false, "DistFuncMomentCalc: Moments must be one of M0, M1i, M2ij, M2, M3i or FiveMoments")
   end
 
   self.momfac = 1.0
   if tbl.momfac then self.momfac = tbl.momfac end
   if tbl.gkfacs then
      self.mass = tbl.gkfacs[1]
      self.bmag = assert(tbl.gkfacs[2], "DistFuncMomentCalc: must provide bmag in gkfacs")
      self.bmagIndexer = self.bmag:genIndexer()
      self.bmagItr = self.bmag:get(1)
   end

   self.onGhosts = xsys.pickBool(tbl.onGhosts, true)

   -- for use in _advance() method
   self.dxv = Lin.Vec(self._pDim) -- cell shape
   self.w = Lin.Vec(self._pDim) -- phase-space cell center
end

-- advance method
function DistFuncMomentCalc:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid
   local distf, mom1 = inFld[1], outFld[1]
   local mom2, mom3
   if self._fivemoments then 
      mom2 = outFld[2]
      mom3 = outFld[3] 
   end

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local localRange = distf:localRange()
   if self.onGhosts then -- extend range to config-space ghosts
      local cdirs = {}
      for dir = 1, cDim do 
	 localRange = localRange:extendDir(dir, distf:lowerGhost(), distf:upperGhost())
      end
   end
   local distfIndexer = distf:genIndexer()
   local mom1Indexer = mom1:genIndexer()
   local mom2Indexer, mom3Indexer
   if self._fivemoments then 
      mom2Indexer = mom2:genIndexer()
      mom3Indexer = mom3:genIndexer() 
   end

   local distfItr, mom1Itr = distf:get(1), mom1:get(1)
   local mom2Itr, mom3Itr
   if self._fivemoments then 
      mom2Itr = mom2:get(1)
      mom3Itr = mom3:get(1) 
   end

   mom1:scale(0.0) -- zero out moments
   if self._fivemoments then 
      mom2:scale(0.0)
      mom3:scale(0.0)
   end

   -- construct ranges for nested loops
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = localRange:selectFirst(cDim), numSplit = grid:numSharedProcs() }
   local velRange = localRange:selectLast(vDim)
   local tId = grid:subGridSharedId() -- local thread ID

   local idx = Lin.IntVec(cDim+vDim)
   -- outer loop is threaded and over configuration space
   for cIdx in confRangeDecomp:colMajorIter(tId) do

      -- inner loop is over velocity space: no threading to avoid race
      -- conditions
      for vIdx in velRange:colMajorIter() do
	 --for d = 1, cDim do idx[d] = cIdx[d] end
         cIdx:copyInto(idx)
         for d = 1, vDim do idx[cDim+d] = vIdx[d] end
         --idx[cDim+1] = vIdx[1]
         --idx[cDim+2] = vIdx[2]
      
	 grid:setIndex(idx)
	 grid:cellCenter(self.w)
         grid:getDx(self.dxv)
      
	 distf:fill(distfIndexer(idx), distfItr)
	 mom1:fill(mom1Indexer(idx), mom1Itr)

	 if self._isGK then
	    self.bmag:fill(self.bmagIndexer(idx), self.bmagItr)
	    if self._fivemoments then
	       mom2:fill(mom2Indexer(idx), mom2Itr)
	       mom3:fill(mom3Indexer(idx), mom3Itr)
	       self._momCalcFun(self.w:data(), self.dxv:data(), self.mass, self.bmagItr:data(), distfItr:data(), 
	        		mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
	    else
	       self._momCalcFun(self.w:data(), self.dxv:data(), self.mass, self.bmagItr:data(), distfItr:data(), mom1Itr:data())
	    end
	 elseif self._fivemoments then 
	    mom2:fill(mom2Indexer(idx), mom2Itr)
	    mom3:fill(mom3Indexer(idx), mom3Itr)
	    self._momCalcFun(self.w:data(), self.dxv:data(), distfItr:data(), mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
	 else
	    self._momCalcFun(self.w:data(), self.dxv:data(), distfItr:data(), mom1Itr:data())
	 end
      end
   end
   if self.momfac ~= 1.0 then mom1:scale(self.momfac) end
end

return DistFuncMomentCalc
