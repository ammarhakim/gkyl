-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments of distribution function on a
-- rectangular (but potentially non-uniform) grid.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase = require "Updater.Base"
local Lin = require "Lib.Linalg"
local Proto = require "Lib.Proto"
local MomDecl = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"
local Mpi = require "Comm.Mpi"
local Range = require "Lib.Range"
local xsys = require "xsys"

-- Moments updater object
local DistFuncMomentCalc = Proto(UpdaterBase)

function DistFuncMomentCalc:isMomentNameGood(nm)
   if nm == "M0" or nm == "M1i" or nm == "M2ij" or nm == "M2" or nm == "M3i" or nm == "FiveMoments" or nm == "StarMoments" then
      return true
   end
   return false
end

function DistFuncMomentCalc:isGkMomentNameGood(nm)
   if nm == "GkM0" or nm == "GkM1" or nm == "GkM2par" or nm == "GkM2perp" or nm == "GkM2" or nm == "GkM3par" or nm == "GkM3perp" or nm == "GkThreeMoments" or nm == "GkStarMoments" then
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

   if mom == "FiveMoments" or mom == "StarMoments" or mom == "GkThreeMoments" or mom == "GkStarMoments" then self._fivemoments = true end

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
      print("DistFuncMomentCalc: Requested moment is", mom)
      assert(false, "DistFuncMomentCalc: Moments must be one of M0, M1i, M2ij, M2, M3i, FiveMoments, or StarMoments")
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
   self.phaseIdx = {}
   self.l, self.u = {}, {}
end

-- advance method
function DistFuncMomentCalc:_advance(tCurr, dt, inFld, outFld)
   local grid = self._onGrid
   local distf, mom1 = inFld[1], outFld[1]
   local mom2, mom3
   if self._fivemoments then 
      mom2 = outFld[2]
      mom3 = outFld[3] 
   end

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local confRange = mom1:localRange()
   local localItrFunc, localItrState = mom1:localRangeIter() -- get ahold of iterator for local region when using shared memory
   local phaseRange = distf:localRange()
   if self.onGhosts then -- extend range to config-space ghosts
      for dir = 1, cDim do
         confRange = confRange:extendDir(dir, mom1:lowerGhost(), mom1:upperGhost())
         localItrFunc, localItrState = mom1:localExtRangeIter() -- use extended shared memory region if onGhosts is true 
         phaseRange = phaseRange:extendDir(dir, distf:lowerGhost(), distf:upperGhost())
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

   for d = 1, vDim do
      self.l[d] = phaseRange:lower(cDim + d)
      self.u[d] = phaseRange:upper(cDim + d)
   end
   local velRange = Range.Range(self.l, self.u)

   -- barrier before shared loop
   local worldComm = self:getWorldComm()
   Mpi.Barrier(worldComm)

   -- loop, computing moments in each cell
   -- The configuration space loop
   for confIdx in localItrFunc, localItrState do
      mom1:fill(mom1Indexer(confIdx), mom1Itr)
      if self._fivemoments then
         mom2:fill(mom2Indexer(confIdx), mom2Itr)
         mom3:fill(mom3Indexer(confIdx), mom3Itr)
      end
      if self._isGK then
         self.bmag:fill(self.bmagIndexer(confIdx), self.bmagItr)
      end

      -- The velocity space loop
      for velIdx in velRange:colMajorIter() do
	 -- Construct the phase space index out of the configuration
	 -- space and velocity space indices
	 for d = 1, cDim do self.phaseIdx[d] = confIdx[d] end
	 for d = 1, vDim do self.phaseIdx[d + cDim] = velIdx[d] end

         -- get ahold of the grid, cell center coordinate, and grid spacing
         grid:setIndex(self.phaseIdx)
         grid:cellCenter(self.w)
         for d = 1, pDim do self.dxv[d] = grid:dx(d) end

         distf:fill(distfIndexer(self.phaseIdx), distfItr)

         if self._isGK then
            if self._fivemoments then
               self._momCalcFun(self.w:data(), self.dxv:data(), self.mass, self.bmagItr:data(), distfItr:data(), 
                             mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
            else
               self._momCalcFun(self.w:data(), self.dxv:data(), self.mass, self.bmagItr:data(), distfItr:data(), mom1Itr:data())
            end
         elseif self._fivemoments then 
            self._momCalcFun(self.w:data(), self.dxv:data(), distfItr:data(), mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
         else
            self._momCalcFun(self.w:data(), self.dxv:data(), distfItr:data(), mom1Itr:data())
         end
      end
   end

   if self.momfac ~= 1.0 then mom1:scale(self.momfac) end
   
   return true, GKYL_MAX_DOUBLE
end

return DistFuncMomentCalc
