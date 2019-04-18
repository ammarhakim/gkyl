-- Gkyl ------------------------------------------------------------------------
--
-- Updater to compute moments of distribution function on a
-- rectangular (but potentially non-uniform) grid.
--
-- If collisions (LBO for now) are included, this updater also computes the
-- boundary corrections and, if using a piecewise polynomial basis, the star
-- moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local Lin          = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local MomDecl      = require "Updater.momentCalcData.DistFuncMomentCalcModDecl"
local Proto        = require "Lib.Proto"
local UpdaterBase  = require "Updater.Base"
local lume         = require "Lib.lume"
local xsys         = require "xsys"

-- Moments updater object.
local DistFuncMomentCalc = Proto(UpdaterBase)

-- Valid moment names for Vlasov and GK equations.
local goodMomNames = {
   "M0", "M1i", "M2ij", "M2", "M3i", "FiveMoments", "FiveMomentsLBO"
}
local goodGkMomNames = {
   "GkM0", "GkM1", "GkM1proj", "GkM2par", "GkM2perp", "GkM2", "GkM3par", "GkM3perp",
   "GkThreeMoments", "GkThreeMomentsLBO"
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
   DistFuncMomentCalc.super.init(self, tbl)    -- Setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.DistFuncMomentCalc: Must provide grid object using 'onGrid'")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.DistFuncMomentCalc: Must provide phase-space basis object using 'phaseBasis'")
   local confBasis = assert(
      tbl.confBasis, "Updater.DistFuncMomentCalc: Must provide configuration-space basis object using 'confBasis'")

   self._basisID   = phaseBasis:id()
   self._polyOrder = phaseBasis:polyOrder()

   -- Dimension of spaces.
   self._pDim = phaseBasis:ndim() 
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Ensure sanity.
   assert(self._polyOrder == confBasis:polyOrder(),
	  "Polynomial orders of phase-space and config-space basis must match")
   assert(self._basisID == confBasis:id(),
	  "Type of phase-space and config-space basis must match")

   local mom = assert(
      tbl.moment, "Updater.DistFuncMomentCalc: Must provide moment to compute using 'moment'.")

   if mom == "FiveMoments" or mom == "GkThreeMoments" then
      self._fivemoments    = true
   elseif mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then
      self._fivemomentsLBO = true
   end

   -- Function to compute specified moment.
   self._isGK = false
   if self:isMomentNameGood(mom) then
      self._kinSpecies = "Vm"
      self._momCalcFun = MomDecl.selectMomCalc(mom, self._basisID, self._cDim, self._vDim, self._polyOrder)
   elseif self:isGkMomentNameGood(mom) then
      self._kinSpecies = "Gk"
      self._momCalcFun = MomDecl.selectGkMomCalc(mom, self._basisID, self._cDim, self._vDim, self._polyOrder)
      self._isGK = true
      assert(tbl.gkfacs, [[DistFuncMomentCalc: must provide a gkfacs table 
                        containing the species mass and the background magnetic field
                        to calculate a Gk moment]])
   else
      assert(false, "DistFuncMomentCalc: Moments must be one of M0, M1i, M2ij, M2, M3i, FiveMoments, or FiveMomentsLBO")
   end

   if self._polyOrder == 1 then
      self._StarM1iM2Calc = MomDecl.selectStarM1iM2Calc(self._kinSpecies, self._basisID, self._cDim, self._vDim)

   end
 
   self.momfac = 1.0
   if tbl.momfac then self.momfac = tbl.momfac end
   if tbl.gkfacs then
      self.mass        = tbl.gkfacs[1]
      self.bmag        = assert(tbl.gkfacs[2], "DistFuncMomentCalc: must provide bmag in gkfacs")
      self.bmagIndexer = self.bmag:genIndexer()
      self.bmagItr     = self.bmag:get(1)
   end

   self.onGhosts = xsys.pickBool(tbl.onGhosts, true)

   -- For use in _advance() method.
   self.dxv = Lin.Vec(self._pDim)    -- Cell shape.
   self.w   = Lin.Vec(self._pDim)    -- Phase-space cell center.
end

-- advance method
function DistFuncMomentCalc:_advance(tCurr, inFld, outFld)
   local grid        = self._onGrid
   local distf, mom1 = inFld[1], outFld[1]
   local mom2, mom3

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local localRange = distf:localRange()
   if self.onGhosts then -- extend range to config-space ghosts
      local cdirs = {}
      for dir = 1, cDim do 
	 localRange = localRange:extendDir(dir, distf:lowerGhost(), distf:upperGhost())
      end
   end

   local distfIndexer      = distf:genIndexer()
   local mom1Indexer       = mom1:genIndexer()
   local distfItr, mom1Itr = distf:get(1), mom1:get(1)

   local mom2Indexer, mom3Indexer
   local cMomBIndexer, cEnergyBIndexer
   local m0StarIndexer, m1StarIndexer, m2StarIndexer
   local mom2Itr, mom3Itr
   local cMomBItr, cEnergyBItr
   local m0StarItr, m1StarItr, m2StarItr
   mom1:scale(0.0) -- Zero out moments.
   if self._fivemoments then 
      mom2 = outFld[2]
      mom3 = outFld[3] 

      mom2Indexer = mom2:genIndexer()
      mom3Indexer = mom3:genIndexer() 

      mom2Itr = mom2:get(1)
      mom3Itr = mom3:get(1) 

      mom2:scale(0.0)
      mom3:scale(0.0)
   elseif self._fivemomentsLBO then 
      mom2      = outFld[2]
      mom3      = outFld[3] 
      cMomB     = outFld[4]
      cEenergyB = outFld[5] 

      mom2Indexer     = mom2:genIndexer()
      mom3Indexer     = mom3:genIndexer() 
      cMomBIndexer    = cMomB:genIndexer()
      cEnergyBIndexer = cEnergyB:genIndexer()

      mom2Itr     = mom2:get(1)
      mom3Itr     = mom3:get(1) 
      cMomBItr    = cMomB:get(1) 
      cEnergyBItr = cEnergyB:get(1) 

      mom2:scale(0.0)
      mom3:scale(0.0)
      cMomB:scale(0.0)
      cEnergyB:scale(0.0)
      if self._polyOrder == 1 then
         m0Star = outFld[6]
         m1Star = outFld[7]
         m2Star = outFld[8] 

         m0StarIndexer = m0Star:genIndexer()
         m1StarIndexer = m1Star:genIndexer()
         m2StarIndexer = m2Star:genIndexer()

         m0StarItr = m0Star:get(1) 
         m1StarItr = m1Star:get(1) 
         m2StarItr = m2Star:get(1) 

         m0Star:scale(0.0)
         m1Star:scale(0.0)
         m2Star:scale(0.0)
      end
   end

   -- construct ranges for nested loops
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = localRange:selectFirst(cDim), numSplit = grid:numSharedProcs() }
   local velRange = localRange:selectLast(vDim)
   local tId = grid:subGridSharedId() -- local thread ID

   local idx = Lin.IntVec(cDim+vDim)

   -- Separate the case with and without LBO collisions to reduce number of if statements.
   if self._fivemomentsLBO then 
      for cIdx in confRangeDecomp:rowMajorIter(tId) do
         grid:setIndex(cIdx)

         mom1:fill(mom1Indexer(cIdx), mom1Itr)
         mom2:fill(mom2Indexer(cIdx), mom2Itr)
         mom3:fill(mom3Indexer(cIdx), mom3Itr)
         cMomB:fill(cMomBIndexer(cIdx), cMomBItr)
         cEnergyB:fill(cEnergyBIndexer(cIdx), cEnergyBItr)

         if (self._isGK) then
            self.bmag:fill(self.bmagIndexer(cIdx), self.bmagItr)
         end

         -- Only when the contributions to m0Star from the first direction
         -- are collected, do we collect contributions to m1Star and m2Star.
         -- Also, since Gk velocities are organized as (vpar,mu) the velocity
         -- correction is only computed for the first velocity direction.
         local firstDir = true
   
         -- isLo=true current cell is the lower boundary cell.
         -- isLo=false current cell is the upper boundary cell.
         local isLo = true

         -- polyOrder=1 and >1 each use separate velocity grid loops to
         -- avoid evaluating (if polyOrder==1) at each velocity coordinate.
         if self._polyOrder > 1 then
            m0Star:fill(m0StarIndexer(cIdx), m0StarItr)
            m1Star:fill(m1StarIndexer(cIdx), m1StarItr)
            m2Star:fill(m2StarIndexer(cIdx), m2StarItr)

            for vDir = 1, vDim do
            end
         end
      end
   else
      -- outer loop is threaded and over configuration space
      for cIdx in confRangeDecomp:rowMajorIter(tId) do

         -- inner loop is over velocity space: no threading to avoid race
         -- conditions
         for vIdx in velRange:rowMajorIter() do
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

   end
   if self.momfac ~= 1.0 then mom1:scale(self.momfac) end
end

return DistFuncMomentCalc
