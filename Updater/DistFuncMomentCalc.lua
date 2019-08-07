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

   -- Number of basis functions.
   self._numBasisP = phaseBasis:numBasis()

   -- Ensure sanity.
   assert(self._polyOrder == confBasis:polyOrder(),
	  "Polynomial orders of phase-space and config-space basis must match")
   assert(self._basisID == confBasis:id(),
	  "Type of phase-space and config-space basis must match")

   local mom = assert(
      tbl.moment, "Updater.DistFuncMomentCalc: Must provide moment to compute using 'moment'.")

   if mom == "FiveMoments" or mom == "GkThreeMoments" or mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then
      self._fiveMoments = true
      if mom == "FiveMomentsLBO" or mom == "GkThreeMomentsLBO" then
         self._fiveMomentsLBO = true
         if mom == "FiveMomentsLBO" then    -- Rename this variable to call the right kernel below.
            mom = "FiveMoments"
         elseif mom == "GkThreeMomentsLBO" then
            mom = "GkThreeMoments"
         end
      end
   end

   -- Function to compute specified moment.
   self._isGk = false
   if self:isMomentNameGood(mom) then
      self._kinSpecies = "Vm"
      self._momCalcFun = MomDecl.selectMomCalc(mom, self._basisID, self._cDim, self._vDim, self._polyOrder)
   elseif self:isGkMomentNameGood(mom) then
      self._kinSpecies = "Gk"
      self._momCalcFun = MomDecl.selectGkMomCalc(mom, self._basisID, self._cDim, self._vDim, self._polyOrder)
      self._isGk       = true
      assert(tbl.gkfacs, [[DistFuncMomentCalc: must provide a gkfacs table 
                        containing the species mass and the background magnetic field
                        to calculate a Gk moment]])
   else
      assert(false, "DistFuncMomentCalc: Moments must be one of M0, M1i, M2ij, M2, M3i, FiveMoments, or FiveMomentsLBO")
   end

   self.momfac = 1.0
   if tbl.momfac then self.momfac = tbl.momfac end
   if tbl.gkfacs then
      self.mass        = tbl.gkfacs[1]
      self.bmag        = assert(tbl.gkfacs[2], "DistFuncMomentCalc: must provide bmag in gkfacs")
      self.bmagItr     = self.bmag:get(1)
   end

   -- Cell index, center and length right of a cell-boundary (also used for current cell for p>1).
   self.idxP = Lin.IntVec(self._pDim)
   self.xcP  = Lin.Vec(self._pDim)
   self.dxP  = Lin.Vec(self._pDim)

   if self._fiveMomentsLBO then
      -- If vDim>1, intFac=2*pi/m or 4*pi/m.
      self._intFac = Lin.Vec(self._vDim)
      for d = 1,self._vDim do
        self._intFac[d] = 1.0
      end
      if self._isGk and (self._vDim > 1) then -- A (vpar,mu) simulation has 3 physical velocity dimensions.
         self._intFac[1] = 2.0*math.pi/self.mass
         self._intFac[2] = 4.0*math.pi/self.mass
      end
      self._isFirst   = true
      self._perpRange = {}    -- Perp ranges in velocity directions.
      if self._polyOrder == 1 then
         self._StarM1iM2Calc = MomDecl.selectStarM1iM2Calc(self._kinSpecies, self._basisID, self._cDim, self._vDim)
         -- Cell index, center and length left of a cell-boundary.
         self.idxM = Lin.IntVec(self._pDim)
         self.xcM  = Lin.Vec(self._pDim)
         self.dxM  = Lin.Vec(self._pDim)
      end
   end

   self.onGhosts = xsys.pickBool(tbl.onGhosts, true)
end

-- advance method
function DistFuncMomentCalc:_advance(tCurr, inFld, outFld)
   local grid        = self._onGrid
   local distf, mom1 = inFld[1], outFld[1]

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim

   local phaseRange = distf:localRange()
   if self.onGhosts then -- extend range to config-space ghosts
      local cdirs = {}
      for dir = 1, cDim do 
	 phaseRange = phaseRange:extendDir(dir, distf:lowerGhost(), distf:upperGhost())
      end
   end

   local phaseIndexer      = distf:genIndexer()
   local confIndexer       = mom1:genIndexer()
   local distfItr, mom1Itr = distf:get(1), mom1:get(1)

   -- Construct ranges for nested loops.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(cDim), numSplit = grid:numSharedProcs() }
   local velRange = phaseRange:selectLast(vDim)
   local tId      = grid:subGridSharedId()    -- Local thread ID.

   local mom2, mom3
   local mom2Itr, mom3Itr
   mom1:scale(0.0) -- Zero out moments.

   local cMomB, cEnergyB
   local m0Star, m1Star, m2Star
   local cMomBItr, cEnergyBItr
   local m0StarItr, m1StarItr, m2StarItr
   local distfItrP, distfItrM
   local uCorrection, vtSqCorrection, StarM0Calc    -- Kernel pointers.
   if self._fiveMoments then 
      mom2 = outFld[2]
      mom3 = outFld[3] 

      mom2Itr = mom2:get(1)
      mom3Itr = mom3:get(1) 

      mom2:scale(0.0)
      mom3:scale(0.0)
      if self._fiveMomentsLBO then 
         cMomB    = outFld[4]
         cEnergyB = outFld[5] 

         cMomBItr    = cMomB:get(1) 
         cEnergyBItr = cEnergyB:get(1) 

         cMomB:scale(0.0)
         cEnergyB:scale(0.0)

         -- Added for corrections and star moments.
         -- Distribution functions left and right of a cell-boundary.
         distfItrP = distf:get(1)
         distfItrM = distf:get(1)

         if self._polyOrder == 1 then
            m0Star = outFld[6]
            m1Star = outFld[7]
            m2Star = outFld[8]

            m0StarItr = m0Star:get(1) 
            m1StarItr = m1Star:get(1) 
            m2StarItr = m2Star:get(1) 

            m0Star:scale(0.0)
            m1Star:scale(0.0)
            m2Star:scale(0.0)
         end
      end
   end

   -- Separate the case with and without LBO collisions to reduce number of if statements.
   if (self._fiveMomentsLBO and (self._polyOrder==1)) then 

      -- Outer loop is threaded and over configuration space.
      for cIdx in confRangeDecomp:rowMajorIter(tId) do

         cIdx:copyInto(self.idxP)

         mom1:fill(confIndexer(cIdx), mom1Itr)
         mom2:fill(confIndexer(cIdx), mom2Itr)
         mom3:fill(confIndexer(cIdx), mom3Itr)

         if self._isGk then
            self.bmag:fill(confIndexer(cIdx), self.bmagItr)
         end

         -- Now loop over velocity space boundary surfaces to compute boundary corrections.
         cMomB:fill(confIndexer(cIdx), cMomBItr)
         cEnergyB:fill(confIndexer(cIdx), cEnergyBItr)

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
      
         -- To have energy conservation with piece-wise linear, we must use
         -- star moments in the second equation of the weak system solved
         -- in SelfPrimMoments.
         m0Star:fill(confIndexer(cIdx), m0StarItr)
         m1Star:fill(confIndexer(cIdx), m1StarItr)
         m2Star:fill(confIndexer(cIdx), m2StarItr)

         for vDir = 1, vDim do
            if (not self._isGk) or (self._isGk and firstDir) then
               StarM0Calc  = MomDecl.selectStarM0Calc(vDir, self._kinSpecies, self._basisID, cDim, vDim)
               uCorrection = MomDecl.selectBoundaryFintegral(vDir, self._kinSpecies, self._basisID, cDim, vDim, self._polyOrder)
            end
            vtSqCorrection = MomDecl.selectBoundaryVFintegral(vDir, self._kinSpecies, self._basisID, cDim, vDim, self._polyOrder)
   
            -- Lower/upper bounds in direction 'vDir': edge indices (including outer edges).
            local dirLoIdx, dirUpIdx = phaseRange:lower(cDim+vDir), phaseRange:upper(cDim+vDir)+1
   
            if self._isFirst then
               -- Restricted velocity range.
               -- Velocity integral in m0Star does not include last cell.
               self._perpRange[vDir] = phaseRange
               for cd = 1, cDim do
                  self._perpRange[vDir] = self._perpRange[vDir]:shorten(cd) -- shorten configuration range.
               end
               self._perpRange[vDir] = self._perpRange[vDir]:shorten(cDim+vDir) -- velocity range orthogonal to 'vDir'.
            end
            local perpRange = self._perpRange[vDir]

            -- Outer loop is over directions orthogonal to 'vDir' and
            -- inner loop is over 1D slice in 'vDir'.
            for vPerpIdx in perpRange:rowMajorIter() do
               vPerpIdx:copyInto(self.idxM); vPerpIdx:copyInto(self.idxP)
               for d = 1, cDim do self.idxM[d] = cIdx[d] end
               for d = 1, cDim do self.idxP[d] = cIdx[d] end
      
               for i = dirLoIdx, dirUpIdx do     -- This loop is over edges.
                  self.idxM[cDim+vDir], self.idxP[cDim+vDir] = i-1, i -- Cell left/right of edge 'i'.
      
                  grid:setIndex(self.idxM)
                  grid:getDx(self.dxM)
                  grid:cellCenter(self.xcM)
      
                  grid:setIndex(self.idxP)
                  grid:getDx(self.dxP)
                  grid:cellCenter(self.xcP)
      
                  distf:fill(phaseIndexer(self.idxM), distfItrM)
                  distf:fill(phaseIndexer(self.idxP), distfItrP)
      
                  if i>dirLoIdx and i<dirUpIdx then
                     if (self._isGk) then
                        if (firstDir) then
                           StarM0Calc(self._intFac[1], self.xcM:data(), self.xcP:data(), self.dxM:data(), self.dxP:data(), distfItrM:data(), distfItrP:data(), m0StarItr:data())
                        end
                     else
                        StarM0Calc(self.xcM:data(), self.xcP:data(), self.dxM:data(), self.dxP:data(), distfItrM:data(), distfItrP:data(), m0StarItr:data())
                     end
                  end
                  if firstDir and i<dirUpIdx then
                     if self._isGk then
                        self._momCalcFun(self.xcP:data(), self.dxP:data(), self.mass, self.bmagItr:data(), distfItrP:data(), 
                         		 mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
                        self._StarM1iM2Calc(self.xcP:data(), self.dxP:data(), self._intFac[1], self.mass, self.bmagItr:data(), distfItrP:data(), m1StarItr:data(), m2StarItr:data())
                     else
                        self._momCalcFun(self.xcP:data(), self.dxP:data(), distfItrP:data(), mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
                        self._StarM1iM2Calc(self.xcP:data(), self.dxP:data(), distfItrP:data(), m1StarItr:data(), m2StarItr:data())
                     end
                  end

                  if i==dirLoIdx or i==dirUpIdx-1 then
                     local vBound = 0.0
                     -- Careful: for vBound below we assume idxP was set after idxM above.
                     if isLo then
                        vBound = grid:cellLowerInDir(cDim + vDir)
                     else
                        vBound = grid:cellUpperInDir(cDim + vDir)
                     end
                     if (self._isGk) then
                        if (firstDir) then
                           uCorrection(isLo, self._intFac[1], vBound, self.dxP:data(), distfItrP:data(), cMomBItr:data())
                        end
                        vtSqCorrection(isLo, self._intFac[vDir], vBound, self.dxP:data(), distfItrP:data(), cEnergyBItr:data())
                     else
                        uCorrection(isLo, vBound, self.dxP:data(), distfItrP:data(), cMomBItr:data())
                        vtSqCorrection(isLo, vBound, self.dxP:data(), distfItrP:data(), cEnergyBItr:data())
                     end
      
                     isLo = not isLo
                  end    -- i==dirLoIdx or i==dirUpIdx-1.
      
               end    -- Loop over edges.
            end    -- Loop over directions perpendicular to vDir.
            firstDir = false

         end    -- vDir loop.
      end    -- Loop over configuration space.


   else    -- if self._fiveMomentsLBO and self._polyOrder=1. 

      -- Outer loop is threaded and over configuration space.
      for cIdx in confRangeDecomp:rowMajorIter(tId) do

         cIdx:copyInto(self.idxP)

         -- Inner loop is over velocity space: no threading to avoid race conditions.
         for vIdx in velRange:rowMajorIter() do
            for d = 1, vDim do self.idxP[cDim+d] = vIdx[d] end
         
            grid:setIndex(self.idxP)
            grid:cellCenter(self.xcP)
            grid:getDx(self.dxP)
         
            distf:fill(phaseIndexer(self.idxP), distfItr)
            mom1:fill(confIndexer(cIdx), mom1Itr)

            if self._isGk then
               self.bmag:fill(confIndexer(cIdx), self.bmagItr)
               if self._fiveMoments then
                  mom2:fill(confIndexer(cIdx), mom2Itr)
                  mom3:fill(confIndexer(cIdx), mom3Itr)
                  self._momCalcFun(self.xcP:data(), self.dxP:data(), self.mass, self.bmagItr:data(), distfItr:data(), 
                   		mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
               else
                  self._momCalcFun(self.xcP:data(), self.dxP:data(), self.mass, self.bmagItr:data(), distfItr:data(), mom1Itr:data())
               end
            elseif self._fiveMoments then 
               mom2:fill(confIndexer(cIdx), mom2Itr)
               mom3:fill(confIndexer(cIdx), mom3Itr)
               self._momCalcFun(self.xcP:data(), self.dxP:data(), distfItr:data(), mom1Itr:data(), mom2Itr:data(), mom3Itr:data())
            else
               self._momCalcFun(self.xcP:data(), self.dxP:data(), distfItr:data(), mom1Itr:data())
            end
         end

         if self._fiveMomentsLBO then  -- and polyOrder>1.
            -- Now loop over velocity space boundary surfaces to compute boundary corrections.
            cMomB:fill(confIndexer(cIdx), cMomBItr)
            cEnergyB:fill(confIndexer(cIdx), cEnergyBItr)

            -- Only when the contributions to m0Star from the first direction
            -- are collected, do we collect contributions to m1Star and m2Star.
            -- Also, since Gk velocities are organized as (vpar,mu) the velocity
            -- correction is only computed for the first velocity direction.
            local firstDir = true
   
            -- isLo=true current cell is the lower boundary cell.
            -- isLo=false current cell is the upper boundary cell.
            local isLo = true

            for vDir = 1, vDim do
               if (not self._isGk) or (self._isGk and firstDir) then
                  uCorrection = MomDecl.selectBoundaryFintegral(vDir, self._kinSpecies, self._basisID, cDim, vDim, self._polyOrder)
               end
               vtSqCorrection = MomDecl.selectBoundaryVFintegral(vDir, self._kinSpecies, self._basisID, cDim, vDim, self._polyOrder)

               -- Lower/upper bounds in direction 'vDir': cell indices.
               local dirLoIdx, dirUpIdx = phaseRange:lower(cDim+vDir), phaseRange:upper(cDim+vDir)

               if self._isFirst then
                  self._perpRange[vDir] = phaseRange
                  for cd = 1, cDim do
                     self._perpRange[vDir] = self._perpRange[vDir]:shorten(cd) -- shorten configuration range.
                  end
                  self._perpRange[vDir] = self._perpRange[vDir]:shorten(cDim+vDir) -- velocity range orthogonal to 'vDir'.
               end
               local perpRange = self._perpRange[vDir]

               for vPerpIdx in perpRange:rowMajorIter() do
                  vPerpIdx:copyInto(self.idxP)
                  for d = 1, cDim do self.idxP[d] = cIdx[d] end
      
                  for _, i in ipairs({dirLoIdx, dirUpIdx}) do     -- This loop is over edges.
                     self.idxP[cDim+vDir] = i
      
                     grid:setIndex(self.idxP)
                     grid:getDx(self.dxP)
                     grid:cellCenter(self.xcP)
      
                     distf:fill(phaseIndexer(self.idxP), distfItrP)
      
                     local vBound = 0.0
                     if isLo then
                        vBound = grid:cellLowerInDir(cDim + vDir)
                     else
                        vBound = grid:cellUpperInDir(cDim + vDir)
                     end
      
                     if (self._isGk) then
                        if (firstDir) then
                          uCorrection(isLo, self._intFac[1], vBound, self.dxP:data(), distfItrP:data(), cMomBItr:data())
                        end
                        vtSqCorrection(isLo, self._intFac[vDir], vBound, self.dxP:data(), distfItrP:data(), cEnergyBItr:data())
                     else
                        uCorrection(isLo, vBound, self.dxP:data(), distfItrP:data(), cMomBItr:data())
                        vtSqCorrection(isLo, vBound, self.dxP:data(), distfItrP:data(), cEnergyBItr:data())
                     end
      
                     isLo = not isLo
                  end
               end    -- vPerpIdx loop.
               firstDir = false
            end    -- vDir loop.
      
         end    -- if self._fiveMomentsLBO.
      end    -- Loop over configuration space.

   end    -- if self._fiveMomentsLBO and polyOrder=1.
   if self.momfac ~= 1.0 then mom1:scale(self.momfac) end
end

return DistFuncMomentCalc
