-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate the primitive moments, u and vtSq=sqrt(T/m), given
-- the distribution function and its moments.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries
local UpdaterBase     = require "Updater.Base"
local Lin             = require "Lib.Linalg"
local Proto           = require "Lib.Proto"
local PrimMomentsDecl = require "Updater.primMomentsCalcData.PrimMomentsModDecl"
local xsys            = require "xsys"
local ffi = require "ffi"
local ffiC = ffi.C

ffi.cdef [[
  void gkylCopyToField(double *f, double *data, unsigned numComponents, unsigned c);
  void gkylCartFieldAssignAll(unsigned s, unsigned nv, double val, double *out);
]]

-- function to check if operator option is correct
local function isOperatorGood(nm)
   if nm == "VmLBO" or nm == "GkLBO" then
      return true
   end
   return false
end

-- Moments updater object.
local SelfPrimMoments = Proto(UpdaterBase)

function SelfPrimMoments:init(tbl)
   SelfPrimMoments.super.init(self, tbl) -- setup base object.

   self._onGrid = assert(
      tbl.onGrid, "Updater.SelfPrimMoments: Must provide grid object using 'onGrid'.")

   self._phaseGrid = assert(
      tbl.phaseGrid, "Updater.SelfPrimMoments: Must provide phase grid object using 'phaseGrid'.")

   local phaseBasis = assert(
      tbl.phaseBasis, "Updater.SelfPrimMoments: Must provide the phase basis object using 'phaseBasis'.")

   local confBasis = assert(
      tbl.confBasis, "Updater.SelfPrimMoments: Must provide the configuration basis object using 'confBasis'.")

   local operator = assert(
      tbl.operator, "Updater.SelfPrimMoments: Must specify the collision operator (VmLBO or GkLBO for now) using 'operator'.")

   -- dimension of phase space.
   self._pDim = phaseBasis:ndim()
   -- Basis name and polynomial order.
   self._basisID   = confBasis:id()
   self._polyOrder = confBasis:polyOrder()

   -- ensure sanity.
   assert(phaseBasis:polyOrder() == self._polyOrder,
          "Polynomial orders of phase and conf basis must match.")
   assert(phaseBasis:id() == self._basisID,
          "Type of phase and conf basis must match.")
   -- determine configuration and velocity space dims.
   self._cDim = confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   -- Number of basis functions. Used to compute number of vector components.
   self._numBasisP = phaseBasis:numBasis()
   self._numBasisC = confBasis:numBasis()

   self._uDim       = self._vDim  -- Dimensionality of flow velocity vector.
   self._kinSpecies = "Vm"        -- Vlasov-Maxwell species.
   if isOperatorGood(operator) then
     self._isGkLBO = false
     if operator == "GkLBO" then
        self._isGkLBO    = true
        self._kinSpecies = "Gk"  -- Gyrokinetic species.

        assert(tbl.gkfacs, [[DistFuncMomentCalc: must provide a gkfacs table
                            containing the species mass and the background magnetic field
                            to calculate a Gk moment]])
        self.mass = tbl.gkfacs[1]
        self.bmag = assert(tbl.gkfacs[2], "DistFuncMomentCalc: must provide bmag in gkfacs")

        self.bmagIndexer = self.bmag:genIndexer()
        self.bmagItr     = self.bmag:get(1)

        -- If vDim>1, intFac=2*pi/m or 4*pi/m.
        self._intFac       = Lin.Vec(self._vDim)
        for d = 1,self._vDim do 
          self._intFac[d]  = 1.0
        end
        if self._vDim > 1 then -- A (vpar,mu) simulation has 3 physical velocity dimensions.
           self._intFac[1] = 2.0*math.pi/self.mass
           self._intFac[2] = 4.0*math.pi/self.mass
           self._uDim      = 1
        end
     end
   else
      assert(false, string.format(
                "SelfPrimMoments: Operator option must be 'VmLBO' or 'GkLBO'. Requested %s instead.", operator))
   end

   self._SelfPrimMomentsCalc = PrimMomentsDecl.selectSelfPrimMomentsCalc(self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)
   self.onGhosts = xsys.pickBool(true, tbl.onGhosts)
   
   self._isFirst   = true
   self._perpRange = {} -- perp ranges in velocity directions.
   if self._polyOrder == 1 then
      self._StarM1iM2Calc = PrimMomentsDecl.selectStarM1iM2Calc(self._kinSpecies, self._basisID, self._cDim, self._vDim) 

   end

   self._binOpData = ffiC.new_binOpData_t(self._numBasisC*(self._uDim+1), 0) 
end

-- advance method
function SelfPrimMoments:_advance(tCurr, inFld, outFld)
   local grid = self._onGrid

   local uOut           = outFld[1]
   local vtSqOut        = outFld[2]

   local uOutIndexer    = uOut:genIndexer()
   local uOutItr        = uOut:get(1)
   local vtSqOutIndexer = vtSqOut:genIndexer()
   local vtSqOutItr     = vtSqOut:get(1)

   local m0fld, m1fld, fIn = inFld[1], inFld[2], inFld[4]

   local m0fldIndexer   = m0fld:genIndexer()
   local m0fldItr       = m0fld:get(1)
   local m1fldIndexer   = m1fld:genIndexer()
   local m1fldItr       = m1fld:get(1)
   local phaseIndexer   = fIn:genIndexer()
   local fInItrP         = fIn:get(1)
   local fInItrM         = fIn:get(1)

   local m0Star, m1Star, m2Star
   local m2fld, m2fldIndexer, m2fldItr
   if self._polyOrder == 1 then
      m0Star = Lin.Vec(self._numBasisC)
      m1Star = Lin.Vec(self._numBasisC*self._uDim)
      m2Star = Lin.Vec(self._numBasisC)
   else
      m2fld         = inFld[3]
      m2fldIndexer  = m2fld:genIndexer()
      m2fldItr      = m2fld:get(1)
   end

   local confRange  = m0fld:localRange()
   if self.onGhosts then confRange = m0fld:localExtRange() end

   local phaseRange = fIn:localRange()

   -- Distribution functions left and right of a cell-boundary.
   local fInP = Lin.Vec(self._numBasisP)
   local fInM = Lin.Vec(self._numBasisP)

   -- These store the corrections to momentum and energy
   -- from the (velocity) boundary surface integrals.
   local cMomB    = Lin.Vec(self._numBasisC*self._uDim)
   local cEnergyB = Lin.Vec(self._numBasisC)

   -- Cell index, center and length left and right of a cell-boundary.
   local idxM, idxP = Lin.IntVec(self._pDim), Lin.IntVec(self._pDim)
   local xcM, xcP   = Lin.Vec(self._pDim), Lin.Vec(self._pDim)
   local dxM, dxP   = Lin.Vec(self._pDim), Lin.Vec(self._pDim)

   for confIdx in confRange:rowMajorIter() do
      grid:setIndex(confIdx)

      m0fld:fill(m0fldIndexer(confIdx), m0fldItr)
      m1fld:fill(m1fldIndexer(confIdx), m1fldItr)
      uOut:fill(uOutIndexer(confIdx), uOutItr)
      vtSqOut:fill(vtSqOutIndexer(confIdx), vtSqOutItr)
      if (self._isGkLBO) then
         self.bmag:fill(self.bmagIndexer(confIdx), self.bmagItr)
      end

      -- Compute the corrections to u and vtSq due to 
      -- finite velocity grid to ensure conservation.
      ffiC.gkylCartFieldAssignAll(0, self._numBasisC, 0.0, cEnergyB:data())
      ffiC.gkylCartFieldAssignAll(0, self._uDim*self._numBasisC, 0.0, cMomB:data())

      -- Only when the contributions to m0Star from the first direction
      -- are collected, do we collect contributions to m1Star and m2Star.
      -- Also, since Gk velocities are organized as (vpar,mu) the velocity
      -- correction is only computed for the first velocity direction.
      local firstDir = true

      -- isLo=true current cell is the lower boundary cell.
      -- isLo=false current cell is the upper boundary cell.
      local isLo = true

      -- MF: Apologies if the code seems lengthy. 
      -- polyOrder=1 and >1 each use separate velocity grid loops to
      -- avoid evaluating (if polyOrder==1) at each velocity coordinate.
      if self._polyOrder > 1 then
         m2fld:fill(m2fldIndexer(confIdx), m2fldItr)

         for vDir = 1, self._vDim do
            if (not self._isGkLBO) or (self._isGkLBO and firstDir) then
               self._uCorrection    = PrimMomentsDecl.selectBoundaryFintegral(vDir, self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)
            end
            self._vtSqCorrection = PrimMomentsDecl.selectBoundaryVFintegral(vDir, self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)

            -- Lower/upper bounds in direction 'vDir': cell indices.
            local dirLoIdx, dirUpIdx = phaseRange:lower(self._cDim+vDir), phaseRange:upper(self._cDim+vDir)

            if self._isFirst then
               self._perpRange[vDir] = phaseRange
               for cd = 1, self._cDim do
                  self._perpRange[vDir] = self._perpRange[vDir]:shorten(cd) -- shorten configuration range.
               end
               self._perpRange[vDir] = self._perpRange[vDir]:shorten(self._cDim+vDir) -- velocity range orthogonal to 'vDir'.
            end
            local perpRange = self._perpRange[vDir]

            for vPerpIdx in perpRange:rowMajorIter() do
               vPerpIdx:copyInto(idxP)
               for d = 1, self._cDim do idxP[d] = confIdx[d] end

               for _, i in ipairs({dirLoIdx, dirUpIdx}) do     -- This loop is over edges.
                  idxP[self._cDim+vDir] = i 

                  self._phaseGrid:setIndex(idxP)
                  self._phaseGrid:getDx(dxP)
                  self._phaseGrid:cellCenter(xcP)

                  fIn:fill(phaseIndexer(idxP), fInItrP)

                  local vBound = 0.0
                  if isLo then
                     vBound = self._phaseGrid:cellLowerInDir(self._cDim + vDir)
                  else
                     vBound = self._phaseGrid:cellUpperInDir(self._cDim + vDir)
                  end

                  if (self._isGkLBO) then
                     if (firstDir) then
                       self._uCorrection(isLo, self._intFac[1], vBound, dxP:data(), fInItrP:data(), cMomB:data())
                     end
                     self._vtSqCorrection(isLo, self._intFac[vDir], vBound, dxP:data(), fInItrP:data(), cEnergyB:data())
                  else
                     self._uCorrection(isLo, vBound, dxP:data(), fInItrP:data(), cMomB:data())
                     self._vtSqCorrection(isLo, vBound, dxP:data(), fInItrP:data(), cEnergyB:data())
                  end

                  isLo = not isLo
               end
            end
            firstDir = false
         end

         self._SelfPrimMomentsCalc(self._binOpData, m0fldItr:data(), m1fldItr:data(), m2fldItr:data(), cMomB:data(), cEnergyB:data(), uOutItr:data(), vtSqOutItr:data())

      else
         -- To have energy conservation with piece-wise linear, we must use
         -- star moments in the second equation of the weak system solved
         -- in SelfPrimMoments.
         ffiC.gkylCartFieldAssignAll(0, self._numBasisC, 0.0, m0Star:data())
         ffiC.gkylCartFieldAssignAll(0, self._uDim*self._numBasisC, 0.0, m1Star:data())
         ffiC.gkylCartFieldAssignAll(0, self._numBasisC, 0.0, m2Star:data())

         for vDir = 1, self._vDim do
            if (not self._isGkLBO) or (self._isGkLBO and firstDir) then
               self._StarM0Calc = PrimMomentsDecl.selectStarM0Calc(vDir, self._kinSpecies, self._basisID, self._cDim, self._vDim)
               self._uCorrection    = PrimMomentsDecl.selectBoundaryFintegral(vDir, self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)
            end
            self._vtSqCorrection = PrimMomentsDecl.selectBoundaryVFintegral(vDir, self._kinSpecies, self._basisID, self._cDim, self._vDim, self._polyOrder)

            -- Lower/upper bounds in direction 'vDir': edge indices (including outer edges).
            local dirLoIdx, dirUpIdx = phaseRange:lower(self._cDim+vDir), phaseRange:upper(self._cDim+vDir)+1

            if self._isFirst then
               -- Restricted velocity range.
               -- Velocity integral in m0Star does not include last cell.
               self._perpRange[vDir] = phaseRange
               for cd = 1, self._cDim do
                  self._perpRange[vDir] = self._perpRange[vDir]:shorten(cd) -- shorten configuration range.
               end
               self._perpRange[vDir] = self._perpRange[vDir]:shorten(self._cDim+vDir) -- velocity range orthogonal to 'vDir'.
            end
            local perpRange = self._perpRange[vDir]

            -- Outer loop is over directions orthogonal to 'vDir' and
            -- inner loop is over 1D slice in 'vDir'.
            for vPerpIdx in perpRange:rowMajorIter() do
               vPerpIdx:copyInto(idxM); vPerpIdx:copyInto(idxP)
               for d = 1, self._cDim do idxM[d] = confIdx[d] end
               for d = 1, self._cDim do idxP[d] = confIdx[d] end

               for i = dirLoIdx, dirUpIdx do     -- This loop is over edges.
                  idxM[self._cDim+vDir], idxP[self._cDim+vDir] = i-1, i -- Cell left/right of edge 'i'.

	          self._phaseGrid:setIndex(idxM)
                  self._phaseGrid:getDx(dxM)
	          self._phaseGrid:cellCenter(xcM)

	          self._phaseGrid:setIndex(idxP)
                  self._phaseGrid:getDx(dxP)
	          self._phaseGrid:cellCenter(xcP)

                  fIn:fill(phaseIndexer(idxM), fInItrM)
                  fIn:fill(phaseIndexer(idxP), fInItrP)

                  if i>dirLoIdx and i<dirUpIdx then  
                     if (self._isGkLBO) then
                        if (firstDir) then
                           self._StarM0Calc(self._intFac[1], xcM:data(), xcP:data(), dxM:data(), dxP:data(), fInItrM:data(), fInItrP:data(), m0Star:data())
                        end
                     else
                        self._StarM0Calc(xcM:data(), xcP:data(), dxM:data(), dxP:data(), fInItrM:data(), fInItrP:data(), m0Star:data())
                     end
                  end
             	  if firstDir and i<dirUpIdx then
                     if self._isGkLBO then
                        self._StarM1iM2Calc(xcP:data(), dxP:data(), self._intFac[1], self.mass, self.bmagItr:data(), fInItrP:data(), m1Star:data(), m2Star:data())
                     else
                        self._StarM1iM2Calc(xcP:data(), dxP:data(), fInItrP:data(), m1Star:data(), m2Star:data())
                     end
                  end

                  if i==dirLoIdx or i==dirUpIdx-1 then
                     local vBound = 0.0
                     -- Careful: for vBound below we assume idxP was set after idxM above.
                     if isLo then
                        vBound = self._phaseGrid:cellLowerInDir(self._cDim + vDir)
                     else
                        vBound = self._phaseGrid:cellUpperInDir(self._cDim + vDir)
                     end
                     if (self._isGkLBO) then
                        if (firstDir) then
                           self._uCorrection(isLo, self._intFac[1], vBound, dxP:data(), fInItrP:data(), cMomB:data())
                        end
                        self._vtSqCorrection(isLo, self._intFac[vDir], vBound, dxP:data(), fInItrP:data(), cEnergyB:data())
                     else
                        self._uCorrection(isLo, vBound, dxP:data(), fInItrP:data(), cMomB:data())
                        self._vtSqCorrection(isLo, vBound, dxP:data(), fInItrP:data(), cEnergyB:data())
                     end

                     isLo = not isLo
                  end

               end
            end
            firstDir = false
         end

      self._SelfPrimMomentsCalc(self._binOpData, m0fldItr:data(), m1fldItr:data(), m0Star:data(), m1Star:data(), m2Star:data(), cMomB:data(), cEnergyB:data(), uOutItr:data(), vtSqOutItr:data())

      end

   end
   self._isFirst = false
end

return SelfPrimMoments
