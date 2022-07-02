-- Gkyl ------------------------------------------------------------------------
--
-- Updater to create the Maxwellian distribution from the conserved
-- moments and project it on basis functions. Uses Gaussian
-- quadrature.
--
-- Note: There's an implementation below that has more calculations in Lua.
--       This is slower but it is kept because as of 11/03/2020 it is the only
--       way to do quadrature with >polyOrder+1 quad points. We could add support
--       for such quadratures in the C implementation by generating more kernels.
--
--    _______     ___
-- + 6 @ |||| # P ||| +
--------------------------------------------------------------------------------

-- Gkyl libraries.
local GaussQuadRules    = require "Lib.GaussQuadRules"
local LinearDecomp      = require "Lib.LinearDecomp"
local Lin               = require "Lib.Linalg"
local Proto             = require "Proto"
local Range             = require "Lib.Range"
local Time              = require "Lib.Time"
local UpdaterBase       = require "Updater.Base"
local MaxwellianModDecl = require "Updater.maxwellianOnBasisData.MaxwellianOnBasisModDecl" 

-- System libraries.
local xsys = require "xsys"

-- Inherit the base Updater from UpdaterBase updater object.
local MaxwellianOnBasis = Proto(UpdaterBase)

function MaxwellianOnBasis:init(tbl)
   MaxwellianOnBasis.super.init(self, tbl) -- setup base object

   self.phaseGrid  = assert(tbl.onGrid,
                            "Updater.MaxwellianOnBasis: Must provide phase space grid object 'onGrid'")
   self.confGrid   = assert(tbl.confGrid,
                            "Updater.MaxwellianOnBasis: Must provide configuration space grid object 'confGrid'")
   self.confBasis  = assert(tbl.confBasis,
                            "Updater.MaxwellianOnBasis: Must provide configuration space basis object 'confBasis'")
   self.phaseBasis = assert(tbl.phaseBasis,
                            "Updater.MaxwellianOnBasis: Must provide phase space basis object 'phaseBasis'")

   -- Number of quadrature points in each direction
   local numQuad1D = tbl.numConfQuad and tbl.numConfQuad or self.confBasis:polyOrder() + 1
   assert(numQuad1D<=8, "Updater.MaxwellianOnBasis: Gaussian quadrature only implemented for numQuad<=8 in each dimension")

   self.quadType = "Gauss"
   self.quadImpl = tbl.implementation and tbl.implementation or "C"

   -- The C implementation only allows p+1 quadrature points (in 1D).
   if numQuad1D ~= self.confBasis:polyOrder()+1 then self.quadImpl="Lua" end

   self.onGhosts = xsys.pickBool(tbl.onGhosts, false)

   -- Mass and magnetic field amplitude are needed to project a gyrokinetic Maxwellian.
   self.mass = tbl.mass

   self._pDim, self._cDim = self.phaseBasis:ndim(), self.confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   if self.mass then
      self.isGK = true
   else
      self.isGK = false
      self.mass = 0.0    -- Not used.
   end

   if self.quadImpl == "C" then

      self.numConfOrds = numQuad1D^self._cDim 

      if self.isGK then
         self.uFlowOrd = Lin.Vec(self.numConfOrds)
      else
         self.uFlowOrd = Lin.Vec(self.numConfOrds*self._vDim)
      end
      self.vtSqOrd    = Lin.Vec(self.numConfOrds)
      self.normFacOrd = Lin.Vec(self.numConfOrds)
      self.bmagOrd    = Lin.Vec(self.numConfOrds)

      local quadFuncs = MaxwellianModDecl.selectQuad(self.confBasis:id(), self._cDim, self._vDim, self.confBasis:polyOrder(), self.quadType, self.isGK)
      self._evAtConfOrds = quadFuncs[1]
      self._phaseQuad    = quadFuncs[2]

   elseif self.quadImpl == "Lua" then

      self._innerLoop = MaxwellianModDecl.selectInnerLoop(self.isGK)
      -- 1D weights and ordinates.
      local ordinates = GaussQuadRules.ordinates[numQuad1D]
      local weights   = GaussQuadRules.weights[numQuad1D]

      -- Configuration space ordinates ----------------------------------
      local l, u = {}, {}
      for d = 1, self._cDim do l[d], u[d] = 1, numQuad1D end
      self.confQuadRange = Range.Range(l, u)
      self.confQuadIndexer = 
         Range.makeColMajorGenIndexer(self.confQuadRange)
      self.numConfOrds = self.confQuadRange:volume()
      self.numConfBasis = self.confBasis:numBasis()
      self.confOrdinates = Lin.Mat(self.numConfOrds, self._cDim)
      self.confBasisAtOrds = Lin.Mat(self.numConfOrds, self.numConfBasis)
      for ordIndexes in self.confQuadRange:colMajorIter() do
         local ordIdx = self.confQuadIndexer(ordIndexes)
         for d = 1, self._cDim do
            self.confOrdinates[ordIdx][d] = ordinates[ordIndexes[d]]
         end
         self.confBasis:evalBasis(self.confOrdinates:getRow(ordIdx),
                                  self.confBasisAtOrds:getRow(ordIdx))
      end

      -- Phase space ordinates and weights ------------------------------
      for d = 1, self._pDim do l[d], u[d] = 1, numQuad1D end
      self.phaseQuadRange = Range.Range(l, u)
      self.phaseQuadIndexer = 
         Range.makeColMajorGenIndexer(self.phaseQuadRange)
      self.numPhaseOrds = self.phaseQuadRange:volume()
      self.numPhaseBasis = self.phaseBasis:numBasis()
      self.phaseBasisAtOrds = Lin.Mat(self.numPhaseOrds,
           			   self.numPhaseBasis)
      self.phaseOrdinates = Lin.Mat(self.numPhaseOrds, self._pDim)
      self.phaseWeights = Lin.Vec(self.numPhaseOrds) -- Needed for integration.
      for ordIndexes in self.phaseQuadRange:colMajorIter() do
         local ordIdx = self.phaseQuadIndexer(ordIndexes)
         self.phaseWeights[ordIdx] = 1.0
         for d = 1, self._pDim do
            self.phaseWeights[ordIdx] =
               self.phaseWeights[ordIdx]*weights[ordIndexes[d]]
            self.phaseOrdinates[ordIdx][d] = ordinates[ordIndexes[d]]
         end
         self.phaseBasis:evalBasis(self.phaseOrdinates:getRow(ordIdx),
                                   self.phaseBasisAtOrds:getRow(ordIdx))
      end

      -- Construct the phase space to conf space ordinate map.
      self.phaseToConfOrdMap = Lin.Vec(self.numPhaseOrds)
      for ordIndexes in self.phaseQuadRange:colMajorIter() do
         local confOrdIdx = self.confQuadIndexer(ordIndexes)
         local phaseOrdIdx = self.phaseQuadIndexer(ordIndexes)
         self.phaseToConfOrdMap[phaseOrdIdx] = confOrdIdx
      end

   end

   -- Cell index, center, and dx.
   self.idx = Lin.IntVec(self._pDim)
   self.xc  = Lin.Vec(self._pDim)
   self.dx  = Lin.Vec(self._pDim)
end

function MaxwellianOnBasis:_advance(tCurr, inFld, outFld)
   -- Get the inputs and outputs.
   local nIn     = assert(inFld[1], "MaxwellianOnBasis: Must specify density 'inFld[1]'")
   local uFlowIn = assert(inFld[2], "MaxwellianOnBasis: Must specify drift speed 'inFld[2]'")
   local vtSqIn  = assert(inFld[3], "MaxwellianOnBasis: Must specify thermal velocity squared 'inFld[3]'")
   local fOut    = assert(outFld[1], "MaxwellianOnBasis: Must specify an output field 'outFld[1]'")

   local pDim, cDim, vDim = self._pDim, self._cDim, self._vDim
   local uDim = uFlowIn:numComponents()/self.confBasis:numBasis()   -- Number of vector components of uFlowIn.

   local nItr, uFlowItr, vtSqItr = nIn:get(1), uFlowIn:get(1), vtSqIn:get(1)
   local fItr = fOut:get(1)

   local uKer = 1   -- This kernel label indicates whether to take in a gyrokinetic or vlasov drift velocity.
   local bmagIn, bmagItr
   if self.isGK then
      if uDim>1 then uKer=2 end
      bmagIn  = assert(inFld[4], "MaxwellianOnBasis: Must specify magnetic field amplitude 'inFld[4]'")
      bmagItr = bmagIn:get(1)
   else
      if uDim~=vDim then uKer=2 end
      bmagItr = nIn:get(1)   -- Not used.
   end

   local confIndexer  = nIn:genIndexer()
   local phaseIndexer = fOut:genIndexer()
   local phaseRange
   if self.onGhosts then
      phaseRange = fOut:localExtRange()
   else
      phaseRange = fOut:localRange()
   end

   -- Construct ranges for nested loops.
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(cDim), numSplit = self.phaseGrid:numSharedProcs() }
   local velRange = phaseRange:selectLast(vDim)
   local tId      = self.phaseGrid:subGridSharedId()   -- Local thread ID.

   if self.quadImpl == "C" then

      -- The configuration space loop
      for cIdx in confRangeDecomp:rowMajorIter(tId) do
         nIn:fill(confIndexer(cIdx), nItr)
         uFlowIn:fill(confIndexer(cIdx), uFlowItr)
         vtSqIn:fill(confIndexer(cIdx), vtSqItr)
         if self.isGK then bmagIn:fill(confIndexer(cIdx), bmagItr) end

         self._evAtConfOrds[uKer](nItr:data(), uFlowItr:data(), vtSqItr:data(), bmagItr:data(),
                                  self.uFlowOrd:data(), self.vtSqOrd:data(), self.normFacOrd:data(), self.bmagOrd:data())

         -- The velocity space loop
         for vIdx in velRange:rowMajorIter() do
            cIdx:copyInto(self.idx)
            for d = 1, vDim do self.idx[cDim+d] = vIdx[d] end

   	    fOut:fill(phaseIndexer(self.idx), fItr)
   
   	    -- Get cell shape, cell center coordinates
   	    self.phaseGrid:setIndex(self.idx)
            self.phaseGrid:getDx(self.dx)
            self.phaseGrid:cellCenter(self.xc)
   
            self._phaseQuad(self.uFlowOrd:data(), self.vtSqOrd:data(), self.normFacOrd:data(),
                            self.bmagOrd:data(), self.mass,
                            self.xc:data(), self.dx:data(), fItr:data())

         end
      end

   elseif self.quadImpl == "Lua" then

      local nOrd     = Lin.Vec(self.numConfOrds)
      local uFlowOrd = self.isGK and Lin.Vec(self.numConfOrds) or Lin.Mat(self.numConfOrds, vDim)
      local vtSqOrd  = Lin.Vec(self.numConfOrds)
      local bmagOrd  = Lin.Vec(self.numConfOrds)
   
      -- Additional preallocated variables
      local ordIdx = nil
   
      -- The configuration space loop
      for cIdx in confRangeDecomp:rowMajorIter(tId) do
         nIn:fill(confIndexer(cIdx), nItr)
         uFlowIn:fill(confIndexer(cIdx), uFlowItr)
         vtSqIn:fill(confIndexer(cIdx), vtSqItr)
         if self.isGK then bmagIn:fill(confIndexer(cIdx), bmagItr) end
   
         -- Evaluate the the primitive variables (given as expansion
         -- coefficients) on the ordinates
         for ordIndexes in self.confQuadRange:rowMajorIter() do
            ordIdx = self.confQuadIndexer(ordIndexes)
            nOrd[ordIdx], vtSqOrd[ordIdx] = 0.0, 0.0
           
            for k = 1, self.numConfBasis do
               nOrd[ordIdx]    = nOrd[ordIdx] + nItr[k]*self.confBasisAtOrds[ordIdx][k]
               vtSqOrd[ordIdx] = vtSqOrd[ordIdx] + vtSqItr[k]*self.confBasisAtOrds[ordIdx][k]
            end
            if self.isGK then 
               uFlowOrd[ordIdx], bmagOrd[ordIdx] = 0.0, 0.0
               for k = 1, self.numConfBasis do
                  bmagOrd[ordIdx] = bmagOrd[ordIdx] + bmagItr[k]*self.confBasisAtOrds[ordIdx][k]
                  if uDim > 1 and cDim > 1 then   -- Get the z-component of fluid velocity.
                     uFlowOrd[ordIdx] = uFlowOrd[ordIdx] + uFlowItr[self.numConfBasis*2+k]*self.confBasisAtOrds[ordIdx][k]
                  else
                     uFlowOrd[ordIdx] = uFlowOrd[ordIdx] + uFlowItr[k]*self.confBasisAtOrds[ordIdx][k]
                  end
               end
            else
               for d = 1, vDim do uFlowOrd[ordIdx][d] = 0.0 end
               if uDim == vDim then
                  for d = 1, vDim do
                     for k = 1, self.numConfBasis do
                        uFlowOrd[ordIdx][d] = uFlowOrd[ordIdx][d] +
                           uFlowItr[self.numConfBasis*(d-1)+k]*self.confBasisAtOrds[ordIdx][k]
                     end
                  end
               elseif uDim == 1 and cDim==1 and vDim==3 then -- if uPar passed from GkSpecies, fill d=1 component of u
                  for k = 1, self.numConfBasis do
                     uFlowOrd[ordIdx][1] = uFlowOrd[ordIdx][1] + uFlowItr[k]*self.confBasisAtOrds[ordIdx][k]
                  end
               elseif uDim == 1 and cDim>1 and vDim==3 then -- if uPar passed from GkSpecies, fill d=3 component of u
                  for k = 1, self.numConfBasis do
                     uFlowOrd[ordIdx][vDim] = uFlowOrd[ordIdx][vDim] + uFlowItr[k]*self.confBasisAtOrds[ordIdx][k]
                  end
               else
                  print("Updater.MaxwellianOnBasis: incorrect uDim")	 
               end
            end
         end
   
         -- The velocity space loop
         for vIdx in velRange:rowMajorIter() do
            -- Construct the phase space index ot of the configuration
            -- space and velocity space indices.
            cIdx:copyInto(self.idx)
            for d = 1, vDim do self.idx[cDim+d] = vIdx[d] end
   	    fOut:fill(phaseIndexer(self.idx), fItr)
   
   	    -- Get cell shape, cell center coordinates
   	    self.phaseGrid:setIndex(self.idx)
            self.phaseGrid:getDx(self.dx)
            self.phaseGrid:cellCenter(self.xc)
   
            self._innerLoop(nOrd:data(), uFlowOrd:data(), vtSqOrd:data(),
                            bmagOrd:data(), self.mass,
                            fItr:data(),
                            self.phaseWeights:data(), self.dx:data(), self.xc:data(),
                            self.phaseOrdinates:data(),
                            self.phaseBasisAtOrds:data(),
                            self.phaseToConfOrdMap:data(),
                            self.numPhaseBasis,
                            self.numConfOrds, self.numPhaseOrds,
                            cDim, pDim)
         end
      end
   
   end
   -- Set id of output to id of projection basis.
   fOut:setBasisId(self.phaseBasis:id())
end

return MaxwellianOnBasis
