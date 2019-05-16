-- Gkyl ------------------------------------------------------------------------
--
-- Updater to create to fix a Maxwellian distribution function so it
-- matches the moments used to create it
--
--------------------------------------------------------------------------------

-- Gkyl libraries
local Lin = require "Lib.Linalg"
local LinearDecomp = require "Lib.LinearDecomp"
local Proto = require "Proto"
local Range = require "Lib.Range"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"
local ffi = require "ffi"

-- System libraries
local xsys = require "xsys"

-- Updater object
local LagrangeFix = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function LagrangeFix:init(tbl)
   LagrangeFix.super.init(self, tbl) -- setup base object

   self.confGrid = assert(tbl.confGrid,
			  "Updater.LagrangeFix: Must provide configuration space grid object using 'confGrid'")
   self.confBasis = assert(tbl.confBasis,
			   "Updater.LagrangeFix: Must provide configuration space basis object using 'confBasis'")
   self.phaseGrid = tbl.onGrid
   self.phaseBasis = assert(tbl.phaseBasis,
			    "Updater.LagrangeFix: Must provide phase space basis object using 'phaseBasis'")

   self.mode = assert(tbl.mode,
		      "Updater.LagrangeFix: Must specify 'mode'; options are: 'Vlasov' and 'Gk' for gyrokinetics")
   assert(self.mode == 'Vlasov' or self.mode == 'Gk',
	  "Updater.LagrangeFix: Supported options of 'mode' are: 'Vlasov' and 'Gk' for gyrokinetics")
   if self.mode == 'Gk' then
      self.mass = assert(tbl.mass, "Updater.LagrangeFix: Must provide 'mass' when in 'Gk' mode")
   end

   -- Dimension of spaces.
   self._pDim = self.phaseBasis:ndim() 
   self._cDim = self.confBasis:ndim()
   self._vDim = self._pDim - self._cDim

   self.basisID = self.confBasis:id()
   self.polyOrder = self.confBasis:polyOrder()
   -- Ensure sanity
   assert(self.phaseBasis:polyOrder() == self.polyOrder,
          "Polynomial orders of phase and conf basis must match")
   assert(self.phaseBasis:id() == self.basisID,
          "Polynomial orders of phase and conf basis must match")

   local LagrangeFixDecl = nil
   if self.mode == 'Vlasov' then
      LagrangeFixDecl = require "Updater.lagrangeFixData.VlasovLagrangeFixDecl"
   else
      LagrangeFixDecl = require "Updater.lagrangeFixData.GkLagrangeFixDecl"
   end
   self.lagrangeFixFn =
      LagrangeFixDecl.selectLagrangeFixFunction(self.basisID,
						self._cDim,
						self._vDim,
						self.polyOrder)

   -- Cell index, center, center velocity and upper velocity.
   self.idxP = Lin.IntVec(self._pDim)
   self.xcP = Lin.Vec(self._pDim)
   self.vcP = ffi.new("double[3]")
   self.Nv = ffi.new("double[3]")

   -- self.L = ffi.new("double[3]", {})
   -- self.lo = ffi.new("double[3]", {})
   self.L = ffi.new("double[3]")
   self.lo = ffi.new("double[3]")
   for d = 1, self._vDim do
      self.L[d-1] = self.phaseGrid:upper(self._cDim + d) - 
	 self.phaseGrid:lower(self._cDim + d) 
      self.lo[d-1] = self.phaseGrid:lower(self._cDim + d) 
   end
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function LagrangeFix:_advance(tCurr, inFld, outFld)
   -- Get the inputs and outputs
   local dm0 = assert(inFld[1],
		      "LagrangeFix.advance: Must specify dm0 as 'inFld[1]'")
   local dm1 = assert(inFld[2],
		      "LagrangeFix.advance: Must specify dm1 as 'inFld[2]'")
   local dm2 = assert(inFld[3],
		      "LagrangeFix.advance: Must specify dm2 as 'inFld[3]'")
   local B = nil
   if self.mode == 'Gk' then
      B = assert(inFld[4],
		 "LagrangeFix.advance: Must specify B as 'inFld[4]'")
   end
   local f = assert(outFld[1], "LagrangeFix.advance: Must specify an output distribution function")

   local cDim, vDim = self._cDim, self._vDim

   local dm0Itr = dm0:get(1)
   local dm1Itr = dm1:get(1)
   local dm2Itr = dm2:get(1)
   local BItr = nil
   if self.mode == 'Gk' then
      BItr = B:get(1)
   end
   local fItr = f:get(1)

   -- Get the Ranges to loop over the domain
   local confRange = dm0:localRange()
   local confIndexer = dm0:genIndexer()
   local phaseRange = f:localRange()
   local phaseIndexer = f:genIndexer()

   for d = 1, vDim do
      self.Nv[d-1] = phaseRange:upper(cDim + d)
   end

   -- construct ranges for nested loops
   local confRangeDecomp = LinearDecomp.LinearDecompRange {
      range = phaseRange:selectFirst(cDim), numSplit = self.phaseGrid:numSharedProcs() }
   local velRange = phaseRange:selectLast(vDim)
   local tId = self.phaseGrid:subGridSharedId() -- local thread ID

   -- The configuration space loop
   for cIdx in confRangeDecomp:rowMajorIter(tId) do
      dm0:fill(confIndexer(cIdx), dm0Itr)
      dm1:fill(confIndexer(cIdx), dm1Itr)
      dm2:fill(confIndexer(cIdx), dm2Itr)
      if self.mode == 'Gk' then
	 B:fill(confIndexer(cIdx), BItr)
      end

      -- The velocity space loop
      for vIdx in velRange:rowMajorIter() do
         cIdx:copyInto(self.idxP)
         for d = 1, vDim do self.idxP[cDim+d] = vIdx[d] end

	 f:fill(phaseIndexer(self.idxP), fItr)
	 self.phaseGrid:setIndex(self.idxP)
	 self.phaseGrid:cellCenter(self.xcP)
	 for d = 1, vDim do
	    self.vcP[d-1] = self.xcP[d + cDim]
	 end

	 if self.mode == 'Vlasov' then
	    self.lagrangeFixFn(dm0Itr:data(), dm1Itr:data(), dm2Itr:data(),
			       self.lo, self.L, self.Nv, self.vcP,
			       fItr:data())
	 else
	    self.lagrangeFixFn(dm0Itr:data(), dm1Itr:data(), dm2Itr:data(),
			       BItr:data(), self.mass, self.lo, self.L, self.Nv, self.vcP,
			       fItr:data())
	 end
      end
   end
end

return LagrangeFix
