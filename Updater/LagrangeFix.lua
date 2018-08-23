-- Gkyl ------------------------------------------------------------------------
--
-- Updater to create to fix a Maxwellian distribution function so it
-- matches the moments used to create it
--
--------------------------------------------------------------------------------

-- Gkyl libraries
local LagrangeFixDecl = require "Updater.lagrangeFixData.lagrangeFixDecl"
local Lin = require "Lib.Linalg"
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
   self.phaseGrid = assert(tbl.phaseGrid,
			   "Updater.LagrangeFix: Must provide phase space grid object using 'phaseGrid'")
   self.phaseBasis = assert(tbl.phaseBasis,
			    "Updater.LagrangeFix: Must provide phase space basis object using 'phaseBasis'")


   self.numConfDims = self.confBasis:ndim()
   self.numPhaseDims = self.phaseBasis:ndim()
   self.numVelDims = self.phaseGrid:ndim() - self.confGrid:ndim()

   self.basisID = self.confBasis:id()
   self.polyOrder = self.confBasis:polyOrder()
   -- Ensure sanity
   assert(self.phaseBasis:polyOrder() == self.polyOrder,
          "Polynomial orders of phase and conf basis must match")
   assert(self.phaseBasis:id() == self.basisID,
          "Polynomial orders of phase and conf basis must match")

   self.lagrangeFixFn =
      LagrangeFixDecl.selectLagrangeFixFunction(self.basisID,
						self.numConfDims,
						self.numVelDims,
						self.polyOrder)

   self.L = ffi.new("double[3]", {})
   for d = 1, self.numVelDims do
      self.L[d-1] = self.phaseGrid:upper(self.numConfDims + d) - 
	 self.phaseGrid:lower(self.numConfDims + d) 
   end
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function LagrangeFix:_advance(tCurr, dt, inFld, outFld)
   -- Get the inputs and outputs
   local dm0 = assert(inFld[1],
		      "LagrangeFix.advance: Must specify dm0 as 'inFld[1]")
   local dm1 = assert(inFld[2],
		      "LagrangeFix.advance: Must specify dm1 as 'inFld[2]")
   local dm2 = assert(inFld[3],
		      "LagrangeFix.advance: Must specify dm2 as 'inFld[3]")
   local f = assert(outFld[1], "LagrangeFix.advance: Must specify an output distribution function")

   local dm0Itr = dm0:get(1)
   local dm1Itr = dm1:get(1)
   local dm2Itr = dm2:get(1)
   local fItr = f:get(1)

   -- Get the Ranges to loop over the domain
   local confRange = dm0:localRange()
   local confIndexer = dm0:genIndexer()
   local phaseRange = f:localRange()
   local phaseIndexer = f:genIndexer()
   local phaseIdx = {}
   local l, u = {}, {}
   local Nv = ffi.new("double[3]")
   for d = 1, self.numVelDims do
      l[d] = phaseRange:lower(self.numConfDims + d)
      u[d] = phaseRange:upper(self.numConfDims + d)
      Nv[d-1] = phaseRange:upper(self.numConfDims + d)
   end
   local velRange = Range.Range(l, u)
   local zc = Lin.Vec(self.numPhaseDims) -- cell center
   local vc = ffi.new("double[3]") -- cell center

   -- The configuration space loop
   for confIdx in confRange:colMajorIter() do
      dm0:fill(confIndexer(confIdx), dm0Itr)
      dm1:fill(confIndexer(confIdx), dm1Itr)
      dm2:fill(confIndexer(confIdx), dm2Itr)

      -- The velocity space loop
      for velIdx in velRange:colMajorIter() do
	 -- Construct the phase space index ot of the configuration
	 -- space a velocity space indices
	 for d = 1, self.numConfDims do phaseIdx[d] = confIdx[d] end
	 for d = 1, self.numVelDims do
	    phaseIdx[d + self.numConfDims] = velIdx[d]
	 end
	 f:fill(phaseIndexer(phaseIdx), fItr)
	 self.phaseGrid:setIndex(phaseIdx)
	 self.phaseGrid:cellCenter(zc)
	 for d = 1, self.numVelDims do
	    vc[d-1] = zc[d + self.numConfDims]
	 end

	 self.lagrangeFixFn(dm0Itr:data(), dm1Itr:data(), dm2Itr:data(),
			    self.L, Nv, vc,
			    fItr:data())
      end
   end
   return true, GKYL_MAX_DOUBLE
end

return LagrangeFix
