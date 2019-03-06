-- Gkyl ------------------------------------------------------------------------
--
-- Updater to create the Maxwellian distribution from the coserved
-- moments and project it on basis functions. Uses Gaussian
-- quadrature.
--
--------------------------------------------------------------------------------

-- Gkyl libraries
local GaussQuadRules = require "Lib.GaussQuadRules"
local Lin = require "Lib.Linalg"
local Proto = require "Proto"
local Range = require "Lib.Range"
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"

-- System libraries
local ffi = require "ffi"
local ffiC = ffi.C
local xsys = require "xsys"

ffi.cdef [[
  void MaxwellianInnerLoop(double * n, double * u, double * vth2,
			   double * fItr,
			   double * weights, double * dz, double * zc,
			   double * ordinates,
                           double * basisAtOrdinates,
                           double * phaseToConfOrdMap,
			   int numPhaseBasis, 
                           int numConfOrds, int numPhaseOrds,
			   int numConfDims, int numPhaseDims);
]]

-- Inherit the base Updater from UpdaterBase updater object
local MaxwellianOnBasis = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function MaxwellianOnBasis:init(tbl)
   MaxwellianOnBasis.super.init(self, tbl) -- setup base object

   self.confGrid = assert(tbl.confGrid,
			  "Updater.MaxwellianOnBasis: Must provide configuration space grid object 'confGrid'")
   self.confBasis = assert(tbl.confBasis,
			   "Updater.MaxwellianOnBasis: Must provide configuration space basis object 'confBasis'")
   self.phaseGrid = assert(tbl.onGrid,
			   "Updater.MaxwellianOnBasis: Must provide phase space grid object 'onGrid'")
   self.phaseBasis = assert(tbl.phaseBasis,
			    "Updater.MaxwellianOnBasis: Must provide phase space basis object 'phaseBasis'")

   -- Number of quadrature points in each direction
   local N = tbl.numConfQuad and tbl.numConfQuad or self.confBasis:polyOrder() + 1
   assert(N<=8, "Updater.MaxwellianOnBasis: Gaussian quadrature only implemented for numQuad<=8 in each dimension")

   self.projectOnGhosts = xsys.pickBool(tbl.projectOnGhosts, false)

   -- 1D weights and ordinates
   local ordinates = GaussQuadRules.ordinates[N]
   local weights = GaussQuadRules.weights[N]

   -- Configuration space ordinates ----------------------------------
   self.numConfDims = self.confBasis:ndim()
   local l, u = {}, {}
   for d = 1, self.numConfDims do l[d], u[d] = 1, N end
   self.confQuadRange = Range.Range(l, u)
   self.confQuadIndexer = 
      Range.makeColMajorGenIndexer(self.confQuadRange)
   self.numConfOrds = self.confQuadRange:volume()
   self.numConfBasis = self.confBasis:numBasis()
   self.confOrdinates = Lin.Mat(self.numConfOrds, self.numConfDims)
   self.confBasisAtOrds = Lin.Mat(self.numConfOrds, self.numConfBasis)
   for ordIndexes in self.confQuadRange:colMajorIter() do
      local ordIdx = self.confQuadIndexer(ordIndexes)
      for d = 1, self.numConfDims do
	 self.confOrdinates[ordIdx][d] = ordinates[ordIndexes[d]]
      end
      self.confBasis:evalBasis(self.confOrdinates[ordIdx],
			       self.confBasisAtOrds[ordIdx])
   end

   -- Phase space ordinates and weights ------------------------------
   self.numPhaseDims = self.phaseBasis:ndim()
   for d = 1, self.numPhaseDims do l[d], u[d] = 1, N end
   self.phaseQuadRange = Range.Range(l, u)
   self.phaseQuadIndexer = 
      Range.makeColMajorGenIndexer(self.phaseQuadRange)
   self.numPhaseOrds = self.phaseQuadRange:volume()
   self.numPhaseBasis = self.phaseBasis:numBasis()
   self.phaseBasisAtOrds = Lin.Mat(self.numPhaseOrds,
				   self.numPhaseBasis)
   self.phaseOrdinates = Lin.Mat(self.numPhaseOrds, self.numPhaseDims)
   self.phaseWeights = Lin.Vec(self.numPhaseOrds) -- Needed for integration
   for ordIndexes in self.phaseQuadRange:colMajorIter() do
      local ordIdx = self.phaseQuadIndexer(ordIndexes)
      self.phaseWeights[ordIdx] = 1.0
      for d = 1, self.numPhaseDims do
	 self.phaseWeights[ordIdx] =
	    self.phaseWeights[ordIdx]*weights[ordIndexes[d]]
	 self.phaseOrdinates[ordIdx][d] = ordinates[ordIndexes[d]]
      end
      self.phaseBasis:evalBasis(self.phaseOrdinates[ordIdx],
				self.phaseBasisAtOrds[ordIdx])
   end

   -- Construct the phase space to conf space ordinate map
   self.phaseToConfOrdMap = Lin.Vec(self.numPhaseOrds)
   for ordIndexes in self.phaseQuadRange:colMajorIter() do
      local confOrdIdx = self.confQuadIndexer(ordIndexes)
      local phaseOrdIdx = self.phaseQuadIndexer(ordIndexes)
      self.phaseToConfOrdMap[phaseOrdIdx] = confOrdIdx
   end

   self.numVelDims = self.numPhaseDims - self.numConfDims
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function MaxwellianOnBasis:_advance(tCurr, inFld, outFld)
   -- Get the inputs and outputs
   local nIn = assert(inFld[1], "MaxwellianOnBasis.advance: Must specify density 'inFld[1]'")
   local uIn = assert(inFld[2], "MaxwellianOnBasis.advance: Must specify drift speed 'inFld[2]'")
   local vth2In = assert(inFld[3], "MaxwellianOnBasis.advance: Must specify thermal velocity squared 'inFld[3]'")
   local fOut = assert(outFld[1], "MaxwellianOnBasis.advance: Must specify an output field 'outFld[1]'")

   local nItr, nOrd = nIn:get(1), Lin.Vec(self.numConfOrds)
   local uItr, uOrd = uIn:get(1), Lin.Mat(self.numConfOrds, self.numVelDims)
   local vth2Itr, vth2Ord = vth2In:get(1), Lin.Vec(self.numConfOrds)
   local fItr = fOut:get(1)

   -- Get the Ranges to loop over the domain
   local confRange = nIn:localRange()
   local confIndexer = nIn:genIndexer()
   local phaseRange = fOut:localRange()
   local phaseIndexer = fOut:genIndexer()
   local l, u = {}, {}
   for d = 1, self.numVelDims do
      l[d] = phaseRange:lower(self.numConfDims + d)
      u[d] = phaseRange:upper(self.numConfDims + d)
   end
   local velRange = Range.Range(l, u)

   -- Additional preallocated variables
   local dz = Lin.Vec(self.numPhaseDims) -- cell shape
   local zc = Lin.Vec(self.numPhaseDims) -- cell center
   local ordIdx = nil
   local phaseIndexes = {}

   -- The conf. space loop
   for confIndexes in confRange:colMajorIter() do
      nIn:fill(confIndexer(confIndexes), nItr)
      uIn:fill(confIndexer(confIndexes), uItr)
      vth2In:fill(confIndexer(confIndexes), vth2Itr)

      -- Evaluate the the primitive variables (given as expansion
      -- coefficiens) on the ordinates
      for ordIndexes in self.confQuadRange:colMajorIter() do
	 ordIdx = self.confQuadIndexer(ordIndexes)
	 nOrd[ordIdx], vth2Ord[ordIdx] = 0.0, 0.0
	 for d = 1, self.numVelDims do uOrd[ordIdx][d] = 0.0 end

	 for k = 1, self.numConfBasis do
	    nOrd[ordIdx] = nOrd[ordIdx] + nItr[k]*self.confBasisAtOrds[ordIdx][k]
	    vth2Ord[ordIdx] = vth2Ord[ordIdx] + vth2Itr[k]*self.confBasisAtOrds[ordIdx][k]
	 end
	 for d = 1, self.numVelDims do
	    for k = 1, self.numConfBasis do
	       uOrd[ordIdx][d] = uOrd[ordIdx][d] +
		  uItr[self.numConfBasis*(d-1)+k]*self.confBasisAtOrds[ordIdx][k]
	    end
	 end
      end

      -- The velocity space loop
      for velIndexes in velRange:colMajorIter() do
	 -- Construct the phase space index ot of the configuration
	 -- space a velocity space indices
	 for d = 1, self.numConfDims do
	    phaseIndexes[d] = confIndexes[d]
	 end
	 for d = 1, self.numVelDims do
	    phaseIndexes[self.numConfDims+d] = velIndexes[d]
	 end
	 fOut:fill(phaseIndexer(phaseIndexes), fItr)

	 -- Get cell shape, cell center coordinates
	 self.phaseGrid:setIndex(phaseIndexes)
	 for d = 1, self.numPhaseDims do dz[d] = self.phaseGrid:dx(d) end
	 self.phaseGrid:cellCenter(zc)

	 ffiC.MaxwellianInnerLoop(nOrd:data(), uOrd:data(), vth2Ord:data(),
				   fItr:data(),
				   self.phaseWeights:data(), dz:data(), zc:data(),
				   self.phaseOrdinates:data(),
				   self.phaseBasisAtOrds:data(),
				   self.phaseToConfOrdMap:data(),
				   self.numPhaseBasis,
				   self.numConfOrds, self.numPhaseOrds,
				   self.numConfDims, self.numPhaseDims)
      end

   end

   -- set id of output to id of projection basis
   fOut:setBasisId(self.phaseBasis:id())
end

return MaxwellianOnBasis
