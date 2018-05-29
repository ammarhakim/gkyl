-- Gkyl ------------------------------------------------------------------------
--
-- Updater to create the Maxwellian distribution from the coserved
-- moments and project it on basis functions. Uses Gaussian
-- quadrature.
--
--------------------------------------------------------------------------------

-- Gkyl libraries
local DataStruct = require "DataStruct"
local GaussQuadRules = require "Lib.GaussQuadRules"
local Lin = require "Lib.Linalg"
local Proto = require "Proto"
local Range = require "Lib.Range"
local Time = require "Lib.Time"
local CartFieldBinOp = require "Updater.CartFieldBinOp"
local UpdaterBase = require "Updater.Base"

-- System libraries
local xsys = require "xsys"

-- Template for function to map computional space -> physical space
local compToPhysTempl = xsys.template([[
return function (eta, dx, xc, xOut)
|for i = 1, NDIM do
   xOut[${i}] = 0.5*dx[${i}]*eta[${i}] + xc[${i}]
|end
end
]])

-- BGK Collisions updater object
local MaxwellianOnBasis = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function MaxwellianOnBasis:init(tbl)
   MaxwellianOnBasis.super.init(self, tbl) -- setup base object

   self.confGrid = assert(tbl.confGrid,
			  "Updater.MaxwellianOnBasis: Must provide configuration space grid object using 'confGrid'")
   self.confBasis = assert(tbl.confBasis,
			   "Updater.MaxwellianOnBasis: Must provide configuration space basis object using 'confBasis'")
   self.phaseGrid = assert(tbl.phaseGrid,
			   "Updater.MaxwellianOnBasis: Must provide phase space grid object using 'phaseGrid'")
   self.phaseBasis = assert(tbl.phaseBasis,
			    "Updater.MaxwellianOnBasis: Must provide phase space basis object using 'phaseBasis'")

   -- Number of quadrature points in each direction
   local N = tbl.numConfQuad and tbl.numConfQuad or self.confBasis:polyOrder() + 1
   assert(N<=8, "Updater.MaxwellianOnBasis: Gaussian quadrature only implemented for numQuad<=8 in each dimension")

   self.projectOnGhosts = xsys.pickBool(tbl.projectOnGhosts, false)
   self.numVelDims = self.phaseGrid:ndim() - self.confGrid:ndim()

   -- 1D weights and ordinates
   local ordinates = GaussQuadRules.ordinates[N]
   local weights = GaussQuadRules.weights[N]

   -- Conf. space ordinates
   self.numConfDims = self.confBasis:ndim()
   local l, u = {}, {}
   for d = 1, self.numConfDims do l[d], u[d] = 1, N end
   self.confQuadRange = Range.Range(l, u)
   self.numConfOrdinates = self.confQuadRange:volume()
   self.confOrdinates = Lin.Mat(self.numConfOrdinates, self.numConfDims)
   local nodeNum = 1
   for idx in self.confQuadRange:colMajorIter() do
      for d = 1, self.numConfDims do
	 self.confOrdinates[nodeNum][d] = ordinates[idx[d]]
      end
      nodeNum = nodeNum + 1
   end
   self.numConfBasis = self.confBasis:numBasis()
   self.confBasisAtOrdinates = Lin.Mat(self.numConfOrdinates, self.numConfBasis)
   -- Pre-compute values of basis functions at quadrature nodes
   for n = 1, self.numConfOrdinates do
      self.confBasis:evalBasis(self.confOrdinates[n],
				self.confBasisAtOrdinates[n])
   end

   -- Phase space ordinates and weights
   self.numPhaseDims = self.phaseBasis:ndim()
   for d = 1, self.numPhaseDims do l[d], u[d] = 1, N end
   self.phaseQuadRange = Range.Range(l, u)
   self.numPhaseOrdinates = self.phaseQuadRange:volume()
   self.phaseOrdinates = Lin.Mat(self.numPhaseOrdinates,
				 self.numPhaseDims)
   self.phaseWeights = Lin.Vec(self.numPhaseOrdinates)
   nodeNum = 1
   for idx in self.phaseQuadRange:colMajorIter() do
      self.phaseWeights[nodeNum] = 1.0
      for d = 1, self.numPhaseDims do
	 self.phaseWeights[nodeNum] =
	    self.phaseWeights[nodeNum] * weights[idx[d]]
	 self.phaseOrdinates[nodeNum][d] = ordinates[idx[d]]
      end
      nodeNum = nodeNum + 1
   end
   self.numPhaseBasis = self.phaseBasis:numBasis()
   self.phaseBasisAtOrdinates = Lin.Mat(self.numPhaseOrdinates,
					self.numPhaseBasis)
   -- Pre-compute values of basis functions at quadrature nodes
   for n = 1, self.numPhaseOrdinates do
      self.phaseBasis:evalBasis(self.phaseOrdinates[n],
				self.phaseBasisAtOrdinates[n])
   end

   -- Construct the logical space to physical space mapping
   self.compToPhys = loadstring(compToPhysTempl {NDIM = self.numPhaseDims} )()

   -- Weak division of two configuration space fields
   self.confDiv = CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "Divide",
   }
   -- Dot product of two configuration space vector fields
   self.confDotProduct = CartFieldBinOp {
      onGrid = self.confGrid,
      weakBasis = self.confBasis,
      operation = "DotProduct",
   }

   -- Flow velocity in vdim directions
   self.bulkVelocity = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.numConfBasis*self.numVelDims,
      ghost = {1, 1},
   }
   -- Thermal velocity squared; vth=sqrt(T/m)
   self.thermVelocity2 = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.numConfBasis,
      ghost = {1, 1},
   }
   -- Magnitude of kinetic energy density vector
   self.kinEnergyDens = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.numConfBasis,
      ghost = {1, 1},
   }
   self.thermEnergyDens = DataStruct.Field {
      onGrid = self.confGrid,
      numComponents = self.numConfBasis,
      ghost = {1, 1},
   }
end

-- Computes primitive variables from the moments
function MaxwellianOnBasis:primVariables(M0, M1i, M2)
   self.confDiv:advance(0., 0., {M0, M1i}, {self.bulkVelocity})
   self.confDotProduct:advance(0., 0., {self.bulkVelocity, M1i},
			       {self.kinEnergyDens})
   self.thermEnergyDens:combine(1.0/self.numVelDims, M2,
				   -1.0/self.numVelDims, self.kinEnergyDens)
   self.confDiv:advance(0., 0., {M0, self.thermEnergyDens},
			{self.thermVelocity2})
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function MaxwellianOnBasis:_advance(tCurr, dt, inFld, outFld)
   -- Define iterrators for looping over ordinates
   local confQuadIndexer = Range.makeColMajorGenIndexer(self.confQuadRange)
   local phaseQuadIndexer = Range.makeColMajorGenIndexer(self.phaseQuadRange)
   local phaseIdx = {}

   -- Variables for evaluating the moments at quadrature points
   local numDensityOrd = Lin.Vec(self.numConfOrdinates)
   local bulkVelocityOrd = Lin.Mat(self.numConfOrdinates, self.numVelDims)
   local thermVelocity2Ord = Lin.Vec(self.numConfOrdinates)

   -- Variables to get the physical coordinates from of the
   -- computational space
   local dz = Lin.Vec(self.numPhaseDims) -- cell shape
   local zc = Lin.Vec(self.numPhaseDims) -- cell center
   local zPhys = Lin.Vec(self.numPhaseDims)

   -- Additional variables which we don't want to allocate inside the
   -- loops
   local maxwellian = nil
   local confMu = nil
   local phaseMu = nil
   local offset = nil

   -- Get the inputs and outputs
   assert(inFld[1], "MaxwellianOnBasis.advance: Must specify M0 input fluid moments field")
   assert(inFld[2], "MaxwellianOnBasis.advance: Must specify M1i input fluid moments field")
   assert(inFld[3], "MaxwellianOnBasis.advance: Must specify M2 input fluid moments field")
   -- load the moments and calculate the primitive variables
   self.numDensity = inFld[1]
   self:primVariables(inFld[1], inFld[2], inFld[3])

   local fOut = assert(outFld[1], "MaxwellianOnBasis.advance: Must specify an output field")
   local numDensityItr = self.numDensity:get(1)
   local bulkVelocityItr = self.bulkVelocity:get(1)
   local thermVelocity2Itr = self.thermVelocity2:get(1)
   local fItr = fOut:get(1)

   -- Get the Ranges to loop over the domain
   local confRange = self.numDensity:localRange()
   local confIndexer = self.numDensity:genIndexer()
   local phaseRange = fOut:localRange()
   local phaseIndexer = fOut:genIndexer()
   local l, u = {}, {}
   for d = 1, self.numVelDims do
      l[d] = phaseRange:lower(self.numConfDims + d)
      u[d] = phaseRange:upper(self.numConfDims + d)
   end
   local velRange = Range.Range(l, u)

   -- The conf. space loop
   for confIdx in confRange:colMajorIter() do
      self.numDensity:fill(confIndexer(confIdx), numDensityItr)
      self.bulkVelocity:fill(confIndexer(confIdx), bulkVelocityItr)
      self.thermVelocity2:fill(confIndexer(confIdx), thermVelocity2Itr)

      -- Evaluate the the primitive variables (given as expansion
      -- coefficiens) on the ordinates
      for muIdx in self.confQuadRange:colMajorIter() do
	 confMu = confQuadIndexer(muIdx)
	 numDensityOrd[confMu] = 0
	 for d = 1, self.numVelDims do bulkVelocityOrd[confMu][d] = 0 end
	 thermVelocity2Ord[confMu] = 0

	 for k = 1, self.numConfBasis do
	    numDensityOrd[confMu] = numDensityOrd[confMu] +
	       numDensityItr[k] * self.confBasisAtOrdinates[confMu][k]
	    thermVelocity2Ord[confMu] = thermVelocity2Ord[confMu] +
		  thermVelocity2Itr[k] * self.confBasisAtOrdinates[confMu][k]
	 end
	 offset = 0
	 for d = 1, self.numVelDims do
	    for k = 1, self.numConfBasis do
	       bulkVelocityOrd[confMu][d] = bulkVelocityOrd[confMu][d] +
		  bulkVelocityItr[offset + k] *
		  self.confBasisAtOrdinates[confMu][k]
	    end
	    offset = offset + self.numConfBasis
	 end
      end

      -- The velocity space loop
      for velIdx in velRange:colMajorIter() do
	 -- Construct the phase space index ot of the configuration
	 -- space a velocity space indices
	 for d = 1, self.numConfDims do phaseIdx[d] = confIdx[d] end
	 for d = 1, self.numVelDims do
	    phaseIdx[d + self.numConfDims] = velIdx[d]
	 end
	 fOut:fill(phaseIndexer(phaseIdx), fItr)

	 -- Get cell shape, cell center coordinates
	 self.phaseGrid:setIndex(phaseIdx)
	 for d = 1, self.numPhaseDims do dz[d] = self.phaseGrid:dx(d) end
	 self.phaseGrid:cellCenter(zc)

	 for k = 1, self.numPhaseBasis do fItr[k] = 0 end
	 for muIdx in self.phaseQuadRange:colMajorIter() do
	    confMu = confQuadIndexer(muIdx)
	    phaseMu = phaseQuadIndexer(muIdx)

	    -- Get the physical velocity coordinates
	    self.compToPhys(self.phaseOrdinates[phaseMu], dz, zc, zPhys)

	    -- Construct the Maxwellian
	    maxwellian = numDensityOrd[confMu]
	    for d = 1, self.numVelDims do
	       maxwellian = maxwellian /
		  math.sqrt(2*math.pi * thermVelocity2Ord[confMu]) *
		  math.exp(-(zPhys[self.numConfDims+d]-bulkVelocityOrd[confMu][d]) *
			      (zPhys[self.numConfDims+d]-bulkVelocityOrd[confMu][d]) /
			      (2 * thermVelocity2Ord[confMu]))
	    end
	    -- Project the Maxwellian on basis
	    for k = 1, self.numPhaseBasis do
	       fItr[k] = fItr[k] +
		  self.phaseWeights[phaseMu] * maxwellian *
		  self.phaseBasisAtOrdinates[phaseMu][k]
	    end
	 end
      end
   end
   return true, GKYL_MAX_DOUBLE
end

return MaxwellianOnBasis
