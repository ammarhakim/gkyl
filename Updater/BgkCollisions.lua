-- Gkyl ------------------------------------------------------------------------
--
-- Updater to calculate BGK collisions for various species
--
--------------------------------------------------------------------------------

-- Gkyl libraries
local GaussQuadRules = require "Lib.GaussQuadRules"
local Lin = require "Lib.Linalg"
local Proto = require "Proto"
local Range = require "Lib.Range"
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
local BgkCollisions = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function BgkCollisions:init(tbl)
   BgkCollisions.super.init(self, tbl) -- setup base object

   self._cfl = assert(tbl.cfl, "Updater.BgkCollisions: Must provide CFL number using 'cfl'")
   self._speciesList = assert(tbl.speciesList,
			      "Updater.BgkCollisions: Must provide a list of species using 'speciesList'")
   self._collFreq = assert(tbl.collFreq,
			   "Updater.BgkCollisions: Must provide a collision frequency using 'collFreq'")
   self._confGrid = assert(tbl.confGrid,
			   "Updater.BgkCollisions: Must provide configuration space grid object using 'confGrid'")
   self._confBasis = assert(tbl.confBasis,
			    "Updater.BgkCollisions: Must provide configuration space basis object using 'confBasis'")
   self._phaseGrid = assert(tbl.phaseGrid,
			    "Updater.BgkCollisions: Must provide phase space grid object using 'phaseGrid'")
   self._phaseBasis = assert(tbl.phaseBasis,
			     "Updater.BgkCollisions: Must provide phase space basis object using 'phaseBasis'")

   assert(#self._speciesList == #self._collFreq,
	  "List of species table and collision frequences table must have the same size")

   -- Number of quadrature points in each direction
   self._N = tbl.numConfQuad and tbl.numConfQuad or self._confBasis:polyOrder() + 1

   -- 1D weights and ordinates
   local ordinates = GaussQuadRules.ordinates[self._N]
   local weights = GaussQuadRules.weights[self._N]

   local numConfDims = self._confBasis:ndim()
   local l, u = {}, {}
   for d = 1, numConfDims do l[d], u[d] = 1, self._N end
   local quadRange = Range.Range(l, u) -- for looping over quadrature nodes
   local numOrdinates = quadRange:volume() -- number of ordinates
   -- construct weights and ordinates for integration in multiple dimensions
   self._confOrdinates = Lin.Mat(numOrdinates, numConfDims)
   local nodeNum = 1
   for idx in quadRange:colMajorIter() do
      for d = 1, numConfDims do
	 self._confOrdinates[nodeNum][d] = ordinates[idx[d]]
      end
      nodeNum = nodeNum + 1
   end
   local numBasis = self._confBasis:numBasis()
   self._confBasisAtOrdinates = Lin.Mat(numOrdinates, numBasis)
   -- pre-compute values of basis functions at quadrature nodes
   for n = 1, numOrdinates do
      self._confBasis:evalBasis(self._confOrdinates[n],
				self._confBasisAtOrdinates[n])
   end

   local numPhaseDims = nil
   for _, nm in ipairs(self._speciesList) do
      if numPhaseDims == nil then
	 numPhaseDims = self._phaseBasis:ndim()
      else
	 assert(numPhaseDims == self._phaseBasis[nm]:ndim(),
		"Phase spaces for different species must have the same number of dimensions")
      end
   end
   for d = 1, numPhaseDims do l[d], u[d] = 1, self._N end
   quadRange = Range.Range(l, u) -- for looping over quadrature nodes
   numOrdinates = quadRange:volume() -- number of ordinates
   -- construct weights and ordinates for integration in multiple dimensions
   self._phaseOrdinates = Lin.Mat(numOrdinates, numPhaseDims)
   self._phaseWeights = Lin.Vec(numOrdinates)
   nodeNum = 1
   for idx in quadRange:colMajorIter() do
      self._phaseWeights[nodeNum] = 1.0
      for d = 1, numPhaseDims do
	 self._phaseWeights[nodeNum] =
	    self._phaseWeights[nodeNum] * weights[idx[d]]
	 self._phaseOrdinates[nodeNum][d] = ordinates[idx[d]]
      end
      nodeNum = nodeNum + 1
   end
   numBasis = self._phaseBasis:numBasis()
   self._phaseBasisAtOrdinates = Lin.Mat(numOrdinates, numBasis)
   -- pre-compute values of basis functions at quadrature nodes
   for n = 1, numOrdinates do
      self._phaseBasis:evalBasis(self._phaseOrdinates[n],
				 self._phaseBasisAtOrdinates[n])
   end

   -- construct various functions from template representations
   self._compToPhys = loadstring(compToPhysTempl {NDIM = numPhaseDims} )()
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function BgkCollisions:_advance(tCurr, dt, inFld, outFld)
   local numConfDims = self._confGrid:ndim()
   local numConfBasis = self._confBasis:numBasis()
   local numPhaseDims = self._phaseGrid:ndim()
   local numPhaseBasis = self._phaseBasis:numBasis()
   local numVelDims = numPhaseDims - numConfDims

   -- Define regions for looping over ordinates
   local l, u = {}, {}
   for d = 1, numConfDims do l[d], u[d] = 1, self._N end
   local confQuadRange = Range.Range(l, u)
   local confQuadIndexer = Range.makeColMajorGenIndexer(confQuadRange)
   for d = 1, numPhaseDims do l[d], u[d] = 1, self._N end
   local phaseQuadRange = Range.Range(l, u)
   local phaseQuadIndexer = Range.makeColMajorGenIndexer(phaseQuadRange)
   local phaseIdx = {}

   -- Variables for evaluating the moments at quadrature points
   local numConfOrdinates = confQuadRange:volume()
   local numDensityOrd = Lin.Vec(numConfOrdinates)
   local momDensityOrd = Lin.Mat(numConfOrdinates, numVelDims)
   local ptclEnergyOrd = Lin.Vec(numConfOrdinates)
   local velBulkOrd = Lin.Mat(numConfOrdinates, numVelDims)
   local velTherm2Ord = Lin.Vec(numConfOrdinates)

   -- Variables to get the physical coordinates from of the
   -- computational space
   local dz = Lin.Vec(numPhaseDims) -- cell shape
   local zc = Lin.Vec(numPhaseDims) -- cell center
   local zPhys = Lin.Vec(numPhaseDims)

   -- Additional variables which we don't want to allocate inside the
   -- loops
   local maxwellian = nil
   local maxwellianModal = Lin.Vec(numPhaseBasis)
   local u2 = nil
   local confMu = nil
   local phaseMu = nil
   local offset = nil

   -- Main loop over the species
   for i, nm in ipairs(self._speciesList) do
      -- Get the inputs and outputs
      local numDensityIn = assert(inFld[nm][1],
				  "BgkCollisions.advance: Must specify an input fluid moments field")
      local momDensityIn = assert(inFld[nm][2],
				  "BgkCollisions.advance: Must specify an input fluid moments field")
      local ptclEnergyIn = assert(inFld[nm][3],
				  "BgkCollisions.advance: Must specify an input fluid moments field")
      local fOut = assert(outFld[nm],
			  "BgkCollisions.advance: Must specify an output field")
      local ndItr = numDensityIn:get(1)
      local mdItr = momDensityIn:get(1)
      local peItr = ptclEnergyIn:get(1)
      local fItr = fOut:get(1)

      -- Get the Ranges to loop over the domain
      local confRange = numDensityIn:localRange()
      local confIndexer = numDensityIn:genIndexer()
      local phaseRange = fOut:localRange()
      local phaseIndexer = fOut:genIndexer()
      l, u = {}, {}
      for d = 1, numVelDims do
	 l[d] = phaseRange:lower(numConfDims + d)
	 u[d] = phaseRange:upper(numConfDims + d)
      end
      local velRange = Range.Range(l, u)

      local nu = self._collFreq[i]  -- Collision frequency function

      -- Configuration space loop
      for confIdx in confRange:colMajorIter() do
	 numDensityIn:fill(confIndexer(confIdx), ndItr)
	 momDensityIn:fill(confIndexer(confIdx), mdItr)
	 ptclEnergyIn:fill(confIndexer(confIdx), peItr)

	 -- Evaluate the the moments (given as expansion coefficiens)
	 -- on the ordinates
	 for muIdx in confQuadRange:colMajorIter() do
	    confMu = confQuadIndexer(muIdx)
	    numDensityOrd[confMu] = 0
	    for d = 1, numVelDims do
	       momDensityOrd[confMu][d] = 0
	    end
	    ptclEnergyOrd[confMu] = 0
	    for k = 1, numConfBasis do
	       numDensityOrd[confMu] = numDensityOrd[confMu] +
		  ndItr[k] * self._confBasisAtOrdinates[confMu][k]
	       ptclEnergyOrd[confMu] = ptclEnergyOrd[confMu] +
		  peItr[k] * self._confBasisAtOrdinates[confMu][k]
	    end
	    offset = 0
	    for d = 1, numVelDims do
	       for k = 1, numConfBasis do
		  momDensityOrd[confMu][d] = momDensityOrd[confMu][d] +
		     mdItr[offset + k] * self._confBasisAtOrdinates[confMu][k]
	       end
	       offset = offset + numConfBasis
	    end

	    -- Calculate the bulk velocity out of the flux and thermal
	    -- velocity out of the energy
	    u2 = 0
	    for d = 1, numVelDims do
	       velBulkOrd[confMu][d] = momDensityOrd[confMu][d] / numDensityOrd[confMu]
	       u2 = u2 + velBulkOrd[confMu][d] * velBulkOrd[confMu][d]
	    end
	    velTherm2Ord[confMu] = ptclEnergyOrd[confMu] / numDensityOrd[confMu] - u2
	 end

	 -- Velocity space loop
	 for velIdx in velRange:colMajorIter() do
	    -- Construct the phase space index ot of the configuration
	    -- space a velocity space indices
	    for d = 1, numConfDims do phaseIdx[d] = confIdx[d] end
	    for d = 1, numVelDims do phaseIdx[d + numConfDims] = velIdx[d] end
	    fOut:fill(phaseIndexer(phaseIdx), fItr)

	    -- Get cell shape, cell center coordinates
	    self._phaseGrid:setIndex(phaseIdx)
	    for d = 1, numPhaseDims do dz[d] = self._phaseGrid:dx(d) end
	    self._phaseGrid:cellCenter(zc)

	    for k = 1, numPhaseBasis do maxwellianModal[k] = 0 end
	    for muIdx in phaseQuadRange:colMajorIter() do
	       confMu = confQuadIndexer(muIdx)
	       phaseMu = phaseQuadIndexer(muIdx)

	       -- Get the physical velocity coordinates
	       self._compToPhys(self._phaseOrdinates[phaseMu], dz, zc, zPhys)

	       -- Constract the Maxwellian
	       maxwellian = numDensityOrd[confMu]
	       for d = 1, numVelDims do
		  maxwellian = maxwellian / math.sqrt(2 * math.pi * velTherm2Ord[confMu]) *
		     math.exp(-(zPhys[numConfDims + d] - velBulkOrd[confMu][d]) *
				 (zPhys[numConfDims + d] - velBulkOrd[confMu][d]) /
				 (2 * velTherm2Ord[confMu]))
	       end
	       -- Project the Maxwellian on basis
	       for k = 1, numPhaseBasis do
		  maxwellianModal[k] = maxwellianModal[k] +
		     self._phaseWeights[phaseMu] * maxwellian * self._phaseBasisAtOrdinates[phaseMu][k]
	       end
	    end
	    -- Modify the solution with BGK collisions
	    for k = 1, numPhaseBasis do
	       fItr[k] = fItr[k] + self._cfl * dt * nu(1, 1) *
		  (maxwellianModal[k] - fItr[k])
	    end
	 end
      end
   end
   return true, GKYL_MAX_DOUBLE
end

return BgkCollisions
