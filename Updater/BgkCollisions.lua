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
local Time = require "Lib.Time"
local UpdaterBase = require "Updater.Base"

-- System libraries
local xsys = require "xsys"

-- BGK Collisions updater object
local BgkCollisions = Proto(UpdaterBase)

----------------------------------------------------------------------
-- Updater Initialization --------------------------------------------
function BgkCollisions:init(tbl)
   BgkCollisions.super.init(self, tbl) -- setup base object

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

   -- -- Number of quadrature points in each direction
   -- self._N = tbl.numConfQuad and tbl.numConfQuad or self._confBasis:polyOrder() + 1

   -- -- 1D weights and ordinates
   -- local ordinates = GaussQuadRules.ordinates[self._N]
   -- local weights = GaussQuadRules.weights[self._N]

   -- local numConfDims = self._confBasis:ndim()
   -- local l, u = {}, {}
   -- for d = 1, numConfDims do l[d], u[d] = 1, self._N end
   -- local quadRange = Range.Range(l, u) -- for looping over quadrature nodes
   -- local numOrdinates = quadRange:volume() -- number of ordinates
   -- -- construct weights and ordinates for integration in multiple dimensions
   -- self._confOrdinates = Lin.Mat(numOrdinates, numConfDims)
   -- local nodeNum = 1
   -- for idx in quadRange:colMajorIter() do
   --    for d = 1, numConfDims do
   -- 	 self._confOrdinates[nodeNum][d] = ordinates[idx[d]]
   --    end
   --    nodeNum = nodeNum + 1
   -- end
   -- local numBasis = self._confBasis:numBasis()
   -- self._confBasisAtOrdinates = Lin.Mat(numOrdinates, numBasis)
   -- -- pre-compute values of basis functions at quadrature nodes
   -- for n = 1, numOrdinates do
   --    self._confBasis:evalBasis(self._confOrdinates[n],
   -- 				self._confBasisAtOrdinates[n])
   -- end

   -- numPhaseDims = self._phaseBasis:ndim()
   -- for d = 1, numPhaseDims do l[d], u[d] = 1, self._N end
   -- quadRange = Range.Range(l, u) -- for looping over quadrature nodes
   -- numOrdinates = quadRange:volume() -- number of ordinates
   -- -- construct weights and ordinates for integration in multiple dimensions
   -- self._phaseOrdinates = Lin.Mat(numOrdinates, numPhaseDims)
   -- self._phaseWeights = Lin.Vec(numOrdinates)
   -- nodeNum = 1
   -- for idx in quadRange:colMajorIter() do
   --    self._phaseWeights[nodeNum] = 1.0
   --    for d = 1, numPhaseDims do
   -- 	 self._phaseWeights[nodeNum] =
   -- 	    self._phaseWeights[nodeNum] * weights[idx[d]]
   -- 	 self._phaseOrdinates[nodeNum][d] = ordinates[idx[d]]
   --    end
   --    nodeNum = nodeNum + 1
   -- end
   -- numBasis = self._phaseBasis:numBasis()
   -- self._phaseBasisAtOrdinates = Lin.Mat(numOrdinates, numBasis)
   -- -- pre-compute values of basis functions at quadrature nodes
   -- for n = 1, numOrdinates do
   --    self._phaseBasis:evalBasis(self._phaseOrdinates[n],
   -- 				 self._phaseBasisAtOrdinates[n])
   -- end

   -- timings
   self._tmEvalMom = 0.0
   self._tmProjectMaxwell = 0.0
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function BgkCollisions:_advance(tCurr, dt, inFld, outFld)
   -- local numConfDims = self._confGrid:ndim()
   -- local numConfBasis = self._confBasis:numBasis()
   local numPhaseDims = self._phaseGrid:ndim()
   local numPhaseBasis = self._phaseBasis:numBasis()

   -- Define regions for looping over ordinates
   -- local l, u = {}, {}
   -- for d = 1, numConfDims do l[d], u[d] = 1, self._N end
   -- local confQuadRange = Range.Range(l, u)
   -- local confQuadIndexer = Range.makeColMajorGenIndexer(confQuadRange)
   -- for d = 1, numPhaseDims do l[d], u[d] = 1, self._N end
   -- local phaseQuadRange = Range.Range(l, u)
   -- local phaseQuadIndexer = Range.makeColMajorGenIndexer(phaseQuadRange)
   -- local phaseIdx = {}

   -- -- Variables for evaluating the moments at quadrature points
   -- local numConfOrdinates = confQuadRange:volume()
   -- local numDensityOrd = Lin.Vec(numConfOrdinates)
   -- local momDensityOrd = Lin.Mat(numConfOrdinates, numVelDims)
   -- local ptclEnergyOrd = Lin.Vec(numConfOrdinates)
   -- local velBulkOrd = Lin.Mat(numConfOrdinates, numVelDims)
   -- local velTherm2Ord = Lin.Vec(numConfOrdinates)

   -- -- Variables to get the physical coordinates from of the
   -- -- computational space
   -- local dz = Lin.Vec(numPhaseDims) -- cell shape
   -- local zc = Lin.Vec(numPhaseDims) -- cell center
   -- local zPhys = Lin.Vec(numPhaseDims)

   -- -- Additional variables which we don't want to allocate inside the
   -- -- loops
   -- local maxwellian = nil
   -- local maxwellianModal = Lin.Vec(numPhaseBasis)
   -- local u2 = nil
   -- local confMu = nil
   -- local phaseMu = nil
   -- local offset = nil

   -- Get the inputs and outputs
   local fIn = assert(inFld[1],
		      "BgkCollisions.advance: Must specify an input distribution function field as input[1]")
   local fMaxwell = assert(inFld[2],
			   "BgkCollisions.advance: Must specify the Maxwellian distribution function field as input[2]")
   local fOut = assert(outFld[1],
		       "BgkCollisions.advance: Must specify an output field")
   
   local fInItr = fIn:get(1)
   local fMaxwellItr = fMaxwell:get(1)
   local fOutItr = fOut:get(1)

   -- Get the Ranges to loop over the domain
   -- local confRange = numDensityIn:localRange()
   -- local confIndexer = numDensityIn:genIndexer()
   local phaseRange = fOut:localRange()
   local phaseIndexer = fOut:genIndexer()

   local nu = self._collFreq  -- Collision frequency function

   -- Phase space loop
   for phaseIdx in phaseRange:colMajorIter() do
      fIn:fill(phaseIndexer(phaseIdx), fInItr)
      fMaxwell:fill(phaseIndexer(phaseIdx), fMaxwellItr)
      fOut:fill(phaseIndexer(phaseIdx), fOutItr)

      for k = 1, numPhaseBasis do
	 fOutItr[k] = fOutItr[k] + dt * nu(1, 1) * -- FIX collision frequency!
	    (fMaxwellItr[k] - fInItr[k])
      end
   end
   return true, GKYL_MAX_DOUBLE
end

function BgkCollisions:evalMomTime() return self._tmEvalMom end
function BgkCollisions:projectMaxwellTime() return self._tmProjectMaxwell end

return BgkCollisions
