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

   self._confGrid = assert(tbl.confGrid,
			   "Updater.BgkCollisions: Must provide configuration space grid object using 'confGrid'")
   self._confBasis = assert(tbl.confBasis,
			    "Updater.BgkCollisions: Must provide configuration space basis object using 'confBasis'")
   self._phaseGrid = assert(tbl.phaseGrid,
			    "Updater.BgkCollisions: Must provide phase space grid object using 'phaseGrid'")
   self._phaseBasis = assert(tbl.phaseBasis,
			     "Updater.BgkCollisions: Must provide phase space basis object using 'phaseBasis'")

   self.collFreq = tbl.collFreq
   if not self.collFreq then
      self.mass = assert(
	 tbl.mass, "Updater.BgkCollisions: Must specify mass with 'mass' ('collFreq' is not specified, so classical \nu is used instead)")
      self.charge = assert(
	 tbl.charge, "Updater.BgkCollisions: Must specify charge with 'charge' ('collFreq' is not specified, so classical \nu is used instead)")
      self.elemCharge = assert(
	 tbl.elemCharge, "Updater.BgkCollisions: Must specify elementary charge with 'elemCharge' ('collFreq' is not specified, so classical \nu is used instead)")
      self.epsilon0 = assert(
	 tbl.epsilon0, "Updater.BgkCollisions: Must specify vacuum permitivity with 'epsilon0' ('collFreq' is not specified, so classical \nu is used instead)")
      -- self.coulombLog = assert(
      -- 	 tbl.coulombLog, "Updater.BgkCollisions: Must specify Coulomb logaritm with 'coulombLog' ('collFreq' is not specified, so classical \nu is used instead)")
   end
end

----------------------------------------------------------------------
-- Updater Advance ---------------------------------------------------
function BgkCollisions:_advance(tCurr, inFld, outFld)
   local numPhaseDims = self._phaseGrid:ndim()
   local numConfDims = self._confGrid:ndim()
   local numPhaseBasis = self._phaseBasis:numBasis()

   -- Get the inputs and outputs
   local fIn = assert(inFld[1],
		      "BgkCollisions.advance: Must specify an input distribution function field as input[1]")
   local fMaxwell = assert(inFld[2],
			   "BgkCollisions.advance: Must specify the Maxwellian distribution function field as input[2]")
   local nuFrac = 1.0
   if inFld[3] then
      nuFrac = inFld[3]
   end
   local nIn = nil
   local vth2In = nil
   if not self.collFreq then
      nIn = assert(inFld[4],
		   "BgkCollisions.advance: Must specify an input number density field as input[4] ('collFreq' is not specified, so classical \nu is used instead)")
      vth2In = assert(inFld[5],
		      "BgkCollisions.advance: Must specify an input thermal velocity squared field as input[5] ('collFreq' is not specified, so classical \nu is used instead)")
   end

   local fRhsOut = assert(outFld[1],
		       "BgkCollisions.advance: Must specify an output field")

   local fInItr = fIn:get(1)
   local fMaxwellItr = fMaxwell:get(1)
   local nInItr = nil
   local vth2InItr = nil
   if not self.collFreq then
      nInItr = nIn:get(1)
      vth2InItr = vth2In:get(1)
   end

   local fRhsOutItr = fRhsOut:get(1)

   -- Get the Ranges to loop over the domain
   -- local confRange = numDensityIn:localRange()
   -- local confIndexer = numDensityIn:genIndexer()
   local phaseIndexer = fRhsOut:genIndexer()
   local confIndexer = nil
   if not self.collFreq then
      confIndexer = nIn:genIndexer()
   end
   -- Phase space loop
   local nu = self.collFreq
   local n, vth2 = nil, nil
   local T, coulombLog = nil, nil
   for phaseIdx in fRhsOut:localRangeIter() do
      fIn:fill(phaseIndexer(phaseIdx), fInItr)
      fMaxwell:fill(phaseIndexer(phaseIdx), fMaxwellItr)
      fRhsOut:fill(phaseIndexer(phaseIdx), fRhsOutItr)
      if not self.collFreq then
	 nIn:fill(confIndexer(phaseIdx), nInItr)
	 vth2In:fill(confIndexer(phaseIdx), vth2InItr)
	 n = nInItr[1]/2^numConfDims
	 vth2 = vth2InItr[1]/2^numConfDims
	 if n <= 0. or vth2 <= 0. then 
	    nu = 0.0
	 else
	    T = self.mass*vth2/self.elemCharge
	    if T < 50 then
	       coulombLog = 23.4-1.15*math.log10(n*1e-6)+3.45*math.log10(T)
	    else
	       coulombLog = 25.3-1.15*math.log10(n*1e-6)+2.3*math.log10(T)
	    end
	    nu = self.charge^4/(6*math.sqrt(2)*math.pi^1.5*self.epsilon0^2*self.mass^2) *
	       coulombLog * n / math.sqrt(vth2)^3
	 end
      end

      for k = 1, numPhaseBasis do
	 if fMaxwellItr[k] == fMaxwellItr[k] then -- NaN check
	    fRhsOutItr[k] = fRhsOutItr[k] + nuFrac * nu *
	       (fMaxwellItr[k] - fInItr[k])
	 end
      end
   end
end

return BgkCollisions
